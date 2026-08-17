//! Throttled progress reporting for the long unattended sweeps (design §12.1 extension, Q19).

use rayon::prelude::*;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::{Duration, Instant};

/// Time-based emit throttle: at most one progress line per this interval.
pub const EMIT_INTERVAL: Duration = Duration::from_secs(15);

/// Whole microseconds in `d`, saturating (a sweep would have to run ~584 000 years to overflow).
fn micros(d: Duration) -> u64 {
    u64::try_from(d.as_micros()).unwrap_or(u64::MAX)
}

/// Render a duration at the scale a human reads it at: `9s`, `16m05s`, `55h04m`.
fn fmt_duration(d: Duration) -> String {
    let secs = d.as_secs();
    match secs {
        0..60 => format!("{secs}s"),
        60..3600 => format!("{}m{:02}s", secs / 60, secs % 60),
        _ => format!("{}h{:02}m", secs / 3600, (secs % 3600) / 60),
    }
}

/// Time remaining, extrapolated linearly from the rate achieved so far:
/// `elapsed × (total − completed) / completed`. `None` before the first case finishes (no rate
/// yet). An *estimate* — cases are not equal-cost, so a sweep that starts on its dilute (fast)
/// rows will overshoot early and settle as the expensive rows land.
fn eta(completed: usize, total: usize, elapsed: Duration) -> Option<Duration> {
    if completed == 0 {
        return None;
    }
    let remaining = total.saturating_sub(completed);
    let secs = elapsed.as_secs_f64() * remaining as f64 / completed as f64;
    Duration::try_from_secs_f64(secs).ok()
}

/// A shared completion counter for a rayon sweep, emitting a throttled progress line to stderr.
#[derive(Debug)]
pub struct Progress {
    label: &'static str,
    total: usize,
    completed: AtomicUsize,
    /// Elapsed micros at the last emitted line; the throttle's only mutable state.
    last_emit: AtomicU64,
    start: Instant,
}

impl Progress {
    pub fn new(label: &'static str, total: usize) -> Self {
        Self {
            label,
            total,
            completed: AtomicUsize::new(0),
            last_emit: AtomicU64::new(0),
            start: Instant::now(),
        }
    }

    /// Record one finished case against the wall clock, writing any throttled line to **stderr**
    /// so it never contaminates a stdout artifact.
    pub fn complete(&self) {
        if let Some(line) = self.complete_at(self.start.elapsed()) {
            eprintln!("  {line}");
        }
    }

    /// The closing line: unconditional, so a finished sweep always states its own wall time.
    pub fn finish(&self) {
        eprintln!(
            "  {}: done {}/{} in {}",
            self.label,
            self.completed(),
            self.total,
            fmt_duration(self.start.elapsed()),
        );
    }

    /// Cases finished so far.
    pub fn completed(&self) -> usize {
        self.completed.load(Ordering::Relaxed)
    }

    /// Record one finished case at `elapsed` since the sweep started, returning the progress line
    /// to emit if the throttle allows one. The clock seam: `complete` supplies the real elapsed.
    pub fn complete_at(&self, elapsed: Duration) -> Option<String> {
        let completed = self.completed.fetch_add(1, Ordering::Relaxed) + 1;
        let now = micros(elapsed);
        let last = self.last_emit.load(Ordering::Relaxed);
        if now < last.saturating_add(micros(EMIT_INTERVAL)) {
            return None;
        }
        // Compare-and-swap so exactly one racing thread wins the emit slot for this interval.
        self.last_emit
            .compare_exchange(last, now, Ordering::Relaxed, Ordering::Relaxed)
            .ok()?;
        let pct = 100.0 * completed as f64 / self.total as f64;
        let remaining =
            eta(completed, self.total, elapsed).map_or_else(|| "?".to_string(), fmt_duration);
        Some(format!(
            "{}: {completed}/{} ({pct:.0}%) elapsed {}, eta ~{remaining}",
            self.label,
            self.total,
            fmt_duration(elapsed),
        ))
    }
}

/// Map `cases` in parallel (rayon) with throttled progress on stderr, **preserving input order**
/// so the JSONL rows land deterministically regardless of the parallel schedule.
pub fn par_map_with_progress<T, R, F>(label: &'static str, cases: &[T], f: F) -> Vec<R>
where
    T: Sync,
    R: Send,
    F: Fn(&T) -> R + Sync + Send,
{
    let progress = Progress::new(label, cases.len());
    eprintln!("  {label}: {} cases", cases.len());
    let out = cases
        .par_iter()
        .map(|case| {
            let r = f(case);
            progress.complete();
            r
        })
        .collect();
    progress.finish();
    out
}

#[cfg(test)]
mod tests {
    use super::{EMIT_INTERVAL, Progress, eta, fmt_duration, par_map_with_progress};
    use rayon::prelude::*;
    use std::time::Duration;

    /// The retrofit wrapper every sweep calls. Determinism of the JSONL artifacts rests on the
    /// rayon map preserving input order (each case is independent, so the schedule must not show
    /// up in the output) — instrumenting it must not disturb that.
    #[test]
    fn par_map_preserves_input_order() {
        let cases: Vec<usize> = (0..500).collect();
        let out = par_map_with_progress("order", &cases, |&i| format!("case-{i}"));

        let expected: Vec<String> = (0..500).map(|i| format!("case-{i}")).collect();
        assert_eq!(out, expected);
    }

    /// The counter is shared across rayon's worker threads, so it must count every case exactly
    /// once and hand the emit slot to exactly one racing thread — a double-count would misreport
    /// progress and a double-emit would interleave two half-written lines on stderr.
    #[test]
    fn counts_and_emits_exactly_once_under_rayon() {
        const N: usize = 2_000;

        // All inside the first interval: every case counted, nothing emitted.
        let quiet = Progress::new("parallel", N);
        (0..N)
            .into_par_iter()
            .for_each(|_| assert!(quiet.complete_at(Duration::ZERO).is_none()));
        assert_eq!(quiet.completed(), N);

        // All past the boundary: still every case counted, and exactly one thread emits.
        let noisy = Progress::new("parallel", N);
        let lines: Vec<String> = (0..N)
            .into_par_iter()
            .filter_map(|_| noisy.complete_at(EMIT_INTERVAL))
            .collect();
        assert_eq!(noisy.completed(), N);
        assert_eq!(lines.len(), 1, "expected one line, got {lines:?}");
    }

    /// The emitted line is what Seth reads off a background run, so its content is the contract:
    /// which sweep, how far along, how long it has been, and how much longer — the last marked
    /// `~` because it is an extrapolation, not a schedule.
    #[test]
    fn line_reports_count_percent_elapsed_and_eta() {
        let p = Progress::new("heavyplate", 100);
        for _ in 0..24 {
            assert_eq!(p.complete_at(Duration::ZERO), None);
        }
        assert_eq!(
            p.complete_at(Duration::from_secs(60)).as_deref(),
            Some("heavyplate: 25/100 (25%) elapsed 1m00s, eta ~3m00s")
        );
    }

    /// Durations render at the scale a human reads them at: a 16-minute sweep in minutes, a
    /// multi-hour one in hours. Worked by hand.
    #[test]
    fn durations_render_at_human_scale() {
        assert_eq!(fmt_duration(Duration::from_secs(9)), "9s");
        assert_eq!(fmt_duration(Duration::from_secs(59)), "59s");
        assert_eq!(fmt_duration(Duration::from_secs(60)), "1m00s");
        assert_eq!(fmt_duration(Duration::from_secs(965)), "16m05s"); // the 348-case sweep
        assert_eq!(fmt_duration(Duration::from_secs(3599)), "59m59s");
        assert_eq!(fmt_duration(Duration::from_secs(3600)), "1h00m");
        assert_eq!(fmt_duration(Duration::from_secs(3304 * 60 + 30)), "55h04m"); // 55h04m CPU
    }

    /// Linear extrapolation from the observed rate. Worked by hand: 25 of 100 cases in 60 s is
    /// 2.4 s/case, so the remaining 75 want 180 s. Cases are not equal-cost (dilute ones finish
    /// faster), which is why the emitted line labels this an estimate.
    #[test]
    fn eta_extrapolates_linearly_from_the_observed_rate() {
        assert_eq!(
            eta(25, 100, Duration::from_secs(60)),
            Some(Duration::from_secs(180))
        );
        // Half done at 10 minutes ⇒ 10 minutes to go.
        assert_eq!(
            eta(174, 348, Duration::from_secs(600)),
            Some(Duration::from_secs(600))
        );
        // Nothing finished yet: no rate to extrapolate from, so no estimate.
        assert_eq!(eta(0, 348, Duration::from_secs(30)), None);
        // Finished: nothing left to wait for.
        assert_eq!(
            eta(348, 348, Duration::from_secs(600)),
            Some(Duration::ZERO)
        );
    }

    /// The throttle is the whole point of the helper: a 348-case sweep must not emit 348 lines.
    /// Nothing before the first interval elapses, then at most one line per interval.
    #[test]
    fn emits_at_most_one_line_per_interval() {
        let p = Progress::new("test", 100);

        // Four cases finish inside the first interval: silence.
        for s in [1, 5, 10, 14] {
            assert_eq!(
                p.complete_at(Duration::from_secs(s)),
                None,
                "emitted before the first interval elapsed (t = {s}s)"
            );
        }

        // The case that crosses the interval boundary emits.
        assert!(
            p.complete_at(EMIT_INTERVAL).is_some(),
            "no line at the interval boundary"
        );

        // ...and the ones right behind it do not.
        for s in [16, 20, 29] {
            assert_eq!(
                p.complete_at(Duration::from_secs(s)),
                None,
                "emitted twice inside one interval (t = {s}s)"
            );
        }

        // The next boundary emits again.
        assert!(
            p.complete_at(2 * EMIT_INTERVAL).is_some(),
            "no line at the second interval boundary"
        );
    }
}
