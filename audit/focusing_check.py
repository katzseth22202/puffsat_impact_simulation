"""Q-D: the contour's survivable ceiling omits the concave focusing factor the discrete
frontier applies. Re-solve the contour with it, as a fixed point (focusing depends on shape,
shape depends on rho, rho depends on focusing)."""
import sys, json; sys.path.insert(0,"python")
from pathlib import Path
from puffsat import contour, heavyplate as hp
from puffsat.analysis import impact_density

rows=[json.loads(l) for f in ("data/results/sweep_geometry_m40.jsonl","data/results/sweep_geometry_wide.jsonl")
      for l in Path(f).read_text().splitlines() if l.strip()]
contour._GEO_CACHE=[{k:float(v) for k,v in r.items() if isinstance(v,(int,float))} for r in rows]
flat={(r["l_over_d"],r["r_foot_over_r"]):r["peak_local_pressure"] for r in rows if r["d_over_d"]==0.0}

def bilin(field,ld,rf,dd):
    """Same bilinear scheme contour.eta_at uses, on an arbitrary field."""
    sub=[r for r in contour._GEO_CACHE if abs(r["d_over_d"]-dd)<1e-9]
    xs=sorted({r["l_over_d"] for r in sub}); ys=sorted({r["r_foot_over_r"] for r in sub})
    node={(r["l_over_d"],r["r_foot_over_r"]):field(r) for r in sub}
    def br(v,t):
        t=min(max(t,v[0]),v[-1])
        for a,b in zip(v,v[1:]):
            if a<=t<=b: return a,b,(0.0 if b==a else (t-a)/(b-a))
        return v[-2],v[-1],1.0
    x0,x1,tx=br(xs,ld); y0,y1,ty=br(ys,rf)
    a=node[(x0,y0)]*(1-ty)+node[(x0,y1)]*ty; b=node[(x1,y0)]*(1-ty)+node[(x1,y1)]*ty
    return a*(1-tx)+b*tx

def focusing_at(ld,rf,dd):
    if dd==0.0: return 1.0
    return bilin(lambda r:r["peak_local_pressure"],ld,rf,dd)/bilin(
        lambda r:flat[(r["l_over_d"],r["r_foot_over_r"])],ld,rf,dd)

sw=hp.read_heavyplate_sweep()
LOD_LO,LOD_HI=contour.L_OVER_D_BOX

def solve(v,mass,radius,rf_hi,dd,with_focus,p_limit=400.0e6,samples=1200):
    c=hp.stagnation_coefficient_at_v(sw,min(v,63000.))
    rho_box=impact_density(LOD_LO,0.3,mass,radius)
    rho=min(contour.rho_ceiling(v,c,p_limit),rho_box)
    for _ in range(60):                      # fixed point on focusing
        best=None
        for i in range(samples):
            rf=0.3+i*(rf_hi-0.3)/(samples-1)
            lod=contour.l_over_d_for(rho,rf,mass,radius)
            if not (LOD_LO<=lod<=LOD_HI): continue
            eta=contour.eta_at(lod,rf,dd)
            if best is None or eta>best[0]: best=(eta,lod,rf)
        if best is None: return None
        if not with_focus: break
        fo=focusing_at(best[1],best[2],dd)
        new=min(contour.rho_ceiling(v,c,p_limit)/fo,rho_box)
        if abs(new-rho)<1e-9: rho=new; break
        rho=0.5*rho+0.5*new                  # damped, focusing is steep in rf
    e=hp.e_eff_interpolator_at_v(sw,min(v,63000.))(rho)
    return dict(rho=rho,eta=best[0],lod=best[1],rf=best[2],f=best[0]*(1+e)/2,
                focus=focusing_at(best[1],best[2],dd))

for tag,mass,radius,rf_hi in (("HEAVY 100kg/15m",100.,15.,0.7),("PAPER 25kg/5m",25.,5.,1.0)):
    print(f"\n{tag}   (box rf/R <= {rf_hi})")
    print(f"{'v':>4} | {'as published: rho':>17} {'f':>6} | {'+focusing, concave: rho':>23} {'foc':>5} {'f':>6}"
          f" | {'+focusing, FLAT: rho':>20} {'f':>6}")
    for v in (28000,34000,45000,55000,63000):
        a=solve(v,mass,radius,rf_hi,0.10,False)
        b=solve(v,mass,radius,rf_hi,0.10,True)
        c=solve(v,mass,radius,rf_hi,0.00,True)
        fa=f"{a['rho']:17.4f} {a['f']:6.3f}" if a else f"{'--':>17} {'--':>6}"
        fb=f"{b['rho']:23.4f} {b['focus']:5.2f} {b['f']:6.3f}" if b else f"{'--':>23} {'--':>5} {'--':>6}"
        fc=f"{c['rho']:20.4f} {c['f']:6.3f}" if c else f"{'--':>20} {'--':>6}"
        print(f"{v/1000:4.0f} | {fa} | {fb} | {fc}")
