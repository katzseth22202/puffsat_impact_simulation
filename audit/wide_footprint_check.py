import sys, math, json; sys.path.insert(0,"python")
from pathlib import Path
from puffsat import contour, heavyplate as hp
from puffsat.analysis import impact_density

# Union the published m40 grid with design 7's remaining r_foot/R nodes.
rows=[json.loads(l) for f in ("data/results/sweep_geometry_m40.jsonl","data/results/sweep_geometry_wide.jsonl")
      for l in Path(f).read_text().splitlines() if l.strip()]
contour._GEO_CACHE=[{k:float(v) for k,v in r.items() if isinstance(v,(int,float))} for r in rows]

print("eta_capture at d/D=0.10, the headline curvature (new nodes marked *):")
print(f"{'rf/R':>6} " + " ".join(f"L/D={x:<5.1f}" for x in (0.3,0.6,1.0)))
for rf in (0.3,0.5,0.7,0.8,0.9,1.0):
    cells=" ".join(f"{contour.eta_at(ld,rf,0.10):9.4f}" for ld in (0.3,0.6,1.0))
    print(f"{rf:6.1f} {cells}{'  *' if rf>0.7 else ''}")

sw=hp.read_heavyplate_sweep()
LOD_LO,LOD_HI=contour.L_OVER_D_BOX

def solve(rho,mass,radius,rf_hi,samples=3000):
    best=None
    for i in range(samples):
        rf=0.3+i*(rf_hi-0.3)/(samples-1)
        lod=contour.l_over_d_for(rho,rf,mass,radius)
        if not (LOD_LO<=lod<=LOD_HI): continue
        eta=contour.eta_at(lod,rf,0.10)
        if best is None or eta>best[0]: best=(eta,lod,rf)
    return best

print("\nPAPER PLATE (25 kg / 5 m): does completing the grid recover the gap?")
hdr=(f"{'v':>5} {'rho':>7} | {'box<=0.7: rf/R':>14} {'L/D':>6} {'eta':>7} {'f':>6}"
     f" | {'box<=1.0: rf/R':>14} {'L/D':>6} {'eta':>7} {'f':>6} | {'gain':>6} {'vs heavy':>8}")
print(hdr); print("-"*len(hdr))
for v in (34000,40000,45000,50000,55000,58000,63000):
    c=hp.stagnation_coefficient_at_v(sw,min(v,63000.))
    e_at=hp.e_eff_interpolator_at_v(sw,min(v,63000.))
    rho=min(contour.rho_ceiling(v,c,400.0e6),impact_density(0.3,0.3,25.,5.))
    e=e_at(rho)
    old=solve(rho,25.,5.,0.7); new=solve(rho,25.,5.,1.0)
    fh=solve(rho,100.,15.,1.0)
    f_old=old[0]*(1+e)/2 if old else None; f_new=new[0]*(1+e)/2 if new else None
    f_h=fh[0]*(1+e)/2 if fh else None
    o=(f"{old[2]:14.3f} {old[1]:6.3f} {old[0]:7.4f} {f_old:6.3f}" if old
       else f"{'--':>14} {'--':>6} {'--':>7} {'INFEA':>6}")
    n=(f"{new[2]:14.3f} {new[1]:6.3f} {new[0]:7.4f} {f_new:6.3f}" if new
       else f"{'--':>14} {'--':>6} {'--':>7} {'INFEA':>6}")
    g=f"{f_new-f_old:+6.3f}" if (old and new) else f"{'--':>6}"
    h=f"{f_new-f_h:+8.3f}" if (new and f_h) else f"{'--':>8}"
    print(f"{v/1000:5.0f} {rho:7.4f} | {o} | {n} | {g} {h}")

lo=impact_density(LOD_HI,1.0,25.,5.)
k=400.0e6/hp.stagnation_coefficient_at_v(sw,55000.)
print(f"\nreach: dilutest cloud now {lo:.4f} kg/m^3 -> contour flyable to {math.sqrt(k/lo)/1000:.1f} km/s")
