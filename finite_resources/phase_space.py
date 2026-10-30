"""Nonnegative f(E) audit of the frozen finite-halo profiles.

An isotropic stationary distribution in a spherical reflecting box may have
positive energy relative to the wall. Its density and radial pressure are
linear energy-bin integrals. Simultaneously fitting both tests a preparation
stronger than assigning independent local speed distributions.
"""
from pathlib import Path
import sys,json,time,hashlib
import numpy as np
from scipy.optimize import lsq_linear

from .profiles import Core,G,C,KMS,M0

def powers_difference(a,b,p):
    """a**p-b**p for a>=b>=0 without subtracting near-equal powers."""
    a=np.maximum(a,0);b=np.maximum(b,0)
    result=np.zeros(np.broadcast_shapes(a.shape,b.shape))
    a,b=np.broadcast_arrays(a,b)
    take=a>0
    result[take]=a[take]**p*(-np.expm1(p*np.log(np.maximum(b[take]/a[take],1e-300))))
    return result


def moments_matrix(phi,edges):
    hi=np.maximum(edges[None,1:]-phi[:,None],0)
    lo=np.maximum(edges[None,:-1]-phi[:,None],0)
    density=8*np.sqrt(2)*np.pi/3*powers_difference(hi,lo,1.5)
    pressure=16*np.sqrt(2)*np.pi/15*powers_difference(hi,lo,2.5)
    return density,pressure


def prepare(case,ne=640,nr=200,output_dir="outputs/preparation"):
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    core=Core(case)
    rr,rho,sig,b,md,vmax=core.state(M0,'fixed')
    acc=G*(M0+md+core.baryons(rr)[0])/rr**2
    from scipy.integrate import cumulative_trapezoid
    phi=cumulative_trapezoid(acc[::-1],rr[::-1],initial=0)[::-1]
    if not np.all(np.isfinite(phi)) or not (phi[0]<0 and phi[-1]==0 and np.all(np.diff(phi)>0)):
        raise ArithmeticError('Potential must increase from negative values to zero at the wall')
    # Phi(wall)=0; central negative energies span approximately 11 decades.
    select=np.geomspace(100*G*M0/C**2,core.rvir,nr)
    interp=lambda a:np.interp(np.log(select),np.log(rr),a)
    psel=interp(phi);rhosel=np.exp(interp(np.log(rho)));psel_pressure=np.exp(interp(np.log(rho*sig**2)))
    vwall2=float(sig[-1]**2)
    emin=float(psel[0])*1.04
    negative=-np.geomspace(-emin,vwall2*1e-4,ne*3//4)
    positive=np.r_[0,np.geomspace(vwall2*1e-4,32*vwall2,ne-len(negative))]
    edges=np.unique(np.r_[negative,positive])
    if not np.isfinite(edges).all() or len(edges)<ne or np.any(np.diff(edges)<=0):
        raise ArithmeticError('Invalid energy grid')
    dr,pr=moments_matrix(psel,edges)
    a=np.vstack([dr/rhosel[:,None],pr/psel_pressure[:,None]])
    scale=np.maximum(np.linalg.norm(a,axis=0),1e-100)
    scaled=a/scale
    # No smoothing penalty: an infeasible equilibrium must remain visible.
    solution=lsq_linear(scaled,np.ones(2*nr),bounds=(0,np.inf),method='bvls',tol=2e-11,max_iter=5000)
    f=np.maximum(solution.x,0)/scale
    rhorec=dr@f;prec=pr@f
    er=rhorec/rhosel-1;ep=prec/psel_pressure-1
    savedmin=float(core.input['r_pc'][0])
    region=select>=savedmin
    fit=dict(density_rms=float(np.sqrt(np.mean(er**2))),pressure_rms=float(np.sqrt(np.mean(ep**2))),
             density_max=float(np.max(np.abs(er))),pressure_max=float(np.max(np.abs(ep))),
             saved_domain_density_rms=float(np.sqrt(np.mean(er[region]**2))),
             saved_domain_pressure_rms=float(np.sqrt(np.mean(ep[region]**2))),
             saved_domain_density_max=float(np.max(np.abs(er[region]))),
             saved_domain_pressure_max=float(np.max(np.abs(ep[region]))))
    plan={'density_relative_rms_tolerance':.03,'pressure_relative_rms_tolerance':.05,'density_local_relative_tolerance':.10,'pressure_local_relative_tolerance':.15}
    passed=fit['density_rms']<plan['density_relative_rms_tolerance'] and fit['pressure_rms']<plan['pressure_relative_rms_tolerance'] and fit['density_max']<plan['density_local_relative_tolerance'] and fit['pressure_max']<plan['pressure_local_relative_tolerance']
    label=f'{case}_E{ne}_R{nr}'
    np.savez_compressed(output_dir/(label+'.npz'),r=select,phi=psel,rho_input=rhosel,pressure_input=psel_pressure,rho_reconstructed=rhorec,pressure_reconstructed=prec,energy_edges=edges,f=f,rho_error=er,pressure_error=ep)
    np.savetxt(output_dir/(case+'_distribution.csv'),np.column_stack([edges[:-1],edges[1:],f]),delimiter=',',header='energy_low,energy_high,f_mass',comments='',fmt='%.17e')
    return dict(case=case,energy_bins=len(edges)-1,radial_samples=nr,optimizer_success=bool(solution.success),
                optimality=float(solution.optimality),iterations=int(solution.nit),minimum_f=float(f.min()),
                active_energy_bins=int(np.count_nonzero(f>0)),passed_initial_state=bool(passed),**fit,
                worst_density_radius_pc=float(select[np.argmax(np.abs(er))]),
                worst_pressure_radius_pc=float(select[np.argmax(np.abs(ep))]),raw_file=label+'.npz')


def main():
    import argparse
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',choices=['compact','extended'],default='compact')
    ap.add_argument('--energy-bins',type=int,default=640)
    ap.add_argument('--radial-samples',type=int,default=200)
    ap.add_argument('--output-dir',default='outputs/preparation')
    a=ap.parse_args();row=prepare(a.case,a.energy_bins,a.radial_samples,a.output_dir)
    (Path(a.output_dir)/(a.case+'_fit.json')).write_text(json.dumps(row,indent=2))
    print(json.dumps(row,indent=2))
    if not row['passed_initial_state']:raise SystemExit(1)

if __name__=='__main__':main()
