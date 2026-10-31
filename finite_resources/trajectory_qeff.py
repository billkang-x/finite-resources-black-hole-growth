"""Matched-cell Sigma and q_eff extraction from an axisymmetric wind snapshot.

Input NPZ keys: state[field,radius,ray], r_edges[cm], mu_edges=cos(theta).
state[0] is density [g/cm^3]; state[1] is radial momentum density [g/cm^2/s].
q_eff = integral dr/(v_r r^2) on outward cells above a declared threshold.
This is a snapshot streamline diagnostic, not a time-resolved labeled history
or a converged production RHD result. No simulation snapshot is distributed.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from .ceers_selection import sha256
from .joint_constraints import MSUN, YEAR


def extract_arrays(state, r_edges, mu_edges, threshold_km_s=1.):
    state=np.asarray(state,dtype=float); edges=np.asarray(r_edges,dtype=float)
    mu=np.asarray(mu_edges,dtype=float)
    if state.ndim!=3 or state.shape[0]<2 or edges.ndim!=1 or mu.ndim!=1:
        raise ValueError('expected state[field,radius,ray] and one-dimensional edges')
    if state.shape[1:]!=(len(edges)-1,len(mu)-1):
        raise ValueError('edge arrays do not match state dimensions')
    if len(edges)<2 or len(mu)<2 or not np.all(np.isfinite(edges)) or np.any(edges<=0) or np.any(np.diff(edges)<=0):
        raise ValueError('positive increasing radial edges required')
    if not np.all(np.isfinite(mu)) or np.any(np.diff(mu)<=0) or mu[0]<-1 or mu[-1]>1:
        raise ValueError('increasing mu edges within [-1,1] required')
    if not math.isfinite(threshold_km_s) or threshold_km_s<=0:
        raise ValueError('outward velocity threshold must be positive')
    rho,momentum=state[:2]
    if np.any(~np.isfinite(rho)) or np.any(~np.isfinite(momentum)) or np.any(rho<0) or np.any((rho==0)&(momentum!=0)):
        raise ValueError('finite nonnegative density and consistent momentum required')
    vr=np.divide(momentum,rho,out=np.zeros_like(rho),where=rho>0)
    mask=vr>threshold_km_s*1e5
    inverse_v=np.divide(1.,vr,out=np.zeros_like(vr),where=mask)
    q_ray=np.sum((1/edges[:-1]-1/edges[1:])[:,None]*inverse_v,axis=0)
    sigma_ray=np.sum(np.where(mask,rho,0.)*np.diff(edges)[:,None],axis=0)
    weights=np.diff(mu)/2
    valid=q_ray>0
    if not np.any(valid):
        raise ValueError('no outward cells above threshold')
    q=float(np.dot(weights,q_ray));sigma=float(np.dot(weights,sigma_ray))
    return dict(radial_cells=state.shape[1],angular_rays=state.shape[2],
                velocity_threshold_km_s=threshold_km_s,r_inner_cm=float(edges[0]),r_outer_cm=float(edges[-1]),
                outward_solid_angle_fraction=float(np.sum(weights[valid])),
                q_eff_weighted_s_cm2=q,sigma_weighted_g_cm2=sigma,
                a_sigma_msun_yr=4*math.pi*sigma/q*YEAR/MSUN,
                q_eff_ray_s_cm2=q_ray.tolist(),sigma_ray_g_cm2=sigma_ray.tolist(),
                angular_weights=weights.tolist(),
                definition='Sigma and q_eff use exactly the same outward cells and finite radial limits',
                scope='Eulerian snapshot streamline reconstruction; not a time-resolved mass-weighted R_Q')


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('snapshot',type=Path)
    ap.add_argument('--threshold-km-s',type=float,default=1.)
    ap.add_argument('--exposure-myr',type=float,default=50.)
    ap.add_argument('--budget-msun',type=float,default=1e9)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/qeff'))
    args=ap.parse_args()
    if not math.isfinite(args.exposure_myr) or args.exposure_myr<0 or not math.isfinite(args.budget_msun) or args.budget_msun<=0:
        ap.error('nonnegative finite exposure and positive finite budget required')
    with np.load(args.snapshot,allow_pickle=False) as data:
        report=extract_arrays(data['state'],data['r_edges'],data['mu_edges'],args.threshold_km_s)
    report.update(snapshot_filename=args.snapshot.name,snapshot_sha256=sha256(args.snapshot),
                  exposure_myr=args.exposure_myr,budget_msun=args.budget_msun,
                  Lambda=report['a_sigma_msun_yr']*args.exposure_myr*1e6/args.budget_msun)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'trajectory_qeff.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('q_eff_weighted_s_cm2','sigma_weighted_g_cm2','a_sigma_msun_yr','Lambda','scope')},indent=2))


if __name__=='__main__':
    main()
