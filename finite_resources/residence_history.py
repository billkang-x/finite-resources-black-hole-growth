"""Integrate all visits of distinct parcels without double-counting supplied mass.

JSON input: {"parcels": [{"id": ..., "mass_g": ..., "source": "initial"
or "external", "segments": [{"time_s": [...], "radius_cm": [...],
"opacity_cm2_g": scalar_or_array}, ...]}], "reference":
{"speed_cm_s": ..., "radius_cm": ..., "opacity_cm2_g": ...}}.

Segments are disjoint time intervals when the parcel contributes to the column.
Repeated segments of one parcel increase exposure but never its distinct mass.
The kernel alone does not establish dynamical support or overlap with luminosity.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def _integral(values, time):
    return float(np.trapezoid(values,time) if hasattr(np,'trapezoid') else np.trapz(values,time))


def integrate_history(history):
    ref=history['reference']
    v,r,k=(float(ref[key]) for key in ('speed_cm_s','radius_cm','opacity_cm2_g'))
    if not all(math.isfinite(x) and x>0 for x in (v,r,k)):
        raise ValueError('positive finite reference speed, radius and opacity required')
    seen=set(); rows=[]; total_q=0.; initial=0.; external=0.
    for parcel in history['parcels']:
        key=str(parcel['id'])
        if key in seen:
            raise ValueError('duplicate parcel id: '+key)
        seen.add(key)
        mass=float(parcel['mass_g'])
        source=parcel['source']
        if not math.isfinite(mass) or mass<=0 or source not in ('initial','external'):
            raise ValueError('each distinct parcel needs positive mass and a declared source')
        previous_end=-math.inf; q=0.
        for seg in parcel['segments']:
            t=np.asarray(seg['time_s'],dtype=float)
            rad=np.asarray(seg['radius_cm'],dtype=float)
            opacity=np.broadcast_to(np.asarray(seg.get('opacity_cm2_g',k),dtype=float),t.shape)
            if t.ndim!=1 or len(t)<2 or rad.shape!=t.shape or not np.all(np.isfinite(t)):
                raise ValueError('segment time/radius arrays must have identical one-dimensional shapes')
            if np.any(np.diff(t)<=0) or t[0]<0 or t[0]<previous_end:
                raise ValueError('segments must be ordered, nonnegative, and nonoverlapping in time')
            if not np.all(np.isfinite(rad)) or np.any(rad<=0) or not np.all(np.isfinite(opacity)) or np.any(opacity<0):
                raise ValueError('radii must be positive; opacities nonnegative and finite')
            previous_end=float(t[-1])
            q+=_integral(opacity/rad**2,t)
        total_q+=mass*q
        if source=='initial': initial+=mass
        else: external+=mass
        rows.append(dict(id=key,source=source,mass_g=mass,segments=len(parcel['segments']),
                         q_s_per_g=q,mass_weighted_q_s=mass*q))
    distinct=initial+external
    if distinct<=0:
        raise ValueError('at least one parcel is required')
    budget=float(history.get('available_budget_g',distinct))
    if not math.isfinite(budget) or budget<distinct:
        raise ValueError('available budget cannot be smaller than represented distinct mass')
    ratio=total_q*v*r/(k*budget)
    return dict(distinct_initial_g=initial,distinct_external_g=external,
                distinct_represented_g=distinct,available_budget_g=budget,
                Q_s=total_q,R_Q=ratio,parcels=rows,
                scope='Quadrature of prescribed labeled histories; no dynamics or g(t) overlap inferred')


def demo_histories():
    """Kinematic controls in cgs; not fitted or production RHD trajectories."""
    r0=3.0856775814913673e18; v=1e8; mass=1e9*1.98847e33
    ref=dict(radius_cm=r0,speed_cm_s=v,opacity_cm2_g=1.)
    crossing=r0/v
    def outward(duration=1., start=0.):
        t=np.linspace(0,duration*crossing,2049)
        return dict(time_s=(start*crossing+t).tolist(),radius_cm=(r0+v*t).tolist())
    def history(segments,source='initial',fraction=1.):
        return dict(reference=ref,parcels=[dict(id='a',mass_g=mass*fraction,source=source,segments=segments)])
    first=outward()
    t=np.linspace(crossing,2*crossing,2049)
    inward=dict(time_s=t.tolist(),radius_cm=(2*r0-v*(t-crossing)).tolist())
    cases=dict(short_residence=history([outward(.1)]),one_outward_pass=history([first]),
               returned_same_parcel=history([first,inward,outward(start=2.)]))
    # A recurrent luminous mask does not restart the outward trajectory.
    long=outward(3.)
    split=[]
    for begin,end in [(0,683),(683,1366),(1366,2048)]:
        split.append({key:value[begin:end+1] for key,value in long.items()})
    cases['outward_episode_subdivision']=history(split)
    cases['continuous_same_outward_path']=history([long])
    cases['external_fresh_parcel']=history([first])
    cases['external_fresh_parcel']['parcels'].append(dict(id='b',source='external',mass_g=mass,segments=[outward(start=2.)]))
    return cases


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',type=Path,help='Labeled-parcel JSON; omit to run kinematic controls')
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/history'))
    args=ap.parse_args()
    histories={'user':json.loads(args.input.read_text(encoding='utf-8'))} if args.input else demo_histories()
    report={name:integrate_history(value) for name,value in histories.items()}
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'residence_histories.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps({name:{k:row[k] for k in ('R_Q','distinct_initial_g','distinct_external_g')} for name,row in report.items()},indent=2))


if __name__=='__main__':
    main()
