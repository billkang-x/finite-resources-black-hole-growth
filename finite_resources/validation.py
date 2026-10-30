"""Independent analytic and adaptive-integration controls for the core kernels."""
from pathlib import Path
import argparse,json
import numpy as np
from scipy.integrate import solve_ivp,quad
from numpy.polynomial.legendre import leggauss
from .inner_transport import (scatter,gamma_of,moller,static_to_pg,pg_to_static,
                             specific_energy,orbit_step,rhs,incoming)
from .orbit_control import bound_critical
from .profiles import G,C,M0
from .bounds import outward_kernel,column_cost


def collision_control():
    rng=np.random.default_rng(261030);ee=0.;pe=0.
    for _ in range(5000):
        p1=rng.normal(size=3)*10**rng.uniform(-3,1.5)
        p2=rng.normal(size=3)*10**rng.uniform(-3,1.5)
        d=rng.normal(size=3);d/=np.linalg.norm(d)
        q1,q2=scatter(p1,p2,d)
        before=gamma_of(p1)+gamma_of(p2)
        ee=max(ee,abs(gamma_of(q1)+gamma_of(q2)-before)/before)
        pe=max(pe,np.linalg.norm(q1+q2-p1-p2)/before)
    v=.4;me=0.
    for u,w,expected in [(np.array([v,0,0]),np.array([-v,0,0]),2*v),
                           (np.array([v,0,0]),np.array([v,0,0]),0.),
                           (np.array([v,0,0]),np.array([0,v,0]),np.sqrt(2*v*v-v**4))]:
        me=max(me,abs(moller(u/np.sqrt(1-u@u),w/np.sqrt(1-w@w))-expected))
    return dict(pairs=5000,energy_error=float(ee),momentum_error=float(pe),moller_error=float(me),passed=bool(max(ee,pe,me)<1e-12))


def coordinate_control():
    rng=np.random.default_rng(261031);worst=0.;energy=0.
    for _ in range(5000):
        x=np.exp(rng.uniform(np.log(2.01),np.log(1000)));v=rng.uniform(0,.99);mu=rng.uniform(-1,1)
        vr=v*mu;vt=v*np.sqrt(1-mu*mu)
        p,ell=static_to_pg(x,vr,vt);rr,tt=pg_to_static(x,p,ell)
        worst=max(worst,abs(vr-rr),abs(vt-tt))
        ref=np.sqrt(1-2/x)/np.sqrt(1-v*v)
        energy=max(energy,abs(specific_energy(x,p,ell)/ref-1))
    return dict(samples=5000,velocity_error=float(worst),killing_energy_error=float(energy),passed=bool(max(worst,energy)<1e-11))


def orbit_control():
    rows=[]
    for v,mu in [(.05,.3),(.1,.3),(.1,.7),(.14,.7),(.14,.95),(.1,.995)]:
        p0,ell=static_to_pg(100.,-v*mu,v*np.sqrt(1-mu*mu))
        e=specific_energy(100.,p0,ell);_,critical=bound_critical(e)
        def horizon(t,y):return y[0]-2
        def outward(t,y):return y[0]-100
        horizon.terminal=True;horizon.direction=-1
        outward.terminal=True;outward.direction=1
        def f(t,y):return rhs(y[0],y[1],ell)
        ref=solve_ivp(f,[0,20000],[100,p0],method='DOP853',rtol=1e-12,atol=1e-13,events=[horizon,outward],max_step=5)
        captured=bool(len(ref.t_events[0]));expect=bool(ell<critical)
        if not (len(ref.t_events[0]) or len(ref.t_events[1])):raise ArithmeticError('Unresolved orbit')
        end=float(ref.t[-1]*.8)
        fine=solve_ivp(f,[0,end],[100,p0],method='DOP853',rtol=2e-13,atol=1e-13).y[:,-1]
        errors=[]
        for dt in [2.,1.,.5]:
            x=100.;p=p0;t=0.
            while t<end:
                h=min(dt,end-t);x,p=orbit_step(x,p,ell,h);t+=h
            errors.append(float(np.linalg.norm((np.array([x,p])-fine)/np.maximum(np.abs(fine),1))))
        rows.append(dict(speed=v,mu=mu,analytic_capture=expect,adaptive_capture=captured,
            fixed_step_errors=errors,passed=captured==expect and errors[-1]<2e-5))
    return dict(rows=rows,passed=all(r['passed'] for r in rows))


def flux_capture_reference(case):
    lo,hi,cdf,flux=incoming(case);weights=np.diff(np.r_[0,cdf])
    z,w=leggauss(48);z=(z+1)/2;w=w/2;answer=0.
    for a,b,prob in zip(lo,hi,weights):
        if prob==0:continue
        v2=np.sqrt(a*a+z*(b*b-a*a));gam=1/np.sqrt(1-v2)
        _,ell=bound_critical(np.sqrt(.98)*gam)
        fraction=np.minimum(ell**2/(10000*gam**2*v2),1)
        answer+=prob*np.dot(w,fraction)
    return dict(case=case,inward_flux_msun_myr=flux,collisionless_capture_fraction=float(answer))


def residence_control():
    # The same finite radial layer must occur in both the column and kernel.
    rate=7.;speed=3.;inner=2.;ratios=[]
    for factor in [1.01,1.5,2,10,1000]:
        outer=inner*factor
        col=quad(lambda r:rate/(4*np.pi*speed*r*r),inner,outer,epsabs=1e-12)[0]
        kernel=outward_kernel(1.,speed,inner,outer)
        ratios.append(4*np.pi*col/kernel/rate)
    error=float(np.max(np.abs(np.array(ratios)-1)))
    # Correct pc/Myr constants: t_g for a 10^5 Msun seed is 0.4926 seconds.
    seconds=G*M0/C**3*1e6*31557600
    return dict(finite_layer_relative_error=error,tg_seconds=seconds,
                passed=error<1e-11 and .492<seconds<.494)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',default='outputs/validation',type=Path)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    checks=dict(collisions=collision_control(),coordinates=coordinate_control(),
                geodesics=orbit_control(),residence_and_units=residence_control())
    out=dict(checks=checks,passed=all(x['passed'] for x in checks.values()),
             flux_controls=[flux_capture_reference(c) for c in ['compact','extended']])
    (a.output_dir/'validation.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
    if not out['passed']:raise SystemExit(1)

if __name__=='__main__':main()
