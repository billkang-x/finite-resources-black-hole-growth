"""Transient inner Boltzmann DSMC in fixed Schwarzschild geometry.

Units: r_g=GM/c^2, t_g=GM/c^3, particle specific momentum p/(mc).
Painleve-Gullstrand slices cross the horizon without a coordinate singularity.
This is an INNER RESPONSE experiment with prescribed incoming f(E), not a
finite-halo 100 Myr growth calculation. No immediate reflection is used.
"""
from pathlib import Path
import sys, json, time, argparse
import numpy as np
from numba import njit
from scipy.integrate import cumulative_trapezoid

from .profiles import Core,G,C,M0,load_distribution
S=1.988409870698051e33/3.085677581491367e18**2


@njit(cache=True)
def gamma_of(p):
    return np.sqrt(1.0+np.dot(p,p))


@njit(cache=True)
def boost(q0,p,beta):
    b2=np.dot(beta,beta)
    if b2==0.0:return q0,p.copy()
    g=1.0/np.sqrt(1.0-b2)
    d=np.dot(beta,p)
    return g*(q0-d),p+(((g-1.0)*d/b2)-g*q0)*beta


@njit(cache=True)
def scatter(p1,p2,direction):
    """Constant-sigma isotropic equal-mass elastic COM map, c=m=1."""
    g1=gamma_of(p1);g2=gamma_of(p2)
    beta=(p1+p2)/(g1+g2)
    e1,_=boost(g1,p1,beta);e2,_=boost(g2,p2,beta)
    star=.5*(e1+e2)
    pstar=np.sqrt(max(0.0,(star-1.0)*(star+1.0)))
    _,out1=boost(star,pstar*direction,-beta)
    _,out2=boost(star,-pstar*direction,-beta)
    return out1,out2


@njit(cache=True)
def moller(p1,p2):
    v1=p1/gamma_of(p1);v2=p2/gamma_of(p2)
    diff=v1-v2;cross=np.cross(v1,v2)
    return np.sqrt(max(0.0,np.dot(diff,diff)-np.dot(cross,cross)))


@njit(cache=True)
def specific_energy(x,p,ell):
    return np.sqrt(1+p*p+ell*ell/(x*x))-np.sqrt(2/x)*p


@njit(cache=True)
def static_to_pg(x,vr,vt):
    gs=1/np.sqrt(1-vr*vr-vt*vt)
    beta=np.sqrt(2/x);alpha=np.sqrt(1-2/x)
    pr=(gs*vr+beta*gs)/alpha
    return pr,x*gs*vt


@njit(cache=True)
def pg_to_static(x,p,ell):
    ge=np.sqrt(1+p*p+ell*ell/(x*x))
    beta=np.sqrt(2/x);alpha=np.sqrt(1-2/x)
    gs=(ge-beta*p)/alpha
    prs=(p-beta*ge)/alpha
    return prs/gs,ell/x/gs


@njit(cache=True)
def rhs(x,p,ell):
    ge=np.sqrt(1+p*p+ell*ell/(x*x));beta=np.sqrt(2/x)
    return p/ge-beta,ell*ell/(ge*x*x*x)-beta*p/(2*x)


@njit(cache=True)
def orbit_step(x,p,ell,dt):
    k1,l1=rhs(x,p,ell)
    k2,l2=rhs(x+.5*dt*k1,p+.5*dt*l1,ell)
    k3,l3=rhs(x+.5*dt*k2,p+.5*dt*l2,ell)
    k4,l4=rhs(x+dt*k3,p+dt*l3,ell)
    return x+dt*(k1+2*k2+2*k3+k4)/6,p+dt*(l1+2*l2+2*l3+l4)/6


def incoming(case,xmatch=100.,distribution=None):
    core=Core(case)
    rr,_,_,_,md,_=core.state(M0,'fixed')
    acc=G*(M0+md+core.baryons(rr)[0])/rr**2
    phi=cumulative_trapezoid(acc[::-1],rr[::-1],initial=0)[::-1]
    rm=xmatch*G*M0/C**2
    ph=float(np.interp(np.log(rm),np.log(rr),phi))
    data=load_distribution(case,distribution)
    edges=data['energy_edges'];f=data['f']
    lo=np.maximum(2*(edges[:-1]-ph),0.);hi=np.maximum(2*(edges[1:]-ph),0.)
    if np.any((hi>=C*C)&(f>0)):
        raise ValueError('Input distribution contains superluminal matching speeds')
    rates=f*(hi*hi-lo*lo)/4
    if rates.sum()<=0:raise ValueError('Input has no inward flux at the matching surface')
    flux=4*np.pi*np.pi*rm*rm*np.sum(rates)
    # Static proper-time flux -> fixed-radius PG/Schwarzschild coordinate time.
    flux*=np.sqrt(1-2/xmatch)
    cdf=np.cumsum(rates);cdf/=cdf[-1]
    return lo/(C*C),hi/(C*C),cdf,float(flux)


@njit(cache=True)
def simulate(lo,hi,cdf,xmatch,inner,cells,dt,total_time,nin,sigma_weight,seed,capacity):
    np.random.seed(seed)
    x=np.zeros(capacity);p=np.zeros(capacity);ell=np.zeros(capacity)
    born=np.zeros(capacity);initial_energy=np.zeros(capacity)
    identity=np.zeros(capacity,np.int64)
    # One immutable label per equal-rest-mass packet. Final state code:
    # 0 = still resident, 1 = horizon captured, 2 = returned through r_match.
    records=np.zeros((nin,11))
    n=0;now=0.;injected=0;captured=0;returned=0;events=0;candidate=0
    ein=0.;ecap=0.;eret=0.;event_drift=0.;step_drift=0.
    cap_late=0;ret_late=0;ein_late=0.;ecap_late=0.;eret_late=0.
    injected_late=0;res_time_sum=0.;res_time_sq=0.;overflow=False
    re=np.exp(np.linspace(np.log(inner),np.log(xmatch),cells+1))
    counts=np.zeros(cells,np.int64)
    indices=np.zeros((cells,capacity),np.int64)
    sample_count=0;population_sum=0.;population_sq=0.
    history=np.zeros((12,5))
    speedmax=0.;event_energy_error=0.;event_momentum_error=0.
    while now<total_time:
        h=min(dt,total_time-now)
        block=min(11,int(now/total_time*12))
        # Deterministic equally spaced arrival times make equal macroparticle
        # weights compatible with the declared coordinate-time mass flux.
        should=min(nin,int((now+h)/total_time*nin))
        while injected<should:
            if n==capacity:
                overflow=True;break
            u=np.random.random();j=np.searchsorted(cdf,u)
            v4=lo[j]*lo[j]+np.random.random()*(hi[j]*hi[j]-lo[j]*lo[j])
            v=np.sqrt(np.sqrt(v4));mu=np.sqrt(np.random.random())
            # The inward flux distribution has density 2 mu, not uniform mu.
            pp,ll=static_to_pg(xmatch,-v*mu,v*np.sqrt(1-mu*mu))
            x[n]=xmatch;p[n]=pp;ell[n]=ll;born[n]=now
            en=specific_energy(xmatch,pp,ll);initial_energy[n]=en
            identity[n]=injected
            records[injected,0]=now;records[injected,1]=en
            records[injected,2]=ll
            ein+=en;injected+=1
            history[block,0]+=1
            if now>=total_time/2:injected_late+=1;ein_late+=en
            n+=1
        if overflow:break
        keep=0
        for i in range(n):
            en0=specific_energy(x[i],p[i],ell[i])
            xn,pn=orbit_step(x[i],p[i],ell[i],h)
            en1=specific_energy(xn,pn,ell[i])
            step_drift+=en1-en0
            if xn<=inner:
                ident=identity[i]
                records[ident,3]=1;records[ident,4]=now+h
                records[ident,5]=en1;records[ident,6]=ell[i]
                records[ident,7]=xn;records[ident,8]=pn
                captured+=1;ecap+=en1
                history[block,1]+=1
                if now>=total_time/2:cap_late+=1;ecap_late+=en1
                residence=now+h-born[i]
                res_time_sum+=residence;res_time_sq+=residence*residence
                continue
            if xn>=xmatch and rhs(xn,pn,ell[i])[0]>0:
                ident=identity[i]
                records[ident,3]=2;records[ident,4]=now+h
                records[ident,5]=en1;records[ident,6]=ell[i]
                records[ident,7]=xn;records[ident,8]=pn
                returned+=1;eret+=en1
                history[block,2]+=1
                if now>=total_time/2:ret_late+=1;eret_late+=en1
                residence=now+h-born[i]
                res_time_sum+=residence;res_time_sq+=residence*residence
                continue
            x[keep]=xn;p[keep]=pn;ell[keep]=ell[i]
            born[keep]=born[i];initial_energy[keep]=initial_energy[i]
            identity[keep]=identity[i]
            keep+=1
        n=keep
        if sigma_weight>0 and n>1:
            counts[:]=0
            for i in range(n):
                k=np.searchsorted(re,x[i])-1
                if k<0 or k>=cells:raise ValueError('particle outside collision domain')
                indices[k,counts[k]]=i;counts[k]+=1
            for k in range(cells):
                count=counts[k]
                if count<2:continue
                volume=4*np.pi*(re[k+1]**3-re[k]**3)/3
                # Moller flux <=2 in c=1, constant cross section. PG slices
                # have Euclidean proper volume and unit lapse.
                majorant=count*(count-1)*sigma_weight/volume
                attempts=np.random.poisson(majorant*h)
                for _ in range(attempts):
                    candidate+=1
                    a=np.random.randint(count);b=np.random.randint(count-1)
                    if b>=a:b+=1
                    i=indices[k,a];j=indices[k,b]
                    az=2*np.pi*np.random.random()
                    pi=np.array([p[i],ell[i]/x[i],0.])
                    pj=np.array([p[j],ell[j]/x[j]*np.cos(az),ell[j]/x[j]*np.sin(az)])
                    vm=moller(pi,pj)
                    if np.random.random()>=vm/2:continue
                    z=2*np.random.random()-1;az=2*np.pi*np.random.random()
                    direction=np.array([z,np.sqrt(1-z*z)*np.cos(az),np.sqrt(1-z*z)*np.sin(az)])
                    qi,qj=scatter(pi,pj,direction)
                    before=gamma_of(pi)+gamma_of(pj);after=gamma_of(qi)+gamma_of(qj)
                    event_energy_error=max(event_energy_error,abs(after-before)/before)
                    event_momentum_error=max(event_momentum_error,np.sqrt(np.sum((qi+qj-pi-pj)**2))/before)
                    old=specific_energy(x[i],p[i],ell[i])+specific_energy(x[j],p[j],ell[j])
                    p[i]=qi[0];ell[i]=x[i]*np.hypot(qi[1],qi[2])
                    p[j]=qj[0];ell[j]=x[j]*np.hypot(qj[1],qj[2])
                    new=specific_energy(x[i],p[i],ell[i])+specific_energy(x[j],p[j],ell[j])
                    event_drift+=new-old;events+=1
                    records[identity[i],9]+=1;records[identity[j],9]+=1
                    history[block,3]+=1
        now+=h
        history[block,4]=n
        if now>=total_time/2:
            sample_count+=1;population_sum+=n;population_sq+=n*n
    eres=0.
    for i in range(n):
        en=specific_energy(x[i],p[i],ell[i]);eres+=en
        ident=identity[i]
        records[ident,4]=now;records[ident,5]=en
        records[ident,6]=ell[i];records[ident,7]=x[i]
        records[ident,8]=p[i]
    records[:injected,10]=1.0
    return (injected,captured,returned,n,events,candidate,ein,ecap,eret,eres,
            event_drift,step_drift,cap_late,ret_late,injected_late,ein_late,
            ecap_late,eret_late,res_time_sum,res_time_sq,overflow,now,
            population_sum/max(sample_count,1),population_sq/max(sample_count,1),
            event_energy_error,event_momentum_error,history,records)


def run(case='compact',cells=64,dt=.5,total_time=12000.,nin=6000,
        sigma=1.,seed=20261030,inner=1.99,export_packets=False,output_dir="outputs/inner",distribution=None):
    start=time.perf_counter()
    if min(cells,dt,total_time,nin)<=0 or sigma<0:raise ValueError("Positive resolution and nonnegative cross section required")
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    xmatch=100.;lo,hi,cdf,flux=incoming(case,xmatch,distribution)
    rg=G*M0/C**2;tg=rg/C
    weight=flux*total_time*tg/nin
    sigma_weight=sigma*S*weight/(rg*rg)
    out=simulate(lo,hi,cdf,xmatch,inner,cells,dt,total_time,nin,sigma_weight,seed,max(nin,2000))
    (ni,nc,nr,nres,ev,ca,ei,ec,er,es,de,ds,ncl,nrl,nil,eil,ecl,erl,
     rt,rt2,over,actual,meanpop,meanpop2,ee,pe,history,records)=out
    row=dict(case=case,cells=cells,dt_rg_over_c=dt,requested_duration_rg_over_c=total_time,
        simulated_duration_rg_over_c=actual,simulated_duration_myr=actual*tg,
        cross_section_cm2_g=sigma,seed=seed,inner_absorber_rg=inner,matching_rg=xmatch,
        injected_packets=ni,captured_packets=nc,returned_packets=nr,resident_packets=nres,
        accepted_collisions=ev,candidate_collisions=ca,packet_mass_msun=weight,
        input_mass_rate_msun_myr=flux,
        captured_rest_msun=nc*weight,returned_rest_msun=nr*weight,resident_rest_msun=nres*weight,
        input_killing_energy_msun=ei*weight,captured_killing_energy_msun=ec*weight,
        returned_killing_energy_msun=er*weight,resident_killing_energy_msun=es*weight,
        rest_packet_count_residual=ni-nc-nr-nres,
        energy_ledger_relative_residual=(ec+er+es-ei)/max(ei,1e-100),
        finite_cell_collision_energy_drift_relative=de/max(ei,1e-100),
        orbit_energy_drift_relative=ds/max(ei,1e-100),
        energy_ledger_unexplained_residual=(ec+er+es-ei-de-ds)/max(ei,1e-100),
        local_event_energy_relative_error=ee,local_event_momentum_relative_error=pe,
        history_blocks=[dict(start_rg_over_c=i*total_time/12,end_rg_over_c=(i+1)*total_time/12,
            injected=int(h[0]),captured=int(h[1]),returned=int(h[2]),collisions=int(h[3]),
            resident_at_end=int(h[4])) for i,h in enumerate(history)],
        late_captured_packets=ncl,late_returned_packets=nrl,late_injected_packets=nil,
        late_capture_per_injected=ncl/max(nil,1),
        late_capture_fraction_among_exits=ncl/max(ncl+nrl,1),
        late_resident_population_mean=meanpop,
        mean_completed_residence_rg_over_c=rt/max(nc+nr,1),
        capacity_exhausted=bool(over),elapsed_seconds=time.perf_counter()-start,
        scope='Fixed-geometry prescribed-input inner response only. No finite-halo growth history, BH feedback, or demonstrated stationary response.',
        physical_branch_selection_computed=False)
    label=f'{case}_c{cells}_dt{dt:g}_T{total_time:g}_N{nin}_x{sigma:g}_seed{seed}'
    if export_packets:
        label+='_packets'
        columns=['injection_time_tg','input_killing_energy','input_angular_momentum_rg_c',
                 'state_code','final_time_tg','final_killing_energy',
                 'final_angular_momentum_rg_c','final_radius_rg',
                 'final_pg_radial_momentum_mc','collision_count','injected_flag']
        packet_file=label+'.npz'
        np.savez_compressed(output_dir/packet_file,records=records[:ni],columns=columns,
            packet_mass_msun=weight,rg_pc=rg,tg_myr=tg,matching_rg=xmatch,
            horizon_absorber_rg=inner,
            scope='Prescribed-input fixed-geometry inner response; resident records are right censored.')
        row['packet_file']=packet_file
    (output_dir/(label+'.json')).write_text(json.dumps(row,indent=2))
    print(json.dumps(row,indent=2),flush=True)
    return row


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--case',choices=['compact','extended'],default='compact')
    ap.add_argument('--output-dir',default='outputs/inner')
    ap.add_argument('--distribution',help='Optional regenerated distribution CSV')
    ap.add_argument('--cells',type=int,default=64)
    ap.add_argument('--dt',type=float,default=.5)
    ap.add_argument('--duration',type=float,default=12000.)
    ap.add_argument('--packets',type=int,default=6000)
    ap.add_argument('--sigma',type=float,default=1.)
    ap.add_argument('--seed',type=int,default=20261030)
    ap.add_argument('--export-packets',action='store_true')
    a=ap.parse_args();run(a.case,a.cells,a.dt,a.duration,a.packets,a.sigma,a.seed,
                         export_packets=a.export_packets,output_dir=a.output_dir,distribution=a.distribution)
