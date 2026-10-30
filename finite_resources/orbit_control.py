"""Finite-reservoir collisionless control from the reconstructed global f(E).

Newtonian spherical orbits outside 100 GM/c^2 are matched to local static
Schwarzschild velocities at that surface. Only inward geodesics that cross
the relativistic potential barrier are absorbed; inward passage through the
matching sphere alone is not capture. The potential is held fixed for this
sigma/m=0 control. No collisional SIDM growth history is inferred here.

Units: pc, Msun, Myr. The stationary distribution is a mass distribution.
"""
from pathlib import Path
import json, time, argparse
from functools import lru_cache
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from numpy.polynomial.legendre import leggauss
from .profiles import Core,G,C,M0,S,load_distribution


@lru_cache(maxsize=16)
def nodes(n):
    return leggauss(n)


def bound_critical(e):
    """Unstable Schwarzschild circular orbit, including bound E<1.

    The matching surface at x=100 has E>=sqrt(.98)>sqrt(8/9),
    so the unstable barrier exists for every admitted local velocity.
    """
    e = np.asarray(e)
    if np.any(e < np.sqrt(8/9)):
        raise ValueError('No unstable circular barrier in supported branch')
    delta = 2/(3*e*e-2 + e*np.sqrt(9*e*e-8))
    x = 3+delta
    return x, x/np.sqrt(delta)


class OrbitControl:
    def __init__(self, name, ne=640, nr=200, distribution=None):
        self.core = Core(name)
        rr, rho, *_ = self.core.state(M0, 'fixed')
        a = G*(M0+self.core.md+self.core.baryons(rr)[0])/rr**2
        phi = cumulative_trapezoid(a[::-1], rr[::-1], initial=0)[::-1]
        self.phi_table = PchipInterpolator(np.log(rr), phi)
        self.rg = G*M0/C**2
        self.match = 100*self.rg
        self.phim = float(self.phi(self.match))
        self.data = load_distribution(name,distribution)
        self.column = cumulative_trapezoid(S*rho, rr, initial=0)
        self.column_table = PchipInterpolator(np.log(rr), self.column)

    def phi(self, r):
        return self.phi_table(np.log(r))

    def capture_limit(self, energy):
        speed2 = 2*(energy-self.phim)
        if np.any(speed2 < 0) or np.any(speed2 >= C*C):
            raise ValueError('Invalid matching speed')
        gamma = 1/np.sqrt(1-speed2/C**2)
        egr = np.sqrt(.98)*gamma
        _, ell = bound_critical(egr)
        j2 = (self.rg*C*ell/gamma)**2
        j2 = np.minimum(j2, self.match**2*speed2)
        return j2, egr, speed2

    def outer(self, energy, j2):
        def radial2(r):
            return 2*(energy-float(self.phi(r)))-j2/r**2
        if radial2(self.core.rvir) >= 0:
            return self.core.rvir
        return brentq(radial2, self.match, self.core.rvir,
                      xtol=self.match*1e-16, rtol=4*np.finfo(float).eps)

    def radial_integral(self, energy, j2, n):
        r, wdr, v2, outer = self.quadrature(energy,j2,n)
        if np.min(v2) <= 0:
            raise ArithmeticError(f'Nonpositive radial speed: E={energy}, j2={j2}, outer={outer}, minv2={v2.min()}, r={r[np.argmin(v2)]}, match={self.match}, phim={self.phim}')
        period = 2*np.sum(wdr/np.sqrt(v2))
        # Jensen's inequality for the initially isotropic background gives
        # <|u-v|> >= |u| >= |u_r|. This is a lower bound on the expected
        # collision count during a first inward leg, not a trajectory history.
        optical_lower = float(self.column_table(np.log(outer))-
                              self.column_table(np.log(self.match)))
        return period, outer, optical_lower

    def quadrature(self, energy, j2, n):
        outer=self.outer(energy,j2)
        y,w=nodes(n);y=(y+1)/2;w=w/2
        distance=(outer-self.match)*(1-y)**2
        r=outer-distance
        # Near a turning point r and outer may round to the same float.
        # Retain the quadrature distance and subtract the local PCHIP cubic
        # analytically; never clip a negative radial kinetic energy.
        delta=-np.log1p(-distance/outer)
        logouter=np.log(outer)
        k=int(np.clip(np.searchsorted(self.phi_table.x,logouter)-1,0,
                      len(self.phi_table.x)-2))
        z0=logouter-self.phi_table.x[k];z=z0-delta
        c=self.phi_table.c[:,k]
        difference=delta*(c[0]*(z0*z0+z0*z+z*z)+c[1]*(z0+z)+c[2])
        crossed=z<0
        difference[crossed]=float(self.phi(outer))-self.phi(r[crossed])
        wall_v2=2*(energy-float(self.phi(outer)))-j2/outer**2 if outer==self.core.rvir else 0.
        v2=2*difference-j2/outer**2*np.expm1(2*delta)+wall_v2
        return r,w*2*(outer-self.match)*(1-y),v2,outer

    def spatial_inventory(self, energy, j2max, n):
        # Independent phase-volume integral over position and velocity angle.
        r,wdr,v2,outer=self.quadrature(energy,0.,n)
        if np.min(v2)<=0:raise ArithmeticError('Invalid spatial phase-volume node')
        v=np.sqrt(v2)
        z=np.minimum(j2max/(r*v)**2, 1.)
        fraction=z/(1+np.sqrt(1-z))
        return 16*np.pi**2*np.sum(wdr*r*r*v*fraction)


def weighted_quantile(x, w, q):
    order=np.argsort(x);x=np.asarray(x)[order];w=np.asarray(w)[order]
    return np.interp(q,np.cumsum(w)/np.sum(w),x)


def run_case(name, ne=640, nr=200, neq=12, nj=16, nt=384, output_dir="outputs/orbits", distribution=None):
    start=time.perf_counter()
    output_dir=Path(output_dir);output_dir.mkdir(parents=True,exist_ok=True)
    model=OrbitControl(name,ne,nr,distribution)
    e_edges=model.data['energy_edges']; f=model.data['f']
    ne=len(f)
    xe,we=nodes(neq);xj,wj=nodes(nj)
    flux=[]; periods=[]; horizons=[]; taus=[]; outers=[]; direct=[]
    energy_list=[];j2_list=[];sphere_flux=0.
    for low, high, fi in zip(e_edges[:-1],e_edges[1:],f):
        low=max(low,model.phim)
        if high<=low or fi<=0:
            continue
        energies=(low+high)/2+(high-low)/2*xe
        eweights=(high-low)/2*we*fi
        jmax,egr,v2=model.capture_limit(energies)
        sphere_flux+=np.sum(4*np.pi**2*eweights*model.match**2*v2)
        for en,ew,jm,eg in zip(energies,eweights,jmax,egr):
            direct.append(ew*model.spatial_inventory(en,jm,nt*2))
            for xx,ww in zip(xj,wj):
                jj=jm*(xx+1)/2
                ff=2*np.pi**2*ew*jm*ww
                period,outer,tau=model.radial_integral(en,jj,nt)
                periods.append(period);flux.append(ff);horizons.append(eg)
                taus.append(tau);outers.append(outer);energy_list.append(en);j2_list.append(jj)
    flux=np.array(flux);periods=np.array(periods);horizons=np.array(horizons)
    taus=np.array(taus);outers=np.array(outers)
    mass=flux*periods
    times=np.r_[0.,np.geomspace(1e-8,100.,180)]
    removed=np.minimum(times[:,None]/periods[None,:],1.)@mass
    remaining=np.maximum(1-times[:,None]/periods[None,:],0.)@mass
    rate=(times[:,None]<periods[None,:])@flux
    # At most one crossing of the absorbing relativistic barrier per element.
    # A uniform initial radial phase drains linearly until one radial period.
    # The negligible geodesic travel delay inside x=100 is omitted explicitly.
    tag=f'{name}_E{ne}_eq{neq}_j{nj}_t{nt}'
    np.savez_compressed(output_dir/f'{tag}.npz',time_myr=times,
        captured_rest_msun=removed,remaining_capture_orbit_rest_msun=remaining,
        capture_rate_msun_myr=rate,orbit_mass_msun=mass,orbit_period_myr=periods,
        orbit_flux_msun_myr=flux,energy_schwarzschild=horizons,
        first_infall_collision_count_lower=taus,outer_turning_radius_pc=outers,
        energy_newtonian=np.array(energy_list),angular_momentum_squared=np.array(j2_list))
    direct_mass=float(np.sum(direct)); inventory=float(mass.sum())
    result=dict(case=name,energy_bins=ne,energy_quadrature=neq,angular_quadrature=nj,
        radial_quadrature=nt,orbits=len(mass),gas_radius_pc=model.core.R,
        initial_capture_rate_rest_msun_myr=float(flux.sum()),
        all_inward_matching_sphere_rate_rest_msun_myr=float(sphere_flux),
        true_geodesic_fraction_of_inward_rate=float(flux.sum()/sphere_flux),
        initial_capture_orbit_inventory_rest_msun=inventory,
        independent_spatial_inventory_rest_msun=direct_mass,
        phase_volume_relative_discrepancy=abs(inventory/direct_mass-1),
        captured_rest_by_100myr_msun=float(removed[-1]),
        captured_seed_mass_fraction=float(removed[-1]/M0),
        time_to_drain_half_inventory_myr=float(weighted_quantile(periods,mass,.5)),
        # The period median is not itself the half-drain time; corrected below.
        mass_weighted_first_infall_collision_lower_median=float(weighted_quantile(taus,mass,.5)),
        mass_fraction_orbits_with_first_infall_collision_lower_ge_one=float(mass[taus>=1].sum()/inventory),
        flux_fraction_orbits_with_first_infall_collision_lower_ge_one=float(flux[taus>=1].sum()/flux.sum()),
        mass_balance_relative_max=float(np.max(abs(remaining+removed-inventory))/inventory),
        raw_file=tag+'.npz',wall_seconds=time.perf_counter()-start)
    result['surface']='geodesic_filtered_x100'
    result['mass_weighted_radial_period_median_myr']=result.pop('time_to_drain_half_inventory_myr')
    result['time_to_drain_half_inventory_myr']=float(brentq(
        lambda t:np.minimum(t/periods,1)@mass-.5*inventory,0,float(periods.max())))
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',choices=['compact','extended'],default='compact')
    ap.add_argument('--energy-order',type=int,default=12)
    ap.add_argument('--angular-order',type=int,default=16)
    ap.add_argument('--radial-order',type=int,default=384)
    ap.add_argument('--distribution')
    ap.add_argument('--output-dir',default='outputs/orbits')
    a=ap.parse_args()
    row=run_case(a.case,neq=a.energy_order,nj=a.angular_order,nt=a.radial_order,output_dir=a.output_dir,distribution=a.distribution)
    (Path(a.output_dir)/(a.case+'_orbit_summary.json')).write_text(json.dumps(row,indent=2))
    print(json.dumps(row,indent=2))

if __name__=='__main__':main()
