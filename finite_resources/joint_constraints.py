"""Count likelihood and explicitly conditional finite-resource comparisons.

Rates are in Msun/yr and time in Myr. The general bound permits reuse between
fuel and column ledgers. The separate disjoint-sector application adds their
costs, assumes covering C during activity, and reproduces the manuscript table.
Neither likelihood identifies SIDM origin from a colour count alone.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom

ROOT = Path(__file__).resolve().parents[1]
MSUN = 1.98847e33
YEAR = 31557600.0
C_LIGHT = 2.99792458e10


def count_summary(k, n):
    if not isinstance(n, (int, np.integer)) or not isinstance(k, (int, np.integer)) or not 0 <= k <= n or n <= 0:
        raise ValueError('counts must be integers with 0 <= k <= n and n > 0')
    interval = [0.0 if k == 0 else float(beta.ppf(.025, k, n-k+1)),
                1.0 if k == n else float(beta.ppf(.975, k+1, n-k))]
    return dict(N=int(n), K=int(k), F_mle=k/n, exact_CP95=interval,
                uniform_prior_beta95=[float(x) for x in beta.ppf([.025, .5, .975], k+1, n-k+1)],
                posterior_prior='uniform on total selected incidence F; not on SIDM occupation')


def cosmic_age_myr(z, h0=67.4, omega_m=.315):
    """Flat matter+Lambda comparison clock; radiation omitted, as in the paper."""
    if z < 0 or h0 <= 0 or not 0 < omega_m < 1:
        raise ValueError('invalid cosmology')
    hubble = h0 * 1e5 / 3.0856775814913673e24
    return 2/(3*hubble*math.sqrt(1-omega_m))*math.asinh(
        math.sqrt((1-omega_m)/omega_m)/(1+z)**1.5)/(YEAR*1e6)


def outflow_models():
    """Liu et al. continuity prescription, NOT its larger printed rho coefficient."""
    rate_edd = (1.09e44/.75)/(.1*C_LIGHT**2)*YEAR/MSUN
    rows = []
    for name, m6, dotm in [('Fiducial',1.,2000), ('GTO-2974',.23,3084),
                            ('JADES-28074',2.8,1933), ('UNCOVER-45924',7.6,1544)]:
        radius = 5.22e16*m6**(10/11)*(dotm/2000)**(51/22)
        speed = 746e5*(dotm/2000)**(-.5)
        rate = (dotm-1)*rate_edd*m6
        rows.append(dict(name=name, M6=m6, dotm=dotm, radius_cm=radius,
                         speed_cm_s=speed, a_sigma_msun_yr=rate,
                         sigma_g_cm2=rate*MSUN/YEAR/(4*math.pi*speed*radius)))
    return rows


def incidence_cap(a_sigma, duration_myr, budget_msun=1e9, residence=1.,
                  eligible_fraction=1., sampling_cap=1., external_msun=0.,
                  fuel_rate=0., covering=.5, accounting='general'):
    vals = (a_sigma,duration_myr,budget_msun,residence,eligible_fraction,
            sampling_cap,external_msun,fuel_rate,covering)
    if not all(math.isfinite(x) for x in vals):
        raise ValueError('parameters must be finite')
    if min(a_sigma,duration_myr,residence) <= 0 or min(budget_msun,external_msun,fuel_rate) < 0:
        raise ValueError('invalid rates, time, residence or budget')
    if not 0 <= eligible_fraction <= 1 or sampling_cap < 1 or not 0 < covering <= 1:
        raise ValueError('invalid population, sampling cap or covering')
    if accounting == 'general':
        cost = max(a_sigma/residence, fuel_rate)
        ceiling = 1.
    elif accounting == 'disjoint_sector':
        cost = a_sigma/residence + fuel_rate/covering
        ceiling = covering
    else:
        raise ValueError('accounting must be general or disjoint_sector')
    return eligible_fraction*min(ceiling, sampling_cap*(budget_msun+external_msun)/
                                  (duration_myr*1e6*cost))


def mixture_profile_loglike(f_sidm, k, n, cap=1.):
    """Analytically maximize over 0 <= F_other <= 1-F_SIDM.

    At fixed F_SIDM the optimum total is max(K/N, F_SIDM). This eliminates
    the discretization error and memory cost of a two-dimensional grid.
    """
    count_summary(k, n)
    fs = np.asarray(f_sidm, dtype=float)
    if not 0 <= cap <= 1:
        raise ValueError('cap must lie in [0,1]')
    valid = np.isfinite(fs) & (fs >= 0) & (fs <= cap)
    q = np.clip(np.maximum(k/n, fs), 0., 1.)
    return np.where(valid, binom.logpmf(k, n, q), -np.inf)


def build_report(selection=None):
    ref = json.loads((ROOT/'inputs/ceers_proxy_reference.json').read_text(encoding='utf-8'))
    parent = (selection or ref)['samples']['parent']
    k, n = int(parent['K']), int(parent['N'])
    counts = count_summary(k,n)
    clock = cosmic_age_myr(5)-cosmic_age_myr(9)
    fuel = 1e43/(.1*C_LIGHT**2)*YEAR/MSUN
    scenarios = []
    for model in outflow_models():
        for column, factor in [('nominal',1.),('conservative_corner',.8**3)]:
            for rq in (.5,1.,2.):
                for external in (0.,1e9):
                    for pg in (1.,.5,.1):
                        for sampling in (1.,2.,5.):
                            rate = model['a_sigma_msun_yr']*factor
                            kwargs = dict(residence=rq,external_msun=external,
                                          eligible_fraction=pg,sampling_cap=sampling,fuel_rate=fuel)
                            cap = incidence_cap(rate,clock,accounting='disjoint_sector',**kwargs)
                            cap_general = incidence_cap(rate,clock,accounting='general',**kwargs)
                            scenarios.append(dict(profile=model['name'],column_branch=column,
                                a_sigma_msun_yr=rate,R_Q=rq,external_msun=external,p_g=pg,
                                K_sampling=sampling,cap_disjoint_sector=cap,cap_general=cap_general,
                                central_value_satisfies_necessary_condition=counts['F_mle'] <= cap,
                                maximum_SIDM_share_at_central_incidence=min(1.,cap/counts['F_mle']) if k else 1.,
                                exact_CP95_contains_cap=counts['exact_CP95'][0] <= cap <= counts['exact_CP95'][1],
                                all_branch_binomial_upper_tail=float(binom.sf(k-1,n,cap)),
                                uniform_total_prior_probability_above_cap=float(beta.sf(cap,k+1,n-k+1))))
    return dict(data=counts,comparison_clock_myr=clock,
                luminous_fuel_rate_msun_yr=fuel,models=outflow_models(),scenarios=scenarios,
                general_bound='p_g min(1, K_sampling B / [T max(a_sigma/R_Q,a_fuel)])',
                disjoint_sector_application='p_g min(C, K_sampling B / [T (a_sigma/R_Q+a_fuel/C)]); C=0.5',
                disjoint_sector_assumptions='Disjoint fuel and column gas, fixed covering during activity, uniform viewing directions and bounded temporal sampling',
                interpretation='Central-value necessary-condition comparison; the fiducial cap lies inside the 95% count interval. This is not a 95% statistical rejection.',
                identifiability='The count likelihood identifies only F_total=F_SIDM+F_other. A bound on one declared branch does not identify SIDM occupation.',
                selection_scope='Public catalogue hard-cut colour proxy; no official injection/recovery completeness correction.',
                reference_source=ref['public_sources'])


def write_report(report, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    (out/'joint_constraints.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    with (out/'branch_scan.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(report['scenarios'][0])); w.writeheader(); w.writerows(report['scenarios'])
    k,n=report['data']['K'],report['data']['N']
    fs=np.linspace(0, .15, 301)
    with (out/'mixture_profile_likelihood.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['profile','column_branch','R_Q','F_SIDM','profile_loglike']
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader()
        for row in report['scenarios']:
            if row['p_g'] != 1 or row['K_sampling'] != 1 or row['external_msun'] != 0:
                continue
            ll=mixture_profile_loglike(fs,k,n,row['cap_disjoint_sector'])
            for x,y in zip(fs,ll):
                w.writerow(dict(profile=row['profile'],column_branch=row['column_branch'],R_Q=row['R_Q'],
                                F_SIDM=float(x),profile_loglike=float(y)))


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--selection-json',type=Path,help='Output of ceers_selection; omit to use shipped aggregate counts')
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/ceers'))
    args=ap.parse_args()
    selection=json.loads(args.selection_json.read_text(encoding='utf-8')) if args.selection_json else None
    report=build_report(selection)
    write_report(report,args.output_dir)
    print(json.dumps(dict(data=report['data'],comparison_clock_myr=report['comparison_clock_myr'],
                         output=str(args.output_dir)),indent=2))


if __name__ == '__main__':
    main()
