"""Finite residence in cgs units; covering-weighted duration is in seconds.

All mass elements have nonnegative mass. For a nonreturning outward history,
opacity <= kappa_max, r >= r_min and dr/dt >= v_min imply
q = integral(kappa dt/r^2) <= kappa_max/(v_min*r_min).
"""
import math
import numpy as np


def column_cost(column_g_cm2, speed_cm_s, radius_cm):
    """Minimum distinct gas mass per covering-weighted second, in g/s."""
    if min(column_g_cm2,speed_cm_s,radius_cm)<=0:
        raise ValueError('Positive column, speed and radius required')
    return 4*math.pi*column_g_cm2*speed_cm_s*radius_cm


def outward_kernel(kappa_cm2_g,speed_cm_s,radius_cm,outer_radius_cm=math.inf):
    """Maximal opacity-weighted column-time per gram for the specified layer."""
    if min(kappa_cm2_g,speed_cm_s,radius_cm)<=0 or outer_radius_cm<=radius_cm:
        raise ValueError('Positive parameters and outer radius > inner radius required')
    return kappa_cm2_g/speed_cm_s*(1/radius_cm-1/outer_radius_cm)


def incidence_bound(eligible_fraction,budget_g,duration_s,column_cost_g_s,
                    fuel_cost_g_s=0.,sampling_cap=1.,residence_enhancement=1.):
    """Joint growth-obscuration incidence bound for homogeneous eligible hosts.

    sampling_cap=K means normalized time weight w <= K/T. Fuel cost is
    retained-growth fuel per qualifying covering-weighted time. The maximum
    permits sequential reuse of the same gas; use an additive cost only when
    distinct wind and fuel streams have separately been established.
    """
    if not 0<=eligible_fraction<=1 or budget_g<0:
        raise ValueError('Invalid host fraction or budget')
    if min(duration_s,column_cost_g_s,residence_enhancement)<=0 or sampling_cap<1 or fuel_cost_g_s<0:
        raise ValueError('Invalid duration, cost, sampling cap or residence enhancement')
    cost=max(column_cost_g_s/residence_enhancement,fuel_cost_g_s)
    return eligible_fraction*min(1.,sampling_cap*budget_g/(duration_s*cost))


def trajectory_kernel(times_s,radii_cm,kappa_cm2_g):
    """Numerical q for one resolved trajectory, including any actual returns.

    Times must increase. This quadrature alone does not certify a speed floor
    between nodes. Check temporal resolution for curved/returning trajectories.
    """
    t=np.asarray(times_s,float);r=np.asarray(radii_cm,float)
    k=np.broadcast_to(np.asarray(kappa_cm2_g,float),t.shape)
    if t.ndim!=1 or len(t)<2 or t.shape!=r.shape or np.any(np.diff(t)<=0) or np.any(r<=0) or np.any(k<0):
        raise ValueError('Invalid trajectory')
    from scipy.integrate import trapezoid
    return float(trapezoid(k/r**2,t))


def residence_enhancement(masses_g,kernels,opacity_ceiling_cm2_g,speed_floor_cm_s,radius_floor_cm):
    """R_Q of distinct mass elements; repeated passages belong in each q."""
    m=np.asarray(masses_g,float);q=np.asarray(kernels,float)
    if m.shape!=q.shape or np.any(m<0) or np.any(q<0) or m.sum()<=0:
        raise ValueError('Nonnegative distinct masses and kernels required')
    if min(opacity_ceiling_cm2_g,speed_floor_cm_s,radius_floor_cm)<=0:
        raise ValueError('Positive comparison floors required')
    return float(np.dot(m,q)*speed_floor_cm_s*radius_floor_cm/(opacity_ceiling_cm2_g*m.sum()))
