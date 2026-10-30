"""Frozen spherical preparations. Units are pc, solar masses, and Myr.

The gas disk is represented by its spherical monopole. The tabulated final
density is continued inward with gamma=3/4, normalized to the stated dark
mass, and supported by the original finite reflecting-wall Jeans pressure.
This specifies an input state; it does not simulate formation of that state.
"""
from pathlib import Path
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid
INPUTS=Path(__file__).resolve().parents[1]/'inputs'
KMS=1.022712165045695
G=4.300917270036279e-3*KMS**2
C=299792.458*KMS
S=1.988409870698051e33/3.085677581491367e18**2
M0=1e5

def canonical_case(name):
    aliases={'finite_disk_compact_slow':'compact','finite_disk_fine':'extended'}
    name=aliases.get(name,name)
    if name not in ('compact','extended'):raise ValueError('case must be compact or extended')
    return name

def load_distribution(case,filename=None):
    path=Path(filename) if filename else INPUTS/(canonical_case(case)+'_distribution.csv')
    a=np.genfromtxt(path,delimiter=',',names=True)
    if not all(np.all(np.isfinite(a[n])) for n in ['energy_low','energy_high','f_mass']):
        raise ValueError('Energy distribution must be finite')
    if np.any(a['f_mass']<0) or not np.all(a['energy_high']>a['energy_low']):
        raise ValueError('Energy distribution must have increasing bins and nonnegative mass density')
    if not np.array_equal(a['energy_high'][:-1],a['energy_low'][1:]):
        raise ValueError('Energy bins must be contiguous')
    return {'energy_edges':np.r_[a['energy_low'],a['energy_high'][-1]],'f':a['f_mass']}

def nfw_fn(x):
    x=np.asarray(x)
    return np.where(x<1e-4,x*x*(.5-x*(2/3-x*.75)),np.log1p(x)-x/(1+x))

class Core:
    def __init__(self,name,n=5000,gamma=.75):
        name=canonical_case(name)
        self.name=name
        self.gamma=gamma
        self.meta=json.loads((INPUTS/f'{name}.json').read_text())
        a=np.genfromtxt(INPUTS/f'{name}_profiles.csv',delimiter=',',names=True)
        a=a[a['t_myr']==a['t_myr'].max()]
        self.input=a
        self.R=self.meta['final']['gas_scale_pc']
        self.rvir=G*5e10/(2*(76.51034218019993*KMS)**2)
        self.rs=self.rvir/3
        self.mdiff=5e10*.049/.315-1e8-M0
        self.r=np.geomspace(1e-12,self.rvir,n)
        lr=np.log(a['r_pc']);rho=np.exp(np.interp(np.log(self.r),lr,np.log(a['rho_msun_pc3'])))
        inside=self.r<a['r_pc'][0]
        rho[inside]=a['rho_msun_pc3'][0]*(self.r[inside]/a['r_pc'][0])**(-gamma)
        dm=4*np.pi*rho*self.r**3
        md=np.r_[0,cumulative_trapezoid(dm,np.log(self.r))]
        md+=4*np.pi*rho[0]*self.r[0]**3/(3-gamma)
        self.normalization=self.meta['final']['M_dark_msun']/md[-1]
        self.rho=rho*self.normalization;self.md=md*self.normalization
        self.wall_pressure=self.rho[-1]*(76.51034218019993*KMS)**2

    def baryons(self,r):
        x=np.asarray(r)/self.R;h=.6;xx=np.minimum(x,np.sqrt(1+h*h))
        mg=1e8*2/h*(np.where(xx<=h,xx**3/3,.5*h*xx**2-h**3/6)-np.maximum(xx*xx-1,0)**1.5/3)
        xxn=np.asarray(r)/self.rs
        mdiff=self.mdiff*nfw_fn(xxn)/nfw_fn(3.)
        low=np.sqrt(np.maximum(1-1/np.maximum(x,1e-300)**2,0))
        frac=np.maximum(np.minimum(1,h/np.maximum(x,1e-300))-low,0)
        deriv=1e8*2/h*x**3*frac+self.mdiff*xxn**2/(1+xxn)**2/nfw_fn(3.)
        return mg+mdiff,deriv

    def state(self,M,response):
        r=self.r;rho=self.rho;md=self.md
        if response!='fixed':
            raise ValueError('Only the frozen preparation is part of this release')
        bn,_=self.baryons(r)
        mt=M+md+bn
        acc=G*mt/r**2
        # Integrate from the outer boundary inward to avoid cancellation.
        ii=-cumulative_trapezoid((rho*acc)[::-1],r[::-1],initial=0)[::-1]
        sig2=(self.wall_pressure+ii)/rho
        ph=-cumulative_trapezoid(acc[::-1],r[::-1],initial=0)[::-1]+G*mt[-1]/r[-1]
        # The Kepler rho~r^-3/4 cusp has beta=1/4 and vmax^2=3.5 sigma^2.
        # A finite hot reflecting halo wall needs the same positive support
        # extension when its Jeans pressure includes formally unbound speeds.
        vmax2=np.maximum(2*ph,3.5*sig2)
        beta=vmax2/(2*sig2)-1.5
        return r,rho,np.sqrt(sig2),beta,md,vmax2

