"""Recreate the revised Letter's continuous residence-enhancement figure."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from .joint_constraints import build_report, incidence_cap


def make_figure(output_dir, report=None):
    report=build_report() if report is None else report
    clock=report['comparison_clock_myr']; fuel=report['luminous_fuel_rate_msun_yr']
    rate=report['models'][0]['a_sigma_msun_yr']; observed=report['data']['F_mle']
    rq=np.linspace(.5,3.,300)
    def cap(r,factor=1.):
        return incidence_cap(rate*factor,clock,residence=float(r),fuel_rate=fuel,accounting='disjoint_sector')
    nominal=np.array([cap(x) for x in rq])
    conservative=np.array([cap(x,.8**3) for x in rq])
    # Exact crossing for the disjoint-sector formula, before its covering cap.
    crossing=rate/(1e9/(observed*clock*1e6)-fuel/.5)
    with plt.rc_context({'font.size':8,'axes.labelsize':9,'legend.fontsize':7}):
        fig,ax=plt.subplots(figsize=(3.35,2.35),layout='constrained')
        ax.plot(rq,100*nominal,color='#2468a0',label='Fiducial')
        ax.plot(rq,100*conservative,color='#49a66a',label='Conservative corner')
        ax.axhline(100*observed,color='black',ls='--',lw=1.,label=fr'CEERS proxy: {100*observed:.3f}%')
        for x in (1.,2.):
            ax.plot(x,100*cap(x),'o',color='#2468a0',ms=3)
        ax.axvline(crossing,color='#2468a0',ls=':',lw=.8)
        box=dict(facecolor='white',edgecolor='none',alpha=.95,pad=1.5)
        ax.annotate(fr'$R_Q\simeq {crossing:.2f}$',xy=(crossing,100*observed),xytext=(1.68,1.15),
                    bbox=box,arrowprops=dict(arrowstyle='-',color='#2468a0',lw=.7))
        ax.annotate(fr'$R_Q=1$: {100*cap(1):.3f}%',xy=(1,100*cap(1)),xytext=(.72,.75),
                    bbox=box,arrowprops=dict(arrowstyle='-',color='#2468a0',lw=.7))
        ax.annotate(fr'$R_Q=2$: {100*cap(2):.3f}%',xy=(2,100*cap(2)),xytext=(1.9,9.15),
                    bbox=box,arrowprops=dict(arrowstyle='-',color='#2468a0',lw=.7))
        ax.set(xlim=(.5,3),ylim=(0,20),xlabel=r'Residence enhancement $R_Q$',ylabel=r'$F_{\max}$ (%)')
        ax.grid(alpha=.18)
        ax.legend(loc='upper left',framealpha=1.,facecolor='white')
        out=Path(output_dir);out.mkdir(parents=True,exist_ok=True)
        for suffix in ('pdf','png'):
            fig.savefig(out/f'residence_scaling.{suffix}',dpi=240)
        plt.close(fig)
    return dict(crossing_R_Q=crossing,nominal_cap_RQ1=cap(1),nominal_cap_RQ2=cap(2),
                conservative_cap_RQ1=cap(1,.8**3),F_proxy=observed)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/revision'))
    args=ap.parse_args();print(make_figure(args.output_dir))


if __name__=='__main__':
    main()
