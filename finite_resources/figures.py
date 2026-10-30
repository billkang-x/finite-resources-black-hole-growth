from pathlib import Path
import json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
ap=argparse.ArgumentParser(description='Plot saved lightweight summary products; no transport rerun.')
ap.add_argument('--output-dir',default='outputs/figures')
a=ap.parse_args();P=Path(a.output_dir);P.mkdir(parents=True,exist_ok=True)
D=Path(__file__).resolve().parents[1]/'summary'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':6.7,'pdf.fonttype':42})
fig,(ax,bx)=plt.subplots(2,1,figsize=(3.42,3.95),gridspec_kw={'height_ratios':[1.1,1]},layout='constrained')
h=np.geomspace(1,120,180)
for name,rate,col in [('Fiducial',51.3005,'#111111'),('GTO-2974',18.20,'#0072B2'),('JADES-28074',138.83,'#D55E00'),('UNCOVER-45924',300.95,'#009E73')]:
 y=rate*h/1000
 ax.plot(h,y,lw=1.35,color=col,label=name)
 ax.fill_between(h,.512*y,y,color=col,alpha=.12,lw=0)
ax.axhline(1,color='.45',ls='--',lw=.8);ax.axvline(50,color='.45',ls=':',lw=.8)
ax.set(xscale='log',yscale='log',xlim=(1,120),ylim=(.01,50),
 xlabel=r'Covering-weighted duration $H$ (Myr)',ylabel=r'Required gas ($10^9\,M_\odot$)')
ax.legend(loc='upper left',frameon=False,ncol=1)
ax.text(.97,.06,'(a)',transform=ax.transAxes,ha='right',weight='bold')
suite=json.loads((D/'inner_suite.json').read_text())
names=['collisionless','baseline','space','time','particles','extended']
labels=['No\nscatter','128\ncells','256\ncells','Half\nstep','Double\npackets','Extended\ncore']
means=[]
for i,n in enumerate(names):
 vals=np.array(suite['groups'][n]['independent_seed_values'])*100;means.append(float(vals.mean()))
 color='#0072B2' if n in ['collisionless','extended'] else '#D55E00'
 bx.bar(i,vals.mean(),width=.65,color=color,alpha=.65,zorder=2)
 bx.scatter(i+np.array([-.12,.12]),vals,color='black',s=9,zorder=3)
bx.axhline(100*suite['compact_collisionless_analytic_fraction'],color='.3',ls='--',lw=.8,zorder=1)
bx.set(xticks=range(6),xticklabels=labels,ylabel='Late capture / injected (%)',ylim=(0,17))
bx.text(.97,.9,'(b)',transform=bx.transAxes,ha='right',weight='bold')
for aa in [ax,bx]:
 aa.spines[['top','right']].set_visible(False)
fig.savefig(P/'residence_capture.pdf');fig.savefig(P/'residence_capture.png',dpi=220)
plt.close(fig)
# Independent, finite-inventory curve: no exponential reset of radial phases.
fig,ax=plt.subplots(figsize=(5.7,2.8),layout='constrained')
for name,label,c in [('compact','Compact','#D55E00'),('extended','Extended','#0072B2')]:
 a=np.genfromtxt(D/f'{name}_inventory.csv',delimiter=',',names=True)
 meta=json.loads((D/f'{name}_inventory_metadata.json').read_text())
 t=a['time_myr'][1:];ax.plot(t,a['captured_rest_msun'][1:],label=label,color=c)
 ax.axhline(meta['inventory_msun'],color=c,lw=.8,ls=':')
ax.set(xscale='log',yscale='log',xlabel='Elapsed time (Myr)',ylabel=r'Collisionless capture ($M_\odot$)',ylim=(1e-8,.003))
ax.legend(frameon=False);ax.spines[['top','right']].set_visible(False)
fig.savefig(P/'finite_inventory.pdf');plt.close(fig)
fig,ax=plt.subplots(figsize=(3.42,1.95),layout='constrained')
for name,label,c in [('compact','Compact','#D55E00'),('extended','Extended','#0072B2')]:
 a=np.genfromtxt(D/f'{name}_inventory.csv',delimiter=',',names=True)
 meta=json.loads((D/f'{name}_inventory_metadata.json').read_text())
 t=a['time_myr'][1:];ax.plot(t,a['captured_rest_msun'][1:],label=label,color=c,lw=1.3)
 ax.axhline(meta['inventory_msun'],color=c,lw=.8,ls=':')
ax.set(xscale='log',yscale='log',xlabel='Elapsed time (Myr)',ylabel=r'Captured mass ($M_\odot$)',ylim=(1e-8,.003))
ax.legend(frameon=False,loc='lower right');ax.spines[['top','right']].set_visible(False)
fig.savefig(P/'finite_inventory_letter.pdf');plt.close(fig)
print('Figures generated from lightweight saved summary tables; transport not rerun.')
