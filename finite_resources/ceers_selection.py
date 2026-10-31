"""Reconstruct a CEERS DR1 colour proxy from a user-supplied public FITS file.

The four detection bands and F444W colour/flux band have different roles. The
Gaussian noise-margin weights are diagnostics, not completeness corrections.
Astropy is only required by the CLI's FITS loader, not by the other modules.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import ndtr

from .joint_constraints import ROOT, count_summary

BANDS = ('F150W','F200W','F277W','F356W','F444W')
FIELDS = ['NUMBER','RA','DEC','STAR_FLAG','BAD_REGION_FLAG','LP_Z_MED','LP_NBANDS','LP_MASS_MED']
FIELDS += [b+s for b in BANDS for s in ('_FLUX','_FLUXERR_EMP')]


def sha256(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1<<20),b''):
            h.update(block)
    return h.hexdigest()


def selection_masks(table):
    """Pure vectorized predicates; parent selection NEVER includes colour cuts."""
    f=np.column_stack([np.asarray(table[b+'_FLUX'],dtype=float) for b in BANDS])
    err=np.column_stack([np.asarray(table[b+'_FLUXERR_EMP'],dtype=float) for b in BANDS])
    finite=np.isfinite(f)&np.isfinite(err)&(f>0)&(err>0)
    snr=np.divide(f,err,out=np.zeros_like(f),where=np.isfinite(err)&(err>0))
    z=np.asarray(table['LP_Z_MED'],dtype=float)
    mass=np.asarray(table['LP_MASS_MED'],dtype=float)
    clean=(np.asarray(table['STAR_FLAG'])==0)&(np.asarray(table['BAD_REGION_FLAG'])==0)
    base=clean & np.isfinite(z) & (np.asarray(table['LP_NBANDS'])>=8) & np.isfinite(mass)
    mag=np.full(len(z),np.inf)
    positive=np.isfinite(f[:,4])&(f[:,4]>0)
    mag[positive]=31.4-2.5*np.log10(f[positive,4])
    quality4=np.all(finite[:,:4]&(snr[:,:4]>=3),axis=1)
    quality5=np.all(finite&(snr>=3),axis=1)
    base=base&(mag<=27.5)
    parent=base&(z>=5)&(z<=9)&quality4
    strict5=base&(z>=5)&(z<=9)&quality5
    broad=base&(z>=4.5)&(z<=9.5)&quality4
    # Flux-ratio tests implement AB colour inequalities without log rounding.
    colour=np.all(np.isfinite(f[:,[0,1,2,4]]),axis=1)&(f[:,4]>0)&(f[:,1]>0)
    colour=colour&(f[:,2]<f[:,4]*10**(-.4))&(f[:,0]>f[:,1]*10**(-.2))
    weights=np.prod(ndtr(snr[:,:4]-3),axis=1)
    weights[~np.all(np.isfinite(snr[:,:4]),axis=1)]=0.
    masks=dict(parent=parent,broad_parent=broad,mass_parent=parent&(mass>9)&(mass<10),
               strict5=strict5,strict5_mass_parent=strict5&(mass>9)&(mass<10))
    return masks,colour,weights,dict(mag444=mag,min_snr=np.min(snr,axis=1))


def reconstruct(table, bootstrap_samples=4096, seed=20261031):
    if bootstrap_samples < 0:
        raise ValueError('bootstrap count cannot be negative')
    masks,colour,weights,diagnostics=selection_masks(table)
    rng=np.random.default_rng(seed)
    samples={}
    for name,mask in masks.items():
        n=int(np.count_nonzero(mask)); k=int(np.count_nonzero(mask&colour))
        if n == 0:
            raise ValueError('empty comparison sample: '+name)
        stats=count_summary(k,n)
        sw=float(np.sum(weights[mask])); sk=float(np.sum(weights[mask&colour]))
        row=dict(N=n,K=k,q_count=k/n,count_CP95=stats['exact_CP95'],
                 weighted_N=sw,weighted_K=sk,q_weighted=sk/sw if sw else None)
        if bootstrap_samples:
            # Row bootstrap of a binary indicator has exactly Binomial(n,k/n)
            # counts. This samples the same law without large index matrices.
            b=rng.binomial(n,k/n,size=bootstrap_samples)/n
            row['row_bootstrap95']=[float(x) for x in np.quantile(b,[.025,.975])]
        samples[name]=row
    ref=json.loads((ROOT/'inputs/ceers_proxy_reference.json').read_text(encoding='utf-8'))
    return dict(catalogue_rows=len(colour),samples=samples,
                public_sources=ref['public_sources'],public_rules=ref['public_rules'],
                bootstrap=dict(samples=bootstrap_samples,seed=seed,
                    scope='Count resampling only; excludes cosmic variance and latent-source incompleteness'),
                weight_scope='Product Phi(SNR_i-3) in four detection bands, conditional on measured fluxes; not inverse-completeness weights',
                scope='Public CEERS catalogue hard-cut colour proxy; no official injection/recovery selection function'),masks,colour,weights,diagnostics


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('catalogue',type=Path,help='Path to public ceers_cat_v1.0.fits.gz')
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/ceers'))
    ap.add_argument('--bootstrap-samples',type=int,default=4096)
    ap.add_argument('--seed',type=int,default=20261031)
    ap.add_argument('--allow-different-catalogue',action='store_true',help='Allow a hash mismatch; record it explicitly')
    args=ap.parse_args()
    try:
        from astropy.io import fits
    except ImportError as exc:
        raise SystemExit('Install FITS support with: python -m pip install -r requirements-observation.txt') from exc
    ref=json.loads((ROOT/'inputs/ceers_proxy_reference.json').read_text(encoding='utf-8'))
    digest=sha256(args.catalogue)
    matched=digest==ref['catalogue_sha256']
    if not matched and not args.allow_different_catalogue:
        raise SystemExit('Catalogue SHA256 differs from the paper input. Use --allow-different-catalogue for a declared new analysis.')
    with fits.open(args.catalogue,memmap=False) as hdul:
        table={key:np.asarray(hdul[1].data[key]) for key in FIELDS}
    report,masks,colour,weights,diagnostics=reconstruct(table,args.bootstrap_samples,args.seed)
    report.update(catalogue_filename=args.catalogue.name,catalogue_sha256=digest,
                  matches_paper_catalogue=matched)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    (args.output_dir/'public_selection_reconstruction.json').write_text(
        json.dumps(report,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    # These generated catalogue-level records stay local under ignored outputs/.
    selected=np.logical_or.reduce(list(masks.values()))
    with (args.output_dir/'public_parent_records.csv').open('w',newline='',encoding='utf-8') as f:
        fields=['source_id','z','logMstar','mag444','min_snr','pdet4_proxy','color_candidate',*masks]
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for j in np.flatnonzero(selected):
            row=dict(source_id=int(table['NUMBER'][j]),z=float(table['LP_Z_MED'][j]),
                     logMstar=float(table['LP_MASS_MED'][j]),mag444=float(diagnostics['mag444'][j]),
                     min_snr=float(diagnostics['min_snr'][j]),pdet4_proxy=float(weights[j]),color_candidate=bool(colour[j]))
            row.update({k:bool(v[j]) for k,v in masks.items()});w.writerow(row)
    print(json.dumps(dict(catalogue_rows=report['catalogue_rows'],samples=report['samples'],
                         matches_paper_catalogue=matched,output=str(args.output_dir)),indent=2))


if __name__ == '__main__':
    main()
