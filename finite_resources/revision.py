"""One-command audit of the additional revised-Letter analysis."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
import scipy

from .joint_constraints import ROOT, build_report, write_report
from .residence_history import demo_histories, integrate_history


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',type=Path,default=Path('outputs/revision'))
    ap.add_argument('--selection-json',type=Path,help='Optional reconstructed FITS selection; otherwise aggregate inputs are used')
    ap.add_argument('--skip-figure',action='store_true')
    args=ap.parse_args()
    manifest=json.loads((ROOT/'INPUT_SHA256.json').read_text(encoding='utf-8'))
    for name,digest in manifest.items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise SystemExit('Input checksum failed: '+name)
    selection=json.loads(args.selection_json.read_text(encoding='utf-8')) if args.selection_json else None
    report=build_report(selection)
    write_report(report,args.output_dir)
    histories={name:integrate_history(h) for name,h in demo_histories().items()}
    (args.output_dir/'residence_histories.json').write_text(json.dumps(histories,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    base=next(row for row in report['scenarios'] if row['profile']=='Fiducial' and row['column_branch']=='nominal'
              and row['R_Q']==1 and row['external_msun']==0 and row['p_g']==1 and row['K_sampling']==1)
    controls={
        'all_input_checksums':True,
        'quadrature_outward_RQ_half':math.isclose(histories['one_outward_pass']['R_Q'],.5,rel_tol=1e-6),
        'subdivision_does_not_create_exposure':math.isclose(histories['outward_episode_subdivision']['Q_s'],
                                                          histories['continuous_same_outward_path']['Q_s'],rel_tol=1e-12),
        'return_increases_RQ_without_new_mass':histories['returned_same_parcel']['R_Q']>1 and
                      histories['returned_same_parcel']['distinct_represented_g']==histories['one_outward_pass']['distinct_represented_g'],
        'new_supply_is_counted_once':histories['external_fresh_parcel']['distinct_external_g']==histories['one_outward_pass']['distinct_initial_g'],
    }
    if selection is None:
        controls.update(reference_count=(report['data']['N']==1242 and report['data']['K']==48),
                        reference_cap=math.isclose(base['cap_disjoint_sector'],.03114623936067968,rel_tol=1e-12))
    figure=None
    if not args.skip_figure:
        from .residence_scaling import make_figure
        figure=make_figure(args.output_dir,report)
    result=dict(passed=all(controls.values()),checks=controls,data=report['data'],fiducial=base,
                figure=figure,software=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
                scope='Aggregate-count quick reproduction, not a rerun of the raw FITS selection or a production RHD calculation')
    (args.output_dir/'revision_audit.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not result['passed']:
        raise SystemExit(1)


if __name__=='__main__':
    main()
