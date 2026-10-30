"""Quick audit: shipped inputs, gas budgets, analytic flux, saved summary curves."""
from pathlib import Path
import argparse,json,hashlib,subprocess,sys
import numpy as np
from .validation import flux_capture_reference
from .profiles import G,C,M0


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output-dir',default='outputs',type=Path)
    a=ap.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
    root=Path(__file__).resolve().parents[1]
    manifest=json.loads((root/'INPUT_SHA256.json').read_text())
    failures=[name for name,h in manifest.items() if hashlib.sha256((root/name).read_bytes()).hexdigest()!=h]
    if failures:raise RuntimeError('Input checksum mismatch: '+str(failures))
    subprocess.run([sys.executable,'-m','finite_resources.gas_budget','--output-dir',str(a.output_dir)],check=True,stdout=subprocess.DEVNULL)
    gas=json.loads((a.output_dir/'gas_budget_audit.json').read_text())
    suite=json.loads((root/'summary/inner_suite.json').read_text())
    energy=json.loads((root/'summary/capture_energy.json').read_text())
    controls=[flux_capture_reference(c) for c in ['compact','extended']]
    assert abs(controls[0]['collisionless_capture_fraction']-suite['compact_collisionless_analytic_fraction'])<1e-12
    assert abs(controls[0]['inward_flux_msun_myr']/energy['inner_records'][0]['rest_flux_msun_myr']-1)<1e-12
    assert abs(gas['model_rows'][0]['bill_msun']/2.565e9-1)<1e-4
    inventories={}
    for name in ['compact','extended']:
        tab=np.genfromtxt(root/f'summary/{name}_inventory.csv',delimiter=',',names=True)
        meta=json.loads((root/f'summary/{name}_inventory_metadata.json').read_text())
        assert np.all(np.diff(tab['captured_rest_msun'])>=0)
        assert np.allclose(tab['captured_rest_msun']+tab['remaining_rest_msun'],meta['inventory_msun'],rtol=3e-14,atol=1e-20)
        inventories[name]=dict(inventory_msun=meta['inventory_msun'],captured_100myr_msun=float(tab['captured_rest_msun'][-1]))
    baseline=suite['groups']['baseline']['mean_late_capture_per_injected']
    collisionless=suite['groups']['collisionless']['mean_late_capture_per_injected']
    eflux=float(np.mean([x['input_energy_flux_msun_myr'] for x in energy['inner_records']]))
    out=dict(passed=True,files_verified=len(manifest),fiducial_gas_bill_msun=gas['model_rows'][0]['bill_msun'],
        retained_fraction_dotm2000=gas['stationary_growth'][1]['retained_fraction'],
        input_distribution_controls=controls,inventory_summary_checks=inventories,
        saved_summary_capture_ratio=baseline/collisionless,
        inner_duration_myr=96000*G*M0/C**3,fixed_input_energy_100myr_msun=100*eflux,
        required_mean_energy_supply_amplification=(1e7-1e5)/(100*eflux),
        scope='Gas and analytic flux recomputed; packet and orbit summaries checked. No transport simulation rerun by this command.')
    (a.output_dir/'summary_check.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
