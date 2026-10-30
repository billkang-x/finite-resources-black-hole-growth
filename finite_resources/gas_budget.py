from pathlib import Path
import json, math
import argparse
ap=argparse.ArgumentParser(description='Gas-column, finite-layer and stationary-retention checks.')
ap.add_argument('--output-dir',default='outputs')
a=ap.parse_args();R=Path(a.output_dir);R.mkdir(parents=True,exist_ok=True)
c=2.99792458e10; msun=1.98847e33; year=31557600; pc=3.0856775814913673e18
Ledd=1.09e44/.75
mdot_edd=Ledd/(.1*c*c)*year/msun
rows=[]
for name,M6,dotm in [('Fiducial',1.,2000),('GTO-2974',.23,3084),('JADES-28074',2.8,1933),('UNCOVER-45924',7.6,1544)]:
    radius=5.22e16*M6**(10/11)*(dotm/2000)**(51/22)
    speed=746e5*(dotm/2000)**(-.5)
    rate=(dotm-1)*mdot_edd*M6
    column=rate*msun/year/(4*math.pi*speed*radius)
    bill=rate*50e6
    rows.append(dict(name=name,M6=M6,dotm=dotm,radius_cm=radius,speed_cm_s=speed,column_g_cm2=column,rate_msun_yr=rate,bill_msun=bill,lower_bill_msun=bill*.8**3))
rho=rows[0]['column_g_cm2']/rows[0]['radius_cm']
retention=[]
for dotm in [1500,2000,5000]:
    fuel=1e6*(math.e-1)/.9; wind=.5*(dotm-1)*fuel
    retention.append(dict(dotm=dotm,retained_fraction=.9/(1+.5*(dotm-1)),mass_for_one_efold=fuel+wind,minimum_reused_fraction=1-(1e9-fuel)/wind))
finite_layer=[]
for factor in [1.01,1.5,2,10,1000]:
    kernel=(1-1/factor)/(rows[0]['speed_cm_s']*rows[0]['radius_cm'])
    column=rows[0]['column_g_cm2']*(1-1/factor)
    cost=4*math.pi*column/kernel*year/msun
    assert math.isclose(cost,rows[0]['rate_msun_yr'],rel_tol=1e-12)
    finite_layer.append(dict(R_over_r=factor,cost_msun_yr=cost))
assert rows[0]['lower_bill_msun']>1e9
assert rows[1]['bill_msun']<1e9
out=dict(constants=dict(c=c,solar_mass_g=msun,year_seconds=year,pc_cm=pc,Ledd_erg_s=Ledd),model_rows=rows,density_normalization=dict(continuity_g_cm3=rho,printed_g_cm3=1.51e-15,printed_over_continuity=1.51e-15/rho),stationary_growth=retention,finite_layer=finite_layer)
(R/'gas_budget_audit.json').write_text(json.dumps(out,indent=2),encoding='utf8')
print(json.dumps(out,indent=2))
