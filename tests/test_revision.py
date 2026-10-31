import copy
import math
import unittest

import numpy as np
from scipy.stats import binom

from finite_resources.bounds import incidence_bound
from finite_resources.ceers_selection import BANDS, reconstruct, selection_masks
from finite_resources.joint_constraints import (MSUN,YEAR,build_report,count_summary,
    incidence_cap,mixture_profile_loglike)
from finite_resources.residence_history import demo_histories, integrate_history
from finite_resources.trajectory_qeff import extract_arrays


class SelectionControls(unittest.TestCase):
    def setUp(self):
        n=10
        self.table={key:np.full(n,value,dtype=float) for key,value in
                    [('STAR_FLAG',0),('BAD_REGION_FLAG',0),('LP_Z_MED',6),('LP_NBANDS',10),('LP_MASS_MED',9.5)]}
        for b,f in zip(BANDS,[100,100,150,150,500]):
            self.table[b+'_FLUX']=np.full(n,f,dtype=float)
            self.table[b+'_FLUXERR_EMP']=np.full(n,10.,dtype=float)
        self.table['F444W_FLUXERR_EMP'][1]=1000. # colour band, not a detection band
        self.table['F150W_FLUXERR_EMP'][2]=1000.
        self.table['STAR_FLAG'][3]=1
        self.table['LP_Z_MED'][4]=4.6
        self.table['LP_NBANDS'][5]=7
        self.table['F277W_FLUX'][6]=500. # stays in parent; fails red colour
        self.table['LP_MASS_MED'][7]=np.nan
        self.table['LP_MASS_MED'][8]=10.5
        self.table['LP_Z_MED'][9]=9.1

    def test_parent_denominator_and_detection_bands(self):
        masks,colour,_,_=selection_masks(self.table)
        self.assertEqual(np.flatnonzero(masks['parent']).tolist(),[0,1,6,8])
        self.assertEqual(np.flatnonzero(masks['strict5']).tolist(),[0,6,8])
        self.assertEqual(np.flatnonzero(masks['broad_parent']).tolist(),[0,1,4,6,8,9])
        self.assertFalse(colour[6])
        result,*_=reconstruct(self.table,bootstrap_samples=0)
        self.assertEqual(result['samples']['parent']['N'],4)
        self.assertEqual(result['samples']['parent']['K'],3)

    def test_empty_sample_is_explicit(self):
        self.table['STAR_FLAG'][:]=1
        with self.assertRaises(ValueError): reconstruct(self.table)


class StatisticalControls(unittest.TestCase):
    def test_exact_interval_endpoints(self):
        self.assertEqual(count_summary(0,10)['exact_CP95'][0],0)
        self.assertEqual(count_summary(10,10)['exact_CP95'][1],1)
        with self.assertRaises(ValueError): count_summary(1,0)

    def test_profile_against_direct_maximization(self):
        k,n=4,100
        fs=np.array([0.,.02,.04,.08,.12])
        actual=mixture_profile_loglike(fs,k,n,cap=.1)
        for i,x in enumerate(fs[:-1]):
            q=np.unique(np.r_[np.linspace(x,1,20001),max(x,k/n)])
            expected=np.max(binom.logpmf(k,n,q))
            self.assertAlmostEqual(actual[i],expected,places=12)
        self.assertTrue(np.isneginf(actual[-1]))
        self.assertEqual(actual[0],actual[1]) # origin mixture is not identified

    def test_general_bound_matches_original_and_scales(self):
        cap=incidence_cap(51.3,625.8,residence=2.,fuel_rate=.002)
        old=incidence_bound(1.,1e9,625.8e6,51.3,fuel_cost_g_s=.002,residence_enhancement=2.)
        self.assertAlmostEqual(cap,old,places=14)
        self.assertAlmostEqual(incidence_cap(51.3,625.8,budget_msun=2e9,residence=2.,fuel_rate=.002),2*cap)
        self.assertLess(incidence_cap(51.3,625.8,residence=2.,fuel_rate=.002,accounting='disjoint_sector'),cap)
        with self.assertRaises(ValueError): incidence_cap(1,1,sampling_cap=.5)

    def test_manuscript_values_and_statistical_scope(self):
        report=build_report()
        row=next(r for r in report['scenarios'] if r['profile']=='Fiducial' and r['column_branch']=='nominal'
                 and r['R_Q']==1 and r['external_msun']==0 and r['p_g']==1 and r['K_sampling']==1)
        self.assertAlmostEqual(report['comparison_clock_myr'],625.8106358755068,places=10)
        self.assertAlmostEqual(row['cap_disjoint_sector'],.03114623936067968,places=14)
        self.assertFalse(row['central_value_satisfies_necessary_condition'])
        self.assertTrue(row['exact_CP95_contains_cap'])
        self.assertGreater(row['all_branch_binomial_upper_tail'],.05)


class HistoryControls(unittest.TestCase):
    def test_analytic_outward_and_distinct_mass(self):
        cases=demo_histories(); rows={k:integrate_history(v) for k,v in cases.items()}
        self.assertAlmostEqual(rows['short_residence']['R_Q'],1-1/1.1,places=7)
        self.assertAlmostEqual(rows['one_outward_pass']['R_Q'],.5,places=7)
        self.assertAlmostEqual(rows['returned_same_parcel']['R_Q'],1.5,places=6)
        self.assertEqual(rows['returned_same_parcel']['distinct_represented_g'],rows['one_outward_pass']['distinct_represented_g'])
        self.assertAlmostEqual(rows['outward_episode_subdivision']['Q_s'],rows['continuous_same_outward_path']['Q_s'],delta=1e-12*rows['continuous_same_outward_path']['Q_s'])
        self.assertAlmostEqual(rows['external_fresh_parcel']['R_Q'],rows['one_outward_pass']['R_Q'])
        self.assertEqual(rows['external_fresh_parcel']['distinct_represented_g'],2*rows['one_outward_pass']['distinct_represented_g'])

    def test_invalid_relabeling_overlap_and_budget(self):
        h=demo_histories()['one_outward_pass']
        duplicate=copy.deepcopy(h); duplicate['parcels'].append(copy.deepcopy(h['parcels'][0]))
        with self.assertRaises(ValueError): integrate_history(duplicate)
        overlap=copy.deepcopy(h); overlap['parcels'][0]['segments']*=2
        with self.assertRaises(ValueError): integrate_history(overlap)
        h['available_budget_g']=1.
        with self.assertRaises(ValueError): integrate_history(h)


class SnapshotControls(unittest.TestCase):
    def test_outward_mask_is_shared_by_both_integrals(self):
        edges=np.array([1.,2.,4.,8.])
        rho=np.array([[2.],[3.],[4.]])
        # Middle cell is stationary; it contributes to neither integral.
        velocity=np.array([[2e5],[0.],[4e5]])
        result=extract_arrays(np.stack([rho,rho*velocity]),edges,[-1.,1.])
        self.assertEqual(result['sigma_weighted_g_cm2'],2.*1.+4.*4.)
        self.assertAlmostEqual(result['q_eff_weighted_s_cm2'],.5/2e5+.125/4e5,places=15)

    def test_steady_finite_layer_and_angular_cancellation(self):
        # Analytic constant-speed rho=A/r^2 cell averages. Both finite integrals
        # cancel for any outer radius or represented solid angle.
        a,v,r0=2e18,7.46e7,5.22e16
        for outer in (1.01,2.,100.):
            edges=np.geomspace(r0,r0*outer,33)
            for mu in (np.array([-1.,0.,1.]),np.array([.2,.6,.9])):
                rho=np.repeat((a/(edges[:-1]*edges[1:]))[:,None],len(mu)-1,axis=1)
                state=np.stack([rho,rho*v])
                result=extract_arrays(state,edges,mu)
                self.assertAlmostEqual(result['a_sigma_msun_yr']/(4*math.pi*a*v*YEAR/MSUN),1.,places=12)
                expected=(mu[-1]-mu[0])/2*(1/r0-1/(r0*outer))/v
                self.assertAlmostEqual(result['q_eff_weighted_s_cm2']/expected,1.,places=12)
        with self.assertRaises(ValueError): extract_arrays(state,edges,mu,threshold_km_s=1000.)


if __name__=='__main__':
    unittest.main()
