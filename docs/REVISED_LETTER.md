# Additional analysis in the revised Letter (v1.1.0)

This release adds code for *A Finite-Residence Bound on Obscured Black-Hole
Growth*. It preserves the v1.0.0 capture calculation and inputs. Manuscript
sources, PDFs, raw survey tables and wind snapshots are not distributed.

## Public catalogue reconstruction

Download `ceers_cat_v1.0.fits.gz` from the public
[CEERS DR1 catalogue directory](https://web.corral.tacc.utexas.edu/ceersdata/DR1/Catalog/).
The exact compressed-input SHA256, release description links, predicates and
aggregate results are in `inputs/ceers_proxy_reference.json`. Source descriptions:

- [Catalogue readme](https://web.corral.tacc.utexas.edu/ceersdata/DR1/Catalog/ceers_cat.readme.md)
- [NIRCam readme](https://web.corral.tacc.utexas.edu/ceersdata/DR1/NIRCam/README.txt)
- [Catalogue paper](https://arxiv.org/abs/2510.08743)
- [Public reduction code](https://github.com/ceers/ceers-nircam)

The analyst-defined parent uses clean star/bad-region flags, 5 <= LP_Z_MED <= 9,
at least eight fitted bands, finite LP mass, empirical S/N >= 3 in F150W,
F200W, F277W and F356W, and F444W <= 27.5 AB. The numerator then applies
F277W-F444W > 1 and F150W-F200W < 0.5. The colour cut does not define its own
denominator. The four-band parent has 1242 rows and 48 candidates; the strict
five-band control has 1241 and 48. AB magnitudes use nJy zero point 31.4.

The optional row bootstrap samples its exact Binomial(N,K/N) count law. This
has the same resampling distribution as bootstrapping binary row indicators,
although random realizations differ from the earlier index-array implementation.
The weighted ratio uses product Phi(SNR_i-3), conditional on measured catalogue
fluxes. It is a noise-margin sensitivity diagnostic, not an inverse-completeness
weight. The catalogue operator is reproducible; official latent-source
injection/recovery completeness, size/crowding corrections and cosmic variance
are not supplied by this implementation.

Use `--allow-different-catalogue` only for an explicitly new input release;
the output records a mismatch. Object-level generated CSVs are local outputs,
never repository inputs. The software license does not relicense CEERS data.
Please follow the data providers' citation instructions.

## Bounds, posterior and branch comparison

`joint_constraints.py` reconstructs the four outflow examples directly from
the published continuity prescription. Their parameters, constants and
threshold factors are contained in the code. It does not read the PRD paper's
local tables or import its manuscript. The density coefficient is the smaller
continuity value, not the printed coefficient that is 19.3% higher.

Rates a_col and a_fuel are in solar masses/year, B is in solar masses, and T
is converted from Myr to years. K_sampling is distinct from candidate count K.
The comparison clock uses flat matter+Lambda with H0=67.4 km/s/Mpc and
Omega_m=0.315, giving T(z=9 to 5)=625.8106359 Myr. The luminosity example uses
L_min=1e43 erg/s and efficiency 0.1. The two explicitly separated formulas are:

1. **General conservative bound:**
   F <= p_g min(1, K_sampling B / [T max(a_col/R_Q, a_fuel)]).
   The fuel and column constraints are independent; this permits their gas
   ledgers to overlap. No fixed covering is required.
2. **Disjoint-sector application:**
   F <= p_g min(C, K_sampling B / [T (a_col/R_Q + a_fuel/C)]), C=0.5.
   This additionally adopts disjoint wind/fuel streams, fixed covering
   during activity, uniform viewing directions and bounded temporal sampling.
   It reproduces the full-precision 3.114623936% coefficient
   used for the paper's rounded 3.115% Fiducial comparison. The general bound
   is reported separately; their small fuel correction is not silently equated.

Reducing Sigma, r and v by 20% simultaneously reduces a_col by 0.8^3; this is a
declared sensitivity corner, not three measured independent errors. External
supply enlarges the distinct budget. R_Q changes residence; it is not a return
fraction or a new gas budget. The branch table scans both these operations,
growth eligibility and bounded sampling.

The selected-count likelihood is Binomial(K|N,F_total), where
F_total=F_SIDM+F_other for disjoint origins. Only total incidence is identified.
The code reports a Beta(K+1,N-K+1) posterior **given a uniform prior on total
incidence**, Clopper-Pearson intervals and an analytically profiled likelihood
over F_other. This is not a posterior for physical SIDM parameters. With
F_SIDM fixed, the optimal total incidence is max(K/N,F_SIDM), subject to the
physical cap. No large two-dimensional likelihood grid is needed.

For N=1242, K=48, the central fraction is 3.8647343%, while its exact 95% count
interval is [2.8630452%,5.0916453%]. The Fiducial cap lies inside that interval.
The central-value necessary-condition failure is therefore **not a 95%
statistical exclusion**. Tail probabilities are reported for clarity; no
universal SIDM exclusion or official completeness-corrected abundance follows.

## Labeled trajectory histories

`residence_history.py` accepts JSON with this structure (cgs units):

```json
{
  "reference": {"speed_cm_s": 100000000, "radius_cm": 1e18, "opacity_cm2_g": 1},
  "parcels": [
    {"id": "gas-a", "mass_g": 1e33, "source": "initial",
     "segments": [
       {"time_s": [0, 1e9, 2e9], "radius_cm": [1e18, 1.1e18, 1.2e18], "opacity_cm2_g": 1}
     ]}
  ]
}
```

Every parcel ID appears once, with all its contributing visits in ordered,
nonoverlapping segments. Gaps correspond to time outside the selected region.
For each parcel q = integral kappa dt/r^2 has units s/g; Q=sum(m q) has units s.
The reported R_Q=Q v_ref r_ref/(kappa_ref B) uses distinct available mass B.
By default B is the sum of represented initial and external masses; optional
`available_budget_g` can include unused material, but cannot undercount parcels.
External parcels are counted once. Returning visits add q, not distinct mass.

The controls compare a short outward passage, a longer passage, returning reuse,
subdivision of an unchanged outward path, and an additional external parcel.
They are prescribed kinematic trajectories. Relabeling a luminous episode
cannot reset an outward radius. The module does not infer support, forces,
thermal balance, or overlap with the qualifying luminosity mask g(t).

## Optional wind-snapshot adapter and the resource scale

`trajectory_qeff.py` reads user-supplied NPZ fields:

- `state`: shape (field, radial_cell, angular_ray), density in field 0 and
  radial momentum density in field 1, in cgs units;
- `r_edges`: positive increasing radial edges in cm;
- `mu_edges`: increasing cos(theta) edges in [-1,1].

The same cells with v_r above the stated positive threshold enter
Sigma_ray=sum(rho dr) and q_ray=sum[(1/r_lo-1/r_hi)/v_r]. The latter is exact
for cellwise constant velocity; dr/r_geometric^2 is algebraically identical.
Angular weights are delta_mu/2, without renormalizing partial sky coverage.
The output reports a_col=4 pi <Sigma>/<q> and
Lambda=a_col H_eff/B_eff, with user-specified exposure and distinct budget.
It rejects invalid arrays and does not unpickle input objects.

This reconstructs radial streamlines of one Eulerian snapshot. It is not the
mass-weighted, time-resolved Q of tracked parcels, nor a converged production
RHD coefficient. The original local states were artificial neutral-wind
benchmarks; they are excluded from this release. In tests, a constant-speed
rho=A/r^2 wind gives a_col=4 pi A v independent of the matched finite outer
radius and represented solid angle. A trajectory-derived physical R_Q should
instead use labeled histories or a justified transport reconstruction.

## Source mapping and portability changes

| Local research implementation | Public implementation |
|---|---|
| `public_selection_joint_reconstruction_v1.py` | `ceers_selection.py` |
| `public_selection_joint_report_v1.py`, `public_joint_likelihood_v1.py` | `joint_constraints.py` |
| `history_kernel_qeff_v1.py` | `residence_history.py`, `trajectory_qeff.py` |
| `make_residence_scaling_figure.py` | `residence_scaling.py` |

Absolute machine paths and manuscript-table dependencies were replaced with
CLI inputs, exact formulas and a small hashed aggregate reference. Snapshot
and history APIs explicitly separate Eulerian and labeled interpretations.
The count mixture is profiled analytically, the row bootstrap is sampled from
its equivalent count law, and the plot computes its coefficients at full
precision. No physical conclusion is inferred from whether a software test passes.
