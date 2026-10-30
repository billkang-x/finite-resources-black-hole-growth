# Provenance and scope

Code by Dong-Biao Kang. This release accompanies the revised Letter
*Finite Resources Limit Obscuration and Capture in Black-Hole Growth*.
The release is a portable extraction of the finite-residence/gas accounting,
nonnegative energy-distribution fit, geodesic-filtered collisionless inventory,
fixed-Schwarzschild inner DSMC and their diagnostics. Changes for portability
replace local paths with packaged CSV inputs and user-chosen output folders;
the transport kernels and orbit quadrature formulae are retained.

Small profiles and energy distributions are the same prescribed preparations
used for the reported results. They are not observed halo data or cosmological
initial conditions. `INPUT_SHA256.json` identifies all shipped inputs and saved
reference summaries. Python/library and platform changes can change stochastic
streams; compare conservation, analytic controls and ensemble diagnostics.

The inner experiment lasts 96000 GM/c^3, about 13.1 hours for the 10^5 Msun seed.
The 100 Myr energy comparison holds input fixed; it is not an evolved halo.
Earlier experimental long-term boundary routines were excluded after audits
identified invalid time-unit and returned-angular-momentum mappings. They are
not supplied and none of their endpoints enter the revised manuscript.

Gas residence and dark-matter capture use distinct material ledgers. These
calculations do not by themselves exclude all SIDM seeds or all LRD models.

Related accepted article (background, not the DOI of this code or the Letter):
Dong-Biao Kang, *Seed-mass dependence of baryon-assisted black-hole growth in
self-interacting dark matter haloes*, The Astrophysical Journal (accepted),
DOI: https://doi.org/10.3847/1538-4357/aeb5e1 (author-supplied accepted metadata).

AI tools assisted with language polishing and code writing. The author takes
responsibility for the scientific content and the final manuscript.
