# Release validation (2026-10-30)

The quick audit, figure reconstruction, nonnegative 640-bin distribution refit,
relativistic pair and coordinate tests,
six analytic-versus-adaptive geodesic controls, finite-layer and physical-time
controls, a 6000-packet local scattering run, and both full 12/16/384 orbit
quadratures were run from this portable repository.

The independent compact and extended orbit inventories reproduce the recorded
1.5185165615102135e-3 and 5.59031485657213e-4 solar masses, respectively.
The compact independent spatial phase volume differs by 3.75e-7 relatively.
The 6000-packet run has zero rest-count residual; its unexplained relative energy
residual after measured drifts is approximately 2.1e-15. This small run is a
runtime check, not a replacement for the reported refined ensemble.

The validation CLI tests 5000 pair collisions (energy/momentum errors below
3.3e-14), 5000 frame round trips (error below 3.8e-14), and six geodesics.
Its finite-layer cancellation error is below 3e-16; t_g=0.492549 seconds for
the adopted 100000-solar-mass seed. All independent controls pass.

Tested library versions:

- Python: 3.13.5
- NumPy: 2.3.1
- SciPy: 1.16.0
- Numba: 0.66.0
- Matplotlib: 3.11.0

## Revised-Letter additions (2026-10-31)

The portable FITS selection was rerun on the full 87,345-row public CEERS DR1
catalogue with the recorded SHA256. It reproduces all five archived sample
counts and weighted ratios, including parent N=1242, K=48 and strict-five-band
N=1241, K=48. The raw catalogue and generated object records remain local.

`python -m unittest discover -s tests -v` passes 10 independent controls:
synthetic parent/numerator separation and detection-band sensitivity; empty
sample handling; exact interval endpoints; mixture profiling against direct
maximization; agreement with the original general-bound implementation and
budget scaling; manuscript values and their count-statistical scope; analytic
outward history, returning reuse and distinct-supply accounting; duplicate-ID,
overlap and underbudget rejection; matched outward-cell masking; and analytic
steady-wind cancellation across finite radial layers and angular coverage.

The revision quick audit and the complete FITS-to-joint-analysis workflow pass.
The comparison clock is 625.8106359 Myr. The Fiducial disjoint-sector cap is
0.03114623936068 and the general maximum-cost cap is 0.03114838352262.
The all-branch binomial upper-tail diagnostic is 0.07834; its central-value
necessary-condition failure is not a 95% statistical rejection. The origin
mixture remains unidentifiable from the count alone.

The optional NPZ adapter was also run on the local saved 64x64 artificial
neutral-wind benchmark at a 1 km/s outward threshold, giving
q_eff=4.948026889e-28 s/cm^2, weighted Sigma=0.1430220032 g/cm^2,
a_col=57.6455243 Msun/yr and Lambda=2.88227621 at 50 Myr and 1e9 Msun.
This checks adapter execution; it does not establish a converged production
RHD coefficient or a time-resolved R_Q. No snapshot is included.

The original `finite_resources.reproduce` command still passes and verifies
all numerical input hashes plus the new aggregate reference. The new scaling
figure was rendered and visually checked. Original capture integrations were
not rerun because their code and inputs have not changed.
