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
