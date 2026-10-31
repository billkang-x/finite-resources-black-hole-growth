# Small fixed inputs

Each case contains its final saved density profile (only the final time slice),
gas radius and total dark mass, and the 640-bin nonnegative isotropic energy
distribution used in the Letter. CSV columns retain 17-digit floating precision.
`model_parameters.json` states units and shared physical assumptions.

`f_mass` is mass per spatial volume per velocity-space volume; energy is
Newtonian specific energy in (pc/Myr)^2, with Phi(outer wall)=0. The wall reflects,
so positive energies are admitted. The density interior to the saved profile is
continued as r^(-3/4); this is part of the declared preparation.

These are small numerical inputs, not particle dumps. Halo/gas formation histories,
full particle records and large binary quadratures are omitted. The distribution
can optionally be refitted by `python -m finite_resources.phase_space`.

`ceers_proxy_reference.json` adds only aggregate audit inputs for the revised
Letter: public data URLs, compressed-catalogue SHA256, exact proxy cuts and
counts in five comparison windows. It contains no individual source records,
coordinates or official completeness weights. Reconstruct the aggregates from
the separately downloaded CEERS FITS file with `finite_resources.ceers_selection`.
