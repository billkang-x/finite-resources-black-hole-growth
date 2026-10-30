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
