# IceShellXtal
This code simulates a freezing intrusion into an ocean world's ice shell:

- Nature of the solids (ice and salts) that precipitate
- Spatial distribution of the ice and salts depending on assumptions about whether the chamber liquid is mixed faster than salts settle
- Composition of the residual fluid, including water activity and solution density and conductivity.

IceShellXtal uses [PHREEQC](https://www.usgs.gov/software/phreeqc-version-3) (Parkhurst & Appelo 2013) to carry out the compositional calculations, with two databases described in my [frezchem](https://github.com/MarcNeveu/frezchem) repository. The PHREEQC output is analyzed in terms of solid volumes and, in the case of no mixing, the complex spatial distribution is plotted using [SDL2](https://www.libsdl.org) libraries. Running IceShellXtal requires installing the IPhreeqc modules, SDL2, and SDL2_image dependencies. IceShellXtal has only been developed and run on a Mac (OS 12.4 Monterey).
