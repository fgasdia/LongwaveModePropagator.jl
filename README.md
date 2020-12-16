# Longwave Mode Propagator

Calculate the electric field produced by a very low frequency (VLF) radio transmitter in the ![Earth-ionosphere waveguide](https://en.wikipedia.org/wiki/Earth%E2%80%93ionosphere_waveguide) using the waveguide mode theory developed by Budden [^Budden1955], [^Budden1962].

`LongwaveModePropagator.jl` is similar in function to the U.S. Navy's Long-Wave Propagation Model in the Long-Wavelength Propagation Capability (LWPC) [^Ferguson1998], but can be more easily extended in part because it is written in the ![Julia](https://julialang.org/) programming language.

## Running the model

To get the most out of this package, interface with it through your own Julia project. (link to examples?)

## File interface

Describe JSON files. Show example from Matlab?

## New to Julia?

- What is julia ("compiles" what it can, multiple dispatch; build in Pkg)
- Download Julia.
- How to install this pkg.
- How to instantiate.
- Will to longer to run the first time because everything needs to be (pre)compiled. Also mention can be convenient to keep it open when running small individual jobs
- Finding help

## References

[^Budden1955]: K. G. Budden, “The numerical solution of differential equations governing reflexion of long radio waves from the ionosphere,” Proc. R. Soc. Lond. A, vol. 227, no. 1171, pp. 516–537, Feb. 1955, doi: 10.1098/rspa.1955.0027.

[^Budden1962]: K. G. Budden, “The influence of the earth’s magnetic field on radio propagation by wave-guide modes,” Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences, vol. 265, no. 1323, pp. 538–553, Feb. 1962, doi: 10.1098/rspa.1962.0041.

[^Ferguson1998]: J. A. Ferguson, “Computer programs for assessment of long-wavelength radio communications, version 2.0: User’s guide and source files,” Space and Naval Warfare Systems Center, San Diego, CA, Technical Document 3030, May 1998. [Online]. Available: http://www.dtic.mil/docs/citations/ADA350375.

## Citing

We encourage you to cite this package if used in scientific work. See the Zenodo
badge above or refer to [CITATION.bib](CITATION.bib).
