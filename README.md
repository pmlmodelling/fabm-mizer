This is an implementation of [mizer](http://dx.doi.org/10.1111/2041-210X.12256), a dynamic size-structured model for aquatic predators, in [FABM](https://fabm.net).
It can be used to model the size spectrum of the higher trophic level community (e.g., all biomass between 1 mg and 1000 kg) or the size structure of individual populations.
It is typically used to describe fish, and accordingly supports several different parameterisatons of fisheries.

## Included modules

This codebase contains separate modules for a single size-structured population ("size_structured_population") and a resource spectrum ("resource_spectrum").
Multiple size-structured populations may be added to fabm.yaml to create a multi-species community. It also includes a multi-element version of
the size structured population ("multi_element_mizer"), which explicitly treats carbon, nitrogen and phosphorus (and, in prey, silicon). This can be coupled to
element-resolving prey models such as [ERSEM](https://ersem.com/) and ensures mass of all different elements is conserved.

## Offline simulation

The `python` subdirectory includes tools for offline simulations with mizer, driven with outputs (environmental conditions and prey)
of lower trophic level or biogeochemical models. Such outputs usually come from simulations with 1D or 3D models, for instance, GOTM or NEMO.

To use the tools for offline simulation, you need pyfabm (the Python interface to FABM).
To build and install pyfabm, please follow the [instructions on the FABM web site](https://github.com/fabm-model/fabm/wiki/python).
 **Note:** You will need to explicitly activate mizer support within pyfabm by providing `-DFABM_INSTITUTES="<OTHER_MODELS>;mizer" -DFABM_MIZER_BASE=<MIZERDIR>` to cmake. 

## Two-way coupled simulation

Through FABM, these modules can be embedded in variety of hydrodynamic models and combined with models for lower trophic levels to build a two-way coupled end-to-end ecosystem model.
For more information on using FABM with various hydrodynamic models, please visit the [FABM wiki](https://fabm.net/wiki).

To include mizer modules during a FABM build, add the following to your call to CMake: `-DFABM_INSTITUTES="<OTHER_MODELS>;mizer" -DFABM_MIZER_BASE=<MIZERDIR>`.
Here, `<OTHER_MODELS>` is a semi-colon-separated list of additional FABM packages you want to build, e.g., ersem.
NB the double quotes are _required_ to prevent your shell (e.g., bash) from interpreting the semi-colon as "end of command".
`<MIZERDIR>` is the root of the mizer repository.
