This is an implementation of [mizer](http://dx.doi.org/10.1111/2041-210X.12256) in [FABM](https://fabm.net).

It contains separate modules for a single size-structured population ("size_structured_population") and a resource spectrum ("resource_spectrum").
Multiple size-structured populations may be added to fabm.yaml to create a multi-species community. It also includes a multi-element version of
the size structured population ("multi_element_mizer"), which explicitly treats carbon, nitrogen and phosphorus (plus, in prey, silicon). This can be coupled to
element-resolving prey models such as ERSEM and ensures mass of all different elements is conserved.

To include mizer modules during a FABM build, add the following to your call to CMake: `-DFABM_INSTITUTES="<OTHER_MODELS>;mizer" -DFABM_MIZER_BASE=<MIZERDIR>`.
Here, `<OTHER_MODELS>` is a semi-colon-separated list of additional FABM packages you want to build, e.g., ersem (quotes are _required_ to prevent your shell (e.g., bash)
from interpreting the semi-colon as "end of command"). `<MIZERDIR>` is the root of the mizer repository.

# pyfabm

The tools to perform mizer simulations driven by existing outputs from lower trophic level (LTL) models
are written in Python and available in the `python` subdirectory. To use these, you first need to build [pyfabm](https://github.com/fabm-model/fabm/wiki/python)
(the Python interface to FABM).

On Linux you can build pyfabm with mizer support as follows:

    mkdir -p ~/build/pyfabm && cd ~/build/pyfabm
    cmake <FABMDIR>/src/drivers/python -DFABM_INSTITUTES="mizer" -DFABM_MIZER_BASE=<MIZERDIR>
    make install

# stand-alone box model

Similarly, you can build FABM's 0d driver with mizer and ERSEM support as follows:

    mkdir -p ~/build/fabm0d && cd ~/build/fabm0d
    cmake <FABMDIR>/src/drivers/0d -DGOTM_BASE=<GOTMDIR> -DFABM_INSTITUTES="ersem;mizer" -DFABM_ERSEM_BASE=<ERSEMDIR> -DFABM_MIZER_BASE=<MIZERDIR>
    make install

This will give you an executable at `~/local/fabm/0d/bin/fabm0d`.

A test case for FABM's 0d driver is incldued at testcases/mizer_0d.
