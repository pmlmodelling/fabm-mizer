This is an implementation of [mizer](http://dx.doi.org/10.1111/2041-210X.12256) in [FABM](http://fabm.net).

It contains separate modules for a single size-structured population ("size_structured_population") and a resource spectrum ("resource_spectrum").
Multiple size-structured populations may be added to fabm.yaml to create a multi-species community.

To include mizer modules during a FABM build, add the following to your call to CMake: `-DFABM_INSTITUTES="<OTHER_MODELS>;mizer" -DFABM_MIZER_BASE=<MIZERDIR>`.
Here, `<OTHER_MODELS>` is a semi-colon-separated list of additional FABM packages you want to build, e.g., ersem (quotes are _required_ to prevent your shell (e.g., bash)
from interpreting the semi-colon as "end of command"). `<MIZERDIR>` is the root of the mizer repository.

For instance, on Linux you can build FABM's 0d driver with mizer and ERSEM support as follows:

    mkdir -p ~/build/fabm0d && cd ~/build/fabm0d
    cmake <FABMDIR>/src/drivers/0d -DGOTM_BASE=<GOTMDIR> -DFABM_INSTITUTES="ersem;mizer" -DFABM_ERSEM_BASE=<ERSEMDIR> -DFABM_MIZER_BASE=<MIZERDIR>
    make install

This will give you an executable at `~/local/fabm/0d/bin/fabm0d`.

A test case for FABM's 0d driver is incldued at testcases/mizer_0d.
