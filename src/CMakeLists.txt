add_library(fabm_models_mizer OBJECT
            model_library.F90
            size_structured_population.F90
            resource_spectrum.F90
           )

add_dependencies(fabm_models_mizer fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2015-2016 - PML
#-----------------------------------------------------------------------
