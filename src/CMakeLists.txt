add_library(fabm_models_mizer OBJECT
            model_library.F90
            size_structured_base.F90
            size_structured_population.F90
            resource_spectrum.F90
            prey.F90
            multi_element_mizer.F90
            multi_element_support.F90
	    multi_element_mizer_demersal.F90
            multi_element_mizer_demersal_pelagic.F90
	    multi_element_mizer_demersal_pelagic_bsdo.F90
           )

add_dependencies(fabm_models_mizer fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2015-2016 - PML
#-----------------------------------------------------------------------
