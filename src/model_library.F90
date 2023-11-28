module mizer_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   use mizer_size_structured_population
   use mizer_resource_spectrum
   use mizer_prey
   use mizer_multi_element_population
   use mizer_multi_element_demersal_pelagic_population
   use mizer_multi_element_dem_pelagic_population
   ! Add new mizer models here

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: mizer_model_factory

contains

   subroutine create(self,name,model)

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('size_structured_population'); allocate(type_size_structured_population::model)
         case ('multi_element_population');   allocate(type_multi_element_population::model)
         case ('resource_spectrum');          allocate(type_resource_spectrum::model)
         case ('prey');                       allocate(type_prey::model)
         case ('multi_element_demersal_pelagic_population');   allocate(type_multi_element_demersal_pelagic_population::model)
         case ('multi_element_dem_pelagic_population');   allocate(type_multi_element_dem_pelagic_population::model)
         ! Add new mizer models here
      end select

   end subroutine

end module
