#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: prey
!
! !INTERFACE:
module mizer_prey
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_prey
      integer :: nclass
      type (type_bottom_state_variable_id) :: id_Nw
contains
      procedure :: initialize
   end type type_prey
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Initialise the module
!
! !INTERFACE:
   subroutine initialize(self, configunit)
!
! !INPUT PARAMETERS:
   class (type_prey), intent(inout), target :: self
   integer,           intent(in )           :: configunit
!
! !LOCAL VARIABLES:
   real(rk) :: w
!EOP
!-----------------------------------------------------------------------
!BOC
   call self%register_state_variable(self%id_Nw, 'Nw', 'g m-3', 'biomass')
   call self%get_parameter(w, 'w', 'g', 'mass per individual')
   call self%set_variable_property(self%id_Nw, 'particle_mass', w)

   end subroutine initialize
!EOC

end module mizer_prey

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2016
!-----------------------------------------------------------------------
