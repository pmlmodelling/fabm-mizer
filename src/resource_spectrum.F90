#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: resource_spectrum
!
! !INTERFACE:
module mizer_resource_spectrum
!
! !DESCRIPTION:
! Size-structured resource spectrum growing to carrying capacity
!
! Based on resource in:
!
! Blanchard, J. L., Andersen, K. H., Scott, F., Hintzen, N. T., Piet, G., & Jennings, S. (2014)
! Evaluating targets and trade-offs among fisheries and conservation objectives using a multispecies size spectrum model
! Journal of Applied Ecology, 51(3), 612–622
! doi:10.1111/1365-2664.12238
!
! Symbols are taken equal to those used in the mizer R package.
!
! !USES:
   use fabm_types

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_resource_spectrum
      ! Variable identifiers
      type (type_bottom_state_variable_id),allocatable :: id_Nw(:)
      type (type_bottom_state_variable_id)             :: id_waste
      integer :: nclass
      real(rk), allocatable :: w(:)
      real(rk), allocatable :: r(:)
      real(rk), allocatable :: K(:)
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type type_resource_spectrum

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable),parameter :: total_mass = type_bulk_standard_variable(name='total_mass',units='g',aggregate_variable=.true.,conserved=.true.)
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
   subroutine initialize(self,configunit)
!
! !INPUT PARAMETERS:
   class (type_resource_spectrum), intent(inout), target :: self
   integer,                        intent(in )           :: configunit
!
! !LOCAL VARIABLES:
   integer              :: iclass
   real(rk)             :: delta_logw, logw
   real(rk)             :: kappa, lambda, r0, w_min, w_max, n
   real(rk),allocatable :: logw_bounds(:)
   real(rk),parameter   :: sec_per_year = 86400 * 365.2425_rk
   character(len=10)    :: strindex
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read parameters
   ! All rate coefficients are converted from d-1 to s-1 ensure source terms are directly compatible with FABM.
   ! Default values taken from Table S5
   call self%get_parameter(self%nclass,'nclass', '-',           'number of size classes', default=30)
   call self%get_parameter(w_min,      'w_min',  'g',           'minimum size', minimum=0.0_rk)
   call self%get_parameter(w_max,      'w_max',  'g',           'maximum size', minimum=0.0_rk)
   call self%get_parameter(r0,         'r0',     'g^(1-n) yr-1','pre-factor for growth rate r', default=10._rk, scale_factor=1/sec_per_year, minimum=0.0_rk)
   call self%get_parameter(n,          'n',      '-',           'exponent of max. consumption', default=2.0_rk/3.0_rk)
   call self%get_parameter(kappa,      'kappa',  'g^(lambda-1) m-3','pre-factor for carrying capacity', default=1.e11_rk, minimum=0.0_rk)
   call self%get_parameter(lambda,     'lambda', '-',           'exponent of background resource spectrum', default=2.8_rk-n)

   delta_logw = (log(w_max) - log(w_min)) / self%nclass       ! Log-weight distance between size classes [constant across spectrum]
   allocate(logw_bounds(self%nclass + 1))
   do iclass = 1, self%nclass + 1
      logw_bounds(iclass) = log(w_min) + delta_logw * (iclass - 1)
   end do

   allocate(self%w(self%nclass))
   self%w(:) = exp(0.5_rk * (logw_bounds(:self%nclass) + logw_bounds(2:)))  ! Weight of each size class at bin centres

   ! Compute size-specific parameters that do not vary in time. It is most efficient to compute these once, during initialization.
   allocate(self%r(self%nclass))
   allocate(self%K(self%nclass))
   self%r(:) = r0*self%w**(n-1)                   ! specific growth rate r
   self%K(:) = kappa*self%w**(-lambda+1) * (exp(logw_bounds(2:)) - exp(logw_bounds(:self%nclass)))  ! carrying capacity K (in Blanchard et al in units # g-1 - here, we premultiply with biomass and bin width!!)

   allocate(self%id_Nw(self%nclass))
   do iclass=1,self%nclass
      ! Postfix for size-class-specific variable names (an integer number)
      write (strindex,'(i0)') iclass

      ! Register state variable, store associated individual mass (used by predators, if any, to determine grazing preference).
      call self%register_state_variable(self%id_Nw(iclass), 'Nw'//trim(strindex), 'g m-3', 'biomass in size class '//trim(strindex), self%K(iclass), minimum=0.0_rk, maximum=self%K(iclass))
      call self%set_variable_property(self%id_Nw(iclass), 'particle_mass', self%w(iclass))
      call self%set_variable_property(self%id_Nw(iclass), 'min_particle_mass', exp(logw_bounds(iclass)))
      call self%set_variable_property(self%id_Nw(iclass), 'max_particle_mass', exp(logw_bounds(iclass + 1)))
      call self%add_to_aggregate_variable(total_mass, self%id_Nw(iclass))
   end do
   call self%register_state_variable(self%id_waste, 'waste', 'g m-3', 'waste')
   call self%add_to_aggregate_variable(total_mass,self%id_waste)

   end subroutine initialize
!EOC

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_resource_spectrum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer                         :: iclass
      real(rk),dimension(self%nclass) :: Nw, dNw

      _HORIZONTAL_LOOP_BEGIN_
         do iclass=1,self%nclass
            _GET_HORIZONTAL_(self%id_Nw(iclass), Nw(iclass))
         end do

         ! Semi-chemostatic growth for all size classes [not logistic!]
         dNw = self%r * (self%K - Nw)

         do iclass=1,self%nclass
            _SET_BOTTOM_ODE_(self%id_Nw(iclass), dNw(iclass))
         end do
         _SET_BOTTOM_ODE_(self%id_waste, -sum(dNw))
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

!-----------------------------------------------------------------------

end module mizer_resource_spectrum

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2016
!-----------------------------------------------------------------------
