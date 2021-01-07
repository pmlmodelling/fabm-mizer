#include "fabm_driver.h"

module mizer_resource_spectrum

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

   use fabm_types
   use size_structured_base

   implicit none

   private

   type, extends(type_size_structured_base),public :: type_resource_spectrum
      ! Variable identifiers
      type (type_bottom_state_variable_id)             :: id_waste
      real(rk), allocatable :: r(:)
      real(rk), allocatable :: K(:)
   contains
      procedure :: initialize
      procedure :: do_bottom
   end type type_resource_spectrum

contains

   subroutine initialize(self,configunit)
      class (type_resource_spectrum), intent(inout), target :: self
      integer,                        intent(in )           :: configunit

      real(rk)             :: kappa, lambda, r0, w_min, w_max, n
      real(rk),parameter   :: sec_per_year = 86400 * 365.2425_rk

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

      call self%configure(self%nclass, w_min, w_max)
      call self%request_coupling(self%id_depth, standard_variables%bottom_depth)

      ! Compute size-specific parameters that do not vary in time. It is most efficient to compute these once, during initialization.
      allocate(self%r(self%nclass))
      allocate(self%K(self%nclass))
      self%r(:) = r0*self%w**(n-1)                   ! specific growth rate r
      self%K(:) = kappa*self%w**(-lambda+1) * (self%w_bounds(2:) - self%w_bounds(:self%nclass))  ! carrying capacity K (in Blanchard et al in units # g-1 - here, we premultiply with biomass and bin width!!)

      call self%register_state_variable(self%id_waste, 'waste', 'g m-3', 'waste')
      call self%add_to_aggregate_variable(total_mass, self%id_waste)
   end subroutine initialize

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_resource_spectrum),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer                          :: iclass
      real(rk), dimension(self%nclass) :: Nw, dNw
      real(rk)                         :: depth

      _HORIZONTAL_LOOP_BEGIN_
         do iclass=1,self%nclass
            _GET_HORIZONTAL_(self%id_Nw(iclass), Nw(iclass))
         end do
         _GET_HORIZONTAL_(self%id_depth, depth)

         ! Semi-chemostatic growth for all size classes [not logistic!]
         ! Units are g m-2 s-1
         dNw = self%r * (self%K * depth - Nw)

         do iclass=1,self%nclass
            _SET_BOTTOM_ODE_(self%id_Nw(iclass), dNw(iclass))
         end do
         _SET_BOTTOM_ODE_(self%id_waste, -sum(dNw))
      _HORIZONTAL_LOOP_END_
   end subroutine do_bottom

end module mizer_resource_spectrum

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2016
!-----------------------------------------------------------------------
