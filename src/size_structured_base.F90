#include "fabm_driver.h"

module size_structured_base

   use fabm_types

   implicit none

   private

   public :: type_size_structured_base, total_mass

   type, extends(type_base_model) :: type_map_density_to_concentration
      type (type_horizontal_dependency_id),          allocatable :: id_density(:)
      type (type_horizontal_dependency_id)                       :: id_depth
      type (type_horizontal_diagnostic_variable_id), allocatable :: id_concentration(:)
   contains
      procedure :: do_bottom => map_density_to_concentration_do_bottom
   end type

   type, extends(type_base_model) :: type_apply_loss
      type (type_horizontal_dependency_id), allocatable :: id_bottom_flux(:)
      type (type_bottom_state_variable_id), allocatable :: id_density(:)
      type (type_horizontal_diagnostic_variable_id), allocatable :: id_mortality(:)
   contains
      procedure :: do_bottom => apply_loss_do_bottom
   end type

   type, extends(type_base_model) :: type_size_structured_base
      type (type_bottom_state_variable_id), allocatable :: id_Nw(:) ! Total mass (g m-2) per size class (sum over all individuals)
      type (type_horizontal_dependency_id)              :: id_depth              ! Predator-prey interaction depth
      integer :: nclass

      ! Size class characteristics
      real(rk), allocatable :: logw(:)         ! log individual mass
      real(rk), allocatable :: w(:)            ! individual mass
      real(rk), allocatable :: w_bounds(:)     ! individualmass at interfaces
      real(rk), allocatable :: delta_w(:)      ! individual mass difference between centers of current and next size class
   contains
      procedure :: configure
   end type

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable), parameter :: total_mass = type_bulk_standard_variable(name='total_mass',units='g',aggregate_variable=.true.,conserved=.true.)

contains

   subroutine configure(self, nclass, w_min, w_max)
      class (type_size_structured_base), intent(inout) :: self
      integer,                           intent(in) :: nclass
      real(rk),                          intent(in) :: w_min, w_max

      integer            :: iclass
      real(rk)           :: delta_logw
      character(len=10)  :: strindex
      class(type_map_density_to_concentration),       pointer :: scale_state
      class(type_apply_loss),                         pointer :: loss_copier

      self%nclass = nclass

      ! Determine size classes (log-spaced between size at birth and infinite size)
      allocate(self%logw(self%nclass))
      allocate(self%w(self%nclass))
      allocate(self%w_bounds(self%nclass + 1))
      allocate(self%delta_w(self%nclass))
      delta_logw = (log(w_max) - log(w_min)) / self%nclass                    ! Log-mass distance between size classes [constant across spectrum]
      self%w_bounds(1) = w_min
      do iclass = 1, self%nclass
         self%logw(iclass) = log(w_min) + delta_logw * (iclass - 0.5_rk)      ! Log mass of at the center of the size class
         self%w_bounds(iclass + 1) = exp(log(w_min) + delta_logw * iclass)    ! Mass of right bound of the size class
      end do
      self%w(:) = exp(self%logw)                                              ! Mass of each size class
      self%delta_w(:) = exp(self%logw + delta_logw) - self%w                  ! Mass difference between consecutive size classes

      ! Depth over which predator-prey interact
      ! If cannibalism is active, the fish per area will be divided by this depth to compute the prey concentration
      call self%register_dependency(self%id_depth, 'interaction_depth', 'm', 'predator-prey interaction depth')

      ! Set up submodel that calculates prey concentrations (g m-3) from internal biomass density (g m-2)
      ! This is used only for our own size classes - not external prey, which are assumed to be already expressed in g m-3
      allocate(scale_state)
      call self%add_child(scale_state, 'biomass_as_prey', configunit=-1)
      allocate(scale_state%id_density(self%nclass))
      allocate(scale_state%id_concentration(self%nclass))
      call scale_state%register_dependency(scale_state%id_depth, 'interaction_depth', 'm', 'predator-prey interaction depth')
      call scale_state%request_coupling(scale_state%id_depth, '../interaction_depth')

      ! Create submodel that applies the bottom flux of our own size class concentrations (g m-3) to the
      ! size class densities (g m-2) and that computes a specific loss rate (d-1) due to predation.
      allocate(loss_copier)
      call scale_state%add_child(loss_copier, 'flux_copier', configunit=-1)
      allocate(loss_copier%id_density(self%nclass))
      allocate(loss_copier%id_bottom_flux(self%nclass))
      allocate(loss_copier%id_mortality(self%nclass))
   
      allocate(self%id_Nw(self%nclass))
      do iclass = 1, self%nclass
         ! Postfix for size-class-specific variable names (an integer number)
         write (strindex,'(i0)') iclass

         ! Register state variable for biomass density (g m-2) for this size class
         ! Add attribute for individual mass (used by predators, if any, to determine grazing preference).
         call self%register_state_variable(self%id_Nw(iclass), 'Nw'//trim(strindex), 'g m-2', 'biomass in size class '//trim(strindex), 1.0_rk, minimum=0.0_rk)
         call self%set_variable_property(self%id_Nw(iclass), 'particle_mass',self%w(iclass))
         call self%set_variable_property(self%id_Nw(iclass), 'min_particle_mass', self%w_bounds(iclass))
         call self%set_variable_property(self%id_Nw(iclass), 'max_particle_mass', self%w_bounds(iclass + 1))

         ! Register this size class' contribution to total mass in the system (for mass conservation checks)
         call self%add_to_aggregate_variable(total_mass,self%id_Nw(iclass))

         ! Convert the biomass density (g m-2) into concentration (g m-3) for use by predators.
         call scale_state%register_dependency(scale_state%id_density(iclass), 'bm'//trim(strindex), 'g m-2', 'biomass '//trim(strindex))
         call scale_state%request_coupling(scale_state%id_density(iclass), '../Nw'//trim(strindex))
         call scale_state%register_diagnostic_variable(scale_state%id_concentration(iclass), 'c'//trim(strindex), 'g m-3', 'biomass concentration', act_as_state_variable=.true., output=output_none, domain=domain_bottom, source=source_do_bottom)
         call self%set_variable_property(scale_state%id_concentration(iclass), 'particle_mass', self%w(iclass))

         ! Apply bottom flux of biomass concentration (g m-2 s-1) to original biomass density
         call loss_copier%register_state_dependency(loss_copier%id_density(iclass), 'bm'//trim(strindex), 'g m-2', 'biomass '//trim(strindex))
         call loss_copier%request_coupling(loss_copier%id_density(iclass), '../../Nw'//trim(strindex))
         call loss_copier%register_dependency(loss_copier%id_bottom_flux(iclass), 'bottom_flux'//trim(strindex), 'g m-2 s-1', 'biomass '//trim(strindex))
         call loss_copier%request_coupling(loss_copier%id_bottom_flux(iclass), '../c'//trim(strindex)//'_sms_tot')
         call loss_copier%register_diagnostic_variable(loss_copier%id_mortality(iclass), 'loss'//trim(strindex), 'd-1', 'mortality '//trim(strindex), source=source_do_bottom)
      end do
   end subroutine

   subroutine map_density_to_concentration_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_map_density_to_concentration), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iclass
      real(rk) :: depth, density

      _HORIZONTAL_LOOP_BEGIN_
         _GET_HORIZONTAL_(self%id_depth, depth)
         do iclass = 1, size(self%id_density)
            _GET_HORIZONTAL_(self%id_density(iclass), density)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_concentration(iclass), density / depth)
         end do
      _HORIZONTAL_LOOP_END_
   end subroutine

   subroutine apply_loss_do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_apply_loss), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iclass
      real(rk) :: bottom_flux, density

      _HORIZONTAL_LOOP_BEGIN_
         do iclass = 1, size(self%id_density)
            _GET_HORIZONTAL_(self%id_bottom_flux(iclass), bottom_flux)
            _GET_HORIZONTAL_(self%id_density(iclass), density)
            _SET_BOTTOM_ODE_(self%id_density(iclass), bottom_flux)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_mortality(iclass), -86400 * bottom_flux / (density + tiny(density)))
         end do
      _HORIZONTAL_LOOP_END_
   end subroutine

end module