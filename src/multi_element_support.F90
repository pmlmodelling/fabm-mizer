#include "fabm_driver.h"

module multi_element_support

   use fabm_types
   use fabm_particle
   use fabm_builtin_models

   implicit none

   type,extends(type_particle_model),public :: type_depth_integrated_sink
      type (type_model_id) :: id_target
      type (type_horizontal_diagnostic_variable_id) :: id_c_int
      type (type_horizontal_diagnostic_variable_id) :: id_n_int
      type (type_horizontal_diagnostic_variable_id) :: id_p_int
      type (type_horizontal_diagnostic_variable_id) :: id_s_int
      type (type_state_variable_id) :: id_c
      type (type_state_variable_id) :: id_n
      type (type_state_variable_id) :: id_p
      type (type_state_variable_id) :: id_s
      type (type_dependency_id) :: id_w
      type (type_horizontal_dependency_id) :: id_w_int
   contains
      procedure :: initialize => depth_integrated_sink_initialize
   end type type_depth_integrated_sink
   
   type,extends(type_particle_model),public :: type_depth_averaged_prey
      type (type_dependency_id) :: id_w
      type (type_dependency_id) :: id_c_pel
      type (type_dependency_id) :: id_n_pel
      type (type_dependency_id) :: id_p_pel
      type (type_dependency_id) :: id_s_pel
      type (type_horizontal_dependency_id) :: id_w_int
      type (type_horizontal_dependency_id) :: id_int_c_w
      type (type_horizontal_dependency_id) :: id_int_n_w
      type (type_horizontal_dependency_id) :: id_int_p_w
      type (type_horizontal_dependency_id) :: id_int_s_w
      type (type_horizontal_diagnostic_variable_id) :: id_c
      type (type_horizontal_diagnostic_variable_id) :: id_n
      type (type_horizontal_diagnostic_variable_id) :: id_p
      type (type_horizontal_diagnostic_variable_id) :: id_s
   contains
      procedure :: initialize => depth_averaged_prey_initialize
      procedure :: do_bottom => depth_averaged_prey_do_bottom
   end type type_depth_averaged_prey

   type,extends(type_particle_model),public :: type_relative_rate_distributor
      type (type_dependency_id) :: id_w
      type (type_horizontal_dependency_id) :: id_w_int
      type (type_horizontal_dependency_id) :: id_target_w_int
      type (type_horizontal_dependency_id) :: id_sms_int
      type (type_diagnostic_variable_id) :: id_r
      type (type_model_id) :: id_target
   contains
      procedure :: initialize => relative_rate_distributor_initialize
      procedure :: do => relative_rate_distributor_do
   end type

   type,extends(type_particle_model),public :: type_absolute_rate_distributor
      type (type_dependency_id) :: id_w
      type (type_horizontal_dependency_id) :: id_w_int
      type (type_horizontal_dependency_id) :: id_sms_int
      type (type_state_variable_id) :: id_target
      type (type_diagnostic_variable_id) :: id_r
   contains
      procedure :: initialize => absolute_rate_distributor_initialize
      procedure :: do => absolute_rate_distributor_do
   end type

   type,extends(type_particle_model),public :: type_depth_averaged_class
      type (type_horizontal_dependency_id) :: id_int_w
      type (type_horizontal_dependency_id) :: id_int_w2
      type (type_horizontal_dependency_id) :: id_int_c
      type (type_horizontal_diagnostic_variable_id) :: id_c
      real(rk) :: qnc, qpc
   contains
      procedure :: initialize => depth_averaged_class_initialize
      procedure :: do_bottom => depth_averaged_class_do_bottom
   end type type_depth_averaged_class

   type,extends(type_particle_model),public :: type_product
      type (type_dependency_id)          :: id_term1
      type (type_dependency_id)          :: id_term2
      type (type_diagnostic_variable_id) :: id_result
   contains
      procedure :: initialize => product_initialize
      procedure :: do => product_do
   end type type_product

   type,extends(type_particle_model),public :: type_square
      type (type_dependency_id)          :: id_source
      type (type_diagnostic_variable_id) :: id_result
      integer :: output = output_none
   contains
      procedure :: initialize => square_initialize
      procedure :: do => square_do
   end type

   type,extends(type_base_model),public :: type_pelagic_size_spectrum
      type (type_dependency_id), allocatable :: id_biomass(:)
      type (type_horizontal_dependency_id), allocatable :: id_biomass_int(:)
      type (type_horizontal_diagnostic_variable_id) :: id_slope
      type (type_horizontal_diagnostic_variable_id) :: id_offset
      real(rk), allocatable :: scale_factors(:, :), mass_grid(:)
   contains
      procedure :: initialize => pelagic_size_spectrum_initialize
      procedure :: do_bottom => pelagic_size_spectrum_do_bottom
   end type

contains

   subroutine depth_averaged_class_initialize(self, configunit)
!
! !INPUT PARAMETERS:
   class (type_depth_averaged_class), intent(inout), target :: self
   integer,                           intent(in)            :: configunit

   call self%register_dependency(self%id_int_w, 'int_w', '', 'depth-integral of weight')
   call self%register_dependency(self%id_int_w2, 'int_w2', '', 'depth-integral of squared weight')
   call self%register_dependency(self%id_int_c, 'int_c', '', 'depth-integral of carbon')
   call self%register_diagnostic_variable(self%id_c, 'c', 'mmol C/m^3', 'depth-averaged carbon', act_as_state_variable=.true., domain=domain_bottom, source=source_do_bottom, output=output_none)
   call self%add_to_aggregate_variable(standard_variables%total_carbon,     self%id_c)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen,   self%id_c, scale_factor=self%qnc)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_c, scale_factor=self%qpc)
   call copy_horizontal_fluxes(self, self%id_c, './int_c')

   end subroutine depth_averaged_class_initialize

   subroutine depth_averaged_class_do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
    class (type_depth_averaged_class), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_BOTTOM_

    real(rk) :: int_w, int_w2, int_c

    _HORIZONTAL_LOOP_BEGIN_
      _GET_HORIZONTAL_(self%id_int_w, int_w)
      _GET_HORIZONTAL_(self%id_int_w2, int_w2)
      _GET_HORIZONTAL_(self%id_int_c, int_c)
      _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c, int_c*int_w2/int_w/int_w)
    _HORIZONTAL_LOOP_END_
   end subroutine depth_averaged_class_do_bottom

   subroutine depth_integrated_sink_initialize(self, configunit)
!
! !INPUT PARAMETERS:
   class (type_depth_integrated_sink), intent(inout), target :: self
   integer,                            intent(in)            :: configunit

   character(len=10) :: composition

   call self%register_model_dependency(self%id_target, 'target')
   call self%register_dependency(self%id_w, 'w', '', 'weights')
   call self%register_dependency(self%id_w_int, 'w_int', '', 'depth integrated weights')
   
   call self%get_parameter(composition, 'composition', '', 'elemental composition')
   if (index(composition, 'c') /= 0)  call add_constituent('c', 'mmol C', 'carbon', self%id_c_int, self%id_c, standard_variables%total_carbon)
   if (index(composition, 'n') /= 0)  call add_constituent('n', 'mmol N', 'nitrogen', self%id_n_int, self%id_n, standard_variables%total_nitrogen)
   if (index(composition, 'p') /= 0)  call add_constituent('p', 'mmol P', 'phosphorus', self%id_p_int, self%id_p, standard_variables%total_phosphorus)
   if (index(composition, 's') /= 0)  call add_constituent('s', 'mmol Si', 'silicate', self%id_s_int, self%id_s, standard_variables%total_silicate)

   contains
   
      subroutine add_constituent(name, units, long_name, id_int, id_pelstate, standard_variable)
         character(len=*), intent(in) :: name, units, long_name
         type (type_horizontal_diagnostic_variable_id), intent(inout), target :: id_int
         type (type_state_variable_id),                 intent(inout), target :: id_pelstate
         type (type_bulk_standard_variable),            intent(in)            :: standard_variable

         class (type_absolute_rate_distributor), pointer :: rate_distributor

         call self%register_diagnostic_variable(id_int, name//'_int', units//'/m^2', 'depth-integrated '//long_name, act_as_state_variable=.true., domain=domain_bottom, source=source_none, output=output_none)
         call self%add_to_aggregate_variable(standard_variable, id_int)
         call self%register_state_dependency(id_pelstate, name, units//'/m^3', long_name)
         call self%request_coupling_to_model(id_pelstate, self%id_target, standard_variable)

         allocate(rate_distributor)
         call self%add_child(rate_distributor, name//'_source_distributor', configunit=-1)
         call rate_distributor%request_coupling(rate_distributor%id_target, '../'//name)
         call rate_distributor%request_coupling(rate_distributor%id_w, '../w')
         call rate_distributor%request_coupling(rate_distributor%id_w_int, '../w_int')
         call rate_distributor%request_coupling(rate_distributor%id_sms_int, '../'//name//'_int_sms_tot')
      end subroutine

   end subroutine depth_integrated_sink_initialize

   subroutine relative_rate_distributor_initialize(self, configunit)
      class (type_relative_rate_distributor),intent(inout),target :: self
      integer,                               intent(in)           :: configunit

      call self%register_model_dependency(self%id_target, 'target')
      call self%register_dependency(self%id_w, 'w', '', 'weights')
      call self%register_dependency(self%id_target_w_int, 'target_w_int', '', 'depth-integrated weights')
      call self%register_dependency(self%id_sms_int,'sms_int', '', 'depth-integrated sources-sinks')
      call self%register_diagnostic_variable(self%id_r, 'rate', 'd-1', 'relative change')
   end subroutine relative_rate_distributor_initialize

   subroutine relative_rate_distributor_do(self,_ARGUMENTS_DO_)
      class (type_relative_rate_distributor),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: target_w_int, sms_int, r, target, w
      integer :: istate

      _LOOP_BEGIN_
         ! First compute relative rate of change of depth-integrated target variable.
         _GET_HORIZONTAL_(self%id_target_w_int, target_w_int)
         _GET_HORIZONTAL_(self%id_sms_int, sms_int)
         if (target_w_int /= 0.0_rk) then
            r = sms_int/target_w_int
         else
            r = 0.0_rk
         end if

         _GET_(self%id_w, w)
         do istate=1,size(self%id_target%state)
            _GET_(self%id_target%state(istate), target)
            _SET_ODE_(self%id_target%state(istate), r * target * w)
         end do
         _SET_DIAGNOSTIC_(self%id_r, r * w * 86400)
      _LOOP_END_
   end subroutine relative_rate_distributor_do
   
   subroutine absolute_rate_distributor_initialize(self, configunit)
      class (type_absolute_rate_distributor), intent(inout),target :: self
      integer,                                intent(in)           :: configunit

      call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
      call self%register_dependency(self%id_w,'w','-','weights')
      call self%register_dependency(self%id_w_int,'w_int','','depth-integrated weights')
      call self%register_dependency(self%id_sms_int,'sms','','depth-integrated sources-sinks')
      call self%register_diagnostic_variable(self%id_r, 'rate', 'TARGET UNITS d-1', 'rate of change')
   end subroutine absolute_rate_distributor_initialize

   subroutine absolute_rate_distributor_do(self,_ARGUMENTS_DO_)
      class (type_absolute_rate_distributor), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_

      real(rk) :: sms_int, w, w_int

      _LOOP_BEGIN_
         ! First compute relative rate of change of depth-integrated target variable.
         _GET_HORIZONTAL_(self%id_sms_int, sms_int)
         _GET_HORIZONTAL_(self%id_w_int, w_int)
         _GET_(self%id_w, w)
         _SET_ODE_(self%id_target, sms_int * w/w_int)
         _SET_DIAGNOSTIC_(self%id_r, sms_int * w/w_int * 86400)
      _LOOP_END_
   end subroutine absolute_rate_distributor_do

   subroutine depth_averaged_prey_initialize(self, configunit)
!
! !INPUT PARAMETERS:
   class (type_depth_averaged_prey), intent(inout), target :: self
   integer,                          intent(in)            :: configunit

   class (type_relative_rate_distributor), pointer :: rate_distributor

   call self%register_model_dependency(name='source')
   call self%register_dependency(self%id_w, 'w', '', 'weight')
   call self%register_dependency(self%id_w_int, 'w_int', '', 'depth-integral of weight')

   call add_constituent('c', 'mmol C', 'carbon', self%id_c_pel, self%id_int_c_w, self%id_c, standard_variables%total_carbon)
   call add_constituent('n', 'mmol N', 'nitrogen', self%id_n_pel, self%id_int_n_w, self%id_n, standard_variables%total_nitrogen)
   call add_constituent('p', 'mmol P', 'phosphorus', self%id_p_pel, self%id_int_p_w, self%id_p, standard_variables%total_phosphorus)
   call add_constituent('s', 'mmol Si', 'silicate', self%id_s_pel, self%id_int_s_w, self%id_s, standard_variables%total_silicate)

   allocate(rate_distributor)
   call self%add_child(rate_distributor, 'loss_distributor', configunit=-1)
   call rate_distributor%request_coupling('w', '../w')
   call rate_distributor%request_coupling('target_w_int', '../int_c_w')
   call rate_distributor%request_coupling('sms_int', '../c_sms_tot')
   call rate_distributor%couplings%set_string('target', '../source')

   contains

      subroutine add_constituent(name, units, long_name, id_pel, id_int_c_w, id_ave, standard_variable)
         character(len=*),                              intent(in)            :: name, units, long_name
         type (type_dependency_id),                     intent(inout), target :: id_pel
         type (type_horizontal_dependency_id),          intent(inout), target :: id_int_c_w
         type (type_horizontal_diagnostic_variable_id), intent(inout), target :: id_ave
         type (type_bulk_standard_variable),            intent(in)            :: standard_variable

         class (type_product),        pointer :: product
         class (type_depth_integral), pointer :: depth_integral

         call self%register_dependency(id_pel, name//'_pel', units//'/m^3', long_name)
         call self%request_coupling_to_model(id_pel, 'source', standard_variable)

         allocate(product)
         call self%add_child(product, name//'_w_calculator', configunit=-1)
         call product%request_coupling('term1', '../w')
         call product%request_coupling('term2', '../'//name//'_pel')

         allocate(depth_integral)
         call self%add_child(depth_integral, name//'_w_integrator', configunit=-1)
         call depth_integral%request_coupling('source', '../'//name//'_w_calculator/result')
         depth_integral%id_output%link%target%output = output_none

         call self%register_dependency(id_int_c_w, 'int_'//name//'_w', '', 'depth-integral of '//long_name//' * weight')
         call self%request_coupling(id_int_c_w, name//'_w_integrator/result')

         call self%register_diagnostic_variable(id_ave, name, units//'/m^3', 'depth-averaged '//long_name, act_as_state_variable=.true., domain=domain_bottom, source=source_do_bottom)
         call self%add_to_aggregate_variable(standard_variable, id_ave)

      end subroutine

   end subroutine depth_averaged_prey_initialize

   subroutine depth_averaged_prey_do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
    class (type_depth_averaged_prey), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_BOTTOM_

    real(rk) :: c_w_int, n_w_int, p_w_int, s_w_int, w_int

    _HORIZONTAL_LOOP_BEGIN_
        _GET_HORIZONTAL_(self%id_w_int, w_int)
        _GET_HORIZONTAL_(self%id_int_c_w, c_w_int)
        _GET_HORIZONTAL_(self%id_int_n_w, n_w_int)
        _GET_HORIZONTAL_(self%id_int_p_w, p_w_int)
        _GET_HORIZONTAL_(self%id_int_s_w, s_w_int)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c, c_w_int/w_int)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_n, n_w_int/w_int)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_p, p_w_int/w_int)
        _SET_HORIZONTAL_DIAGNOSTIC_(self%id_s, s_w_int/w_int)
    _HORIZONTAL_LOOP_END_
   end subroutine depth_averaged_prey_do_bottom

   subroutine square_initialize(self, configunit)
!
! !INPUT PARAMETERS:
   class (type_square), intent(inout), target :: self
   integer,             intent(in)            :: configunit

   call self%register_dependency(self%id_source, 'source', '', 'source')
   call self%register_diagnostic_variable(self%id_result, 'result', '', 'result', output=self%output)

   end subroutine square_initialize

   subroutine square_do(self, _ARGUMENTS_DO_)
    class (type_square), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) :: value

    _LOOP_BEGIN_
        _GET_(self%id_source, value)
        _SET_DIAGNOSTIC_(self%id_result, value*value)
    _LOOP_END_
   end subroutine square_do

   subroutine product_initialize(self, configunit)
!
! !INPUT PARAMETERS:
   class (type_product), intent(inout), target :: self
   integer,              intent(in)            :: configunit

   call self%register_dependency(self%id_term1, 'term1', '', 'term1')
   call self%register_dependency(self%id_term2, 'term2', '', 'term2')
   call self%register_diagnostic_variable(self%id_result, 'result', '', 'term1 * term2', output=output_none)

   end subroutine product_initialize

   subroutine product_do(self, _ARGUMENTS_DO_)
    class (type_product), intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    real(rk) :: term1, term2

    _LOOP_BEGIN_
        _GET_(self%id_term1, term1)
        _GET_(self%id_term2, term2)
        _SET_DIAGNOSTIC_(self%id_result, term1 * term2)
    _LOOP_END_
   end subroutine product_do

   subroutine pelagic_size_spectrum_initialize(self, configunit)
      class (type_pelagic_size_spectrum), intent(inout), target :: self
      integer,                            intent(in)            :: configunit

      integer  :: ngrid, nsource, i, isource
      real(rk) :: w_min, w_max, dlog_w
      character(len=10)  :: strindex
      class (type_depth_integral), pointer :: depth_integral

      call self%get_parameter(ngrid, 'ngrid', '', 'grid size', default=100)
      call self%get_parameter(nsource, 'nsource', '', 'number of source pools')
      call self%get_parameter(w_min, 'w_min', '', 'minimum mass')
      call self%get_parameter(w_max, 'w_max', '', 'maximum mass')

      allocate(self%mass_grid(ngrid))
      allocate(self%scale_factors(ngrid, nsource))
      allocate(self%id_biomass(nsource))
      allocate(self%id_biomass_int(nsource))
      dlog_w = (log(w_max) - log(w_min)) / ngrid
      do i = 1, ngrid
         self%mass_grid(i) =  log(w_min) + (i - 0.5_rk) * dlog_w
      end do
      self%scale_factors = 0

      do isource = 1, nsource
        write (strindex,'(i0)') isource
        call self%register_dependency(self%id_biomass(isource), 'source'//trim(strindex), '<SOURCE UNITS>', 'source '//trim(strindex))
        call self%get_parameter(w_min, 'w_source'//trim(strindex)//'_min', 'g', 'minimum mass of source '//trim(strindex))
        call self%get_parameter(w_max, 'w_source'//trim(strindex)//'_max', 'g', 'maximum mass of source '//trim(strindex))
        do i = 1, ngrid
           if (self%mass_grid(i) >= log(w_min) .and. self%mass_grid(i) <= log(w_max)) &
              self%scale_factors(i, isource) = (min(log(w_max), self%mass_grid(i) + dlog_w/2) - max(log(w_min), self%mass_grid(i) - dlog_w/2))/(log(w_max) - log(w_min))
        end do
        allocate(depth_integral)
        call self%add_child(depth_integral, 'source'//trim(strindex)//'_integrator', configunit=-1)
        call depth_integral%request_coupling('source', '../source'//trim(strindex))
        call self%register_dependency(self%id_biomass_int(isource), 'source'//trim(strindex)//'_int', '<SOURCE UNITS> m', 'depth integrated source '//trim(strindex))
        call self%request_coupling(self%id_biomass_int(isource), './source'//trim(strindex)//'_integrator/result')
      end do
      do i = 1, size(self%mass_grid)
         self%scale_factors(i, :) = self%scale_factors(i, :) / (exp(self%mass_grid(i) + dlog_w/2) - exp(self%mass_grid(i) - dlog_w/2))
      end do
      call self%register_diagnostic_variable(self%id_slope, 'slope', '-', 'slope of size spectrum', source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_offset, 'offset', 'ln <SOURCE UNITS> m', 'offset of size spectrum', source=source_do_bottom)
   end subroutine pelagic_size_spectrum_initialize

   subroutine pelagic_size_spectrum_do_bottom(self, _ARGUMENTS_DO_BOTTOM_)
      class (type_pelagic_size_spectrum), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer  :: isource
      real(rk) :: spectrum(size(self%mass_grid)), bm_int, cov, var, slope, offset
      real(rk), allocatable :: x(:),  y(:)
      integer  :: i, n

      _HORIZONTAL_LOOP_BEGIN_
         spectrum = 0
         do isource = 1, size(self%scale_factors, 2)
            _GET_HORIZONTAL_(self%id_biomass_int(isource), bm_int)
            spectrum = spectrum + self%scale_factors(:, isource) * bm_int
         end do
         n = 0
         do i = 1, size(self%mass_grid)
            if (spectrum(i) > 0) n = n + 1
         end do
         allocate(x(n))
         allocate(y(n))
         n = 0
         do i = 1, size(self%mass_grid)
            if (spectrum(i) > 0) then
               n = n + 1
               x(n) = self%mass_grid(i)
               y(n) = log(spectrum(i))
            end if
         end do
         cov = sum(x*y)/n - sum(x)*sum(y)/n/n
         var = sum(x**2)/n - (sum(x)/n)**2
         slope = cov/var
         offset = sum(y)/n - sum(slope*x)/n
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_slope, slope)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_offset, offset)
         deallocate(x)
         deallocate(y)
      _HORIZONTAL_LOOP_END_

   end subroutine pelagic_size_spectrum_do_bottom

end module multi_element_support