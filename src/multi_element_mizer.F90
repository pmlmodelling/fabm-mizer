#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mizer
!
! !INTERFACE:
module mizer_multi_element_population
!
! !DESCRIPTION:
! Size-structured population for 2D fauna (e.g., depth-integrated fish or benthic fauna), based on:
!
! Blanchard, J. L., Andersen, K. H., Scott, F., Hintzen, N. T., Piet, G., & Jennings, S. (2014)
! Evaluating targets and trade-offs among fisheries and conservation objectives using a multispecies size spectrum model
! Journal of Applied Ecology, 51(3), 612-622
! doi:10.1111/1365-2664.12238
!
! Notation mostly follows the mizer R package.
!
! The biomass/abundance per species is represented by the biomass in the respective size bin (g),
! as in Sheldon's original postulate. This quantity can be converted to a biomass density (g g-1) by dividing
! by the bin width (in g), and to abundance density (# g-1) by further dividing by the individual biomass (g)
! at the centre of the bin. Note that each division roughly corresponds to a decrease of 1 in spectrum slope.
!
! !USES:
   use fabm_types
   use fabm_particle

   use multi_element_support

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_particle_model),public :: type_multi_element_population
      ! Variable identifiers
      type (type_bottom_state_variable_id),         allocatable :: id_c(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_c(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_n(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_p(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_s(:)
      type (type_bottom_state_variable_id)                      :: id_waste_c
      type (type_bottom_state_variable_id)                      :: id_waste_n
      type (type_bottom_state_variable_id)                      :: id_waste_p
      type (type_bottom_state_variable_id)                      :: id_waste_s
      type (type_bottom_state_variable_id)                      :: id_o2
      type (type_bottom_state_variable_id)                      :: id_dic
      type (type_bottom_state_variable_id)                      :: id_din
      type (type_bottom_state_variable_id)                      :: id_dip
      type (type_bottom_state_variable_id)                      :: id_discard_c
      type (type_bottom_state_variable_id)                      :: id_discard_n
      type (type_bottom_state_variable_id)                      :: id_discard_p
      type (type_bottom_state_variable_id)                      :: id_landings           ! State variable that will serve as sink for all landed biomass
      type (type_horizontal_diagnostic_variable_id)             :: id_total_reproduction ! Total reproduction
      type (type_horizontal_diagnostic_variable_id)             :: id_R_p                ! Density-independent recruitment
      type (type_horizontal_diagnostic_variable_id)             :: id_R                  ! Density-dependent recruitment
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_reproduction(:)    ! Reproduction per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_f(:)               ! Functional response per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_g(:)               ! Specific growth rate per size class
      type (type_dependency_id),                    allocatable :: id_pelprey_c(:)

      type (type_horizontal_dependency_id) :: id_slope, id_offset

      type (type_horizontal_dependency_id)                      :: id_T_w_int
      type (type_horizontal_dependency_id)                      :: id_w_int
      type (type_horizontal_diagnostic_variable_id)             :: id_T_ave
      type (type_horizontal_diagnostic_variable_id)             :: id_c_tot
      type (type_horizontal_diagnostic_variable_id)             :: id_c_size1
      type (type_horizontal_diagnostic_variable_id)             :: id_c_size2
      type (type_horizontal_diagnostic_variable_id)             :: id_c_size3
      type (type_horizontal_diagnostic_variable_id)             :: id_c_lfi
      real(rk)                                                  :: w_threshold
      real(rk)                                                  :: w_threshold2
      real(rk)                                                  :: w_threshold3
      real(rk)                                                  :: lfi_w_threshold


      ! Number of size classes and prey
      integer :: nclass
      integer :: nprey

      ! Size class characteristics
      real(rk),allocatable :: logw(:)     ! log mass
      real(rk),allocatable :: w(:)        ! mass
      real(rk),allocatable :: delta_w(:)  ! mass difference between consecutive size classes

      ! Size-class-independent parameters
      real(rk) :: w_min       ! egg mass
      real(rk) :: alpha       ! assimilation efficiency
      real(rk) :: beta        ! preferred predator:prey mass ratio
      real(rk) :: sigma       ! s.d. of lognormal prey size selection function
      real(rk) :: erepro      ! efficiency of egg production
      real(rk) :: xi          ! fraction of mass consisting of lipid reserve (which can fuel maintenance)
      integer  :: SRR         ! type of stock-recruitment relationship (SRR)
      real(rk) :: recruitment ! constant recruitment flux (SSR=0)
      real(rk) :: R_max       ! maximum recruitment flux (SSR=2)
      real(rk) :: R_relax     ! rate of relaxation towards expected egg density - derived from prey spectrum (SSR=3)
      real(rk) :: qnc         ! nitrogen:carbon ratio (mol:mol, constant)
      real(rk) :: qpc         ! phosphorus:carbon ratio (mol:mol, constant)
      real(rk) :: alpha_eg    ! fraction of food that is egested

      integer  :: T_dependence ! Type of temperature dependence (0: none, 1: Arrhenius)
      real(rk) :: c1           ! Reference constant in Arrhenius equation = E_a/k/(T_ref+Kelvin)
      real(rk) :: E_a          ! Activation energy (eV)
 !     real(rk) :: resp_o2C      ! oxygen produced per carbon respired (mol:mol,constant)

      logical :: feedback

      ! Size-class-dependent parameters that will be precomputed during initialization
      real(rk), allocatable :: V(:)            ! volumetric search rate (Eq M2)
      real(rk), allocatable :: I_max(:)        ! maximum ingestion rate (Eq M4)
      real(rk), allocatable :: std_metab(:)    ! standard metabolism (k*w^p in Eq M7)
      real(rk), allocatable :: mu_b(:)         ! background mortality (temperature dependent)
      real(rk), allocatable :: mu_s(:)         ! senescence mortality (temperature independent)
      real(rk), allocatable :: F(:)            ! fishing mortality
      real(rk), allocatable :: psi(:)          ! allocation to reproduction
      real(rk), allocatable :: phi(:,:)        ! prey preference

      real(rk) :: c0_proxy_width
   contains
      procedure :: initialize
      procedure :: do_bottom
      procedure :: after_coupling
   end type type_multi_element_population

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable),parameter :: total_mass = type_bulk_standard_variable(name='total_mass',units='g',aggregate_variable=.true.,conserved=.true.)

   real(rk) :: Kelvin = 273.15_rk      ! offset of Celsius temperature scale (K)
   real(rk) :: Boltzmann = 8.62e-5_rk  ! Boltzmann constant
   real(rk) :: g_per_mmol_carbon = 0.12_rk ! assume 10% of wet mass is carbon. 1 mmol carbon = 0.012 g, so that is equivalent to 0.12 g wet mass
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
   class (type_multi_element_population), intent(inout),target :: self
   integer,                               intent(in )          :: configunit
!
! !LOCAL VARIABLES:
   integer            :: iclass, iprey
   logical            :: cannibalism, biomass_has_prey_unit
   real(rk)           :: delta_logw
   real(rk)           :: k_vb,n,q,p,w_mat,w_inf,gamma,h,ks,f0,z0
   real(rk)           :: z0pre,z0exp,w_s,z_s,z_spre
   real(rk)           :: kappa,lambda
   real(rk)           :: T_ref
   real(rk)           :: S1,S2,F,w_minF,F_a,F_b
   integer            :: z0_type
   integer            :: fishing_type
   real(rk)           :: w_prey_min, w_prey_max
   real(rk)           :: c_ini
   character(len=10)  :: strindex, strindex2
   real(rk),parameter :: pi = 4*atan(1.0_rk)
   real(rk),parameter :: sec_per_year = 86400*365.2425_rk
   class (type_weighted_sum), pointer :: total_pelprey_calculator
   class (type_depth_averaged_prey),  pointer :: depth_averaged_prey
   class (type_depth_averaged_class), pointer :: depth_averaged_class
   class (type_square),               pointer :: square
   class (type_depth_integral),       pointer :: depth_integral
   class (type_product),              pointer :: product
   class (type_pelagic_size_spectrum),pointer :: pelagic_size_spectrum
   logical, parameter :: report_statistics = .false.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read parameters
   ! All rate coefficients are converted from d-1 to s-1 ensure source terms are directly compatible with FABM.
   ! Default values taken from Table S5
   call self%get_parameter(c_ini,'c_ini', 'mmol C/m2',     'initial density per size class', default=0.0_rk)
   call self%get_parameter(self%nclass,'nclass', '',     'number of size classes', default=100)
   call self%get_parameter(self%nprey, 'nprey',  '',     'number of prey')
   call self%get_parameter(self%alpha, 'alpha',  '-',    'assimilation efficiency',            default=0.6_rk,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%alpha_eg, 'alpha_eg',  '-',    'fraction of food egested', default=1-self%alpha,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%erepro,'erepro', '-',    'reproductive efficiency',            default=1.0_rk,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%w_min, 'w_min',  'g',    'egg mass',                           default=0.001_rk, minimum=0.0_rk)
   call self%get_parameter(n,          'n',      '-',    'exponent of max. consumption',       default=2.0_rk/3.0_rk)
   call self%get_parameter(q,          'q',      '-',    'exponent of search volume',          default=0.8_rk)
   call self%get_parameter(p,          'p',      '-',    'exponent of standard metabolism',    default=0.7_rk)
   call self%get_parameter(z0_type,    'z0_type','',     'type of background mortality (0: constant, 1: allometric function of size)', default=0)
   call self%get_parameter(z0pre,      'z0pre',  'yr-1', 'pre-factor for background mortality (= mortality at 1 g)',default=0.6_rk,   minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(z0exp,      'z0exp',  '-',    'exponent of background mortality',   default=n-1)
   call self%get_parameter(w_s,        'w_s',    'g',    'start mass for senescence mortality',default=0._rk, minimum=0.0_rk)
   call self%get_parameter(z_s,        'z_s',    '-',    'exponent for senescence mortality',default=0.3_rk, minimum=0.0_rk)
   call self%get_parameter(z_spre,     'z_spre', 'yr-1', 'pre-factor for senescence mortality (= mortality at w_s g)',default=0.2_rk, minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(w_mat,      'w_mat',  'g',    'maturation mass', default=0.0_rk, minimum=0.0_rk)
   call self%get_parameter(w_inf,      'w_inf',  'g',    'asymptotic mass', default=1e3_rk, minimum=0.0_rk)
   call self%get_parameter(self%w_threshold,'w_threshold','g','mass separating size1 and large individuals (for diagnostic output only)', default=10._rk, minimum=0.0_rk)
   call self%get_parameter(self%w_threshold2,'w_threshold2','g','mass separating size2 and large individuals (for diagnostic output only)', default=80._rk, minimum=0.0_rk)
   call self%get_parameter(self%w_threshold3,'w_threshold3','g','mass separating size3 and large individuals (for diagnostic output only)', default=1000._rk, minimum=0.0_rk)
   call self%get_parameter(self%lfi_w_threshold,'lfi_w_threshold','g','mass for lfi index separating small and large individuals (for diagnostic output only)', default=582._rk, minimum=0.0_rk)
   call self%get_parameter(self%beta,  'beta',   '-',    'preferred predator:prey mass ratio', default=100.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%sigma, 'sigma',  '-',    'width of prey size preference (sd in ln mass units)', default=1.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%xi,    'xi',     '-',    'fraction of mass consisting of lipid reserve', default=0.1_rk, minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(k_vb,       'k_vb',   'yr-1', 'von Bertalanffy growth rate', minimum=0.0_rk, default=0.0_rk, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(lambda,     'lambda', '-',        'exponent of background resource spectrum', default=2+q-n)
   call self%get_parameter(kappa,      'kappa',  'g^(lambda-1)','carrying capacity of background resource spectrum', default=1e11_rk, minimum=0.0_rk)
   call self%get_parameter(f0,         'f0',     '-',        'background feeding level', default=0.6_rk, minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%SRR,   'SRR',    '',     'stock-recruitment relationship (0: constant recruitment, 1: density-independent recruitment, 2: Beverton-Holt)')
   select case (self%SRR)
   case (0)
      call self%get_parameter(self%recruitment, 'recruitment', '# m-2 yr-1', 'constant recruitment flux', minimum=0.0_rk, default=kappa*self%w_min**(-lambda)*sec_per_year, scale_factor=1._rk/sec_per_year)
   case (2)
      call self%get_parameter(self%R_max, 'R_max','# m-2 yr-1','maximum recruitment flux', minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   case (3)
      call self%get_parameter(self%R_relax, 'R_relax', 'd-1', 'relaxation towards expected egg density', minimum=0.0_rk, default=1.0_rk, scale_factor=1._rk/86400)
   end select

   call self%get_parameter(cannibalism,'cannibalism','',         'whether to enable intraspecific predation', default=.true.)
   if (cannibalism) call self%get_parameter(biomass_has_prey_unit, 'biomass_has_prey_unit', '', 'biomass has the same unit as prey', default=.true.)
   call self%get_parameter(self%qnc,   'qnc',        'mol mol-1','nitrogen to carbon ratio', default=16.0_rk/106.0_rk)
   call self%get_parameter(self%qpc,   'qpc',        'mol mol-1','phosphorus to carbon ratio', default=1.0_rk/106.0_rk)
!      call self%get_parameter(self%resp_o2C,   'resp_o2C',        'mol mol-1','oxygen consumed per carbon respired', default=127_rk/106.0_rk)
   call self%get_parameter(self%feedback, 'feedback', '', 'feedback from fish to ecosystem (prey, waste, O2, CO2)', default=.true.)

   call self%get_parameter(self%T_dependence, 'T_dependence', '', 'temperature dependence (0: none, 1: Arrhenius)', default=0)
   select case (self%T_dependence)
   case (1)
      call self%get_parameter(self%E_a, 'E_a',   'eV',              'activation energy',     default=0.63_rk)
      call self%get_parameter(T_ref,    'T_ref', 'degrees_Celsius', 'reference temperature', default=13._rk)
      self%c1 = self%E_a/Boltzmann/(T_ref+Kelvin)
      call self%register_dependency(self%id_T_w_int, 'T_w_int', 'degrees_Celsius', 'average experienced temperature')
      call self%register_dependency(self%id_w_int, 'w_int', '', 'depth-integrated weight')

      allocate(product)
      call self%add_child(product, 'T_w_calculator', configunit=-1)
      call product%request_coupling(product%id_term1, standard_variables%temperature)
      call product%request_coupling(product%id_term2, '../total_pelprey_calculator/result')

      allocate(depth_integral)
      call self%add_child(depth_integral, 'T_w_integrator', configunit=-1)
      call depth_integral%request_coupling('source', '../T_w_calculator/result')
      depth_integral%id_output%link%target%output = output_none

      call self%request_coupling(self%id_w_int, './w_integrator/result')
      call self%request_coupling(self%id_T_w_int, './T_w_integrator/result')
      call self%register_diagnostic_variable(self%id_T_ave, 'T_ave' ,'degrees_Celsius', 'average experienced temperature', source=source_do_bottom)
   end select

   ! Pre-factor for maximum ingestion rate derived from von Bertalanffy growth rate and background feeding level as described in appendix B
   h = 3*k_vb*w_inf**(1.0_rk/3.0_rk)/self%alpha/f0
   call self%get_parameter(h, 'h', 'yr-1 g^(-n)', 'pre-factor for maximum ingestion rate', minimum=0.0_rk, default=h*sec_per_year, scale_factor=1._rk/sec_per_year)

   ! Pre-factor for volumetric search rate (Eq M2)
   gamma = f0*h*self%beta**(2-lambda)/((1-f0)*sqrt(2*pi)*kappa*self%sigma)
   !gamma = f0*h*self%beta**(2-lambda)/((1-f0)*sqrt(2*pi)*kappa*self%sigma*exp((lambda-2)**2 * self%sigma**2 / 2)) ! add exp term taken from actual R code
   call self%get_parameter(gamma, 'gamma', 'PV yr-1 g^(-q)', 'pre-factor for volumetric search rate', minimum=0.0_rk, default=gamma*sec_per_year, scale_factor=1._rk/sec_per_year)

   ! Allow user override of standard metabolism pre-factor (e.g., Blanchard community size spectrum model has ks=0)
   call self%get_parameter(ks, 'ks', 'yr-1 g^(-p)', 'pre-factor for standard metabolism', minimum=0.0_rk, default=0.2_rk*h*sec_per_year, scale_factor=1._rk/sec_per_year)

   call self%get_parameter(z0, 'z0', 'yr-1', 'background mortality', minimum=0.0_rk, default=z0pre*w_inf**z0exp*sec_per_year, scale_factor=1._rk/sec_per_year)

   ! Determine size classes (log-spaced between size at birth and infinite size)
   allocate(self%logw(self%nclass))
   allocate(self%w(self%nclass))
   allocate(self%delta_w(self%nclass))
   delta_logw = (log(w_inf)-log(self%w_min))/(self%nclass-1)       ! Log-mass distance between size classes [constant across spectrum]
   do iclass=1,self%nclass
      self%logw(iclass) = log(self%w_min)+delta_logw*(iclass-1)    ! Log mass of each size class
   end do
   self%w = exp(self%logw)                                         ! Mass of each size class
   self%delta_w = exp(self%logw+delta_logw) - self%w               ! Mass difference between consecutive size classes

   ! Compute size-class-specific parameters that do not vary in time. It is most efficient to compute these once, during initialization.
   allocate(self%I_max(self%nclass))
   allocate(self%std_metab(self%nclass))
   allocate(self%mu_b(self%nclass))
   allocate(self%mu_s(self%nclass))
   allocate(self%F(self%nclass))
   allocate(self%psi(self%nclass))
   allocate(self%V(self%nclass))
   self%V(:) = gamma*self%w**(q-1)      ! specific volumetric search rate [m3 s-1 g-1] (mass-specific, hence the -1!)
   self%V = self%V * g_per_mmol_carbon  ! express volumetric search rate per mmol carbon (instead of per g)
   self%I_max(:) = h*self%w**(n-1)      ! specific maximum ingestion rate [s-1]; Eq M4, but specific, hence the -1!
   self%std_metab(:) = ks*self%w**(p-1) ! specific metabolism [s-1]; second term in Eq M7, but specific, hence the -1!
   select case (z0_type)
   case (0)
      self%mu_b(:) = z0                    ! background mortality [s-1]; Eq M11
   case (1)
      self%mu_b(:) = z0pre*self%w**z0exp
   end select
   if (w_s > 0.0_rk) then
      self%mu_s = z_spre*(self%w/w_s)**z_s  ! Blanchard et al. 10.1098/rstb.2012.0231 Table S1
   else
      self%mu_s = 0
   end if

   ! Fishing mortality
   self%F = 0.0_rk
   call self%get_parameter(fishing_type,'fishing_type', '', 'fishing regime (0: none, 1: constant/knife-edge, 2: logistic)',default=0, minimum=0, maximum=3)
   if (fishing_type > 0) call self%get_parameter(w_minF, 'w_minF', 'g', 'minimum mass for fishing selectivity', default=0.0_rk, minimum=0.0_rk)
   select case (fishing_type)
   case (1)
      ! constant
      call self%get_parameter(F, 'F', 'yr-1', 'fishing effort', default=0.0_rk, minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
      do iclass=1,self%nclass
         if (self%w(iclass) > w_minF) self%F(iclass) = F
      end do
   case (2)
      ! mizer fishing mortality [s-1]; Eqs M13 and M14 combined
      call self%get_parameter(F,  'F',  'yr-1', 'maximum fishing effort',                        default=0.0_rk, minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
      call self%get_parameter(S1, 'S1', '-',    'offset for fishing selectivity exponent',       default=0.0_rk, minimum=0.0_rk)
      call self%get_parameter(S2, 'S2', 'g-1',  'scale factor for fishing selectivity exponent', default=0.0_rk, minimum=0.0_rk)
      do iclass=1, self%nclass
         if (self%w(iclass) > w_minF) self%F(iclass) = F/(1+exp(S1-S2*self%w(iclass)))
      end do
   case (3)
      ! linearly increasing mortality as in Blanchard et al 2009 J Anim Ecol
      call self%get_parameter(F_a, 'F_a', 'yr-1 (log10 g)-1', 'scale factor for fishing mortality as function of log10 mass', default=0.0_rk, minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
      call self%get_parameter(F_b, 'F_b', 'yr-1', 'offset for fishing mortality as function of log10 mass', default=0.0_rk, minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
      do iclass=1, self%nclass
         if (self%w(iclass) > w_minF) self%F(iclass) = F_a * log10(self%w(iclass)) + F_b
      end do
   end select

   if (w_mat==0.0_rk) then
      ! No explicit reproduction as in original Blanchard community size spectrum model. Recruitment will be constant.
      self%psi(:) = 0
   else
      self%psi(:) = (self%w/w_inf)**(1-n) / (1+(self%w/w_mat)**(-10)) ! allocation to reproduction [-]; Eqs M13 and M14 combined
   end if
   if (report_statistics) then
      write (*,*) 'Specific search volume (yr-1):'
      write (*,*) '  @ mass = 1:',gamma*sec_per_year,'(gamma)'
      write (*,*) '  @ minimum mass:',self%V(1)*sec_per_year
      write (*,*) '  @ maximum mass:',self%V(self%nclass)*sec_per_year
      write (*,*) 'Specific ingestion rate (yr-1)'
      write (*,*) '  @ mass = 1:',h*sec_per_year,'(h)'
      write (*,*) '  @ minimum mass:',self%I_max(1)*sec_per_year
      write (*,*) '  @ maximum mass:',self%I_max(self%nclass)*sec_per_year
      write (*,*) 'Specific metabolism (yr-1):'
      write (*,*) '  @ mass = 1:',self%std_metab(1)*sec_per_year,'(ks)'
      write (*,*) '  @ minimum mass:',self%std_metab(1)*sec_per_year
      write (*,*) '  @ maximum mass:',self%std_metab(self%nclass)*sec_per_year
      write (*,*) 'Background mortality (yr-1):'
      select case (z0_type)
      case (0)
      write (*,*) '  @ mass = 1:',z0*sec_per_year,'(z0)'
      case (1)
         write (*,*) '  @ mass = 1:',z0pre*sec_per_year,'(z0pre)'
      end select
      write (*,*) '  @ minimum mass:',self%mu_b(1)*sec_per_year
      write (*,*) '  @ maximum mass:',self%mu_b(self%nclass)*sec_per_year
      write (*,*) 'Senescence mortality (yr-1):'
      write (*,*) '  @ minimum mass:',self%mu_s(1)*sec_per_year
      write (*,*) '  @ maximum mass:',self%mu_s(self%nclass)*sec_per_year
      write (*,*) 'Fishing mortality at minimum size:',self%F(1)*sec_per_year,'yr-1'
      write (*,*) 'Fishing mortality at maximum size:',self%F(self%nclass)*sec_per_year,'yr-1'
   end if

   allocate(pelagic_size_spectrum)
   call pelagic_size_spectrum%parameters%set_integer('nsource', self%nprey)
   call pelagic_size_spectrum%parameters%set_real('w_min', self%w_min/self%beta**2)
   call pelagic_size_spectrum%parameters%set_real('w_max', self%w_min)

   ! Register dependencies for all prey.
   ! If the population is cannibalistic, autoamtically add all our size classes to the set of prey types.
   if (cannibalism) self%nprey = self%nprey + self%nclass
   allocate(self%id_prey_c(self%nprey))
   allocate(self%id_prey_n(self%nprey))
   allocate(self%id_prey_p(self%nprey))
   allocate(self%id_prey_s(self%nprey))
   allocate(self%id_pelprey_c(self%nprey))
   allocate(total_pelprey_calculator)
   do iprey=1,self%nprey
      write (strindex,'(i0)') iprey

      call self%register_bottom_state_dependency(self%id_prey_c(iprey), 'prey'//trim(strindex)//'_c', 'mmol C m-3', 'average encountered concentration of carbon in prey '//trim(strindex))
      call self%register_bottom_state_dependency(self%id_prey_n(iprey), 'prey'//trim(strindex)//'_n', 'mmol N m-3', 'average encountered concentration of nitrogen in prey '//trim(strindex))
      call self%register_bottom_state_dependency(self%id_prey_p(iprey), 'prey'//trim(strindex)//'_p', 'mmol P m-3', 'average encountered concentration of phosphorus in prey '//trim(strindex))
      call self%register_bottom_state_dependency(self%id_prey_s(iprey), 'prey'//trim(strindex)//'_s', 'mmol Si m-3', 'average encountered concentration of silicate in prey '//trim(strindex))

      call self%request_coupling_to_model(self%id_prey_c(iprey), 'ave_prey'//trim(strindex), standard_variables%total_carbon)
      call self%request_coupling_to_model(self%id_prey_n(iprey), 'ave_prey'//trim(strindex), standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_prey_p(iprey), 'ave_prey'//trim(strindex), standard_variables%total_phosphorus)
      call self%request_coupling_to_model(self%id_prey_s(iprey), 'ave_prey'//trim(strindex), standard_variables%total_silicate)

      if (iprey <= self%nprey - self%nclass .or. .not. cannibalism) then
        ! Prey is pelagic
        call self%register_dependency(self%id_pelprey_c(iprey), 'pelprey_c'//trim(strindex), 'mmol C m-3', 'carbon in pelagic prey '//trim(strindex))
        call self%request_coupling_to_model(self%id_pelprey_c(iprey), 'prey'//trim(strindex), standard_variables%total_carbon)
        call total_pelprey_calculator%add_component('pelprey_c'//trim(strindex))
        allocate(depth_averaged_prey)
        call self%add_child(depth_averaged_prey, 'pelprey'//trim(strindex)//'_depth_average', configunit=-1)
        call depth_averaged_prey%couplings%set_string('source', '../prey'//trim(strindex))
        call depth_averaged_prey%request_coupling('w', '../total_pelprey_calculator/result')
        call depth_averaged_prey%request_coupling('w_int', '../w_integrator/result')
        call self%request_coupling(self%id_prey_c(iprey), 'pelprey'//trim(strindex)//'_depth_average/c')
        call self%request_coupling(self%id_prey_n(iprey), 'pelprey'//trim(strindex)//'_depth_average/n')
        call self%request_coupling(self%id_prey_p(iprey), 'pelprey'//trim(strindex)//'_depth_average/p')
        call self%request_coupling(self%id_prey_s(iprey), 'pelprey'//trim(strindex)//'_depth_average/s')
        call self%couplings%set_string('ave_prey'//trim(strindex), 'pelprey'//trim(strindex)//'_depth_average')
        call self%get_parameter(w_prey_min, 'w_prey'//trim(strindex)//'_min', 'g', 'minimum mass of prey '//trim(strindex))
        call self%get_parameter(w_prey_max, 'w_prey'//trim(strindex)//'_max', 'g', 'maximum mass of prey '//trim(strindex))
        call depth_averaged_prey%set_variable_property(depth_averaged_prey%id_c, 'min_particle_mass', w_prey_min)
        call depth_averaged_prey%set_variable_property(depth_averaged_prey%id_c, 'max_particle_mass', w_prey_max)
        call pelagic_size_spectrum%parameters%set_real('w_source'//trim(strindex)//'_min', w_prey_min)
        call pelagic_size_spectrum%parameters%set_real('w_source'//trim(strindex)//'_max', w_prey_max)
        call pelagic_size_spectrum%couplings%set_string('source'//trim(strindex), '../pelprey_c'//trim(strindex))
      end if
   end do
   call self%add_child(pelagic_size_spectrum, 'pelagic_size_spectrum', configunit=-1)

   call self%add_child(total_pelprey_calculator,'total_pelprey_calculator',configunit=-1)

   allocate(depth_integral)
   call self%add_child(depth_integral, 'w_integrator', configunit=-1)
   call depth_integral%request_coupling('source', '../total_pelprey_calculator/result')

   allocate(square)
   call self%add_child(square, 'w_squarer', configunit=-1)
   call square%request_coupling('source', '../total_pelprey_calculator/result')

   allocate(depth_integral)
   call self%add_child(depth_integral, 'w2_integrator', configunit=-1)
   call depth_integral%request_coupling('source', '../w_squarer/result')

   ! Allocate size-class-specific identifiers for abundance state variable and diagnostics.
   allocate(self%id_c(self%nclass))
   if (self%SRR == 1 .or. self%SRR == 2)  allocate(self%id_reproduction(self%nclass))
   allocate(self%id_f(self%nclass))
   allocate(self%id_g(self%nclass))
   do iclass=1, self%nclass
      ! Postfix for size-class-specific variable names (an integer number)
      write (strindex,'(i0)') iclass

      ! Register state variable, store associated individual mass (used by predators, if any, to determine grazing preference).
      call self%register_state_variable(self%id_c(iclass), 'c'//trim(strindex), 'mmol m-2', 'carbon in size class '//trim(strindex), c_ini, minimum=0.0_rk)
      call self%set_variable_property(self%id_c(iclass), 'particle_mass', self%w(iclass))

      ! Register this size class' contribution to total mass in the system (for mass conservation checks)
      call self%add_to_aggregate_variable(standard_variables%total_carbon,self%id_c(iclass))
      call self%add_to_aggregate_variable(standard_variables%total_nitrogen,self%id_c(iclass), scale_factor=self%qnc)
      call self%add_to_aggregate_variable(standard_variables%total_phosphorus,self%id_c(iclass), scale_factor=self%qpc)

      allocate(depth_averaged_class)
      depth_averaged_class%qnc = self%qnc
      depth_averaged_class%qpc = self%qpc
      call self%add_child(depth_averaged_class, 'class'//trim(strindex), configunit=-1)
      call depth_averaged_class%request_coupling('int_w', '../w_integrator/result')
      call depth_averaged_class%request_coupling('int_w2', '../w2_integrator/result')
      call depth_averaged_class%request_coupling('int_c', '../c'//trim(strindex))
      call depth_averaged_class%set_variable_property(depth_averaged_class%id_c,'particle_mass',self%w(iclass))

      ! If population is cannibalistic, add this size class as one of the prey (after the user-specified prey set).
      if (cannibalism) then
         write (strindex2,'(i0)') self%nprey - self%nclass + iclass
         call self%couplings%set_string('ave_prey'//trim(strindex2), './class'//trim(strindex))
      end if

      ! Register size-class-specific diagnostics
      if (self%SRR == 1 .or. self%SRR == 2) call self%register_diagnostic_variable(self%id_reproduction(iclass),'reproduction'//trim(strindex),'g m-2 d-1','allocation to reproduction in size class '//trim(strindex),         source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_f(iclass),           'f'//trim(strindex),           '-',        'functional response of size class '//trim(strindex),                source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_g(iclass),           'g'//trim(strindex),           'd-1',      'specific growth rate of individuals in size class '//trim(strindex),source=source_do_bottom)
   end do

   allocate(self%phi(self%nprey,self%nclass))

   ! Register a state variable for waste (faeces, maintenance, dead matter resulting from non-predation mortality, fraction of offspring that does not survive)
   call register_waste('egested_matter', 'cnps', self%id_waste_c, self%id_waste_n, self%id_waste_p, self%id_waste_s)
 !  call register_waste('respired_oxygen', 'c', id_c=self%id_o2)
   call register_waste('respired_carbon', 'c', id_c=self%id_dic)
   call register_waste('excreted_nitrogen', 'n', id_n=self%id_din)
   call register_waste('excreted_phosphorus', 'p', id_p=self%id_dip)

   call self%register_bottom_state_dependency(self%id_discard_c, 'discard_c', 'mmol C m-2', 'organic carbon discards')
   call self%register_bottom_state_dependency(self%id_discard_n, 'discard_n', 'mmol N m-2', 'organic nitrogen discards')
   call self%register_bottom_state_dependency(self%id_discard_p, 'discard_p', 'mmol P m-2', 'organic phosphorus discards')
   call self%request_coupling_to_model(self%id_discard_c, 'discards', standard_variables%total_carbon)
   call self%request_coupling_to_model(self%id_discard_n, 'discards', standard_variables%total_nitrogen)
   call self%request_coupling_to_model(self%id_discard_p, 'discards', standard_variables%total_phosphorus)
   call self%couplings%set_string('discards', './egested_matter')

   call self%register_state_variable(self%id_landings, 'landings', 'g m-2', 'landed biomass')
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_landings, scale_factor=1.0_rk/g_per_mmol_carbon)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_landings, scale_factor=self%qnc/g_per_mmol_carbon)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_landings, scale_factor=self%qpc/g_per_mmol_carbon)

   ! Register diagnostic for total offspring production across population.
   if (self%SRR == 1 .or. self%SRR == 2) then
      call self%register_diagnostic_variable(self%id_total_reproduction,'total_reproduction','mmol C m-2 d-1','total reproduction',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_p,'R_p','# m-2 d-1','density-independent recruitment',source=source_do_bottom)
   elseif (self%SRR == 3) then
      ! Infer biomass of lowest size class by extending spectrum of (small) prey
      call self%register_dependency(self%id_offset, 'prey_spectrum_offset', '-', 'offset of pelagic prey spectrum')
      call self%register_dependency(self%id_slope, 'prey_spectrum_slope', '-', 'slope of pelagic prey spectrum')
      call self%request_coupling(self%id_offset, './pelagic_size_spectrum/offset')
      call self%request_coupling(self%id_slope, './pelagic_size_spectrum/slope')
   end if
   call self%register_diagnostic_variable(self%id_R,'R','# m-2 d-1','recruitment',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_tot, 'c_tot', 'g m-2', 'total biomass', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_size1, 'c_size1', 'g m-2', 'fish biomass smaller than threshold1', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_size2, 'c_size2', 'g m-2', 'fish biomass smaller than threshold2', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_size3, 'c_size3', 'g m-2', 'fish biomass smaller than threshold3', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_lfi, 'c_lfi', 'g m-2', 'Large fish index', source=source_do_bottom)
   contains
   
   subroutine register_waste(name, composition, id_c, id_n, id_p, id_s)
      character(len=*), intent(in) :: name, composition
      type (type_bottom_state_variable_id), intent(inout), target, optional :: id_c, id_n, id_p, id_s

      class (type_depth_integrated_sink),pointer :: waste

      allocate(waste)
      call waste%parameters%set_string('composition', composition)
      call self%register_model_dependency(name=name)
      call self%add_child(waste, name//'_int', configunit=-1)
      if (present(id_c)) then
         call self%register_bottom_state_dependency(id_c, name//'_c', 'mmol C m-2', 'particulate organic carbon waste')
         call self%request_coupling(id_c, name//'_int/c_int')
      end if
      if (present(id_n)) then
         call self%register_bottom_state_dependency(id_n, name//'_n', 'mmol N m-2', 'particulate organic nitrogen waste')
         call self%request_coupling(id_n, name//'_int/n_int')
      end if
      if (present(id_p)) then
         call self%register_bottom_state_dependency(id_p, name//'_p', 'mmol P m-2', 'particulate organic phosphorus waste')
         call self%request_coupling(id_p, name//'_int/p_int')
      end if
      if (present(id_s)) then
         call self%register_bottom_state_dependency(id_s, name//'_s', 'mmol Si m-2','particulate silicate waste')
         call self%request_coupling(id_s, name//'_int/s_int')
      end if
      call waste%request_coupling('w', '../total_pelprey_calculator/result')
      call waste%request_coupling('w_int', '../w_integrator/result')
      call waste%couplings%set_string('target', '../'//name)
   end subroutine

   end subroutine initialize
!EOC

   subroutine after_coupling(self)
      class (type_multi_element_population),intent(inout) :: self

      integer           :: iprey, iclass, i
      character(len=10) :: strindex
      integer, parameter :: n = 100
      real(rk)          :: w_p, w_p_min, w_p_max, log_w_ps(n)

      ! Coupling with prey has completed.
      ! Now we can query all prey for their wet mass per individual. From that we precompute predator-prey preferences.
      do iprey=1, self%nprey
         ! First retrieve individual mass of prey (throw fatal error if not available)
         w_p = self%id_prey_c(iprey)%link%target%properties%get_real('particle_mass',-1._rk)
         w_p_min = self%id_prey_c(iprey)%link%target%properties%get_real('min_particle_mass', -1._rk)
         w_p_max = self%id_prey_c(iprey)%link%target%properties%get_real('max_particle_mass', -1._rk)
         if (w_p_min > 0 .and. w_p_max > 0) then
            ! A size range (minimum wet mass - maximum wet mass) for this prey is given.
            ! We assume that within this range, biomass is uniformly distributed in log wet mass - log total biomass space (= a Sheldon size spectrum)
            ! The effective preference is computed by averaging the (Gaussian) preference curve over the relevant interval.
            ! That is done below by discretizing the wet mass interval.

            ! Create log-spaced grid for prey size range
            do i = 1, n
               log_w_ps(i) = log(w_p_min) + (i - 0.5_rk) * (log(w_p_max) - log(w_p_min)) / n
            end do

            ! Compute size-class-specific preference for current prey; Eq 4 in Hartvig et al. 2011 JTB, but note sigma typo confirmed by KH Andersen
            ! This is a log-normal distribution of prey mass, scaled such that at optimum prey mass (=predator mass/beta), the preference equals 1.
            ! sigma is the standard deviation in ln mass units.
            do iclass=1, self%nclass
               self%phi(iprey, iclass) = sum(exp(-(log_w_ps - self%logw(iclass) + log(self%beta))**2 / self%sigma**2 / 2)) / n
            end do
         elseif (w_p > 0) then
            ! A single size (wet mass) for this prey is given.
            ! Compute size-class-specific preference for current prey; Eq 4 in Hartvig et al. 2011 JTB, but note sigma typo confirmed by KH Andersen
            ! This is a log-normal distribution of prey mass, scaled such that at optimum prey mass (=predator mass/beta), the preference equals 1.
            ! sigma is the standard deviation in ln mass units.
            do iclass=1, self%nclass
               self%phi(iprey, iclass) = exp(-(log(w_p) - self%logw(iclass) + log(self%beta))**2 / self%sigma**2 / 2)
            end do
         else
            write (strindex, '(i0)') iprey
            call self%fatal_error('after_coupling', 'prey '//trim(strindex)//' does not have attribute "particle_mass" or the combination "min_particle_mass", "max_particle_mass".')
         end if
      end do
      if (any(self%phi < 0)) call self%fatal_error('after_coupling', 'one or more entries of phi are < 0')

   end subroutine after_coupling

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_multi_element_population),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iclass,iprey,istate
      real(rk) :: c_lfi, c_size1,c_size2, c_size3, slope, offset, endpoint, expected_eggs
      real(rk) :: E_e,E_a_c,E_a_n,E_a_p,E_a_s,f,total_reproduction,T_lim,temp,T_w_int,w_int,g_tot,R,R_p,nflux(0:self%nclass),prey_state,g_tot_c,g_tot_n,g_tot_p
      real(rk),dimension(self%nprey)  :: prey_c,prey_n,prey_p,prey_s,prey_loss
      real(rk),dimension(self%nclass) :: Nw,I_c,I_n,I_p,I_s,maintenance,g,mu,reproduction
      real(rk), parameter :: delta_t = 900

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve size-class-specific abundances
         c_size1 = 0
         c_size2 = 0
         c_size3 = 0
         c_lfi = 0
         do iclass=1,self%nclass
            _GET_HORIZONTAL_(self%id_c(iclass), Nw(iclass))
            if (self%w(iclass) < self%w_threshold) c_size1 = c_size1 + Nw(iclass)
            if (self%w(iclass) < self%w_threshold2) c_size2 = c_size2 + Nw(iclass)
            if (self%w(iclass) < self%w_threshold3) c_size3 = c_size3 + Nw(iclass)
            if (self%w(iclass) > self%lfi_w_threshold) c_lfi = c_lfi + Nw(iclass)
         end do
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c_tot, sum(Nw)*g_per_mmol_carbon)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c_size1, c_size1*g_per_mmol_carbon)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c_size2, c_size2*g_per_mmol_carbon)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c_size3, c_size3*g_per_mmol_carbon)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_c_lfi, (c_lfi/sum(Nw))*g_per_mmol_carbon)

         ! Retrieve prey abundances
         do iprey=1,self%nprey
            _GET_HORIZONTAL_(self%id_prey_c(iprey), prey_c(iprey))
            _GET_HORIZONTAL_(self%id_prey_n(iprey), prey_n(iprey))
            _GET_HORIZONTAL_(self%id_prey_p(iprey), prey_p(iprey))
            _GET_HORIZONTAL_(self%id_prey_s(iprey), prey_s(iprey))
         end do

         ! Temperature limitation factor affecting all rates (not in Blanchard et al.)
         if (self%T_dependence==1) then
            _GET_HORIZONTAL_(self%id_T_w_int, T_w_int)
            _GET_HORIZONTAL_(self%id_w_int, w_int)
            temp = T_w_int/w_int
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_T_ave, temp)
            T_lim = exp(self%c1-self%E_A/Boltzmann/(temp+Kelvin))
         else
            T_lim = 1
         end if

         ! Food uptake (all size classes, all prey types)
         ! This computes total ingestion per size class (over all prey), and total loss per prey type (over all size classes)
         prey_loss = 0.0_rk
         do iclass=1,self%nclass
            ! Compute total prey availability (concentration summed over all prey, scaled with prey-specific preference)
            E_a_c = sum(self%phi(:,iclass)*prey_c)
            E_a_n = sum(self%phi(:,iclass)*prey_n)
            E_a_p = sum(self%phi(:,iclass)*prey_p)
            E_a_s = sum(self%phi(:,iclass)*prey_s)
#ifndef NDEBUG
            if (isnan(E_a_c)) &
               call self%fatal_error('do_bottom','E_a_c is nan')
#endif

            ! Compute actual encounter in mmol prey s-1 (mmol C in predator)-1
            ! availability per volume (mmol prey m-3) times mass-specific volumetric search rate (m3 mmol predator-1 s-1)
            E_e = self%V(iclass)*E_a_c ! Eq M3

            ! Compute ingestion rate (s-1) - mmol C of prey per mmol C in predator per time!
            f = E_e/(E_e+self%I_max(iclass))   ! Eq M5
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f(iclass),f)
            I_c(iclass) = T_lim*self%I_max(iclass)*f ! ingestion part of M7
#ifndef NDEBUG
            if (isnan(I_c(iclass))) &
               call self%fatal_error('do_bottom','ingestion is nan')
#endif
            ! Predator-specific ingestion of nitrogen, phosphorus, silicate (mmol prey s-1 mmol-1)
            I_n(iclass) = I_c(iclass)/E_a_c*E_a_n
            I_p(iclass) = I_c(iclass)/E_a_c*E_a_p
            I_s(iclass) = I_c(iclass)/E_a_c*E_a_s

            ! Account for this size class' ingestion in specific loss rate of all prey
            ! Units go from (mmol prey s-1 g-1) to (m s-1) - we divide by prey concentration and multiply by predator areal density
            prey_loss(:) = prey_loss(:) + I_c(iclass)/E_a_c*self%phi(:,iclass)*Nw(iclass)
         end do

#ifndef NDEBUG
         if (any(prey_loss<0)) &
            call self%fatal_error('do_bottom','prey_loss is negative')
         if (any(prey_loss>1/delta_t)) &
            call self%fatal_error('do_bottom','prey_loss is high (prey will be <0 within time step)')
#endif

         ! Initialize size-class-specific mortality (s-1) with precomputed size-dependent background value.
         mu = self%mu_b*T_lim + self%mu_s

         ! Individual physiology (per size class)
         do iclass=1,self%nclass
            ! Specific maintenance rate (s-1)
            maintenance(iclass) = T_lim*self%std_metab(iclass)

            ! Net specific element availability (mmol prey s-1 mmol-1)
            g_tot_c = self%alpha*I_c(iclass) - maintenance(iclass)
            g_tot_n = self%alpha*I_n(iclass)
            g_tot_p = self%alpha*I_p(iclass)

            ! Specific growth rate (s-1) is minimum supported by different resources.
            g_tot = min(g_tot_c, g_tot_n/self%qnc, g_tot_p/self%qpc)

            ! Avoid shrinking: limit maintenance to maximum sustainable value and increase starvation mortality.
            maintenance(iclass) = min(maintenance(iclass),self%alpha*I_c(iclass))
            mu(iclass) = mu(iclass) + max(0.0_rk,-g_tot/self%w(iclass)/self%xi)
            g_tot = max(0.0_rk,g_tot)

            ! Individual growth (s-1)
            g(iclass) = (1-self%psi(iclass))*g_tot ! Eq M7

            ! Mass flux towards reproduction (mmol C m-2 s-1) - sum over all individuals in this size class
            reproduction(iclass) = self%psi(iclass)*g_tot*Nw(iclass)
         end do

         ! Compute number of individuals moving from each size class to the next (units: # s-1)
         nflux(1:self%nclass) = Nw*g/(self%delta_w/g_per_mmol_carbon)

         ! Sum reproductive output of entire population in # m-2 s-1 (Eq 10 of Hartvig et al. 2011 JTB)
         ! Note: division by 2 is the result of the fact that reproductive output applies to females only,
         ! which are assumed to be 50% of the population.
         total_reproduction = sum(reproduction)
         R_p = self%erepro/2 * total_reproduction / (self%w_min/g_per_mmol_carbon)

         ! Use stock-recruitment relationship to translate density-independent recruitment into actual recruitment (units: # s-1)
         if (self%SRR==0) then
            ! Constant recruitment
            R = self%recruitment
         elseif (self%SRR==1) then
            ! Density-independent recruitment
            R = R_p
         elseif (self%SRR==2) then
            ! Beverton-Holt recruitment
            R = self%R_max*R_p/(R_p + self%R_max)
         else
            _GET_HORIZONTAL_(self%id_offset, offset)
            _GET_HORIZONTAL_(self%id_slope, slope)
            endpoint = offset + slope * log(self%w_min)
            expected_eggs = exp(endpoint) * self%delta_w(1)
            R = max(expected_eggs - Nw(1), 0._rk)/(self%w_min/g_per_mmol_carbon) * self%R_relax
         end if

         ! Use recruitment as number of incoming individuals for the first size class.
         nflux(0) = R

         ! Destroy all prey constituents that we are aware of (we only need to destroy carbon
         ! to get the correct impact on prey, but by destroying all we enable conservation checks)
         do iprey = 1, self%nprey
            if (iprey > self%nprey - self%nclass) then
               ! Prey is one of our own size classes
               _SET_BOTTOM_ODE_(self%id_prey_c(iprey), -prey_loss(iprey)*prey_c(iprey))
            elseif (self%feedback) then
               ! Prey is an external variable (not one of our size classes)
               _SET_BOTTOM_ODE_(self%id_prey_c(iprey), -prey_loss(iprey)*prey_c(iprey))
               _SET_BOTTOM_ODE_(self%id_prey_n(iprey), -prey_loss(iprey)*prey_n(iprey))
               _SET_BOTTOM_ODE_(self%id_prey_p(iprey), -prey_loss(iprey)*prey_p(iprey))
               _SET_BOTTOM_ODE_(self%id_prey_s(iprey), -prey_loss(iprey)*prey_s(iprey))
            end if
         end do

         ! Transfer size-class-specific source terms and diagnostics to FABM
         do iclass=1,self%nclass
            ! Apply specific mortality (s-1) to size-class-specific abundances and apply upwind advection - this is a time-explicit version of Eq G.1 of Hartvig et al.
            _SET_BOTTOM_ODE_(self%id_c(iclass),-(mu(iclass) + self%F(iclass))*Nw(iclass) + (nflux(iclass-1)-nflux(iclass))*self%w(iclass)/g_per_mmol_carbon)

            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_g(iclass),g(iclass)*86400)
            if (self%SRR == 1 .or. self%SRR == 2) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_reproduction(iclass),reproduction(iclass)*86400)
            end if
         end do

         if (self%SRR == 1 .or. self%SRR == 2) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_total_reproduction,total_reproduction*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_p,R_p*86400)
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R,R*86400)

         ! Compute waste fluxes: total ingestion plus mortality, minus mass used in growth, minus recruitment, plus growth over right edge of resolved size range.
         if (self%feedback) then
  !          _SET_BOTTOM_ODE_(self%id_o2,-self%resp_o2C*(1-self%alpha-self%alpha_eg)*sum(I_c*Nw))
            _SET_BOTTOM_ODE_(self%id_dic,(1-self%alpha-self%alpha_eg)*sum(I_c*Nw))
            _SET_BOTTOM_ODE_(self%id_din,(1-self%alpha-self%alpha_eg)*sum(I_n*Nw))
            _SET_BOTTOM_ODE_(self%id_dip,(1-self%alpha-self%alpha_eg)*sum(I_p*Nw))
            _SET_BOTTOM_ODE_(self%id_waste_c,sum(((self%alpha+self%alpha_eg)*I_c +  mu - g/(1-self%psi)          )*Nw) - R*self%w_min/g_per_mmol_carbon          + nflux(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon)
            _SET_BOTTOM_ODE_(self%id_waste_n,sum(((self%alpha+self%alpha_eg)*I_n + (mu - g/(1-self%psi))*self%qnc)*Nw) - R*self%w_min/g_per_mmol_carbon*self%qnc + nflux(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qnc)
            _SET_BOTTOM_ODE_(self%id_waste_p,sum(((self%alpha+self%alpha_eg)*I_p + (mu - g/(1-self%psi))*self%qpc)*Nw) - R*self%w_min/g_per_mmol_carbon*self%qpc + nflux(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qpc)
            _SET_BOTTOM_ODE_(self%id_waste_s,sum(I_s*Nw))
         end if
         _SET_BOTTOM_ODE_(self%id_landings,sum(self%F*Nw*g_per_mmol_carbon))
      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

!-----------------------------------------------------------------------

end module mizer_multi_element_population

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2017
!-----------------------------------------------------------------------
