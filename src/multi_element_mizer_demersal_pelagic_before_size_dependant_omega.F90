#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mizer
!
! !INTERFACE:
module mizer_multi_element_demersal_pelagic_population
       
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

!Hp note: Here MIZER is run in units of mmol C m-2 s-1. Model outputs for fish (fish_c_tot, fish_c_size1,..., fish_landings) are given in units of wet weight. Diagnostics are converted individually to per day either here or in multi element support. Unit conversion between mmolC and g wet weight are given below. 
!
! !USES:
   use fabm_types
   use fabm_particle

   use multi_element_support

   implicit none

!  default: all is private.
   private
   
   type type_prey
      type (type_model_id)                          :: id_source
      type (type_horizontal_dependency_id)          :: id_c,id_n,id_p,id_s
      type (type_horizontal_dependency_id)          :: id_c_an
      type (type_dependency_id)                     :: id_c_pel,id_n_pel,id_p_pel,id_s_pel
    !  type (type_horizontal_diagnostic_variable_id) :: id_fc,id_fn,id_fp,id_fs

      ! To achieve compatibility with legacy ERSEM, we need to be able to decouple the variable
      ! from which food availability is derived from the variable that absorbs the loss due to
      ! gross food uptake. The following variables absorb the loss due to food uptake - by default
      ! they are coupled to the same variable from which available food is derived.
      type (type_model_id) :: id_loss_source

      logical  :: isben
   end type
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_particle_model),public :: type_multi_element_demersal_pelagic_population
      ! Variable identifiers
      type (type_bottom_state_variable_id),         allocatable :: id_c(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_c(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_n(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_p(:)
      type (type_bottom_state_variable_id),         allocatable :: id_prey_s(:)
      
      type (type_prey), allocatable :: prey(:)
      
      
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
      
      type (type_state_variable_id)                      :: id_beno2
      type (type_state_variable_id)                      :: id_bendic
      type (type_state_variable_id)                      :: id_bendip
      type (type_state_variable_id)                      :: id_bendin
      type (type_state_variable_id)                      :: id_bendiscard_c
      type (type_state_variable_id)                      :: id_bendiscard_n
      type (type_state_variable_id)                      :: id_bendiscard_p
      type (type_state_variable_id)                      :: id_bendiscard_s
      type (type_dependency_id)                          :: id_ETW          
      
      
      
      
      type (type_horizontal_diagnostic_variable_id)             :: id_total_reproduction ! Total reproduction
      type (type_horizontal_diagnostic_variable_id)             :: id_R_p                ! Density-independent recruitment
      type (type_horizontal_diagnostic_variable_id)             :: id_R                ! Density-dependent recruitment
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_reproduction(:)    ! Reproduction per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_f_pel(:)               ! Functional response per size class    
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_g_pel(:)               ! Specific growth rate per size class
      type (type_dependency_id),                    allocatable :: id_pelprey_c(:)
      type (type_horizontal_dependency_id),                    allocatable :: id_benprey_c(:)
      type (type_horizontal_dependency_id),                    allocatable :: id_benprey_n(:)
      type (type_horizontal_dependency_id),                    allocatable :: id_benprey_p(:)
      type (type_horizontal_dependency_id),                    allocatable :: id_benprey_s(:)
               
    !  type (type_horizontal_diagnostic_variable_id)             :: id_total_reproduction_ben ! Total reproduction
     ! type (type_horizontal_diagnostic_variable_id)             :: id_R_p_ben                ! Density-independent recruitment
     ! type (type_horizontal_diagnostic_variable_id)             :: id_R_ben                 ! Density-dependent recruitment
     ! type (type_horizontal_diagnostic_variable_id),allocatable :: id_reproduction_ben(:)    ! Reproduction per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_f_ben(:)               ! Functional response per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_g_ben(:)               ! Specific growth rate per size class


      type (type_horizontal_dependency_id) :: id_slope, id_offset
      type (type_horizontal_dependency_id) :: id_benslope, id_benoffset

      type (type_horizontal_dependency_id)                      :: id_T_w_int
      type (type_horizontal_dependency_id)                      :: id_w_int
      type (type_horizontal_diagnostic_variable_id)             :: id_T_ave, id_omega
      type (type_horizontal_diagnostic_variable_id)             :: id_c_tot
      type (type_horizontal_diagnostic_variable_id)             :: id_c_size1
      type (type_horizontal_diagnostic_variable_id)             :: id_c_size2
      type (type_horizontal_diagnostic_variable_id)             :: id_c_size3
      type (type_horizontal_diagnostic_variable_id)             :: id_c_lfi
      type (type_horizontal_diagnostic_variable_id)             :: id_test
      type (type_horizontal_diagnostic_variable_id), allocatable:: id_bprey_c(:), id_bprey_n(:),id_bprey_p(:),id_bprey_s(:)
      type (type_horizontal_diagnostic_variable_id)             :: id_fish_benDIP,id_fish_benDIC, id_fish_benDIN, id_benO2_fish
      type (type_horizontal_diagnostic_variable_id)             :: id_fish_benPOC, id_fish_benPON, id_fish_benPOP, id_fish_benPOS

      real(rk)                                                  :: w_threshold
      real(rk)                                                  :: w_threshold2
      real(rk)                                                  :: w_threshold3
      real(rk)                                                  :: lfi_w_threshold


      ! Number of size classes and prey
      integer :: nclass
      integer :: nprey, npelprey,nbenprey

      ! Size class characteristics
      real(rk),allocatable :: logw(:)     ! log mass
      real(rk),allocatable :: w(:)        ! mass
      real(rk),allocatable :: delta_w(:)  ! mass difference between consecutive size classes

      ! Size-class-independent parameters
      real(rk) :: w_min, w_min_ben       ! egg mass
      real(rk) :: alpha       ! assimilation efficiency
      real(rk) :: beta        ! preferred predator:prey mass ratio
      real(rk) :: sigma       ! s.d. of lognormal prey size selection function
      real(rk) :: erepro      ! efficiency of egg production
      real(rk) :: xi          ! fraction of mass consisting of lipid reserve (which can fuel maintenance)
      integer  :: SRR         ! type of stock-recruitment relationship (SRR)
      real(rk) :: recruitment   ! constant recruitment flux (SSR=0)
      real(rk) :: R_max       ! maximum recruitment flux (SSR=2)
      real(rk) :: R_relax     ! rate of relaxation towards expected egg density - derived from prey spectrum (SSR=3)
      real(rk) :: qnc         ! nitrogen:carbon ratio (mol:mol, constant)
      real(rk) :: qpc         ! phosphorus:carbon ratio (mol:mol, constant)
      real(rk) :: alpha_eg    ! fraction of food that is egested
      real(rk) :: omega    ! fraction of time that fish spend in pelagic

      integer  :: T_dependence ! Type of temperature dependence (0: none, 1: Arrhenius)
      real(rk) :: c1           ! Reference constant in Arrhenius equation = E_a/k/(T_ref+Kelvin)
      real(rk) :: E_a          ! Activation energy (eV)
      real(rk) :: resp_o2C      ! oxygen used per carbon respired (mol:mol,constant)

      logical :: feedback
      logical :: isben        ! Argument to define whether resource is pelagic or benthic
      logical :: emergent_omega

      ! Size-class-dependent parameters that will be precomputed during initialization
      real(rk), allocatable :: V(:)            ! volumetric search rate (Eq M2)
      real(rk), allocatable :: VB(:)            ! benthic volumetric search rate (Eq M2)
      real(rk), allocatable :: I_max(:)        ! maximum ingestion rate (Eq M4)
      real(rk), allocatable :: std_metab(:)    ! standard metabolism (k*w^p in Eq M7)
      real(rk), allocatable :: mu_b(:)         ! background mortality (temperature dependent)
      real(rk), allocatable :: mu_s(:)         ! senescence mortality (temperature independent)
      real(rk), allocatable :: F(:)            ! fishing mortality
      real(rk), allocatable :: psi(:)          ! allocation to reproduction
      real(rk), allocatable :: phi(:,:)        ! prey preference
      real(rk), allocatable :: pelspec(:)        ! prey available to pelagic fish
      real(rk), allocatable :: benspec(:)        ! prey available to benthic fish
            

      real(rk) :: c0_proxy_width
   contains
      procedure :: initialize
      procedure :: do_bottom
      procedure :: after_coupling
   end type type_multi_element_demersal_pelagic_population

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable),parameter :: total_mass = type_bulk_standard_variable(name='total_mass',units='g',aggregate_variable=.true.,conserved=.true.)

   real(rk) :: Kelvin = 273.15_rk      ! offset of Celsius temperature scale (K)
   real(rk) :: Boltzmann = 8.62e-5_rk  ! Boltzmann constant
   real(rk) :: g_per_mmol_carbon = 0.12_rk ! assume 10% of wet mass is carbon. 1 mmol carbon = 0.012 g, so that is equivalent to 0.12 g wet mass
   real(rk) :: CMass = 12.011_rk !Taken from ERSEM shared.F90
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
   class (type_multi_element_demersal_pelagic_population), intent(inout),target :: self
   integer,                               intent(in )          :: configunit
!
! !LOCAL VARIABLES:
   integer            :: b 
   integer            :: iclass, iprey
   logical            :: cannibalism, biomass_has_prey_unit
   real(rk)           :: delta_logw
   real(rk)           :: k_vb,n,q,qb, p,w_mat,w_inf,gamma,gammaB,h,ks,f0,z0,omega
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
   !class (type_demersal_size_spectrum),pointer :: demersal_size_spectrum
   logical, parameter :: report_statistics = .false.
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read parameters
   ! All rate coefficients are converted from d-1 to s-1 ensure source terms are directly compatible with FABM.
   ! Default values taken from Table S5
   call self%get_parameter(c_ini,'c_ini', 'mmol C/m2',     'initial density per size class', default=0.0_rk)
   call self%get_parameter(self%nclass,'nclass', '',     'number of size classes', default=100)
   call self%get_parameter(self%npelprey, 'npelprey',  '',     'number of pelagic prey')
   call self%get_parameter(self%nbenprey, 'nbenprey',  '',     'number of benthic prey')
   call self%get_parameter(self%nprey, 'nprey',  '',     'number of total prey') 
   call self%get_parameter(self%alpha, 'alpha',  '-',    'assimilation efficiency',            default=0.6_rk,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%alpha_eg, 'alpha_eg',  '-',    'fraction of food egested', default=1-self%alpha,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%erepro,'erepro', '-',    'reproductive efficiency',            default=1.0_rk,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%w_min, 'w_min',  'g',    'egg mass',                           default=0.001_rk, minimum=0.0_rk)
   call self%get_parameter(self%w_min_ben, 'w_min_ben',  'g',    'egg mass',                           default=0.001_rk, minimum=0.0_rk)
   call self%get_parameter(n,          'n',      '-',    'exponent of max. consumption',       default=2.0_rk/3.0_rk)
   call self%get_parameter(q,          'q',      '-',    'exponent of search volume',          default=0.8_rk)
   call self%get_parameter(qB,          'qB',      '-',    'benthic exponent of search volume',          default=q)
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
    !  call self%get_parameter(self%recruitment_ben, 'recruitment_ben', '# m-2 yr-1', 'constant recruitment flux', minimum=0.0_rk, default=kappa*self%w_min_ben**(-lambda)*sec_per_year, scale_factor=1._rk/sec_per_year)
   case (2)
      call self%get_parameter(self%R_max, 'R_max','# m-2 yr-1','maximum recruitment flux', minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   case (3)
      call self%get_parameter(self%R_relax, 'R_relax', 'd-1', 'relaxation towards expected egg density', minimum=0.0_rk, default=1.0_rk, scale_factor=1._rk/86400)
   end select

   call self%get_parameter(cannibalism,'cannibalism','',         'whether to enable intraspecific predation', default=.true.)
   if (cannibalism) call self%get_parameter(biomass_has_prey_unit, 'biomass_has_prey_unit', '', 'biomass has the same unit as prey', default=.true.)
   call self%get_parameter(self%qnc,   'qnc',        'mol mol-1','nitrogen to carbon ratio', default=16.0_rk/106.0_rk)
   call self%get_parameter(self%qpc,   'qpc',        'mol mol-1','phosphorus to carbon ratio', default=1.0_rk/106.0_rk)
   call self%get_parameter(self%resp_o2C,   'resp_o2C',        'mol mol-1','oxygen consumed per carbon respired', default=127_rk/106.0_rk)
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
      


      call self%register_dependency(self%id_ETW,standard_variables%temperature)

      
   end select

   ! Pre-factor for maximum ingestion rate derived from von Bertalanffy growth rate and background feeding level as described in appendix B
   h = 3*k_vb*w_inf**(1.0_rk/3.0_rk)/self%alpha/f0
   call self%get_parameter(h, 'h', 'yr-1 g^(-n)', 'pre-factor for maximum ingestion rate', minimum=0.0_rk, default=h*sec_per_year, scale_factor=1._rk/sec_per_year)

   ! Pre-factor for volumetric search rate (Eq M2)
   gamma = f0*h*self%beta**(2-lambda)/((1-f0)*sqrt(2*pi)*kappa*self%sigma)
   !gamma = f0*h*self%beta**(2-lambda)/((1-f0)*sqrt(2*pi)*kappa*self%sigma*exp((lambda-2)**2 * self%sigma**2 / 2)) ! add exp term taken from actual R code
   call self%get_parameter(gamma, 'gamma', 'PV yr-1 g^(-q)', 'pre-factor for volumetric search rate', minimum=0.0_rk, default=gamma*sec_per_year, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(gammaB, 'gammaB', 'PV yr-1 g^(-q)', 'pre-factor for volumetric search rate', minimum=0.0_rk, default=gamma*sec_per_year, scale_factor=1._rk/sec_per_year)

   ! Allow user override of standard metabolism pre-factor (e.g., Blanchard community size spectrum model has ks=0)
   call self%get_parameter(ks, 'ks', 'yr-1 g^(-p)', 'pre-factor for standard metabolism', minimum=0.0_rk, default=0.2_rk*h*sec_per_year, scale_factor=1._rk/sec_per_year)

   call self%get_parameter(z0, 'z0', 'yr-1', 'background mortality', minimum=0.0_rk, default=z0pre*w_inf**z0exp*sec_per_year, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(self%emergent_omega,'emergent_omega','-',          'use emergent omega', default=.false.)
   if (self%emergent_omega .eqv. .false. )    call self%get_parameter(self%omega,      'omega',  '-',    'fraction of time that fish spend in the pelagic',  minimum=0.0_rk, maximum=1.0_rk,default=1.0_rk)
   
   
    call self%register_diagnostic_variable(self%id_omega, 'omega' ,'-', 'time fish spend in the pelagic', source=source_do_bottom)
    
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
   allocate(self%VB(self%nclass))
   self%V(:) = gamma*self%w**(q-1)      ! specific volumetric search rate [m3 s-1 g-1] (mass-specific, hence the -1!)
   self%V = self%V * g_per_mmol_carbon  ! express volumetric search rate per mmol carbon (instead of per g)
   self%VB(:) = gammaB*self%w**(qB-1)      ! specific volumetric search rate [m3 s-1 g-1] (mass-specific, hence the -1!)
   self%VB = self%VB * g_per_mmol_carbon  ! express volumetric search rate per mmol carbon (instead of per g)
   
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
   call pelagic_size_spectrum%parameters%set_integer('nsource', self%npelprey)
   call pelagic_size_spectrum%parameters%set_real('w_min', self%w_min/self%beta**2)
   call pelagic_size_spectrum%parameters%set_real('w_max', self%w_min)
    
  ! allocate (demersal_size_spectrum)
  ! call demersal_size_spectrum%parameters%set_integer('nsource', self%nbenprey)
  ! call demersal_size_spectrum%parameters%set_integer('nsource_start', self%npelprey)
  ! call demersal_size_spectrum%parameters%set_real('w_min', self%w_min_ben/self%beta**2)
  ! call demersal_size_spectrum%parameters%set_real('w_max', self%w_min_ben)
   
   allocate (self%id_bprey_c(self%nbenprey ))
   allocate (self%id_bprey_n(self%nbenprey ))
   allocate (self%id_bprey_p(self%nbenprey))
   allocate (self%id_bprey_s(self%nbenprey))

   ! Register dependencies for all prey.
   ! If the population is cannibalistic, autoamtically add all our size classes to the set of prey types.
   allocate(self%prey(self%nprey))
   if (cannibalism) self%nprey = self%nprey + self%nclass
  ! allocate(self%id_prey(self%nprey))
   allocate(self%id_prey_c(self%nprey))
   allocate(self%id_prey_n(self%nprey))
   allocate(self%id_prey_p(self%nprey))
   allocate(self%id_prey_s(self%nprey))
   allocate(self%id_pelprey_c(self%npelprey))
   allocate (self%id_benprey_c(self%nbenprey))
   allocate (self%id_benprey_n(self%nbenprey))
   allocate (self%id_benprey_p(self%nbenprey))
   allocate (self%id_benprey_s(self%nbenprey))
   allocate(total_pelprey_calculator)

   allocate (self%pelspec(self%nprey))
   allocate (self%benspec(self%nprey))
   
   b=0._rk
   do iprey=1,self%nprey
      write (strindex,'(i0)') iprey
         !Compute whether prey is avaiable to pelagic and benthic fish.  Only plankton/benthos will not be available.
      self%pelspec(iprey)=1.0_rk
      self%benspec(iprey)=1.0_rk
      call self%register_bottom_state_dependency(self%id_prey_c(iprey), 'prey'//trim(strindex)//'_c', 'mmol C m-3', 'average encountered concentration of carbon in prey '//trim(strindex))
      call self%register_bottom_state_dependency(self%id_prey_n(iprey), 'prey'//trim(strindex)//'_n', 'mmol N m-3', 'average encountered concentration of nitrogen in prey '//trim(strindex))
      call self%register_bottom_state_dependency(self%id_prey_p(iprey), 'prey'//trim(strindex)//'_p', 'mmol P m-3', 'average encountered concentration of phosphorus in prey '//trim(strindex))
      call self%register_bottom_state_dependency(self%id_prey_s(iprey), 'prey'//trim(strindex)//'_s', 'mmol Si m-3', 'average encountered concentration of silicate in prey '//trim(strindex))

      call self%request_coupling_to_model(self%id_prey_c(iprey), 'ave_prey'//trim(strindex), standard_variables%total_carbon)
      call self%request_coupling_to_model(self%id_prey_n(iprey), 'ave_prey'//trim(strindex), standard_variables%total_nitrogen)
      call self%request_coupling_to_model(self%id_prey_p(iprey), 'ave_prey'//trim(strindex), standard_variables%total_phosphorus)
      call self%request_coupling_to_model(self%id_prey_s(iprey), 'ave_prey'//trim(strindex), standard_variables%total_silicate)
      
     

      if (iprey <= self%nprey - self%nclass .or. .not. cannibalism) then
        call self%get_parameter(self%prey(iprey)%isben,'prey'//trim(strindex)//'isben','','prey '//trim(strindex)//' is benthic',default=.false.)

        if (self%prey(iprey)%isben) then
           !b=self%nprey - self%nclass-self%npelprey
           b=b+1
           write (strindex2,'(i0)') b
      ! prey is benthic
        !   call self%register_dependency(self%id_pelprey_c(iprey), 'pelprey_c'//trim(strindex), 'mmol C m-2', 'carbon in demersal prey '//trim(strindex))
      !     call self%request_coupling_to_model(self%id_pelprey_c(iprey), 'prey'//trim(strindex), standard_variables%total_carbon)    
          ! call total_pelprey_calculator%add_component('pelprey_c'//trim(strindex)) 
         !  call self%register_horizontal_dependency(self%id_benprey_c(b), 'benprey_c'//trim(strindex2), 'mmol C m-2', 'carbon in benthic prey '//trim(strindex2)) 
         !  call self%register_horizontal_dependency(self%id_benprey_n(b), 'benprey_n'//trim(strindex2), 'mmol C m-2', 'carbon in benthic prey '//trim(strindex2)) 
         !  call self%register_horizontal_dependency(self%id_benprey_p(b), 'benprey_p'//trim(strindex2), 'mmol C m-2', 'carbon in benthic prey '//trim(strindex2)) 
         !  call self%register_horizontal_dependency(self%id_benprey_s(b), 'benprey_s'//trim(strindex2), 'mmol C m-2', 'carbon in benthic prey '//trim(strindex2)) 
          ! call self%request_coupling_to_model(self%id_benprey_c(b), 'prey'//trim(strindex), standard_variables%total_carbon)
        !   call self%request_coupling_to_model(self%id_benprey_n(b), 'prey'//trim(strindex), standard_variables%total_nitrogen)
        !   call self%request_coupling_to_model(self%id_benprey_p(b), 'prey'//trim(strindex), standard_variables%total_phosphorus)
        !   call self%request_coupling_to_model(self%id_benprey_s(b), 'prey'//trim(strindex), standard_variables%total_silicate)
           !call self%request_coupling(self%id_prey_c(iprey), 'benprey'//trim(strindex2)//'_depth_average/c')
           !call self%request_coupling(self%id_prey_n(iprey), 'benprey'//trim(strindex2)//'_depth_average/n')
           !call self%request_coupling(self%id_prey_p(iprey), 'benprey'//trim(strindex2)//'_depth_average/p')
           !call self%request_coupling(self%id_prey_s(iprey), 'benprey'//trim(strindex2)//'_depth_average/s')
           call self%couplings%set_string('ave_prey'//trim(strindex), 'prey'//trim(strindex))
           call self%request_coupling_to_model(self%id_prey_c(iprey), 'prey'//trim(strindex), standard_variables%total_carbon)
           call self%request_coupling_to_model(self%id_prey_n(iprey), 'prey'//trim(strindex), standard_variables%total_nitrogen)
           call self%request_coupling_to_model(self%id_prey_p(iprey), 'prey'//trim(strindex), standard_variables%total_phosphorus)
           call self%request_coupling_to_model(self%id_prey_s(iprey), 'prey'//trim(strindex), standard_variables%total_silicate)
         

         
         
         
           call self%get_parameter(w_prey_min, 'w_prey'//trim(strindex)//'_min', 'g', 'minimum mass of prey '//trim(strindex))
           call self%get_parameter(w_prey_max, 'w_prey'//trim(strindex)//'_max', 'g', 'maximum mass of prey '//trim(strindex))   
  
          
           call self%set_variable_property(self%id_prey_c(iprey), 'min_particle_mass', w_prey_min)
           call self%set_variable_property(self%id_prey_c(iprey), 'max_particle_mass', w_prey_max) 
 !      call depth_averaged_prey%set_variable_property(depth_averaged_prey%id_c, 'min_particle_mass', w_prey_min)
 !      call depth_averaged_prey%set_variable_property(depth_averaged_prey%id_c, 'max_particle_mass', w_prey_max)          
      !     call demersal_size_spectrum%parameters%set_real('w_source'//trim(strindex)//'_min', w_prey_min)
       !    call demersal_size_spectrum%parameters%set_real('w_source'//trim(strindex)//'_max', w_prey_max)
       !    call demersal_size_spectrum%couplings%set_string('source'//trim(strindex), '../benprey_c'//trim(strindex2))      

           call self%register_diagnostic_variable(self%id_bprey_c(b),'bprey'//trim(strindex)//'c','mg C/m^2/d',   'uptake of carbon in food source '//trim(strindex),    source=source_do_bottom)
           call self%register_diagnostic_variable(self%id_bprey_n(b),'bprey'//trim(strindex)//'n','mmol C/m^2/d',   'uptake of nitrogen in food source '//trim(strindex),    source=source_do_bottom)
           call self%register_diagnostic_variable(self%id_bprey_p(b),'bprey'//trim(strindex)//'p','mmol C/m^2/d',   'uptake of phosphorus in food source '//trim(strindex),    source=source_do_bottom)
           call self%register_diagnostic_variable(self%id_bprey_s(b),'bprey'//trim(strindex)//'s','mmol C/m^2/d',   'uptake of silicon in food source '//trim(strindex),    source=source_do_bottom)
            self%pelspec(iprey)=0.0_rk             


     !     call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_bprey_c(b))
      !    call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_bprey_n(b))
      !    call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_bprey_p(b))
      !    call self%add_to_aggregate_variable(standard_variables%total_silicate, self%id_bprey_s(b))
      else
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
            
            self%benspec(iprey)=0.0_rk      
      end if
     end if
   end do
   call self%add_child(pelagic_size_spectrum, 'pelagic_size_spectrum', configunit=-1)
  ! call self%add_child(demersal_size_spectrum, 'demersal_size_spectrum', configunit=-1)

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
   if (self%SRR == 1 .or. self%SRR == 2)  then 
       allocate(self%id_reproduction(self%nclass))
       !allocate(self%id_reproduction_ben(self%nclass))
   end if
   
   
   allocate(self%id_f_pel(self%nclass))
   allocate(self%id_f_ben(self%nclass))
   
   allocate(self%id_g_pel(self%nclass))
   allocate(self%id_g_ben(self%nclass))
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
      
!      allocate(depth_averaged_class_demersal)
!      depth_averaged_class_demersal%qnc = self%qnc
!      depth_averaged_class_demersal%qpc = self%qpc
!      call self%add_child(depth_averaged_class_demersal, 'class'//trim(strindex), configunit=-1)
!      call depth_averaged_class_demersal%request_coupling('int_c', '../c'//trim(strindex))
!      call depth_averaged_class_demersal%set_variable_property(depth_averaged_class_demersal%id_c,'particle_mass',self%w(iclass))

      ! If population is cannibalistic, add this size class as one of the prey (after the user-specified prey set).
      if (cannibalism) then
         write (strindex2,'(i0)') self%nprey - self%nclass + iclass
         call self%couplings%set_string('ave_prey'//trim(strindex2), './class'//trim(strindex))
      end if

      ! Register size-class-specific diagnostics
      if (self%SRR == 1 .or. self%SRR == 2) call self%register_diagnostic_variable(self%id_reproduction(iclass),'reproduction'//trim(strindex),'g m-2 d-1','allocation to reproduction in size class '//trim(strindex),         source=source_do_bottom)
   !    if (self%SRR == 1 .or. self%SRR == 2) call self%register_diagnostic_variable(self%id_reproduction_ben(iclass),'reproduction_ben'//trim(strindex),'g m-2 d-1','allocation to reproduction in size class '//trim(strindex),         source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_f_pel(iclass),           'f_pel'//trim(strindex),           '-',        'functional response of size class '//trim(strindex),                source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_f_ben(iclass),           'f_ben'//trim(strindex),           '-',        'functional response of size class '//trim(strindex),                source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_g_pel(iclass),           'g_pel'//trim(strindex),           'd-1',      'specific growth rate of individuals in size class '//trim(strindex),source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_g_ben(iclass),           'g_ben'//trim(strindex),           'd-1',      'specific growth rate of individuals in size class '//trim(strindex),source=source_do_bottom)
   end do

   allocate(self%phi(self%nprey,self%nclass))


   ! Register a state variable for waste (faeces, maintenance, dead matter resulting from non-predation mortality, fraction of offspring that does not survive)
   call register_waste('egested_matter', 'cnps', self%id_waste_c, self%id_waste_n, self%id_waste_p, self%id_waste_s)
   call register_waste('oxygen_consumed', 'o', id_o2=self%id_o2)
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
   
   !Add in links so can add the demersal pool
   call self%register_state_dependency(self%id_Bendip,'dipp','mmol P/m^2','disssolved inorganic phosphorus')
   call self%register_state_dependency(self%id_Bendin,'dinn','mmol N/m^2','dissolved inorganic nitrogen')
   call self%register_state_dependency(self%id_Bendic,'O3c','mmol C/m^2','carbon dioxide')
   call self%register_state_dependency(self%id_Beno2,'O2o','mmol O_2/m^2','oxygen')
   call self%register_state_dependency(self%id_bendiscard_c, 'bendiscard_c', 'mmol C m-2', 'organic carbon discards')
   call self%register_state_dependency(self%id_bendiscard_n, 'bendiscard_n', 'mmol N m-2', 'organic nitrogen discards')
   call self%register_state_dependency(self%id_bendiscard_p, 'bendiscard_p', 'mmol P m-2', 'organic phosphorus discards')
   call self%register_state_dependency(self%id_bendiscard_s, 'bendiscard_s', 'mmol S m-2', 'silicate discards')
   
   
   call self%request_coupling_to_model(self%id_Bendip,   'Benthic_excreted_phosphorus',   standard_variables%total_phosphorus)
   call self%request_coupling_to_model(self%id_Bendin,   'Benthic_excreted_nitrogen',   standard_variables%total_nitrogen)
   call self%request_coupling_to_model(self%id_Bendic,   'Benthic_respired_carbon',   standard_variables%total_carbon)
   call self%request_coupling_to_model(self%id_Beno2,   'Benthic_oxygen_consumed', 'o')
   call self%request_coupling_to_model(self%id_bendiscard_c, 'bendiscards', standard_variables%total_carbon)
   call self%request_coupling_to_model(self%id_bendiscard_n, 'bendiscards', standard_variables%total_nitrogen)
   call self%request_coupling_to_model(self%id_bendiscard_p, 'bendiscards', standard_variables%total_phosphorus)
   call self%request_coupling_to_model(self%id_bendiscard_s, 'bendiscards', standard_variables%total_silicate)
   
   call self%register_diagnostic_variable(self%id_fish_benDIP,'fish_benDIP','mmol m-2 d-1','benthic fish excretion to DIP',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fish_benDIN,'fish_benDIN','mmol m-2 d-1','benthic fish excretion to DIN',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fish_benDIC,'fish_benDIC','mmol m-2 d-1','benthic fish production of DIC from respiration',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_benO2_fish,'benO2_fish','mmol m-2 d-1','benthic fish consumption of O2 during repisration',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fish_benPOC,'fish_benPOC','mmol m-2 d-1','benthic fish egestion of POC',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fish_benPON,'fish_benPON','mmol m-2 d-1','benthic fish egestion of PON',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fish_benPOP,'fish_benPOP','mmol m-2 d-1','benthic fish egestion of POP',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_fish_benPOS,'fish_benPOS','mmol m-2 d-1','benthic fish egestion of POS',source=source_do_bottom)

   call self%register_state_variable(self%id_landings, 'landings', 'g m-2', 'landed biomass')
   call self%add_to_aggregate_variable(standard_variables%total_carbon, self%id_landings, scale_factor=1.0_rk/g_per_mmol_carbon)
   call self%add_to_aggregate_variable(standard_variables%total_nitrogen, self%id_landings, scale_factor=self%qnc/g_per_mmol_carbon)
   call self%add_to_aggregate_variable(standard_variables%total_phosphorus, self%id_landings, scale_factor=self%qpc/g_per_mmol_carbon)

   ! Register diagnostic for total offspring production across population.
   if (self%SRR == 1 .or. self%SRR == 2) then
      call self%register_diagnostic_variable(self%id_total_reproduction,'total_reproduction','mmol C m-2 d-1','total pelagic reproduction',source=source_do_bottom)
     ! call self%register_diagnostic_variable(self%id_total_reproduction_ben,'total_reproduction_ben','mmol C m-2 d-1','total benthic reproduction',source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_R_p,'R_p','# m-2 d-1','density-independent recruitment',source=source_do_bottom)
      !call self%register_diagnostic_variable(self%id_R_p_ben,'R_p_ben','# m-2 d-1','density-independent recruitment',source=source_do_bottom)
   elseif (self%SRR == 3) then
      ! Infer biomass of lowest size class by extending spectrum of (small) prey
      call self%register_dependency(self%id_offset, 'prey_spectrum_offset', '-', 'offset of pelagic prey spectrum')
      call self%register_dependency(self%id_slope, 'prey_spectrum_slope', '-', 'slope of pelagic prey spectrum')
      call self%request_coupling(self%id_offset, './pelagic_size_spectrum/offset')
      call self%request_coupling(self%id_slope, './pelagic_size_spectrum/slope')
     ! call self%register_dependency(self%id_benoffset, 'benthic_prey_spectrum_offset', '-', 'offset of benthic prey spectrum')
     ! call self%register_dependency(self%id_benslope, 'benthic_prey_spectrum_slope', '-', 'slope of benthic prey spectrum')
     ! call self%request_coupling(self%id_benoffset, './demersal_size_spectrum/benoffset')
     ! call self%request_coupling(self%id_benslope, './demersal_size_spectrum/benslope')
      
   end if
   call self%register_diagnostic_variable(self%id_R,'R','# m-2 d-1','recruitment',source=source_do_bottom)
  ! call self%register_diagnostic_variable(self%id_R_ben,'R_ben','# m-2 d-1','recruitment',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_tot, 'c_tot', 'g m-2', 'total biomass', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_size1, 'c_size1', 'g m-2', 'fish biomass smaller than threshold1', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_size2, 'c_size2', 'g m-2', 'fish biomass smaller than threshold2', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_size3, 'c_size3', 'g m-2', 'fish biomass smaller than threshold3', source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_c_lfi, 'c_lfi', 'g m-2', 'Large fish index', source=source_do_bottom)
   
   call self%register_diagnostic_variable(self%id_test, 'test', 'g m-2', 'test', source=source_do_bottom)
   contains
   
   subroutine register_waste(name, composition, id_c, id_n, id_p, id_s, id_o2)
      character(len=*), intent(in) :: name, composition
      type (type_bottom_state_variable_id), intent(inout), target, optional :: id_c, id_n, id_p, id_s, id_o2

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
      if (present(id_o2)) then
         call self%register_bottom_state_dependency(id_o2, name//'_o', 'mmol O_2 m-2','oxygen source')
         call self%request_coupling(id_o2, name//'_int/o2_int')
      end if
      call waste%request_coupling('w', '../total_pelprey_calculator/result')
      call waste%request_coupling('w_int', '../w_integrator/result')
      call waste%couplings%set_string('target', '../'//name)
   end subroutine

   end subroutine initialize
!EOC

   subroutine after_coupling(self)
      class (type_multi_element_demersal_pelagic_population),intent(inout) :: self

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
   !         if (self%prey(iprey)%ispel) then
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
            ! Compute size-class-specific preference for current prey;  Eq 4 in Hartvig et al. 2011 JTB, but note sigma typo confirmed by KH Andersen
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
      class (type_multi_element_demersal_pelagic_population),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iclass,iprey,istate
      real(rk) :: c_lfi, c_size1,c_size2, c_size3, slope, offset, endpoint, expected_eggs
      real(rk) :: total_reproduction,T_lim,T_lim_bot,temp,bot_temp,T_w_int,w_int,g_tot_pel, g_tot_ben
      real(rk) :: R,R_p
      real(rk) :: nflux_pel(0:self%nclass),nflux_ben(0:self%nclass),prey_state
      real(rk) :: E_e_pel,E_a_c_pel,E_a_n_pel,E_a_p_pel,E_a_s_pel,E_e_ben,E_a_c_ben,E_a_n_ben,E_a_p_ben,E_a_s_ben
      real(rk) :: f_pel,f_ben
      real(rk) :: ETW,b, omega,total_ben_prey
      real(rk) :: g_tot_c_pel,g_tot_n_pel,g_tot_p_pel,g_tot_c_ben,g_tot_n_ben,g_tot_p_ben
      real(rk),dimension(self%nprey)  :: prey_c,prey_n,prey_p,prey_s,prey_loss_pel, prey_loss_ben
      real(rk),dimension(self%nclass) :: Nw,I_c_pel, I_c_ben,I_n_pel,I_p_pel,I_s_pel,I_n_ben,I_p_ben,I_s_ben
      real(rk),dimension(self%nclass) :: mu_pel,reproduction,maintenance_pel,g_pel,maintenance_ben,g_ben,mu_ben
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
            _GET_(self%id_ETW,ETW)
            temp = T_w_int/w_int
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_T_ave, temp)
            T_lim = exp(self%c1-self%E_A/Boltzmann/(temp+Kelvin))
            bot_temp=ETW
            T_lim_bot=exp(self%c1-self%E_A/Boltzmann/(bot_temp+Kelvin))
         else
            T_lim = 1._rk
            T_lim_bot=1._rk
         end if

         !Determine fraction of time fish spend in pelagic based on ratio of pelagic to demersal food
         if (self%emergent_omega) then
             total_ben_prey=0._rk
             do iprey = 1, self%nprey
                  if (self%prey(iprey)%isben) then
                      total_ben_prey=total_ben_prey+prey_c(iprey)
                  end if
             end do
             omega=w_int/(w_int+total_ben_prey)
              
         else
             omega=self%omega
          end if
          _SET_HORIZONTAL_DIAGNOSTIC_(self%id_omega, omega)
          
         ! Food uptake (all size classes, all prey types)
         ! This computes total ingestion per size class (over all prey), and total loss per prey type (over all size classes)
         prey_loss_pel = 0.0_rk
         prey_loss_ben = 0.0_rk
         do iclass=1,self%nclass
            ! Compute total prey availability (concentration summed over all prey, scaled with prey-specific preference)
            E_a_c_pel = sum(self%phi(:,iclass)*prey_c*self%pelspec(:))
            E_a_n_pel = sum(self%phi(:,iclass)*prey_n*self%pelspec(:))
            E_a_p_pel = sum(self%phi(:,iclass)*prey_p*self%pelspec(:))
            E_a_s_pel = sum(self%phi(:,iclass)*prey_s*self%pelspec(:))
            E_a_c_ben = sum(self%phi(:,iclass)*prey_c*self%benspec(:))
            E_a_n_ben = sum(self%phi(:,iclass)*prey_n*self%benspec(:))
            E_a_p_ben = sum(self%phi(:,iclass)*prey_p*self%benspec(:))
            E_a_s_ben = sum(self%phi(:,iclass)*prey_s*self%benspec(:))
            
#ifndef NDEBUG
            if (isnan(E_a_c_pel)) &
               call self%fatal_error('do_bottom','E_a_c_pel is nan')
            if (isnan(E_a_c_ben)) &
               call self%fatal_error('do_bottom','E_a_c_ben is nan')
#endif

            ! Compute actual encounter in mmol prey s-1 (mmol C in predator)-1
            ! availability per volume (mmol prey m-3) times mass-specific volumetric search rate (m3 mmol predator-1 s-1)
           ! E_e = sum(self%V(:,iclass)*self%phi(:,iclass)*prey_c) ! Eq M3
             E_e_pel= self%V(iclass)*E_a_c_pel! Eq M3
             E_e_ben= self%VB(iclass)*E_a_c_ben! Eq M3
             
             
            ! Compute ingestion rate (s-1) - mmol C of prey per mmol C in predator per time!
            f_pel = E_e_pel/(E_e_pel+self%I_max(iclass))   ! Eq M5
            f_ben = E_e_ben/(E_e_ben+self%I_max(iclass))   ! Eq M5
                                   
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_pel(iclass),f_pel)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f_ben(iclass),f_ben)
                        
            I_c_pel(iclass) = T_lim*self%I_max(iclass)*f_pel ! ingestion part of M7
            I_c_ben(iclass) = T_lim_bot*self%I_max(iclass)*f_ben ! ingestion part of M7
#ifndef NDEBUG
            if (isnan(I_c_pel(iclass))) &
               call self%fatal_error('do_bottom','ingestion is nan')
               
            if (isnan(I_c_ben(iclass))) &
               call self%fatal_error('do_bottom','ingestion is nan')
#endif
            ! Predator-specific ingestion of nitrogen, phosphorus, silicate (mmol prey s-1 mmol-1)
            I_n_pel(iclass) = I_c_pel(iclass)/E_a_c_pel*E_a_n_pel
            I_p_pel(iclass) = I_c_pel(iclass)/E_a_c_pel*E_a_p_pel
            I_s_pel(iclass) = I_c_pel(iclass)/E_a_c_pel*E_a_s_pel
            
            I_n_ben(iclass) = I_c_ben(iclass)/E_a_c_ben*E_a_n_ben
            I_p_ben(iclass) = I_c_ben(iclass)/E_a_c_ben*E_a_p_ben
            I_s_ben(iclass) = I_c_ben(iclass)/E_a_c_ben*E_a_s_ben

            ! Account for this size class' ingestion in specific loss rate of all prey
            ! Units go from (mmol prey s-1 g-1) to (m s-1) - we divide by prey concentration and multiply by predator areal density
            !HP: Add array for whether the resource is pelagic or benthic. Stops benthic fish eating pelagic plankton and pelagic fish eating benthos. Loss still needs to be normalised to percentage of time that fish is in pelagic vs benthic systems
  
            
            prey_loss_pel(:) = prey_loss_pel(:) + I_c_pel(iclass)/E_a_c_pel*self%phi(:,iclass)*Nw(iclass)
            prey_loss_ben(:) = prey_loss_ben(:) + I_c_ben(iclass)/E_a_c_ben*self%phi(:,iclass)*Nw(iclass)
         end do

#ifndef NDEBUG
         if (any(prey_loss_pel<0)) &
            call self%fatal_error('do_bottom','prey_loss is negative')
         if (any(prey_loss_ben<0)) &
            call self%fatal_error('do_bottom','prey_loss is negative')
         if (any(prey_loss_pel>1/delta_t)) &
            call self%fatal_error('do_bottom','prey_loss is high (prey will be <0 within time step)')
         if (any(prey_loss_ben>1/delta_t)) &
            call self%fatal_error('do_bottom','prey_loss is high (prey will be <0 within time step)')
#endif

         ! Initialize size-class-specific mortality (s-1) with precomputed size-dependent background value.
         mu_pel = self%mu_b*T_lim + self%mu_s
         mu_ben = self%mu_b*T_lim_bot + self%mu_s
         
         ! Individual physiology (per size class)
         do iclass=1,self%nclass
            ! Specific maintenance rate (s-1)
            maintenance_pel(iclass) = T_lim*self%std_metab(iclass)
            maintenance_ben(iclass) = T_lim_bot*self%std_metab(iclass)

            ! Net specific element availability (mmol prey s-1 mmol-1)
            g_tot_c_pel = self%alpha*I_c_pel(iclass) - maintenance_pel(iclass)
            g_tot_c_ben = self%alpha*I_c_ben(iclass) - maintenance_ben(iclass)
            g_tot_n_pel = self%alpha*I_n_pel(iclass)
            g_tot_p_pel = self%alpha*I_p_pel(iclass)
            g_tot_n_ben = self%alpha*I_n_ben(iclass)
            g_tot_p_ben = self%alpha*I_p_ben(iclass)

            ! Specific growth rate (s-1) is minimum supported by different resources.
            g_tot_pel = min(g_tot_c_pel, g_tot_n_pel/self%qnc, g_tot_p_pel/self%qpc)
            
            ! Specific growth rate (s-1) is minimum supported by different resources.
            g_tot_ben = min(g_tot_c_ben, g_tot_n_ben/self%qnc, g_tot_p_ben/self%qpc)

            ! Avoid shrinking: limit maintenance to maximum sustainable value and increase starvation mortality.
            maintenance_pel(iclass) = min(maintenance_pel(iclass),self%alpha*I_c_pel(iclass))
            maintenance_ben(iclass) = min(maintenance_ben(iclass),self%alpha*I_c_ben(iclass))
            mu_pel(iclass) = mu_pel(iclass) + max(0.0_rk,-g_tot_pel/self%w(iclass)/self%xi)
            mu_ben(iclass) = mu_ben(iclass) + max(0.0_rk,-g_tot_ben/self%w(iclass)/self%xi)
            g_tot_pel = max(0.0_rk,g_tot_pel)
            g_tot_ben = max(0.0_rk,g_tot_ben)

            ! Individual growth (s-1)
            g_pel(iclass) = (1-self%psi(iclass))*g_tot_pel ! Eq M7
            g_ben(iclass) = (1-self%psi(iclass))*g_tot_ben ! Eq M7

            ! Mass flux towards reproduction (mmol C m-2 s-1) - sum over all individuals in this size class
            reproduction(iclass) = (omega*self%psi(iclass)*g_tot_pel*Nw(iclass)) + (1._rk-omega)*self%psi(iclass)*g_tot_ben*Nw(iclass)
         end do





         ! Compute number of individuals moving from each size class to the next (units: # s-1)
         nflux_pel(1:self%nclass) = Nw*g_pel/(self%delta_w/g_per_mmol_carbon)
         nflux_ben(1:self%nclass) = Nw*g_ben/(self%delta_w/g_per_mmol_carbon)

         ! Sum reproductive output of entire population in # m-2 s-1 (Eq 10 of Hartvig et al. 2011 JTB)
         ! Note: division by 2 is the result of the fact that reproductive output applies to females only,
         ! which are assumed to be 50% of the population.
         total_reproduction = sum(reproduction)
        ! total_reproduction_ben = sum(reproduction_ben)
         R_p = self%erepro/2 * total_reproduction / (self%w_min/g_per_mmol_carbon)
       !  R_p_ben = self%erepro/2 * total_reproduction_ben / (self%w_min_ben/g_per_mmol_carbon)

         ! Use stock-recruitment relationship to translate density-independent recruitment into actual recruitment (units: # s-1)
         if (self%SRR==0) then
            ! Constant recruitment
            R = self%recruitment
           ! R_ben = self%recruitment_ben
         elseif (self%SRR==1) then
            ! Density-independent recruitment
            R = R_p
          !  R_ben = R_p_ben
         elseif (self%SRR==2) then
            ! Beverton-Holt recruitment
            R = self%R_max*R_p/(R_p + self%R_max)
            !R_ben = self%R_max*R_p_ben/(R_p_ben + self%R_max)
         else
            _GET_HORIZONTAL_(self%id_offset, offset)
            _GET_HORIZONTAL_(self%id_slope, slope)
          !  _GET_HORIZONTAL_(self%id_benoffset, benoffset)
          !  _GET_HORIZONTAL_(self%id_benslope, benslope)
            endpoint = offset + slope * log(self%w_min)
          !  endpoint_ben= benoffset + benslope * log(self%w_min_ben)
            expected_eggs= exp(endpoint) * self%delta_w(1)
          !  expected_eggs_ben= exp(endpoint_ben) * self%delta_w(1)
            R = max(expected_eggs - Nw(1), 0._rk)/(self%w_min/g_per_mmol_carbon) * self%R_relax
            
            !HP:need to think about whether want to use first size class for benthic fish
           ! R_ben = max(expected_eggs_ben - Nw(1), 0._rk)/(self%w_min_ben/g_per_mmol_carbon) * self%R_relax
         end if

         ! Use recruitment as number of incoming individuals for the first size class.
         !HP: need to figure put how to get recruitment so it is not impacted by omega!
         nflux_pel(0) = R
      !   nflux_ben(0) = R

         ! Destroy all prey constituents that we are aware of (we only need to destroy carbon
         ! to get the correct impact on prey, but by destroying all we enable conservation checks)
         b=0._rk      
                  
         
         do iprey = 1, self%nprey
            if (iprey > self%nprey - self%nclass) then
               ! Prey is one of our own size classes
               _SET_BOTTOM_ODE_(self%id_prey_c(iprey), -omega*(prey_loss_pel(iprey)*prey_c(iprey)) -(1._rk-omega)*prey_loss_ben(iprey)*prey_c(iprey))
               
               
               !HP: NEED TO CHANGE
               !_SET_BOTTOM_ODE_(self%id_prey_c(iprey), -(1._rk-self%omega)*prey_loss_ben(iprey)*prey_c(iprey))
            elseif (self%feedback) then
               ! Prey is an external variable (not one of our size classes)
                   if (self%prey(iprey)%isben) then
                        b=b+1._rk
                  !     do istate=1,size(self%id_prey(iprey)%state)
                   !      _GET_(self%id_prey(iprey)%state(istate),preyP) 
                   !      _SET_ODE_(self%id_prey(iprey)%state(istate),-(1._rk-self%omega)*prey_loss_ben(iprey)*preyP)
                   !    end do
                   ! HP: only set ODE for carbon because that is only variable associated with YX (and Hx) - P and N are also automatically changed due to it having a constant stoichometry
                       _SET_BOTTOM_ODE_(self%id_prey_c(iprey),  -(1._rk-omega)*prey_loss_ben(iprey)*prey_c(iprey))
                    !   _SET_BOTTOM_ODE_(self%id_prey_n(iprey),  -(1._rk-self%omega)*prey_loss_ben(iprey)*prey_n(iprey))
                    !   _SET_BOTTOM_ODE_(self%id_prey_p(iprey),  -(1._rk-self%omega)*prey_loss_ben(iprey)*prey_p(iprey))
                    !  _SET_BOTTOM_ODE_(self%id_prey_s(iprey),  -(1._rk-self%omega)*prey_loss_ben(iprey)*prey_s(iprey))
                       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bprey_c(b),  -(1._rk-omega)*prey_loss_ben(iprey)*prey_c(iprey)*86400*CMass) 
                       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bprey_n(b),  -(1._rk-omega)*prey_loss_ben(iprey)*prey_n(iprey)*86400)
                       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bprey_p(b),  -(1._rk-omega)*prey_loss_ben(iprey)*prey_p(iprey)*86400)
                       _SET_HORIZONTAL_DIAGNOSTIC_(self%id_bprey_s(b),  -(1._rk-omega)*prey_loss_ben(iprey)*prey_s(iprey)*86400)
                   else
                       _SET_BOTTOM_ODE_(self%id_prey_c(iprey), -omega*prey_loss_pel(iprey)*prey_c(iprey))
                       _SET_BOTTOM_ODE_(self%id_prey_n(iprey), -omega*prey_loss_pel(iprey)*prey_n(iprey))
                       _SET_BOTTOM_ODE_(self%id_prey_p(iprey), -omega*prey_loss_pel(iprey)*prey_p(iprey))
                       _SET_BOTTOM_ODE_(self%id_prey_s(iprey), -omega*prey_loss_pel(iprey)*prey_s(iprey))                                      
                   end if
            end if
         end do
          !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_test,-sum((1._rk-self%omega)*prey_loss_ben*prey_c))
       !   PRINT*, self%omega
       !   PRINT*, omega
        !  _SET_HORIZONTAL_DIAGNOSTIC_(self%id_test,(omega))
         ! Transfer size-class-specific source terms and diagnostics to FABM
         do iclass=1,self%nclass
            ! Apply specific mortality (s-1) to size-class-specific abundances and apply upwind advection - this is a time-explicit version of Eq G.1 of Hartvig et al.
            _SET_BOTTOM_ODE_(self%id_c(iclass),(-(mu_pel(iclass) + self%F(iclass))*Nw(iclass) + (nflux_pel(iclass-1)-nflux_pel(iclass))*self%w(iclass)/g_per_mmol_carbon)*omega)
            _SET_BOTTOM_ODE_(self%id_c(iclass),(-(mu_ben(iclass) + self%F(iclass))*Nw(iclass) + (nflux_ben(iclass-1)-nflux_ben(iclass))*self%w(iclass)/g_per_mmol_carbon)*(1._rk-omega))
            
            
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_g_pel(iclass),g_pel(iclass)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_g_ben(iclass),g_ben(iclass)*86400)
            if (self%SRR == 1 .or. self%SRR == 2) then
               _SET_HORIZONTAL_DIAGNOSTIC_(self%id_reproduction(iclass),reproduction(iclass)*86400)
               !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_reproduction_ben(iclass),reproduction_ben(iclass)*86400)
            end if
         end do

         if (self%SRR == 1 .or. self%SRR == 2) then
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_total_reproduction,total_reproduction*86400)
            !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_total_reproduction_ben,total_reproduction_ben*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_p,R_p*86400)
           ! _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_p_ben,R_p_ben*86400)
         end if
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R,R*86400)
        ! _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_ben,R_ben*86400)

         ! Compute waste fluxes: total ingestion plus mortality, minus mass used in growth, minus recruitment, plus growth over right edge of resolved size range.
         if (self%feedback) then
            _SET_BOTTOM_ODE_(self%id_o2,-omega*self%resp_o2C*(1-self%alpha-self%alpha_eg)*sum(I_c_pel*Nw))            
            _SET_BOTTOM_ODE_(self%id_dic,omega*(1-self%alpha-self%alpha_eg)*sum(I_c_pel*Nw))
            _SET_BOTTOM_ODE_(self%id_din,omega*(1-self%alpha-self%alpha_eg)*sum(I_n_pel*Nw))
            _SET_BOTTOM_ODE_(self%id_dip,omega*(1-self%alpha-self%alpha_eg)*sum(I_p_pel*Nw))
            _SET_BOTTOM_ODE_(self%id_waste_c,(sum(((self%alpha+self%alpha_eg)*I_c_pel +  mu_pel - g_pel/(1-self%psi)          )*Nw) - R*self%w_min/g_per_mmol_carbon          + nflux_pel(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon)*omega)
            _SET_BOTTOM_ODE_(self%id_waste_n,(sum(((self%alpha+self%alpha_eg)*I_n_pel + (mu_pel - g_pel/(1-self%psi))*self%qnc)*Nw) - R*self%w_min/g_per_mmol_carbon*self%qnc + nflux_pel(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qnc)*omega)
            _SET_BOTTOM_ODE_(self%id_waste_p,(sum(((self%alpha+self%alpha_eg)*I_p_pel + (mu_pel - g_pel/(1-self%psi))*self%qpc)*Nw) - R*self%w_min/g_per_mmol_carbon*self%qpc + nflux_pel(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qpc)*omega)
            _SET_BOTTOM_ODE_(self%id_waste_s,sum(I_s_pel*Nw)*omega)
            
            _SET_BOTTOM_EXCHANGE_(self%id_beno2,-(1._rk-omega)*self%resp_o2C*(1-self%alpha-self%alpha_eg)*sum(I_c_ben*Nw))            
            _SET_BOTTOM_EXCHANGE_(self%id_bendic,(1._rk-omega)*(1-self%alpha-self%alpha_eg)*sum(I_c_ben*Nw))
            _SET_BOTTOM_EXCHANGE_(self%id_bendin,(1._rk-omega)*(1-self%alpha-self%alpha_eg)*sum(I_n_ben*Nw))
            _SET_BOTTOM_EXCHANGE_(self%id_bendip,(1._rk-omega)*(1-self%alpha-self%alpha_eg)*sum(I_p_ben*Nw))
            _SET_BOTTOM_EXCHANGE_(self%id_bendiscard_c, (sum(((self%alpha+self%alpha_eg)*I_c_ben +  mu_ben - g_ben/(1-self%psi)          )*Nw)          + nflux_ben(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon)*(1._rk-omega))
            
            _SET_BOTTOM_EXCHANGE_(self%id_bendiscard_n,(sum(((self%alpha+self%alpha_eg)*I_n_ben + (mu_ben - g_ben/(1-self%psi))*self%qnc)*Nw) + nflux_ben(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qnc)*(1._rk-omega))
            
            _SET_BOTTOM_EXCHANGE_(self%id_bendiscard_p,(sum(((self%alpha+self%alpha_eg)*I_p_ben + (mu_ben - g_ben/(1-self%psi))*self%qpc)*Nw)  + nflux_ben(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qpc)*(1._rk-omega))
            
            _SET_BOTTOM_EXCHANGE_(self%id_bendiscard_s,sum(I_s_ben*Nw)*(1._rk-omega))
            
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benDIP,  (1._rk-omega)*(1-self%alpha-self%alpha_eg)*sum(I_p_ben*Nw)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_benO2_fish,  (1._rk-omega)*self%resp_o2C*(1-self%alpha-self%alpha_eg)*sum(I_c_ben*Nw)*86400) 
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benDIC,  (1._rk-omega)*(1-self%alpha-self%alpha_eg)*sum(I_c_ben*Nw)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benDIN, (1._rk-omega)*(1-self%alpha-self%alpha_eg)*sum(I_n_ben*Nw)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benPOC, (sum(((self%alpha+self%alpha_eg)*I_c_ben +  mu_ben - g_ben/(1-self%psi)          )*Nw)          + nflux_ben(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon)*(1._rk-omega)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benPON, (sum(((self%alpha+self%alpha_eg)*I_n_ben + (mu_ben - g_ben/(1-self%psi))*self%qnc)*Nw) + nflux_ben(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qnc)*(1._rk-omega)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benPOP,(sum(((self%alpha+self%alpha_eg)*I_p_ben + (mu_ben - g_ben/(1-self%psi))*self%qpc)*Nw)  + nflux_ben(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass))/g_per_mmol_carbon*self%qpc)*(1._rk-omega)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_fish_benPOS,sum(I_s_ben*Nw)*(1._rk-omega))
  
            
         end if
         _SET_BOTTOM_ODE_(self%id_landings,sum(self%F*Nw*g_per_mmol_carbon))
      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

!-----------------------------------------------------------------------

end module mizer_multi_element_demersal_pelagic_population

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2017
!-----------------------------------------------------------------------
