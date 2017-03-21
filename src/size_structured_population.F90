#include "fabm_driver.h"

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mizer
!
! !INTERFACE:
module mizer_size_structured_population
!
! !DESCRIPTION:
! Size-structured population for 2D fauna (e.g., depth-integrated fish or benthic fauna), based on:
!
! Blanchard, J. L., Andersen, K. H., Scott, F., Hintzen, N. T., Piet, G., & Jennings, S. (2014)
! Evaluating targets and trade-offs among fisheries and conservation objectives using a multispecies size spectrum model
! Journal of Applied Ecology, 51(3), 612–622
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

   implicit none

!  default: all is private.
   private
!
! !PUBLIC DERIVED TYPES:
   type,extends(type_base_model),public :: type_size_structured_population
      ! Variable identifiers
      type (type_bottom_state_variable_id),         allocatable :: id_Nw(:)              ! Total weight per size class (sum over all individuals)
      type (type_bottom_state_variable_id),         allocatable :: id_Nw_prey(:)         ! Total weight per prey (sum over all individuals)
      type (type_bottom_state_variable_id)                      :: id_waste              ! State variable that will serve as sink for all waste
      type (type_horizontal_diagnostic_variable_id)             :: id_total_reproduction ! Total reproduction
      type (type_horizontal_diagnostic_variable_id)             :: id_R_p                ! Density-independent recruitment
      type (type_horizontal_diagnostic_variable_id)             :: id_R                  ! Density-dependent recruitment
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_reproduction(:)    ! Reproduction per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_f(:)               ! Functional response per size class
      type (type_horizontal_diagnostic_variable_id),allocatable :: id_g(:)               ! Specific growth rate per size class
      type (type_dependency_id)                                 :: id_T                  ! Temperature

      ! Number of size classes and prey
      integer :: nclass
      integer :: nprey

      ! Size class characteristics
      real(rk),allocatable :: logw(:)     ! log mass
      real(rk),allocatable :: w(:)        ! mass
      real(rk),allocatable :: delta_w(:)  ! mass difference between consecutive size classes

      ! Size-class-independent parameters
      real(rk) :: w_min       ! egg weight
      real(rk) :: alpha       ! assimilation efficiency
      real(rk) :: beta        ! preferred predator:prey mass ratio
      real(rk) :: sigma       ! s.d. of lognormal prey size selection function
      real(rk) :: erepro      ! efficiency of egg production
      real(rk) :: xi          ! fraction of mass consiting of lipid reserve (which can fuel maintenance)
      integer  :: SRR         ! type of stock-recruitment relationship (SRR)
      real(rk) :: recruitment ! constant recruitment flux (SSR=0)
      real(rk) :: R_max       ! maximum recruitment flux (SSR=2)

      integer  :: T_dependence ! Type of temperature dependence (0: none, 1: Arrhenius)
      real(rk) :: c1           ! Reference constant in Arrhenius equation = E_a/k/(T_ref+Kelvin)
      real(rk) :: E_a          ! Activation energy (eV)

      ! Size-class-dependent parameters that will be precomputed during initialization
      real(rk), allocatable :: V(:)            ! volumetric search rate (Eq M2)
      real(rk), allocatable :: I_max(:)        ! maximum ingestion rate (Eq M4)
      real(rk), allocatable :: std_metab(:)    ! standard metabolism (k*w^p in Eq M7)
      real(rk), allocatable :: mu_b(:)         ! background mortality
      real(rk), allocatable :: F(:)            ! fishing mortality
      real(rk), allocatable :: psi(:)          ! allocation to reproduction
      real(rk), allocatable :: phi(:,:)        ! prey preference
   contains
      procedure :: initialize
      procedure :: do_bottom
      procedure :: after_coupling
   end type type_size_structured_population

   ! Standard variable ("total mass") used for mass conservation checking
   type (type_bulk_standard_variable),parameter :: total_mass = type_bulk_standard_variable(name='total_mass',units='g',aggregate_variable=.true.,conserved=.true.)

   real(rk) :: Kelvin = 273.15_rk      ! offset of Celsius temperature scale (K)
   real(rk) :: Boltzmann = 8.62e-5_rk  ! Boltzmann constant
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
   class (type_size_structured_population), intent(inout),target :: self
   integer,                                 intent(in )          :: configunit
!
! !LOCAL VARIABLES:
   integer            :: iclass, iprey
   logical            :: cannibalism
   real(rk)           :: delta_logw
   real(rk)           :: k_vb,n,q,p,w_mat,w_inf,gamma,h,ks,f0,z0
   real(rk)           :: z0pre,z0exp,w_s,z_s
   real(rk)           :: kappa,lambda
   real(rk)           :: T_ref
   real(rk)           :: S1,S2,F
   integer            :: z0_type
   character(len=10)  :: strindex
   real(rk),parameter :: pi = 4*atan(1.0_rk)
   real(rk),parameter :: sec_per_year = 86400*365.2425_rk
!EOP
!-----------------------------------------------------------------------
!BOC
   ! Read parameters
   ! All rate coefficients are converted from d-1 to s-1 ensure source terms are directly compatible with FABM.
   ! Default values taken from Table S5
   call self%get_parameter(self%nclass,'nclass', '-',    'number of size classes', default=100)
   call self%get_parameter(self%nprey, 'nprey',  '-',    'number of prey')
   call self%get_parameter(self%alpha, 'alpha',  '-',    'assimilation efficiency',            default=0.6_rk,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%erepro,'erepro', '-',    'reproductive efficiency',            default=1.0_rk,   minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%w_min, 'w_min',  'g',    'egg weight',                         default=0.001_rk, minimum=0.0_rk)
   call self%get_parameter(n,          'n',      '-',    'exponent of max. consumption',       default=2.0_rk/3.0_rk)
   call self%get_parameter(q,          'q',      '-',    'exponent of search volume',          default=0.8_rk)
   call self%get_parameter(p,          'p',      '-',    'exponent of standard metabolism',    default=0.7_rk)
   call self%get_parameter(z0_type,    'z0_type','',     'type of background mortality (0: constant, 1: allometric function of size)', default=0)
   call self%get_parameter(z0pre,      'z0pre',  'yr-1', 'pre-factor for background mortality',default=0.6_rk,   minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(z0exp,      'z0exp',  '-',    'exponent of background mortality',   default=n-1)
   call self%get_parameter(w_s,        'w_s',    'g',    'start weight for senescence mortality',default=0._rk, minimum=0.0_rk)
   call self%get_parameter(z_s,        'z_s',    '-',    'exponent for senescence mortality',default=0.3_rk, minimum=0.0_rk)
   call self%get_parameter(w_mat,      'w_mat',  'g',    'maturation weight', default=0.0_rk, minimum=0.0_rk)
   call self%get_parameter(w_inf,      'w_inf',  'g',    'asymptotic weight', default=1e3_rk, minimum=0.0_rk)
   call self%get_parameter(self%beta,  'beta',   '-',    'preferred predator:prey mass ratio', default=100.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%sigma, 'sigma',  '-',    'width of prey size preference (sd in ln weight units)', default=1.0_rk, minimum=0.0_rk)
   call self%get_parameter(self%xi,    'xi',     '-',    'fraction of weight consisting of lipid reserve', default=0.1_rk, minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(k_vb,       'k_vb',   'yr-1', 'von Bertalanffy growth rate', minimum=0.0_rk, default=0.0_rk, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(lambda,     'lambda', '-',        'exponent of background resource spectrum', default=2+q-n)
   call self%get_parameter(kappa,      'kappa',  'g^(lambda-1)','carrying capacity of background resource spectrum', default=1e11_rk, minimum=0.0_rk)
   call self%get_parameter(f0,         'f0',     '-',        'background feeding level', default=0.6_rk, minimum=0.0_rk, maximum=1.0_rk)
   call self%get_parameter(self%SRR,   'SRR',    '',     'stock-recruitment relationship (0: constant recruitment, 1: density-independent recruitment, 2: Beverton-Holt)')
   select case (self%SRR)
   case (0)
      call self%get_parameter(self%recruitment, 'recruitment', '# yr-1', 'constant recruitment flux', minimum=0.0_rk, default=kappa*self%w_min**(-lambda)*sec_per_year, scale_factor=1._rk/sec_per_year)
   case (2)
      call self%get_parameter(self%R_max, 'R_max','# yr-1','maximum recruitment flux', minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   end select
   call self%get_parameter(S1,         'S1',         '-',        'offset for fishing selectivity exponent',       default=0.0_rk, minimum=0.0_rk)
   call self%get_parameter(S2,         'S2',         'g-1',      'scale factor for fishing selectivity exponent', default=0.0_rk, minimum=0.0_rk)
   call self%get_parameter(F,          'F',          'yr-1',     'fishing effort',                                default=0.0_rk, minimum=0.0_rk, scale_factor=1._rk/sec_per_year)
   call self%get_parameter(cannibalism,'cannibalism','',         'whether to enable intraspecific predation', default=.true.)

   call self%get_parameter(self%T_dependence, 'T_dependence', '', 'temperature dependence (0: none, 1: Arrhenius)', default=0)
   select case (self%T_dependence)
   case (1)
      call self%get_parameter(self%E_a, 'E_a',   'eV',              'activation energy',     default=0.63_rk)
      call self%get_parameter(T_ref,    'T_ref', 'degrees_Celsius', 'reference temperature', default=13._rk)
      self%c1 = self%E_a/Boltzmann/(T_ref+Kelvin)
      call self%register_dependency(self%id_T, standard_variables%temperature)
   end select

   ! Pre-factor for maximum ingestion rate derived from von Bertalanffy growth rate and background feeding level as described in appendix B
   h = 3*k_vb*w_inf**(1.0_rk/3.0_rk)/self%alpha/f0
   call self%get_parameter(h, 'h', 'yr-1 g^(-n)', 'pre-factor for maximum ingestion rate', minimum=0.0_rk, default=h*sec_per_year, scale_factor=1._rk/sec_per_year)

   ! Pre-factor for volumetric search rate (Eq M2)
   gamma = f0*h*self%beta**(2-lambda)/((1-f0)*sqrt(2*pi)*kappa*self%sigma)
   !gamma = f0*h*self%beta**(2-lambda)/((1-f0)*sqrt(2*pi)*kappa*self%sigma*exp((lambda-2)**2 * self%sigma**2 / 2)) ! add exp term taken from actual R code
   call self%get_parameter(gamma, 'gamma', 'yr-1 g^(-q)', 'pre-factor for volumetric search rate', minimum=0.0_rk, default=gamma*sec_per_year, scale_factor=1._rk/sec_per_year)

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
   allocate(self%F(self%nclass))
   allocate(self%psi(self%nclass))
   allocate(self%V(self%nclass))
   self%V(:) = gamma*self%w**(q-1)      ! specific volumetric search rate (specific, hence the -1!)
   self%I_max(:) = h*self%w**(n-1)      ! specific maximum ingestion rate [s-1]; Eq M4, but specific, hence the -1!
   self%std_metab(:) = ks*self%w**(p-1) ! specific metabolism [s-1]; second term in Eq M7, but specific, hence the -1!
   select case (z0_type)
   case (0)
      self%mu_b(:) = z0                    ! background mortality [s-1]; Eq M11
   case (1)
      self%mu_b(:) = z0pre*self%w**z0exp
   end select
   if (w_s>0.0_rk) self%mu_b = self%mu_b + 0.2/sec_per_year*(self%w/w_s)**z_s  ! Blanchard et al. 10.1098/rstb.2012.0231 Table S1
   self%F(:) = F/(1+exp(S1-S2*self%w))  ! fishing mortality [s-1]; Eqs M13 and M14 combined
   if (w_mat==0.0_rk) then
      ! No explicit reproduction as in original Blanchard community size spectrum model. Recruitment will be constant.
      self%psi(:) = 0
   else
      self%psi(:) = (self%w/w_inf)**(1-n) / (1+(self%w/w_mat)**(-10)) ! allocation to reproduction [-]; Eqs M13 and M14 combined
   end if
   write (*,*) 'Specific search volume (yr-1):'
   write (*,*) '  @ weight = 1:',gamma*sec_per_year,'(gamma)'
   write (*,*) '  @ minimum weight:',self%V(1)*sec_per_year
   write (*,*) '  @ maximum weight:',self%V(self%nclass)*sec_per_year
   write (*,*) 'Specific ingestion rate (yr-1)'
   write (*,*) '  @ weight = 1:',h*sec_per_year,'(h)'
   write (*,*) '  @ minimum weight:',self%I_max(1)*sec_per_year
   write (*,*) '  @ maximum weight:',self%I_max(self%nclass)*sec_per_year
   write (*,*) 'Specific metabolism (yr-1):'
   write (*,*) '  @ weight = 1:',self%std_metab(1)*sec_per_year,'(ks)'
   write (*,*) '  @ minimum weight:',self%std_metab(1)*sec_per_year
   write (*,*) '  @ maximum weight:',self%std_metab(self%nclass)*sec_per_year
   write (*,*) 'Background mortality (yr-1):'
   write (*,*) '  @ weight = 1:',z0*sec_per_year,'(z0)'
   write (*,*) '  @ minimum weight:',self%mu_b(1)*sec_per_year
   write (*,*) '  @ maximum weight:',self%mu_b(self%nclass)*sec_per_year
   write (*,*) 'Fishing mortality at minimum size:',self%F(1)*sec_per_year,'yr-1'
   write (*,*) 'Fishing mortality at maximum size:',self%F(self%nclass)*sec_per_year,'yr-1'

   ! Register dependencies for all prey.
   ! If the population is cannibalistic, autoamtically add all our size classes to the set of prey types.
   if (cannibalism) self%nprey = self%nprey + self%nclass
   allocate(self%id_Nw_prey(self%nprey))
   do iprey=1,self%nprey
      write (strindex,'(i0)') iprey
      call self%register_bottom_state_dependency(self%id_Nw_prey(iprey),'Nw_prey'//trim(strindex),'g m-2','biomass of prey '//trim(strindex))
   end do

   ! Allocate size-class-specific identifiers for abundance state variable and diagnostics.
   allocate(self%id_Nw(self%nclass))
   allocate(self%id_reproduction(self%nclass))
   allocate(self%id_f(self%nclass))
   allocate(self%id_g(self%nclass))
   do iclass=1,self%nclass
      ! Postfix for size-class-specific variable names (an integer number)
      write (strindex,'(i0)') iclass

      ! Register state variable, store associated individual mass (used by predators, if any, to determine grazing preference).
      call self%register_state_variable(self%id_Nw(iclass),'Nw'//trim(strindex),'g m-2','biomass in size class '//trim(strindex), 1.0_rk, minimum=0.0_rk)
      call self%set_variable_property(self%id_Nw(iclass),'particle_mass',self%w(iclass))

      ! Register this size class' contribution to total mass in the system (for mass conservation checks)
      call self%add_to_aggregate_variable(total_mass,self%id_Nw(iclass))

      ! If population is cannibalistic, add this size class as one of the prey (after the user-specified prey set).
      if (cannibalism) call self%request_coupling(self%id_Nw_prey(self%nprey-self%nclass+iclass),'Nw'//trim(strindex))

      ! Register size-class-specific diagnostics
      call self%register_diagnostic_variable(self%id_reproduction(iclass),'reproduction'//trim(strindex),'g m-2 d-1','allocation to reproduction in size class '//trim(strindex),         source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_f(iclass),           'f'//trim(strindex),           '-',        'functional response of size class '//trim(strindex),                source=source_do_bottom)
      call self%register_diagnostic_variable(self%id_g(iclass),           'g'//trim(strindex),           'd-1',      'specific growth rate of individuals in size class '//trim(strindex),source=source_do_bottom)
   end do

   allocate(self%phi(self%nprey,self%nclass))

   ! Register a state variable for waste (faeces, maintenance, dead matter resulting from non-predation mortality, fraction of offspring that does not survive)
   call self%register_bottom_state_variable(self%id_waste,'waste','g m-2','target variable for faeces')
   call self%add_to_aggregate_variable(total_mass,self%id_waste)

   ! Register diagnostic for total offspring production across population.
   call self%register_diagnostic_variable(self%id_total_reproduction,'total_reproduction','g m-2 d-1','total reproduction',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_R_p,'R_p','# d-1','density-independent recruitment',source=source_do_bottom)
   call self%register_diagnostic_variable(self%id_R,'R','# d-1','recruitment',source=source_do_bottom)

   end subroutine initialize
!EOC

   subroutine after_coupling(self)
      class (type_size_structured_population),intent(inout) :: self

      integer           :: iprey,iclass
      character(len=10) :: strindex
      real(rk)          :: w_p

      ! Coupling with prey has completed.
      ! Now we can query all prey for their weight per individual. From that we precompute predator-prey preferences.

      do iprey=1,self%nprey
         ! First retrieve individual mass of prey (throw fatal error if not available)
         w_p = self%id_Nw_prey(iprey)%link%target%properties%get_real('particle_mass',-1._rk)
         if (w_p<0) then
            write (strindex,'(i0)') iprey
            call self%fatal_error('after_coupling','prey '//trim(strindex)//' does not have attribute "particle_mass" set.')
         end if

         ! Compute size-class-specific preference for current prey; Eq 4 in Hartvig et al. 2011 JTB, but note sigma typo confirmed by KH Andersen
         ! This is a log-normal distribution of prey mass, scaled such that at optimum prey mass (=predator mass/beta), the preference equals 1.
         ! sigma is the standard deviation in ln mass units.
         do iclass=1,self%nclass
            self%phi(iprey,iclass) = exp(-(log(w_p)-self%logw(iclass)+log(self%beta))**2/self%sigma**2/2)
         end do
      end do
      if (any(self%phi<0)) call self%fatal_error('after_coupling','one or more entries of phi are < 0')

   end subroutine after_coupling

   subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
      class (type_size_structured_population),intent(in) :: self
      _DECLARE_ARGUMENTS_DO_BOTTOM_

      integer :: iclass,iprey
      real(rk) :: E_e,E_a,f,total_reproduction,T_lim,temp,g_tot,R,R_p,nflux(0:self%nclass)
      real(rk),dimension(self%nprey)  :: Nw_prey,prey_loss
      real(rk),dimension(self%nclass) :: Nw,I,maintenance,g,mu,reproduction
      real(rk), parameter :: delta_t = 12._rk/86400

      _HORIZONTAL_LOOP_BEGIN_

         ! Retrieve size-class-specific abundances
         do iclass=1,self%nclass
            _GET_HORIZONTAL_(self%id_Nw(iclass),Nw(iclass))
         end do

         ! Retrieve prey abundances
         do iprey=1,self%nprey
            _GET_HORIZONTAL_(self%id_Nw_prey(iprey),Nw_prey(iprey))
         end do
         Nw_prey = max(Nw_prey,0._rk)

         ! Temperature limitation factor affecting all rates (not in Blanchard et al.)
         if (self%T_dependence==1) then
            _GET_(self%id_T, temp)
            T_lim = exp(self%c1-self%E_A/Boltzmann/(temp+Kelvin))
         else
            T_lim = 1
         end if

         ! Food uptake (all size classes, all prey types)
         ! This computes total ingestion per size class (over all prey), and total loss per prey type (over all size classes)
         prey_loss = 0.0_rk
         do iclass=1,self%nclass
            ! Compute total prey availability (mass summed over all prey, scaled with prey-specific preference)
            E_a = sum(self%phi(:,iclass)*Nw_prey)
#ifndef NDEBUG
            if (isnan(E_a)) &
               call self%fatal_error('do_bottom','E_a is nan')
#endif

            ! Compute actual encounter (availability per volume times volumetric search rate)
            E_e = self%V(iclass)*E_a ! Eq M3

            ! Compute ingestion rate (g prey s-1 g-1) - per predator biomass!
            f = E_e/(E_e+self%I_max(iclass))   ! Eq M5
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_f(iclass),f)
            I(iclass) = T_lim*self%I_max(iclass)*f ! ingestion part of M7
#ifndef NDEBUG
            if (isnan(I(iclass))) &
               call self%fatal_error('do_bottom','ingestion is nan')
#endif

            ! Account for this size class' ingestion in specific loss rate of all prey: from (g s-1 #-1) to (s-1)
            prey_loss(:) = prey_loss(:) + I(iclass)/E_a*self%phi(:,iclass)*Nw(iclass)
         end do

#ifndef NDEBUG
         if (any(prey_loss<0)) &
            call self%fatal_error('do_bottom','prey_loss is negative')
         if (any(prey_loss>delta_t)) &
            call self%fatal_error('do_bottom','prey_loss is high (prey will be <0 within time step)')
#endif

         ! Initialize size-class-specific mortality (s-1) with precomputed size-dependent background value.
         mu = self%mu_b + self%F

         ! Individual physiology (per size class)
         do iclass=1,self%nclass
            ! Specific maintenance rate (s-1)
            maintenance(iclass) = T_lim*self%std_metab(iclass)

            ! Net specific energy availability (s-1)
            g_tot = self%alpha*I(iclass) - maintenance(iclass)

            ! Avoid shrinking: limit maintenance to maximum sustainable value and increase starvation mortality.
            maintenance(iclass) = min(maintenance(iclass),self%alpha*I(iclass))
            mu(iclass) = mu(iclass) + max(0.0_rk,-g_tot/self%w(iclass)/self%xi)
            g_tot = max(0.0_rk,g_tot)

            ! Individual growth (s-1)
            g(iclass) = (1-self%psi(iclass))*g_tot ! Eq M7

            ! Mass flux towards reproduction (g s-1) - sum over all individuals in this size class
            reproduction(iclass) = self%psi(iclass)*g_tot*Nw(iclass)
         end do

         ! Compute number of individuals moving from each size class to the next (units: # s-1)
         nflux(1:self%nclass) = Nw*g/self%delta_w

         ! Sum reproductive output of entire population in g s-1 (Eq 10 of Hartvig et al. 2011 JTB)
         ! Note: division by 2 is the result of the fact that reproductive output applies to females only,
         ! which are assumed to be 50% of the population.
         total_reproduction = sum(reproduction)
         R_p = self%erepro/2*total_reproduction/self%w_min

         ! Use stock-recruitment relationship to translate density-independent recruitment into actual recruitment (units: # s-1)
         if (self%SRR==0) then
            ! Constant recruitment
            R = self%recruitment
         elseif (self%SRR==1) then
            ! Density-independent recruitment
            R = R_p
         else
            ! Beverton-Holt recruitment
            R = self%R_max*R_p/(R_p + self%R_max)
         end if

         ! Use recruitment as number of incoming individuals for the first size class.
         nflux(0) = R

         ! Send prey consumption to FABM (combined impact of all size classes, per prey type)
         do iprey=1,self%nprey
            _SET_BOTTOM_ODE_(self%id_Nw_prey(iprey),-prey_loss(iprey)*Nw_prey(iprey))
         end do

         ! Transfer size-class-specific source terms and diagnostics to FABM
         do iclass=1,self%nclass
            ! Apply specific mortality (s-1) to size-class-specific abundances and apply upwind advection - this is a time-explicit version of Eq G.1 of Hartvig et al.
            _SET_BOTTOM_ODE_(self%id_Nw(iclass),-mu(iclass)*Nw(iclass) + (nflux(iclass-1)-nflux(iclass))*self%w(iclass))

            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_g(iclass),g(iclass)*86400)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_reproduction(iclass),reproduction(iclass)*86400)
         end do

         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_total_reproduction,total_reproduction*86400)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R_p,R_p*86400)
         _SET_HORIZONTAL_DIAGNOSTIC_(self%id_R,R*86400)
         _SET_BOTTOM_ODE_(self%id_waste,sum(((1-self%alpha)*I + maintenance + mu)*Nw) + total_reproduction - R*self%w_min + nflux(self%nclass)*(self%w(self%nclass)+self%delta_w(self%nclass)))
      _HORIZONTAL_LOOP_END_

   end subroutine do_bottom

!-----------------------------------------------------------------------

end module mizer_size_structured_population

!-----------------------------------------------------------------------
! Copyright Jorn Bruggeman/PML 2015-2016
!-----------------------------------------------------------------------
