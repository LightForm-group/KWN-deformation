!--------------------------------------------------------------------------------------------------
!> @author Madeleine Bignon, University of Manchester
!> @author Pratheek Shanthraj, University of Manchester
! KWN precipitation model including the effect of deformation via excess vacancy concentration
!---------------------------------------------------------------------------------------------------
!Useful references:
! [1] Robson, J. D. (2020). Metallurgical and Materials Transactions A: Physical Metallurgy and Materials Science, 51(10), 5401–5413. https://doi.org/10.1007/s11661-020-05960-5
! [2] Deschamps, A., Livet, F., & Bréchet, Y. (1998). . Acta Materialia, 47(1), 281–292. https://doi.org/10.1016/S1359-6454(98)00293-6
! [3] Bignon, M., Shanthraj, P., & Robson, J. D. (2022). Acta Materialia, 234, 118036. https://doi.org/10.1016/J.ACTAMAT.2022.118036 
program KWN

implicit none
integer, parameter :: pReal = selected_real_kind(25) 

real(pReal), parameter :: &
     
    kB = 1.38064852e-23_pReal, &
    R  = 8.3145_pReal, &
	PI = 3.14159265359_pReal, & 
	na = 6.02214076e23_pReal

    type :: tParameters

    	integer :: &
      		kwn_nSteps              ! discretization in r-space
    	real(pReal) :: &
      		kwn_stepsize, &         ! discretization in r-space
      		kwn_step0               ! minimum radius
    	real(pReal) :: &
      		lattice_param, &        ! lattice parameter in meter
      		atomic_volume, &        ! atomic volume in meter^3
      		molar_volume, &         ! molar volume in m^3/mol
      		misfit_energy, &        ! normalized precipitate misfit energy in J/m^3
      		gamma_coherent, &       ! coherent precipitate surface energy in J/m^2
      		vacancy_generation, &   ! vacancy generation rate coefficient
      		vacancy_sink_spacing, & ! vacancy sink spacing
      		vacancy_energy, &       ! normalized vacancy formation energy (Q/kB) in 1/k   
      		vacancy_migration_energy, & ! solute migration energy in J/at
      		vacancy_diffusion0, &   ! vacancy diffusivity in m^2/s 
      		mean_radius_initial,&  ! average radius initial distribution in meters
	  		standard_deviation,& ! standard deviation initial distribution (log normal law assumed)
	  		volume_fraction_initial, &! initial total precipitate distribution
	  		rho_0, & !initial dislocation density
	  		rho_s, & !saturation dislocation density
	  		strain_rate, & ! strain rate in /s
	  		dislocation_arrangement, & ! constant related to the dislocation density in the vacancy annihilation term, cf [1]
	  		burgers, & !matrice burgers vector
	  		jog_formation_energy, & ! formation energy for jogs
	  		q_dislocation, & ! activation energy for diffusion at dislocation (pipe diffusion) in J/at - not used yet but to be updated
	  		solute_strength, & ! constant related to the solid solution hardening- cf [2]
	  		enthalpy, & ! formation enthalpy of precipitates - cf [3] 
	  		entropy !formation entropy of precipitates - cf [3] 
	  
		! the following variables are allocatable to allow for precipitates with multiple elements (only situations with 2 elements are used here)
		real(pReal), dimension(:), allocatable :: &
	  		c0_matrix, &            ! initial matrix solute composition in mol fraction : [Mg, Zn]
      		ceq_matrix, &           ! equilibrium matrix composition for a flat interface (infinite precipitate radius) in mol fraction : [Mg, Zn]
      		ceq_precipitate, &      ! equilibrium precipitate composition in mol fraction : [Mg, Zn]
      		diffusion0, &           ! solute diffusivity in m^2/s : [Mg, Zn] - in the present version of the code - the diffusion coefficient is taken as identical for both solute and equal to the slowest diffuser (Mg)
      		migration_energy        !  solute migration energy in J/at
       
         
    	real(pReal), dimension(:),   allocatable :: &
     		bins                    ! Bins for class sizes in KWN model
    	character(len=15), allocatable, dimension(:) :: &
      		output
    end type tParameters



    type :: tKwnpowerlawState
     	real(pReal), pointer, dimension(:,:) :: &
      		precipitate_density ! table with precipitate density number in each class size [/m^4]
    	real(pReal),  dimension(  :), allocatable :: &
      		c_vacancy, & ! concentration in excess vacancy
      		time ! time [s]
    end type tKwnpowerlawState

!corresponds to the name if the type in Damask
 type :: tKwnpowerlawMicrostructure
   real(pReal),                  dimension(:),   allocatable :: &
     total_precipitate_density, &                                                             
     avg_precipitate_radius, &
     precipitate_volume_frac
    
  real(pReal), dimension(:,:), allocatable :: c_matrix
 end type tKwnpowerlawMicrostructure

!--------------------------------------------------------------------------------------------------
! containers for parameters and state ! changed compared to damask
  type(tParameters) :: prm
  type(tKwnpowerlawState) :: dot, &
    					     stt
  type(tKwnpowerlawMicrostructure) :: dst
  integer :: &
    ph, i, &
    Nmembers
!     sizeState, sizeDotState, &
!     startIndex, endIndex
    
  	
  logical ::decrease_time_step !1 if the time step needs to be reduced                                                                                
  ! test Madeleine, en will be the time increment
  integer :: en, &
  			 bin, &
  			 k, &
  			 Ntime, &
  			 N_elements!number of different elements in the precipitate	 	
  			 !flux_equality, &
  			 !size_q !to calculate the Guinier radius, number of steps in the values of q
  			 
  integer, dimension(:), allocatable :: stoechiometry !precipitate stoechiometry 			 
  
  real(pReal), dimension(:,:), allocatable :: &
  			results, &
  			 normalized_distribution_function
  			 
  	 real(pReal), dimension(:), allocatable :: &		 
  			 growth_rate_array, &
  			 x_eq_interface 
   !      
  real(pReal) :: T !temperature in K  
  real(pReal) :: &
    deltaGv, &
    interface_energy, &
    radius_crit, &
    nucleation_site_density, &
    zeldovich_factor, &
    beta_star, &
    incubation_time, &
    nucleation_rate,&
    radiusL, radiusR, radiusC, &
    growth_rate, flux, &
    interface_c, &
    vacancy_generation, &
    vacancy_annihilation, &
    vacancy_sink_spacing, &
    vacancy_diffusion, &
!     debug_min, &
!     debug_max, &
    N_0=0.0 ,&!number by which the normalized precipitate distribution should be multiplied when there is an initial one	  
    max_growth_rate, &
    time_record, &
    time_record_step =0.5, &
    beta, & ! parameter for the evolution of the dislocation density
	flow_stress, &
	c_thermal_vacancy, &
	dislocation_density, &
	production_rate, &
	annihilation_rate, &
	k1_km, &
	k2_km,  &
	alpha, &
	mu, &
	c_j, &
	guinier_radius, &
	time_cooling, &
	strain, &
	shape_parameter, &
	radius_transition, &
	mean_particle_strength, &
	tau_kwn, &
	gr_mean_radius, &
	x_eq_interface_mean_radius
	

  real(pReal), dimension(:), allocatable ::   diffusion_coefficient  
  

   
   ! the 'temp' variables are to store the previous step and adapt the time step at each iteration
  real(pReal) :: & 
		temp_total_precipitate_density =0.0, &
		temp_avg_precipitate_radius= 0.0, &
        temp_precipitate_volume_frac=0.0, &
        temp_radius_crit=0.0, &
        temp_c_vacancy=0.0, &
        delta_volume_fraction =1.0, & !to measure the difference between one step and the next one
        delta_total_precipitation_density=1.0, &   !to measure the difference between one step and the next one     
        temp_diffusion_coefficient, &
        delta_max=0.0, & 
        temp_dislocation_density, &
        T_start, &
        sigma_r, &
        A, &
        Q_stress, &
        n
 	
  
        
  real(pReal), dimension(:), allocatable ::  temp_c_matrix, &
        precipitate_atomic_fraction, &
        temp_x_eq_interface, &
        temp_x_eq_matrix, &
        	q,	 &
		intensity 
        
    real(pReal), dimension(:), allocatable :: 	temp_precipitate_density ,&
    temp_dot_precipitate_density, &
    k1,k2,k3,k4
	  
    real(pReal) :: dt, dt_temp, dt_max, total_time, delta_num, T_max, heating_rate, T_temp, h
	
	real(pReal) :: temperature_array(6)
	real(pReal) :: strain_rate_array(8)
	integer :: i_temp, k_temp
	character*100 :: filename
	character*100 :: filesuffix
	
	INTEGER :: nvals = 0 ! Number of values read in
	INTEGER :: status ! I/O status
	
	! --------------------------
	! --------------------------	 
	! --------------------------
	! --------------------------	 
	 N_elements=2 ! number of solute elements in the precipitate
	
	  allocate(prm%migration_energy(N_elements), source=0.0_pReal)	
	  allocate(prm%diffusion0(N_elements), source=0.0_pReal)	
	  allocate(prm%c0_matrix(N_elements), source=0.0_pReal)
	  allocate(prm%ceq_matrix(N_elements), source=0.0_pReal)
	  allocate(prm%ceq_precipitate(N_elements), source=0.0_pReal)
	  allocate(stoechiometry(N_elements+1))
	  
      

      

      Nmembers=1
    
    
    
    !!!!!!!!!! 
	! Read the inputs from the input.dat fil
	!!!!!!!!!!!	  
    OPEN (UNIT=1, FILE='input.dat', STATUS='OLD', ACTION='READ', IOSTAT=status)	   
    print*, status
    READ(1,*, IOSTAT=status) prm%kwn_step0! starting bin radius  (m) 
    print*, status
    READ(1,*) prm%kwn_stepsize ! spacing between bins (m) 
    READ(1,*) prm%kwn_nSteps ! no. of radius bins   
    READ(1,*) T
    READ(1,*) prm%strain_rate   
    READ(1,*) prm%lattice_param  ! lattice parameter of the product phase (m)
    READ(1,*) prm%atomic_volume ! atomic volume of matrix
    READ(1,*) prm%molar_volume  ! molar volume of precipitating phase      
	READ(1,*) prm%misfit_energy   ! elastic strain energy of the precipitate J/m^3
	READ(1,*) prm%gamma_coherent ! interfacial energy  J/m^2

    

      !----------- allocate variables
            
      en=1
	  allocate(dst%total_precipitate_density(Nmembers), source=0.0_pReal)
      allocate(dst%avg_precipitate_radius   (Nmembers), source=0.0_pReal)
      allocate(dst%precipitate_volume_frac  (Nmembers), source=0.0_pReal)
      allocate(dst%c_matrix                 (N_elements,Nmembers), source=0.0_pReal)
      
      
      allocate(dot%precipitate_density(prm%kwn_nSteps+1,1), source=0.0_pReal)        
      allocate(stt%precipitate_density(prm%kwn_nSteps+1,1), source=0.0_pReal)
      allocate( normalized_distribution_function(prm%kwn_nSteps+1,1), source=0.0_pReal)
      allocate(growth_rate_array(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(x_eq_interface(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(temp_x_eq_interface(prm%kwn_nSteps+1), source=0.0_pReal)
       
     
      allocate(temp_precipitate_density(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(temp_dot_precipitate_density(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(k1(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(k2(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(k3(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(k4(prm%kwn_nSteps+1), source=0.0_pReal)
      allocate(stt%time (Nmembers), source=0.0_pReal)
      allocate(stt%c_vacancy (Nmembers), source=0.0_pReal)
      allocate(dot%c_vacancy (Nmembers), source=0.0_pReal)
	  
	  	  ! define the bin in linear scale to match the paper on Co-Cu
       allocate(prm%bins(0:prm%kwn_nSteps), source=0.0_pReal)     
      	
      	
      	allocate(temp_c_matrix(N_elements), source=0.0_pReal)	    
      	allocate(temp_x_eq_matrix(N_elements), source=0.0_pReal)	 
        allocate(results(1,7))  

	  ph=1 
    
!!!!!!!!!! 
! Read the inputs from the input.dat fil
!!!!!!!!!!!
	
	 READ(1,*) prm%diffusion0(1) ! m^2/s
	 READ(1,*) prm%diffusion0(2)! m^2/s
	 READ(1,*) prm%migration_energy(1) !J/mol
	 READ(1,*) prm%migration_energy(2) !J/mol
	 READ(1,*) prm%c0_matrix(1) !at fraction
	 READ(1,*) prm%c0_matrix(2) ! at fraction
	 READ(1,*) stoechiometry(1)
	 READ(1,*) stoechiometry(2)
	 READ(1,*) stoechiometry(3)
	 READ(1,*) prm%enthalpy ! [J/mol]
	 READ(1,*) prm%entropy ![J/mol/K]
	 READ(1,*) prm%vacancy_sink_spacing ![m]
	 READ(1,*) prm%vacancy_diffusion0 ![m^2/s]
	 READ(1,*) prm%jog_formation_energy ![ev]
	 READ(1,*) prm%q_dislocation ! [J/mol]
	 READ(1,*) prm%vacancy_energy ! [ev]
	 READ(1,*) prm%vacancy_migration_energy ![ev]
	 READ(1,*) prm%dislocation_arrangement ![]
	 READ(1,*) prm%vacancy_generation ![]
	 READ(1,*) prm%rho_0 ![/m^2]
	 READ(1,*) prm%rho_s ![/m^2]
	 READ(1,*) prm%burgers ! [m]
	 READ(1,*) prm%mean_radius_initial ![m]
	 READ(1,*) prm%volume_fraction_initial ![]
	 READ(1,*) shape_parameter
	 READ(1,*) total_time ![s]
	 READ(1,*) dt_max ![s]
	 READ(1,*) time_record_step ![s]
	 READ(1,*) sigma_r ![MPa] - sinepower law for stress
	 READ(1,*) A ![/s] - sinepower law for stress 
	 READ(1,*) Q_stress ![J/mol] - activation energy in flow stress law
	 READ(1,*) n !exponent in sinepower law for stress
	 
	CLOSE(1)
		

		!!! some conversions
	
	prm%migration_energy=prm%migration_energy/na ! convert in J/at
	prm%jog_formation_energy = prm%jog_formation_energy*1.602176634e-19  ! convert from ev to  J/at
	prm%q_dislocation=prm%q_dislocation/na ! convert in J/at
	prm%vacancy_energy=prm%vacancy_energy*1.602176634e-19  ! convert from ev to  J/at
	prm%vacancy_migration_energy=prm%vacancy_migration_energy*1.602176634e-19  ! convert from ev to  J/at
	  
	  	  
	prm%solute_strength=840.0e6
    radius_transition=3.3e-9 !Fribourg2011
	prm%ceq_precipitate =real(stoechiometry(1:2))/real(sum(stoechiometry))
      
	dt=0.001
	dt=min(dt_max, dt)




	write(filesuffix, '(I3,"C_strain_rate",ES9.3, ".txt")'), int(T)-273, prm%strain_rate

 	
 	
 	!! initialise the bins for the size distribution
 	kwnbins_init: do i = 0, prm%kwn_nSteps
      prm%bins(i) = 10.0_pReal**(real(i,pReal)*prm%kwn_stepsize + prm%kwn_step0)      
    enddo kwnbins_init   
            
        !initialize some outputs
        growth_rate_array=0.0*growth_rate_array        
        stt%precipitate_density=0.0*stt%precipitate_density
        dst%total_precipitate_density=0.0_pReal      	
        dst%avg_precipitate_radius(en)=prm%mean_radius_initial
        dst%precipitate_volume_frac(en)=prm%volume_fraction_initial
        stt%c_vacancy(en) =0.0_pReal
        dislocation_density=prm%rho_0
      
      
      print*, 'Initial distribution:'        
      print*, 'standard deviation', prm%standard_deviation
      print*, 'step size:', prm%kwn_stepsize
      print*, 'number of steps', prm%kwn_nSteps
      
       if (dst%avg_precipitate_radius(en)>0) then
      
       !!!!!!!!initialize precipitate density for an already existing distribution
         
     
     distribution_function : do i=0, prm%kwn_nSteps
      		!log normal distribution	    			
     	normalized_distribution_function(i,en)		=	1.0_pReal/sqrt(PI*2.0_pReal) &
      																/shape_parameter/prm%bins(i)	&
      																*exp(-1.0/2.0*(log(prm%bins(i)/prm%mean_radius_initial)+shape_parameter**2/2)**2/shape_parameter**2)					
    enddo distribution_function

	!initialize the number density distribution
		 !the normalized distribution function gives the shape of the distribution, it needs to be multiplied by N0 to match the initial precipitate fraction
    
        normalized_distribution_function(:,en)=normalized_distribution_function(:,en)/sum(normalized_distribution_function(:,en))
    N_0=0.0_pReal
   do bin=1,prm%kwn_nSteps
   		radiusL = prm%bins(bin-1)
      	radiusR = prm%bins(bin  )
      	radiusC=(radiusL+radiusR)/2
     	N_0=N_0+ normalized_distribution_function(bin,en)*radiusC**3*(radiusR-RadiusL)*4.0/3.0*PI    
   enddo
       	N_0=dst%precipitate_volume_frac(en)/N_0
	  	stt%precipitate_density(:,en)=normalized_distribution_function(:,en)*N_0
        dst%total_precipitate_density(en)=sum(stt%precipitate_density)*prm%kwn_stepsize 
        
        stt%precipitate_density(1,en)=0
     
  
         do bin=1,prm%kwn_nSteps       	
         	if (stt%precipitate_density(bin,en)<1.0e-50_pReal) then         	
         	stt%precipitate_density(bin,en)=0.0_pReal
         	endif          
         enddo
         dot%precipitate_density=stt%precipitate_density/dt
        !recompute the initial radius, volume fraction etc 
         dst%avg_precipitate_radius(en) =0.0_pReal
         dst%total_precipitate_density =0.0_pReal
         dst%precipitate_volume_frac(en)=0.0_pReal 
          
        
        
        do bin=1,prm%kwn_nSteps
   			radiusL = prm%bins(bin-1)
      		radiusR = prm%bins(bin  )
    
   	 	!update precipitate density  
   		 dst%total_precipitate_density = dst%total_precipitate_density &
   								   + stt%precipitate_density(bin,en) &
   								   *(radiusR - radiusL) 								   
   		!update average radius  								   
    	dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     + stt%precipitate_density(bin,en) &
                                     * (radiusR**2.0_pReal - radiusL**2.0_pReal)/2.0_pReal ! at that stage in m/m^3 
                                     
        !update volume fraction
   		 dst%precipitate_volume_frac(en) = dst%precipitate_volume_frac(en) &
                                      + 1.0_pReal/6.0_pReal*PI &
                                      * (radiusR+ radiusL)**3.0_pReal &
                                      * (radiusR - radiusL) &
                                      * stt%precipitate_density(bin,en) 
                                				      																   
		enddo 
		      if (dst%total_precipitate_density(en) > 0.0_pReal) then
     			 dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     / dst%total_precipitate_density(en)
	 		  endif
        endif  
      
     !Write the initial precipitate distribution in a textfile
     filename='initial_precipitation_distribution_'    
     filename=trim(filename)//trim(filesuffix) 
     open(1, file = filename,  ACTION="write", STATUS="replace") 
     write(1,*) ' # Bin [m], Precipitate density distribution [/m^4]'
     
     do bin=1,prm%kwn_nSteps
     	write(1, 901) prm%bins(bin), stt%precipitate_density(bin,en)/sum(stt%precipitate_density(:,en))
     enddo
     close(1)
     
     
    
    !initialize the concentration in the matrix
     dst%c_matrix(:,en) = (prm%c0_matrix(:) - dst%precipitate_volume_frac(en)*prm%ceq_precipitate(:))&
     						/(1.0-dst%precipitate_volume_frac(en))

     diffusion_coefficient = 	prm%diffusion0*exp(-(prm%migration_energy )/T/kb) 	 +2*(dislocation_density)*prm%atomic_volume/prm%burgers&
  	 						 *prm%diffusion0*exp(-(prm%q_dislocation )/T/kb)  ! include pipe diffusion
     temp_diffusion_coefficient=diffusion_coefficient(1)

   	 call get_equilibrium(T, stoechiometry(:), dst%c_matrix(:,en), prm%ceq_matrix(:), &
   	 						N_elements, diffusion_coefficient, dst%precipitate_volume_frac(en), prm%entropy, prm%enthalpy)
     	                     
   	 call  interface_composition(T,  N_elements, prm%kwn_nSteps, stoechiometry, prm%c0_matrix,prm%ceq_matrix, &
   	 							prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, prm%bins, prm%gamma_coherent, &
   	 							 R, x_eq_interface, diffusion_coefficient, dst%precipitate_volume_frac(en), prm%entropy, prm%enthalpy)
    
     call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c, x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,  diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit ) 
	 print*, 'radius crit:', radius_crit*1,0e9, ' nm'               
    


     
           filename='temperature_'
           filename=trim(filename)//trim(filesuffix)
      	   open(1,file = filename,  ACTION="write", STATUS="replace" ) 
      	   write(1,*) '# Time [s], Temperature [K]'
      	   close(1)
      	             
           filename='diffusion_coefficient_'
           filename=trim(filename)//trim(filesuffix)
           open(1, file = filename,  ACTION="write", STATUS="replace")    
           write(1,*) '# Time [s], Diffusion coefficient [m^2/s] '       
           close(1)
                      
           filename='vacancies_'
           filename=trim(filename)//trim(filesuffix)	
	  	   open(1, file = filename,  ACTION="write", STATUS="replace")  
	  	   write(1,*) '# Time [s], c_{ex}/c_{th}, total number of produced vacancies/c_{th}, total number of annihilated vacancies /c_{th}'
	  	   close(1)
	  	   		 
		   filename='strain_'
           filename=trim(filename)//trim(filesuffix)
      	   open(1,file= filename,  ACTION="write", STATUS="replace")
      	   write(1,*) '# Time [s], Strain []'
      	   close(1) 
		
       	   filename='dislocation_density_'
       	   filename=trim(filename)//trim(filesuffix)
       	   open(1, file = filename, ACTION="write", STATUS="replace")
       	   write(1,*) '# Time [s], dislocation density [/m^2]'
       	   close(1)
       	   

       	  
            
      ! Write all the input parameters in a file
         filename='KWN_parameters_' 
         filename=trim(filename)//trim(filesuffix)
         open(201,file= filename,  ACTION="write", STATUS="replace") 
        WRITE (201,*) ' '
      	WRITE (201,*) 'KWN parameters'
        WRITE (201,*) ' '
		WRITE (201,100) prm%gamma_coherent
		100 FORMAT ('Interfacial energy: ', F7.3, ' J/m^2') 
      	WRITE (201,200) prm%migration_energy(1)*na
		200 FORMAT ('Migration energy: ', E30.6, ' J/mol') 
		WRITE (201,300) prm%diffusion0(1)
		300 FORMAT ('D0: ', E30.6, ' m^2/s') 	
		WRITE(201, 202) prm%enthalpy
		202 FORMAT ('Formation enthalpy:',E30.6, ' J/mol' )
		WRITE(201, 203) prm%entropy
		203 FORMAT ('Formation entropy:',E30.6, ' J/mol/K' )
		WRITE (201,*) ' '
		WRITE (201,*) 'Initial distribution'		
		WRITE (201,302) prm%mean_radius_initial	 
		302 FORMAT ('Initial mean radius: ', E30.6, ' m')    	
     	WRITE (201,303) prm%volume_fraction_initial
     	303 FORMAT ('Initial volume fraction: ', F7.3)     	 	
     	WRITE (201,304) prm%standard_deviation
     	304 FORMAT ('Standard deviation: ', E30.6)    	 
     	WRITE (201,320) shape_parameter
     	320 FORMAT ('Shape parameter: ', F7.3)      
		WRITE (201,*) ' '
		WRITE (201,315) prm%kwn_step0 
		315 FORMAT ('Starting bin : ', E30.6, ' m') 
		WRITE (201,316) prm%kwn_stepsize
		316 FORMAT ('Bin width : ', E30.6, ' m') 
		WRITE (201,317) prm%kwn_nSteps
		317 FORMAT ('Number of steps : ', I4) 	
		WRITE (201,*) ' '
		WRITE (201,*) 'Vacancy model parameters '
        WRITE (201,307) prm%vacancy_energy/(1.602176634e-19)
		307 FORMAT ('Vacancy formation energy: ', F7.3, ' [eV]')
        WRITE (201,308) prm%vacancy_migration_energy/(1.602176634e-19)
		308 FORMAT ('Vacancy migration energy: ',F7.3, ' [eV]')
	    WRITE (201,309) prm%vacancy_diffusion0
		309 FORMAT ('Pre-factor diffusion vacancy: ', E30.6, ' [m^2/s]')
		WRITE (201,305) prm%vacancy_generation 
		305 FORMAT ('Mechanical vacancy production constant: ', F7.3, ' []')
        WRITE (201,306) prm%vacancy_sink_spacing
		306 FORMAT ('Vacancy sink spacing (grain size): ', E30.6, ' [m]')  
        WRITE (201,310) prm%jog_formation_energy/(1.602176634e-19)
		310 FORMAT ('Jog formation energy: ', F7.3, ' [eV]')
		WRITE (201,311) prm%dislocation_arrangement
		311 FORMAT ('Dislocation arrangement parameter: ', F7.3, ' []')
		WRITE (201,312) prm%rho_0
		312 FORMAT ('Initial dislocation density: ', E30.6, ' [/m^2]')
		WRITE (201,318) prm%rho_s
		318 FORMAT ('Saturation dislocation density: ', E30.6, ' [/m^2]')
		WRITE (201,*) ''
		WRITE (201,*) 'Deformation conditions'	 
		WRITE (201,313) prm%strain_rate
		313 FORMAT ('Strain rate: ', E30.6, ' [/m^2]')
		WRITE (201,314) T-273
		314 FORMAT ('Temperature: ', F7.3, ' [°C]')
       close(201)
     
     
     
      !the following are used to store the outputs from the previous iteration  
      
    	temp_total_precipitate_density =dst%total_precipitate_density(en)
		temp_avg_precipitate_radius= dst%avg_precipitate_radius(en)
        temp_precipitate_volume_frac=dst%precipitate_volume_frac(en)	
        temp_c_matrix(:)=dst%c_matrix(:,en)
        temp_precipitate_density=stt%precipitate_density(:,en)
        temp_dot_precipitate_density=dot%precipitate_density(:,en)
        temp_radius_crit=radius_crit
        temp_x_eq_matrix=prm%ceq_matrix
		temp_x_eq_interface=x_eq_interface
		temp_c_vacancy=stt%c_vacancy(en)
		temp_dislocation_density=prm%rho_0
        T_temp=T
      	stt%time(en)=0.0_pReal 
      	time_record=-dt
      	nucleation_rate=0.0_pReal
      	h=dt


		   filename='kinetics_data_'
           filename=trim(filename)//trim(filesuffix)

          open(1, file = filename,  ACTION="write", STATUS="replace")  
          write(1,*), '#Time, [s], Average Radius [nm], Total precipitate density [/micron^3], Volume fraction [], Concentration in the matrix [at %]'
     

       		      		
       		!record initial state
       		     				   	results(1,1)=stt%time(en)
 									results(1,2)=dst%avg_precipitate_radius(en)*1.0e9 !in nm
    								results(1,3)=dst%total_precipitate_density(en)*1.0e-18 !in microns/m3
    		
    								if (results(1,3)<1.0e-30_pReal) then
    									results(1,3)=0.0
    								endif
    								
    								results(1,4)=dst%precipitate_volume_frac(en)
    								if (results(1,4)<1.0e-30_pReal) then
    									results(1,4)=0.0
    								endif
   										results(1,7:8)=dst%c_matrix(:,en)
   										results(1,6)=nucleation_rate*1.0e-18
   									if (results(1,6)<1.0e-30_pReal) then
    									results(1,6)=0.0
    								endif
   			
        							results(1,5)=radius_crit*1.0e9
        	
        	
    	    	 					WRITE(1,14) (results(1,i), i=1,8)
        							14 FORMAT(F40.6,F40.6, 6E40.6)
        									close(1)
        							filename='precipitation_distribution_'
           				filename=trim(filename)//trim(filesuffix)
           				
           	
           
        							open(2, file = filename,  ACTION="write", STATUS="replace")          
       			 					WRITE(2,'(E40.15)'), stt%time(en), stt%precipitate_density(:,en)
        							close(2)
        								
        								
        			    filename='diffusion_coefficient_'
           				filename=trim(filename)//trim(filesuffix)
						open(1, file = filename,  ACTION="write", position="append")  
		 				write(1, 601), stt%time(en), diffusion_coefficient(1)
		 	
		 				close(1)	 							
		 				filename='vacancies_'
		 				filename=trim(filename)//trim(filesuffix)	
						open(1, file = filename,  ACTION="write", position="append") 
						write(1, 1001), stt%time(en), stt%c_vacancy(en)/c_thermal_vacancy, production_rate/c_thermal_vacancy, annihilation_rate/c_thermal_vacancy
						
						close(1)
						filename='dislocation_density_'
           				filename=trim(filename)//trim(filesuffix)
		 				open(1, file = filename,  ACTION="write", position="append")  
		 				write(1, 901), stt%time(en), dislocation_density
		 					
		 				close(1)	
		 								     							
        							               
       print*, 'Initial mean radius :', dst%avg_precipitate_radius  
       print*, 'Initial precipitate density :',dst%total_precipitate_density(en)
	   print*, 'Initial volume fraction :', dst%precipitate_volume_frac(en)	
	   print*, 'Bulk composition', prm%c0_matrix(:)
       print*, 'Initial matrix composition :', dst%c_matrix(:,en)
	   print*, 'Equilibrium composition, precipitate :', prm%ceq_precipitate
	   print*, 'Equilibrium composition matrix :', prm%ceq_matrix	
	   

stt%time(en)=0				
k=0  
	
if (1==1) then
	loop_time : do while  (stt%time(en).LE. total_time)
  	k=k+1 
     		 print*, "dt:", dt
    		 print*, 'Time:', stt%time(en)  
    		 print*, 'Temperature', T
    	     print*, 'Mean radius : ', dst%avg_precipitate_radius(en)*1e9, 'nm'		
	diffusion_coefficient = prm%diffusion0*exp(-(prm%migration_energy )/T/kb) 	
	mu=(27.0+(21.5-27.0)/(650.0-273.0)*(T-273.0))*1.0e9; !McLellan 1987 MPa

	if (prm%strain_rate>0.0) then	
	flow_stress=sigma_r*asinh(((prm%strain_rate/(A))*exp(Q_stress/8.314/T))**(1/n))
	c_thermal_vacancy=23.0*exp(-prm%vacancy_energy/kB/T)
	c_j=exp(-prm%jog_formation_energy/kB/T)
	strain=prm%strain_rate*stt%time(en)	
	!from Detemple 1995
	dislocation_density=prm%rho_s*(1- (sqrt(prm%rho_s)-sqrt(prm%rho_0))/sqrt(prm%rho_s)*exp(-1.0/2.0*86*(strain))  )**2	
	endif				     
          
    ! calculate production and annihilation rate of excess vacancies					 
	production_rate = prm%vacancy_generation*flow_stress*prm%atomic_volume/prm%vacancy_energy*prm%strain_rate &
	                + 0.5*c_j*prm%atomic_volume/4.0/prm%burgers**3*prm%strain_rate
	                
	annihilation_rate=prm%vacancy_diffusion0*exp(-prm%vacancy_migration_energy/kB/T) &
					*(dislocation_density/prm%dislocation_arrangement**2+1.0/prm%vacancy_sink_spacing**2)*stt%c_vacancy(en)    
				
  	dot%c_vacancy(en) =production_rate-annihilation_rate
  	stt%c_vacancy(en) =stt%c_vacancy(en)+dot%c_vacancy(en)*dt
	
  	diffusion_coefficient = prm%diffusion0*exp(-(prm%migration_energy )/T/kb)&
  	 						*(1.0+stt%c_vacancy(en)/c_thermal_vacancy  ) &!don't include pipe diffusion
  	 						+2*(dislocation_density)*prm%atomic_volume/prm%burgers&
  	 						 *prm%diffusion0*exp(-(prm%q_dislocation )/T/kb) 
  	
  	
	   
    ! calculate nucleation rate
    nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume

    zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kB/T) &
                     / 2.0_pReal/PI/radius_crit/radius_crit

    beta_star = 4.0_pReal*PI&
              * radius_crit*radius_crit/(prm%lattice_param**4.0) &
              *1/(sum(1/diffusion_coefficient(:)*1/dst%c_matrix(:,en) ))

    incubation_time =  2.0/PI/zeldovich_factor/zeldovich_factor/beta_star
         
    if (stt%time(en) > 0.0_pReal) then
      nucleation_rate = nucleation_site_density*zeldovich_factor*beta_star &
                      * exp( &
                             - 4.0_pReal*PI*prm%gamma_coherent*radius_crit*radius_crit/3.0_pReal/kB/T &
                             - incubation_time/stt%time(en) &
                            )                                
    else
      nucleation_rate = 0.0_pReal
    endif  
    
    dot%precipitate_density=0.0*dot%precipitate_density
    	
    call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c,&
    		 x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    		 stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate, &
    		  diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit ) 

   
    	dot%precipitate_density(0,en)=0.0_pReal
    	dot%precipitate_density(1,en)=0.0_pReal	
    	stt%precipitate_density(0,en)=0.0_pReal
    	stt%precipitate_density(1,en)=0.0_pReal		
	
  
  ! Runge Kutta 4th order 
  ! https://en.wikipedia.org/wiki/Runge–Kutta_methods
  h=dt 
  k1=dot%precipitate_density(:,en) 
  stt%time(en)=stt%time(en)+h/2.0
  stt%precipitate_density(:,en)=temp_precipitate_density+h/2.0*k1


   call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c,&
    		 x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    		 stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,&
    		   diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit ) 
	
    nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume
    zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kB/T) &
                     / 2.0_pReal/PI/radius_crit/radius_crit
    !print*, 'Zeldovich factor:',zeldovich_factor
    beta_star = 4.0_pReal*PI&
              * radius_crit*radius_crit/(prm%lattice_param**4.0) &
              !*diffusion_coefficient(1)
              *1/(sum(1/diffusion_coefficient(:)*1/dst%c_matrix(:,en) ))
   ! print*, 'Beta star:', beta_star


	!slight change compared to the version in Damask
    incubation_time =  2.0/PI/zeldovich_factor/zeldovich_factor/beta_star
    
    !incubation_time=0.0
    !print*, 'Incubation time:', incubation_time
          
    if (stt%time(en) > 0.0_pReal) then
      nucleation_rate = nucleation_site_density*zeldovich_factor*beta_star &
                      * exp( &
                             - 4.0_pReal*PI*prm%gamma_coherent*radius_crit*radius_crit/3.0_pReal/kB/T &
                             - incubation_time/stt%time(en) &
                            )              
                            
                           ! print*, 'gamma', prm%gamma_coherent
                          !  print*, 'critical radius', radius_crit
                            !print*, 'kb', kB
                  
    else
      nucleation_rate = 0.0_pReal
    endif  
    
    
    dot%precipitate_density=0.0*dot%precipitate_density		
    call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c,&
    		 x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    		 	stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,  diffusion_coefficient, &
    		 	 dst%c_matrix(:,en), growth_rate_array, radius_crit ) 
    	
	
    !print*, 'critical radius method 2', radius_crit
    	dot%precipitate_density(0,en)=0.0_pReal
    	dot%precipitate_density(1,en)=0.0_pReal	
    	stt%precipitate_density(0,en)=0.0_pReal
    	stt%precipitate_density(1,en)=0.0_pReal			

 	
 	k2=dot%precipitate_density(:,en)
  
   stt%precipitate_density(:,en)=temp_precipitate_density+h/2.0*k2


    dot%precipitate_density=0.0*dot%precipitate_density
    
    call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c, &
    	 x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    	 stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate, &
    	  diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit ) 
		
	
    
    
    	dot%precipitate_density(0,en)=0.0_pReal
    	dot%precipitate_density(1,en)=0.0_pReal	
    	stt%precipitate_density(0,en)=0.0_pReal
    	stt%precipitate_density(1,en)=0.0_pReal		
	
 	
 	
 	k3=dot%precipitate_density(:,en)
  
  
  	stt%precipitate_density(:,en)=temp_precipitate_density+h*k3
	stt%time(en)= stt%time(en) +h/2.0

           

    
    nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume
    zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kB/T) &
                     / 2.0_pReal/PI/radius_crit/radius_crit
    !print*, 'Zeldovich factor:',zeldovich_factor
    beta_star = 4.0_pReal*PI&
              * radius_crit*radius_crit/(prm%lattice_param**4.0) &
              !*diffusion_coefficient(1)
              *1/(sum(1/diffusion_coefficient(:)*1/dst%c_matrix(:,en) ))
    incubation_time =  2.0/PI/zeldovich_factor/zeldovich_factor/beta_star
    
    incubation_time=0.0
    !print*, 'Incubation time:', incubation_time
          
    if (stt%time(en) > 0.0_pReal) then
      nucleation_rate = nucleation_site_density*zeldovich_factor*beta_star &
                      * exp( &
                             - 4.0_pReal*PI*prm%gamma_coherent*radius_crit*radius_crit/3.0_pReal/kB/T &
                             - incubation_time/stt%time(en) &
                            )              
                            

                  
    else
      nucleation_rate = 0.0_pReal
    endif  

    dot%precipitate_density=0.0*dot%precipitate_density
     call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c, &
     		 x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
     		  stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,  &
     		  diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit ) 
		!print*, 'growth rate array', growth_rate_array
		

    dot%precipitate_density(0,en)=0.0_pReal
    dot%precipitate_density(1,en)=0.0_pReal	
    stt%precipitate_density(0,en)=0.0_pReal
    stt%precipitate_density(1,en)=0.0_pReal		
    	
    k4=dot%precipitate_density(:,en)
    
    stt%precipitate_density(:,en)=temp_precipitate_density + h/6.0*(k1+2.0*k2+2.0*k3+k4)      
    dst%precipitate_volume_frac(en) = 0.0_pReal
    dst%total_precipitate_density(en) = 0.0_pReal
    dst%avg_precipitate_radius(en) = 0.0_pReal

! update radius, total precipitate density and volume fraction

kwnbins: do bin=1,prm%kwn_nSteps
   	radiusL = prm%bins(bin-1)
      radiusR = prm%bins(bin  )
    
    !update precipitate density  
    dst%total_precipitate_density = dst%total_precipitate_density &
   								   + stt%precipitate_density(bin,en) &
   								   *(radiusR - radiusL) 								   
   	!update average radius  								   
    dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     + stt%precipitate_density(bin,en) &
                                     * (radiusR**2.0_pReal - radiusL**2.0_pReal)/2.0_pReal ! at that stage in m/m^3 
    !update volume fraction
    !small change compared to the expression in Damask
    dst%precipitate_volume_frac(en) = dst%precipitate_volume_frac(en) &
                                      + 1.0_pReal/6.0_pReal*PI &
                                      * (radiusR+ radiusL)**3.0_pReal &
                                      * (radiusR - radiusL) &
                                      * stt%precipitate_density(bin,en) 	
                                      
                                                                      					      																   
enddo kwnbins

      if (dst%total_precipitate_density(en) > 0.0_pReal) then
      dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     / dst%total_precipitate_density(en)
	 endif
		
	!update matrix composition
	print*, 'precipitate composition :', prm%ceq_precipitate(:)
	
	!print*, 'initial composition', prm%c0_matrix(:)
	!print*, 'precipitate composition', prm%ceq_precipitate(:)
	 dst%c_matrix(:,en) = (prm%c0_matrix(:) - dst%precipitate_volume_frac(en)*prm%ceq_precipitate(:))&
	 						/(1-dst%precipitate_volume_frac(en))
		
		

	
   ! print*, ''
    print*, 'Total precipitate density : ' , dst%total_precipitate_density*1e-18 , '/micron^3'
    print*, 'Precipitate volume fraction :',  dst%precipitate_volume_frac(en)
    print*, 'Solute concentration in the matrix' , dst%c_matrix(1,en)
	print*, 'Nucleation rate :/micron^3 ', nucleation_rate*1.0e-18
   ! print*, ''
 
   ! Adapt time step so that the outputs do not vary to much between too time steps
        
    !if  either:
    ! - the precipitate distribution in one class is negative,
    ! or - the concentration in any element is negative in the matrix
    ! or - the concentration changes very fast
    ! then go back to the previous step and decrease the time step
  
if (any(isnan(stt%precipitate_density(:,en))))then
print*, 'time :', stt%time(en)
endif 

    if  ( ((stt%c_vacancy(en)) <0.0_pReal)  .OR.(minval(stt%precipitate_density(:,en)) <0.0_pReal).OR.(  minval(dst%c_matrix(:,en))<0.0_pReal).OR.any(isnan(stt%precipitate_density(:,en))))  then	       
    	!reduce the time step
    	print*,'delta c vacancy', abs(temp_c_vacancy-stt%c_vacancy(en))/stt%c_vacancy(en)
    	stt%time(en)=stt%time(en)-dt
    	dst%total_precipitate_density= temp_total_precipitate_density 
		stt%precipitate_density(:,en)=temp_precipitate_density(:) 
		dot%precipitate_density(:,en)=temp_dot_precipitate_density(:) 
		dst%avg_precipitate_radius(en) = temp_avg_precipitate_radius 
        dst%precipitate_volume_frac(en) = temp_precipitate_volume_frac  
        dst%c_matrix(:,en)  = temp_c_matrix
        stt%c_vacancy(en) = temp_c_vacancy
        dislocation_density=temp_dislocation_density
        T=T_temp
        radius_crit=temp_radius_crit
        prm%ceq_matrix= temp_x_eq_matrix
		x_eq_interface = temp_x_eq_interface
		diffusion_coefficient(1)=temp_diffusion_coefficient
       !decrease factor arbitrary chosen, ask about that
        dt=0.5*dt


    else
    !this is to make sure that a precipitate cannot grow from more than the space between two adjacent classes
    
    	    	!not necessary to adapt the time step to the growth rate if there is no precipitate
    	    	if (dst%total_precipitate_density(en)>1.0_pReal) then
    	  
    	    	dt=min(dt_max,(prm%bins(1)-prm%bins(0))/maxval(abs(growth_rate_array)))
    	    	print*,'dt growth rate', (prm%bins(1)-prm%bins(0))/maxval(abs(growth_rate_array))
    	    	else
    	    	dt=min(dt_max,dt*1.2) !increase slightly the time step by an arbitrary factor as long as there are no precipitates
				endif
    	
    	temp_dot_precipitate_density=dot%precipitate_density(:,en)   	
    	temp_total_precipitate_density =dst%total_precipitate_density(en)
		temp_precipitate_density = stt%precipitate_density(:,en)
		temp_avg_precipitate_radius = dst%avg_precipitate_radius(en)
        temp_precipitate_volume_frac = dst%precipitate_volume_frac(en)
        temp_c_matrix = dst%c_matrix(:,en)
        temp_radius_crit=radius_crit
        temp_x_eq_matrix=prm%ceq_matrix
		temp_x_eq_interface=x_eq_interface
		temp_c_vacancy=stt%c_vacancy(en)
		temp_dislocation_density= dislocation_density
        T_temp=T
        temp_diffusion_coefficient=diffusion_coefficient(1)

!         	print*, 'precipitate distribution :', stt%precipitate_density(:,en)
            if (time_record<stt%time(en)) then !record this every 'time_record' seconds
       	
     				   results(1,1)=stt%time(en)
 					   results(1,2)=dst%avg_precipitate_radius(en)*1.0e9 
    				   results(1,3)=dst%total_precipitate_density(en)*1.0e-18
    		
    				   if (results(1,3)<1.0e-30_pReal) then
    									results(1,3)=0.0
    					endif
    								
    				  results(1,4)=dst%precipitate_volume_frac(en)
    				  if (results(1,4)<1.0e-30_pReal) then
    						results(1,4)=0.0
    				  endif
   					  results(1,7:8)=dst%c_matrix(:,en)
   					  results(1,6)=nucleation_rate*1.0e-18
   						if (results(1,6)<1.0e-30_pReal) then
    						results(1,6)=0.0
    					endif
        			  results(1,5)=radius_crit*1.0e9
        	
        			  filename='kinetics_data_'
           			  filename=trim(filename)//trim(filesuffix)

          			  open(1, file = filename,  ACTION="write", position="append")  	
    	    	 	  WRITE(1,13) (results(1,i), i=1,8)
        			  13 FORMAT(F40.6,F40.6,E40.6,E40.6,E40.6, E40.6, 2E40.6 )
        			  close(1)
        							
        			 ! writes the current distribution
        				filename='precipitation_distribution_'
           				filename=trim(filename)//trim(filesuffix)
           
        							open(2, file = filename,  ACTION="write", STATUS="replace")          
       			 					WRITE(2,'(E40.15)'), stt%time(en), stt%precipitate_density(:,en)
        							close(2)
        								
        								
        			    filename='diffusion_coefficient_'
           				filename=trim(filename)//trim(filesuffix)
						open(1, file = filename,  ACTION="write", position="append")  
		 				write(1, 601), stt%time(en), diffusion_coefficient(1)
		 				601 FORMAT(2E40.6)
		 				close(1)	 							
		 				filename='vacancies_'
		 				filename=trim(filename)//trim(filesuffix)	
						open(1, file = filename,  ACTION="write", position="append") 
						write(1, 1001), stt%time(en), stt%c_vacancy(en)/c_thermal_vacancy, production_rate/c_thermal_vacancy, annihilation_rate/c_thermal_vacancy
						1001 FORMAT(4E40.6)
						close(1)
						filename='dislocation_density_'
           				filename=trim(filename)//trim(filesuffix)
		 				open(1, file = filename,  ACTION="write", position="append")  
		 				write(1, 901), stt%time(en), dislocation_density
		 				901 FORMAT(3E40.6)			
		 				close(1)	
	       				
        				time_record=time_record+time_record_step
    	   	endif   	      	   	                   
 	endif           
end do loop_time   
endif





 end program KWN
 
 
subroutine get_equilibrium(T, stoechiometry, c_matrix, ceq_matrix, N_elements, diffusion_coefficient, volume_fraction, entropy, enthalpy)
	! get equilibrium for a flat interface
	implicit none
	integer, parameter:: pReal = selected_real_kind(25) 

	integer, intent(in):: N_elements
	integer, intent(in),  dimension(N_elements+1) :: stoechiometry
	real(pReal), intent(in), dimension(N_elements) ::  c_matrix, diffusion_coefficient
	real(pReal), intent(inout), dimension(N_elements) :: ceq_matrix
	real(pReal), intent(in) :: T, volume_fraction, entropy, enthalpy
	real(pReal) :: x_min, x_max, solubility_product, delta
	
	
	 if (N_elements>1) then   


	solubility_product=exp((entropy/8.314_pReal-enthalpy/8.314_pReal/T))  ! for 7449
		x_min=0.0_pReal
		x_max=1.0_pReal!

	!  find the intersection between stoichiometric line and solubility line by dichotomy
	do while (abs(x_min-x_max)/x_max>0.001_pReal.AND.x_max>1.0e-8_pReal)	

	ceq_matrix(1)=(x_min+x_max)/2.0_pReal
	ceq_matrix(2)=(c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))&
				*diffusion_coefficient(1)/diffusion_coefficient(2)&
				*(ceq_matrix(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction)
	if (ceq_matrix(2)<0.0_pReal) then
	
	x_min =ceq_matrix(1)
	
	else

	delta=ceq_matrix(1)**stoechiometry(1)*((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))&
		*diffusion_coefficient(1)/diffusion_coefficient(2)*(ceq_matrix(1)&
		*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2)&
			-solubility_product
	
	!print*, 'solubility product:', solubility_product
	delta=delta
	essai: if (delta > 0.0) then

		x_max=ceq_matrix(1)

		
	else 

		x_min=ceq_matrix(1)

	endif essai
	
	endif

	enddo
	!print*, 'c_matrix(2)', c_matrix(2)
	ceq_matrix(2)=(c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))&
				*diffusion_coefficient(1)/diffusion_coefficient(2)&
				*(ceq_matrix(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction)
	endif
	return
end subroutine


subroutine interface_composition(T,  N_elements, N_steps, stoechiometry, &
			c_matrix,ceq_matrix, atomic_volume, na, molar_volume, ceq_precipitate, &
			bins, gamma_coherent, R,  x_eq_interface, diffusion_coefficient, volume_fraction, entropy, enthalpy)
	implicit none
	integer, parameter :: pReal = selected_real_kind(25) 
	integer, intent(in) :: N_elements, N_steps  
	integer, intent(in), dimension(N_elements+1) :: stoechiometry
	real(pReal), intent(in), dimension(N_elements) :: c_matrix, ceq_precipitate, diffusion_coefficient, ceq_matrix
	real(pReal), intent(in) :: T,  atomic_volume, na, molar_volume, gamma_coherent, R, volume_fraction, entropy, enthalpy
	real(pReal), intent(inout), dimension(N_steps+1) :: x_eq_interface
    real(pReal), intent(in), dimension(N_steps+1) :: bins
	real(pReal) :: xmin, xmax, solubility_product, delta
	integer :: i
	 
		xmin=0.0_pReal
		xmax=1.0_pReal! no dissolution to see
	
   interface_equilibrium: do i = 0, N_steps
   if (N_elements>1) then   	
   		solubility_product=exp((entropy/8.314_pReal-enthalpy/8.314_pReal/T))
	
	
		xmin=ceq_matrix(1)
		xmax=atomic_volume*na/molar_volume*ceq_precipitate(1)! to be changed
		do while (abs(xmin-xmax)/xmax >0.0001.AND.xmax>1.0e-8)

		!print*,'x eq infinity', ceq_matrix(1)
		x_eq_interface(i)=(xmin+xmax)/2.0
	
		!  find the intersection between stoichiometric line and solubility line by dichotomy
		delta=x_eq_interface(i)**stoechiometry(1)*&
		((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))*diffusion_coefficient(1)/diffusion_coefficient(2)*&
		(x_eq_interface(i)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2)&
			-solubility_product*exp(2.0*molar_volume*gamma_coherent/R/T/bins(i)*real(sum(stoechiometry)) )
	
		
		if (delta<0.0_pReal) then
			xmin=x_eq_interface(i)
		else
			xmax=x_eq_interface(i)
		endif	
		enddo
	
	else 	
	x_eq_interface(i) =ceq_matrix(1)*exp(2.0*molar_volume*gamma_coherent/R/T/bins(i)/ceq_precipitate(1))
	
	endif

    enddo    interface_equilibrium
     			
return
end subroutine

subroutine growth_precipitate(N_elements, N_steps, bins, interface_c, &
x_eq_interface,atomic_volume, na, molar_volume, ceq_precipitate, precipitate_density, &
dot_precipitate_density, nucleation_rate, diffusion_coefficient, c_matrix, growth_rate_array , radius_crit) 
implicit none
	integer, parameter :: pReal = selected_real_kind(25) 
	integer, intent(in) :: N_Steps, N_elements
	real(pReal), intent(in), dimension(N_steps) :: bins
	real(pReal), intent(in), dimension(N_steps) :: x_eq_interface, precipitate_density
	real(pReal), intent(in), dimension(N_elements) :: ceq_precipitate, diffusion_coefficient, c_matrix
	real(pReal), intent(in) :: atomic_volume, na, molar_volume,  nucleation_rate
	real(pReal), intent(inout), dimension(N_steps) :: dot_precipitate_density(N_steps)
	real(pReal), intent(inout), dimension(N_steps-1) :: growth_rate_array
	real(pReal), intent(out)::  radius_crit
	real(pReal) :: radiusC, radiusL, radiusR, interface_c, growth_rate, flux
	integer :: bin
	!print*, 'routine nucleation rate:', nucleation_rate
    kwnbins_growth: do bin = 1, N_Steps-1
    
      radiusC = bins(bin  )
      radiusL = bins(bin-1)
      radiusR = bins(bin+1)
      
     interface_c=x_eq_interface(bin)
	  growth_rate = diffusion_coefficient(1)/radiusC &
                  * (c_matrix(1)    - interface_c) &
                  / (atomic_volume*na/molar_volume*ceq_precipitate(1) - interface_c)
       
     !the growth rate is stored to change the time step
     growth_rate_array(bin)=growth_rate

      if (growth_rate > 0.0_pReal) then
        flux = precipitate_density(bin)*growth_rate        
      else
        flux = precipitate_density(bin+1)*growth_rate     
      endif   
      
      dot_precipitate_density(bin  ) = dot_precipitate_density(bin ) - flux/(radiusC - radiusL)
      dot_precipitate_density(bin+1) = dot_precipitate_density(bin+1) + flux/(radiusR - radiusC)   
    enddo kwnbins_growth  
    
    radius_crit=bins(minloc(abs(growth_rate_array(1: N_Steps-1-2)),1))	
    
    nucleation: do bin = 1, N_Steps-1
    
      radiusC = bins(bin  )
      radiusL = bins(bin-1)
      radiusR = bins(bin+1)
        if (radius_crit > radiusL .and. radius_crit < radiusC) then
                    dot_precipitate_density(bin) = dot_precipitate_density(bin) &
                                        + nucleation_rate/(radiusC - radiusL)                        
        endif   
        
            enddo nucleation  
          
end subroutine growth_precipitate



