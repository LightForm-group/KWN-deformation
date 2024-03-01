module KWN_io

    use KWN_precision
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_parameters

contains
    
subroutine read_configuration( &
                            prm)
    

    implicit none

    type(tParameters), intent(inout) :: prm
        

    ! local variables
    INTEGER :: status ! I/O status
	integer :: &
			kwn_nSteps, &           ! discretization in r-space
            incubation              ! for cases where nucleation is considered, set to 1 to consider incubation time
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
			volume_fraction_initial, &! initial total precipitate distribution
            shape_parameter, &
			rho_0, & !initial dislocation density
			rho_s, & !saturation dislocation density
			strain_rate, & ! strain rate in /s
            Temperature, &
            total_time, & !total heat treatment time for the simulation
			dislocation_arrangement, & ! constant related to the dislocation density in the vacancy annihilation term, cf [1]
			burgers, & !matrice burgers vector
			jog_formation_energy, & ! formation energy for jogs
			q_dislocation, & ! activation energy for diffusion at dislocation (pipe diffusion) in J/at - not used yet but to be updated
            enthalpy, & ! enthalpy of precipitation
            entropy, &  ! entropy of precipitation
            shear_modulus, & ! matrix shear modulus
            ! option 1 for flow stress calculation
            sigma_r, &  ![MPa] - sinepower law for stress
            A, & ![/s] - sinepower law for stress
            Q_stress, &  ![J/mol] - activation energy in flow stress law
            n, &  !exponent in sinepower law for stress
            ! option 2 for flow stress calculation
            k_s, & !constant parameter in regard to solute strength
            k_p, & !constant parameter in regard to precipitate strength
            transition_radius, & ! Transition radius between bypassing and shearing
            M, & ! Taylor Factor
            dt_max, &
            time_record_step


	! the following variables are allocatable to allow for precipitates with multiple elements (only situations with 2 elements are used here)
	real(pReal), dimension(:), allocatable :: &
			c0_matrix, &            ! initial matrix solute composition in mol fraction : [Mg, Zn]
			ceq_matrix, &           ! equilibrium matrix composition for a flat interface (infinite precipitate radius) in mol fraction : [Mg, Zn]
			diffusion0, &           ! solute diffusivity in m^2/s : [Mg, Zn] - in the present version of the code - the diffusion coefficient is taken as identical for both solute and equal to the slowest diffuser (Mg)
			migration_energy        !  solute migration energy in J/at

    integer, dimension(:), allocatable :: stoechiometry !precipitate stoechiometry in the following order : Mg Zn Al

    character*100 :: testfolder !folder where the input file is


    ! define namelist for reading model configuration
    namelist /config/ kwn_nSteps, kwn_stepsize, kwn_step0, lattice_param, atomic_volume, &
                      molar_volume, misfit_energy, gamma_coherent, vacancy_generation, &
                      vacancy_sink_spacing, vacancy_energy, vacancy_migration_energy, &
                      vacancy_diffusion0, mean_radius_initial, volume_fraction_initial, &
                      rho_0, rho_s, strain_rate, dislocation_arrangement, burgers, &
                      jog_formation_energy, q_dislocation, c0_matrix, ceq_matrix, &
                      diffusion0, migration_energy, & 
                      testfolder, Temperature, stoechiometry, shape_parameter, &
                      total_time, dt_max, time_record_step, sigma_r, A, Q_stress, n, &
                      incubation, enthalpy, entropy, k_s, shear_modulus, & !constant parameter in regard to solute strength
                	  k_p, & !constant parameter in regard to precipitate strength
               		  transition_radius, & ! Transition radius between bypassing and shearing
               		  M ! Taylor Factor for yield stress calculation




    ! set default values for parameters in case the user does not define them
    ! set to 1 to consider incubation time
    incubation=0
    enthalpy=0.0_pReal ! if enthalpy and entropy are not given, they are set to 0 and ignored
    entropy=0.0_pReal ! if enthaly and entropy are given, equilibrium concentration is calculated from the solubility product

    ! the following have been set for aluminium and allow to calculate flow stress 
    sigma_r = 0.0 ! 100MPa for aluminium, by default no external stress so we set this to zero
    A = 5.2140e-06 
    Q_stress = 6.0000e+04 
    n = 6.6831e+00 
	!ToDO - if these parameters are given, they allow to calculate the flow stress as a function of solid solution hardening, precipitation and dislocation
    ! if used they should be calibrated for each temperature - possibly each strain rate... 
	k_p=0.035_pReal ! 
	k_s=683.0e+06_pReal
	M=2.0_pReal
	transition_radius=3.3e-9_pReal
    shear_modulus=0.0_pReal

    ! no deformation parameters given -> no deformation 
    strain_rate=0.0_pReal
! ! if nothing about vacancies is specified, ignore them

    vacancy_generation=0.0_pReal
    ! random values for this to allow user not to specify them in case no deformation is needed
    vacancy_energy = 1.0_pReal ! arbitrarily high 
    vacancy_migration_energy=1.0_pReal ! arbitrarily high
    vacancy_diffusion0 = 1.000e-05 
    jog_formation_energy = 3.000e-01 
    vacancy_sink_spacing = 5.000e-05 
    dislocation_arrangement = 1.000e+01 
    rho_0 = 1.000e+14 
    rho_s = 1.000e+14 
    ! activation energy for pipe diffusion - not considered so far but might be useful in the future
    q_dislocation = 1.083e+05 

    ! if nothing is said about the initial precipitation state, consider there are no precipitates
    volume_fraction_initial=0.0_pReal
    mean_radius_initial=0.0_pReal
    shape_parameter = 2.0000e-01 
    ! no elastic strain energy specified ==> considered as negligible
    misfit_energy=0.0_pReal
    ! default value for max time step for integration
    dt_max=0.5
    !default value for period to store the outputs
    time_record_step=1.0_pReal

    ! ensure allocatable arrays are allocated to same size as prm arrays
    allocate(migration_energy(N_elements), source=0.0_pReal)
    allocate(diffusion0(N_elements), source=0.0_pReal)
    allocate(c0_matrix(N_elements), source=0.0_pReal)
    allocate(ceq_matrix(N_elements), source=0.0_pReal)
    allocate(stoechiometry(N_elements+1))

    print*, 'Reading input file...'

    ! Read the inputs from the input.namelist file
    !!!!!!!!!!!
    OPEN (UNIT=1, FILE='namelist.input', STATUS='OLD', ACTION='READ', IOSTAT=status)
    print*,'' 
    read(1, config )
    CLOSE(1)

   

    prm%kwn_nSteps = kwn_nSteps
    prm%kwn_stepsize = kwn_stepsize
    prm%kwn_step0 = kwn_step0
    prm%lattice_param = lattice_param
    prm%atomic_volume = atomic_volume
    prm%molar_volume = molar_volume
    prm%misfit_energy = misfit_energy
    prm%gamma_coherent = gamma_coherent
    prm%vacancy_generation = vacancy_generation
    prm%vacancy_sink_spacing = vacancy_sink_spacing
    prm%vacancy_energy = vacancy_energy
    prm%vacancy_migration_energy = vacancy_migration_energy
    prm%vacancy_diffusion0 = vacancy_diffusion0
    prm%mean_radius_initial = mean_radius_initial
    prm%volume_fraction_initial = volume_fraction_initial
    prm%shape_parameter = shape_parameter ! mean_radius*shape_parameter= standard devitation
    prm%rho_0 = rho_0
    prm%rho_s = rho_s
    prm%strain_rate = strain_rate
    prm%Temperature=Temperature
    prm%shear_modulus=shear_modulus
    prm%dislocation_arrangement = dislocation_arrangement
    prm%burgers = burgers
    prm%jog_formation_energy = jog_formation_energy
    prm%q_dislocation = q_dislocation
    prm%c0_matrix = c0_matrix
    prm%ceq_matrix = ceq_matrix
    prm%diffusion0 = diffusion0
    prm%migration_energy = migration_energy
    prm%enthalpy = enthalpy 
    prm%entropy = entropy
    prm%sigma_r = sigma_r 
    prm%A = A
    prm%Q_stress = Q_stress
    prm%n = n
    prm%k_p=k_p
    prm%k_s=k_s
    prm%M=M
    prm%transition_radius=transition_radius
    prm%stoechiometry=stoechiometry
    prm%total_time=total_time ! heat treatment time for the simulation 
    prm%dt_max=dt_max
    prm%time_record_step=time_record_step
    prm%testfolder=testfolder
    prm%incubation=incubation
    !print*, 'Writing output parameter file...'
    ! Write the namelist to our test folder, for record keeping
     !open (unit=2, file=trim(testfolder)//'/namelist.output', status='replace', iostat=status)
     !print*, ''
     !write(2, config)
     !close(2)
    !print*, 'Output file written'
    !write (*, config)
    
end subroutine read_configuration


subroutine output_results(testfolder, filesuffix, stt, dst, &
                         en)

    type(tKwnpowerlawState), intent(in) :: stt
    type(tKwnpowerlawMicrostructure), intent(in) :: dst



    
    integer, intent(in) :: en
    
    character*100, intent(in) :: filesuffix !the file suffix contains the temperature and strain rate used for the simulation
    character*100, intent(in) :: testfolder !folder where the input file is

    ! local variables
    real(pReal), dimension(:,:), allocatable :: &
        results !variable to store the results
    character*100 :: filename !name of the gile where the outputs will be written
    integer :: bin, i

    
    allocate(results(1,8)) ! the results are stored in this array

    ! write outputs in textfiles
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
    results(1,6)=dst%nucleation_rate*1.0e-18
    if (results(1,6)<1.0e-30_pReal) then
        results(1,6)=0.0
    endif
    results(1,5)=dst%radius_crit*1.0e9

    filename='results/kinetics_data_'
    filename=trim(testfolder)//trim(filename)//trim(filesuffix)

    open(1, file = filename,  ACTION="write", position="append")
        WRITE(1,13) (results(1,i), i=1,8)
        13 FORMAT(F40.6,F40.6,E40.6,E40.6,E40.6, E40.6, 2E40.6 )
    close(1)

    ! writes the current distribution
    filename='results/precipitation_distribution_'
    filename=trim(testfolder)//trim(filename)//trim(filesuffix)

    open(2, file = filename,  ACTION="write", STATUS="replace")
        WRITE(2,'(E40.15)') stt%time(en), stt%precipitate_density(:,en)
    close(2)

    filename = 'results/diffusion_coefficient_'
    filename = trim(testfolder)//trim(filename)//trim(filesuffix)

    open(1, file = filename,  ACTION="write", position="append")
        write(1, 601) stt%time(en), dst%diffusion_coefficient(1,en)
    close(1)

    filename = 'results/stress_'
    filename = trim(testfolder)//trim(filename)//trim(filesuffix)

    open(1, file = filename,  ACTION="write", position="append")
        write(1, 601) stt%time(en), dst%yield_stress(en)
    close(1)


    filename = 'results/vacancies_'
    filename = trim(testfolder)//trim(filename)//trim(filesuffix)
    open(1, file = filename,  ACTION="write", position="append")
        write(1, 901) stt%time(en), stt%c_vacancy(en),  dst%c_thermal_vacancy 
    close(1)

    filename = 'results/dislocation_density_'
    filename = trim(testfolder)//trim(filename)//trim(filesuffix)
    open(1, file = filename,  ACTION="write", position="append")
        write(1, 601) stt%time(en), dst%dislocation_density
    close(1)

601 FORMAT(2E40.6)
901 FORMAT(3E40.6)
1001 FORMAT(4E40.6)

end subroutine output_results


subroutine print_results(prm, stt, dst, en)
!to display results on the commandline
        
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawMicrostructure), intent(in) :: dst
    type(tKwnpowerlawState), intent(in) ::  stt
    integer, intent(in) :: en
    ! print the current system state
    print*, ' '
    print*, 'Time:', stt%time(en)
    print*, 'Temperature', prm%Temperature
    print*, 'Mean radius : ', dst%avg_precipitate_radius(en)*1e9, 'nm'
    print*, 'Total precipitate density : ' , dst%total_precipitate_density*1e-18 , '/micron^3'
    print*, 'Precipitate volume fraction :',  dst%precipitate_volume_frac(en)
    print*, 'Solute concentration in the matrix' , dst%c_matrix(:,en)
    print*, 'Equilibrium concentration in the matrix' , prm%ceq_matrix(:)
    print*, 'Equilibrium volume fraction', (prm%c0_matrix(1)-prm%ceq_matrix(1))/(prm%ceq_precipitate(1)-prm%ceq_matrix(1))
    print*, 'Nucleation rate :part/micron^3/s ', dst%nucleation_rate*1.0e-18
    print*, 'Critical Radius : ', dst%radius_crit*1e9, 'nm'
    print*, 'Yield stress:', dst%yield_stress*1e-6, 'MPa'
        


end subroutine print_results

end module KWN_io