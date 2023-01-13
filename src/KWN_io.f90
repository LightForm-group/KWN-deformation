module KWN_io

    use KWN_precision
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    
contains
    
subroutine read_configuration( &
                            testfolder, &
                            prm, &
                            Temperature, &
                            shape_parameter, &  ! the initial distribution is defined by mean radius, volume fraction and shape parameter of a log normal distribution - see e.g. ref [4]
                            total_time, &  ![s]
                            dt_max, &  ![s]
                            time_record_step, &  ![s]
                            sigma_r, &  ![MPa] - sinepower law for stress
                            A, & ![/s] - sinepower law for stress
                            Q_stress, &  ![J/mol] - activation energy in flow stress law
                            n, &  !exponent in sinepower law for stress
                            incubation,  & !incubation prefactor, either 0 or 1)
                            stoechiometry, &
                            N_elements &
                            )
    

    implicit none

    character*100, intent(out) :: testfolder
    type(tParameters), intent(inout) :: prm
    integer, dimension(:), allocatable, intent(inout) :: stoechiometry !precipitate stoechiometry in the following order : Mg Zn Al
    integer, intent(in) :: N_elements
    real(pReal), intent(out) :: &
        Temperature, &
        shape_parameter, &
        total_time, &
        dt_max, &
        time_record_step, &
        sigma_r, &
        A, &
        Q_stress, &
        n, &
        incubation

    ! local variables
    INTEGER :: status ! I/O status
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
			volume_fraction_initial, &! initial total precipitate distribution
			rho_0, & !initial dislocation density
			rho_s, & !saturation dislocation density
			strain_rate, & ! strain rate in /s
			dislocation_arrangement, & ! constant related to the dislocation density in the vacancy annihilation term, cf [1]
			burgers, & !matrice burgers vector
			jog_formation_energy, & ! formation energy for jogs
			q_dislocation ! activation energy for diffusion at dislocation (pipe diffusion) in J/at - not used yet but to be updated

	! the following variables are allocatable to allow for precipitates with multiple elements (only situations with 2 elements are used here)
	real(pReal), dimension(:), allocatable :: &
			c0_matrix, &            ! initial matrix solute composition in mol fraction : [Mg, Zn]
			ceq_matrix, &           ! equilibrium matrix composition for a flat interface (infinite precipitate radius) in mol fraction : [Mg, Zn]
			diffusion0, &           ! solute diffusivity in m^2/s : [Mg, Zn] - in the present version of the code - the diffusion coefficient is taken as identical for both solute and equal to the slowest diffuser (Mg)
			migration_energy        !  solute migration energy in J/at

    
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
                      incubation

    ! ensure allocatable arrays are allocated to same size as prm arrays
    allocate(migration_energy(N_elements), source=0.0_pReal)
    allocate(diffusion0(N_elements), source=0.0_pReal)
    allocate(c0_matrix(N_elements), source=0.0_pReal)
    allocate(ceq_matrix(N_elements), source=0.0_pReal)




    ! Read the inputs from the input.namelist file
    !!!!!!!!!!!
    OPEN (UNIT=1, FILE='namelist.input', STATUS='OLD', ACTION='READ', IOSTAT=status)
    print*, status
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
    prm%rho_0 = rho_0
    prm%rho_s = rho_s
    prm%strain_rate = strain_rate
    prm%dislocation_arrangement = dislocation_arrangement
    prm%burgers = burgers
    prm%jog_formation_energy = jog_formation_energy
    prm%q_dislocation = q_dislocation
    prm%c0_matrix = c0_matrix
    prm%ceq_matrix = ceq_matrix
    prm%diffusion0 = diffusion0
    prm%migration_energy = migration_energy

    ! Write the namelist to our test folder, for record keeping
    open (unit=1, file=trim(testfolder)//'/namelist.output', status='replace', iostat=status)
    print*, status
    write(1, config)
    close(1)
    
end subroutine read_configuration


subroutine output_results(testfolder, filesuffix, stt, dst, diffusion_coefficient, c_thermal_vacancy, &
                        nucleation_rate, production_rate, annihilation_rate, dislocation_density, &
                        radius_crit, en)

    type(tKwnpowerlawState), intent(in) :: stt
    type(tKwnpowerlawMicrostructure), intent(in) :: dst

    real(pReal), dimension(:), allocatable, intent(in) ::   &
        diffusion_coefficient  ! diffusion coefficient for Mg and Zn

    real(pReal), intent(in) :: &
        c_thermal_vacancy, & ! concentration in thermal vacancies
        production_rate, & ! production rate for excess vacancies
        annihilation_rate, & !annihilation rate for excess vacancies
        nucleation_rate, & ! part/m^3/s
        dislocation_density, & ![/m^2]
        radius_crit !critical radius, [m]
    
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
    results(1,6)=nucleation_rate*1.0e-18
    if (results(1,6)<1.0e-30_pReal) then
        results(1,6)=0.0
    endif
    results(1,5)=radius_crit*1.0e9

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
        write(1, 601) stt%time(en), diffusion_coefficient(1)
    close(1)

    filename = 'results/vacancies_'
    filename = trim(testfolder)//trim(filename)//trim(filesuffix)
    open(1, file = filename,  ACTION="write", position="append")
        write(1, 1001) stt%time(en), stt%c_vacancy(en)/c_thermal_vacancy, &
                       production_rate/c_thermal_vacancy, annihilation_rate/c_thermal_vacancy
    close(1)

    filename = 'results/dislocation_density_'
    filename = trim(testfolder)//trim(filename)//trim(filesuffix)
    open(1, file = filename,  ACTION="write", position="append")
        write(1, 901) stt%time(en), dislocation_density
    close(1)

601 FORMAT(2E40.6)
901 FORMAT(3E40.6)
1001 FORMAT(4E40.6)

end subroutine output_results

end module KWN_io