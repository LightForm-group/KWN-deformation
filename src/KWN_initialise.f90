module KWN_initialise

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_model_routines, only: interface_composition, growth_precipitate , equilibrium_flat_interface
    use KWN_model_functions, only: calculate_binary_alloy_critical_radius
    use KWN_io, only: read_configuration, output_results

contains

subroutine initialise_model_state(prm, dot, stt, dst, &
                                Nmembers, en, &
                                interface_c, &
                                dt, &
                                x_eq_interface &
                                )


    implicit none

    type(tParameters), intent(out) :: prm
    type(tKwnpowerlawState), intent(out) ::  dot, stt
    type(tKwnpowerlawMicrostructure), intent(out) :: dst

    integer, intent(out) :: &
        Nmembers, &
        en


    real(pReal), dimension(:), allocatable, intent(out) :: &
        x_eq_interface  !array with equilibrium concentrations at the interface between matrix and precipitates of each bin

    real(pReal), intent(out) :: &
        interface_c
    

    real(pReal), intent(out) :: &
        dt  !time step for integration [s]



    ! local variables
    integer ::  bin
    real(pReal) :: & 
        radiusL, radiusR, radiusC, & ! used for the calculation of the growth rate in the different bins
        N_0 = 0.0, & !parameter used
        dislocation_density, & ![/m^2]
        production_rate, & ! production rate for excess vacancies
        annihilation_rate, & !annihilation rate for excess vacancies
        integral_dist_function, & ! used to calculate the initial distribution
        mean_particle_strength, & !particle strength for precipitation hardening effect calculation - ref[2]
        nucleation_rate, & ! part/m^3/s
        deltaGv, & ! chemical driving force [J/mol]
        vol_int


  
    Nmembers = 1 ! Nmembers is the number of point in the simulation - here 1 as it's a single point model

    en = 1

 


    ! --------------------------
    ! allocating variables needed for reading configuration file
    allocate(prm%migration_energy(N_elements), source=0.0_pReal)
    allocate(prm%diffusion0(N_elements), source=0.0_pReal)
    allocate(prm%c0_matrix(N_elements), source=0.0_pReal)
    allocate(prm%ceq_matrix(N_elements), source=0.0_pReal)
    allocate(prm%ceq_precipitate(N_elements), source=0.0_pReal)
    allocate(prm%stoechiometry(N_elements+1))






    !!! read the configuration file, using data arrays allocated above
    call read_configuration( & 
                            prm &
                            )

    
    
    !-------- add a backslash to the folder path
    prm%testfolder = trim(prm%testfolder)//'/'


    !----------- allocate variables not needed in configuration file
    ! scalar values
    allocate(dst%total_precipitate_density(Nmembers), source=0.0_pReal)
    allocate(dst%avg_precipitate_radius   (Nmembers), source=0.0_pReal)
    allocate(dst%precipitate_volume_frac  (Nmembers), source=0.0_pReal)
    ! size of the array = number of considered elements
    allocate(dst%c_matrix                 (N_elements,Nmembers), source=0.0_pReal)
    allocate(dst%diffusion_coefficient     (N_elements,Nmembers), source=0.0_pReal)

    ! size of these arrays: number of bins
    allocate(prm%bins(0:prm%kwn_nSteps), source=0.0_pReal)
    allocate(dot%precipitate_density(prm%kwn_nSteps,1), source=0.0_pReal)  ! time derivative of the precipitate density in each bin
    allocate(stt%precipitate_density(prm%kwn_nSteps,1), source=0.0_pReal)  ! precipitate density in each bin
    allocate(stt%normalized_distribution_function(prm%kwn_nSteps,1), source=0.0_pReal) ! distribution function for precipitate density [/m^4]
    allocate(stt%growth_rate_array(prm%kwn_nSteps-1), source=0.0_pReal) ! array containing the growth rate in each bin
    allocate(x_eq_interface(0:prm%kwn_nSteps), source=0.0_pReal) ! equilibrium concentration at the interface taking into account Gibbs Thomson effect (one equilibrium concentration for each bin)
    allocate(stt%time (Nmembers), source=0.0_pReal) ! Time array
    allocate(stt%c_vacancy (Nmembers), source=0.0_pReal) ! Number of excess vacancies
    allocate(stt%yield_stress (Nmembers), source=0.0_pReal) ! Yield stress that will be calculated from solid solution, dislocation and precipitates
    allocate(dot%c_vacancy (Nmembers), source=0.0_pReal) !Time derivative of excess vacancies



    !!! some conversions
    prm%migration_energy = prm%migration_energy / na ! convert form J/mol to J/at
    prm%jog_formation_energy = prm%jog_formation_energy * ev_to_Jat  ! convert from ev to  J/at
    prm%q_dislocation = prm%q_dislocation / na ! convert to J/at
    prm%vacancy_energy = prm%vacancy_energy * ev_to_Jat  ! convert from ev to  J/at
    prm%vacancy_migration_energy = prm%vacancy_migration_energy * ev_to_Jat  ! convert from ev to  J/at


    prm%ceq_precipitate = real(prm%stoechiometry(1:2)) / real(sum(prm%stoechiometry)) ! calculate the concentration of the precipitate from the stoichiometry

      !calculate initial diffusion coefficient
    dst%diffusion_coefficient(:,en) = prm%diffusion0 * exp( -(prm%migration_energy) / prm%Temperature / kb ) ! +2*(dislocation_density)*prm%atomic_volume/prm%burgers&
                            !   *prm%diffusion0*exp(-(prm%q_dislocation )/Temperature/kb)  ! include pipe diffusion


    ! if the enthalpy and entropy are provided, then the equilibrium concentration should be calculated, otherwise take the input value for equilibrium concentration
    if (prm%entropy>0.0_pReal) then
		call 			equilibrium_flat_interface(prm%Temperature,  N_elements,  prm%stoechiometry, &
												   prm%c0_matrix,prm%ceq_matrix, prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
												   dst%diffusion_coefficient, dst%precipitate_volume_frac(en), prm%enthalpy, prm%entropy)
	endif


    ! initial value for the time step
    dt = 0.001
    ! if the initial value for the time step is higher than dt max replace by dt max
    dt = min(prm%dt_max, dt)

    ! define the output file suffix (contains Temperature in °C and strain rate in /s)
    if (prm%Temperature<1273) then
        write(prm%filesuffix, '(I3,"C_strain_rate",ES9.3, ".txt")') int(prm%Temperature)-273, prm%strain_rate
    else
        write(prm%filesuffix, '(I4,"C_strain_rate",ES9.3, ".txt")') int(prm%Temperature)-273, prm%strain_rate
    endif
    !! initialise the bins for the size distribution
    ! SAM: Adjusted binning
    !---------------------------------------------------------------------------------------------------------------------------------
    kwnbins_init: do bin = 0, prm%kwn_nSteps
        if (prm%kwn_step0 < 0.0) then
            ! if initial step is smaller than 0 then make log bins
            prm%bins(bin) = 10**(real(bin,pReal) * prm%kwn_stepsize + prm%kwn_step0)
        else
            ! otherwise linear bins are used
            prm%bins(bin) = real(bin,pReal) * prm%kwn_stepsize + prm%kwn_step0
        endif
    enddo kwnbins_init
    !---------------------------------------------------------------------------------------------------------------------------------

    !initialize some outputs
    stt%growth_rate_array = 0.0_pReal
    stt%precipitate_density = 0.0_pReal
    dst%total_precipitate_density(en) = 0.0_pReal
    dst%avg_precipitate_radius(en) = prm%mean_radius_initial
    dst%precipitate_volume_frac(en) = prm%volume_fraction_initial
    stt%c_vacancy(en) = 0.0_pReal
    dislocation_density = prm%rho_0


    !Calculated the precipitate density from the initial mean radius and volume fraction for an already existing distribution
    if (dst%avg_precipitate_radius(en) > 0) then

        distribution_function : do bin = 1, prm%kwn_nSteps
                                !definition of a log normal distribution
                                radiusL = prm%bins(bin-1)
                                radiusR = prm%bins(bin  )
                                radiusC = (radiusL+radiusR)/2
                                stt%normalized_distribution_function(bin,en) =  1.0_pReal / sqrt(PI * 2.0_pReal) &
                                                                                / prm%shape_parameter / radiusC &
                                                                                * exp( -1.0 / 2.0 * &
                                                                                ( log( radiusC / prm%mean_radius_initial )&
                                                                                + prm%shape_parameter**2 / 2 )**2 / prm%shape_parameter**2 )
                                enddo distribution_function

        !print*, normalized_distribution_function


        ! calculate the integral of the distribution function
        integral_dist_function=0.0_pReal
        do bin = 1,prm%kwn_nSteps
            radiusL = prm%bins(bin-1)
            radiusR = prm%bins(bin  )
            integral_dist_function = integral_dist_function + (radiusR - radiusL) * stt%normalized_distribution_function(bin,en)
        enddo

        !print*, 'integral dist', integral_dist_function


        ! normalized_distribution is not normalised yet
        stt%normalized_distribution_function(:,en) = stt%normalized_distribution_function(:,en) / integral_dist_function
        ! now it is normalised

        !the normalized distribution function gives the shape of the distribution, it needs to be multiplied by the number density N0 such that the initial precipitate fraction is the one given as input
        ! calculate the volume integral of the normalised precipitate distribution
        vol_int = 0.0
        do bin = 1, prm%kwn_nSteps
            radiusL = prm%bins(bin-1)
            radiusR = prm%bins(bin  )
            radiusC = (radiusL + radiusR)/2
            vol_int = vol_int + stt%normalized_distribution_function(bin,en) * radiusC**3 * (radiusR-RadiusL) * 4.0/3.0 * PI
        enddo
        !number density * total volume of precipitates = volume fraction
        N_0 = dst%precipitate_volume_frac(en) / vol_int

        ! now the normalised distribution function is multiplied by the total number density to get the number density per bin size [m^{-4}]
        stt%precipitate_density(:,en) = stt%normalized_distribution_function(:,en) * N_0
        dst%total_precipitate_density(en) = N_0

        stt%precipitate_density(1,en) = 0

        ! to avoid some problems when writing the outputs
        do bin = 1, prm%kwn_nSteps
            if (stt%precipitate_density(bin,en) < 1.0e-50_pReal) then
                stt%precipitate_density(bin,en) = 0.0_pReal
            endif
        enddo


        ! initialise time derivative
        dot%precipitate_density = stt%precipitate_density / dt


        !recompute the initial radius, volume fraction etc to avoid discontinuities between first and next steps


        dst%avg_precipitate_radius(en)  = 0.0_pReal
        dst%total_precipitate_density   = 0.0_pReal
        dst%precipitate_volume_frac(en) = 0.0_pReal



        do bin = 1, prm%kwn_nSteps
            radiusL = prm%bins(bin-1)
            radiusR = prm%bins(bin  )

            !update precipitate density
            dst%total_precipitate_density(en)   =   dst%total_precipitate_density(en) &
                                                + stt%precipitate_density(bin,en) &
                                                * (radiusR - radiusL)
            !update average radius
            dst%avg_precipitate_radius(en)  =   dst%avg_precipitate_radius(en) &
                                                + stt%precipitate_density(bin,en) &
                                                * (radiusR**2.0_pReal - radiusL**2.0_pReal) / 2.0_pReal ! at that stage in m/m^3

            !update volume fraction
             dst%precipitate_volume_frac(en) =  dst%precipitate_volume_frac(en) &
                                                + 1.0_pReal / 6.0_pReal * PI &
                                                * (radiusR+ radiusL)**3.0_pReal &
                                                * (radiusR - radiusL) &
                                                * stt%precipitate_density(bin,en)

        enddo

        !convert from m/m^3 to m
        if (dst%total_precipitate_density(en) > 0.0_pReal) then
            dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                            / dst%total_precipitate_density(en)
        endif
    endif



    !initialize the concentration in the matrix as a function of the volume fraction and composition of precipitate

    dst%c_matrix(:,en) = ( prm%c0_matrix(:) - dst%precipitate_volume_frac(en) * prm%ceq_precipitate(:) ) &
                            / ( 1.0 - dst%precipitate_volume_frac(en) )

  




    !calculate the equilibrium composition at the interface between precipitates and matrix as a function of their size (Gibbs Thomson effect)
    call interface_composition( prm%Temperature,  N_elements, prm%kwn_nSteps, prm%stoechiometry, prm%c0_matrix,prm%ceq_matrix, &
                                prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, prm%bins, prm%gamma_coherent, &
                                R, x_eq_interface, dst%diffusion_coefficient, dst%precipitate_volume_frac(en), prm%misfit_energy)

    !TODO: Have users set N_elements, and test for N_elements==1 here to define a binary alloy
    ! calculate critical radius in the case of a binary alloy
    if (dst%c_matrix(2,en)==0) then
        stt%radius_crit = calculate_binary_alloy_critical_radius(dst, prm, en)
    end if

    !calculate the initial growth rate of precipitates of all sizes
    call growth_precipitate( N_elements, prm%kwn_nSteps, prm%bins, interface_c, x_eq_interface,prm%atomic_volume, &
                            na, prm%molar_volume, prm%ceq_precipitate, stt%precipitate_density, &
                            dot%precipitate_density(:,en), nucleation_rate,  dst%diffusion_coefficient, &
                            dst%c_matrix(:,en), stt%growth_rate_array,stt%radius_crit )


    !the critical radius for dissolution if calculated from the precipitate growth rate array - display it
   print*, ''
   print*, 'Critical radius:', stt%radius_crit*1.0e9, ' nm'


    stt%c_thermal_vacancy = 1.0
    production_rate = 0.0
    annihilation_rate = 0.0

    print*, 'Initial mean radius :', dst%avg_precipitate_radius
    print*, 'Initial precipitate density :', dst%total_precipitate_density(en)
    print*, 'Initial volume fraction :', dst%precipitate_volume_frac(en)
    print*, 'Bulk composition', prm%c0_matrix(:)
    print*, 'Initial matrix composition :', dst%c_matrix(:,en)
    print*, 'Equilibrium composition, precipitate :', prm%ceq_precipitate
    print*, 'Equilibrium composition matrix :', prm%ceq_matrix

    print*, 'Initialising outputs'
    

    call initialise_outputs(prm, stt, dst, nucleation_rate,  &
                                dt, stt%growth_rate_array, &
                               mean_particle_strength, &
                               x_eq_interface, en)
    
    print*, 'Writing outputs'

    call output_results(prm%testfolder, prm%filesuffix, stt, dst,  &
                        nucleation_rate, production_rate, annihilation_rate, dislocation_density, &
                         en)
    print*, 'End initialisation'

end subroutine initialise_model_state


subroutine initialise_outputs( prm, stt, dst, nucleation_rate, &
                               dt, growth_rate_array, &
                               mean_particle_strength, &
                               x_eq_interface, en)
    implicit none

    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawState), intent(in) :: stt
    type(tKwnpowerlawMicrostructure), intent(in) :: dst

    real(pReal), dimension(:), allocatable, intent(in) :: &
        growth_rate_array, &!array that contains the precipitate growth rate of each bin
        x_eq_interface  !array with equilibrium concentrations at the interface between matrix and precipitates of each bin

    real(pReal), intent(in) :: &
        dt, & !time step for integration [s]
        mean_particle_strength, & !particle strength for precipitation hardening effect calculation - ref[2]
        nucleation_rate ! part/m^3/s
    


    integer, intent(in) :: en

    ! local variables
    character*200 :: filename !name of the gile where the outputs will be written
    integer :: bin, i, status

    ! print*, 'file_suffix : ', filesuffix

    !Write the initial precipitate distribution in a textfile
    filename = 'results/initial_precipitation_distribution_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
   !  print*, 'File name:', filename 
    print*, 'Creating initial distribution output file... '
  
    open(unit=1, file = filename,  STATUS='replace', ACTION='WRITE', IOSTAT=status)
 
        print*, 'File opened'
        write(1,*) ' # Bin [m], Precipitate density distribution [/m^4]'

        do bin = 1, prm%kwn_nSteps
            print*, 'bin : ', bin
            if (sum(stt%precipitate_density(:,en)) > 0.0_pReal) then
                write(1, 901) prm%bins(bin), stt%precipitate_density(bin,en) / sum(stt%precipitate_density(:,en))
            else
                write(1, 901) prm%bins(bin), stt%precipitate_density(bin,en)
            endif
        enddo
    close(1)
    print*, 'Creating temperature file... '
    ! record the temperature (for versions where there would be a temperature ramp for example)
    filename = 'results/temperature_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
    open(1,file = filename,  ACTION="write", STATUS="replace" )
        write(1,*) '# Time [s], Temperature [K]'
    close(1)

    ! record the diffusion coefficient
    filename = 'results/diffusion_coefficient_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
    open(1, file = filename,  ACTION="write", STATUS="replace")
        write(1,*) '# Time [s], Diffusion coefficient [m^2/s] '
    close(1)

    ! record the number of excess vacancies
    filename = 'results/vacancies_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
    open(1, file = filename,  ACTION="write", STATUS="replace")
        write(1,*) '# Time [s], c_{ex}/c_{th}, total number of produced vacancies/c_{th}, total number of annihilated vacancies /c_{th}'
    close(1)

    ! record the dislocation density
    filename = 'results/dislocation_density_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
    open(1, file = filename, ACTION="write", STATUS="replace")
        write(1,*) '# Time [s], dislocation density [/m^2]'
    close(1)

    ! record calculated flow stress
    filename = 'results/stress_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
    open(1, file = filename, ACTION="write", STATUS="replace")
        write(1,*) '# Time [s], Stress [Pa]'
    close(1)


    ! this file will be used to store most of the results
    filename = 'results/kinetics_data_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
    open(1, file = filename,  ACTION="write", STATUS="replace")
        write(1,*) '#Time, [s], Average Radius [nm], Total precipitate density [/micron^3], Volume fraction [], Concentration in the matrix [at %]'
    close(1)

    ! Write all the input parameters in a file
    filename = 'results/KWN_parameters_'
    filename = trim(prm%testfolder)//trim(filename)//trim(prm%filesuffix)
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
            WRITE (201,*) ' '
            WRITE (201,*) 'Initial distribution'
            WRITE (201,302) prm%mean_radius_initial
            302 FORMAT ('Initial mean radius: ', E30.6, ' m')
            WRITE (201,303) prm%volume_fraction_initial
            303 FORMAT ('Initial volume fraction: ', F7.3)
            WRITE (201,320) prm%shape_parameter
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
            WRITE (201,314) prm%Temperature-273
            314 FORMAT ('Temperature: ', F7.3, ' [°C]')
        close(201)




601 FORMAT(2E40.6)
901 FORMAT(3E40.6)
1001 FORMAT(4E40.6)



end subroutine initialise_outputs


end module KWN_initialise