module KWN_model

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_model_routines, only: interface_composition, growth_precipitate, &
                                  update_precipate_properties, update_diffusion_coefficient
    use KWN_model_functions, only: calculate_binary_alloy_critical_radius, &
                                   calculate_beta_star, calculate_nucleation_rate, calculate_shear_modulus, &
                                   calculate_yield_stress
    use KWN_io, only: output_results
    
contains

subroutine run_model(prm, dot, stt, dst, &
                    prm_temp,  stt_temp, dst_temp, &
                    Nmembers, en, &
                    dt &
                    )

    implicit none

    type(tParameters), intent(inout) :: prm, prm_temp
    type(tKwnpowerlawState), intent(inout) ::  dot,  stt, stt_temp
    type(tKwnpowerlawMicrostructure), intent(inout) :: dst, dst_temp
    



    integer, intent(in) :: &
        Nmembers, &
        en
   
    real(pReal), intent(inout) :: &
        dt !time step for integration [s]

    !!local variables
    ! the temp variables will be used to store the results of previous time step
    ! ideally, the prm should all stay constant but because the temperature migt change, they can change right now
    ! to do remove temperature from prm and put it somewhere else

    integer ::  bin, k, i

    real(pReal) :: &
        deltaGv, & ! chemical driving force [J/mol]
        radiusL, radiusR, radiusC, & ! used for the calculation of the growth rate in the different bins
        growth_rate, flux, & ! growth rate and flux between different bins for the precipitates
        time_record ! used to record the outputs in files


    real(pReal), dimension(:,:), allocatable :: &
        results !variable to store the results

    ! the 'temp' variables are to store the previous step and adapt the time step at each iteration
    real(pReal), dimension(:), allocatable ::   &
        temp_c_matrix, &
        temp_diffusion_coefficient, &
        temp_x_eq_interface, &
        temp_x_eq_matrix, &
        temp_precipitate_density ,&
        temp_dot_precipitate_density, &
        k1,k2,k3,k4 ! variables used for Runge Kutta integration

    real(pReal) :: &
        temp_total_precipitate_density, &
        temp_avg_precipitate_radius, &
        temp_precipitate_volume_frac, &
        temp_radius_crit, &
        temp_c_vacancy, &
        temp_dislocation_density, &
        dt_temp, &
        Temperature_temp, &
        h !used for Runge Kutta integration

    INTEGER :: status ! I/O status

    ! allocate arrays for Runga Kutta and temporary work
    allocate(k1(prm%kwn_nSteps), source=0.0_pReal) ! Runge Kutta
    allocate(k2(prm%kwn_nSteps), source=0.0_pReal) !
    allocate(k3(prm%kwn_nSteps), source=0.0_pReal)
    allocate(k4(prm%kwn_nSteps), source=0.0_pReal)
    allocate(temp_x_eq_interface(0:prm%kwn_nSteps), source=0.0_pReal)
    allocate(temp_precipitate_density(prm%kwn_nSteps), source=0.0_pReal)
    allocate(temp_dot_precipitate_density(prm%kwn_nSteps), source=0.0_pReal)
    allocate(temp_c_matrix(N_elements), source=0.0_pReal)
    allocate(temp_diffusion_coefficient(N_elements), source=0.0_pReal)
    allocate(temp_x_eq_matrix(N_elements), source=0.0_pReal)
    allocate(results(1,8)) ! the results are stored in this array

    !the following are used to store the outputs from the previous iteration
    temp_diffusion_coefficient(:) = dst%diffusion_coefficient(:,en)
    dst_temp%total_precipitate_density = dst%total_precipitate_density
    temp_avg_precipitate_radius = dst%avg_precipitate_radius(en)
    temp_precipitate_volume_frac = dst%precipitate_volume_frac(en)
    temp_c_matrix(:) = dst%c_matrix(:,en)
    temp_precipitate_density = stt%precipitate_density(:,en)
    temp_dot_precipitate_density = dot%precipitate_density(:,en)
    temp_radius_crit = stt%radius_crit
    temp_x_eq_matrix = prm%ceq_matrix
    temp_x_eq_interface = stt%x_eq_interface
    temp_c_vacancy = stt%c_vacancy(en)
    temp_dislocation_density = prm%rho_0
    Temperature_temp = prm%Temperature
    stt%time(en) = 0.0_pReal
    ! time_record is used to record the results in textfiles
    time_record = -dt

    h = dt



    k = 0
    loop_time : do while  (stt%time(en).LE. prm%total_time)
        k=k+1
        
        ! print*, "dt:", dt
        print*, ' '
        print*, 'Time:', stt%time(en)
        print*, 'Temperature', prm%Temperature
        print*, 'Mean radius : ', dst%avg_precipitate_radius(en)*1e9, 'nm'
        print*, 'Total precipitate density : ' , dst%total_precipitate_density*1e-18 , '/micron^3'
        print*, 'Precipitate volume fraction :',  dst%precipitate_volume_frac(en)
        print*, 'Solute concentration in the matrix' , dst%c_matrix(:,en)
        print*, 'Equilibrium concentration in the matrix' , prm%ceq_matrix(:)
        print*, 'Equilibrium volume fraction', (prm%c0_matrix(1)-prm%ceq_matrix(1))/(prm%ceq_precipitate(1)-prm%ceq_matrix(1))
        print*, 'Nucleation rate :part/micron^3/s ', stt%nucleation_rate*1.0e-18
        print*, 'Critical Radius : ', stt%radius_crit*1e9, 'nm'
        print*, 'Yield stress:', stt%yield_stress*1e-6, 'MPa'
        
        
        
        ! update diffusion coefficient taking into account strain induced vacancies
        call update_diffusion_coefficient(prm, stt, dst, dot, dt, en)                                        
        
        ! calculate nucleation rate
        if (stt%time(en) > 0.0_pReal) then
            stt%nucleation_rate = calculate_nucleation_rate(prm, stt, &
                                                            dst, en)
        else
            stt%nucleation_rate = 0.0_pReal
        endif

        !calculate the precipitate growth in all bins dot%precipitate_density
        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, &
                                    stt%x_eq_interface,prm%atomic_volume,  prm%molar_volume, prm%ceq_precipitate, &
                                    stt%precipitate_density, dot%precipitate_density(:,en), stt%nucleation_rate, &
                                    dst%diffusion_coefficient(:,en), dst%c_matrix(:,en), stt%growth_rate_array, stt%radius_crit )


        ! Runge Kutta 4th order to calculate the derivatives

        ! https://en.wikipedia.org/wiki/Runge–Kutta_methods


        ! Runge Kutta k2 calculation
        ! repeat the calculations above for t = t+dt/2

        h = dt
        k1 = dot%precipitate_density(:,en)

        stt%time(en) = stt%time(en) + h / 2.0
        stt%precipitate_density(:,en) = temp_precipitate_density + h / 2.0 * k1


        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins,&
                                stt%x_eq_interface,prm%atomic_volume, prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), stt%nucleation_rate,&
                                dst%diffusion_coefficient(:,en), dst%c_matrix(:,en), stt%growth_rate_array, stt%radius_crit )

        

        if (stt%time(en) > 0.0_pReal) then
            stt%nucleation_rate = calculate_nucleation_rate(prm, stt, &
                                                            dst, en)
        else
            stt%nucleation_rate = 0.0_pReal
        endif


        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, &
                                stt%x_eq_interface,prm%atomic_volume, prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), stt%nucleation_rate,  dst%diffusion_coefficient(:,en), &
                                dst%c_matrix(:,en), stt%growth_rate_array, stt%radius_crit )


        k2 = dot%precipitate_density(:,en)

        ! Runge Kutta k3 calculation
        stt%precipitate_density(:,en) = temp_precipitate_density + h / 2.0 * k2

        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins,&
                                stt%x_eq_interface,prm%atomic_volume,  prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), stt%nucleation_rate, &
                                dst%diffusion_coefficient, dst%c_matrix(:,en), stt%growth_rate_array, stt%radius_crit )


        k3 = dot%precipitate_density(:,en)

        ! Runge Kutta k4 calculation
        stt%precipitate_density(:,en) = temp_precipitate_density + h * k3
        stt%time(en) = stt%time(en) + h / 2.0


        if (stt%time(en) > 0.0_pReal) then
            stt%nucleation_rate = calculate_nucleation_rate(prm, stt, &
                                                            dst, en)
        else
                stt%nucleation_rate = 0.0_pReal
        endif

        !calculate precipitate growth rate in all bins
        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins,  &
                                stt%x_eq_interface,prm%atomic_volume,  prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), stt%nucleation_rate,  &
                                dst%diffusion_coefficient, dst%c_matrix(:,en), stt%growth_rate_array, stt%radius_crit )



        k4 = dot%precipitate_density(:,en)


        !Runge Kutta, calculate precipitate density in all bins (/m^4)
        stt%precipitate_density(:,en) = temp_precipitate_density + h / 6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)


        ! update precipate (dst) volume frac, density, avg radius, and matrix composition
        call update_precipate_properties(prm, dst, stt, en)


        ! print*, ''

        
            

        ! Adapt time step so that the outputs do not vary too much between two time steps
        !if  either:
        !    - the precipitate distribution in one class/ the vacancy concentration / the concentration in the matrix becomes negative,
        ! or - the growth rate is sufficiently fast for some precipitates to be able to jump from two size classes or more during one time step
        ! then go back to the previous step and decrease the time step


        if  (      (stt%c_vacancy(en) < 0.0_pReal) &
              .OR. (minval(stt%precipitate_density(:,en)) < 0.0_pReal) &
              .OR. (minval(dst%c_matrix(:,en)) < 0.0_pReal) &
              .OR. any(isnan(stt%precipitate_density(:,en))) &
            )  then
        ! go back one step before

            stt%time(en) = stt%time(en) - dt
            dst%total_precipitate_density(en) = temp_total_precipitate_density
            stt%precipitate_density(:,en) = temp_precipitate_density(:)
            dot%precipitate_density(:,en) = temp_dot_precipitate_density(:)
            dst%avg_precipitate_radius(en) = temp_avg_precipitate_radius
            dst%precipitate_volume_frac(en) = temp_precipitate_volume_frac
            dst%c_matrix(:,en) = temp_c_matrix
            !stt%c_vacancy(en) = temp_c_vacancy
            !dislocation_density = temp_dislocation_density
            prm%Temperature = Temperature_temp
            !radius_crit = temp_radius_crit
            prm%ceq_matrix = temp_x_eq_matrix
            stt%x_eq_interface = temp_x_eq_interface
            !diffusion_coefficient(1) = temp_diffusion_coefficient

            !decrease time step by a factor arbitrarily chose (2)
            dt = 0.5 * dt


        else
        !set the time step so that a precipitate cannot grow from more than the space between two adjacent classes

            !not necessary to adapt the time step to the growth rate if there is no precipitate
            if ( dst%total_precipitate_density(en) > 1.0_pReal ) then

                dt = min( prm%dt_max, &
                          (prm%bins(1)-prm%bins(0)) / maxval(abs(stt%growth_rate_array)) &
                        )
                ! print*,'dt growth rate', (prm%bins(1) - prm%bins(0)) / maxval(abs(growth_rate_array))

            else
                dt = min(prm%dt_max, dt*1.2) !increase slightly the time step by an arbitrary factor as long as there are no precipitates
            endif


            ! store the new values of the outputs in the temporary variables
            temp_dot_precipitate_density = dot%precipitate_density(:,en)
            dst_temp%total_precipitate_density = dst%total_precipitate_density
            temp_precipitate_density = stt%precipitate_density(:,en)
            temp_avg_precipitate_radius = dst%avg_precipitate_radius(en)
            temp_precipitate_volume_frac = dst%precipitate_volume_frac(en)
            temp_c_matrix = dst%c_matrix(:,en)
            !temp_radius_crit = radius_crit
            temp_x_eq_matrix = prm%ceq_matrix
            temp_x_eq_interface = stt%x_eq_interface
            !temp_c_vacancy = stt%c_vacancy(en)
            !temp_dislocation_density = dislocation_density
            Temperature_temp = prm%Temperature
            !temp_diffusion_coefficient = diffusion_coefficient(1)
            !stt%yield_stress=calculate_yield_stress(calculate_shear_modulus(Temperature),dislocation_density,dst,prm,en)
          

            if (time_record < stt%time(en)) then !record the outputs every 'time_record' seconds


                call output_results(prm%testfolder, prm%filesuffix, stt, dst, &
                                     en)


                ! next time for which the outputs should be written
                        if (prm%time_record_step > 0) then
                            !Save data linearly
                            time_record = time_record + prm%time_record_step

                        else
                            !save data logarithimically
                            time_record = time_record + 10**(INT(LOG10(stt%time(en)))-1)

                  endif

            endif
        endif
        
    
    end do loop_time


	

end subroutine run_model


end module KWN_model