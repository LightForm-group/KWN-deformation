module KWN_model

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_model_routines, only: interface_composition, growth_precipitate, &
                                  update_precipate_properties, update_diffusion_coefficient, next_time_increment
    use KWN_model_functions, only: calculate_binary_alloy_critical_radius, &
                                   calculate_beta_star, calculate_nucleation_rate, calculate_shear_modulus, &
                                   calculate_yield_stress
    use KWN_io, only: output_results, print_results
    
contains

subroutine run_model(prm, dot, stt, dst, &
                    prm_temp, dot_temp,  stt_temp, dst_temp, &
                    Nmembers, en, &
                    dt &
                    )

    implicit none

    type(tParameters), intent(inout) :: prm, prm_temp
    type(tKwnpowerlawState), intent(inout) ::  dot, dot_temp, stt, stt_temp
    type(tKwnpowerlawMicrostructure), intent(inout) :: dst, dst_temp
    

    integer, intent(in) :: &
        Nmembers, &
        en
   
    real(pReal), intent(inout) :: &
        dt !time step for integration [s]

    integer ::  k

    real(pReal) :: &
        time_record  ! used to record the outputs in files


    INTEGER :: status ! I/O status

    ! time_record is used to record the results in textfiles
    time_record = -dt


    k = 0
    loop_time : do while  (stt%time(en).LE. prm%total_time)
        k=k+1
        
        ! display current state on terminal
        call print_results(prm, stt, dst, en)
        
        ! go to next time increment and update the variables
        ! this routine calculate new state for t = t+dt - main program
        call next_time_increment(prm, dst, dst_temp, dot, dot_temp, stt, stt_temp, dt, en)


        ! Adapt time step - and possibly go back to previous step -  so that the outputs do not vary too much between two time steps
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
            dst = dst_temp
            stt = stt_temp
            prm = prm_temp
            dot = dot_temp

            !because it didn't go well with this time step, decrease the time step by a factor arbitrarily chosen (2)
            dt = 0.5 * dt

        else
        !set the time step (for next step) so that a precipitate cannot grow from more than the space between two adjacent classes

            !not necessary to adapt the time step to the growth rate if there is no precipitate
            if ( dst%total_precipitate_density(en) > 1.0_pReal ) then
                ! condition that guarantees that precipitates move from one class maximum
                dt = min( prm%dt_max, &
                          (prm%bins(1)-prm%bins(0)) / maxval(abs(dst%growth_rate_array)) &
                        )

            else
                dt = min(prm%dt_max, dt*1.2) !increase slightly the time step by an arbitrary factor as long as there are no precipitates
            endif


            ! store the updated values of the outputs in the temporary variables
            dst_temp=dst
            prm_temp=prm
            stt_temp=stt
            dot_temp=dot

          

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