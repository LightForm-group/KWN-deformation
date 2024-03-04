module KWN_model_routines

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_model_functions, only : calculate_shear_modulus, calculate_dislocation_density, calculate_yield_stress, calculate_nucleation_rate,&
                                    calculate_temperature, calculate_misfit_energy


contains


subroutine update_diffusion_coefficient(prm, stt, dst, dot, dt, en)

    implicit none
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawMicrostructure), intent(inout) :: dst
    type(tKwnpowerlawState), intent(inout) :: dot, stt
    real(pReal), intent(in) :: &
        dt !time step for integration [s]
    integer, intent(in) :: en



    real(pReal) :: &
        mu, & !shear modulus [Pa]
        flow_stress, & ! flow stress in the material [Pa]
        c_j, & ! jog concentration - ref [1]
        strain, & !macroscopic strain
        production_rate, & ! production rate for excess vacancies
        annihilation_rate!annihilation rate for excess vacancies
    

	dst%diffusion_coefficient(:,en) = prm%diffusion0 * exp( -(prm%migration_energy) / (prm%Temperature * kb) )
	mu = calculate_shear_modulus(prm)



    dst%dislocation_density = calculate_dislocation_density(prm%rho_0, prm%rho_s, strain)
    	! if there is deformation, calculate the vacancy related parameters

    ! two situations: if the user defines parameters for precipitation hardening and solid solution hardening, use them
    ! otherwise; use the asinh function for the flow stress 
    if(prm%sigma_r>0.0_pReal) then
        dst%yield_stress = prm%sigma_r * asinh(((prm%strain_rate / (prm%A)) * exp(prm%Q_stress / ( 8.314 * prm%Temperature) )) ** (1/prm%n))    
    else    
        dst%yield_stress=calculate_yield_stress(dst,prm,stt, en)
    endif
    dst%c_thermal_vacancy = 23.0 * exp(-prm%vacancy_energy / (kB * prm%Temperature) ) !TODO change this 23 to 2.3
	c_j = exp(-prm%jog_formation_energy / (kB * prm%Temperature) )
	strain = prm%strain_rate * stt%time(en)
	

	! calculate production and annihilation rate of excess vacancies as described in ref [1] and [3]
	production_rate =   prm%vacancy_generation * dst%yield_stress(en) * prm%atomic_volume * prm%strain_rate &
							 / prm%vacancy_energy &
						+ 0.5 * c_j * prm%atomic_volume * prm%strain_rate / (4.0 * prm%burgers**3)

	annihilation_rate = prm%vacancy_diffusion0 * exp( -prm%vacancy_migration_energy / (kB * prm%Temperature) ) &
						* ( dst%dislocation_density / prm%dislocation_arrangement**2 &
							+ 1.0 / prm%vacancy_sink_spacing**2 ) &
						* stt%c_vacancy(en)
	
	! variation in vacancy concentration
	dot%c_vacancy(en) = production_rate - annihilation_rate

	! total number of vacancies
	stt%c_vacancy(en) = stt%c_vacancy(en) + dot%c_vacancy(en) * dt



	!update the diffusion coefficient as a function of the vacancy concentration
	! the first term adds the contribution of excess vacancies,the second adds the contribution of dislocation pipe diffusion
	dst%diffusion_coefficient(:,en) = prm%diffusion0 * exp( -prm%migration_energy / (prm%Temperature * kb) ) &
							* (1.0 + stt%c_vacancy(en) / dst%c_thermal_vacancy )! &
						!   +2*(dislocation_density)*prm%atomic_volume/prm%burgers&
						!   *prm%diffusion0*exp(-(prm%q_dislocation )/Temperature/kb)

end subroutine update_diffusion_coefficient


subroutine update_precipate_properties(prm, dst, stt, en)

    implicit none
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawState), intent(inout) ::  stt
    type(tKwnpowerlawMicrostructure), intent(inout) :: dst
    integer, intent(in) :: en

    integer :: bin
    real(pReal) :: radiusL, radiusR

	!stt%precipitate density contains all information to calculate mean radius, volume fraction and avg radius - calculate them now
	dst%precipitate_volume_frac(en) = 0.0_pReal
	dst%total_precipitate_density(en) = 0.0_pReal
	dst%avg_precipitate_radius(en) = 0.0_pReal


	! update radius, total precipitate density and volume fraction

	kwnbins:    do bin=1,prm%kwn_nSteps
					radiusL = prm%bins(bin-1)
					radiusR = prm%bins(bin  )

					!update precipitate density
					dst%total_precipitate_density = dst%total_precipitate_density &
													+ stt%precipitate_density(bin,en) &
													* (radiusR - radiusL)
					!update average radius
					dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
													+ stt%precipitate_density(bin,en) &
													* (radiusR**2.0_pReal - radiusL**2.0_pReal) / 2.0_pReal ! at that stage in m/m^3
					!update volume fraction

					dst%precipitate_volume_frac(en) = dst%precipitate_volume_frac(en) &
													+ 1.0_pReal / 6.0_pReal * PI &
													* (radiusR + radiusL)**3.0_pReal &
													* (radiusR - radiusL) &
													* stt%precipitate_density(bin,en)


				enddo kwnbins
	! mean radius from m/m^3 to m
	if (dst%total_precipitate_density(en) > 0.0_pReal) then
				dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
												/ dst%total_precipitate_density(en)
				
	endif


	!update matrix composition

	 dst%c_matrix(:,en) = (prm%c0_matrix(:) - dst%precipitate_volume_frac(en) * prm%ceq_precipitate(:)) &
							/ (1 - dst%precipitate_volume_frac(en))


end subroutine update_precipate_properties

subroutine equilibrium_flat_interface(T,  N_elements, stoechiometry, &
									 c_matrix,x_eq, atomic_volume, na, molar_volume, ceq_precipitate, &
									 diffusion_coefficient, volume_fraction, enthalpy, entropy)

	!TODO change this to take prm, stt etc as inputs rather than the variables separately
    !  find the intersection between stoichiometric line and solubility line for precipitates of different sizes by dichotomy - more information in ref [3] or [6] + [5]
	implicit none
	integer, parameter :: pReal = selected_real_kind(25)
	integer, intent(in) :: N_elements
	integer, intent(in), dimension(N_elements+1) :: stoechiometry
	real(pReal), intent(in), dimension(N_elements) :: c_matrix, ceq_precipitate, diffusion_coefficient
	real(pReal), intent(in) :: T,  atomic_volume, na, molar_volume,  enthalpy, entropy, volume_fraction
	real(pReal), intent(inout), dimension(N_elements+1) :: x_eq
	real(pReal) :: xmin, xmax, solubility_product, delta
	integer :: i

	xmin=0.0_pReal
	xmax=1.0_pReal

  


   								solubility_product=exp(entropy/8.314-enthalpy/8.314/T)


								xmin=0.0! 
								xmax=atomic_volume*na/molar_volume*ceq_precipitate(1) ! the equilibrium concentration at the interface cannot be higher than the concentration in the precipitate

								do while (abs(xmin-xmax)/xmax >0.0001.AND.xmax>1.0e-8)


									x_eq(1)=(xmin+xmax)/2.0

									!  find the intersection between stoichiometric line and solubility line by dichotomy
									! delta=0 ==> intersection between stoichiometry and solubility line - please see e.g. Nicolas et al Acta Mat 2003 or ref 3 for more explanations 

                                    ! "sufficient" approximation if the precipitate is mainly made of alloying elements (the main element is ignored in the solubility product)
                                          !  delta = x_eq(1)**stoechiometry(1)*&
										!	((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))*diffusion_coefficient(1)/diffusion_coefficient(2)*&
										!	(x_eq(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2) &
										!	-solubility_product
                                    ! uncomment these lines if the precipitate is mainly made of principal elements (e.g. Ni3Al in Ni alloys)
                                    ! ToDO : change that as the solubility product function might not be adapted for all stroichiometries (not necessarily a monotonic function of c_eq)
                                    		delta = x_eq(1)**stoechiometry(1)*&
											((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))&
                                            *diffusion_coefficient(1)/diffusion_coefficient(2)*&
											(x_eq(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2)*&
                                            (1-x_eq(1)-((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))&
                                            *diffusion_coefficient(1)/diffusion_coefficient(2)*&
											(x_eq(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction)))**stoechiometry(3)&
											-solubility_product

									if (delta<0.0_pReal) then
										xmin=x_eq(1)
									else
										xmax=x_eq(1)
									endif

								enddo

								! stoichiometry line to find x_eq
								x_eq(2)=c_matrix(2)-real(stoechiometry(2))/real(stoechiometry(1))*(c_matrix(1)-x_eq(1))
								



return
end subroutine

subroutine interface_composition(Temperature,  N_elements, N_steps, stoechiometry, &
                                c_matrix,ceq_matrix, atomic_volume, na, molar_volume, ceq_precipitate, &
                                bins, gamma_coherent, R,  x_eq_interface, diffusion_coefficient, volume_fraction, misfit_energy)

    !  find the intersection between stoichiometric line and solubility line for precipitates of different sizes by dichotomy - more information in ref [3] or [6] + [5]
    implicit none
    integer, intent(in) :: N_elements, N_steps
    integer, intent(in), dimension(N_elements+1) :: stoechiometry
    real(pReal), intent(in), dimension(N_elements) :: c_matrix, ceq_precipitate, diffusion_coefficient, ceq_matrix
    real(pReal), intent(in) :: Temperature,  atomic_volume, na, molar_volume, gamma_coherent, R, volume_fraction, misfit_energy
    real(pReal), intent(inout), dimension(0:N_steps) :: x_eq_interface
    real(pReal), intent(in), dimension(0:N_steps) :: bins
    ! local variables
    real(pReal) :: xmin, xmax, solubility_product, delta
    integer :: i


    xmin=0.0_pReal
    xmax=1.0_pReal

   interface_equilibrium:   do i = 0, N_steps

                                ! the solubility product is only necessary in a ternary alloy as the interface energy has a simple expression for a binary alloy
                                if (stoechiometry(2)>0) then
                                !if (1==1) then

                                    solubility_product=ceq_matrix(1)**stoechiometry(1)*ceq_matrix(2)**stoechiometry(2)


                                    xmin=ceq_matrix(1)! the equilibrium concentration at the interface for a precipitate of size r cannot be lower than that of a flat interface
                                    xmax=atomic_volume*na/molar_volume*ceq_precipitate(1) ! the equilibrium concentration at the interface cannot be higher than the concentration in the precipitate

                                    do while (abs(xmin-xmax)/xmax >0.0001.AND.xmax>1.0e-8)


                                        x_eq_interface(i)=(xmin+xmax)/2.0

                                        !  find the intersection between stoichiometric line and solubility line by dichotomy
                                        ! delta=0 ==> intersection between stoichiometry and solubility line
                                        delta = x_eq_interface(i)**stoechiometry(1)*&
                                                ((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))&
                                                *diffusion_coefficient(1)/diffusion_coefficient(2)*&
                                                (x_eq_interface(i)*(1-volume_fraction)-c_matrix(1)))/&
                                                (1-volume_fraction))**stoechiometry(2)&
                                                -solubility_product*exp(2.0*molar_volume*gamma_coherent/R/Temperature/bins(i)&
                                                *real(sum(stoechiometry)) )

                                        if (delta<0.0_pReal) then
                                            xmin=x_eq_interface(i)
                                        else
                                            xmax=x_eq_interface(i)
                                        endif

                                    enddo

                                else
                                    x_eq_interface(i) =ceq_matrix(1)&
                                    *exp((2.0*molar_volume*gamma_coherent/(R*Temperature*bins(i)*ceq_precipitate(1)))&
                                    +molar_volume*misfit_energy/(R*Temperature)) ! Gibbs Thomson effect for a precipitate with a single alloying element

                                endif

                            enddo    interface_equilibrium

end subroutine interface_composition


!TODO replace this function to take as entry only dot, stt, etc. 
subroutine  growth_precipitate(N_elements, N_steps, bins,  &
            x_eq_interface,atomic_volume, molar_volume, ceq_precipitate, precipitate_density, &
            dot_precipitate_density, nucleation_rate, diffusion_coefficient, c_matrix, growth_rate_array , radius_crit)

    ! calculate precipitate growth rate for all class sizes
    implicit none

  

    integer, intent(in) :: N_Steps, N_elements
    real(pReal), intent(in), dimension(0:N_steps) :: bins
    real(pReal), intent(inout), dimension(N_steps) :: precipitate_density
    real(pReal), intent(in), dimension(N_elements) :: ceq_precipitate, diffusion_coefficient, c_matrix
    real(pReal), intent(in) :: atomic_volume,  molar_volume,  nucleation_rate
    real(pReal), intent(inout), dimension(N_steps) :: dot_precipitate_density
    real(pReal), intent(inout), dimension(0:N_steps) :: growth_rate_array, x_eq_interface
    real(pReal), intent(inout)::  radius_crit
    
    !local variables
    real(pReal) :: radiusC, radiusL, radiusR, interface_c, growth_rate, flux
    integer :: bin

  
    dot_precipitate_density=0.0_pReal

   ! the growth rate is stored to change the time step in the main program
    growth_rate_array = diffusion_coefficient(1)/bins&
    * (c_matrix(1)    - x_eq_interface) &
    / (atomic_volume*na/molar_volume*ceq_precipitate(1) - x_eq_interface)



    kwnbins_growth: do bin = 1, N_Steps-1

                        ! consider two classes and their interface (radius_c)
                        radiusC = bins(bin  ) ! center
                        radiusL = bins(bin-1) ! left
                        radiusR = bins(bin+1) ! right

                        ! concentration at the interface between matrix and precipitate of the considered bin
                        interface_c =x_eq_interface(bin)
                        ! classical growth rate equation
                        growth_rate = growth_rate_array(bin)
                        ! if the growth rate is positive, precipitates grow (smaller class -> bigger class) and the flux is positive
                        if (growth_rate > 0.0_pReal) then
                            flux = precipitate_density(bin)*growth_rate

                        ! if the growth rate is positive, precipitates shrink (bigger class -> smaller class) and the flux is negative
                        else
                            flux = precipitate_density(bin+1)*growth_rate

                        endif

                        dot_precipitate_density(bin  ) = dot_precipitate_density(bin ) - flux/(radiusC - radiusL)
                        dot_precipitate_density(bin+1) = dot_precipitate_density(bin+1) + flux/(radiusR - radiusC)

                    ! populate the class of the critical radius with nucleating particle
                        ! in binary alloys, the critical radius can be explicitely calculated and it's made in the main program
                        ! for ternary alloys, the critical radius is calculated as the bin at which the growth rate is zero
                        if (c_matrix(2)>0) then
                        
                            if (growth_rate_array(bin-1)<0 .and. growth_rate_array(bin+1)>0) then
                                radius_crit=radiusC
                                    dot_precipitate_density(bin+1) = dot_precipitate_density(bin+1) &
                                                            + nucleation_rate/(radiusR - radiusC)
                            endif
                        else
                            if (radiusL<=radius_crit.and.radiusC>radius_crit) then
                                dot_precipitate_density(bin+1) = dot_precipitate_density(bin+1) &
                                                            + nucleation_rate/(radiusR - radiusC)
                            endif

                        endif
                    
                    
                    enddo kwnbins_growth
    ! empty the first bin to avoid precipitate accumulation
    dot_precipitate_density(1) = 0.0_pReal
    precipitate_density(1) = 0.0_pReal      
      


end subroutine growth_precipitate


subroutine next_time_increment(prm, dst, dst_temp, dot, dot_temp, stt, stt_temp, dt, en)
    !this is the routine in the main loop, calculate all variables for t=t+dt
    implicit none

    type(tParameters), intent(inout) :: prm
    type(tKwnpowerlawMicrostructure), intent(inout) :: dst, dst_temp
    type(tKwnpowerlawState), intent(inout) :: dot, dot_temp, stt, stt_temp
    real(pReal), intent(in) :: &
        dt !time step for integration [s]
    integer, intent(in) :: en
    !local variables 
    real(pReal), dimension(:), allocatable ::   &
        k1,k2,k3,k4 ! variables used for Runge Kutta integration

    real(pReal) :: &
        h !used for Runge Kutta integration

    INTEGER :: status ! I/O status

    ! allocate arrays for Runge Kutta 
    allocate(k1(prm%kwn_nSteps), source=0.0_pReal) ! Runge Kutta
    allocate(k2(prm%kwn_nSteps), source=0.0_pReal) !
    allocate(k3(prm%kwn_nSteps), source=0.0_pReal)
    allocate(k4(prm%kwn_nSteps), source=0.0_pReal)

        ! update the temperature if considering cyclic heating.
        if (prm%heating_freq > 0.0_pReal) then
            prm%Temperature = calculate_temperature(stt,prm,en)

            ! update the equilibrium Gibbs-Thomson effect if the temperature changes
            call interface_composition( prm%Temperature,  N_elements, prm%kwn_nSteps, prm%stoechiometry, prm%c0_matrix,prm%ceq_matrix, &
            prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, prm%bins, prm%gamma_coherent, &
            R, dst%x_eq_interface, dst%diffusion_coefficient, dst%precipitate_volume_frac(en), prm%misfit_energy)

            ! update the equilibrium concentration if entropy is provided, otherwise provide a warning.
            if (prm%entropy>0.0_pReal) then
                call 			equilibrium_flat_interface(prm%Temperature,  N_elements,  prm%stoechiometry, &
                                                           prm%c0_matrix,prm%ceq_matrix, prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
                                                           dst%diffusion_coefficient, dst%precipitate_volume_frac(en), prm%enthalpy, prm%entropy)
            else
                print*,'WARNING: Entropy not provided for cyclic heating.'
            endif

        endif

        ! SAM: Hard coded misfit energy parameter.
        prm%misfit_energy = calculate_misfit_energy(prm)

        print*,'Misfit Energy',prm%misfit_energy

        ! update diffusion coefficient taking into account strain induced vacancies
        call update_diffusion_coefficient(prm, stt, dst, dot, dt, en)     
        
        ! calculate nucleation rate
        if (stt%time(en) > 0.0_pReal) then
            dst%nucleation_rate = calculate_nucleation_rate(prm, stt, &
                                                            dst, en)
        else
            dst%nucleation_rate = 0.0_pReal
        endif

        !calculate the precipitate growth in all bins dot%precipitate_density
        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, &
                                    dst%x_eq_interface,prm%atomic_volume,  prm%molar_volume, prm%ceq_precipitate, &
                                    stt%precipitate_density, dot%precipitate_density(:,en), dst%nucleation_rate, &
                                    dst%diffusion_coefficient(:,en), dst%c_matrix(:,en), dst%growth_rate_array, dst%radius_crit )


        ! Runge Kutta 4th order to calculate the derivatives

        ! https://en.wikipedia.org/wiki/Runge–Kutta_methods


        ! Runge Kutta k2 calculation
        ! repeat the calculations above for t = t+dt/2

        h = dt
        k1 = dot%precipitate_density(:,en)

        stt%time(en) = stt%time(en) + h / 2.0
        stt%precipitate_density(:,en) =  stt_temp%precipitate_density(:,en) + h / 2.0 * k1


        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins,&
                                dst%x_eq_interface,prm%atomic_volume, prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), dst%nucleation_rate,&
                                dst%diffusion_coefficient(:,en), dst%c_matrix(:,en), dst%growth_rate_array, dst%radius_crit )

        

        if (stt%time(en) > 0.0_pReal) then
            dst%nucleation_rate = calculate_nucleation_rate(prm, stt, &
                                                            dst, en)
        else
            dst%nucleation_rate = 0.0_pReal
        endif


        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, &
                                dst%x_eq_interface,prm%atomic_volume, prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), dst%nucleation_rate,  dst%diffusion_coefficient(:,en), &
                                dst%c_matrix(:,en), dst%growth_rate_array, dst%radius_crit )


        k2 = dot%precipitate_density(:,en)

        ! Runge Kutta k3 calculation
        stt%precipitate_density(:,en) = stt_temp%precipitate_density(:,en) + h / 2.0 * k2

        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins,&
                                dst%x_eq_interface,prm%atomic_volume,  prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), dst%nucleation_rate, &
                                dst%diffusion_coefficient, dst%c_matrix(:,en), dst%growth_rate_array, dst%radius_crit )


        k3 = dot%precipitate_density(:,en)

        ! Runge Kutta k4 calculation
        stt%precipitate_density(:,en) = stt_temp%precipitate_density(:,en) + h * k3
        stt%time(en) = stt%time(en) + h / 2.0


        if (stt%time(en) > 0.0_pReal) then
            dst%nucleation_rate = calculate_nucleation_rate(prm, stt, &
                                                            dst, en)
        else
                dst%nucleation_rate = 0.0_pReal
        endif

        !calculate precipitate growth rate in all bins
        call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins,  &
                                dst%x_eq_interface,prm%atomic_volume,  prm%molar_volume, prm%ceq_precipitate, &
                                stt%precipitate_density, dot%precipitate_density(:,en), dst%nucleation_rate,  &
                                dst%diffusion_coefficient, dst%c_matrix(:,en), dst%growth_rate_array, dst%radius_crit )



        k4 = dot%precipitate_density(:,en)


        !Runge Kutta, calculate precipitate density in all bins (/m^4)
        stt%precipitate_density(:,en) = stt_temp%precipitate_density(:,en) + h / 6.0 * (k1 + 2.0*k2 + 2.0*k3 + k4)


        ! update precipate (dst) volume frac, density, avg radius, and matrix composition
        call update_precipate_properties(prm, dst, stt, en)




end subroutine next_time_increment



end module KWN_model_routines