module KWN_model_routines

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_model_functions, only : calculate_shear_modulus, calculate_dislocation_density


contains


subroutine set_initial_timestep_constants(prm, stt, dot, Temperature, sigma_r, A, Q_stress, n, dt, en, &
                                          diffusion_coefficient, c_thermal_vacancy, dislocation_density, &
                                          production_rate, annihilation_rate)

    implicit none
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawState), intent(inout) :: dot, stt
    real(pReal), intent(in) :: &
        Temperature, & !temperature in K
        sigma_r, & ! constant in the sinepowerlaw for flow stress [MPa]
        A, &  ! constant in the sinepowerlaw for flow stress  [/s]
        Q_stress, &  ! activation energy in the sinepowerlaw for flow stress [J/mol]
        n, & ! stress exponent in the sinepower law for flow stress
        dt !time step for integration [s]
    integer, intent(in) :: en

    real(pReal), dimension(:), allocatable, intent(out) ::   &
        diffusion_coefficient  ! diffusion coefficient for Mg and Zn

    real(pReal), intent(out) :: &
        c_thermal_vacancy, & ! concentration in thermal vacancies
        production_rate, & ! production rate for excess vacancies
        annihilation_rate, & !annihilation rate for excess vacancies
        dislocation_density ![/m^2]

    real(pReal) :: &
        mu, & !shear modulus [Pa]
        flow_stress, & ! flow stress in the material [Pa]
        c_j, & ! jog concentration - ref [1]
        strain !macroscopic strain
    

	diffusion_coefficient = prm%diffusion0 * exp( -(prm%migration_energy) / (Temperature * kb) )
	mu = calculate_shear_modulus(Temperature)


	! if there is deformation, calculate the vacancy related parameters


	flow_stress = sigma_r * asinh(((prm%strain_rate / (A)) * exp(Q_stress / ( 8.314 * Temperature) )) ** (1/n))
	c_thermal_vacancy = 23.0 * exp(-prm%vacancy_energy / (kB * Temperature) )
	c_j = exp(-prm%jog_formation_energy / (kB * Temperature) )
	strain = prm%strain_rate * stt%time(en)
	dislocation_density = calculate_dislocation_density(prm%rho_0, prm%rho_s, strain)

	! calculate production and annihilation rate of excess vacancies as described in ref [1] and [3]
	production_rate =   prm%vacancy_generation * flow_stress * prm%atomic_volume * prm%strain_rate &
							 / prm%vacancy_energy &
						+ 0.5 * c_j * prm%atomic_volume * prm%strain_rate / (4.0 * prm%burgers**3)

	annihilation_rate = prm%vacancy_diffusion0 * exp( -prm%vacancy_migration_energy / (kB * Temperature) ) &
						* ( dislocation_density / prm%dislocation_arrangement**2 &
							+ 1.0 / prm%vacancy_sink_spacing**2 ) &
						* stt%c_vacancy(en)
	
	! variation in vacancy concentration
	dot%c_vacancy(en) = production_rate - annihilation_rate

	! total number of vacancies
	stt%c_vacancy(en) = stt%c_vacancy(en) + dot%c_vacancy(en) * dt



	!update the diffusion coefficient as a function of the vacancy concentration
	! the first term adds the contribution of excess vacancies,the second adds the contribution of dislocation pipe diffusion
	diffusion_coefficient = prm%diffusion0 * exp( -prm%migration_energy / (Temperature * kb) ) &
							* (1.0 + stt%c_vacancy(en) / c_thermal_vacancy )! &
						!   +2*(dislocation_density)*prm%atomic_volume/prm%burgers&
						!   *prm%diffusion0*exp(-(prm%q_dislocation )/Temperature/kb)

end subroutine set_initial_timestep_constants


subroutine update_precipate_properties(prm, dst, stt, en)

    implicit none
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawState), intent(in) ::  stt
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
                                            delta = x_eq(1)**stoechiometry(1)*&
											((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))*diffusion_coefficient(1)/diffusion_coefficient(2)*&
											(x_eq(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2) &
											-solubility_product
                                    ! uncomment these lines if the precipitate is mainly made of principal elements (e.g. Ni3Al in Ni alloys)
                                    ! ToDO : change that as the solubility product function might not be adapted for all stroichiometries (not necessarily a monotonic function of c_eq)
                                    		!delta = x_eq(1)**stoechiometry(1)*&
											!((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))*diffusion_coefficient(1)/diffusion_coefficient(2)*&
											!(x_eq(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2)*(1-x_eq(1)-((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))*diffusion_coefficient(1)/diffusion_coefficient(2)*&
											!(x_eq(1)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction)))**stoechiometry(3)&
											!-solubility_product

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
                                                ((c_matrix(2)+real(stoechiometry(2))/real(stoechiometry(1))*diffusion_coefficient(1)/diffusion_coefficient(2)*&
                                                (x_eq_interface(i)*(1-volume_fraction)-c_matrix(1)))/(1-volume_fraction))**stoechiometry(2)&
                                                -solubility_product*exp(2.0*molar_volume*gamma_coherent/R/Temperature/bins(i)*real(sum(stoechiometry)) )

                                        if (delta<0.0_pReal) then
                                            xmin=x_eq_interface(i)
                                        else
                                            xmax=x_eq_interface(i)
                                        endif

                                    enddo

                                else
                                    x_eq_interface(i) =ceq_matrix(1)*exp((2.0*molar_volume*gamma_coherent/(R*Temperature*bins(i)*ceq_precipitate(1)))+molar_volume*misfit_energy/(R*Temperature)) ! Gibbs Thomson effect for a precipitate with a single alloying element

                                endif

                            enddo    interface_equilibrium

end subroutine interface_composition



subroutine  growth_precipitate(N_elements, N_steps, bins, interface_c, &
            x_eq_interface,atomic_volume, na, molar_volume, ceq_precipitate, precipitate_density, &
            dot_precipitate_density, nucleation_rate, diffusion_coefficient, c_matrix, growth_rate_array , radius_crit)

    ! calculate precipitate growth rate for all class sizes
    implicit none


    integer, intent(in) :: N_Steps, N_elements
    real(pReal), intent(in), dimension(0:N_steps) :: bins, x_eq_interface
    real(pReal), intent(in), dimension(N_steps) :: precipitate_density
    real(pReal), intent(in), dimension(N_elements) :: ceq_precipitate, diffusion_coefficient, c_matrix
    real(pReal), intent(in) :: atomic_volume, na, molar_volume,  nucleation_rate
    real(pReal), intent(inout), dimension(N_steps) :: dot_precipitate_density
    real(pReal), intent(inout), dimension(0:N_steps) :: growth_rate_array
    real(pReal), intent(inout)::  radius_crit
    real(pReal) :: radiusC, radiusL, radiusR, interface_c, growth_rate, flux
    integer :: bin



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

end subroutine growth_precipitate


end module KWN_model_routines