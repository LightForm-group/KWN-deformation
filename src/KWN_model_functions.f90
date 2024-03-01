module KWN_model_functions

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawMicrostructure, tKwnpowerlawState


contains

function calculate_shear_modulus(prm)
    implicit none
    type(tParameters), intent(in) :: prm
    real(pReal) :: calculate_shear_modulus !shear modulus [Pa]

    if (prm%shear_modulus > 0.0_pReal) then
        calculate_shear_modulus = prm%shear_modulus
    else
        ! Aluminium temperature dependant shear modulus from McLellan 1987 - https://doi.org/10.1016/0022-3697(87)90147-8
        calculate_shear_modulus = ( 27.0 + (21.5 - 27.0) / (650.0 - 273.0) &
                                  * (prm%Temperature - 273.0) ) * 1.0e9
    end if

end function calculate_shear_modulus


function calculate_dislocation_density(rho_0, rho_s, strain)
    implicit none
    !from Detemple 1995 Physical Review B - Condensed Matter and Materials Physics, 52(1), 125–133.
    real(pReal), intent(in) :: &
                            rho_0, & !initial dislocation density
                            rho_s, & !saturation dislocation density
                            strain   !macroscopic strain
    real(pReal) :: calculate_dislocation_density

    calculate_dislocation_density = rho_s &
                                    * ( 1 & 
                                        - (sqrt(rho_s) - sqrt(rho_0)) &
                                          / sqrt(rho_s) &
                                          * exp( (-1.0 / 2.0) * 86 * strain) &
                                      ) ** 2

end function calculate_dislocation_density


function calculate_binary_alloy_critical_radius(dst, prm, en)
    implicit none
    ! calculate the critical radius for binary alloys
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawMicrostructure), intent(in) :: dst
    integer, intent(in) :: en
    real(pReal) :: deltaGv
    real(pReal) :: calculate_binary_alloy_critical_radius
    
    ! SAM: Added method to calculate the  explicitly
    deltaGv = -R * prm%Temperature &
                 * log( dst%c_matrix(1,en) / prm%ceq_matrix(1) ) &
                 / prm%molar_volume & 
              + prm%misfit_energy

    calculate_binary_alloy_critical_radius = -2.0_pReal * prm%gamma_coherent / deltaGv

end function calculate_binary_alloy_critical_radius


function calculate_beta_star(radius_crit, lattice_param, en, dst)
    !! Function for calculating beta_star for binary or ternary mixtures.
    real(pReal), intent(in) :: radius_crit
    real(pReal), intent(in) :: lattice_param
    type(tKwnpowerlawMicrostructure), intent(in) :: dst
    integer, intent(in) :: en
    
    real(pReal) :: calculate_beta_star

    !TODO: Eventually the sizes of diffusion coefficient and c_matrix should be adjusted
    !      in the main code, so that there is not empty space in these when working with 
    !      binary alloys. Once that happens we can use the ternary alloy solution for 
    !      all calculations. 

    ! expression of beta star for ternary alloys
    if (dst%c_matrix(2,en) > 0) then
        calculate_beta_star = 4.0_pReal * PI &
                    * radius_crit ** 2.0 / (lattice_param ** 4.0) &
                    * 1 / ( sum( 1 / (dst%diffusion_coefficient(:,en) * dst%c_matrix(:,en)) ) )
    ! expression of beta star for binary alloys
    else
        calculate_beta_star = 4.0_pReal * PI &
                    * radius_crit ** 2.0 / (lattice_param ** 4.0) &
                    * 1 / ( ( 1 / (dst%diffusion_coefficient(1,en) * dst%c_matrix(1,en)) ) )
                                        
    endif

end function calculate_beta_star


function calculate_nucleation_rate(prm, stt, dst, &
                                   en)

    type(tKwnpowerlawState), intent(inout) :: stt
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawMicrostructure), intent(inout) :: dst
    integer, intent(in) :: en
    !local variables
    real(pReal) :: calculate_nucleation_rate,&
                   nucleation_site_density, &
                   zeldovich_factor, &
                   beta_star, &
                   incubation_time

    zeldovich_factor = prm%atomic_volume * sqrt(prm%gamma_coherent / ( kB * prm%Temperature ) ) &
                                                / ( 2.0_pReal * PI * dst%radius_crit**2.0 )


    ! expression of beta star for all alloys
    beta_star = calculate_beta_star(dst%radius_crit, prm%lattice_param, &
                                         en, dst)

        !TODO: Have users set N_elements, and test for N_elements==1 here to define a binary alloy
        !TODO: Doug: I think this should be calculated before beta_star in each timestep,
        !            in the setting of initial timestep constants
        !            (it will converge towards the same answer either way, but with slightly
        !             different strain rates early in the simulation)
        ! calculate critical radius in the case of a binary alloy
    if (dst%c_matrix(2,en)==0) then
            dst%radius_crit = calculate_binary_alloy_critical_radius( dst, prm, en)
    end if


    incubation_time = REAL(prm%incubation) * 2.0 &
                        / ( PI * zeldovich_factor**2.0 * beta_star )

 
    nucleation_site_density = sum(dst%c_matrix(:,en)) / prm%atomic_volume

    calculate_nucleation_rate = nucleation_site_density * zeldovich_factor * beta_star &
                                * exp( &
                                  - 4.0_pReal * PI * prm%gamma_coherent * dst%radius_crit ** 2.0 &
                                     / ( 3.0_pReal * kB * prm%Temperature ) &
                                  - incubation_time / stt%time(en) &
                                  )

end function calculate_nucleation_rate

function calculate_yield_stress(dst,prm,stt,en)
    implicit none
    type(tKwnpowerlawState), intent(inout) :: stt
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawMicrostructure), intent(in) :: dst
    integer, intent(in) :: en
    integer :: bin
    real(pReal) :: tau_s,tau_d,tau_p, mu, calculate_yield_stress, &
                    line_tension, obstacle_strength,radiusC, strength_increment,precipitate_density_fraction, bin_width

    !calculate yield stress
    mu = calculate_shear_modulus(prm)

    tau_s=prm%k_s*sum(dst%c_matrix(:,en))**(2.0/3.0)
    print*, 'Solid solution contribution', tau_s*1e-6, 'MPa'
    !print*, 'sum c', sum(dst%c_matrix(:,en))
    !Taylor relation for dislocation contribution
    tau_d=0.3*mu*prm%burgers*sqrt(dst%dislocation_density)
    print*, 'Dislocation contribution', tau_d*1e-6, 'MPa'

    ! Strength model comes from Myhr et al - https://doi.org/10.1016/S1359-6454(00)00301-3 - which is modified from Deschamps and Brechet.

    line_tension = 2*prm%k_p*mu*prm%burgers**2

    obstacle_strength=0.0_pReal
    
    ! Loop through all the bins
    kwnbins: do bin=1,prm%kwn_nSteps
    
        radiusC = prm%bins(bin)
        bin_width = prm%bins(bin) - prm%bins(bin - 1)

        ! Quick if statement prevents division by 0 or nan's from appearing
        if (dst%total_precipitate_density(en) > 0.0_pReal) then
            precipitate_density_fraction = bin_width * stt%precipitate_density(bin,en)/dst%total_precipitate_density(en)
        else
            precipitate_density_fraction = 0.0_pReal
        end if

        if (radiusC < prm%transition_radius) then
            strength_increment = line_tension * (radiusC/prm%transition_radius) * precipitate_density_fraction
        else
            strength_increment = line_tension * precipitate_density_fraction
        end if

        obstacle_strength = obstacle_strength + strength_increment

    enddo kwnbins

    if (dst%avg_precipitate_radius(en) > 0.0_pReal) then    
        tau_p = (obstacle_strength**(3.0_pReal/2.0_pReal))*sqrt(3*dst%precipitate_volume_frac(en)/(2*PI))/(prm%burgers*dst%avg_precipitate_radius(en)*sqrt(line_tension))
    else
        tau_p = 0.0_pReal
    endif

    print*, 'Precipitate contribution', tau_p*1e-6, 'MPa'
    calculate_yield_stress = prm%M*(tau_s + sqrt(tau_p**2+tau_d**2))
    !  print*, 'Yield stress function:', calculate_yield_stress*1e-6, 'MPa'

end function calculate_yield_stress

end module KWN_model_functions




