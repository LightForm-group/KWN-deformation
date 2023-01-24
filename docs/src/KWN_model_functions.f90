module KWN_model_functions

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawMicrostructure


contains

function calculate_shear_modulus(Temperature)
    implicit none
    ! Calculate Shear Modulus after McLellan 1987 MPa
    real(pReal), intent(in) :: Temperature !temperature in K
    real(pReal) :: calculate_shear_modulus !shear modulus [Pa]

        calculate_shear_modulus = ( 27.0 & 
                                  + (21.5 - 27.0) / (650.0 - 273.0) &
                                  * (Temperature - 273.0) &
                                  ) * 1.0e9

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


function calculate_binary_alloy_critical_radius(Temperature, dst, prm, en)
    implicit none
    ! calculate the critical radius for binary alloys
    real(pReal), intent(in) :: Temperature !temperature in K
    type(tParameters), intent(in) :: prm
    type(tKwnpowerlawMicrostructure), intent(in) :: dst
    integer, intent(in) :: en
    real(pReal) :: deltaGv
    real(pReal) :: calculate_binary_alloy_critical_radius
    
    ! SAM: Added method to calculate the  explicitly
    deltaGv = -R * Temperature &
                 * log( dst%c_matrix(1,en) / prm%ceq_matrix(1) ) &
                 / prm%molar_volume & 
              + prm%misfit_energy

    calculate_binary_alloy_critical_radius = -2.0_pReal * prm%gamma_coherent / deltaGv

end function calculate_binary_alloy_critical_radius


function calculate_beta_star(radius_crit, lattice_param, &
                              diffusion_coefficient, c_matrix, en)
    !! Function for calculating beta_star for binary or ternary mixtures.
    real(pReal), intent(in) :: radius_crit
    real(pReal), intent(in) :: lattice_param
    real(pReal), dimension(:), allocatable, intent(in) ::   &
        diffusion_coefficient  ! diffusion coefficient for Mg and Zn
    real(pReal), dimension(:,:), allocatable, intent(in) :: c_matrix
    integer, intent(in) :: en
    
    real(pReal) :: calculate_beta_star

    !TODO: Eventually the sizes of diffusion coefficient and c_matrix should be adjusted
    !      in the main code, so that there is not empty space in these when working with 
    !      binary alloys. Once that happens we can use the ternary alloy solution for 
    !      all calculations. 

    ! expression of beta star for ternary alloys
    if (c_matrix(2,en) > 0) then
        calculate_beta_star = 4.0_pReal * PI &
                    * radius_crit ** 2.0 / (lattice_param ** 4.0) &
                    * 1 / ( sum( 1 / (diffusion_coefficient(:) * c_matrix(:,en)) ) )
    ! expression of beta star for binary alloys
    else
        calculate_beta_star = 4.0_pReal * PI &
                    * radius_crit ** 2.0 / (lattice_param ** 4.0) &
                    * 1 / ( ( 1 / (diffusion_coefficient(1) * c_matrix(1,en)) ) )
                                        
    endif

end function calculate_beta_star


function calculate_nucleation_rate(nucleation_site_density, zeldovich_factor, beta_star, &
                                   gamma_coherent, radius_crit, Temperature, incubation_time, &
                                   time, en)
    real(pReal), intent(in) :: &
        nucleation_site_density, & !nucleation density [pr/m^3]
        zeldovich_factor, & !Zeldovich factor
        beta_star, & ! in the nucleation rate expression
        gamma_coherent, &       ! coherent precipitate surface energy in J/m^2
        radius_crit, &
        Temperature, & !temperature in K
        incubation_time ! in the nucleation rate expression
	real(pReal),  dimension(:), allocatable, intent(in) :: time ! time [s]
    integer, intent(in) :: en

    real(pReal) :: calculate_nucleation_rate


    calculate_nucleation_rate = nucleation_site_density * zeldovich_factor * beta_star &
                                * exp( &
                                  - 4.0_pReal * PI * gamma_coherent * radius_crit ** 2.0 &
                                     / ( 3.0_pReal * kB * Temperature ) &
                                  - incubation_time / time(en) &
                                  )

end function calculate_nucleation_rate

end module KWN_model_functions




