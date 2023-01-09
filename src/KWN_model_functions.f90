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



end module KWN_model_functions