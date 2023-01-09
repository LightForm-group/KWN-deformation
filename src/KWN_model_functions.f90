module KWN_model_functions

use KWN_precision

contains

function calculate_shear_modulus(T)
    ! Calculate Shear Modulus after McLellan 1987 MPa
    real(pReal), intent(in) :: T !temperature in K
    real(pReal) :: calculate_shear_modulus !shear modulus [Pa]

		calculate_shear_modulus = ( &
		                    27.0 & 
		                    + (21.5 - 27.0) / (650.0 - 273.0) &
		                    * (T - 273.0) &
		                 ) * 1.0e9

end function calculate_shear_modulus

end module KWN_model_functions