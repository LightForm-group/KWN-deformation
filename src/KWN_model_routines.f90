module KWN_model_routines

use KWN_precision

contains

subroutine interface_composition(T,  N_elements, N_steps, stoechiometry, &
								c_matrix,ceq_matrix, atomic_volume, na, molar_volume, ceq_precipitate, &
								bins, gamma_coherent, R,  x_eq_interface, diffusion_coefficient, volume_fraction, misfit_energy)
    

	!  find the intersection between stoichiometric line and solubility line for precipitates of different sizes by dichotomy - more information in ref [3] or [6] + [5]
	implicit none
	integer, intent(in) :: N_elements, N_steps
	integer, intent(in), dimension(N_elements+1) :: stoechiometry
	real(pReal), intent(in), dimension(N_elements) :: c_matrix, ceq_precipitate, diffusion_coefficient, ceq_matrix
	real(pReal), intent(in) :: T,  atomic_volume, na, molar_volume, gamma_coherent, R, volume_fraction, misfit_energy
	real(pReal), intent(inout), dimension(0:N_steps) :: x_eq_interface
    real(pReal), intent(in), dimension(0:N_steps) :: bins
	real(pReal) :: xmin, xmax, solubility_product, delta
	integer :: i

	xmin=0.0_pReal
	xmax=1.0_pReal

   interface_equilibrium: 	do i = 0, N_steps

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
												-solubility_product*exp(2.0*molar_volume*gamma_coherent/R/T/bins(i)*real(sum(stoechiometry)) )

										if (delta<0.0_pReal) then
											xmin=x_eq_interface(i)
										else
											xmax=x_eq_interface(i)
										endif

									enddo

								else
									x_eq_interface(i) =ceq_matrix(1)*exp((2.0*molar_volume*gamma_coherent/(R*T*bins(i)*ceq_precipitate(1)))+molar_volume*misfit_energy/(R*T)) ! Gibbs Thomson effect for a precipitate with a single alloying element

								endif

    						enddo    interface_equilibrium

end subroutine interface_composition



subroutine 	growth_precipitate(N_elements, N_steps, bins, interface_c, &
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