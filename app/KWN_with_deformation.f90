!--------------------------------------------------------------------------------------------------
!> @author Madeleine Bignon, University of Manchester
!> @author Pratheek Shanthraj, University of Manchester
! KWN precipitation model including the effect of deformation via excess vacancy concentration
!---------------------------------------------------------------------------------------------------
!References:
! [1] Robson, J. D. (2020). Metallurgical and Materials Transactions A: Physical Metallurgy and Materials Science, 51(10), 5401–5413. https://doi.org/10.1007/s11661-020-05960-5
! [2] Deschamps, A., Livet, F., & Bréchet, Y. (1998). . Acta Materialia, 47(1), 281–292. https://doi.org/10.1016/S1359-6454(98)00293-6
! [3] Bignon, M., Shanthraj, P., & Robson, J. D. (2022). Acta Materialia, 234, 118036. https://doi.org/10.1016/J.ACTAMAT.2022.118036
! [4] Deschamps, A., & De Geuser, F. (2011). Journal of Applied Crystallography, 44(2), 343–352. https://doi.org/10.1107/S0021889811003049
! [5] Perez, M. (2005). Scripta Materialia, 52(8), 709–712. https://doi.org/10.1016/j.scriptamat.2004.12.026
! [6] Nicolas, M., & Deschamps, A. (2003).  Acta Materialia, 51(20), 6077–6094. https://doi.org/10.1016/S1359-6454(03)00429-4
program KWN

    use KWN_parameters
    use KWN_data_types, only: tParameters, tKwnpowerlawState, tKwnpowerlawMicrostructure
    use KWN_model_routines, only: interface_composition, growth_precipitate
    use KWN_initialise, only: initialise_model_state
    
	implicit none


	!--------------------------------------------------------------------------------------------------
	! containers for parameters and state
	type(tParameters) :: prm
	type(tKwnpowerlawState) ::  dot, &
								stt
	type(tKwnpowerlawMicrostructure) :: dst

	integer :: &
		Nmembers, &
		N_elements, & ! number of different elements in the precipitate
		en

	integer :: 	bin, k, i

	integer, dimension(:), allocatable :: stoechiometry !precipitate stoechiometry in the following order : Mg Zn Al

	real(pReal), dimension(:,:), allocatable :: &
		results, & !variable to store the results
		normalized_distribution_function !normalised distribution for the precipitate size
	real(pReal), dimension(:), allocatable :: &
		growth_rate_array, &!array that contains the precipitate growth rate of each bin
		x_eq_interface  !array with equilibrium concentrations at the interface between matrix and precipitates of each bin


	real(pReal) :: &
		T, & !temperature in K
		deltaGv, & ! chemical driving force [J/mol]
		interface_energy, & ![J/m^2]
		radius_crit, & !critical radius, [m]
		nucleation_site_density, & !nucleation density [pr/m^3]
		zeldovich_factor, & !Zeldovich factor
		beta_star, & ! in the nucleation rate expression
		incubation_time, & ! in the nucleation rate expression
		nucleation_rate,& ! part/m^3/s
		radiusL, radiusR, radiusC, & ! used for the calculation of the growth rate in the different bins
		growth_rate, flux, & ! growth rate and flux between different bins for the precipitates
		interface_c, & !interface composition between matrix and a precipitate
		time_record, & ! used to record the outputs in files
		time_record_step, & ! time step for the output [s]
		flow_stress, & ! flow stress in the material [Pa]
		c_thermal_vacancy, & ! concentration in thermal vacancies
		dislocation_density, & ![/m^2]
		production_rate, & ! production rate for excess vacancies
		annihilation_rate, & !annihilation rate for excess vacancies
		mu, & !shear modulus [Pa]
		c_j, & ! jog concentration - ref [1]
		strain, & !macroscopic strain
		shape_parameter, & !shape parameter in the log normal distribution of the precipitates - ref [4]
		sigma_r, & ! constant in the sinepowerlaw for flow stress [MPa]
		A, &  ! constant in the sinepowerlaw for flow stress  [/s]
		incubation, & ! incubation prefactor either 0 or 1
		Q_stress, &  ! activation energy in the sinepowerlaw for flow stress [J/mol]
		n ! stress exponent in the sinepower law for flow stress

	! the 'temp' variables are to store the previous step and adapt the time step at each iteration
	real(pReal), dimension(:), allocatable ::   &
		diffusion_coefficient, &  ! diffusion coefficient for Mg and Zn
		temp_c_matrix, &
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
		temp_diffusion_coefficient, &
		temp_dislocation_density, &
		dt, & !time step for integration [s]
		dt_temp, &
		dt_max, & ! max time step for integration [s]
		total_time, & ![s]
		T_temp, &
		h !used for Runge Kutta integration


	character*100 :: filename !name of the gile where the outputs will be written
	character*100 :: filesuffix !the file suffix contains the temperature and strain rate used for the simulation
	character*100 :: testfolder !folder where the input file is

	INTEGER :: status ! I/O status


    call initialise_model_state(prm, dot, stt, dst, &
                                Nmembers, N_elements, en, &
                                stoechiometry, normalized_distribution_function, &
                                T, radius_crit, interface_c, time_record_step, &
                                c_thermal_vacancy, shape_parameter, sigma_r, A, &
                                incubation, Q_stress, n, diffusion_coefficient, &
                                dt, dt_max, total_time, growth_rate_array, &
                                x_eq_interface, &
                                filesuffix, testfolder &
                                )

    ! allocate arrays for Runga Kutta and temporary work
	allocate(k1(prm%kwn_nSteps), source=0.0_pReal) ! Runge Kutta
	allocate(k2(prm%kwn_nSteps), source=0.0_pReal) !
	allocate(k3(prm%kwn_nSteps), source=0.0_pReal)
	allocate(k4(prm%kwn_nSteps), source=0.0_pReal)
	allocate(temp_x_eq_interface(0:prm%kwn_nSteps), source=0.0_pReal)
	allocate(temp_precipitate_density(prm%kwn_nSteps), source=0.0_pReal)
	allocate(temp_dot_precipitate_density(prm%kwn_nSteps), source=0.0_pReal)
	allocate(temp_c_matrix(N_elements), source=0.0_pReal)
	allocate(temp_x_eq_matrix(N_elements), source=0.0_pReal)
	allocate(results(1,8)) ! the results are stored in this array

	!the following are used to store the outputs from the previous iteration
	temp_diffusion_coefficient = diffusion_coefficient(en)
	temp_total_precipitate_density = dst%total_precipitate_density(en)
	temp_avg_precipitate_radius = dst%avg_precipitate_radius(en)
	temp_precipitate_volume_frac = dst%precipitate_volume_frac(en)
	temp_c_matrix(:) = dst%c_matrix(:,en)
	temp_precipitate_density = stt%precipitate_density(:,en)
	temp_dot_precipitate_density = dot%precipitate_density(:,en)
	temp_radius_crit = radius_crit
	temp_x_eq_matrix = prm%ceq_matrix
	temp_x_eq_interface = x_eq_interface
	temp_c_vacancy = stt%c_vacancy(en)
	temp_dislocation_density = prm%rho_0
	T_temp = T
	stt%time(en) = 0.0_pReal
	! time_record is used to record the results in textfiles
	time_record = -dt
	nucleation_rate = 0.0_pReal
	h = dt


    k = 0
	loop_time : do while  (stt%time(en).LE. total_time)


  					k=k+1
     		 		print*, "dt:", dt
    		 		print*, 'Time:', stt%time(en)
    		 		print*, 'Temperature', T
    	     		print*, 'Mean radius : ', dst%avg_precipitate_radius(en)*1e9, 'nm'
					diffusion_coefficient = prm%diffusion0*exp(-(prm%migration_energy )/T/kb)
					mu=(27.0+(21.5-27.0)/(650.0-273.0)*(T-273.0))*1.0e9; !McLellan 1987 MPa


					! if there is deformation, calculate the vacancy related parameters


					flow_stress=sigma_r*asinh(((prm%strain_rate/(A))*exp(Q_stress/8.314/T))**(1/n))
					c_thermal_vacancy=23.0*exp(-prm%vacancy_energy/kB/T)
					c_j=exp(-prm%jog_formation_energy/kB/T)
					strain=prm%strain_rate*stt%time(en)
					!from Detemple 1995 Physical Review B - Condensed Matter and Materials Physics, 52(1), 125–133.
					dislocation_density=prm%rho_s*(1- (sqrt(prm%rho_s)-sqrt(prm%rho_0))/sqrt(prm%rho_s)*exp(-1.0/2.0*86*(strain))  )**2


    				! calculate production and annihilation rate of excess vacancies as described in ref [1] and [3]
					production_rate = 	prm%vacancy_generation*flow_stress*prm%atomic_volume/prm%vacancy_energy*prm%strain_rate &
	                				  	+ 0.5*c_j*prm%atomic_volume/4.0/prm%burgers**3*prm%strain_rate

					annihilation_rate =	prm%vacancy_diffusion0*exp(-prm%vacancy_migration_energy/kB/T) &
										*(dislocation_density/prm%dislocation_arrangement**2+1.0/prm%vacancy_sink_spacing**2)*stt%c_vacancy(en)
					


  					! variation in vacancy concentration
  					dot%c_vacancy(en) = production_rate-annihilation_rate

  					! total number of vacancies
  					stt%c_vacancy(en) = stt%c_vacancy(en)+dot%c_vacancy(en)*dt

					!update the diffusion coefficient as a function of the vacancy concentration
					! the first term adds the contribution of excess vacancies,the second adds the contribution of dislocation pipe diffusion
  					diffusion_coefficient = prm%diffusion0*exp(-(prm%migration_energy )/T/kb)&
  	 										*(1.0+stt%c_vacancy(en)/c_thermal_vacancy  )! &
  	 									!	+2*(dislocation_density)*prm%atomic_volume/prm%burgers&
  	 						 			!	*prm%diffusion0*exp(-(prm%q_dislocation )/T/kb)


				

    				! calculate nucleation rate
    				nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume

    				zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kB/T) &
                     / 2.0_pReal/PI/radius_crit/radius_crit


					! expression of beta star for ternary alloys
					if (dst%c_matrix(2,en)>0) then
   						beta_star = 4.0_pReal*PI&
              						* radius_crit*radius_crit/(prm%lattice_param**4.0) &
              						*1/(sum(1/diffusion_coefficient(:)*1/dst%c_matrix(:,en) ))
              		! expression of beta star for binary alloys
              		else
              			beta_star = 4.0_pReal*PI&
              						* radius_crit*radius_crit/(prm%lattice_param**4.0) &
              						*1/((1/diffusion_coefficient(1)*1/dst%c_matrix(1,en) ))
!-----------------------------------------------------------------------------------------------------------------------------------
										! SAM: Added method to calculate the  explicitly
										deltaGv = -R*T/prm%molar_volume*log(dst%c_matrix(1,en)/prm%ceq_matrix(1)) + prm%misfit_energy

										radius_crit = -2.0_pReal*prm%gamma_coherent / (deltaGv)
									
										
!-----------------------------------------------------------------------------------------------------------------------------------
              		endif


    				incubation_time =  incubation*2.0/PI/zeldovich_factor/zeldovich_factor/beta_star
    				print*, 'Incubation time', incubation_time

    				if (stt%time(en) > 0.0_pReal) then
      					nucleation_rate =   nucleation_site_density*zeldovich_factor*beta_star &
                      						* exp( &
                             				- 4.0_pReal*PI*prm%gamma_coherent*radius_crit*radius_crit/3.0_pReal/kB/T &
                             				- incubation_time/stt%time(en) )
						print*, 'nucleation rate', nucleation_rate*1e-6, '/cm^3'
						

    				else
      					nucleation_rate = 0.0_pReal
    				endif

    				dot%precipitate_density=0.0*dot%precipitate_density

    				!calculate the precipitate growth in all bins dot%precipitate_density
    				call 	growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c,&
    		 									x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    		 									stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate, &
    		  									diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit )


   					! empty the first bin to avoid precipitate accumulation
    				!dot%precipitate_density(0,en)=0.0_pReal
    				dot%precipitate_density(1,en)=0.0_pReal
    				!stt%precipitate_density(0,en)=0.0_pReal
    				stt%precipitate_density(1,en)=0.0_pReal


  					! Runge Kutta 4th order to calculate the derivatives

  					! https://en.wikipedia.org/wiki/Runge–Kutta_methods


  					! Runge Kutta k2 calculation
  					! repeat the calculations above for t = t+dt/2

  					h=dt
  					k1=dot%precipitate_density(:,en)

  					stt%time(en)=stt%time(en)+h/2.0
  					stt%precipitate_density(:,en)=temp_precipitate_density+h/2.0*k1


   					call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c,&
    		 								x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    		 								stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,&
    		   								diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit )

    				nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume
    				zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kB/T) &
                     					/ 2.0_pReal/PI/radius_crit/radius_crit

					! expression of beta star for ternary alloys
					if (dst%c_matrix(2,en)>0) then
   						beta_star = 4.0_pReal*PI&
              						* radius_crit*radius_crit/(prm%lattice_param**4.0) &
              						*1/(sum(1/diffusion_coefficient(:)*1/dst%c_matrix(:,en) ))
              		! expression of beta star for binary alloys
              		else
              			beta_star = 4.0_pReal*PI&
              						* radius_crit*radius_crit/(prm%lattice_param**4.0) &
              						*1/((1/diffusion_coefficient(1)*1/dst%c_matrix(1,en) ))
              		endif

    				incubation_time = incubation*2.0/PI/zeldovich_factor/zeldovich_factor/beta_star

    				if (stt%time(en) > 0.0_pReal) then
      					nucleation_rate = nucleation_site_density*zeldovich_factor*beta_star &
                      	* exp( &
                        - 4.0_pReal*PI*prm%gamma_coherent*radius_crit*radius_crit/3.0_pReal/kB/T &
                        - incubation_time/stt%time(en) )


   					else
      					nucleation_rate = 0.0_pReal
    				endif


    				dot%precipitate_density=0.0*dot%precipitate_density

    				call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c,&
    		 								x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    		 								stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,  diffusion_coefficient, &
    		 	 							dst%c_matrix(:,en), growth_rate_array, radius_crit )



    				!dot%precipitate_density(0,en)=0.0_pReal
    				dot%precipitate_density(1,en)=0.0_pReal
    				!stt%precipitate_density(0,en)=0.0_pReal
    				stt%precipitate_density(1,en)=0.0_pReal
 					k2=dot%precipitate_density(:,en)

 					! Runge Kutta k3 calculation
   					stt%precipitate_density(:,en)=temp_precipitate_density+h/2.0*k2
    				dot%precipitate_density=0.0*dot%precipitate_density

    				call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c, &
    	 									x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
    	 									stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate, &
    	  									diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit )

    				! empty the first bin to avoid precipitate accumulation
    				!dot%precipitate_density(0,en)=0.0_pReal
    				dot%precipitate_density(1,en)=0.0_pReal
    				!stt%precipitate_density(0,en)=0.0_pReal
    				stt%precipitate_density(1,en)=0.0_pReal

					k3=dot%precipitate_density(:,en)

  					! Runge Kutta k4 calculation
  					stt%precipitate_density(:,en)=temp_precipitate_density+h*k3
					stt%time(en)= stt%time(en) +h/2.0

    				nucleation_site_density = sum(dst%c_matrix(:,en))/prm%atomic_volume
    				zeldovich_factor = prm%atomic_volume*sqrt(prm%gamma_coherent/kB/T) &
                     				/ 2.0_pReal/PI/radius_crit/radius_crit
					! expression of beta star for ternary alloys
					if (dst%c_matrix(2,en)>0) then
   						beta_star = 4.0_pReal*PI&
              						* radius_crit*radius_crit/(prm%lattice_param**4.0) &
              						*1/(sum(1/diffusion_coefficient(:)*1/dst%c_matrix(:,en) ))
              		! expression of beta star for binary alloys
              		else
              			beta_star = 4.0_pReal*PI&
              						* radius_crit*radius_crit/(prm%lattice_param**4.0) &
              						*1/((1/diffusion_coefficient(1)*1/dst%c_matrix(1,en) ))
              		endif
   					incubation_time =  incubation*2.0/PI/zeldovich_factor/zeldovich_factor/beta_star





    				if (stt%time(en) > 0.0_pReal) then
      						nucleation_rate = nucleation_site_density*zeldovich_factor*beta_star &
                      		* exp( &
                             - 4.0_pReal*PI*prm%gamma_coherent*radius_crit*radius_crit/3.0_pReal/kB/T &
                             - incubation_time/stt%time(en) &
                            )
    				else
      						nucleation_rate = 0.0_pReal
    				endif

    				dot%precipitate_density=0.0*dot%precipitate_density
    				!calculate precipitate growth rate in all bins
     				call growth_precipitate(N_elements, prm%kwn_nSteps, prm%bins, interface_c, &
     		 								x_eq_interface,prm%atomic_volume, na, prm%molar_volume, prm%ceq_precipitate, &
     		  								stt%precipitate_density, dot%precipitate_density(:,en), nucleation_rate,  &
     		  								diffusion_coefficient, dst%c_matrix(:,en), growth_rate_array, radius_crit )



    				!dot%precipitate_density(0,en)=0.0_pReal
    				dot%precipitate_density(1,en)=0.0_pReal
    				!stt%precipitate_density(0,en)=0.0_pReal
    				stt%precipitate_density(1,en)=0.0_pReal

    				k4=dot%precipitate_density(:,en)


    				!Runge Kutta, calculate precipitate density in all bins (/m^4)
    				stt%precipitate_density(:,en)=temp_precipitate_density + h/6.0*(k1+2.0*k2+2.0*k3+k4)

    				!stt%precipitate density contains all information to calculate mean radius, volume fraction and avg radius - calculate them now
    				dst%precipitate_volume_frac(en) = 0.0_pReal
    				dst%total_precipitate_density(en) = 0.0_pReal
    				dst%avg_precipitate_radius(en) = 0.0_pReal

					! update radius, total precipitate density and volume fraction

					kwnbins:	do bin=1,prm%kwn_nSteps
   									radiusL = prm%bins(bin-1)
      								radiusR = prm%bins(bin  )

    								!update precipitate density
    								dst%total_precipitate_density = dst%total_precipitate_density &
   								   									+ stt%precipitate_density(bin,en) &
   								   									*(radiusR - radiusL)
   									!update average radius
    								dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     								+ stt%precipitate_density(bin,en) &
                                     								* (radiusR**2.0_pReal - radiusL**2.0_pReal)/2.0_pReal ! at that stage in m/m^3
   									!update volume fraction

    								dst%precipitate_volume_frac(en) = dst%precipitate_volume_frac(en) &
                                      								+ 1.0_pReal/6.0_pReal*PI &
                                      								* (radiusR+ radiusL)**3.0_pReal &
                                      								* (radiusR - radiusL) &
                                      								* stt%precipitate_density(bin,en)


								enddo kwnbins
					! mean radius from m/m^3 to m
      				if (dst%total_precipitate_density(en) > 0.0_pReal) then
      							dst%avg_precipitate_radius(en) = dst%avg_precipitate_radius(en) &
                                     							/ dst%total_precipitate_density(en)
								
	 				endif


					!update matrix composition

					 dst%c_matrix(:,en) = (prm%c0_matrix(:) - dst%precipitate_volume_frac(en)*prm%ceq_precipitate(:))&
	 										/(1-dst%precipitate_volume_frac(en))




   					! print*, ''
    				print*, 'Total precipitate density : ' , dst%total_precipitate_density*1e-18 , '/micron^3'
   					print*, 'Precipitate volume fraction :',  dst%precipitate_volume_frac(en)
    				print*, 'Solute concentration in the matrix' , dst%c_matrix(1,en)
						print*, 'Nucleation rate :part/micron^3/s ', nucleation_rate*1.0e-18
						print*, 'Critical Radius : ', radius_crit*1e9, 'nm'
						

   					! Adapt time step so that the outputs do not vary to much between too time steps
    				!if  either:
    				!    - the precipitate distribution in one class/ the vacancy concentration / the concentration in the matrix becomes negative,
    				! or - at least one size of precipitates grow sufficiently fast to move from two size classes during one time step
    				! then go back to the previous step and decrease the time step


    				if  ( ((stt%c_vacancy(en)) <0.0_pReal)  .OR.(minval(stt%precipitate_density(:,en)) <0.0_pReal).OR.(  minval(dst%c_matrix(:,en))<0.0_pReal).OR.any(isnan(stt%precipitate_density(:,en))))  then
    				! go back one step before

    					stt%time(en)=stt%time(en)-dt
    					dst%total_precipitate_density(en)= temp_total_precipitate_density
						stt%precipitate_density(:,en)=temp_precipitate_density(:)
						dot%precipitate_density(:,en)=temp_dot_precipitate_density(:)
						dst%avg_precipitate_radius(en) = temp_avg_precipitate_radius
       					dst%precipitate_volume_frac(en) = temp_precipitate_volume_frac
        				dst%c_matrix(:,en)  = temp_c_matrix
        				stt%c_vacancy(en) = temp_c_vacancy
        				dislocation_density=temp_dislocation_density
        				T=T_temp
        				radius_crit=temp_radius_crit
        				prm%ceq_matrix= temp_x_eq_matrix
						x_eq_interface = temp_x_eq_interface
						diffusion_coefficient(1)=temp_diffusion_coefficient

       					!decrease time step by a factor arbitrarily chose (2)
        				dt=0.5*dt


    				else
    				!set the time step so that a precipitate cannot grow from more than the space between two adjacent classes

    	    			!not necessary to adapt the time step to the growth rate if there is no precipitate
    	    			if (dst%total_precipitate_density(en)>1.0_pReal) then

    	    				dt=min(dt_max,(prm%bins(1)-prm%bins(0))/maxval(abs(growth_rate_array)))
    	    				print*,'dt growth rate', (prm%bins(1)-prm%bins(0))/maxval(abs(growth_rate_array))

    	    			else
    	    				dt=min(dt_max,dt*1.2) !increase slightly the time step by an arbitrary factor as long as there are no precipitates
						endif


    					! store the new values of the outputs in the temporary variables
    					temp_dot_precipitate_density=dot%precipitate_density(:,en)
    					temp_total_precipitate_density =dst%total_precipitate_density(en)
						temp_precipitate_density = stt%precipitate_density(:,en)
						temp_avg_precipitate_radius = dst%avg_precipitate_radius(en)
        				temp_precipitate_volume_frac = dst%precipitate_volume_frac(en)
        				temp_c_matrix = dst%c_matrix(:,en)
        				temp_radius_crit=radius_crit
        				temp_x_eq_matrix=prm%ceq_matrix
						temp_x_eq_interface=x_eq_interface
						temp_c_vacancy=stt%c_vacancy(en)
						temp_dislocation_density= dislocation_density
        				T_temp=T
        				temp_diffusion_coefficient=diffusion_coefficient(1)


            			if (time_record<stt%time(en)) then !record the outputs every 'time_record' seconds


       						! write outputs in textfiles
     				   		results(1,1)=stt%time(en)
 					   		results(1,2)=dst%avg_precipitate_radius(en)*1.0e9
    				   		results(1,3)=dst%total_precipitate_density(en)*1.0e-18

    				   		if (results(1,3)<1.0e-30_pReal) then
    							results(1,3)=0.0
    						endif

    				  		results(1,4)=dst%precipitate_volume_frac(en)
    				  		if (results(1,4)<1.0e-30_pReal) then
    							results(1,4)=0.0
    				  		endif
   					  		results(1,7:8)=dst%c_matrix(:,en)
   					  		results(1,6)=nucleation_rate*1.0e-18
   							if (results(1,6)<1.0e-30_pReal) then
    							results(1,6)=0.0
    						endif
        			  		results(1,5)=radius_crit*1.0e9

        			  		filename='results/kinetics_data_'
           			  		filename=trim(testfolder)//trim(filename)//trim(filesuffix)

          			  		open(1, file = filename,  ACTION="write", position="append")
    	    	 	  			WRITE(1,13) (results(1,i), i=1,8)
        			  			13 FORMAT(F40.6,F40.6,E40.6,E40.6,E40.6, E40.6, 2E40.6 )
        			  		close(1)

        			 		! writes the current distribution
        					filename='results/precipitation_distribution_'
           					filename=trim(testfolder)//trim(filename)//trim(filesuffix)

        					open(2, file = filename,  ACTION="write", STATUS="replace")
       			 				WRITE(2,'(E40.15)') stt%time(en), stt%precipitate_density(:,en)
        					close(2)


        			    	filename='results/diffusion_coefficient_'
           					filename=trim(testfolder)//trim(filename)//trim(filesuffix)
							open(1, file = filename,  ACTION="write", position="append")
		 						write(1, 601) stt%time(en), diffusion_coefficient(1)
		 						601 FORMAT(2E40.6)
		 					close(1)

		 					filename='results/vacancies_'
		 					filename=trim(testfolder)//trim(filename)//trim(filesuffix)
							open(1, file = filename,  ACTION="write", position="append")
								write(1, 1001) stt%time(en), stt%c_vacancy(en)/c_thermal_vacancy, production_rate/c_thermal_vacancy, annihilation_rate/c_thermal_vacancy
								1001 FORMAT(4E40.6)
							close(1)

							filename='results/dislocation_density_'
           					filename=trim(testfolder)//trim(filename)//trim(filesuffix)
		 					open(1, file = filename,  ACTION="write", position="append")
		 						write(1, 901) stt%time(en), dislocation_density
		 						901 FORMAT(3E40.6)
		 					close(1)


	       					! next time for which the outputs should be written
									if (time_record_step > 0) then
										!Save data linearly
										time_record=time_record+time_record_step

									else
										!save data logarithimically
										time_record =time_record+10**(INT(LOG10(stt%time(en)))-1)

						      endif

    	   				endif
 					endif
	end do loop_time

end program KWN

