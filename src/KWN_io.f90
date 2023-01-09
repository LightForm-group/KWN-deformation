module KWN_io

    use KWN_precision
    use KWN_data_types, only: tParameters, tKwnpowerlawState
    
contains
    
subroutine read_configuration( &
							testfolder, &
							prm, &
							T, &
							shape_parameter, &  ! the initial distribution is defined by mean radius, volume fraction and shape parameter of a log normal distribution - see e.g. ref [4]
							total_time, &  ![s]
							dt_max, &  ![s]
							time_record_step, &  ![s]
							sigma_r, &  ![MPa] - sinepower law for stress
							A, & ![/s] - sinepower law for stress
							Q_stress, &  ![J/mol] - activation energy in flow stress law
							n, &  !exponent in sinepower law for stress
							incubation,  & !incubation prefactor, either 0 or 1)
							stoechiometry &
							)
    

    implicit none

    character*100, intent(out) :: testfolder
    type(tParameters), intent(inout) :: prm
	integer, dimension(:), allocatable, intent(inout) :: stoechiometry !precipitate stoechiometry in the following order : Mg Zn Al
    real(pReal), intent(out) :: &
        T, &
        shape_parameter, &
        total_time, &
        dt_max, &
        time_record_step, &
        sigma_r, &
        A, &
        Q_stress, &
        n, &
        incubation
    INTEGER :: status ! I/O status


	! Read the inputs from the input.dat fil
	!!!!!!!!!!!
	OPEN (UNIT=1, FILE='input.dat', STATUS='OLD', ACTION='READ', IOSTAT=status)
	print*, status
	READ(1,*) testfolder
	READ(1,*, IOSTAT=status) prm%kwn_step0! starting bin radius  (m)
	print*, status
	READ(1,*) prm%kwn_stepsize ! spacing between bins (m)
	READ(1,*) prm%kwn_nSteps ! no. of radius bins
	READ(1,*) T
	READ(1,*) prm%strain_rate
	READ(1,*) prm%lattice_param  ! lattice parameter of the product phase (m)
	READ(1,*) prm%atomic_volume ! atomic volume of matrix
	READ(1,*) prm%molar_volume  ! molar volume of precipitating phase
	READ(1,*) prm%misfit_energy   ! elastic strain energy of the precipitate J/m^3
	READ(1,*) prm%gamma_coherent ! interfacial energy  J/m^2
	READ(1,*) prm%diffusion0(1) ! m^2/s
	READ(1,*) prm%diffusion0(2)! m^2/s
	READ(1,*) prm%migration_energy(1) !J/mol
	READ(1,*) prm%migration_energy(2) !J/mol
	READ(1,*) prm%c0_matrix(1) !at fraction
	READ(1,*) prm%c0_matrix(2) ! at fraction
	READ(1,*) stoechiometry(1)
	READ(1,*) stoechiometry(2)
	READ(1,*) stoechiometry(3)
	READ(1,*) prm%vacancy_sink_spacing ![m]
	READ(1,*) prm%vacancy_diffusion0 ![m^2/s]
	READ(1,*) prm%jog_formation_energy ![ev]
	READ(1,*) prm%q_dislocation ! [J/mol]
	READ(1,*) prm%vacancy_energy ! [ev]
	READ(1,*) prm%vacancy_migration_energy ![ev]
	READ(1,*) prm%dislocation_arrangement ![]
	READ(1,*) prm%vacancy_generation ![]
	READ(1,*) prm%rho_0 ![/m^2]
	READ(1,*) prm%rho_s ![/m^2]
	READ(1,*) prm%burgers ! [m]
	READ(1,*) prm%mean_radius_initial ![m]
	READ(1,*) prm%volume_fraction_initial ![]
	READ(1,*) shape_parameter ! the initial distribution is defined by mean radius, volume fraction and shape parameter of a log normal distribution - see e.g. ref [4]
	READ(1,*) total_time ![s]
	READ(1,*) dt_max ![s]
	READ(1,*) time_record_step ![s]
	READ(1,*) sigma_r ![MPa] - sinepower law for stress
	READ(1,*) A ![/s] - sinepower law for stress
	READ(1,*) Q_stress ![J/mol] - activation energy in flow stress law
	READ(1,*) n !exponent in sinepower law for stress
	READ(1,*) prm%ceq_matrix(1) !equilibrium concentration in the matrix
	READ(1,*) prm%ceq_matrix(2) !equilibrium concentration in the matrix
	READ(1,*) incubation !incubation prefactor, either 0 or 1
	CLOSE(1)

end subroutine read_configuration


subroutine output_results(testfolder, filesuffix, stt, diffusion_coefficient, c_thermal_vacancy, &
                        production_rate, annihilation_rate, dislocation_density, en)

	type(tKwnpowerlawState), intent(in) :: stt

	real(pReal), dimension(:), allocatable, intent(in) ::   &
		diffusion_coefficient  ! diffusion coefficient for Mg and Zn

    real(pReal), intent(in) :: &
		c_thermal_vacancy, & ! concentration in thermal vacancies
		production_rate, & ! production rate for excess vacancies
		annihilation_rate, & !annihilation rate for excess vacancies
		dislocation_density ![/m^2]
    
    integer, intent(in) :: en
    
	character*100, intent(in) :: filesuffix !the file suffix contains the temperature and strain rate used for the simulation
	character*100, intent(in) :: testfolder !folder where the input file is

    ! local variables
	character*100 :: filename !name of the gile where the outputs will be written

	filename = 'results/diffusion_coefficient_'
	filename = trim(testfolder)//trim(filename)//trim(filesuffix)

	open(1, file = filename,  ACTION="write", position="append")
		write(1, 601) stt%time(en), diffusion_coefficient(1)
	close(1)

	filename = 'results/vacancies_'
	filename = trim(testfolder)//trim(filename)//trim(filesuffix)
	open(1, file = filename,  ACTION="write", position="append")
		write(1, 1001) stt%time(en), stt%c_vacancy(en)/c_thermal_vacancy, &
		               production_rate/c_thermal_vacancy, annihilation_rate/c_thermal_vacancy
	close(1)

	filename = 'results/dislocation_density_'
	filename = trim(testfolder)//trim(filename)//trim(filesuffix)
	open(1, file = filename,  ACTION="write", position="append")
		write(1, 901) stt%time(en), dislocation_density
	close(1)

601 FORMAT(2E40.6)
901 FORMAT(3E40.6)
1001 FORMAT(4E40.6)

end subroutine output_results

end module KWN_io