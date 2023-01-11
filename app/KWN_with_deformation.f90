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
    use KWN_initialise, only: initialise_model_state
    use KWN_model, only: run_model
    
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

    integer, dimension(:), allocatable :: stoechiometry !precipitate stoechiometry in the following order : Mg Zn Al

    real(pReal), dimension(:,:), allocatable :: &
        normalized_distribution_function !normalised distribution for the precipitate size
    real(pReal), dimension(:), allocatable :: &
        growth_rate_array, &!array that contains the precipitate growth rate of each bin
        x_eq_interface  !array with equilibrium concentrations at the interface between matrix and precipitates of each bin


    real(pReal) :: &
        Temperature, & !temperature in K
        radius_crit, & !critical radius, [m]
        interface_c, & !interface composition between matrix and a precipitate
        time_record_step, & ! time step for the output [s]
        c_thermal_vacancy, & ! concentration in thermal vacancies
        shape_parameter, & !shape parameter in the log normal distribution of the precipitates - ref [4]
        sigma_r, & ! constant in the sinepowerlaw for flow stress [MPa]
        A, &  ! constant in the sinepowerlaw for flow stress  [/s]
        incubation, & ! incubation prefactor either 0 or 1
        Q_stress, &  ! activation energy in the sinepowerlaw for flow stress [J/mol]
        n ! stress exponent in the sinepower law for flow stress

    ! the 'temp' variables are to store the previous step and adapt the time step at each iteration
    real(pReal), dimension(:), allocatable ::   &
        diffusion_coefficient  ! diffusion coefficient for Mg and Zn

    real(pReal) :: &
        dt, & !time step for integration [s]
        dt_max, & ! max time step for integration [s]
        total_time ![s]

    character*100 :: filesuffix !the file suffix contains the temperature and strain rate used for the simulation
    character*100 :: testfolder !folder where the input file is


    call initialise_model_state(prm, dot, stt, dst, &
                                Nmembers, N_elements, en, &
                                stoechiometry, normalized_distribution_function, &
                                Temperature, radius_crit, interface_c, time_record_step, &
                                c_thermal_vacancy, shape_parameter, sigma_r, A, &
                                incubation, Q_stress, n, diffusion_coefficient, &
                                dt, dt_max, total_time, growth_rate_array, &
                                x_eq_interface, &
                                filesuffix, testfolder &
                                )



    call run_model(prm, dot, stt, dst, &
                   Nmembers, N_elements, en, &
                   stoechiometry, normalized_distribution_function, &
                   Temperature, radius_crit, interface_c, time_record_step, &
                   c_thermal_vacancy, shape_parameter, sigma_r, A, &
                   incubation, Q_stress, n, diffusion_coefficient, &
                   dt, dt_max, total_time, growth_rate_array, &
                   x_eq_interface, &
                   filesuffix, testfolder &
                   )


end program KWN

