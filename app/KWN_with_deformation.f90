!--------------------------------------------------------------------------------------------------
!> @author Madeleine Bignon, University of Manchester
!> @author Pratheek Shanthraj, University of Manchester
! KWN precipitation model including the effect of deformation via excess vacancy concentration - model described in [3]
! Bignon, M., Shanthraj, P., & Robson, J. D. (2022). Acta Materialia, 234, 118036. https://doi.org/10.1016/J.ACTAMAT.2022.118036
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
        en

    real(pReal), dimension(:), allocatable :: &
        x_eq_interface  !array with equilibrium concentrations at the interface between matrix and precipitates of each bin


    real(pReal) :: &
        interface_c!interface composition between matrix and a precipitate

    real(pReal) :: &
        dt !time step for integration [s]


    call initialise_model_state(prm, dot, stt, dst, &
                                Nmembers, en, &
                                interface_c, &
                                dt, &
                                x_eq_interface &
                                )



    call run_model(prm, dot, stt, dst, &
                   Nmembers, en, &
                   interface_c,  &
                   dt,  &
                   x_eq_interface &
                   )


end program KWN

