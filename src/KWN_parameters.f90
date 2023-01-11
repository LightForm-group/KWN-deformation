module KWN_parameters

    use KWN_precision
    
    real(pReal), parameter :: &

        kB = 1.38064852e-23_pReal, &
        R  = 8.3145_pReal, &
        PI = 3.14159265359_pReal, &
        na = 6.02214076e23_pReal

    real, parameter :: &
        ev_to_Jat = 1.602176634e-19 ! convert from ev to  J/at

end module KWN_parameters