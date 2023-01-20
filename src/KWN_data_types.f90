module KWN_data_types

    use KWN_precision
    
    type :: tParameters

        integer :: &
                kwn_nSteps              ! discretization in r-space
        real(pReal) :: &
                kwn_stepsize, &         ! discretization in r-space
                kwn_step0               ! minimum radius
        real(pReal) :: &
                lattice_param, &        ! lattice parameter in meter
                atomic_volume, &        ! atomic volume in meter^3
                molar_volume, &         ! molar volume in m^3/mol
                misfit_energy, &        ! normalized precipitate misfit energy in J/m^3
                gamma_coherent, &       ! coherent precipitate surface energy in J/m^2
                vacancy_generation, &   ! vacancy generation rate coefficient
                vacancy_sink_spacing, & ! vacancy sink spacing
                vacancy_energy, &       ! normalized vacancy formation energy (Q/kB) in 1/k
                vacancy_migration_energy, & ! solute migration energy in J/at
                vacancy_diffusion0, &   ! vacancy diffusivity in m^2/s
                mean_radius_initial,&  ! average radius initial distribution in meters
                standard_deviation,& ! standard deviation initial distribution (log normal law assumed)
                volume_fraction_initial, &! initial total precipitate distribution
                rho_0, & !initial dislocation density
                rho_s, & !saturation dislocation density
                strain_rate, & ! strain rate in /s
                dislocation_arrangement, & ! constant related to the dislocation density in the vacancy annihilation term, cf [1]
                burgers, & !matrice burgers vector
                jog_formation_energy, & ! formation energy for jogs
                q_dislocation, & ! activation energy for diffusion at dislocation (pipe diffusion) in J/at - not used yet but to be updated
                solute_strength ! constant related to the solid solution hardening- cf [2]


        ! the following variables are allocatable to allow for precipitates with multiple elements (only situations with 2 elements are used here)
        real(pReal), dimension(:), allocatable :: &
                c0_matrix, &            ! initial matrix solute composition in mol fraction : [Mg, Zn]
                ceq_matrix, &           ! equilibrium matrix composition for a flat interface (infinite precipitate radius) in mol fraction : [Mg, Zn]
                ceq_precipitate, &      ! equilibrium precipitate composition in mol fraction : [Mg, Zn]
                diffusion0, &           ! solute diffusivity in m^2/s : [Mg, Zn] - in the present version of the code - the diffusion coefficient is taken as identical for both solute and equal to the slowest diffuser (Mg)
                migration_energy        !  solute migration energy in J/at


        real(pReal), dimension(:),   allocatable :: &
                bins                    ! Bins for class sizes in KWN model
            character(len=15), allocatable, dimension(:) :: &
                output
    end type tParameters


    type :: tKwnpowerlawState
        real(pReal), pointer, dimension(:,:) :: &
                precipitate_density ! table with precipitate density number in each class size [/m^4]
        real(pReal),  dimension(  :), allocatable :: &
                c_vacancy, & ! concentration in excess vacancy
                time ! time [s]
    end type tKwnpowerlawState


    type :: tKwnpowerlawMicrostructure
        real(pReal),                  dimension(:),   allocatable :: &
                total_precipitate_density, &
                avg_precipitate_radius, &
                precipitate_volume_frac

        real(pReal), dimension(:,:), allocatable :: c_matrix
     end type tKwnpowerlawMicrostructure
    


end module KWN_data_types