module KWN_data_types

    use KWN_precision
    
    type :: tParameters

        integer :: &
                kwn_nSteps, &           ! discretization in r-space
                incubation              ! for cases where nucleation is considered, set to 1 to consider incubation time
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
                shape_parameter,& ! shape_parameter*mean_radius=standard deviation initial distribution (log normal law assumed)
                volume_fraction_initial, &! initial total precipitate distribution
                rho_0, & !initial dislocation density
                rho_s, & !saturation dislocation density
                strain_rate, & ! strain rate in /s
                Temperature, &
                total_time, & ! heat treatment time for the simulation in s
                dislocation_arrangement, & ! constant related to the dislocation density in the vacancy annihilation term, cf [1]
                burgers, & !matrice burgers vector
                jog_formation_energy, & ! formation energy for jogs
                q_dislocation, & ! activation energy for diffusion at dislocation (pipe diffusion) in J/at - not used yet but to be updated
                solute_strength, & ! constant related to the solid solution hardening- cf [2]
                enthalpy, & ! enthalpy of precipitation 
                entropy, & ! entropy of precipitation
                ! the following parameters are to calculate the flow stress during deformation - there are two options detailed below
                ! option 1 - use a sinepower law that depends on temperature and strain rate for the flow stress, that is then taken as constant during deformation
                sigma_r, & ! constant in the sinepowerlaw for flow stress [MPa]
                A, &  ! constant in the sinepowerlaw for flow stress  [/s]
                Q_stress, &  ! activation energy in the sinepowerlaw for flow stress [J/mol]
                n, & ! stress exponent in the sinepower law for flow stress
                ! option 2 - use an expression for the flow stress that accounts for solid solution hardening, precipitation and dislocations
                k_s, & !constant parameter in regard to solute strength
                k_p, & !constant parameter in regard to precipitate strength
                transition_radius, & ! Transition radius between bypassing and shearing
                M,& ! Taylor Factor
                dt_max,& ! max time step for numerical integration
                time_record_step! defines the frequency for the output files
        

        ! the following variables are allocatable to allow for precipitates with multiple elements (only situations with 2 elements are used here)
        real(pReal), dimension(:), allocatable :: &
                c0_matrix, &            ! initial matrix solute composition in mol fraction : [Mg, Zn]
                ceq_matrix, &           ! equilibrium matrix composition for a flat interface (infinite precipitate radius) in mol fraction : [Mg, Zn]
                ceq_precipitate, &      ! equilibrium precipitate composition in mol fraction : [Mg, Zn]
                diffusion0, &           ! solute diffusivity in m^2/s : [Mg, Zn] - in the present version of the code - the diffusion coefficient is taken as identical for both solute and equal to the slowest diffuser (Mg)
                migration_energy        !  solute migration energy in J/at
        integer, dimension(:), allocatable :: &
                stoechiometry

        real(pReal), dimension(:),   allocatable :: &
                bins                    ! Bins for class sizes in KWN model
        character(len=15), allocatable, dimension(:) :: &
                output
        character*100 :: filesuffix !the file suffix contains the temperature and strain rate used for the simulation
        character*100 :: testfolder !folder where the input file is

    end type tParameters


    type :: tKwnpowerlawState
        real(pReal), pointer, dimension(:,:) :: &
                precipitate_density, & ! table with precipitate density number in each class size [/m^4]
                normalized_distribution_function! table with probability distribution for precipitates in each bin class size 
                
        real(pReal),  dimension(  :), allocatable :: &
                c_vacancy, & ! concentration in excess vacancy
                time, & ! time [s]
                yield_stress, & ! yield stress [MPa]
                growth_rate_array, & ! table with growth rate in each bin
                x_eq_interface ! equilibrium concentration at the interface taking into account Gibbs Thomson effect (one equilibrium concentration for each bin)
        
        real(pReal) :: &
                radius_crit, &
                c_thermal_vacancy, & ! equilibrium concentration of vacancies
                nucleation_rate, & ! part/m^3/s
                dislocation_density ! m^(-2)

    end type tKwnpowerlawState


    type :: tKwnpowerlawMicrostructure
        real(pReal),                  dimension(:),   allocatable :: &
                total_precipitate_density, &
                avg_precipitate_radius, &
                precipitate_volume_frac

        real(pReal), dimension(:,:), allocatable :: &
                c_matrix, &
                diffusion_coefficient

     end type tKwnpowerlawMicrostructure
    


end module KWN_data_types