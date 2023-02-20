The namelist.input options required are described below, with example values.

&config
! output settings
testfolder = 'tests/test_1'             ! path to output folder

! Deformation conditions
Temperature = 423.000000                ! [K] temperature
strain_rate = 1.000e-04                 ! [/s] strain rate - set to zero for static precipitation
incubation = 0.0000e+00                 ! 1 or 0, set to 0 for no precipitate incubation time
total_time = 1.0000e+04                 ! [s] total deformation time 
time_record_step = 1.0000e-01           ! [s] the outputs will be generated for this time step

! Initial distribution
mean_radius_initial = 9.0000e-10        ! [m] initial distribution mean radius
volume_fraction_initial = 7.0000e-03    ! initial distribution volume fraction
shape_parameter = 2.0000e-01            ! s for a log normal distribution

! Composition
c0_matrix = 2.8890e-02, 2.4060e-02      ! bulk composition in [Mg, Zn] in atomic fraction - here A7075
stoechiometry = 5, 7, 6                 ! stoichiometry of the precipitate  Mg, Zn, Al (here it comes from Thermocalc for eta' phase )

! Kinetics model inputs
lattice_param = 4.1e-10                 ! [m] lattice parameter of the product phase
atomic_volume = 1.7e-29                 ! [m^3] atomic volume of matrix
molar_volume = 1.0e-05                  ! [m^3/mol] molar volume of precipitating phase 
misfit_energy = 0.0e+00                 ! [J/m^3] elastic misfit energy
gamma_coherent = 0.265                  ! [J/m^2] coherent surface energy - calibrated with Deschamps 2012 Acta Mat
diffusion0 = 1.4900e-05, 1.4900e-05     ! prefactor diffusion coefficient Mg, Zn [m^2/s]
migration_energy = 1.2650e+05, 1.2650e+05 ! activation energy for solute diffusion (here taken identical for Mg and Zn)
ceq_matrix = 1.2200e-02, 5.0400e-04     ! equilibrium concentration in Mg and Zn in the matrix
!!! alternatively, ceq_matrix can be omitted and replaced by value for enthalpy and entropy of precipitation - in that case ce_matrix will be recalculated
 enthalpy = 290e+03 
 entropy = 62.0e+00

! Vacancy model inputs
vacancy_sink_spacing = 5.000e-05        ! sink spacing, here grain size [m]
vacancy_diffusion0 = 1.000e-05          ! prefactor diffusion coefficient for vacancies - set to zero if there is no deformation
jog_formation_energy = 3.000e-01        ! [eV]
q_dislocation = 1.083e+05               ! [J/mol] -ref Legros et al. 2008
vacancy_energy = 5.200e-01              ! [eV]
vacancy_migration_energy= 9.300e-01     ! [eV]
dislocation_arrangement = 1.000e+01     ! empirical factor from 1 to 10 in Militzer model
vacancy_generation = 3.500e-02          ! mechanical vacancy production constant in Militzer model -  - set to zero if there is no deformation
rho_0 = 1.000e+14                       ! [/m^2] initial dislocation density
rho_s = 1.000e+14                       ! [/m^2] saturation dislocation density
burgers = 2.9000e-10                    ! [m] burgers vector

! sinepower law for stress
sigma_r = 1.0000e+08                    ! [MPa]
A = 5.2140e-06                          ! [/s]
Q_stress = 6.0000e+04                   ! [J/mol]
n = 6.6831e+00                          !

! KWN discretisation
kwn_step0 = -9.300e+00                  ! [m] starting bin radius in logscale
kwn_stepsize = 1.000e-02                ! [m] spacing between bins in logscale
kwn_nsteps = 100                        ! no. of radius bins
dt_max = 1.0000e-02                     ! [s] maximum time increment

/
