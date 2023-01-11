      INTEGER :: unit
      OPEN(UNIT=unit, FILE='input.nml', STATUS='OLD', ACTION='READ')
      INQUIRE(FILE='input.nml',EXIST=exist)
      IF (.not.exist) THEN
          PRINT *,'The namelist file does not exist'
          STOP
      ELSE
          READ(unit, '(A)',ERR=10) line
          READ(line, '(A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4)',ERR=20) &
          T,strain_rate,incubation,time,output_time_step
          READ(line, '(A,F7.4,A,F7.4,A,F7.4)',ERR=20)  initial_mean_radius,initial_volume_fraction,dispersion_parameter
          READ(line, '(A,2F7.4,A,3I3)',ERR=20)  c0_matrix,stoichiometry
          READ(line, '(A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,2F7.4,A,2F7.4)',ERR=20)  &
          lattice_parameter,atomic_volume,molar_volume,misfit_energy,gamma_coherent,solute_diffusion0,solute_migration_energy,c_eq
          READ(line, '(A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,F7.4,A,I3,A,F7.4,A,F7.4,A,F7.4,A,F7.4)',ERR=20) &
          vacancy_sink_spacing,vacancy_diffusion0,jog_formation_energy,activation_energy_pipe_diffusion,vacancy_formation_energy,vacancy_migration_energy,dislocation_arrangement,vacancy_generation,initial_dislocation_density,saturation_dislocation_density,burgers_vector
          READ(line, '(A,F7.4,A,F7.4,A,F7.4,A,F7.4)',ERR=20)  sigma_r,A,Q_stress,n
          READ(line, '(A,F7.4,A,F7.4,A,I4,A,F7.4)',ERR=20)  kwn_step0,kwn_stepsize,kwn_nsteps,dt_max    
      END IF
      CLOSE(unit)
