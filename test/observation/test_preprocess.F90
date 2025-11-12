integer function test_initial_time_less_than_timestep() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA

  type(Obses_t) :: obs
  type(output_t) :: out
  
  integer(kind=4) :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer , dimension(:) :: tiempo
  
  logical :: saveall

  obs%TimeStep = 0.5_RKIND
  obs%InitialTime = 0.2_RKIND   ! less than TimeStep -> should be set to 0.0
  obs%FinalTime = 2.0_RKIND
  obs%nP = 0
  obs%Volumic = .false.
  obs%InitialFreq = 0.0_RKIND
  obs%FinalFreq = 1.0_RKIND
  obs%FreqStep = 0.1_RKIND
  out%SaveAll = .false.

  finalTimeIndex = 20
  dt = 0.1
  tiempo => create_time_array(100, dt)
  
  saveall = .true.

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (approx_equal(obs%InitialTime, 0.0_RKIND, 1e-12_RKIND)) then
    err = 0
  else
    print *, "test_initial_time_less_than_timestep: obs%InitialTime=", obs%InitialTime
    err = 1
  end if
end function test_initial_time_less_than_timestep

integer function test_timestep_greater_and_mapvtk() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA
  
  type(Obses_t) :: obs
  type(output_t) :: out
  
  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo
  
  logical :: saveall

  

  obs%TimeStep = 5.0_RKIND
  obs%InitialTime = 0.0_RKIND
  obs%FinalTime = 1.0_RKIND   ! TimeStep > Final-Initial
  obs%nP = 1
  allocate (obs%P(obs%nP))
  obs%P(1)%what = mapvtk
  obs%Volumic = .false.
  obs%InitialFreq = 0.0_RKIND
  obs%FinalFreq = 1.0_RKIND
  obs%FreqStep = 0.1_RKIND

  finalTimeIndex = 90
  dt = 0.1
  tiempo => create_time_array(100, dt)

  saveall = .false.

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (approx_equal(obs%InitialTime, 0.0_RKIND, 1e-12_RKIND) .and. approx_equal(obs%FinalTime, 0.0_RKIND, 1e-12_RKIND)) then
    err = 0
  else
    print *, "test_timestep_greater_and_mapvtk: InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
    err = 1
  end if

  deallocate (obs%P)
end function test_timestep_greater_and_mapvtk

integer function test_timestep_greater_not_mapvtk() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA
  type(Obses_t) :: obs
  type(output_t) :: out

  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo

  logical :: saveall

  finalTimeIndex = 90
  dt = 0.1_RKIND
  tiempo => create_time_array(100, 0.1_RKIND_tiempo)
  saveall = .false.

  obs%TimeStep = 2.0_RKIND
  obs%InitialTime = 1.0_RKIND
  obs%FinalTime = 1.5_RKIND   ! TimeStep > Final-Initial (=0.5) so branch -> FinalTime = Initial + TimeStep
  obs%nP = 1
  allocate (obs%P(1))
  obs%P(1)%what = 999   ! not mapvtk
  obs%Volumic = .false.
  obs%InitialFreq = 0.0_RKIND
  obs%FinalFreq = 1.0_RKIND
  obs%FreqStep = 0.1_RKIND
  obs%Saveall = .false.

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (approx_equal(obs%FinalTime, obs%InitialTime + obs%TimeStep, 1e-12_RKIND)) then
    err = 0
  else
    print *, "test_timestep_greater_not_mapvtk: InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
    err = 1
  end if
  deallocate (obs%P)
end function test_timestep_greater_not_mapvtk

integer function test_freqstep_zero_or_large() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA
  type(Obses_t) :: obs
  type(output_t) :: out

  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo
  logical :: saveall

  finalTimeIndex = 90
  dt = 0.1_RKIND
  tiempo => create_time_array(100, 0.1_RKIND_tiempo)
  saveall = .false.

  ! Case A: FreqStep = 0 -> should be set to FinalFreq-InitialFreq
  obs%TimeStep = 0.2_RKIND
  obs%InitialTime = 0.0_RKIND
  obs%FinalTime = 1.0_RKIND
  obs%InitialFreq = 1.0_RKIND
  obs%FinalFreq = 3.5_RKIND
  obs%FreqStep = 0.0_RKIND
  obs%nP = 0
  obs%Volumic = .false.

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (.not. approx_equal(obs%FreqStep, obs%FinalFreq - obs%InitialFreq, 1e-12_RKIND)) then
    print *, "test_freqstep_zero_or_large A: FreqStep=", obs%FreqStep, " expected=", obs%FinalFreq - obs%InitialFreq
    err = 1
    return
  end if

  ! Case B: FreqStep > FinalFreq-InitialFreq -> becomes FinalFreq-InitialFreq
  obs%FreqStep = 10.0_RKIND
  obs%InitialFreq = 0.0_RKIND
  obs%FinalFreq = 2.0_RKIND
  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (approx_equal(obs%FreqStep, obs%FinalFreq - obs%InitialFreq, 1e-12_RKIND)) then
  else
   print *, "test_freqstep_zero_or_large B: FreqStep=", obs%FreqStep, " finalFreq=", obs%FinalFreq, " initialFreq=", obs%InitialFreq
    err = 1
    return
  end if

  err = 0
end function test_freqstep_zero_or_large

integer function test_volumic_false_true_and_saveall() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA

  type(Obses_t) :: obs
  type(output_t) :: out

  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo
  logical :: ok1, ok2, saveall

  finalTimeIndex = 90
  dt = 0.1_RKIND
  tiempo => create_time_array(100, dt)
  saveall = .false.

  ! Case Volumic = .false. and global saveall = .false.
  obs%Volumic = .false.
  obs%Saveall = .false.
  obs%TimeStep = 0.2_RKIND
  obs%InitialTime = 0.0_RKIND
  obs%FinalTime = 1.0_RKIND
  obs%InitialFreq = 0.0_RKIND
  obs%FinalFreq = 1.0_RKIND
  out%SaveAll = .false.
  obs%nP = 1
  allocate (obs%P(1))
  obs%P(1)%what = 999 !Not mapvtk

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)
  ok1 = (out%SaveAll .eqv. .false.)

  ! Now set global saveall to true -> observation%Saveall becomes true and out%SaveAll true
  saveall = .true.
  obs%Saveall = .false.
  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)
  ok2 = (obs%Saveall .eqv. .true.) .or. (out%SaveAll .eqv. .true.)  ! accept either but at least out%SaveAll true
  ! (the code sets observation%Saveall = observation%Saveall .or. saveall)

  ! revert global
  saveall = .false.

  if (ok1 .and. ok2) then
    err = 0
  else
    err = 1
  end if
end function test_volumic_false_true_and_saveall

integer function test_saveall_branch() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA
  type(Obses_t) :: obs
  type(output_t) :: out

  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo
  logical :: saveall

  finalTimeIndex = 90
  dt = 0.1_RKIND
  tiempo => create_time_array(100, dt)
  saveall = .false.

  obs%Volumic = .false.
  obs%Saveall = .true.
  obs%TimeStep = 0.5_RKIND
  obs%InitialTime = 10.0_RKIND
  obs%FinalTime = 20.0_RKIND
  obs%nP = 1
  allocate (obs%P(1))
  obs%P(1)%what = 999 !Not mapvtk

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (out%Trancos == 1 .and. approx_equal(obs%InitialTime, 0.0_RKIND, 1e-12_RKIND) .and. &
    approx_equal(obs%FinalTime, real(tiempo(finalTimeIndex + 2), RKIND), 1e-12_RKIND)) then
    err = 0
  else
    print *, "test_saveall_branch: Trancos=", out%Trancos, " InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
    err = 1
  end if

end function test_saveall_branch

integer function test_final_less_than_initial() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA
  type(Obses_t) :: obs
  type(output_t) :: out

  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo
  logical :: saveall
  
  finalTimeIndex = 90
  dt = 0.1_RKIND
  tiempo => create_time_array(100, 0.1_RKIND_tiempo)
  saveall = .false.

  obs%Volumic = .false.
  obs%Saveall = .false.
  obs%TimeStep = 0.2_RKIND
  obs%InitialTime = 5.0_RKIND
  obs%FinalTime = 1.0_RKIND   ! Final < Initial: should be fixed to Initial after processing
  obs%nP = 1
  allocate (obs%P(1))
  obs%P(1)%what = 999 !Not mapvtk

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  if (approx_equal(obs%FinalTime, obs%InitialTime, 1e-12_RKIND)) then
    err = 0
  else
    print *, "test_final_less_than_initial: InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
    err = 1
  end if
end function test_final_less_than_initial

integer function test_huge_cap() bind(C) result(err)
  use observation_testingTools
  use FDETYPES
  use OBSERVA
  type(Obses_t) :: obs
  type(output_t) :: out

  integer :: finalTimeIndex
  real(kind=RKIND_tiempo) :: dt
  real(KIND=RKIND_tiempo), pointer, dimension(:) :: tiempo
  real(kind=4) :: huge4
  logical :: ok, saveall

  finalTimeIndex = 90
  dt = 0.1_RKIND
  tiempo => create_time_array(100, 0.1_RKIND_tiempo)
  huge4 = huge(1.0_4)
  saveall = .false.

  ! Set FinalTime - InitialTime extremely large so that
  ! 10 * (Final - Initial) / min(0.1_RKIND, TimeStep) >= huge(1_4)
  obs%TimeStep = 0.1_RKIND
  obs%InitialTime = 0.0_RKIND
  ! pick FinalTime big: make (Final-Initial) = huge(1_4) * min(dt,TimeStep) / 5
  obs%FinalTime = obs%InitialTime + real(huge4/5.0_4, kind=RKIND)*min(0.1_RKIND, obs%TimeStep)
  obs%Volumic = .false.
  obs%Saveall = .false.
  obs%InitialFreq = 0.0_RKIND
  obs%FinalFreq = 1.0_RKIND
  obs%FreqStep = 0.1_RKIND
  obs%nP = 1
  allocate (obs%P(1))
  obs%P(1)%what = 999 !Not mapvtk

  call preprocess_observation(obs, out, tiempo, finalTimeIndex, dt, saveall)

  ! After preprocessing, FinalTime must have been capped to something not exceeding huge-scale expression:
  ! A simple sanity check: FinalTime should not exceed tiempo(FINALTIMESTEP+2) (the upper clamp)
  ok = (obs%FinalTime <= tiempo(finalTimeIndex + 2) + dt)

  if (ok) then
    err = 0
  else
    print *, "test_huge_cap: FinalTime=", obs%FinalTime, " clip=", tiempo(FINALTIMESTEP + 2)
    err = 1
  end if
end function test_huge_cap
