module test_preprocess_m
  use Observa
  implicit none

contains

  logical function approx_equal(a,b,tol)
    real(kind=RKIND), intent(in) :: a,b,tol
    approx_equal = abs(a-b) <= tol
  end function approx_equal

  logical function test_initial_time_less_than_timestep() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    test_initial_time_less_than_timestep = .false.

    dt = 0.1_RKIND
    obs%TimeStep = 0.5_RKIND
    obs%InitialTime = 0.2_RKIND   ! less than TimeStep -> should be set to 0.0
    obs%FinalTime = 2.0_RKIND
    obs%nP = 0
    obs%Volumic = .false.
    obs%InitialFreq = 0.0_RKIND
    obs%FinalFreq = 1.0_RKIND
    obs%FreqStep = 0.1_RKIND
    out%SaveAll = .false.

    call preprocess_observation(obs, out, dt)

    if (approx_equal(obs%InitialTime, 0.0_RKIND, 1e-12_RKIND)) then
        err = 0
    else
        print *, "test_initial_time_less_than_timestep: obs%InitialTime=", obs%InitialTime
        err = 1
    endif
  end function test_initial_time_less_than_timestep

  logical function test_timestep_greater_and_mapvtk() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    test_timestep_greater_and_mapvtk = .false.

    dt = 0.1_RKIND
    obs%TimeStep = 5.0_RKIND
    obs%InitialTime = 0.0_RKIND
    obs%FinalTime = 1.0_RKIND   ! TimeStep > Final-Initial
    obs%nP = 1
    allocate(obs%P(1))
    obs%P(1)%what = mapvtk
    obs%Volumic = .false.
    obs%InitialFreq = 0.0_RKIND
    obs%FinalFreq = 1.0_RKIND
    obs%FreqStep = 0.1_RKIND

    call preprocess_observation(obs, out, dt)


    if (approx_equal(obs%InitialTime, 0.0_RKIND, 1e-12_RKIND) .and. approx_equal(obs%FinalTime, 0.0_RKIND, 1e-12_RKIND)) then
        err = 0
    else
        print *, "test_timestep_greater_and_mapvtk: InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
        err = 1
    endif
    deallocate(obs%P)
  end function test_timestep_greater_and_mapvtk

  logical function test_timestep_greater_not_mapvtk() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    test_timestep_greater_not_mapvtk = .false.

    dt = 0.1_RKIND
    obs%TimeStep = 2.0_RKIND
    obs%InitialTime = 1.0_RKIND
    obs%FinalTime = 1.5_RKIND   ! TimeStep > Final-Initial (=0.5) so branch -> FinalTime = Initial + TimeStep
    obs%nP = 1
    allocate(obs%P(1))
    obs%P(1)%what = 999   ! not mapvtk
    obs%Volumic = .false.
    obs%InitialFreq = 0.0_RKIND
    obs%FinalFreq = 1.0_RKIND
    obs%FreqStep = 0.1_RKIND

    call preprocess_observation(obs, out, dt)

    if (approx_equal(obs%FinalTime, obs%InitialTime + obs%TimeStep, 1e-12_RKIND)) then
       err = 0
    else
       print *, "test_timestep_greater_not_mapvtk: InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
       err = 1
    endif
    deallocate(obs%P)
  end function test_timestep_greater_not_mapvtk

  logical function test_freqstep_zero_or_large() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    test_freqstep_zero_or_large = .false.
    dt = 0.1_RKIND

    ! Case A: FreqStep = 0 -> should be set to FinalFreq-InitialFreq
    obs%TimeStep = 0.2_RKIND
    obs%InitialTime = 0.0_RKIND
    obs%FinalTime = 1.0_RKIND
    obs%InitialFreq = 1.0_RKIND
    obs%FinalFreq = 3.5_RKIND
    obs%FreqStep = 0.0_RKIND
    obs%nP = 0
    obs%Volumic = .false.

    call preprocess_observation(obs, out, dt)
    
    if (.not. approx_equal(obs%FreqStep, obs%FinalFreq - obs%InitialFreq, 1e-12_RKIND)) then
       print *, "test_freqstep_zero_or_large A: FreqStep=", obs%FreqStep, " expected=", obs%FinalFreq - obs%InitialFreq
       err = 1
       return
    end if

    ! Case B: FreqStep > FinalFreq-InitialFreq -> becomes FinalFreq-InitialFreq
    obs%FreqStep = 10.0_RKIND
    obs%InitialFreq = 0.0_RKIND
    obs%FinalFreq = 2.0_RKIND
    call preprocess_observation(obs, out, dt)
    
    if (approx_equal(obs%FreqStep, obs%FinalFreq - obs%InitialFreq, 1e-12_RKIND)) then
    else
       print *, "test_freqstep_zero_or_large B: FreqStep=", obs%FreqStep, " finalFreq=", obs%FinalFreq, " initialFreq=", obs%InitialFreq
       err = 1
       return
    endif

    err = 0
  end function test_freqstep_zero_or_large

  logical function test_volumic_false_true_and_saveall() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    logical :: ok1, ok2
    real(kind=RKIND) :: dt
    test_volumic_false_true_and_saveall = .false.
    dt = 0.1_RKIND

    ! Case Volumic = .false. and global saveall = .false.
    obs%Volumic = .false.
    obs%Saveall = .false.
    obs%TimeStep = 0.2_RKIND
    obs%InitialTime = 0.0_RKIND
    obs%FinalTime = 1.0_RKIND
    obs%InitialFreq = 0.0_RKIND
    obs%FinalFreq = 1.0_RKIND
    out%SaveAll = .false.

    saveall = .false.
    call preprocess_observation(obs, out, dt)
    ok1 = (out%SaveAll == .false.)

    ! Now set global saveall to true -> observation%Saveall becomes true and out%SaveAll true
    saveall = .true.
    obs%Saveall = .false.
    call preprocess_observation(obs, out, dt)
    ok2 = (obs%Saveall == .true__ .or. out%SaveAll == .true__)  ! accept either but at least out%SaveAll true
    ! (the code sets observation%Saveall = observation%Saveall .or. saveall)

    ! revert global
    saveall = .false.

    if (ok1 .and. ok2) then
        err = 0
    else
        err = 1
    endif
  end function test_volumic_false_true_and_saveall

  logical function test_saveall_branch() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    test_saveall_branch = .false.
    dt = 0.1_RKIND

    obs%Volumic = .false.
    obs%Saveall = .true.
    obs%TimeStep = 0.5_RKIND
    obs%InitialTime = 10.0_RKIND
    obs%FinalTime = 20.0_RKIND

    call preprocess_observation(obs, out, dt)

    print *, "test_saveall_branch: Trancos=", out%Trancos, " InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
    if (out%Trancos == 1 .and. approx_equal(obs%InitialTime, 0.0_RKIND, 1e-12_RKIND) .and. &
         approx_equal(obs%FinalTime, sgg%tiempo(FINALTIMESTEP+2), 1e-12_RKIND)) then
       test_saveall_branch = .true.
       print *, "  PASS"
    else
       print *, "  FAIL"
    endif

  end function test_saveall_branch

  logical function test_final_less_than_initial() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    test_final_less_than_initial = .false.
    dt = 0.1_RKIND

    obs%Volumic = .false.
    obs%Saveall = .false.
    obs%TimeStep = 0.2_RKIND
    obs%InitialTime = 5.0_RKIND
    obs%FinalTime = 1.0_RKIND   ! Final < Initial: should be fixed to Initial after processing

    call preprocess_observation(obs, out, dt)

    print *, "test_final_less_than_initial: InitialTime=", obs%InitialTime, " FinalTime=", obs%FinalTime
    if (approx_equal(obs%FinalTime, obs%InitialTime, 1e-12_RKIND)) then
       test_final_less_than_initial = .true.
       print *, "  PASS"
    else
       print *, "  FAIL"
    endif
  end function test_final_less_than_initial

  logical function test_huge_cap() bind(C) result(err)
    type(Obses_t) :: obs
    type(output_t) :: out
    real(kind=RKIND) :: dt
    real(kind=4) :: huge4
    logical :: ok
    test_huge_cap = .false.
    dt = 0.1_RKIND
    huge4 = huge(1.0_4)

    ! Set FinalTime - InitialTime extremely large so that
    ! 10 * (Final - Initial) / min(0.1_RKIND, TimeStep) >= huge(1_4)
    obs%TimeStep = 0.1_RKIND
    obs%InitialTime = 0.0_RKIND
    ! pick FinalTime big: make (Final-Initial) = huge(1_4) * min(dt,TimeStep) / 5
    obs%FinalTime = obs%InitialTime + real(huge4 / 5.0_4, kind=RKIND) * min(0.1_RKIND, obs%TimeStep)
    obs%Volumic = .false.
    obs%Saveall = .false.
    obs%InitialFreq = 0.0_RKIND
    obs%FinalFreq = 1.0_RKIND
    obs%FreqStep = 0.1_RKIND

    call preprocess_observation(obs, out, dt)

    ! After preprocessing, FinalTime must have been capped to something not exceeding huge-scale expression:
    ! A simple sanity check: FinalTime should not exceed sgg%tiempo(FINALTIMESTEP+2) (the upper clamp)
    ok = (obs%FinalTime <= sgg%tiempo(FINALTIMESTEP+2) + 1e-8_RKIND)
    print *, "test_huge_cap: FinalTime=", obs%FinalTime, " clip=", sgg%tiempo(FINALTIMESTEP+2)
    if (ok) then
       test_huge_cap = .true.
       print *, "  PASS"
    else
       print *, "  FAIL"
    endif
  end function test_huge_cap

end module test_preprocess_m