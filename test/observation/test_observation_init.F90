integer function test_init_time_movie_observation() bind(C) result(err)
    use FDETYPES
    use FDETYPES_TOOLS
    use Observa
    use observation_testingTools

    type(SGGFDTDINFO) ::  sgg
    type(media_matrices_t) :: media
    type(taglist_t) :: tag_numbers
    logical :: ThereAreObservation, ThereAreWires, ThereAreFarFields
    integer :: initialtimestep
    real(kind=RKIND_tiempo) :: lastexecutedtime
    type(limit_t), dimension(1:6) :: SINPML_fullsize
    type(bounds_t)  ::  bounds
    type(sim_control_t) :: control
    type(nf2ff_t) :: facesNF2FF
    real(kind=RKIND) :: eps = EPSILON_VACUUM, mu = MU_VACUUM
    character(LEN=BUFSIZE)  ::  chari, charj, chark, chari2, charj2, chark2, expectedOutputPath

    type(output_t), pointer, dimension(:) :: output

    sgg = create_base_sgg(dt=0.1_RKIND_tiempo, time_steps=100)
    call set_sgg_data(sgg)

    media = create_media(sgg%Alloc)
    tag_numbers = create_tag_list(sgg%Alloc)

    ThereAreObservation = .false.
    ThereAreWires = .false.
    ThereAreFarFields = .false.

    initialtimestep = 1
    lastexecutedtime = 0.0_RKIND_tiempo

    SINPML_fullsize = create_limit_t(0,4,0,4,0,4,3,3,3)

    facesNF2FF = create_facesNF2FF(.false., .false., .false., .false., .false., .false.)
    control = create_control_flags(nEntradaRoot="entryRoot", wiresflavor="wiresflavour",facesNF2FF=facesNF2FF)

    call InitObservation(sgg, media, tag_numbers, &
                            ThereAreObservation, ThereAreWires, ThereAreFarFields,&
                            initialtimestep, lastexecutedtime, &
                            SINPML_fullsize, eps, mu, bounds, control)
    output => GetOutput()
    
    !Check output status
    err = 0
    if (output(1)%timeswritten /= 0) then 
        print *, "File was written: output(1)%timeswritten=", output(1)%timeswritten
        err = err + 1
    end if
    if (output(1)%item(1)%unit /= 1001) then 
        print *, "Unexpected file unit: (output(1)%item(1)%unit=", output(1)%item(1)%unit
        err = err + 1
    end if

    write (chari, '(i7)') sgg%observation(1)%P(1)%XI
    write (charj, '(i7)') sgg%observation(1)%P(1)%YI
    write (chark, '(i7)') max(sgg%Sweep(fieldo(sgg%observation(1)%P(1)%What, 'Z'))%ZI, sgg%observation(1)%P(1)%ZI)
    write (chari2, '(i7)') sgg%observation(1)%P(1)%XE
    write (charj2, '(i7)') sgg%observation(1)%P(1)%YE
    write (chark2, '(i7)') min(sgg%Sweep(fieldo(sgg%observation(1)%P(1)%What, 'Z'))%ZE, sgg%observation(1)%P(1)%ZE)

    expectedOutputPath = trim(adjustl(control%nEntradaRoot))//'_'//trim(adjustl(sgg%observation(1)%outputrequest))//'_ME_'//&
                            trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'//&
                            trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//&
                            '.bin'

    if (output(1)%item(1)%path /= expectedOutputPath) then
        print *, "Unexpected path for output output(1)%item(1)%path:", output(1)%item(1)%path
        print *, "Expected path:", expectedOutputPath
        err = err + 1
    end if

    contains
        subroutine set_sgg_data(baseSGG)
            type(SGGFDTDINFO), intent(inout) :: baseSGG
            allocate(baseSGG%Observation(1))
            baseSGG%Observation(1) = define_time_movie_observation()
        end subroutine 

        function define_time_movie_observation() result(obs)
            type(Obses_t) :: obs
            obs%nP = 1
            allocate(obs%P(obs%nP))
            obs%P(1) = create_time_movie_observable(1_4,1_4,1_4,6_4,6_4,6_4)

            obs%InitialTime = 0.0_RKIND_tiempo
            obs%FinalTime = 1.0_RKIND_tiempo
            obs%TimeStep = 0.1_RKIND_tiempo

            obs%InitialFreq = 0.0_RKIND
            obs%FinalFreq = 0.0_RKIND
            obs%FreqStep = 0.0_RKIND

            obs%outputrequest = 'timeMovie'

            obs%FreqDomain = .false. 
            obs%TimeDomain = .true. 
            obs%Saveall = .false.
            obs%TransFer = .false.
            obs%Volumic = .true.
            obs%Done = .false.
            obs%Begun = .false.
            obs%Flushed = .false.

        end function define_time_movie_observation

        function create_time_movie_observable(XI,YI,ZI,XE,YE,ZE) result(observable)
            type(observable_t) :: observable
            integer (kind=4)  ::  XI,YI,ZI,XE,YE,ZE

            observable%XI = XI
            observable%YI = YI
            observable%ZI = ZI

            observable%XE = XE
            observable%YE = YE
            observable%ZE = ZE

            observable%Xtrancos = 1
            observable%Ytrancos = 1
            observable%Ztrancos = 1

            observable%What = iMEC
        end function create_time_movie_observable
end function test_init_time_movie_observation