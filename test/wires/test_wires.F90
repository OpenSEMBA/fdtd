integer function test_getHwires_returns_associated() bind(C, name="test_getHwires_returns_associated") result(status)
    use HollandWires_m, only: GetHwires
    use wiresHolland_constants_m, only: ThinWires_t
    implicit none
    type(ThinWires_t), pointer :: p

    status = 0
    p => GetHwires()
    if (.not. associated(p)) then
        print *, "test_getHwires_returns_associated FAILED: returned pointer is not associated"
        status = 1
    end if
end function test_getHwires_returns_associated

integer function test_destroyWires_deallocates_wires_data() bind(C, name="test_destroyWires_deallocates_wires_data") result(status)
    use HollandWires_m, only: DestroyWires, GetHwires
    use FDETYPES_m, only: SGGFDTDINFO_t, MediaData_t
    use wiresHolland_constants_m, only: ThinWires_t
    implicit none
    type(SGGFDTDINFO_t) :: sgg
    type(MediaData_t) :: basic_media
    type(ThinWires_t), pointer :: p

    status = 0

    sgg%NumMedia = 1
    allocate(sgg%Med(0:sgg%NumMedia))
    sgg%Med = basic_media

    sgg%Med(1)%Is%ThinWire = .true.
    allocate(sgg%Med(1)%wire(1))
    allocate(sgg%Med(1)%wire(1)%Vsource(1))
    allocate(sgg%Med(1)%wire(1)%Isource(1))

    p => GetHwires()
    allocate(p%WireTipoMedio(1))
    allocate(p%CurrentSegment(1))
    allocate(p%ChargeNode(1))

    call DestroyWires(sgg)

    if (associated(sgg%Med(1)%wire)) then
        print *, "test_destroyWires_deallocates_wires_data FAILED: media wire pointer still associated"
        status = 1
    end if
    if (associated(p%WireTipoMedio)) then
        print *, "test_destroyWires_deallocates_wires_data FAILED: HWires%WireTipoMedio still associated"
        status = 1
    end if
    if (associated(p%CurrentSegment)) then
        print *, "test_destroyWires_deallocates_wires_data FAILED: HWires%CurrentSegment still associated"
        status = 1
    end if
    if (associated(p%ChargeNode)) then
        print *, "test_destroyWires_deallocates_wires_data FAILED: HWires%ChargeNode still associated"
        status = 1
    end if

    if (associated(sgg%Med)) deallocate(sgg%Med)
end function test_destroyWires_deallocates_wires_data

integer function test_destroyWires_keeps_non_thinwire_media() bind(C, name="test_destroyWires_keeps_non_thinwire_media") result(status)
    use HollandWires_m, only: DestroyWires
    use FDETYPES_m, only: SGGFDTDINFO_t, MediaData_t
    implicit none
    type(SGGFDTDINFO_t) :: sgg
    type(MediaData_t) :: basic_media

    status = 0

    sgg%NumMedia = 1
    allocate(sgg%Med(0:sgg%NumMedia))
    sgg%Med = basic_media

    sgg%Med(1)%Is%ThinWire = .false.
    allocate(sgg%Med(1)%wire(1))
    allocate(sgg%Med(1)%wire(1)%Vsource(1))
    allocate(sgg%Med(1)%wire(1)%Isource(1))

    call DestroyWires(sgg)

    if (.not. associated(sgg%Med(1)%wire)) then
        print *, "test_destroyWires_keeps_non_thinwire_media FAILED: non-thinwire media was deallocated"
        status = 1
    end if

    if (associated(sgg%Med(1)%wire(1)%Vsource)) deallocate(sgg%Med(1)%wire(1)%Vsource)
    if (associated(sgg%Med(1)%wire(1)%Isource)) deallocate(sgg%Med(1)%wire(1)%Isource)
    if (associated(sgg%Med(1)%wire)) deallocate(sgg%Med(1)%wire)
    if (associated(sgg%Med)) deallocate(sgg%Med)
end function test_destroyWires_keeps_non_thinwire_media

integer function test_destroyWireMedia_thinwire() bind(C, name="test_destroyWireMedia_thinwire") result(status)
    use HollandWires_m, only: DestroyWireMedia
    use FDETYPES_m, only: MediaData_t
    implicit none
    type(MediaData_t) :: media

    status = 0

    media%Is%ThinWire = .true.
    allocate(media%wire(1))
    allocate(media%wire(1)%Vsource(1))
    allocate(media%wire(1)%Isource(1))

    call DestroyWireMedia(media)

    if (associated(media%wire)) then
        print *, "test_destroyWireMedia_thinwire FAILED: wire storage still associated"
        status = 1
    end if
end function test_destroyWireMedia_thinwire

integer function test_evolucion_out_of_range_low() bind(C, name="test_evolucion_out_of_range_low") result(status)
    use HollandWires_m, only: evolucion
    use FDETYPES_m, only: RKIND_wires
    implicit none
    integer(kind=4) :: numus
    real(kind=RKIND_wires), dimension(0:3) :: evol
    real(kind=RKIND_wires) :: deltaevol, t, result

    status = 0
    numus = 3
    deltaevol = 1.0_RKIND_wires
    evol = [0.0_RKIND_wires, 1.0_RKIND_wires, 2.0_RKIND_wires, 3.0_RKIND_wires]
    ! t = -1.5 -> nprev = int(-1.5) = -1 -> nprev+1 = 0 <= 0 -> returns 0
    t = -1.5_RKIND_wires
    result = evolucion(t, evol, deltaevol, numus)
    if (abs(result) > 1.0e-12_RKIND_wires) then
        print *, "test_evolucion_out_of_range_low FAILED: result=", result
        status = 1
    end if
end function test_evolucion_out_of_range_low

integer function test_evolucion_out_of_range_high() bind(C, name="test_evolucion_out_of_range_high") result(status)
    use HollandWires_m, only: evolucion
    use FDETYPES_m, only: RKIND_wires
    implicit none
    integer(kind=4) :: numus
    real(kind=RKIND_wires), dimension(0:3) :: evol
    real(kind=RKIND_wires) :: deltaevol, t, result

    status = 0
    numus = 3
    deltaevol = 1.0_RKIND_wires
    evol = [0.0_RKIND_wires, 1.0_RKIND_wires, 2.0_RKIND_wires, 3.0_RKIND_wires]
    ! t = 3.0 -> nprev = int(3.0) = 3 = numus -> nprev+1 = 4 > 3 = numus -> returns 0
    t = 3.0_RKIND_wires
    result = evolucion(t, evol, deltaevol, numus)
    if (abs(result) > 1.0e-12_RKIND_wires) then
        print *, "test_evolucion_out_of_range_high FAILED: result=", result
        status = 1
    end if
end function test_evolucion_out_of_range_high

integer function test_evolucion_at_zero() bind(C, name="test_evolucion_at_zero") result(status)
    use HollandWires_m, only: evolucion
    use FDETYPES_m, only: RKIND_wires
    implicit none
    integer(kind=4) :: numus
    real(kind=RKIND_wires), dimension(0:3) :: evol
    real(kind=RKIND_wires) :: deltaevol, t, result

    status = 0
    numus = 3
    deltaevol = 1.0_RKIND_wires
    evol = [5.0_RKIND_wires, 1.0_RKIND_wires, 2.0_RKIND_wires, 3.0_RKIND_wires]
    ! t = 0.0 -> nprev = 0, returns evol(0) = 5.0
    t = 0.0_RKIND_wires
    result = evolucion(t, evol, deltaevol, numus)
    if (abs(result - 5.0_RKIND_wires) > 1.0e-12_RKIND_wires) then
        print *, "test_evolucion_at_zero FAILED: result=", result
        status = 1
    end if
end function test_evolucion_at_zero

integer function test_evolucion_midpoint_interp() bind(C, name="test_evolucion_midpoint_interp") result(status)
    use HollandWires_m, only: evolucion
    use FDETYPES_m, only: RKIND_wires
    implicit none
    integer(kind=4) :: numus
    real(kind=RKIND_wires), dimension(0:3) :: evol
    real(kind=RKIND_wires) :: deltaevol, t, result, expected

    status = 0
    numus = 3
    deltaevol = 1.0_RKIND_wires
    evol = [0.0_RKIND_wires, 2.0_RKIND_wires, 4.0_RKIND_wires, 6.0_RKIND_wires]
    ! t = 0.5 -> nprev = 0, result = (2-0)/1*0.5 + 0 = 1.0
    t = 0.5_RKIND_wires
    expected = 1.0_RKIND_wires
    result = evolucion(t, evol, deltaevol, numus)
    if (abs(result - expected) > 1.0e-10_RKIND_wires) then
        print *, "test_evolucion_midpoint_interp FAILED: result=", result, " expected=", expected
        status = 1
    end if
end function test_evolucion_midpoint_interp

integer function test_evolucion_exact_sample() bind(C, name="test_evolucion_exact_sample") result(status)
    use HollandWires_m, only: evolucion
    use FDETYPES_m, only: RKIND_wires
    implicit none
    integer(kind=4) :: numus
    real(kind=RKIND_wires), dimension(0:3) :: evol
    real(kind=RKIND_wires) :: deltaevol, t, result

    status = 0
    numus = 3
    deltaevol = 1.0_RKIND_wires
    evol = [0.0_RKIND_wires, 7.0_RKIND_wires, 4.0_RKIND_wires, 6.0_RKIND_wires]
    ! t = 1.0 -> nprev = 1, result = (evol(2)-evol(1))/1*0 + evol(1) = 7.0
    t = 1.0_RKIND_wires
    result = evolucion(t, evol, deltaevol, numus)
    if (abs(result - 7.0_RKIND_wires) > 1.0e-12_RKIND_wires) then
        print *, "test_evolucion_exact_sample FAILED: result=", result
        status = 1
    end if
end function test_evolucion_exact_sample

integer function test_destroyWireMedia_non_thinwire() bind(C, name="test_destroyWireMedia_non_thinwire") result(status)
    use HollandWires_m, only: DestroyWireMedia
    use FDETYPES_m, only: MediaData_t
    implicit none
    type(MediaData_t) :: media

    status = 0

    media%Is%ThinWire = .false.
    allocate(media%wire(1))
    allocate(media%wire(1)%Vsource(1))
    allocate(media%wire(1)%Isource(1))

    call DestroyWireMedia(media)

    if (.not. associated(media%wire)) then
        print *, "test_destroyWireMedia_non_thinwire FAILED: non-thinwire storage was deallocated"
        status = 1
    end if

    if (associated(media%wire(1)%Vsource)) deallocate(media%wire(1)%Vsource)
    if (associated(media%wire(1)%Isource)) deallocate(media%wire(1)%Isource)
    if (associated(media%wire)) deallocate(media%wire)
end function test_destroyWireMedia_non_thinwire
