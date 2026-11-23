! Unit tests for startup memory estimation and budget formatting.

integer function test_memory_amount_format() bind(C, name="test_memory_amount_format") result(status)
    use memoryBudget_m
    implicit none

    status = 0
    call expect_text(format_memory_amount(1073741824_8), '1.0 GB')
    call expect_text(format_memory_amount(1048576_8), '1.0 MB')
    call expect_text(format_memory_amount(0_8), '0.0 KB')

contains

    subroutine expect_text(actual, expected)
        character(len=*), intent(in) :: actual, expected
        if (trim(actual) /= expected) then
            print *, 'test_memory_amount_format FAILED: got "', trim(actual), '" expected "', expected, '"'
            status = 1
        end if
    end subroutine expect_text

end function test_memory_amount_format

integer function test_solver_memory_estimate() bind(C, name="test_solver_memory_estimate") result(status)
    use FDETYPES_m
    use memoryBudget_m
    implicit none

    type(SGGFDTDINFO_t) :: sgg
    type(media_matrices_t) :: media
    type(taglist_t) :: tags
    type(memory_estimate_t) :: estimate
    integer(kind=4) :: field

    status = 0
    do field = IEX, IHZ
        sgg%Alloc(field)%XI = 1
        sgg%Alloc(field)%XE = 2
        sgg%Alloc(field)%YI = 1
        sgg%Alloc(field)%YE = 2
        sgg%Alloc(field)%ZI = 1
        sgg%Alloc(field)%ZE = 2
    end do

    allocate(media%sggMiNo(1:2, 1:2, 1:2), media%sggMiEx(1:2, 1:2, 1:2), &
             media%sggMiEy(1:2, 1:2, 1:2), media%sggMiEz(1:2, 1:2, 1:2), &
             media%sggMiHx(1:2, 1:2, 1:2), media%sggMiHy(1:2, 1:2, 1:2), &
             media%sggMiHz(1:2, 1:2, 1:2), media%sggMtag(1:2, 1:2, 1:2))
    media%sggMiNo = 1; media%sggMiEx = 1; media%sggMiEy = 1; media%sggMiEz = 1
    media%sggMiHx = 1; media%sggMiHy = 1; media%sggMiHz = 1; media%sggMtag = 0

    allocate(tags%edge%x(1:2, 1:2, 1:2), tags%edge%y(1:2, 1:2, 1:2), tags%edge%z(1:2, 1:2, 1:2))
    allocate(tags%face%x(1:2, 1:2, 1:2), tags%face%y(1:2, 1:2, 1:2), tags%face%z(1:2, 1:2, 1:2))
    tags%edge%x = 0; tags%edge%y = 0; tags%edge%z = 0
    tags%face%x = 0; tags%face%y = 0; tags%face%z = 0

    call estimate_solver_memory(sgg, media, tags, .true., estimate)
    if (estimate%fieldBytes <= 0_8) status = 1
    if (estimate%distanceBytes <= 0_8) status = 1
    if (estimate%mediaCopyBytes <= 0_8) status = 1
    if (estimate%tagCopyBytes <= 0_8) status = 1
    if (estimate%totalBytes /= estimate%fieldBytes + estimate%distanceBytes + &
                                 estimate%mediaCopyBytes + estimate%tagCopyBytes) status = 1

    deallocate(media%sggMiNo, media%sggMiEx, media%sggMiEy, media%sggMiEz, &
               media%sggMiHx, media%sggMiHy, media%sggMiHz, media%sggMtag)
    deallocate(tags%edge%x, tags%edge%y, tags%edge%z, tags%face%x, tags%face%y, tags%face%z)
    if (status /= 0) print *, 'test_solver_memory_estimate FAILED'
end function test_solver_memory_estimate

integer function test_system_memory_probe() bind(C, name="test_system_memory_probe") result(status)
    use systemMemory_m
    implicit none

    type(system_memory_t) :: status_memory

    status = 0
    call read_system_memory(status_memory)
    if (status_memory%totalKnown) then
        if (status_memory%totalBytes <= 0_8) status = 1
    end if
    if (status_memory%availableKnown) then
        if (status_memory%availableBytes < 0_8) status = 1
        if (status_memory%totalKnown .and. status_memory%availableBytes > status_memory%totalBytes) status = 1
    end if
    if (status_memory%residentKnown .and. status_memory%residentBytes <= 0_8) status = 1
    if (status /= 0) print *, 'test_system_memory_probe FAILED'
end function test_system_memory_probe

integer function test_memory_failure_message() bind(C, name="test_memory_failure_message") result(status)
    use memoryBudget_m
    implicit none

    type(memory_verdict_t) :: verdict
    character(len=1024) :: message

    status = 0
    verdict%commSize = 1
    verdict%requiredBytes = 51754345988_8
    verdict%availableBytes = 33715493274_8
    message = format_memory_verdict_message(verdict, 'before time stepping')
    if (index(message, 'Not enough memory to run the simulation.') == 0) status = 1
    if (index(message, 'Required:') == 0) status = 1
    if (index(message, 'Available:') == 0) status = 1

    verdict%commSize = 4
    verdict%offenderRank = 3
    verdict%offenderNode = 'node03'
    message = format_memory_verdict_message(verdict, 'before SGBC initialization')
    if (index(message, 'rank 3 (node03)') == 0) status = 1
    if (index(message, 'Required by process:') == 0) status = 1
    if (index(message, 'Available to process:') == 0) status = 1
    if (status /= 0) print *, 'test_memory_failure_message FAILED'
end function test_memory_failure_message

integer function test_cpml_estimate() bind(C, name="test_cpml_estimate") result(status)
    use FDETYPES_m
    use BORDERS_CPML_m, only: estimate_cpml_memory
    implicit none

    type(SGGFDTDINFO_t) :: sgg
    type(limit_t) :: sinpml(6)
    integer(kind=8) :: noPmlBytes, withPmlBytes
    integer(kind=4) :: field

    status = 0
    do field = IEX, IHZ
        sgg%Alloc(field)%XI = 1; sgg%Alloc(field)%XE = 4
        sgg%Alloc(field)%YI = 1; sgg%Alloc(field)%YE = 4
        sgg%Alloc(field)%ZI = 1; sgg%Alloc(field)%ZE = 4
        sgg%Sweep(field) = sgg%Alloc(field)
        sinpml(field)%XI = 2; sinpml(field)%XE = 3
        sinpml(field)%YI = 2; sinpml(field)%YE = 3
        sinpml(field)%ZI = 1; sinpml(field)%ZE = 4
    end do
    sgg%Border%IsBackPML = .false.; sgg%Border%IsFrontPML = .false.
    sgg%Border%IsLeftPML = .false.; sgg%Border%IsRightPML = .false.
    sgg%Border%IsUpPML = .false.; sgg%Border%IsDownPML = .false.

    call estimate_cpml_memory(sgg, sinpml, noPmlBytes)
    sgg%Border%IsBackPML = .true.; sgg%Border%IsFrontPML = .true.
    sgg%Border%IsLeftPML = .true.; sgg%Border%IsRightPML = .true.
    sgg%Border%IsUpPML = .true.; sgg%Border%IsDownPML = .true.
    call estimate_cpml_memory(sgg, sinpml, withPmlBytes)
    if (noPmlBytes <= 0_8) status = 1
    if (withPmlBytes <= noPmlBytes) status = 1
    if (status /= 0) print *, 'test_cpml_estimate FAILED'
end function test_cpml_estimate

#ifndef CompileWithStochastic
integer function test_sgbc_estimate() bind(C, name="test_sgbc_estimate") result(status)
    use FDETYPES_m
    use SGBC_nostoch_m, only: estimate_sgbc_memory
    implicit none

    type(SGGFDTDINFO_t) :: sgg
    type(media_matrices_t) :: media
    integer(kind=8) :: bytes
    integer(kind=4) :: field

    status = 0
    sgg%NumMedia = 1
    allocate(sgg%Med(0:1))
    sgg%Med(1)%Is%SGBC = .true.
    sgg%Med(1)%Is%SGBCDispersive = .false.
    allocate(sgg%Med(1)%Multiport(1))
    sgg%Med(1)%Multiport(1)%numLayers = 1
    allocate(sgg%Med(1)%Multiport(1)%width(1), sgg%Med(1)%Multiport(1)%epr(1), &
             sgg%Med(1)%Multiport(1)%mur(1), sgg%Med(1)%Multiport(1)%sigma(1), &
             sgg%Med(1)%Multiport(1)%sigmam(1))
    sgg%Med(1)%Multiport(1)%width = 1.0_RKIND
    sgg%Med(1)%Multiport(1)%epr = 2.0_RKIND
    sgg%Med(1)%Multiport(1)%mur = 1.0_RKIND
    sgg%Med(1)%Multiport(1)%sigma = 0.0_RKIND
    sgg%Med(1)%Multiport(1)%sigmam = 0.0_RKIND
    do field = IEX, IHZ
        sgg%SINPMLSweep(field)%XI = 1; sgg%SINPMLSweep(field)%XE = 1
        sgg%SINPMLSweep(field)%YI = 1; sgg%SINPMLSweep(field)%YE = 1
        sgg%SINPMLSweep(field)%ZI = 1; sgg%SINPMLSweep(field)%ZE = 1
    end do
    allocate(media%sggMiEx(1:1, 1:1, 1:1), media%sggMiEy(1:1, 1:1, 1:1), media%sggMiEz(1:1, 1:1, 1:1))
    media%sggMiEx = 1; media%sggMiEy = 1; media%sggMiEz = 1

    call estimate_sgbc_memory(sgg, media, .true., .true., 1.0e9_RKIND, 1.0_RKIND, 2_4, .false., &
                              EPSILON_VACUUM, MU_VACUUM, bytes)
    if (bytes <= 0_8) status = 1

    deallocate(media%sggMiEx, media%sggMiEy, media%sggMiEz)
    deallocate(sgg%Med(1)%Multiport(1)%width, sgg%Med(1)%Multiport(1)%epr, &
               sgg%Med(1)%Multiport(1)%mur, sgg%Med(1)%Multiport(1)%sigma, &
               sgg%Med(1)%Multiport(1)%sigmam)
    deallocate(sgg%Med(1)%Multiport)
    deallocate(sgg%Med)
    if (status /= 0) print *, 'test_sgbc_estimate FAILED'
end function test_sgbc_estimate
#endif

integer function test_holland_wire_estimate() bind(C, name="test_holland_wire_estimate") result(status)
    use FDETYPES_m
    use HollandWires_m, only: estimate_holland_wires_memory
    implicit none

    type(SGGFDTDINFO_t) :: sgg
    integer(kind=8) :: bytes
    integer(kind=4) :: field

    status = 0
    sgg%NumMedia = 1
    allocate(sgg%Med(0:1))
    sgg%Med(1)%Is%ThinWire = .true.
    allocate(sgg%Med(1)%Wire(1))
    allocate(sgg%Med(1)%Wire(1)%segm(1))
    sgg%Med(1)%Wire(1)%numsegmentos = 1
    sgg%Med(1)%Wire(1)%segm(1)%i = 1
    sgg%Med(1)%Wire(1)%segm(1)%j = 1
    sgg%Med(1)%Wire(1)%segm(1)%k = 1
    sgg%Med(1)%Wire(1)%segm(1)%ori = IEZ
    sgg%Med(1)%Wire(1)%segm(1)%multirabo = .false.
    do field = IEX, IHZ
        sgg%Alloc(field)%XI = 1; sgg%Alloc(field)%XE = 2
        sgg%Alloc(field)%YI = 1; sgg%Alloc(field)%YE = 2
        sgg%Alloc(field)%ZI = 1; sgg%Alloc(field)%ZE = 2
    end do

    call estimate_holland_wires_memory(sgg, bytes)
    if (bytes <= 0_8) status = 1

    deallocate(sgg%Med(1)%Wire(1)%segm)
    deallocate(sgg%Med(1)%Wire)
    deallocate(sgg%Med)
    if (status /= 0) print *, 'test_holland_wire_estimate FAILED'
end function test_holland_wire_estimate

integer function test_conformal_estimate() bind(C, name="test_conformal_estimate") result(status)
    use FDETYPES_m
    use conformal_timestepping_m, only: estimate_conformal_memory
    implicit none

    type(SGGFDTDINFO_t) :: sgg
    integer(kind=8) :: bytes

    status = 0
    sgg%NumMedia = 1
    allocate(sgg%Med(0:1))
    sgg%Med(1)%Is%ConformalPEC = .true.
    sgg%Med(1)%Is%Surface = .true.
    sgg%Med(1)%Is%PML = .false.
    allocate(sgg%Med(1)%ConformalFace(1))
    allocate(sgg%Med(1)%ConformalFace(1)%contour_inside_edges(0))
    allocate(sgg%Med(1)%ConformalEdge(1))

    call estimate_conformal_memory(sgg, bytes)
    if (bytes <= 0_8) status = 1

    deallocate(sgg%Med(1)%ConformalFace(1)%contour_inside_edges)
    deallocate(sgg%Med(1)%ConformalFace)
    deallocate(sgg%Med(1)%ConformalEdge)
    deallocate(sgg%Med)
    if (status /= 0) print *, 'test_conformal_estimate FAILED'
end function test_conformal_estimate
