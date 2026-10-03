! Unit tests for src_main_pub/maloney_thin_sheet.F90
! Checks the in-sheet normal electric field update (eq. 8) and the tangential
! magnetic field correction (eqs. 15-16) of the Maloney-Smith model.

! E update: Eintern = g1*Eintern + g2*((H1-H2)*deltaA - (H3-H4)*deltaB)
integer function test_maloneysheet_e_update() bind(C) result(status)
    use maloney_thin_sheet_m, only: maloneySheets_t, GetMaloneySheets, AdvanceMaloneySheetE
    use FDETYPES_m, only: RKIND
    implicit none
    type(maloneySheets_t), pointer :: mal
    real(kind=RKIND), target :: h1, h2, h3, h4, eo
    real(kind=RKIND) :: curl, expected
    real(kind=RKIND), parameter :: tol = 1.0e-6_RKIND

    status = 0
    mal => GetMaloneySheets()
    if (allocated(mal%nodes)) deallocate(mal%nodes)
    allocate(mal%nodes(1))
    mal%numNodes = 1

    h1 = 1.0_RKIND;  h2 = 0.5_RKIND;  h3 = -2.0_RKIND;  h4 = 0.25_RKIND
    eo = 0.0_RKIND

    mal%nodes(1)%Efield => eo
    mal%nodes(1)%Eintern = 1.0_RKIND
    mal%nodes(1)%deltaA = 2.0_RKIND
    mal%nodes(1)%deltaB = 4.0_RKIND
    mal%nodes(1)%g1 = 0.5_RKIND
    mal%nodes(1)%g2 = 0.1_RKIND
    mal%nodes(1)%H1 => h1
    mal%nodes(1)%H2 => h2
    mal%nodes(1)%H3 => h3
    mal%nodes(1)%H4 => h4

    curl = (h1 - h2) * mal%nodes(1)%deltaA - (h3 - h4) * mal%nodes(1)%deltaB
    expected = 0.5_RKIND * 1.0_RKIND + 0.1_RKIND * curl

    call AdvanceMaloneySheetE()

    if (abs(mal%nodes(1)%Eintern - expected) > tol) then
        print *, "test_maloneysheet_e_update FAILED: ", mal%nodes(1)%Eintern, " expected ", expected
        status = 1
    end if

    deallocate(mal%nodes)
    mal%numNodes = 0
end function test_maloneysheet_e_update

! H correction: each H component is shifted by c*beta*(Eintern - Efield)
integer function test_maloneysheet_h_correction() bind(C) result(status)
    use maloney_thin_sheet_m, only: maloneySheets_t, GetMaloneySheets, AdvanceMaloneySheetH
    use FDETYPES_m, only: RKIND
    implicit none
    type(maloneySheets_t), pointer :: mal
    real(kind=RKIND), target :: h1, h2, h3, h4, eo
    real(kind=RKIND) :: dE
    real(kind=RKIND), parameter :: tol = 1.0e-6_RKIND

    status = 0
    mal => GetMaloneySheets()
    if (allocated(mal%nodes)) deallocate(mal%nodes)
    allocate(mal%nodes(1))
    mal%numNodes = 1

    h1 = 1.0_RKIND;  h2 = 2.0_RKIND;  h3 = 3.0_RKIND;  h4 = 4.0_RKIND
    eo = 1.0_RKIND

    mal%nodes(1)%Efield => eo
    mal%nodes(1)%Eintern = 2.0_RKIND
    mal%nodes(1)%beta = 0.25_RKIND
    mal%nodes(1)%c1 = -1.0_RKIND
    mal%nodes(1)%c2 = 0.5_RKIND
    mal%nodes(1)%c3 = 2.0_RKIND
    mal%nodes(1)%c4 = -3.0_RKIND
    mal%nodes(1)%H1 => h1
    mal%nodes(1)%H2 => h2
    mal%nodes(1)%H3 => h3
    mal%nodes(1)%H4 => h4

    dE = 0.25_RKIND * (2.0_RKIND - 1.0_RKIND)

    call AdvanceMaloneySheetH()

    if (abs(h1 - (1.0_RKIND - dE)) > tol) then
        print *, "test_maloneysheet_h_correction FAILED: h1=", h1
        status = 1
    end if
    if (abs(h2 - (2.0_RKIND + 0.5_RKIND*dE)) > tol) then
        print *, "test_maloneysheet_h_correction FAILED: h2=", h2
        status = 1
    end if
    if (abs(h3 - (3.0_RKIND + 2.0_RKIND*dE)) > tol) then
        print *, "test_maloneysheet_h_correction FAILED: h3=", h3
        status = 1
    end if
    if (abs(h4 - (4.0_RKIND - 3.0_RKIND*dE)) > tol) then
        print *, "test_maloneysheet_h_correction FAILED: h4=", h4
        status = 1
    end if

    deallocate(mal%nodes)
    mal%numNodes = 0
end function test_maloneysheet_h_correction
