integer function test_rotate_generate_anisotropics() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0

    ! Test case 1: Verify that rotations are not supported for anisotropic materials
    mpidir = 2
    allocate(this%ANIMATS)
    this%ANIMATS%nvols = 1
    this%ANIMATS%nsurfs = 1
    this%ANIMATS%nlins = 1
    call rotate_generateANISOTROPICs(this, mpidir)
    call expect_eq_int(test_err, 1, this%ANIMATS%nvols, "rotate_generateANISOTROPICs: nvols should remain unchanged")
    call expect_eq_int(test_err, 1, this%ANIMATS%nsurfs, "rotate_generateANISOTROPICs: nsurfs should remain unchanged")
    call expect_eq_int(test_err, 1, this%ANIMATS%nlins, "rotate_generateANISOTROPICs: nlins should remain unchanged")
    deallocate(this%ANIMATS)

    ! Test case 2: Verify that rotations are not supported for mpidir=1
    mpidir = 1
    allocate(this%ANIMATS)
    this%ANIMATS%nvols = 1
    this%ANIMATS%nsurfs = 1
    this%ANIMATS%nlins = 1
    call rotate_generateANISOTROPICs(this, mpidir)
    call expect_eq_int(test_err, 1, this%ANIMATS%nvols, "rotate_generateANISOTROPICs: nvols should remain unchanged")
    call expect_eq_int(test_err, 1, this%ANIMATS%nsurfs, "rotate_generateANISOTROPICs: nsurfs should remain unchanged")
    call expect_eq_int(test_err, 1, this%ANIMATS%nlins, "rotate_generateANISOTROPICs: nlins should remain unchanged")
    deallocate(this%ANIMATS)

    ! Test case 3: Verify that no warning is printed when there are no anisotropic materials
    mpidir = 2
    allocate(this%ANIMATS)
    this%ANIMATS%nvols = 0
    this%ANIMATS%nsurfs = 0
    this%ANIMATS%nlins = 0
    call rotate_generateANISOTROPICs(this, mpidir)
    call expect_eq_int(test_err, 0, this%ANIMATS%nvols, "rotate_generateANISOTROPICs: nvols should remain 0")
    call expect_eq_int(test_err, 0, this%ANIMATS%nsurfs, "rotate_generateANISOTROPICs: nsurfs should remain 0")
    call expect_eq_int(test_err, 0, this%ANIMATS%nlins, "rotate_generateANISOTROPICs: nlins should remain 0")
    deallocate(this%ANIMATS)

    err = test_err
end function test_rotate_generate_anisotropics