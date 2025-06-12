module test_rotate_generateSONDAs_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_sondas() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    mpidir = 2
    call setup_sonda_test(this)
    call rotate_generatesondas(this, mpidir)
    call verify_sonda_rotation(test_err, this%oldSONDA, mpidir)
    call cleanup_sonda_test(this)

    mpidir = 1
    call setup_sonda_test(this)
    call rotate_generatesondas(this, mpidir)

    call verify_sonda_rotation(test_err, this%oldSONDA, mpidir)
    call cleanup_sonda_test(this)
    err = test_err
end function test_rotate_generate_sondas

subroutine setup_sonda_test(this)
    type(Parseador), intent(inout) :: this
    
    allocate(this%oldSONDA)
    allocate(this%oldSONDA%probes(1))
    allocate(this%oldSONDA%probes(1)%Electric(1))
    allocate(this%oldSONDA%probes(1)%Magnetic(1))
    allocate(this%oldSONDA%probes(1)%NormalElectric(1))
    allocate(this%oldSONDA%probes(1)%NormalMagnetic(1))
    allocate(this%oldSONDA%probes(1)%SurfaceElectricCurrent(1))
    allocate(this%oldSONDA%probes(1)%SurfaceMagneticCurrent(1))
    allocate(this%oldSONDA%probes(1)%FarField(1))

    this%oldSONDA%n_probes=1
    this%oldSONDA%probes(1)%n_Electric = 1
    this%oldSONDA%probes(1)%n_Magnetic = 1
    this%oldSONDA%probes(1)%n_FarField = 1

    call setup_basic_sonda(this%oldSONDA%probes(1)%Electric(1)%probe)
    call setup_basic_sonda(this%oldSONDA%probes(1)%Magnetic(1)%probe)

    
    this%oldSONDA%probes(1)%FarField(1)%probe%thetastart = 1.5
    this%oldSONDA%probes(1)%FarField(1)%probe%phistart = 2
    this%oldSONDA%probes(1)%FarField(1)%probe%thetastop = 1.5
    this%oldSONDA%probes(1)%FarField(1)%probe%phistop = 2

end subroutine setup_sonda_test

subroutine setup_basic_sonda(probe)
    type(Sonda), intent(inout) :: probe
    allocate(probe%i(1), source = 1)
    allocate(probe%j(1), source = 2)
    allocate(probe%k(1), source = 3)
    probe%n_cord = 1

end subroutine setup_basic_sonda

subroutine assert_basic_sonda(test_err, mpidir, probe)
    type(Sonda), intent(inout) :: probe
    integer :: mpidir, test_err
    if (mpidir==2) then
        call expect_eq_int(test_err, 3, probe%i(1), "Assertion failed for i mpidir 2")
        call expect_eq_int(test_err, 1, probe%j(1), "Assertion failed for j mpidir 2")
        call expect_eq_int(test_err, 2, probe%k(1), "Assertion failed for k mpidir 2")
    end if
    if (mpidir==1) then
        call expect_eq_int(test_err, 2, probe%i(1), "Assertion failed for i mpidir 1")
        call expect_eq_int(test_err, 3, probe%j(1), "Assertion failed for j mpidir 1")
        call expect_eq_int(test_err, 1, probe%k(1), "Assertion failed for k mpidir 1")
    end if
end subroutine assert_basic_sonda

subroutine verify_sonda_rotation(test_err, oldSonda, mpidir)
    integer, intent(inout) :: test_err, mpidir
    type(Sondas), intent(in) :: oldSonda

    call assert_basic_sonda(test_err, mpidir, oldSONDA%probes(1)%Electric(1)%probe)
    call assert_basic_sonda(test_err, mpidir, oldSONDA%probes(1)%Magnetic(1)%probe)

    if (mpidir==2) then
        call expect_eq_real(test_err, -1.40200949 ,oldSONDA%probes(1)%FarField(1)%probe%thetastart, "Assertion error for mpidir 2")
        call expect_eq_real(test_err, 2.00000000 ,oldSONDA%probes(1)%FarField(1)%probe%phistart, "Assertion error for mpidir 2")
        call expect_eq_real(test_err, 0.434644908 ,oldSONDA%probes(1)%FarField(1)%probe%thetastop, "Assertion error for mpidir 2")
        call expect_eq_real(test_err, -1.40200949 ,oldSONDA%probes(1)%FarField(1)%probe%phistop, "Assertion error for mpidir 2")
    end if

    if (mpidir==1) then
        call expect_eq_real(test_err, 7.78310671E-02 ,oldSONDA%probes(1)%FarField(1)%probe%thetastart, "Assertion error for mpidir 1")
        call expect_eq_real(test_err, 2.00000000 ,oldSONDA%probes(1)%FarField(1)%probe%phistart, "Assertion error for mpidir 1")
        call expect_eq_real(test_err, 1.99885380 ,oldSONDA%probes(1)%FarField(1)%probe%thetastop, "Assertion error for mpidir 1")
        call expect_eq_real(test_err, 7.78310671E-02 ,oldSONDA%probes(1)%FarField(1)%probe%phistop, "Assertion error for mpidir 1")
    end if
end subroutine verify_sonda_rotation

subroutine cleanup_sonda_test(this)
    type(Parseador), intent(inout) :: this
    deallocate(this%oldSONDA)

end subroutine cleanup_sonda_test

end module test_rotate_generateSONDAs_m 