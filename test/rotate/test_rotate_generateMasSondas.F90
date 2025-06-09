module test_rotate_generateMasSondas_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_mas_sondas() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    ! Test case 1: Rotate MasSondas in mpidir=2 direction
    mpidir = 2
    call setup_massonda_test(this, 2)
    
    ! Set up test data
    call init_massonda_data(this, 1, 1, 4, 2, 5, 3, 6, 1, 1, 1, 1, 1, "probe1")
    call init_massonda_data(this, 2, 5, 8, 6, 9, 7, 10, 1, 1, 1, 2, 2, "probe2")
    
    ! Call rotation and verify results
    call rotate_generateMasSondas(this, mpidir)
    call verify_massonda_rotation_y(test_err, this)
    
    ! Clean up
    call cleanup_massonda_test(this)
    
    ! Test case 2: Rotate MasSondas in mpidir=1 direction
    mpidir = 1
    call setup_massonda_test(this, 2)
    
    ! Set up test data (same as before)
    call init_massonda_data(this, 1, 1, 4, 2, 5, 3, 6, 1, 1, 1, 1, 1, "probe1")
    call init_massonda_data(this, 2, 5, 8, 6, 9, 7, 10, 1, 1, 1, 2, 2, "probe2")
    
    ! Call rotation and verify results
    call rotate_generateMasSondas(this, mpidir)
    call verify_massonda_rotation_x(test_err, this)
    
    ! Clean up
    call cleanup_massonda_test(this)
    
    ! Test case 3: Verify behavior with no MasSondas
    mpidir = 2
    call setup_massonda_test(this, 0)
    call rotate_generateMasSondas(this, mpidir)
    call expect_eq_int(test_err, 0, this%Sonda%length, "rotate_generateMasSondas: length should remain 0")
    call cleanup_massonda_test(this)
    
    err = test_err
end function test_rotate_generate_mas_sondas

subroutine setup_massonda_test(this, n_probes)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_probes
    integer :: i
    
    allocate(this%Sonda)
    this%Sonda%length = n_probes
    this%Sonda%length_max = n_probes
    this%Sonda%len_cor_max = 1  ! Each probe has one coordinate set
    if (n_probes > 0) then
        allocate(this%Sonda%collection(n_probes))
        ! Initialize default values for each probe
        do i = 1, n_probes
            this%Sonda%collection(i)%tstart = 0.0d0
            this%Sonda%collection(i)%tstop = 1.0d0
            this%Sonda%collection(i)%tstep = 0.1d0
            this%Sonda%collection(i)%fstart = 0.0d0
            this%Sonda%collection(i)%fstop = 1.0d0
            this%Sonda%collection(i)%fstep = 0.1d0
            this%Sonda%collection(i)%filename = "probe.dat"
            this%Sonda%collection(i)%outputrequest = "output.txt"
            this%Sonda%collection(i)%len_cor = 1
            this%Sonda%collection(i)%type1 = 1
            this%Sonda%collection(i)%type2 = 1
            ! Allocate coordinates array for each probe
            allocate(this%Sonda%collection(i)%cordinates(1))
        end do
    end if
end subroutine setup_massonda_test

subroutine init_massonda_data(this, idx, xi, xe, yi, ye, zi, ze, xtrancos, ytrancos, ztrancos, type1, type2, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, xi, xe, yi, ye, zi, ze, xtrancos, ytrancos, ztrancos, type1, type2
    character(len=*), intent(in) :: tag
    
    ! Initialize coordinates for the probe
    this%Sonda%collection(idx)%cordinates(1)%Xi = xi
    this%Sonda%collection(idx)%cordinates(1)%Xe = xe
    this%Sonda%collection(idx)%cordinates(1)%Yi = yi
    this%Sonda%collection(idx)%cordinates(1)%Ye = ye
    this%Sonda%collection(idx)%cordinates(1)%Zi = zi
    this%Sonda%collection(idx)%cordinates(1)%Ze = ze
    this%Sonda%collection(idx)%cordinates(1)%Xtrancos = xtrancos
    this%Sonda%collection(idx)%cordinates(1)%Ytrancos = ytrancos
    this%Sonda%collection(idx)%cordinates(1)%Ztrancos = ztrancos
    this%Sonda%collection(idx)%cordinates(1)%Or = type1  ! Using type1 as orientation
    this%Sonda%collection(idx)%cordinates(1)%tag = tag
    this%Sonda%collection(idx)%type1 = type1
    this%Sonda%collection(idx)%type2 = type2
end subroutine init_massonda_data

subroutine verify_massonda_rotation_y(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify MasSonda 1 (rotated around y-axis)
    call expect_eq_int(test_err, -3, this%Sonda%collection(1)%cordinates(1)%Xi, "rotate_generateMasSondas: Xi(1) should be rotated")
    call expect_eq_int(test_err, -6, this%Sonda%collection(1)%cordinates(1)%Xe, "rotate_generateMasSondas: Xe(1) should be rotated")
    call expect_eq_int(test_err, 2, this%Sonda%collection(1)%cordinates(1)%Yi, "rotate_generateMasSondas: Yi(1) should remain unchanged")
    call expect_eq_int(test_err, 5, this%Sonda%collection(1)%cordinates(1)%Ye, "rotate_generateMasSondas: Ye(1) should remain unchanged")
    call expect_eq_int(test_err, -1, this%Sonda%collection(1)%cordinates(1)%Zi, "rotate_generateMasSondas: Zi(1) should be rotated")
    call expect_eq_int(test_err, -4, this%Sonda%collection(1)%cordinates(1)%Ze, "rotate_generateMasSondas: Ze(1) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Xtrancos, "rotate_generateMasSondas: Xtrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Ytrancos, "rotate_generateMasSondas: Ytrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Ztrancos, "rotate_generateMasSondas: Ztrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Or, "rotate_generateMasSondas: Or(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%type1, "rotate_generateMasSondas: type1(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%type2, "rotate_generateMasSondas: type2(1) should remain unchanged")
    
    ! Verify MasSonda 2 (rotated around y-axis)
    call expect_eq_int(test_err, -7, this%Sonda%collection(2)%cordinates(1)%Xi, "rotate_generateMasSondas: Xi(2) should be rotated")
    call expect_eq_int(test_err, -10, this%Sonda%collection(2)%cordinates(1)%Xe, "rotate_generateMasSondas: Xe(2) should be rotated")
    call expect_eq_int(test_err, 6, this%Sonda%collection(2)%cordinates(1)%Yi, "rotate_generateMasSondas: Yi(2) should remain unchanged")
    call expect_eq_int(test_err, 9, this%Sonda%collection(2)%cordinates(1)%Ye, "rotate_generateMasSondas: Ye(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%Sonda%collection(2)%cordinates(1)%Zi, "rotate_generateMasSondas: Zi(2) should be rotated")
    call expect_eq_int(test_err, -8, this%Sonda%collection(2)%cordinates(1)%Ze, "rotate_generateMasSondas: Ze(2) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Xtrancos, "rotate_generateMasSondas: Xtrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Ytrancos, "rotate_generateMasSondas: Ytrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Ztrancos, "rotate_generateMasSondas: Ztrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%cordinates(1)%Or, "rotate_generateMasSondas: Or(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%type1, "rotate_generateMasSondas: type1(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%type2, "rotate_generateMasSondas: type2(2) should remain unchanged")
end subroutine verify_massonda_rotation_y

subroutine verify_massonda_rotation_x(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify MasSonda 1 (rotated around x-axis)
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Xi, "rotate_generateMasSondas: Xi(1) should remain unchanged")
    call expect_eq_int(test_err, 4, this%Sonda%collection(1)%cordinates(1)%Xe, "rotate_generateMasSondas: Xe(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%Sonda%collection(1)%cordinates(1)%Yi, "rotate_generateMasSondas: Yi(1) should be rotated")
    call expect_eq_int(test_err, -6, this%Sonda%collection(1)%cordinates(1)%Ye, "rotate_generateMasSondas: Ye(1) should be rotated")
    call expect_eq_int(test_err, 2, this%Sonda%collection(1)%cordinates(1)%Zi, "rotate_generateMasSondas: Zi(1) should be rotated")
    call expect_eq_int(test_err, 5, this%Sonda%collection(1)%cordinates(1)%Ze, "rotate_generateMasSondas: Ze(1) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Xtrancos, "rotate_generateMasSondas: Xtrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Ytrancos, "rotate_generateMasSondas: Ytrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Ztrancos, "rotate_generateMasSondas: Ztrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Or, "rotate_generateMasSondas: Or(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%type1, "rotate_generateMasSondas: type1(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%type2, "rotate_generateMasSondas: type2(1) should remain unchanged")
    
    ! Verify MasSonda 2 (rotated around x-axis)
    call expect_eq_int(test_err, 5, this%Sonda%collection(2)%cordinates(1)%Xi, "rotate_generateMasSondas: Xi(2) should remain unchanged")
    call expect_eq_int(test_err, 8, this%Sonda%collection(2)%cordinates(1)%Xe, "rotate_generateMasSondas: Xe(2) should remain unchanged")
    call expect_eq_int(test_err, -7, this%Sonda%collection(2)%cordinates(1)%Yi, "rotate_generateMasSondas: Yi(2) should be rotated")
    call expect_eq_int(test_err, -10, this%Sonda%collection(2)%cordinates(1)%Ye, "rotate_generateMasSondas: Ye(2) should be rotated")
    call expect_eq_int(test_err, 6, this%Sonda%collection(2)%cordinates(1)%Zi, "rotate_generateMasSondas: Zi(2) should be rotated")
    call expect_eq_int(test_err, 9, this%Sonda%collection(2)%cordinates(1)%Ze, "rotate_generateMasSondas: Ze(2) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Xtrancos, "rotate_generateMasSondas: Xtrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Ytrancos, "rotate_generateMasSondas: Ytrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Ztrancos, "rotate_generateMasSondas: Ztrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%cordinates(1)%Or, "rotate_generateMasSondas: Or(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%type1, "rotate_generateMasSondas: type1(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%type2, "rotate_generateMasSondas: type2(2) should remain unchanged")
end subroutine verify_massonda_rotation_x

subroutine cleanup_massonda_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%Sonda%length > 0) then
        do i = 1, this%Sonda%length
            deallocate(this%Sonda%collection(i)%cordinates)
        end do
        deallocate(this%Sonda%collection)
    end if
    deallocate(this%Sonda)
end subroutine cleanup_massonda_test

end module test_rotate_generateMasSondas_m 