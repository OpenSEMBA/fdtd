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
    
    ! Test case 1: Rotate Sondas in mpidir=2 direction
    mpidir = 2
    call setup_sonda_test(this, 2)
    
    ! Set up test data - using Electric_Sonda type for testing
    call init_sonda_data(this, 1, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 3.0d0, 0.0d0, 0.0d0, 1.0d0, 1)
    call init_sonda_data(this, 2, 5.0d0, 6.0d0, 7.0d0, 5.0d0, 8.0d0, 9.0d0, 1.0d0, 0.0d0, 0.0d0, 2)
    
    ! Call rotation and verify results
    call rotate_generateSondas(this, mpidir)
    call verify_sonda_rotation_y(test_err, this)
    
    ! Clean up
    call cleanup_sonda_test(this)
    
    ! Test case 2: Rotate Sondas in mpidir=1 direction
    mpidir = 1
    call setup_sonda_test(this, 2)
    
    ! Set up test data (same as before)
    call init_sonda_data(this, 1, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 3.0d0, 0.0d0, 0.0d0, 1.0d0, 1)
    call init_sonda_data(this, 2, 5.0d0, 6.0d0, 7.0d0, 5.0d0, 8.0d0, 9.0d0, 1.0d0, 0.0d0, 0.0d0, 2)
    
    ! Call rotation and verify results
    call rotate_generateSondas(this, mpidir)
    call verify_sonda_rotation_x(test_err, this)
    
    ! Clean up
    call cleanup_sonda_test(this)
    
    ! Test case 3: Verify behavior with no Sondas
    mpidir = 2
    call setup_sonda_test(this, 0)
    call rotate_generateSondas(this, mpidir)
    call expect_eq_int(test_err, 0, this%Sonda%length, "rotate_generateSondas: length should remain 0")
    call cleanup_sonda_test(this)
    
    err = test_err
end function test_rotate_generate_sondas

subroutine setup_sonda_test(this, n_probes)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_probes
    integer :: i
    
    allocate(this%Sonda)
    this%Sonda%length = n_probes
    this%Sonda%length_max = n_probes
    if (n_probes > 0) then
        allocate(this%Sonda%collection(n_probes))
        ! Initialize each probe
        do i = 1, n_probes
            this%Sonda%collection(i)%len_cor = 1
            this%Sonda%len_cor_max = 1
            allocate(this%Sonda%collection(i)%cordinates(1))
        end do
    end if
end subroutine setup_sonda_test

subroutine init_sonda_data(this, idx, x1, y1, z1, x2, y2, z2, dirx, diry, dirz, type_val)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, type_val
    real(8), intent(in) :: x1, y1, z1, x2, y2, z2, dirx, diry, dirz
    
    ! Initialize MasSonda data
    this%Sonda%collection(idx)%cordinates(1)%Xi = int(x1)
    this%Sonda%collection(idx)%cordinates(1)%Yi = int(y1)
    this%Sonda%collection(idx)%cordinates(1)%Zi = int(z1)
    this%Sonda%collection(idx)%cordinates(1)%Xe = int(x2)
    this%Sonda%collection(idx)%cordinates(1)%Ye = int(y2)
    this%Sonda%collection(idx)%cordinates(1)%Ze = int(z2)
    this%Sonda%collection(idx)%cordinates(1)%Xtrancos = int(dirx)
    this%Sonda%collection(idx)%cordinates(1)%Ytrancos = int(diry)
    this%Sonda%collection(idx)%cordinates(1)%Ztrancos = int(dirz)
    this%Sonda%collection(idx)%cordinates(1)%Or = type_val
    this%Sonda%collection(idx)%cordinates(1)%tag = "test_probe"
end subroutine init_sonda_data

subroutine verify_sonda_rotation_y(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify Sonda 1 (rotated around y-axis)
    call expect_eq_int(test_err, -3, this%Sonda%collection(1)%cordinates(1)%Xi, "rotate_generateSondas: Xi(1) should be rotated")
    call expect_eq_int(test_err, 2, this%Sonda%collection(1)%cordinates(1)%Yi, "rotate_generateSondas: Yi(1) should remain unchanged")
    call expect_eq_int(test_err, -1, this%Sonda%collection(1)%cordinates(1)%Zi, "rotate_generateSondas: Zi(1) should be rotated")
    call expect_eq_int(test_err, -3, this%Sonda%collection(1)%cordinates(1)%Xe, "rotate_generateSondas: Xe(1) should be rotated")
    call expect_eq_int(test_err, 5, this%Sonda%collection(1)%cordinates(1)%Ye, "rotate_generateSondas: Ye(1) should remain unchanged")
    call expect_eq_int(test_err, -4, this%Sonda%collection(1)%cordinates(1)%Ze, "rotate_generateSondas: Ze(1) should be rotated")
    call expect_eq_int(test_err, 0, this%Sonda%collection(1)%cordinates(1)%Xtrancos, "rotate_generateSondas: Xtrancos(1) should be rotated")
    call expect_eq_int(test_err, 0, this%Sonda%collection(1)%cordinates(1)%Ytrancos, "rotate_generateSondas: Ytrancos(1) should remain unchanged")
    call expect_eq_int(test_err, -1, this%Sonda%collection(1)%cordinates(1)%Ztrancos, "rotate_generateSondas: Ztrancos(1) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Or, "rotate_generateSondas: Or(1) should remain unchanged")
    
    ! Verify Sonda 2 (rotated around y-axis)
    call expect_eq_int(test_err, -7, this%Sonda%collection(2)%cordinates(1)%Xi, "rotate_generateSondas: Xi(2) should be rotated")
    call expect_eq_int(test_err, 6, this%Sonda%collection(2)%cordinates(1)%Yi, "rotate_generateSondas: Yi(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%Sonda%collection(2)%cordinates(1)%Zi, "rotate_generateSondas: Zi(2) should be rotated")
    call expect_eq_int(test_err, -9, this%Sonda%collection(2)%cordinates(1)%Xe, "rotate_generateSondas: Xe(2) should be rotated")
    call expect_eq_int(test_err, 8, this%Sonda%collection(2)%cordinates(1)%Ye, "rotate_generateSondas: Ye(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%Sonda%collection(2)%cordinates(1)%Ze, "rotate_generateSondas: Ze(2) should be rotated")
    call expect_eq_int(test_err, -1, this%Sonda%collection(2)%cordinates(1)%Xtrancos, "rotate_generateSondas: Xtrancos(2) should be rotated")
    call expect_eq_int(test_err, 0, this%Sonda%collection(2)%cordinates(1)%Ytrancos, "rotate_generateSondas: Ytrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 0, this%Sonda%collection(2)%cordinates(1)%Ztrancos, "rotate_generateSondas: Ztrancos(2) should be rotated")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%cordinates(1)%Or, "rotate_generateSondas: Or(2) should remain unchanged")
end subroutine verify_sonda_rotation_y

subroutine verify_sonda_rotation_x(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify Sonda 1 (rotated around x-axis)
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Xi, "rotate_generateSondas: Xi(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%Sonda%collection(1)%cordinates(1)%Yi, "rotate_generateSondas: Yi(1) should be rotated")
    call expect_eq_int(test_err, 2, this%Sonda%collection(1)%cordinates(1)%Zi, "rotate_generateSondas: Zi(1) should be rotated")
    call expect_eq_int(test_err, 4, this%Sonda%collection(1)%cordinates(1)%Xe, "rotate_generateSondas: Xe(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%Sonda%collection(1)%cordinates(1)%Ye, "rotate_generateSondas: Ye(1) should be rotated")
    call expect_eq_int(test_err, 5, this%Sonda%collection(1)%cordinates(1)%Ze, "rotate_generateSondas: Ze(1) should be rotated")
    call expect_eq_int(test_err, 0, this%Sonda%collection(1)%cordinates(1)%Xtrancos, "rotate_generateSondas: Xtrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 0, this%Sonda%collection(1)%cordinates(1)%Ytrancos, "rotate_generateSondas: Ytrancos(1) should be rotated")
    call expect_eq_int(test_err, -1, this%Sonda%collection(1)%cordinates(1)%Ztrancos, "rotate_generateSondas: Ztrancos(1) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(1)%cordinates(1)%Or, "rotate_generateSondas: Or(1) should remain unchanged")
    
    ! Verify Sonda 2 (rotated around x-axis)
    call expect_eq_int(test_err, 5, this%Sonda%collection(2)%cordinates(1)%Xi, "rotate_generateSondas: Xi(2) should remain unchanged")
    call expect_eq_int(test_err, -7, this%Sonda%collection(2)%cordinates(1)%Yi, "rotate_generateSondas: Yi(2) should be rotated")
    call expect_eq_int(test_err, 6, this%Sonda%collection(2)%cordinates(1)%Zi, "rotate_generateSondas: Zi(2) should be rotated")
    call expect_eq_int(test_err, 5, this%Sonda%collection(2)%cordinates(1)%Xe, "rotate_generateSondas: Xe(2) should remain unchanged")
    call expect_eq_int(test_err, -9, this%Sonda%collection(2)%cordinates(1)%Ye, "rotate_generateSondas: Ye(2) should be rotated")
    call expect_eq_int(test_err, 8, this%Sonda%collection(2)%cordinates(1)%Ze, "rotate_generateSondas: Ze(2) should be rotated")
    call expect_eq_int(test_err, 1, this%Sonda%collection(2)%cordinates(1)%Xtrancos, "rotate_generateSondas: Xtrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 0, this%Sonda%collection(2)%cordinates(1)%Ytrancos, "rotate_generateSondas: Ytrancos(2) should be rotated")
    call expect_eq_int(test_err, 0, this%Sonda%collection(2)%cordinates(1)%Ztrancos, "rotate_generateSondas: Ztrancos(2) should be rotated")
    call expect_eq_int(test_err, 2, this%Sonda%collection(2)%cordinates(1)%Or, "rotate_generateSondas: Or(2) should remain unchanged")
end subroutine verify_sonda_rotation_x

subroutine cleanup_sonda_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%Sonda%length > 0) then
        do i = 1, this%Sonda%length
            if (this%Sonda%collection(i)%len_cor > 0) then
                deallocate(this%Sonda%collection(i)%cordinates)
            end if
        end do
        deallocate(this%Sonda%collection)
    end if
    deallocate(this%Sonda)
end subroutine cleanup_sonda_test

end module test_rotate_generateSONDAs_m 