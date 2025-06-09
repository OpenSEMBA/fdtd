module test_rotate_generateSlantedWires_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_slanted_wires() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    ! Test case 1: Rotate slanted wires in mpidir=2 direction
    mpidir = 2
    call setup_slanted_wires_test(this, 2)
    
    ! Set up test data
    call init_slanted_wire_data(this, 1, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 3.0d0, 1, 0.5d0, "wire1")
    call init_slanted_wire_data(this, 2, 5.0d0, 6.0d0, 7.0d0, 5.0d0, 8.0d0, 10.0d0, 2, 0.8d0, "wire2")
    
    ! Call rotation and verify results
    call rotate_generateSlantedWires(this, mpidir)
    call verify_slanted_wire_rotation_y(test_err, this)
    
    ! Clean up
    call cleanup_slanted_wires_test(this)
    
    ! Test case 2: Rotate slanted wires in mpidir=1 direction
    mpidir = 1
    call setup_slanted_wires_test(this, 2)
    
    ! Set up test data (same as before)
    call init_slanted_wire_data(this, 1, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 3.0d0, 1, 0.5d0, "wire1")
    call init_slanted_wire_data(this, 2, 5.0d0, 6.0d0, 7.0d0, 5.0d0, 8.0d0, 10.0d0, 2, 0.8d0, "wire2")
    
    ! Call rotation and verify results
    call rotate_generateSlantedWires(this, mpidir)
    call verify_slanted_wire_rotation_x(test_err, this)
    
    ! Clean up
    call cleanup_slanted_wires_test(this)
    
    ! Test case 3: Verify behavior with no wires
    mpidir = 2
    call setup_slanted_wires_test(this, 0)
    call rotate_generateSlantedWires(this, mpidir)
    call expect_eq_int(test_err, 0, this%sWires%n_sw, "rotate_generateSlantedWires: n_sw should remain 0")
    call cleanup_slanted_wires_test(this)
    
    err = test_err
end function test_rotate_generate_slanted_wires

subroutine setup_slanted_wires_test(this, n_wires)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_wires
    integer :: i
    
    allocate(this%sWires)
    this%sWires%n_sw = n_wires
    this%sWires%n_sw_max = n_wires
    if (n_wires > 0) then
        allocate(this%sWires%sw(n_wires))
        ! Initialize default values for each wire
        do i = 1, n_wires
            this%sWires%sw(i)%rad = 0.1d0
            this%sWires%sw(i)%disp = .false.
            this%sWires%sw(i)%res = 0.0d0
            this%sWires%sw(i)%ind = 0.0d0
            this%sWires%sw(i)%cap = 0.0d0
            this%sWires%sw(i)%n_swc = 1
            this%sWires%sw(i)%n_swc_max = 1
            ! Allocate components array for each wire
            allocate(this%sWires%sw(i)%swc(1))
        end do
    end if
end subroutine setup_slanted_wires_test

subroutine init_slanted_wire_data(this, idx, x1, y1, z1, x2, y2, z2, nd, m, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, nd
    real(8), intent(in) :: x1, y1, z1, x2, y2, z2, m
    character(len=*), intent(in) :: tag
    
    ! Initialize coordinates for the wire component
    this%sWires%sw(idx)%swc(1)%x = x1
    this%sWires%sw(idx)%swc(1)%y = y1
    this%sWires%sw(idx)%swc(1)%z = z1
    this%sWires%sw(idx)%swc(1)%nd = nd
    this%sWires%sw(idx)%swc(1)%m = m
    this%sWires%sw(idx)%swc(1)%tag = tag
    this%sWires%sw(idx)%swc(1)%srctype = "type1"
    this%sWires%sw(idx)%swc(1)%srcfile = "source.dat"
end subroutine init_slanted_wire_data

subroutine verify_slanted_wire_rotation_y(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify SlantedWire 1 (rotated around y-axis)
    call expect_eq_real(test_err, -3.0_RKIND, this%sWires%sw(1)%swc(1)%x, "rotate_generateSlantedWires: x(1) should be rotated")
    call expect_eq_real(test_err, 2.0_RKIND, this%sWires%sw(1)%swc(1)%y, "rotate_generateSlantedWires: y(1) should remain unchanged")
    call expect_eq_real(test_err, -1.0_RKIND, this%sWires%sw(1)%swc(1)%z, "rotate_generateSlantedWires: z(1) should be rotated")
    call expect_eq_int(test_err, 1, this%sWires%sw(1)%swc(1)%nd, "rotate_generateSlantedWires: nd(1) should remain unchanged")
    call expect_eq_real(test_err, 0.5_RKIND, this%sWires%sw(1)%swc(1)%m, "rotate_generateSlantedWires: m(1) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%sWires%sw(1)%rad, "rotate_generateSlantedWires: rad(1) should remain unchanged")
    
    ! Verify SlantedWire 2 (rotated around y-axis)
    call expect_eq_real(test_err, -7.0_RKIND, this%sWires%sw(2)%swc(1)%x, "rotate_generateSlantedWires: x(2) should be rotated")
    call expect_eq_real(test_err, 6.0_RKIND, this%sWires%sw(2)%swc(1)%y, "rotate_generateSlantedWires: y(2) should remain unchanged")
    call expect_eq_real(test_err, -5.0_RKIND, this%sWires%sw(2)%swc(1)%z, "rotate_generateSlantedWires: z(2) should be rotated")
    call expect_eq_int(test_err, 2, this%sWires%sw(2)%swc(1)%nd, "rotate_generateSlantedWires: nd(2) should remain unchanged")
    call expect_eq_real(test_err, 0.8_RKIND, this%sWires%sw(2)%swc(1)%m, "rotate_generateSlantedWires: m(2) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%sWires%sw(2)%rad, "rotate_generateSlantedWires: rad(2) should remain unchanged")
end subroutine verify_slanted_wire_rotation_y

subroutine verify_slanted_wire_rotation_x(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify SlantedWire 1 (rotated around x-axis)
    call expect_eq_real(test_err, 1.0_RKIND, this%sWires%sw(1)%swc(1)%x, "rotate_generateSlantedWires: x(1) should remain unchanged")
    call expect_eq_real(test_err, -3.0_RKIND, this%sWires%sw(1)%swc(1)%y, "rotate_generateSlantedWires: y(1) should be rotated")
    call expect_eq_real(test_err, 2.0_RKIND, this%sWires%sw(1)%swc(1)%z, "rotate_generateSlantedWires: z(1) should be rotated")
    call expect_eq_int(test_err, 1, this%sWires%sw(1)%swc(1)%nd, "rotate_generateSlantedWires: nd(1) should remain unchanged")
    call expect_eq_real(test_err, 0.5_RKIND, this%sWires%sw(1)%swc(1)%m, "rotate_generateSlantedWires: m(1) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%sWires%sw(1)%rad, "rotate_generateSlantedWires: rad(1) should remain unchanged")
    
    ! Verify SlantedWire 2 (rotated around x-axis)
    call expect_eq_real(test_err, 5.0_RKIND, this%sWires%sw(2)%swc(1)%x, "rotate_generateSlantedWires: x(2) should remain unchanged")
    call expect_eq_real(test_err, -7.0_RKIND, this%sWires%sw(2)%swc(1)%y, "rotate_generateSlantedWires: y(2) should be rotated")
    call expect_eq_real(test_err, 6.0_RKIND, this%sWires%sw(2)%swc(1)%z, "rotate_generateSlantedWires: z(2) should be rotated")
    call expect_eq_int(test_err, 2, this%sWires%sw(2)%swc(1)%nd, "rotate_generateSlantedWires: nd(2) should remain unchanged")
    call expect_eq_real(test_err, 0.8_RKIND, this%sWires%sw(2)%swc(1)%m, "rotate_generateSlantedWires: m(2) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%sWires%sw(2)%rad, "rotate_generateSlantedWires: rad(2) should remain unchanged")
end subroutine verify_slanted_wire_rotation_x

subroutine cleanup_slanted_wires_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%sWires%n_sw > 0) then
        do i = 1, this%sWires%n_sw
            deallocate(this%sWires%sw(i)%swc)
        end do
        deallocate(this%sWires%sw)
    end if
    deallocate(this%sWires)
end subroutine cleanup_slanted_wires_test

end module test_rotate_generateSlantedWires_m 