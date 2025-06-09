module test_rotate_generateThinWires_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_thin_wires() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    ! Test case 1: Rotate thin wires in mpidir=2 direction
    mpidir = 2
    call setup_thin_wires_test(this, 2)
    
    ! Set up test data
    call init_thin_wire_data(this, 1, 1, 2, 3, 4, 2, 3, 1, 1, 0.5d0, "wire1")
    call init_thin_wire_data(this, 2, 5, 6, 7, 5, 6, 10, 2, 3, 0.8d0, "wire2")
    
    ! Call rotation and verify results
    call rotate_generateThinWires(this, mpidir)
    call verify_thin_wire_rotation_y(test_err, this)
    
    ! Clean up
    call cleanup_thin_wires_test(this)
    
    ! Test case 2: Rotate thin wires in mpidir=1 direction
    mpidir = 1
    call setup_thin_wires_test(this, 2)
    
    ! Set up test data (same as before)
    call init_thin_wire_data(this, 1, 1, 2, 3, 4, 2, 3, 1, 1, 0.5d0, "wire1")
    call init_thin_wire_data(this, 2, 5, 6, 7, 5, 6, 10, 2, 3, 0.8d0, "wire2")
    
    ! Call rotation and verify results
    call rotate_generateThinWires(this, mpidir)
    call verify_thin_wire_rotation_x(test_err, this)
    
    ! Clean up
    call cleanup_thin_wires_test(this)
    
    ! Test case 3: Verify behavior with no wires
    mpidir = 2
    call setup_thin_wires_test(this, 0)
    call rotate_generateThinWires(this, mpidir)
    call expect_eq_int(test_err, 0, this%tWires%n_tw, "rotate_generateThinWires: n_tw should remain 0")
    call cleanup_thin_wires_test(this)
    
    err = test_err
end function test_rotate_generate_thin_wires

subroutine setup_thin_wires_test(this, n_wires)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_wires
    integer :: i
    
    allocate(this%tWires)
    this%tWires%n_tw = n_wires
    this%tWires%n_tw_max = n_wires
    if (n_wires > 0) then
        allocate(this%tWires%tw(n_wires))
        ! Initialize default values for each wire
        do i = 1, n_wires
            this%tWires%tw(i)%rad = 0.1d0
            this%tWires%tw(i)%disp = .false.
            this%tWires%tw(i)%res = 0.0d0
            this%tWires%tw(i)%ind = 0.0d0
            this%tWires%tw(i)%cap = 0.0d0
            this%tWires%tw(i)%n_twc = 1
            this%tWires%tw(i)%n_twc_max = 1
            ! Allocate components array for each wire
            allocate(this%tWires%tw(i)%twc(1))
        end do
    end if
end subroutine setup_thin_wires_test

subroutine init_thin_wire_data(this, idx, i1, j1, k1, i2, j2, k2, nd, d, m, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, i1, j1, k1, i2, j2, k2, nd, d
    real(8), intent(in) :: m
    character(len=*), intent(in) :: tag
    
    ! Initialize grid indices for the wire component
    this%tWires%tw(idx)%twc(1)%i = i1
    this%tWires%tw(idx)%twc(1)%j = j1
    this%tWires%tw(idx)%twc(1)%k = k1
    this%tWires%tw(idx)%twc(1)%nd = nd
    this%tWires%tw(idx)%twc(1)%d = d
    this%tWires%tw(idx)%twc(1)%m = m
    this%tWires%tw(idx)%twc(1)%tag = tag
    this%tWires%tw(idx)%twc(1)%srctype = "type1"
    this%tWires%tw(idx)%twc(1)%srcfile = "source.dat"
end subroutine init_thin_wire_data

subroutine verify_thin_wire_rotation_y(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify ThinWire 1 (rotated around y-axis)
    call expect_eq_int(test_err, -3, this%tWires%tw(1)%twc(1)%i, "rotate_generateThinWires: i(1) should be rotated")
    call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(1)%j, "rotate_generateThinWires: j(1) should remain unchanged")
    call expect_eq_int(test_err, -1, this%tWires%tw(1)%twc(1)%k, "rotate_generateThinWires: k(1) should be rotated")
    call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%nd, "rotate_generateThinWires: nd(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%d, "rotate_generateThinWires: d(1) should remain unchanged")
    call expect_eq_real(test_err, 0.5_RKIND, this%tWires%tw(1)%twc(1)%m, "rotate_generateThinWires: m(1) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tWires%tw(1)%rad, "rotate_generateThinWires: rad(1) should remain unchanged")
    
    ! Verify ThinWire 2 (rotated around y-axis)
    call expect_eq_int(test_err, -7, this%tWires%tw(2)%twc(1)%i, "rotate_generateThinWires: i(2) should be rotated")
    call expect_eq_int(test_err, 6, this%tWires%tw(2)%twc(1)%j, "rotate_generateThinWires: j(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%tWires%tw(2)%twc(1)%k, "rotate_generateThinWires: k(2) should be rotated")
    call expect_eq_int(test_err, 2, this%tWires%tw(2)%twc(1)%nd, "rotate_generateThinWires: nd(2) should remain unchanged")
    call expect_eq_int(test_err, 3, this%tWires%tw(2)%twc(1)%d, "rotate_generateThinWires: d(2) should remain unchanged")
    call expect_eq_real(test_err, 0.8_RKIND, this%tWires%tw(2)%twc(1)%m, "rotate_generateThinWires: m(2) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tWires%tw(2)%rad, "rotate_generateThinWires: rad(2) should remain unchanged")
end subroutine verify_thin_wire_rotation_y

subroutine verify_thin_wire_rotation_x(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify ThinWire 1 (rotated around x-axis)
    call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%i, "rotate_generateThinWires: i(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%tWires%tw(1)%twc(1)%j, "rotate_generateThinWires: j(1) should be rotated")
    call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(1)%k, "rotate_generateThinWires: k(1) should be rotated")
    call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%nd, "rotate_generateThinWires: nd(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%d, "rotate_generateThinWires: d(1) should remain unchanged")
    call expect_eq_real(test_err, 0.5_RKIND, this%tWires%tw(1)%twc(1)%m, "rotate_generateThinWires: m(1) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tWires%tw(1)%rad, "rotate_generateThinWires: rad(1) should remain unchanged")
    
    ! Verify ThinWire 2 (rotated around x-axis)
    call expect_eq_int(test_err, 5, this%tWires%tw(2)%twc(1)%i, "rotate_generateThinWires: i(2) should remain unchanged")
    call expect_eq_int(test_err, -7, this%tWires%tw(2)%twc(1)%j, "rotate_generateThinWires: j(2) should be rotated")
    call expect_eq_int(test_err, 6, this%tWires%tw(2)%twc(1)%k, "rotate_generateThinWires: k(2) should be rotated")
    call expect_eq_int(test_err, 2, this%tWires%tw(2)%twc(1)%nd, "rotate_generateThinWires: nd(2) should remain unchanged")
    call expect_eq_int(test_err, 3, this%tWires%tw(2)%twc(1)%d, "rotate_generateThinWires: d(2) should remain unchanged")
    call expect_eq_real(test_err, 0.8_RKIND, this%tWires%tw(2)%twc(1)%m, "rotate_generateThinWires: m(2) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tWires%tw(2)%rad, "rotate_generateThinWires: rad(2) should remain unchanged")
end subroutine verify_thin_wire_rotation_x

subroutine cleanup_thin_wires_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%tWires%n_tw > 0) then
        do i = 1, this%tWires%n_tw
            deallocate(this%tWires%tw(i)%twc)
        end do
        deallocate(this%tWires%tw)
    end if
    deallocate(this%tWires)
end subroutine cleanup_thin_wires_test

end module test_rotate_generateThinWires_m 