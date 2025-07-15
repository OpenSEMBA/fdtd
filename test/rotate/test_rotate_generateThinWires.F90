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
    
    mpidir = 2
    call setup_thin_wires_test(this)
    call init_thin_wire_component_data(this,1, 1, 1, 2, 3, 1, 1, 0.5d0, "wire1")
    call init_thin_wire_component_data(this,1, 2, 5, 6, 7, 2, 3, 0.8d0, "wire2")

    call rotate_generateThinWires(this, mpidir)

    call verify_thin_wire_rotation(test_err, this, mpidir)
    call cleanup_thin_wires_test(this)
    
    mpidir = 1
    call setup_thin_wires_test(this)
    call init_thin_wire_component_data(this,1, 1, 1, 2, 3, 1, 1, 0.5d0, "wire1")
    call init_thin_wire_component_data(this,1, 2, 5, 6, 7, 2, 3, 0.8d0, "wire2")
    
    call rotate_generateThinWires(this, mpidir)

    call verify_thin_wire_rotation(test_err, this, mpidir)
    call cleanup_thin_wires_test(this)
    
    err = test_err
end function test_rotate_generate_thin_wires

subroutine setup_thin_wires_test(this)
    type(Parseador), intent(inout) :: this
    
    allocate(this%tWires)
    this%tWires%n_tw = 1
    allocate(this%tWires%tw(1))
    this%tWires%tw(1)%rad = 0.1_RKIND
    this%tWires%tw(1)%disp = .false.
    this%tWires%tw(1)%res = 0.0_RKIND
    this%tWires%tw(1)%ind = 0.0_RKIND
    this%tWires%tw(1)%cap = 0.0_RKIND
    this%tWires%tw(1)%n_twc = 2
    this%tWires%tw(1)%n_twc_max = 2    
    allocate(this%tWires%tw(1)%twc(1))
    allocate(this%tWires%tw(1)%twc(2))
end subroutine setup_thin_wires_test

subroutine init_thin_wire_component_data(this, twidx, twcidx, i1, j1, k1, nd, d, m, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: twidx, twcidx, i1, j1, k1, nd, d
    real(KIND=8), intent(in) :: m
    character(len=*), intent(in) :: tag

    this%tWires%tw(twidx)%twc(twcidx)%i = i1
    this%tWires%tw(twidx)%twc(twcidx)%j = j1
    this%tWires%tw(twidx)%twc(twcidx)%k = k1
    this%tWires%tw(twidx)%twc(twcidx)%nd = nd
    this%tWires%tw(twidx)%twc(twcidx)%d = d
    this%tWires%tw(twidx)%twc(twcidx)%m = m
    this%tWires%tw(twidx)%twc(twcidx)%tag = tag
    this%tWires%tw(twidx)%twc(twcidx)%srctype = "type1"
    this%tWires%tw(twidx)%twc(twcidx)%srcfile = "source.dat"
end subroutine init_thin_wire_component_data

subroutine verify_thin_wire_rotation(test_err, this, mpidir)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    integer, intent(in) :: mpidir
    if (mpidir==2) then
        call expect_eq_int(test_err, 3, this%tWires%tw(1)%twc(1)%i, "rotate_generateThinWires: i(1) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%j, "rotate_generateThinWires: j(1) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(1)%k, "rotate_generateThinWires: k(1) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%nd, "rotate_generateThinWires: nd(1) should remain unchanged")
        call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(1)%d, "rotate_generateThinWires: d(1) should change to iEy")
        call expect_eq_real(test_err, 0.5_RKIND, this%tWires%tw(1)%twc(1)%m, "rotate_generateThinWires: m(1) should remain unchanged")
        !
        call expect_eq_int(test_err, 7, this%tWires%tw(1)%twc(2)%i, "rotate_generateThinWires: i(2) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 5, this%tWires%tw(1)%twc(2)%j, "rotate_generateThinWires: j(2) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 6, this%tWires%tw(1)%twc(2)%k, "rotate_generateThinWires: k(2) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(2)%nd, "rotate_generateThinWires: nd(2) should remain unchanged")
        call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(2)%d, "rotate_generateThinWires: d(2) should change to iEx")
        call expect_eq_real(test_err, 0.8_RKIND, this%tWires%tw(1)%twc(2)%m, "rotate_generateThinWires: m(2) should remain unchanged")

    elseif(mpidir==1) then
        call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(1)%i, "rotate_generateThinWires: i(1) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 3, this%tWires%tw(1)%twc(1)%j, "rotate_generateThinWires: j(1) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%k, "rotate_generateThinWires: k(1) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 1, this%tWires%tw(1)%twc(1)%nd, "rotate_generateThinWires: nd(1) should remain unchanged")
        call expect_eq_int(test_err, 3, this%tWires%tw(1)%twc(1)%d, "rotate_generateThinWires: d(1) should change to iEz")
        call expect_eq_real(test_err, 0.5_RKIND, this%tWires%tw(1)%twc(1)%m, "rotate_generateThinWires: m(1) should remain unchanged")
        !
        call expect_eq_int(test_err, 6, this%tWires%tw(1)%twc(2)%i, "rotate_generateThinWires: i(2) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 7, this%tWires%tw(1)%twc(2)%j, "rotate_generateThinWires: j(2) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 5, this%tWires%tw(1)%twc(2)%k, "rotate_generateThinWires: k(2) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(2)%nd, "rotate_generateThinWires: nd(2) should remain unchanged")
        call expect_eq_int(test_err, 2, this%tWires%tw(1)%twc(2)%d, "rotate_generateThinWires: d(2) should change to iEy")
        call expect_eq_real(test_err, 0.8_RKIND, this%tWires%tw(1)%twc(2)%m, "rotate_generateThinWires: m(2) should remain unchanged")
    endif
end subroutine verify_thin_wire_rotation

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