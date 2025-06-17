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
    
    mpidir = 2
    call setup_slanted_wires_case(this)

    call rotate_generateSlantedWires(this, mpidir)
    call verify_slanted_wire_rotation(test_err, this, mpidir)
    
    call cleanup_slanted_wires_test(this)
    
    mpidir = 1
    call setup_slanted_wires_case(this)
    
    call rotate_generateSlantedWires(this, mpidir)
    call verify_slanted_wire_rotation(test_err, this, mpidir)
    
    call cleanup_slanted_wires_test(this)
    
    err = test_err
end function test_rotate_generate_slanted_wires

subroutine setup_slanted_wires_case(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    integer :: n_swires = 1, n_swc = 2
    
    allocate(this%sWires)
    this%sWires%n_sw = n_swires
    this%sWires%n_sw_max = n_swires
    allocate(this%sWires%sw(n_swires))
    do i=1, n_swires
        this%sWires%sw(i)%rad = 0.1d0
        this%sWires%sw(i)%disp = .false.
        this%sWires%sw(i)%res = 0.0d0
        this%sWires%sw(i)%ind = 0.0d0
        this%sWires%sw(i)%cap = 0.0d0
        this%sWires%sw(i)%n_swc = n_swc
        this%sWires%sw(i)%n_swc_max = n_swc
        allocate(this%sWires%sw(i)%swc(n_swc))
    end do

    call init_slanted_wire_data(this,1, 1, 1.0d0, 2.0d0, 3.0d0, 1, 0.5d0, "wire1")
    call init_slanted_wire_data(this,1, 2, 5.0d0, 6.0d0, 7.0d0, 2, 0.8d0, "wire2")

end subroutine setup_slanted_wires_case

subroutine init_slanted_wire_data(this, swidx, swcidx, x, y, z, nd, m, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: swidx, swcidx, nd
    real(8), intent(in) :: x, y, z, m
    character(len=*), intent(in) :: tag
    
    this%sWires%sw(swidx)%swc(swcidx)%x = x
    this%sWires%sw(swidx)%swc(swcidx)%y = y
    this%sWires%sw(swidx)%swc(swcidx)%z = z
    this%sWires%sw(swidx)%swc(swcidx)%nd = nd
    this%sWires%sw(swidx)%swc(swcidx)%m = m
    this%sWires%sw(swidx)%swc(swcidx)%tag = tag
    this%sWires%sw(swidx)%swc(swcidx)%srctype = "type1"
    this%sWires%sw(swidx)%swc(swcidx)%srcfile = "source.dat"
end subroutine init_slanted_wire_data

subroutine verify_slanted_wire_rotation(test_err, this, mpidir)
    integer, intent(inout) :: test_err, mpidir
    type(Parseador), intent(in) :: this
    
    if (mpidir==2) then
        call expect_eq_real(test_err, 3.0_RKIND, this%sWires%sw(1)%swc(1)%x, "rotate_generateSlantedWires: x(1) should be rotated with mpidir 2")
        call expect_eq_real(test_err, 1.0_RKIND, this%sWires%sw(1)%swc(1)%y, "rotate_generateSlantedWires: y(1) should be rotated with mpidir 2")
        call expect_eq_real(test_err, 2.0_RKIND, this%sWires%sw(1)%swc(1)%z, "rotate_generateSlantedWires: z(1) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 1, this%sWires%sw(1)%swc(1)%nd, "rotate_generateSlantedWires: nd(1) should remain unchanged")
        call expect_eq_real(test_err, 0.5_RKIND, this%sWires%sw(1)%swc(1)%m, "rotate_generateSlantedWires: m(1) should remain unchanged")
        
        call expect_eq_real(test_err, 7.0_RKIND, this%sWires%sw(1)%swc(2)%x, "rotate_generateSlantedWires: x(2) should be rotated with mpidir 2")
        call expect_eq_real(test_err, 5.0_RKIND, this%sWires%sw(1)%swc(2)%y, "rotate_generateSlantedWires: y(2) should be rotated with mpidir 2")
        call expect_eq_real(test_err, 6.0_RKIND, this%sWires%sw(1)%swc(2)%z, "rotate_generateSlantedWires: z(2) should be rotated with mpidir 2")
        call expect_eq_int(test_err, 2, this%sWires%sw(1)%swc(2)%nd, "rotate_generateSlantedWires: nd(2) should remain unchanged")
        call expect_eq_real(test_err, 0.8_RKIND, this%sWires%sw(1)%swc(2)%m, "rotate_generateSlantedWires: m(2) should remain unchanged")

    else if (mpidir==1) then
        call expect_eq_real(test_err, 2.0_RKIND, this%sWires%sw(1)%swc(1)%x, "rotate_generateSlantedWires: x(1) should be rotated with mpidir 1")
        call expect_eq_real(test_err, 3.0_RKIND, this%sWires%sw(1)%swc(1)%y, "rotate_generateSlantedWires: y(1) should be rotated with mpidir 1")
        call expect_eq_real(test_err, 1.0_RKIND, this%sWires%sw(1)%swc(1)%z, "rotate_generateSlantedWires: z(1) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 1, this%sWires%sw(1)%swc(1)%nd, "rotate_generateSlantedWires: nd(1) should remain unchanged")
        call expect_eq_real(test_err, 0.5_RKIND, this%sWires%sw(1)%swc(1)%m, "rotate_generateSlantedWires: m(1) should remain unchanged")

        call expect_eq_real(test_err, 6.0_RKIND, this%sWires%sw(1)%swc(2)%x, "rotate_generateSlantedWires: x(2) should be rotated with mpidir 1")
        call expect_eq_real(test_err, 7.0_RKIND, this%sWires%sw(1)%swc(2)%y, "rotate_generateSlantedWires: y(2) should be rotated with mpidir 1")
        call expect_eq_real(test_err, 5.0_RKIND, this%sWires%sw(1)%swc(2)%z, "rotate_generateSlantedWires: z(2) should be rotated with mpidir 1")
        call expect_eq_int(test_err, 2, this%sWires%sw(1)%swc(2)%nd, "rotate_generateSlantedWires: nd(2) should remain unchanged")
        call expect_eq_real(test_err, 0.8_RKIND, this%sWires%sw(1)%swc(2)%m, "rotate_generateSlantedWires: m(2) should remain unchanged")
    end if
end subroutine verify_slanted_wire_rotation

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