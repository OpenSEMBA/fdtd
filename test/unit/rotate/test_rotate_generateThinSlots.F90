module test_rotate_generateThinSlots_m
    use smbjson_m
    use nfde_rotate_m
    use rotate_testingTools_m
    implicit none

contains

integer function test_rotate_generate_thin_slots() bind(C) result(err)
    type(Parseador_t) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    mpidir = 2
    call setup_thin_slot_test(this)
    
    call rotate_generateThinSlots(this, mpidir)
    call verify_thin_slot_rotation(test_err, this, mpidir)
    
    call cleanup_thin_slot_test(this)
    
    mpidir = 1
    call setup_thin_slot_test(this)
    
    call rotate_generateThinSlots(this, mpidir)
    call verify_thin_slot_rotation(test_err, this, mpidir)
    
    call cleanup_thin_slot_test(this)
    
    err = test_err
end function test_rotate_generate_thin_slots

subroutine setup_thin_slot_test(this)
    type(Parseador_t), intent(inout) :: this
    integer :: n_slots=1
    integer :: i
    
    allocate(this%tslots)
    this%tslots%n_tg = n_slots
    this%tslots%n_tg_max = n_slots
    
    allocate(this%tslots%tg(n_slots))
    allocate(this%tslots%tg(n_slots)%tgc(2))
    this%tslots%tg(1)%n_tgc = 2
    this%tslots%tg(1)%n_tgc_max = 2
    this%tslots%tg(1)%width = 0.1d0
    call init_thin_slot_component(this, 1, 1, 1, 2, 3, 1, "slot1")
    call init_thin_slot_component(this, 1, 2, 4, 5, 6, 2, "slot1")
end subroutine setup_thin_slot_test

subroutine init_thin_slot_component(this, tgidx, tgcidx, i, j, k, dir, tag)
    type(Parseador_t), intent(inout) :: this
    integer, intent(in) :: tgidx, tgcidx, i, j, k, dir
    character(len=*), intent(in) :: tag
    
    
    this%tslots%tg(tgidx)%tgc(tgcidx)%i = i
    this%tslots%tg(tgidx)%tgc(tgcidx)%j = j
    this%tslots%tg(tgidx)%tgc(tgcidx)%k = k
    this%tslots%tg(tgidx)%tgc(tgcidx)%dir = dir
    this%tslots%tg(tgidx)%tgc(tgcidx)%Or = dir
    this%tslots%tg(tgidx)%tgc(tgcidx)%tag = tag
end subroutine init_thin_slot_component

subroutine verify_thin_slot_rotation(test_err, this, mpidir)
    integer, intent(inout) :: test_err, mpidir
    type(Parseador_t), intent(in) :: this
    if (mpidir==2) then
        call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(1)%i, "rotate_generateThinSlots: i(1) mpidir2 error")
        call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%j, "rotate_generateThinSlots: j(1) mpidir2 error")
        call expect_eq_int(test_err, 2, this%tslots%tg(1)%tgc(1)%k, "rotate_generateThinSlots: k(1) mpidir2 error")
        call expect_eq_int(test_err, 2, this%tslots%tg(1)%tgc(1)%dir, "rotate_generateThinSlots: dir(1) mpidir2 error")
        call expect_eq_int(test_err, 2, this%tslots%tg(1)%tgc(1)%Or, "rotate_generateThinSlots: Or(1) mpidir2 error")

        call expect_eq_int(test_err, 6, this%tslots%tg(1)%tgc(2)%i, "rotate_generateThinSlots: i(2) mpidir2 error")
        call expect_eq_int(test_err, 4, this%tslots%tg(1)%tgc(2)%j, "rotate_generateThinSlots: j(2) mpidir2 error")
        call expect_eq_int(test_err, 5, this%tslots%tg(1)%tgc(2)%k, "rotate_generateThinSlots: k(2) mpidir2 error")
        call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(2)%dir, "rotate_generateThinSlots: dir(2) mpidir2 error")
        call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(2)%Or, "rotate_generateThinSlots: Or(2) mpidir2 error")
    else if (mpidir==1) then
        call expect_eq_int(test_err, 2, this%tslots%tg(1)%tgc(1)%i, "rotate_generateThinSlots: i(1) mpidir1 error")
        call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(1)%j, "rotate_generateThinSlots: j(1) mpidir1 error")
        call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%k, "rotate_generateThinSlots: k(1) mpidir1 error")
        call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(1)%dir, "rotate_generateThinSlots: dir(1) mpidir1 error")
        call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(1)%Or, "rotate_generateThinSlots: Or(1) mpidir1 error")

        call expect_eq_int(test_err, 5, this%tslots%tg(1)%tgc(2)%i, "rotate_generateThinSlots: i(2) mpidir1 error")
        call expect_eq_int(test_err, 6, this%tslots%tg(1)%tgc(2)%j, "rotate_generateThinSlots: j(2) mpidir1 error")
        call expect_eq_int(test_err, 4, this%tslots%tg(1)%tgc(2)%k, "rotate_generateThinSlots: k(2) mpidir1 error")
        call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(2)%dir, "rotate_generateThinSlots: dir(2) mpidir1 error")
        call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(2)%Or, "rotate_generateThinSlots: Or(2) mpidir1 error")
    end if
end subroutine verify_thin_slot_rotation

subroutine cleanup_thin_slot_test(this)
    type(Parseador_t), intent(inout) :: this
    integer :: i
    
    if (this%tslots%n_tg > 0) then
        do i = 1, this%tslots%n_tg
            if (associated(this%tslots%tg(i)%tgc)) deallocate(this%tslots%tg(i)%tgc)
            if (allocated(this%tslots%tg(i)%surfaceCoords)) deallocate(this%tslots%tg(i)%surfaceCoords)
        end do
        deallocate(this%tslots%tg)
    end if
    deallocate(this%tslots)
end subroutine cleanup_thin_slot_test

integer function test_rotate_generate_thin_slots_surface() bind(C) result(err)
    type(Parseador_t) :: this
    integer :: test_err = 0

    call setup_surface_thin_slot_test(this)
    call rotate_generateThinSlots(this, 2)
    call expect_eq_int(test_err, 3, this%tslots%tg(1)%surfaceCoords(1)%Xi, "surface rotate: i mpidir2 error")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%surfaceCoords(1)%Yi, "surface rotate: j mpidir2 error")
    call expect_eq_int(test_err, 2, this%tslots%tg(1)%surfaceCoords(1)%Zi, "surface rotate: k mpidir2 error")
    call expect_eq_int(test_err, +IEX, this%tslots%tg(1)%surfaceCoords(1)%Or, "surface rotate: Or mpidir2 error")
    call cleanup_thin_slot_test(this)

    call setup_surface_thin_slot_test(this)
    call rotate_generateThinSlots(this, 1)
    call expect_eq_int(test_err, 2, this%tslots%tg(1)%surfaceCoords(1)%Xi, "surface rotate: i mpidir1 error")
    call expect_eq_int(test_err, 3, this%tslots%tg(1)%surfaceCoords(1)%Yi, "surface rotate: j mpidir1 error")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%surfaceCoords(1)%Zi, "surface rotate: k mpidir1 error")
    call expect_eq_int(test_err, +IEY, this%tslots%tg(1)%surfaceCoords(1)%Or, "surface rotate: Or mpidir1 error")
    call cleanup_thin_slot_test(this)

    err = test_err
end function test_rotate_generate_thin_slots_surface

subroutine setup_surface_thin_slot_test(this)
    type(Parseador_t), intent(inout) :: this

    allocate(this%tslots)
    this%tslots%n_tg = 1
    this%tslots%n_tg_max = 1
    allocate(this%tslots%tg(1))
    this%tslots%tg(1)%n_tgc = 0
    this%tslots%tg(1)%n_tgc_max = 0
    this%tslots%tg(1)%width = 0.1d0
    allocate(this%tslots%tg(1)%surfaceCoords(1))
    this%tslots%tg(1)%surfaceCoords(1)%Xi = 1
    this%tslots%tg(1)%surfaceCoords(1)%Xe = 1
    this%tslots%tg(1)%surfaceCoords(1)%Yi = 2
    this%tslots%tg(1)%surfaceCoords(1)%Ye = 2
    this%tslots%tg(1)%surfaceCoords(1)%Zi = 3
    this%tslots%tg(1)%surfaceCoords(1)%Ze = 3
    this%tslots%tg(1)%surfaceCoords(1)%Or = +IEZ
    this%tslots%tg(1)%surfaceCoords(1)%tag = 'surface-slot'
end subroutine setup_surface_thin_slot_test

end module test_rotate_generateThinSlots_m 