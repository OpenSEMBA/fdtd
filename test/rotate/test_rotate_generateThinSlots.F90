module test_rotate_generateThinSlots_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_thin_slots() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    ! Test case 1: Rotate thin slots in mpidir=2 direction
    mpidir = 2
    call setup_thin_slot_test(this, 2)
    
    ! Set up test data
    call init_thin_slot_data(this, 1, 1, 2, 3, 1, 1, "slot1")
    call init_thin_slot_data(this, 2, 4, 5, 6, 2, 2, "slot2")
    
    ! Call rotation and verify results
    call rotate_generateThinSlots(this, mpidir)
    call verify_thin_slot_rotation_y(test_err, this)
    
    ! Clean up
    call cleanup_thin_slot_test(this)
    
    ! Test case 2: Rotate thin slots in mpidir=1 direction
    mpidir = 1
    call setup_thin_slot_test(this, 2)
    
    ! Set up test data (same as before)
    call init_thin_slot_data(this, 1, 1, 2, 3, 1, 1, "slot1")
    call init_thin_slot_data(this, 2, 4, 5, 6, 2, 2, "slot2")
    
    ! Call rotation and verify results
    call rotate_generateThinSlots(this, mpidir)
    call verify_thin_slot_rotation_x(test_err, this)
    
    ! Clean up
    call cleanup_thin_slot_test(this)
    
    ! Test case 3: Verify behavior with no thin slots
    mpidir = 2
    call setup_thin_slot_test(this, 0)
    call rotate_generateThinSlots(this, mpidir)
    call expect_eq_int(test_err, 0, this%tslots%n_tg, "rotate_generateThinSlots: n_tg should remain 0")
    call cleanup_thin_slot_test(this)
    
    err = test_err
end function test_rotate_generate_thin_slots

subroutine setup_thin_slot_test(this, n_slots)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_slots
    integer :: i
    
    allocate(this%tslots)
    this%tslots%n_tg = n_slots
    this%tslots%n_tg_max = n_slots
    if (n_slots > 0) then
        allocate(this%tslots%tg(n_slots))
        do i = 1, n_slots
            this%tslots%tg(i)%n_tgc = 1
            this%tslots%tg(i)%n_tgc_max = 1
            allocate(this%tslots%tg(i)%tgc(1))
            this%tslots%tg(i)%width = 0.1d0  ! Example width
        end do
    end if
end subroutine setup_thin_slot_test

subroutine init_thin_slot_data(this, idx, i, j, k, dir, node, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, i, j, k, dir, node
    character(len=*), intent(in) :: tag
    
    ! Initialize ThinSlotComp attributes
    this%tslots%tg(idx)%tgc(1)%i = i
    this%tslots%tg(idx)%tgc(1)%j = j
    this%tslots%tg(idx)%tgc(1)%k = k
    this%tslots%tg(idx)%tgc(1)%dir = dir
    this%tslots%tg(idx)%tgc(1)%node = node
    this%tslots%tg(idx)%tgc(1)%Or = 1  ! Default orientation
    this%tslots%tg(idx)%tgc(1)%tag = tag
end subroutine init_thin_slot_data

subroutine verify_thin_slot_rotation_y(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify ThinSlot 1 (rotated around y-axis)
    call expect_eq_int(test_err, -1, this%tslots%tg(1)%tgc(1)%i, "rotate_generateThinSlots: i(1) should be rotated")
    call expect_eq_int(test_err, 2, this%tslots%tg(1)%tgc(1)%j, "rotate_generateThinSlots: j(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%tslots%tg(1)%tgc(1)%k, "rotate_generateThinSlots: k(1) should be rotated")
    call expect_eq_int(test_err, -1, this%tslots%tg(1)%tgc(1)%dir, "rotate_generateThinSlots: dir(1) should be rotated")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%node, "rotate_generateThinSlots: node(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%Or, "rotate_generateThinSlots: Or(1) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tslots%tg(1)%width, "rotate_generateThinSlots: width(1) should remain unchanged")
    
    ! Verify ThinSlot 2 (rotated around y-axis)
    call expect_eq_int(test_err, -4, this%tslots%tg(2)%tgc(1)%i, "rotate_generateThinSlots: i(2) should be rotated")
    call expect_eq_int(test_err, 5, this%tslots%tg(2)%tgc(1)%j, "rotate_generateThinSlots: j(2) should remain unchanged")
    call expect_eq_int(test_err, -6, this%tslots%tg(2)%tgc(1)%k, "rotate_generateThinSlots: k(2) should be rotated")
    call expect_eq_int(test_err, -2, this%tslots%tg(2)%tgc(1)%dir, "rotate_generateThinSlots: dir(2) should be rotated")
    call expect_eq_int(test_err, 2, this%tslots%tg(2)%tgc(1)%node, "rotate_generateThinSlots: node(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tslots%tg(2)%tgc(1)%Or, "rotate_generateThinSlots: Or(2) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tslots%tg(2)%width, "rotate_generateThinSlots: width(2) should remain unchanged")
end subroutine verify_thin_slot_rotation_y

subroutine verify_thin_slot_rotation_x(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify ThinSlot 1 (rotated around x-axis)
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%i, "rotate_generateThinSlots: i(1) should remain unchanged")
    call expect_eq_int(test_err, -2, this%tslots%tg(1)%tgc(1)%j, "rotate_generateThinSlots: j(1) should be rotated")
    call expect_eq_int(test_err, 3, this%tslots%tg(1)%tgc(1)%k, "rotate_generateThinSlots: k(1) should be rotated")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%dir, "rotate_generateThinSlots: dir(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%node, "rotate_generateThinSlots: node(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tslots%tg(1)%tgc(1)%Or, "rotate_generateThinSlots: Or(1) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tslots%tg(1)%width, "rotate_generateThinSlots: width(1) should remain unchanged")
    
    ! Verify ThinSlot 2 (rotated around x-axis)
    call expect_eq_int(test_err, 4, this%tslots%tg(2)%tgc(1)%i, "rotate_generateThinSlots: i(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%tslots%tg(2)%tgc(1)%j, "rotate_generateThinSlots: j(2) should be rotated")
    call expect_eq_int(test_err, 6, this%tslots%tg(2)%tgc(1)%k, "rotate_generateThinSlots: k(2) should be rotated")
    call expect_eq_int(test_err, 2, this%tslots%tg(2)%tgc(1)%dir, "rotate_generateThinSlots: dir(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%tslots%tg(2)%tgc(1)%node, "rotate_generateThinSlots: node(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%tslots%tg(2)%tgc(1)%Or, "rotate_generateThinSlots: Or(2) should remain unchanged")
    call expect_eq_real(test_err, 0.1_RKIND, this%tslots%tg(2)%width, "rotate_generateThinSlots: width(2) should remain unchanged")
end subroutine verify_thin_slot_rotation_x

subroutine cleanup_thin_slot_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%tslots%n_tg > 0) then
        do i = 1, this%tslots%n_tg
            deallocate(this%tslots%tg(i)%tgc)
        end do
        deallocate(this%tslots%tg)
    end if
    deallocate(this%tslots)
end subroutine cleanup_thin_slot_test

end module test_rotate_generateThinSlots_m 