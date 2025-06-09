integer function test_rotate_generate_current_field_sources() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0  ! Error counter for expect_eq calls
    
    ! Test case 1: X->Y->Z->X rotation (mpidir=2)
    mpidir = 2
    allocate(this%nodsrc)
    this%nodsrc%n_nodSrc = 1
    allocate(this%nodsrc%NodalSource(1))
    
    ! Initialize test data
    this%nodsrc%NodalSource(1)%n_c1P = 1
    this%nodsrc%NodalSource(1)%n_c2P = 1
    allocate(this%nodsrc%NodalSource(1)%c1P(1))
    allocate(this%nodsrc%NodalSource(1)%c2P(1))
    
    this%nodsrc%NodalSource(1)%c1P(1)%XI = 1
    this%nodsrc%NodalSource(1)%c1P(1)%YI = 2
    this%nodsrc%NodalSource(1)%c1P(1)%ZI = 3
    this%nodsrc%NodalSource(1)%c1P(1)%OR = iEx
    
    this%nodsrc%NodalSource(1)%c2P(1)%XI = 4
    this%nodsrc%NodalSource(1)%c2P(1)%YI = 5
    this%nodsrc%NodalSource(1)%c2P(1)%ZI = 6
    this%nodsrc%NodalSource(1)%c2P(1)%OR = iEy
    
    ! Call the routine
    call rotate_generateCurrent_Field_Sources(this, mpidir)
    
    ! Verify results
    call expect_eq_int(test_err, 3, this%nodsrc%NodalSource(1)%c1P(1)%XI, "rotate_generateCurrent_Field_Sources: c1P XI should be 3")
    call expect_eq_int(test_err, 1, this%nodsrc%NodalSource(1)%c1P(1)%YI, "rotate_generateCurrent_Field_Sources: c1P YI should be 1")
    call expect_eq_int(test_err, 2, this%nodsrc%NodalSource(1)%c1P(1)%ZI, "rotate_generateCurrent_Field_Sources: c1P ZI should be 2")
    call expect_eq_int(test_err, iEy, this%nodsrc%NodalSource(1)%c1P(1)%OR, "rotate_generateCurrent_Field_Sources: c1P OR should be iEy")
    
    call expect_eq_int(test_err, 6, this%nodsrc%NodalSource(1)%c2P(1)%XI, "rotate_generateCurrent_Field_Sources: c2P XI should be 6")
    call expect_eq_int(test_err, 4, this%nodsrc%NodalSource(1)%c2P(1)%YI, "rotate_generateCurrent_Field_Sources: c2P YI should be 4")
    call expect_eq_int(test_err, 5, this%nodsrc%NodalSource(1)%c2P(1)%ZI, "rotate_generateCurrent_Field_Sources: c2P ZI should be 5")
    call expect_eq_int(test_err, iEz, this%nodsrc%NodalSource(1)%c2P(1)%OR, "rotate_generateCurrent_Field_Sources: c2P OR should be iEz")
    
    deallocate(this%nodsrc%NodalSource(1)%c1P)
    deallocate(this%nodsrc%NodalSource(1)%c2P)
    deallocate(this%nodsrc%NodalSource)
    deallocate(this%nodsrc)
    
    err = test_err  ! Set the function result to the accumulated error count
end function