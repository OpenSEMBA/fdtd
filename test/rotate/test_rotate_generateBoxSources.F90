integer function test_rotate_generate_box_sources() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    mpidir = 2
    allocate(this%boxSrc)
    this%boxSrc%nvols = 1
    allocate(this%boxSrc%vols(1))
    
    this%boxSrc%vols(1)%coor1(1) = 1
    this%boxSrc%vols(1)%coor1(2) = 2
    this%boxSrc%vols(1)%coor1(3) = 3
    this%boxSrc%vols(1)%coor2(1) = 4
    this%boxSrc%vols(1)%coor2(2) = 5
    this%boxSrc%vols(1)%coor2(3) = 6
    
    call rotate_generateBoxSources(this, mpidir)
    
    call expect_eq_int(test_err, 3, this%boxSrc%vols(1)%coor1(1), "rotate_generateBoxSources: coor1(1) should be 3")
    call expect_eq_int(test_err, 1, this%boxSrc%vols(1)%coor1(2), "rotate_generateBoxSources: coor1(2) should be 1")
    call expect_eq_int(test_err, 2, this%boxSrc%vols(1)%coor1(3), "rotate_generateBoxSources: coor1(3) should be 2")
    
    call expect_eq_int(test_err, 6, this%boxSrc%vols(1)%coor2(1), "rotate_generateBoxSources: coor2(1) should be 6")
    call expect_eq_int(test_err, 4, this%boxSrc%vols(1)%coor2(2), "rotate_generateBoxSources: coor2(2) should be 4")
    call expect_eq_int(test_err, 5, this%boxSrc%vols(1)%coor2(3), "rotate_generateBoxSources: coor2(3) should be 5")
    
    deallocate(this%boxSrc%vols)
    deallocate(this%boxSrc)
    
    mpidir = 1
    allocate(this%boxSrc)
    this%boxSrc%nvols = 1
    allocate(this%boxSrc%vols(1))
    
    this%boxSrc%vols(1)%coor1(1) = 1
    this%boxSrc%vols(1)%coor1(2) = 2
    this%boxSrc%vols(1)%coor1(3) = 3
    this%boxSrc%vols(1)%coor2(1) = 4
    this%boxSrc%vols(1)%coor2(2) = 5
    this%boxSrc%vols(1)%coor2(3) = 6
    
    call rotate_generateBoxSources(this, mpidir)
    
    call expect_eq_int(test_err, 2, this%boxSrc%vols(1)%coor1(1), "rotate_generateBoxSources: coor1(1) should be 2")
    call expect_eq_int(test_err, 3, this%boxSrc%vols(1)%coor1(2), "rotate_generateBoxSources: coor1(2) should be 3")
    call expect_eq_int(test_err, 1, this%boxSrc%vols(1)%coor1(3), "rotate_generateBoxSources: coor1(3) should be 1")
    
    call expect_eq_int(test_err, 5, this%boxSrc%vols(1)%coor2(1), "rotate_generateBoxSources: coor2(1) should be 5")
    call expect_eq_int(test_err, 6, this%boxSrc%vols(1)%coor2(2), "rotate_generateBoxSources: coor2(2) should be 6")
    call expect_eq_int(test_err, 4, this%boxSrc%vols(1)%coor2(3), "rotate_generateBoxSources: coor2(3) should be 4")
    
    deallocate(this%boxSrc%vols)
    deallocate(this%boxSrc)
    
    err = test_err
end function