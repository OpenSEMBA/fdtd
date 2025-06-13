integer function test_rotate_generate_space_steps() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    type(Desplazamiento), pointer :: old_despl => null()
    type(MatrizMedios), pointer :: old_matriz => null()
    integer :: test_err = 0
    
    mpidir = 2
    allocate(this%despl)
    allocate(this%matriz)
    
    this%matriz%totalX = 10
    this%matriz%totalY = 20
    this%matriz%totalZ = 30
    
    this%despl%nX = 5
    this%despl%nY = 15
    this%despl%nZ = 25
    
    this%despl%mX1 = 1
    this%despl%mY1 = 11
    this%despl%mZ1 = 21
    
    this%despl%mX2 = 6
    this%despl%mY2 = 16
    this%despl%mZ2 = 26
    
    this%despl%originX = 0.0_RKIND
    this%despl%originY = 0.0_RKIND
    this%despl%originZ = 0.0_RKIND
    
    call rotate_generateSpaceSteps(this, mpidir)
    
    call expect_eq_int(test_err, 30, this%matriz%totalX, "rotate_generateSpaceSteps: totalX should be 30")
    call expect_eq_int(test_err, 10, this%matriz%totalY, "rotate_generateSpaceSteps: totalY should be 10")
    call expect_eq_int(test_err, 20, this%matriz%totalZ, "rotate_generateSpaceSteps: totalZ should be 20")
    
    call expect_eq_int(test_err, 25, this%despl%nX, "rotate_generateSpaceSteps: nX should be 25")
    call expect_eq_int(test_err, 5, this%despl%nY, "rotate_generateSpaceSteps: nY should be 5")
    call expect_eq_int(test_err, 15, this%despl%nZ, "rotate_generateSpaceSteps: nZ should be 15")
    
    call expect_eq_int(test_err, 21, this%despl%mX1, "rotate_generateSpaceSteps: mX1 should be 21")
    call expect_eq_int(test_err, 1, this%despl%mY1, "rotate_generateSpaceSteps: mY1 should be 1")
    call expect_eq_int(test_err, 11, this%despl%mZ1, "rotate_generateSpaceSteps: mZ1 should be 11")
    
    call expect_eq_int(test_err, 26, this%despl%mX2, "rotate_generateSpaceSteps: mX2 should be 26")
    call expect_eq_int(test_err, 6, this%despl%mY2, "rotate_generateSpaceSteps: mY2 should be 6")
    call expect_eq_int(test_err, 16, this%despl%mZ2, "rotate_generateSpaceSteps: mZ2 should be 16")
    
    deallocate(this%despl)
    deallocate(this%matriz)
    
    mpidir = 1
    allocate(this%despl)
    allocate(this%matriz)
    
    this%matriz%totalX = 10
    this%matriz%totalY = 20
    this%matriz%totalZ = 30
    
    this%despl%nX = 5
    this%despl%nY = 15
    this%despl%nZ = 25
    
    this%despl%mX1 = 1
    this%despl%mY1 = 11
    this%despl%mZ1 = 21
    
    this%despl%mX2 = 6
    this%despl%mY2 = 16
    this%despl%mZ2 = 26
    
    this%despl%originX = 0.0_RKIND
    this%despl%originY = 0.0_RKIND
    this%despl%originZ = 0.0_RKIND
    
    call rotate_generateSpaceSteps(this, mpidir)
    
    call expect_eq_int(test_err, 20, this%matriz%totalX, "rotate_generateSpaceSteps: totalX should be 20")
    call expect_eq_int(test_err, 30, this%matriz%totalY, "rotate_generateSpaceSteps: totalY should be 30")
    call expect_eq_int(test_err, 10, this%matriz%totalZ, "rotate_generateSpaceSteps: totalZ should be 10")
    
    call expect_eq_int(test_err, 15, this%despl%nX, "rotate_generateSpaceSteps: nX should be 15")
    call expect_eq_int(test_err, 25, this%despl%nY, "rotate_generateSpaceSteps: nY should be 25")
    call expect_eq_int(test_err, 5, this%despl%nZ, "rotate_generateSpaceSteps: nZ should be 5")
    
    call expect_eq_int(test_err, 11, this%despl%mX1, "rotate_generateSpaceSteps: mX1 should be 11")
    call expect_eq_int(test_err, 21, this%despl%mY1, "rotate_generateSpaceSteps: mY1 should be 21")
    call expect_eq_int(test_err, 1, this%despl%mZ1, "rotate_generateSpaceSteps: mZ1 should be 1")
    
    call expect_eq_int(test_err, 16, this%despl%mX2, "rotate_generateSpaceSteps: mX2 should be 16")
    call expect_eq_int(test_err, 26, this%despl%mY2, "rotate_generateSpaceSteps: mY2 should be 26")
    call expect_eq_int(test_err, 6, this%despl%mZ2, "rotate_generateSpaceSteps: mZ2 should be 6")
    
    deallocate(this%despl)
    deallocate(this%matriz)
    
    err = test_err
end function