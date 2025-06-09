integer function test_rotate_generate_plane_waves() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools    
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    real(kind=RKIND) :: theta, phi, alpha, beta
    integer :: test_err = 0  ! Error counter for expect_eq calls
    
    ! Test case 1: X->Y->Z->X rotation (mpidir=2)
    mpidir = 2
    allocate(this%plnSrc)
    this%plnSrc%nc = 1
    allocate(this%plnSrc%collection(1))
    
    ! Initialize test data
    this%plnSrc%collection(1)%coor1(1) = 1
    this%plnSrc%collection(1)%coor1(2) = 2
    this%plnSrc%collection(1)%coor1(3) = 3
    this%plnSrc%collection(1)%coor2(1) = 4
    this%plnSrc%collection(1)%coor2(2) = 5
    this%plnSrc%collection(1)%coor2(3) = 6
    
    ! Set angles in radians
    theta = 0.5_RKIND
    phi = 0.3_RKIND
    alpha = 0.7_RKIND
    beta = 0.2_RKIND
    
    this%plnSrc%collection(1)%theta = theta
    this%plnSrc%collection(1)%phi = phi
    this%plnSrc%collection(1)%alpha = alpha
    this%plnSrc%collection(1)%beta = beta
    
    ! Call the routine
    call rotate_generatePlaneWaves(this, mpidir)
    
    ! Verify results
    call expect_eq_int(test_err, 3, this%plnSrc%collection(1)%coor1(1), "rotate_generatePlaneWaves: coor1(1) should be 3")
    call expect_eq_int(test_err, 1, this%plnSrc%collection(1)%coor1(2), "rotate_generatePlaneWaves: coor1(2) should be 1")
    call expect_eq_int(test_err, 2, this%plnSrc%collection(1)%coor1(3), "rotate_generatePlaneWaves: coor1(3) should be 2")
    
    call expect_eq_int(test_err, 6, this%plnSrc%collection(1)%coor2(1), "rotate_generatePlaneWaves: coor2(1) should be 6")
    call expect_eq_int(test_err, 4, this%plnSrc%collection(1)%coor2(2), "rotate_generatePlaneWaves: coor2(2) should be 4")
    call expect_eq_int(test_err, 5, this%plnSrc%collection(1)%coor2(3), "rotate_generatePlaneWaves: coor2(3) should be 5")
    
    ! Verify angles (using approximate equality due to floating point)
    call expect_eq_real(test_err, &
        atan2(sqrt(cos(theta)**2.0_RKIND + cos(phi)**2*sin(theta)**2), sin(phi)*sin(theta)), &
        this%plnSrc%collection(1)%theta, &
        "rotate_generatePlaneWaves: theta rotation incorrect")
        
    call expect_eq_real(test_err, &
        atan2(cos(phi)*sin(theta), cos(theta)), &
        this%plnSrc%collection(1)%phi, &
        "rotate_generatePlaneWaves: phi rotation incorrect")
        
    call expect_eq_real(test_err, &
        atan2(sqrt(cos(alpha)**2.0_RKIND + cos(beta)**2*sin(alpha)**2), sin(beta)*sin(alpha)), &
        this%plnSrc%collection(1)%alpha, &
        "rotate_generatePlaneWaves: alpha rotation incorrect")
        
    call expect_eq_real(test_err, &
        atan2(cos(beta)*sin(alpha), cos(alpha)), &
        this%plnSrc%collection(1)%beta, &
        "rotate_generatePlaneWaves: beta rotation incorrect")
    
    deallocate(this%plnSrc%collection)
    deallocate(this%plnSrc)
    
    ! Test case 2: X->Z->Y->X rotation (mpidir=1)
    mpidir = 1
    allocate(this%plnSrc)
    this%plnSrc%nc = 1
    allocate(this%plnSrc%collection(1))
    
    ! Initialize test data
    this%plnSrc%collection(1)%coor1(1) = 1
    this%plnSrc%collection(1)%coor1(2) = 2
    this%plnSrc%collection(1)%coor1(3) = 3
    this%plnSrc%collection(1)%coor2(1) = 4
    this%plnSrc%collection(1)%coor2(2) = 5
    this%plnSrc%collection(1)%coor2(3) = 6
    
    this%plnSrc%collection(1)%theta = theta
    this%plnSrc%collection(1)%phi = phi
    this%plnSrc%collection(1)%alpha = alpha
    this%plnSrc%collection(1)%beta = beta
    
    ! Call the routine
    call rotate_generatePlaneWaves(this, mpidir)
    
    ! Verify results
    call expect_eq_int(test_err, 2, this%plnSrc%collection(1)%coor1(1), "rotate_generatePlaneWaves: coor1(1) should be 2")
    call expect_eq_int(test_err, 3, this%plnSrc%collection(1)%coor1(2), "rotate_generatePlaneWaves: coor1(2) should be 3")
    call expect_eq_int(test_err, 1, this%plnSrc%collection(1)%coor1(3), "rotate_generatePlaneWaves: coor1(3) should be 1")
    
    call expect_eq_int(test_err, 5, this%plnSrc%collection(1)%coor2(1), "rotate_generatePlaneWaves: coor2(1) should be 5")
    call expect_eq_int(test_err, 6, this%plnSrc%collection(1)%coor2(2), "rotate_generatePlaneWaves: coor2(2) should be 6")
    call expect_eq_int(test_err, 4, this%plnSrc%collection(1)%coor2(3), "rotate_generatePlaneWaves: coor2(3) should be 4")
    
    ! Verify angles
    call expect_eq_real(test_err, &
        atan2(sqrt(cos(theta)**2.0_RKIND + sin(phi)**2*sin(theta)**2), cos(phi)*sin(theta)), &
        this%plnSrc%collection(1)%theta, &
        "rotate_generatePlaneWaves: theta rotation incorrect")
        
    call expect_eq_real(test_err, &
        atan2(cos(theta), sin(phi)*sin(theta)), &
        this%plnSrc%collection(1)%phi, &
        "rotate_generatePlaneWaves: phi rotation incorrect")
        
    call expect_eq_real(test_err, &
        atan2(sqrt(cos(alpha)**2.0_RKIND + sin(beta)**2*sin(alpha)**2), cos(beta)*sin(alpha)), &
        this%plnSrc%collection(1)%alpha, &
        "rotate_generatePlaneWaves: alpha rotation incorrect")
        
    call expect_eq_real(test_err, &
        atan2(cos(alpha), sin(beta)*sin(alpha)), &
        this%plnSrc%collection(1)%beta, &
        "rotate_generatePlaneWaves: beta rotation incorrect")
    
    deallocate(this%plnSrc%collection)
    deallocate(this%plnSrc)
    
    err = test_err  ! Set the function result to the accumulated error count
end function