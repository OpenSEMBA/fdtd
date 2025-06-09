integer function test_rotate_generate_fronteras() bind(C) result(err)
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0  ! Error counter for expect_eq calls
    
    ! Test case 1: X->Y->Z->X rotation (mpidir=2)
    mpidir = 2
    allocate(this%front)
    
    ! Initialize test data
    this%front%tipofrontera(1) = 1  ! XL
    this%front%tipofrontera(2) = 2  ! XU
    this%front%tipofrontera(3) = 3  ! YL
    this%front%tipofrontera(4) = 4  ! YU
    this%front%tipofrontera(5) = 5  ! ZL
    this%front%tipofrontera(6) = 6  ! ZU
    
    ! Initialize PML properties
    this%front%propiedadesPML(1)%orden = 1.0_RKIND    ! XL
    this%front%propiedadesPML(2)%orden = 2.0_RKIND    ! XU
    this%front%propiedadesPML(3)%orden = 3.0_RKIND    ! YL
    this%front%propiedadesPML(4)%orden = 4.0_RKIND    ! YU
    this%front%propiedadesPML(5)%orden = 5.0_RKIND    ! ZL
    this%front%propiedadesPML(6)%orden = 6.0_RKIND    ! ZU
    
    this%front%propiedadesPML(1)%refl = 0.1_RKIND  ! XL
    this%front%propiedadesPML(2)%refl = 0.2_RKIND  ! XU
    this%front%propiedadesPML(3)%refl = 0.3_RKIND  ! YL
    this%front%propiedadesPML(4)%refl = 0.4_RKIND  ! YU
    this%front%propiedadesPML(5)%refl = 0.5_RKIND  ! ZL
    this%front%propiedadesPML(6)%refl = 0.6_RKIND  ! ZU
    
    this%front%propiedadesPML(1)%numCapas = 10  ! XL
    this%front%propiedadesPML(2)%numCapas = 20  ! XU
    this%front%propiedadesPML(3)%numCapas = 30  ! YL
    this%front%propiedadesPML(4)%numCapas = 40  ! YU
    this%front%propiedadesPML(5)%numCapas = 50  ! ZL
    this%front%propiedadesPML(6)%numCapas = 60  ! ZU
    
    ! Call the routine
    call rotate_generateFronteras(this, mpidir)
    
    ! Verify results
    call expect_eq_int(test_err, this%front%tipofrontera(1), 5, "rotate_generateFronteras: XL should be 5")
    call expect_eq_int(test_err, this%front%tipofrontera(2), 6, "rotate_generateFronteras: XU should be 6")
    call expect_eq_int(test_err, this%front%tipofrontera(3), 1, "rotate_generateFronteras: YL should be 1")
    call expect_eq_int(test_err, this%front%tipofrontera(4), 2, "rotate_generateFronteras: YU should be 2")
    call expect_eq_int(test_err, this%front%tipofrontera(5), 3, "rotate_generateFronteras: ZL should be 3")
    call expect_eq_int(test_err, this%front%tipofrontera(6), 4, "rotate_generateFronteras: ZU should be 4")
    
    ! Verify PML properties
    call expect_eq_real(test_err, this%front%propiedadesPML(1)%orden, 5.0_RKIND, "rotate_generateFronteras: XL orden should be 5")
    call expect_eq_real(test_err, this%front%propiedadesPML(2)%orden, 6.0_RKIND, "rotate_generateFronteras: XU orden should be 6")
    call expect_eq_real(test_err, this%front%propiedadesPML(3)%orden, 1.0_RKIND, "rotate_generateFronteras: YL orden should be 1")
    call expect_eq_real(test_err, this%front%propiedadesPML(4)%orden, 2.0_RKIND, "rotate_generateFronteras: YU orden should be 2")
    call expect_eq_real(test_err, this%front%propiedadesPML(5)%orden, 3.0_RKIND, "rotate_generateFronteras: ZL orden should be 3")
    call expect_eq_real(test_err, this%front%propiedadesPML(6)%orden, 4.0_RKIND, "rotate_generateFronteras: ZU orden should be 4")
    
    call expect_eq_real(test_err, this%front%propiedadesPML(1)%refl, 0.5_RKIND, "rotate_generateFronteras: XL refl should be 0.5")
    call expect_eq_real(test_err, this%front%propiedadesPML(2)%refl, 0.6_RKIND, "rotate_generateFronteras: XU refl should be 0.6")
    call expect_eq_real(test_err, this%front%propiedadesPML(3)%refl, 0.1_RKIND, "rotate_generateFronteras: YL refl should be 0.1")
    call expect_eq_real(test_err, this%front%propiedadesPML(4)%refl, 0.2_RKIND, "rotate_generateFronteras: YU refl should be 0.2")
    call expect_eq_real(test_err, this%front%propiedadesPML(5)%refl, 0.3_RKIND, "rotate_generateFronteras: ZL refl should be 0.3")
    call expect_eq_real(test_err, this%front%propiedadesPML(6)%refl, 0.4_RKIND, "rotate_generateFronteras: ZU refl should be 0.4")
    
    call expect_eq_int(test_err, this%front%propiedadesPML(1)%numCapas, 50, "rotate_generateFronteras: XL numCapas should be 50")
    call expect_eq_int(test_err, this%front%propiedadesPML(2)%numCapas, 60, "rotate_generateFronteras: XU numCapas should be 60")
    call expect_eq_int(test_err, this%front%propiedadesPML(3)%numCapas, 10, "rotate_generateFronteras: YL numCapas should be 10")
    call expect_eq_int(test_err, this%front%propiedadesPML(4)%numCapas, 20, "rotate_generateFronteras: YU numCapas should be 20")
    call expect_eq_int(test_err, this%front%propiedadesPML(5)%numCapas, 30, "rotate_generateFronteras: ZL numCapas should be 30")
    call expect_eq_int(test_err, this%front%propiedadesPML(6)%numCapas, 40, "rotate_generateFronteras: ZU numCapas should be 40")
    
    deallocate(this%front)
    
    ! Test case 2: X->Z->Y->X rotation (mpidir=1)
    mpidir = 1
    allocate(this%front)
    
    ! Initialize test data (same as before)
    this%front%tipofrontera(1) = 1  ! XL
    this%front%tipofrontera(2) = 2  ! XU
    this%front%tipofrontera(3) = 3  ! YL
    this%front%tipofrontera(4) = 4  ! YU
    this%front%tipofrontera(5) = 5  ! ZL
    this%front%tipofrontera(6) = 6  ! ZU
    
    this%front%propiedadesPML(1)%orden = 1.0_RKIND    ! XL
    this%front%propiedadesPML(2)%orden = 2.0_RKIND    ! XU
    this%front%propiedadesPML(3)%orden = 3.0_RKIND    ! YL
    this%front%propiedadesPML(4)%orden = 4.0_RKIND    ! YU
    this%front%propiedadesPML(5)%orden = 5.0_RKIND    ! ZL
    this%front%propiedadesPML(6)%orden = 6.0_RKIND    ! ZU
    
    this%front%propiedadesPML(1)%refl = 0.1_RKIND  ! XL
    this%front%propiedadesPML(2)%refl = 0.2_RKIND  ! XU
    this%front%propiedadesPML(3)%refl = 0.3_RKIND  ! YL
    this%front%propiedadesPML(4)%refl = 0.4_RKIND  ! YU
    this%front%propiedadesPML(5)%refl = 0.5_RKIND  ! ZL
    this%front%propiedadesPML(6)%refl = 0.6_RKIND  ! ZU
    
    this%front%propiedadesPML(1)%numCapas = 10  ! XL
    this%front%propiedadesPML(2)%numCapas = 20  ! XU
    this%front%propiedadesPML(3)%numCapas = 30  ! YL
    this%front%propiedadesPML(4)%numCapas = 40  ! YU
    this%front%propiedadesPML(5)%numCapas = 50  ! ZL
    this%front%propiedadesPML(6)%numCapas = 60  ! ZU
    
    ! Call the routine
    call rotate_generateFronteras(this, mpidir)
    
    ! Verify results
    call expect_eq_int(test_err, this%front%tipofrontera(1), 3, "rotate_generateFronteras: XL should be 3")
    call expect_eq_int(test_err, this%front%tipofrontera(2), 4, "rotate_generateFronteras: XU should be 4")
    call expect_eq_int(test_err, this%front%tipofrontera(3), 5, "rotate_generateFronteras: YL should be 5")
    call expect_eq_int(test_err, this%front%tipofrontera(4), 6, "rotate_generateFronteras: YU should be 6")
    call expect_eq_int(test_err, this%front%tipofrontera(5), 1, "rotate_generateFronteras: ZL should be 1")
    call expect_eq_int(test_err, this%front%tipofrontera(6), 2, "rotate_generateFronteras: ZU should be 2")
    
    ! Verify PML properties
    call expect_eq_real(test_err, this%front%propiedadesPML(1)%orden, 3.0_RKIND, "rotate_generateFronteras: XL orden should be 3")
    call expect_eq_real(test_err, this%front%propiedadesPML(2)%orden, 4.0_RKIND, "rotate_generateFronteras: XU orden should be 4")
    call expect_eq_real(test_err, this%front%propiedadesPML(3)%orden, 5.0_RKIND, "rotate_generateFronteras: YL orden should be 5")
    call expect_eq_real(test_err, this%front%propiedadesPML(4)%orden, 6.0_RKIND, "rotate_generateFronteras: YU orden should be 6")
    call expect_eq_real(test_err, this%front%propiedadesPML(5)%orden, 1.0_RKIND, "rotate_generateFronteras: ZL orden should be 1")
    call expect_eq_real(test_err, this%front%propiedadesPML(6)%orden, 2.0_RKIND, "rotate_generateFronteras: ZU orden should be 2")
    
    call expect_eq_real(test_err, this%front%propiedadesPML(1)%refl, 0.3_RKIND, "rotate_generateFronteras: XL refl should be 0.3")
    call expect_eq_real(test_err, this%front%propiedadesPML(2)%refl, 0.4_RKIND, "rotate_generateFronteras: XU refl should be 0.4")
    call expect_eq_real(test_err, this%front%propiedadesPML(3)%refl, 0.5_RKIND, "rotate_generateFronteras: YL refl should be 0.5")
    call expect_eq_real(test_err, this%front%propiedadesPML(4)%refl, 0.6_RKIND, "rotate_generateFronteras: YU refl should be 0.6")
    call expect_eq_real(test_err, this%front%propiedadesPML(5)%refl, 0.1_RKIND, "rotate_generateFronteras: ZL refl should be 0.1")
    call expect_eq_real(test_err, this%front%propiedadesPML(6)%refl, 0.2_RKIND, "rotate_generateFronteras: ZU refl should be 0.2")
    
    call expect_eq_int(test_err, this%front%propiedadesPML(1)%numCapas, 30, "rotate_generateFronteras: XL numCapas should be 30")
    call expect_eq_int(test_err, this%front%propiedadesPML(2)%numCapas, 40, "rotate_generateFronteras: XU numCapas should be 40")
    call expect_eq_int(test_err, this%front%propiedadesPML(3)%numCapas, 50, "rotate_generateFronteras: YL numCapas should be 50")
    call expect_eq_int(test_err, this%front%propiedadesPML(4)%numCapas, 60, "rotate_generateFronteras: YU numCapas should be 60")
    call expect_eq_int(test_err, this%front%propiedadesPML(5)%numCapas, 10, "rotate_generateFronteras: ZL numCapas should be 10")
    call expect_eq_int(test_err, this%front%propiedadesPML(6)%numCapas, 20, "rotate_generateFronteras: ZU numCapas should be 20")
    
    deallocate(this%front)
    
    err = test_err  ! Set the function result to the accumulated error count
end function