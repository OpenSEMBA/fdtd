module test_rotate_generateLossyThinSurface_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_lossy_thin_surface() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    ! Test case 1: Rotate lossy thin surface in mpidir=2 direction
    mpidir = 2
    call setup_lossy_thin_surface_test(this, 2)
    
    ! Initialize test data
    call init_lossy_thin_surface_data(this, 1, 1, 2, 3, 4, 5, 3, 0, 0, 1, 1.0d0)
    call init_lossy_thin_surface_data(this, 2, 5, 6, 7, 5, 8, 10, 1, 0, 0, 2.0d0)
    
    ! Call rotation and verify results
    call rotate_generateLossyThinSurface(this, mpidir)
    call verify_lossy_thin_surface_rotation_y(test_err, this)
    
    ! Clean up
    call cleanup_lossy_thin_surface_test(this)
    
    ! Test case 2: Rotate lossy thin surface in mpidir=1 direction
    mpidir = 1
    call setup_lossy_thin_surface_test(this, 2)
    
    ! Set up test data (same as before)
    call init_lossy_thin_surface_data(this, 1, 1, 2, 3, 4, 5, 3, 0, 0, 1, 1.0d0)
    call init_lossy_thin_surface_data(this, 2, 5, 6, 7, 5, 8, 10, 1, 0, 0, 2.0d0)
    
    ! Call rotation and verify results
    call rotate_generateLossyThinSurface(this, mpidir)
    call verify_lossy_thin_surface_rotation_x(test_err, this)
    
    ! Clean up
    call cleanup_lossy_thin_surface_test(this)
    
    ! Test case 3: Verify behavior with no surfaces
    mpidir = 2
    call setup_lossy_thin_surface_test(this, 0)
    call rotate_generateLossyThinSurface(this, mpidir)
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%length, "rotate_generateLossyThinSurface: length should remain 0")
    call cleanup_lossy_thin_surface_test(this)
    
    err = test_err
end function test_rotate_generate_lossy_thin_surface

subroutine setup_lossy_thin_surface_test(this, n_surfaces)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_surfaces
    integer :: i
    
    allocate(this%LossyThinSurfs)
    this%LossyThinSurfs%length = n_surfaces
    this%LossyThinSurfs%length_max = n_surfaces
    if (n_surfaces > 0) then
        allocate(this%LossyThinSurfs%cs(n_surfaces))
        ! Initialize each surface
        do i = 1, n_surfaces
            this%LossyThinSurfs%cs(i)%nc = 1
            allocate(this%LossyThinSurfs%cs(i)%c(1))
            allocate(this%LossyThinSurfs%cs(i)%sigma(1))
            allocate(this%LossyThinSurfs%cs(i)%eps(1))
            allocate(this%LossyThinSurfs%cs(i)%mu(1))
            allocate(this%LossyThinSurfs%cs(i)%sigmam(1))
            allocate(this%LossyThinSurfs%cs(i)%thk(1))
            allocate(this%LossyThinSurfs%cs(i)%sigma_devia(1))
            allocate(this%LossyThinSurfs%cs(i)%eps_devia(1))
            allocate(this%LossyThinSurfs%cs(i)%mu_devia(1))
            allocate(this%LossyThinSurfs%cs(i)%sigmam_devia(1))
            allocate(this%LossyThinSurfs%cs(i)%thk_devia(1))
            ! Initialize default values
            this%LossyThinSurfs%cs(i)%eps(1) = 1.0_RKIND
            this%LossyThinSurfs%cs(i)%mu(1) = 1.0_RKIND
            this%LossyThinSurfs%cs(i)%sigmam(1) = 0.0_RKIND
            this%LossyThinSurfs%cs(i)%thk(1) = 0.0_RKIND
            this%LossyThinSurfs%cs(i)%eps_devia(1) = 0.0_RKIND
            this%LossyThinSurfs%cs(i)%mu_devia(1) = 0.0_RKIND
            this%LossyThinSurfs%cs(i)%sigmam_devia(1) = 0.0_RKIND
            this%LossyThinSurfs%cs(i)%thk_devia(1) = 0.0_RKIND
            this%LossyThinSurfs%cs(i)%numcapas = 1
        end do
    end if
end subroutine setup_lossy_thin_surface_test

subroutine init_lossy_thin_surface_data(this, idx, x1, y1, z1, x2, y2, z2, dirx, diry, dirz, sigma_val)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, x1, y1, z1, x2, y2, z2, dirx, diry, dirz
    real(8), intent(in) :: sigma_val
    
    ! Initialize coordinates and properties
    this%LossyThinSurfs%cs(idx)%c(1)%Xi = x1
    this%LossyThinSurfs%cs(idx)%c(1)%Yi = y1
    this%LossyThinSurfs%cs(idx)%c(1)%Zi = z1
    this%LossyThinSurfs%cs(idx)%c(1)%Xe = x2
    this%LossyThinSurfs%cs(idx)%c(1)%Ye = y2
    this%LossyThinSurfs%cs(idx)%c(1)%Ze = z2
    this%LossyThinSurfs%cs(idx)%c(1)%Xtrancos = dirx
    this%LossyThinSurfs%cs(idx)%c(1)%Ytrancos = diry
    this%LossyThinSurfs%cs(idx)%c(1)%Ztrancos = dirz
    this%LossyThinSurfs%cs(idx)%c(1)%Or = 0
    this%LossyThinSurfs%cs(idx)%c(1)%tag = "test_surface"
    this%LossyThinSurfs%cs(idx)%sigma(1) = sigma_val
end subroutine init_lossy_thin_surface_data

subroutine verify_lossy_thin_surface_rotation_y(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify Surface 1 (rotated around y-axis)
    call expect_eq_int(test_err, -3, this%LossyThinSurfs%cs(1)%c(1)%Xi, "rotate_generateLossyThinSurface: Xi(1) should be rotated")
    call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Yi, "rotate_generateLossyThinSurface: Yi(1) should remain unchanged")
    call expect_eq_int(test_err, -1, this%LossyThinSurfs%cs(1)%c(1)%Zi, "rotate_generateLossyThinSurface: Zi(1) should be rotated")
    call expect_eq_int(test_err, -3, this%LossyThinSurfs%cs(1)%c(1)%Xe, "rotate_generateLossyThinSurface: Xe(1) should be rotated")
    call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(1)%c(1)%Ye, "rotate_generateLossyThinSurface: Ye(1) should remain unchanged")
    call expect_eq_int(test_err, -4, this%LossyThinSurfs%cs(1)%c(1)%Ze, "rotate_generateLossyThinSurface: Ze(1) should be rotated")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(1)%c(1)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(1) should be rotated")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(1)%c(1)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(1) should remain unchanged")
    call expect_eq_int(test_err, -1, this%LossyThinSurfs%cs(1)%c(1)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(1) should be rotated")
    call expect_eq_real(test_err, 1.0_RKIND, this%LossyThinSurfs%cs(1)%sigma(1), "rotate_generateLossyThinSurface: sigma(1) should remain unchanged")
    
    ! Verify Surface 2 (rotated around y-axis)
    call expect_eq_int(test_err, -7, this%LossyThinSurfs%cs(2)%c(1)%Xi, "rotate_generateLossyThinSurface: Xi(2) should be rotated")
    call expect_eq_int(test_err, 6, this%LossyThinSurfs%cs(2)%c(1)%Yi, "rotate_generateLossyThinSurface: Yi(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%LossyThinSurfs%cs(2)%c(1)%Zi, "rotate_generateLossyThinSurface: Zi(2) should be rotated")
    call expect_eq_int(test_err, -10, this%LossyThinSurfs%cs(2)%c(1)%Xe, "rotate_generateLossyThinSurface: Xe(2) should be rotated")
    call expect_eq_int(test_err, 8, this%LossyThinSurfs%cs(2)%c(1)%Ye, "rotate_generateLossyThinSurface: Ye(2) should remain unchanged")
    call expect_eq_int(test_err, -5, this%LossyThinSurfs%cs(2)%c(1)%Ze, "rotate_generateLossyThinSurface: Ze(2) should be rotated")
    call expect_eq_int(test_err, -1, this%LossyThinSurfs%cs(2)%c(1)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(2) should be rotated")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(2)%c(1)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(2)%c(1)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(2) should be rotated")
    call expect_eq_real(test_err, 2.0_RKIND, this%LossyThinSurfs%cs(2)%sigma(1), "rotate_generateLossyThinSurface: sigma(2) should remain unchanged")
end subroutine verify_lossy_thin_surface_rotation_y

subroutine verify_lossy_thin_surface_rotation_x(test_err, this)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    
    ! Verify Surface 1 (rotated around x-axis)
    call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(1)%Xi, "rotate_generateLossyThinSurface: Xi(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%LossyThinSurfs%cs(1)%c(1)%Yi, "rotate_generateLossyThinSurface: Yi(1) should be rotated")
    call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Zi, "rotate_generateLossyThinSurface: Zi(1) should be rotated")
    call expect_eq_int(test_err, 4, this%LossyThinSurfs%cs(1)%c(1)%Xe, "rotate_generateLossyThinSurface: Xe(1) should remain unchanged")
    call expect_eq_int(test_err, -3, this%LossyThinSurfs%cs(1)%c(1)%Ye, "rotate_generateLossyThinSurface: Ye(1) should be rotated")
    call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(1)%c(1)%Ze, "rotate_generateLossyThinSurface: Ze(1) should be rotated")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(1)%c(1)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(1)%c(1)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(1) should be rotated")
    call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(1)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(1) should be rotated")
    call expect_eq_real(test_err, 1.0_RKIND, this%LossyThinSurfs%cs(1)%sigma(1), "rotate_generateLossyThinSurface: sigma(1) should remain unchanged")
    
    ! Verify Surface 2 (rotated around x-axis)
    call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(2)%c(1)%Xi, "rotate_generateLossyThinSurface: Xi(2) should remain unchanged")
    call expect_eq_int(test_err, -7, this%LossyThinSurfs%cs(2)%c(1)%Yi, "rotate_generateLossyThinSurface: Yi(2) should be rotated")
    call expect_eq_int(test_err, 6, this%LossyThinSurfs%cs(2)%c(1)%Zi, "rotate_generateLossyThinSurface: Zi(2) should be rotated")
    call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(2)%c(1)%Xe, "rotate_generateLossyThinSurface: Xe(2) should remain unchanged")
    call expect_eq_int(test_err, -10, this%LossyThinSurfs%cs(2)%c(1)%Ye, "rotate_generateLossyThinSurface: Ye(2) should be rotated")
    call expect_eq_int(test_err, 8, this%LossyThinSurfs%cs(2)%c(1)%Ze, "rotate_generateLossyThinSurface: Ze(2) should be rotated")
    call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(2)%c(1)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(2)%c(1)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(2) should be rotated")
    call expect_eq_int(test_err, 0, this%LossyThinSurfs%cs(2)%c(1)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(2) should be rotated")
    call expect_eq_real(test_err, 2.0_RKIND, this%LossyThinSurfs%cs(2)%sigma(1), "rotate_generateLossyThinSurface: sigma(2) should remain unchanged")
end subroutine verify_lossy_thin_surface_rotation_x

subroutine cleanup_lossy_thin_surface_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%LossyThinSurfs%length > 0) then
        do i = 1, this%LossyThinSurfs%length
            if (this%LossyThinSurfs%cs(i)%nc > 0) then
                deallocate(this%LossyThinSurfs%cs(i)%c)
                deallocate(this%LossyThinSurfs%cs(i)%sigma)
                deallocate(this%LossyThinSurfs%cs(i)%eps)
                deallocate(this%LossyThinSurfs%cs(i)%mu)
                deallocate(this%LossyThinSurfs%cs(i)%sigmam)
                deallocate(this%LossyThinSurfs%cs(i)%thk)
                deallocate(this%LossyThinSurfs%cs(i)%sigma_devia)
                deallocate(this%LossyThinSurfs%cs(i)%eps_devia)
                deallocate(this%LossyThinSurfs%cs(i)%mu_devia)
                deallocate(this%LossyThinSurfs%cs(i)%sigmam_devia)
                deallocate(this%LossyThinSurfs%cs(i)%thk_devia)
            end if
        end do
        deallocate(this%LossyThinSurfs%cs)
    end if
    deallocate(this%LossyThinSurfs)
end subroutine cleanup_lossy_thin_surface_test

end module test_rotate_generateLossyThinSurface_m 