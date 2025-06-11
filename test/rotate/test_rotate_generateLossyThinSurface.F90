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
    
    mpidir = 2
    call setup_lossy_thin_surface_test(this)
    call rotate_generateLossyThinSurface(this, mpidir)
    call verify_lossy_thin_surface_rotation(test_err, this, mpidir)
    call cleanup_lossy_thin_surface_test(this)
    mpidir = 1
    call setup_lossy_thin_surface_test(this)
    call rotate_generateLossyThinSurface(this, mpidir)
    call verify_lossy_thin_surface_rotation(test_err, this, mpidir)
    call cleanup_lossy_thin_surface_test(this)
    err = test_err
end function test_rotate_generate_lossy_thin_surface

subroutine setup_lossy_thin_surface_test(this)
    type(Parseador), intent(inout) :: this
    integer(kind=4) :: n_lts = 1, n_ltsc = 2
    allocate(this%LossyThinSurfs)
    allocate(this%LossyThinSurfs%cs(n_lts))
    this%LossyThinSurfs%length = n_lts
    this%LossyThinSurfs%length_max = n_lts
    allocate(this%LossyThinSurfs%cs(1)%c(n_ltsc))
    this%LossyThinSurfs%cs(1)%nc = n_ltsc    
    call init_lossy_thin_surface_data(this, 1, 1, 1, 2, 3, 4, 5, 6, 1, 2 , 3, 1,"tag1")
    call init_lossy_thin_surface_data(this, 1, 2, 7, 8, 9, 10, 11, 12, 4, 5 ,6, 2,"tag2")
end subroutine setup_lossy_thin_surface_test

subroutine init_lossy_thin_surface_data(this, ltsidx, ltscidx, x1, y1, z1, x2, y2, z2, xt, yt, zt, or, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: ltsidx, ltscidx, x1, y1, z1, x2, y2, z2, xt, yt, zt, or
    character(len=*), intent(in) :: tag
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Xi = X1
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Xe = X2
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Yi = Y1
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Ye = Y2
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Zi = Z1
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Ze = Z2
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Xtrancos = xt
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Ytrancos = yt
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Ztrancos = zt
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%Or = or
    this%LossyThinSurfs%cs(ltsidx)%c(ltscidx)%tag  = tag
end subroutine init_lossy_thin_surface_data

subroutine verify_lossy_thin_surface_rotation(test_err, this, mpidir)
    integer, intent(inout) :: test_err, mpidir
    type(Parseador), intent(in) :: this

    if (mpidir==2) then
        call expect_eq_int(test_err, 3, this%LossyThinSurfs%cs(1)%c(1)%Xi, "rotate_generateLossyThinSurface: Xi(1) mpidir2 error")
        call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(1)%Yi, "rotate_generateLossyThinSurface: Yi(1) mpidir2 error")
        call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Zi, "rotate_generateLossyThinSurface: Zi(1) mpidir2 error")
        call expect_eq_int(test_err, 6, this%LossyThinSurfs%cs(1)%c(1)%Xe, "rotate_generateLossyThinSurface: Xe(1) mpidir2 error")
        call expect_eq_int(test_err, 4, this%LossyThinSurfs%cs(1)%c(1)%Ye, "rotate_generateLossyThinSurface: Ye(1) mpidir2 error")
        call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(1)%c(1)%Ze, "rotate_generateLossyThinSurface: Ze(1) mpidir2 error")
        call expect_eq_int(test_err, 3, this%LossyThinSurfs%cs(1)%c(1)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(1) mpidir2 error")
        call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(1)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(1) mpidir2 error")
        call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(1) mpidir2 error")
        call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Or, "rotate_generateLossyThinSurface: Or(1) mpidir2 error")

        call expect_eq_int(test_err, 9, this%LossyThinSurfs%cs(1)%c(2)%Xi, "rotate_generateLossyThinSurface: Xi(2) mpidir2 error")
        call expect_eq_int(test_err, 7, this%LossyThinSurfs%cs(1)%c(2)%Yi, "rotate_generateLossyThinSurface: Yi(2) mpidir2 error")
        call expect_eq_int(test_err, 8, this%LossyThinSurfs%cs(1)%c(2)%Zi, "rotate_generateLossyThinSurface: Zi(2) mpidir2 error")
        call expect_eq_int(test_err, 12, this%LossyThinSurfs%cs(1)%c(2)%Xe, "rotate_generateLossyThinSurface: Xe(2) mpidir2 error")
        call expect_eq_int(test_err, 10, this%LossyThinSurfs%cs(1)%c(2)%Ye, "rotate_generateLossyThinSurface: Ye(2) mpidir2 error")
        call expect_eq_int(test_err, 11, this%LossyThinSurfs%cs(1)%c(2)%Ze, "rotate_generateLossyThinSurface: Ze(2) mpidir2 error")
        call expect_eq_int(test_err, 6, this%LossyThinSurfs%cs(1)%c(2)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(2) mpidir2 error")
        call expect_eq_int(test_err, 4, this%LossyThinSurfs%cs(1)%c(2)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(2) mpidir2 error")
        call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(1)%c(2)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(2) mpidir2 error")
        call expect_eq_int(test_err, 3, this%LossyThinSurfs%cs(1)%c(2)%Or, "rotate_generateLossyThinSurface: Or(2) mpidir2 error")

    else if (mpidir==1) then
        call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Xi, "rotate_generateLossyThinSurface: Xi(1) mpidir1 error")
        call expect_eq_int(test_err, 3, this%LossyThinSurfs%cs(1)%c(1)%Yi, "rotate_generateLossyThinSurface: Yi(1) mpidir1 error")
        call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(1)%Zi, "rotate_generateLossyThinSurface: Zi(1) mpidir1 error")
        call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(1)%c(1)%Xe, "rotate_generateLossyThinSurface: Xe(1) mpidir1 error")
        call expect_eq_int(test_err, 6, this%LossyThinSurfs%cs(1)%c(1)%Ye, "rotate_generateLossyThinSurface: Ye(1) mpidir1 error")
        call expect_eq_int(test_err, 4, this%LossyThinSurfs%cs(1)%c(1)%Ze, "rotate_generateLossyThinSurface: Ze(1) mpidir1 error")
        call expect_eq_int(test_err, 2, this%LossyThinSurfs%cs(1)%c(1)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(1) mpidir1 error")
        call expect_eq_int(test_err, 3, this%LossyThinSurfs%cs(1)%c(1)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(1) mpidir1 error")
        call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(1)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(1) mpidir1 error")
        call expect_eq_int(test_err, 3, this%LossyThinSurfs%cs(1)%c(1)%Or, "rotate_generateLossyThinSurface: Or(1) mpidir1 error")

        call expect_eq_int(test_err, 8, this%LossyThinSurfs%cs(1)%c(2)%Xi, "rotate_generateLossyThinSurface: Xi(2) mpidir1 error")
        call expect_eq_int(test_err, 9, this%LossyThinSurfs%cs(1)%c(2)%Yi, "rotate_generateLossyThinSurface: Yi(2) mpidir1 error")
        call expect_eq_int(test_err, 7, this%LossyThinSurfs%cs(1)%c(2)%Zi, "rotate_generateLossyThinSurface: Zi(2) mpidir1 error")
        call expect_eq_int(test_err, 11, this%LossyThinSurfs%cs(1)%c(2)%Xe, "rotate_generateLossyThinSurface: Xe(2) mpidir1 error")
        call expect_eq_int(test_err, 12, this%LossyThinSurfs%cs(1)%c(2)%Ye, "rotate_generateLossyThinSurface: Ye(2) mpidir1 error")
        call expect_eq_int(test_err, 10, this%LossyThinSurfs%cs(1)%c(2)%Ze, "rotate_generateLossyThinSurface: Ze(2) mpidir1 error")
        call expect_eq_int(test_err, 5, this%LossyThinSurfs%cs(1)%c(2)%Xtrancos, "rotate_generateLossyThinSurface: Xtrancos(2) mpidir1 error")
        call expect_eq_int(test_err, 6, this%LossyThinSurfs%cs(1)%c(2)%Ytrancos, "rotate_generateLossyThinSurface: Ytrancos(2) mpidir1 error")
        call expect_eq_int(test_err, 4, this%LossyThinSurfs%cs(1)%c(2)%Ztrancos, "rotate_generateLossyThinSurface: Ztrancos(2) mpidir1 error")
        call expect_eq_int(test_err, 1, this%LossyThinSurfs%cs(1)%c(2)%Or, "rotate_generateLossyThinSurface: Or(2) mpidir1 error")
    end if
end subroutine verify_lossy_thin_surface_rotation

subroutine cleanup_lossy_thin_surface_test(this)
    type(Parseador), intent(inout) :: this

    deallocate(this%LossyThinSurfs%cs)
    deallocate(this%LossyThinSurfs)
end subroutine cleanup_lossy_thin_surface_test

end module test_rotate_generateLossyThinSurface_m 