module test_rotate_generateVolumicProbes_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_volumic_probes() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0

    mpidir = 2
    call setup_volumic_probe_test(this, 2)

    call init_volumic_probe_data(this, 1, 1, 4, 2, 5, 3, 6, 1, 2, 3, 1, "probe1")
    call init_volumic_probe_data(this, 2, 5, 8, 6, 9, 7, 10, 1, 2, 3, 2, "probe2")
    
    call rotate_generateVolumicProbes(this, mpidir)
    call verify_volumic_probe_rotation(test_err, this, mpidir)
    
    call cleanup_volumic_probe_test(this)
    
    mpidir = 1
    call setup_volumic_probe_test(this, 2)
    
    call init_volumic_probe_data(this, 1, 1, 4, 2, 5, 3, 6, 1, 2, 3, 1, "probe1")
    call init_volumic_probe_data(this, 2, 5, 8, 6, 9, 7, 10, 1, 2, 3, 2, "probe2")
    
    call rotate_generateVolumicProbes(this, mpidir)
    call verify_volumic_probe_rotation(test_err, this, mpidir)
    
    call cleanup_volumic_probe_test(this)    
    err = test_err
end function test_rotate_generate_volumic_probes

subroutine setup_volumic_probe_test(this, n_probes)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_probes
    integer :: i
    
    allocate(this%VolPrb)
    this%VolPrb%length = n_probes
    this%VolPrb%length_max = n_probes
    this%VolPrb%len_cor_max = 1
    allocate(this%VolPrb%collection(n_probes))
    if (n_probes > 0) then
        do i = 1, n_probes
            this%VolPrb%collection(i)%tstart = 0.0d0
            this%VolPrb%collection(i)%tstop = 1.0d0
            this%VolPrb%collection(i)%tstep = 0.1d0
            this%VolPrb%collection(i)%fstart = 0.0d0
            this%VolPrb%collection(i)%fstop = 1.0d0
            this%VolPrb%collection(i)%fstep = 0.1d0
            this%VolPrb%collection(i)%outputrequest = "output.txt"
            this%VolPrb%collection(i)%filename = "probe.dat"
            this%VolPrb%collection(i)%len_cor = 1
            this%VolPrb%collection(i)%type2 = 1
            allocate(this%VolPrb%collection(i)%cordinates(1))
        end do
    end if
end subroutine setup_volumic_probe_test

subroutine init_volumic_probe_data(this, idx, xi, xe, yi, ye, zi, ze, xtrancos, ytrancos, ztrancos, or, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, xi, xe, yi, ye, zi, ze, xtrancos, ytrancos, ztrancos, or
    character(len=*), intent(in) :: tag
    this%VolPrb%collection(idx)%len_cor = 1
    this%VolPrb%collection(idx)%cordinates(1)%Xi = xi
    this%VolPrb%collection(idx)%cordinates(1)%Xi = xi
    this%VolPrb%collection(idx)%cordinates(1)%Xe = xe
    this%VolPrb%collection(idx)%cordinates(1)%Yi = yi
    this%VolPrb%collection(idx)%cordinates(1)%Ye = ye
    this%VolPrb%collection(idx)%cordinates(1)%Zi = zi
    this%VolPrb%collection(idx)%cordinates(1)%Ze = ze
    this%VolPrb%collection(idx)%cordinates(1)%Xtrancos = xtrancos
    this%VolPrb%collection(idx)%cordinates(1)%Ytrancos = ytrancos
    this%VolPrb%collection(idx)%cordinates(1)%Ztrancos = ztrancos
    this%VolPrb%collection(idx)%cordinates(1)%Or = or
    this%VolPrb%collection(idx)%cordinates(1)%tag = tag
end subroutine init_volumic_probe_data

subroutine verify_volumic_probe_rotation(test_err, this, mpidir)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    integer, intent(in) :: mpidir
    if (mpidir==2) then
    call expect_eq_int(test_err, 1, this%VolPrb%collection(1)%cordinates(1)%Xi, "rotate_generateVolumicProbes: Xi(1) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(1)%cordinates(1)%Yi, "rotate_generateVolumicProbes: Yi(1) should remain unchanged")
    call expect_eq_int(test_err, 3, this%VolPrb%collection(1)%cordinates(1)%Zi, "rotate_generateVolumicProbes: Zi(1) should remain unchanged")
    call expect_eq_int(test_err, 4, this%VolPrb%collection(1)%cordinates(1)%Xe, "rotate_generateVolumicProbes: Xe(1) should remain unchanged")
    call expect_eq_int(test_err, 5, this%VolPrb%collection(1)%cordinates(1)%Ye, "rotate_generateVolumicProbes: Ye(1) should remain unchanged")
    call expect_eq_int(test_err, 6, this%VolPrb%collection(1)%cordinates(1)%Ze, "rotate_generateVolumicProbes: Ze(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%VolPrb%collection(1)%cordinates(1)%Xtrancos, "rotate_generateVolumicProbes: Xtrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(1)%cordinates(1)%Ytrancos, "rotate_generateVolumicProbes: Ytrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 3, this%VolPrb%collection(1)%cordinates(1)%Ztrancos, "rotate_generateVolumicProbes: Ztrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%VolPrb%collection(1)%cordinates(1)%Or, "rotate_generateVolumicProbes: Or(1) should remain unchanged")
    
    call expect_eq_int(test_err, 5, this%VolPrb%collection(2)%cordinates(1)%Xi, "rotate_generateVolumicProbes: Xi(2) should remain unchanged")
    call expect_eq_int(test_err, 6, this%VolPrb%collection(2)%cordinates(1)%Yi, "rotate_generateVolumicProbes: Yi(2) should remain unchanged")
    call expect_eq_int(test_err, 7, this%VolPrb%collection(2)%cordinates(1)%Zi, "rotate_generateVolumicProbes: Zi(2) should remain unchanged")
    call expect_eq_int(test_err, 8, this%VolPrb%collection(2)%cordinates(1)%Xe, "rotate_generateVolumicProbes: Xe(2) should remain unchanged")
    call expect_eq_int(test_err, 9, this%VolPrb%collection(2)%cordinates(1)%Ye, "rotate_generateVolumicProbes: Ye(2) should remain unchanged")
    call expect_eq_int(test_err, 10, this%VolPrb%collection(2)%cordinates(1)%Ze, "rotate_generateVolumicProbes: Ze(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%VolPrb%collection(2)%cordinates(1)%Xtrancos, "rotate_generateVolumicProbes: Xtrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(2)%cordinates(1)%Ytrancos, "rotate_generateVolumicProbes: Ytrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 3, this%VolPrb%collection(2)%cordinates(1)%Ztrancos, "rotate_generateVolumicProbes: Ztrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(2)%cordinates(1)%Or, "rotate_generateVolumicProbes: Or(2) should remain unchanged")
    end if
    if (mpidir==1) then
    call expect_eq_int(test_err, 1, this%VolPrb%collection(1)%cordinates(1)%Xi, "rotate_generateVolumicProbes: Xi(1) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(1)%cordinates(1)%Yi, "rotate_generateVolumicProbes: Yi(1) should remain unchanged")
    call expect_eq_int(test_err, 3, this%VolPrb%collection(1)%cordinates(1)%Zi, "rotate_generateVolumicProbes: Zi(1) should remain unchanged")
    call expect_eq_int(test_err, 4, this%VolPrb%collection(1)%cordinates(1)%Xe, "rotate_generateVolumicProbes: Xe(1) should remain unchanged")
    call expect_eq_int(test_err, 5, this%VolPrb%collection(1)%cordinates(1)%Ye, "rotate_generateVolumicProbes: Ye(1) should remain unchanged")
    call expect_eq_int(test_err, 6, this%VolPrb%collection(1)%cordinates(1)%Ze, "rotate_generateVolumicProbes: Ze(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%VolPrb%collection(1)%cordinates(1)%Xtrancos, "rotate_generateVolumicProbes: Xtrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(1)%cordinates(1)%Ytrancos, "rotate_generateVolumicProbes: Ytrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 3, this%VolPrb%collection(1)%cordinates(1)%Ztrancos, "rotate_generateVolumicProbes: Ztrancos(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%VolPrb%collection(1)%cordinates(1)%Or, "rotate_generateVolumicProbes: Or(1) should remain unchanged")
    
    call expect_eq_int(test_err, 5, this%VolPrb%collection(2)%cordinates(1)%Xi, "rotate_generateVolumicProbes: Xi(2) should remain unchanged")
    call expect_eq_int(test_err, 6, this%VolPrb%collection(2)%cordinates(1)%Yi, "rotate_generateVolumicProbes: Yi(2) should remain unchanged")
    call expect_eq_int(test_err, 7, this%VolPrb%collection(2)%cordinates(1)%Zi, "rotate_generateVolumicProbes: Zi(2) should remain unchanged")
    call expect_eq_int(test_err, 8, this%VolPrb%collection(2)%cordinates(1)%Xe, "rotate_generateVolumicProbes: Xe(2) should remain unchanged")
    call expect_eq_int(test_err, 9, this%VolPrb%collection(2)%cordinates(1)%Ye, "rotate_generateVolumicProbes: Ye(2) should remain unchanged")
    call expect_eq_int(test_err, 10, this%VolPrb%collection(2)%cordinates(1)%Ze, "rotate_generateVolumicProbes: Ze(2) should remain unchanged")
    call expect_eq_int(test_err, 1, this%VolPrb%collection(2)%cordinates(1)%Xtrancos, "rotate_generateVolumicProbes: Xtrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(2)%cordinates(1)%Ytrancos, "rotate_generateVolumicProbes: Ytrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 3, this%VolPrb%collection(2)%cordinates(1)%Ztrancos, "rotate_generateVolumicProbes: Ztrancos(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%VolPrb%collection(2)%cordinates(1)%Or, "rotate_generateVolumicProbes: Or(2) should remain unchanged")
    end if
end subroutine verify_volumic_probe_rotation

subroutine cleanup_volumic_probe_test(this)
    type(Parseador), intent(inout) :: this
    integer :: i
    
    if (this%VolPrb%length > 0) then
        do i = 1, this%VolPrb%length
            deallocate(this%VolPrb%collection(i)%cordinates)
        end do
        deallocate(this%VolPrb%collection)
    end if
    deallocate(this%VolPrb)
end subroutine cleanup_volumic_probe_test

end module test_rotate_generateVolumicProbes_m 