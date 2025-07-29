module test_rotate_generateBloqueProbes_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_bloque_probes() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    mpidir = 2
    call setup_bloque_probe_test(this, 2)
    
    call init_bloque_probe_data(this, 1, 1, 4, 2, 5, 3, 6, 1, 1, "probe1")
    call init_bloque_probe_data(this, 2, 5, 8, 6, 9, 7, 10, 2, 2, "probe2")
    
    call rotate_generateBloqueProbes(this, mpidir)
    call verify_bloque_probe_rotation(test_err, this, mpidir)
    
    call cleanup_bloque_probe_test(this)
    
    mpidir = 1
    call setup_bloque_probe_test(this, 2)
    
    call init_bloque_probe_data(this, 1, 1, 4, 2, 5, 3, 6, 1, 1, "probe1")
    call init_bloque_probe_data(this, 2, 5, 8, 6, 9, 7, 10, 2, 2, "probe2")
    
    call rotate_generateBloqueProbes(this, mpidir)
    call verify_bloque_probe_rotation(test_err, this, mpidir)
    
    call cleanup_bloque_probe_test(this)
    
    err = test_err
end function test_rotate_generate_bloque_probes

subroutine setup_bloque_probe_test(this, n_probes)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: n_probes
    integer :: i
    
    allocate(this%BloquePrb)
    this%BloquePrb%n_bp = n_probes
    this%BloquePrb%n_bp_max = n_probes
    if (n_probes > 0) then
        allocate(this%BloquePrb%bp(n_probes))

        do i = 1, n_probes
            this%BloquePrb%bp(i)%tstart = 0.0d0
            this%BloquePrb%bp(i)%tstop = 1.0d0
            this%BloquePrb%bp(i)%tstep = 0.1d0
            this%BloquePrb%bp(i)%fstart = 0.0d0
            this%BloquePrb%bp(i)%fstop = 1.0d0
            this%BloquePrb%bp(i)%fstep = 0.1d0
            this%BloquePrb%bp(i)%FileNormalize = "normalize.txt"
            this%BloquePrb%bp(i)%skip = 1
            this%BloquePrb%bp(i)%nml = 1
            this%BloquePrb%bp(i)%t = .true.
            this%BloquePrb%bp(i)%outputrequest = "output.txt"
        end do
    end if
end subroutine setup_bloque_probe_test

subroutine init_bloque_probe_data(this, idx, i1, i2, j1, j2, k1, k2, type2, skip, tag)
    type(Parseador), intent(inout) :: this
    integer, intent(in) :: idx, i1, i2, j1, j2, k1, k2, type2, skip
    character(len=*), intent(in) :: tag
    
    this%BloquePrb%bp(idx)%i1 = i1
    this%BloquePrb%bp(idx)%i2 = i2
    this%BloquePrb%bp(idx)%j1 = j1
    this%BloquePrb%bp(idx)%j2 = j2
    this%BloquePrb%bp(idx)%k1 = k1
    this%BloquePrb%bp(idx)%k2 = k2
    this%BloquePrb%bp(idx)%type2 = type2
    this%BloquePrb%bp(idx)%skip = skip
    this%BloquePrb%bp(idx)%tag = tag
end subroutine init_bloque_probe_data

subroutine verify_bloque_probe_rotation(test_err, this, mpidir)
    integer, intent(inout) :: test_err
    type(Parseador), intent(in) :: this
    integer, intent(in) :: mpidir
    
    if (mpidir == 2) then
    call expect_eq_int(test_err, 3, this%BloquePrb%bp(1)%i1, "rotate_generateBloqueProbes: i1(1) should be rotated")
    call expect_eq_int(test_err, 1, this%BloquePrb%bp(1)%j1, "rotate_generateBloqueProbes: j1(1) should remain unchanged")
    call expect_eq_int(test_err, 2, this%BloquePrb%bp(1)%k1, "rotate_generateBloqueProbes: k1(1) should be rotated")
    call expect_eq_int(test_err, 6, this%BloquePrb%bp(1)%i2, "rotate_generateBloqueProbes: i2(1) should be rotated")
    call expect_eq_int(test_err, 4, this%BloquePrb%bp(1)%j2, "rotate_generateBloqueProbes: j2(1) should remain unchanged")
    call expect_eq_int(test_err, 5, this%BloquePrb%bp(1)%k2, "rotate_generateBloqueProbes: k2(1) should be rotated")
    call expect_eq_int(test_err, 1, this%BloquePrb%bp(1)%type2, "rotate_generateBloqueProbes: type2(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%BloquePrb%bp(1)%skip, "rotate_generateBloqueProbes: skip(1) should remain unchanged")
    
    call expect_eq_int(test_err, 7, this%BloquePrb%bp(2)%i1, "rotate_generateBloqueProbes: i1(2) should be rotated")
    call expect_eq_int(test_err, 5, this%BloquePrb%bp(2)%j1, "rotate_generateBloqueProbes: j1(2) should remain unchanged")
    call expect_eq_int(test_err, 6, this%BloquePrb%bp(2)%k1, "rotate_generateBloqueProbes: k1(2) should be rotated")
    call expect_eq_int(test_err, 10, this%BloquePrb%bp(2)%i2, "rotate_generateBloqueProbes: i2(2) should be rotated")
    call expect_eq_int(test_err, 8, this%BloquePrb%bp(2)%j2, "rotate_generateBloqueProbes: j2(2) should remain unchanged")
    call expect_eq_int(test_err, 9, this%BloquePrb%bp(2)%k2, "rotate_generateBloqueProbes: k2(2) should be rotated")
    call expect_eq_int(test_err, 2, this%BloquePrb%bp(2)%type2, "rotate_generateBloqueProbes: type2(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%BloquePrb%bp(2)%skip, "rotate_generateBloqueProbes: skip(2) should remain unchanged")
    end if

    if (mpidir == 1) then 
    call expect_eq_int(test_err, 2, this%BloquePrb%bp(1)%i1, "rotate_generateBloqueProbes: i1(1) should remain unchanged")
    call expect_eq_int(test_err, 3, this%BloquePrb%bp(1)%j1, "rotate_generateBloqueProbes: j1(1) should be rotated")
    call expect_eq_int(test_err, 1, this%BloquePrb%bp(1)%k1, "rotate_generateBloqueProbes: k1(1) should be rotated")
    call expect_eq_int(test_err, 5, this%BloquePrb%bp(1)%i2, "rotate_generateBloqueProbes: i2(1) should remain unchanged")
    call expect_eq_int(test_err, 6, this%BloquePrb%bp(1)%j2, "rotate_generateBloqueProbes: j2(1) should be rotated")
    call expect_eq_int(test_err, 4, this%BloquePrb%bp(1)%k2, "rotate_generateBloqueProbes: k2(1) should be rotated")
    call expect_eq_int(test_err, 1, this%BloquePrb%bp(1)%type2, "rotate_generateBloqueProbes: type2(1) should remain unchanged")
    call expect_eq_int(test_err, 1, this%BloquePrb%bp(1)%skip, "rotate_generateBloqueProbes: skip(1) should remain unchanged")
    
    call expect_eq_int(test_err, 6, this%BloquePrb%bp(2)%i1, "rotate_generateBloqueProbes: i1(2) should remain unchanged")
    call expect_eq_int(test_err, 7, this%BloquePrb%bp(2)%j1, "rotate_generateBloqueProbes: j1(2) should be rotated")
    call expect_eq_int(test_err, 5, this%BloquePrb%bp(2)%k1, "rotate_generateBloqueProbes: k1(2) should be rotated")
    call expect_eq_int(test_err, 9, this%BloquePrb%bp(2)%i2, "rotate_generateBloqueProbes: i2(2) should remain unchanged")
    call expect_eq_int(test_err, 10, this%BloquePrb%bp(2)%j2, "rotate_generateBloqueProbes: j2(2) should be rotated")
    call expect_eq_int(test_err, 8, this%BloquePrb%bp(2)%k2, "rotate_generateBloqueProbes: k2(2) should be rotated")
    call expect_eq_int(test_err, 2, this%BloquePrb%bp(2)%type2, "rotate_generateBloqueProbes: type2(2) should remain unchanged")
    call expect_eq_int(test_err, 2, this%BloquePrb%bp(2)%skip, "rotate_generateBloqueProbes: skip(2) should remain unchanged")
    end if
end subroutine verify_bloque_probe_rotation

subroutine cleanup_bloque_probe_test(this)
    type(Parseador), intent(inout) :: this
    
    if (this%BloquePrb%n_bp > 0) then
        deallocate(this%BloquePrb%bp)
    end if
    deallocate(this%BloquePrb)
end subroutine cleanup_bloque_probe_test

end module test_rotate_generateBloqueProbes_m 