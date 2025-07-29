module test_rotate_generateFDMs_m
    use smbjson
    use nfde_rotate_m
    use rotate_testingTools
    implicit none

contains

integer function test_rotate_generate_fdms() bind(C) result(err)
    type(Parseador) :: this
    integer(kind=4) :: mpidir
    integer :: test_err = 0
    
    mpidir = 2
    call setup_fdm_test(this)
    call rotate_generateFDMs(this, mpidir)

    call verify_fdm_rotation(test_err, this%frqDepMats, mpidir)
    call cleanup_fdm_test(this)

    mpidir = 1
    call setup_fdm_test(this)
    call rotate_generateFDMs(this, mpidir)

    call verify_fdm_rotation(test_err, this%frqDepMats, mpidir)
    call cleanup_fdm_test(this)
    err = test_err
end function test_rotate_generate_fdms

subroutine setup_fdm_test(this)
    type(Parseador), intent(inout) :: this
    allocate(this%frqDepMats)
    allocate(this%frqDepMats%Vols(1))
    allocate(this%frqDepMats%Vols(1)%a11(1), source=(11,0))
    allocate(this%frqDepMats%Vols(1)%b11(1), source=(11,0))
    allocate(this%frqDepMats%Vols(1)%am11(1), source=(11,0))
    allocate(this%frqDepMats%Vols(1)%bm11(1), source=(11,0))
    allocate(this%frqDepMats%Vols(1)%a12(1), source=(12,0))
    allocate(this%frqDepMats%Vols(1)%b12(1), source=(12,0))
    allocate(this%frqDepMats%Vols(1)%am12(1), source=(12,0))
    allocate(this%frqDepMats%Vols(1)%bm12(1), source=(12,0))
    allocate(this%frqDepMats%Vols(1)%a13(1), source=(13,0))
    allocate(this%frqDepMats%Vols(1)%b13(1), source=(13,0))
    allocate(this%frqDepMats%Vols(1)%am13(1), source=(13,0))
    allocate(this%frqDepMats%Vols(1)%bm13(1), source=(13,0))
    allocate(this%frqDepMats%Vols(1)%a22(1), source=(22,0))
    allocate(this%frqDepMats%Vols(1)%b22(1), source=(22,0))
    allocate(this%frqDepMats%Vols(1)%am22(1), source=(22,0))
    allocate(this%frqDepMats%Vols(1)%bm22(1), source=(22,0))
    allocate(this%frqDepMats%Vols(1)%a23(1), source=(23,0))
    allocate(this%frqDepMats%Vols(1)%b23(1), source=(23,0))
    allocate(this%frqDepMats%Vols(1)%am23(1), source=(23,0))
    allocate(this%frqDepMats%Vols(1)%bm23(1), source=(23,0))
    allocate(this%frqDepMats%Vols(1)%a33(1), source=(33,0))
    allocate(this%frqDepMats%Vols(1)%b33(1), source=(33,0))
    allocate(this%frqDepMats%Vols(1)%am33(1), source=(33,0))
    allocate(this%frqDepMats%Vols(1)%bm33(1), source=(33,0))

    this%frqDepMats%nVols=1

    this%frqDepMats%Vols(1)%eps11 = 11
    this%frqDepMats%Vols(1)%eps12 = 12
    this%frqDepMats%Vols(1)%eps13 = 13
    this%frqDepMats%Vols(1)%eps22 = 22
    this%frqDepMats%Vols(1)%eps23 = 23
    this%frqDepMats%Vols(1)%eps33 = 33
    
    this%frqDepMats%Vols(1)%mu11 = 11
    this%frqDepMats%Vols(1)%mu12 = 12
    this%frqDepMats%Vols(1)%mu13 = 13
    this%frqDepMats%Vols(1)%mu22 = 22
    this%frqDepMats%Vols(1)%mu23 = 23
    this%frqDepMats%Vols(1)%mu33 = 33

    this%frqDepMats%Vols(1)%sigma11 = 11
    this%frqDepMats%Vols(1)%sigma12 = 12
    this%frqDepMats%Vols(1)%sigma13 = 13
    this%frqDepMats%Vols(1)%sigma22 = 22
    this%frqDepMats%Vols(1)%sigma23 = 23
    this%frqDepMats%Vols(1)%sigma33 = 33

    this%frqDepMats%Vols(1)%sigmam11 = 11
    this%frqDepMats%Vols(1)%sigmam12 = 12
    this%frqDepMats%Vols(1)%sigmam13 = 13
    this%frqDepMats%Vols(1)%sigmam22 = 22
    this%frqDepMats%Vols(1)%sigmam23 = 23
    this%frqDepMats%Vols(1)%sigmam33 = 33

    this%frqDepMats%Vols(1)%K11 = 11
    this%frqDepMats%Vols(1)%K12 = 12
    this%frqDepMats%Vols(1)%K13 = 13
    this%frqDepMats%Vols(1)%K22 = 22
    this%frqDepMats%Vols(1)%K23 = 23
    this%frqDepMats%Vols(1)%K33 = 33

    this%frqDepMats%Vols(1)%Km11 = 11
    this%frqDepMats%Vols(1)%Km12 = 12
    this%frqDepMats%Vols(1)%Km13 = 13
    this%frqDepMats%Vols(1)%Km22 = 22
    this%frqDepMats%Vols(1)%Km23 = 23
    this%frqDepMats%Vols(1)%Km33 = 33

end subroutine setup_fdm_test

subroutine verify_fdm_rotation(test_err, freqDepMat, mpidir)
    integer, intent(inout) :: test_err, mpidir
    type(FreqDepenMaterials), intent(in) :: freqDepMat

    if (mpidir==2) then
        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%a11(1), "Error parsing a11")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%a12(1), "Error parsing a12")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%a13(1), "Error parsing a13")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%a22(1), "Error parsing a22")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%a23(1), "Error parsing a23")
        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%a33(1), "Error parsing a33")

        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%am11(1), "Error parsing am11")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%am12(1), "Error parsing am12")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%am13(1), "Error parsing am13")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%am22(1), "Error parsing am22")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%am23(1), "Error parsing am23")
        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%am33(1), "Error parsing am33")

        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%b11(1), "Error parsing b11")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%b12(1), "Error parsing b12")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%b13(1), "Error parsing b13")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%b22(1), "Error parsing b22")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%b23(1), "Error parsing b23")
        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%b33(1), "Error parsing b33")

        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%bm11(1), "Error parsing bm11")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%bm12(1), "Error parsing bm12")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%bm13(1), "Error parsing bm13")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%bm22(1), "Error parsing bm22")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%bm23(1), "Error parsing bm23")
        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%bm33(1), "Error parsing bm33")

        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%eps11, "Error parsing eps11")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%eps12, "Error parsing eps12")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%eps13, "Error parsing eps13")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%eps22, "Error parsing eps22")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%eps23, "Error parsing eps23")
        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%eps33, "Error parsing eps33")

        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%mu11, "Error parsing mu11")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%mu12, "Error parsing mu12")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%mu13, "Error parsing mu13")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%mu22, "Error parsing mu22")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%mu23, "Error parsing mu23")
        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%mu33, "Error parsing mu33")

        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%sigma11, "Error parsing sigma11")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%sigma12, "Error parsing sigma12")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%sigma13, "Error parsing sigma13")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%sigma22, "Error parsing sigma22")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%sigma23, "Error parsing sigma23")
        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%sigma33, "Error parsing sigma33")

        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%sigmam11, "Error parsing sigmam11")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%sigmam12, "Error parsing sigmam12")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%sigmam13, "Error parsing sigmam13")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%sigmam22, "Error parsing sigmam22")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%sigmam23, "Error parsing sigmam23")
        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%sigmam33, "Error parsing sigmam33")

        call expect_eq_int(test_err, 33, freqDepMat%Vols(1)%k11, "Error parsing k11")
        call expect_eq_int(test_err, 23, freqDepMat%Vols(1)%k12, "Error parsing k12")
        call expect_eq_int(test_err, 12, freqDepMat%Vols(1)%k13, "Error parsing k13")
        call expect_eq_int(test_err, 11, freqDepMat%Vols(1)%k22, "Error parsing k22")
        call expect_eq_int(test_err, 13, freqDepMat%Vols(1)%k23, "Error parsing k23")
        call expect_eq_int(test_err, 22, freqDepMat%Vols(1)%k33, "Error parsing k33")

        call expect_eq_int(test_err, 33, freqDepMat%Vols(1)%km11, "Error parsing km11")
        call expect_eq_int(test_err, 23, freqDepMat%Vols(1)%km12, "Error parsing km12")
        call expect_eq_int(test_err, 12, freqDepMat%Vols(1)%km13, "Error parsing km13")
        call expect_eq_int(test_err, 11, freqDepMat%Vols(1)%km22, "Error parsing km22")
        call expect_eq_int(test_err, 13, freqDepMat%Vols(1)%km23, "Error parsing km23")
        call expect_eq_int(test_err, 22, freqDepMat%Vols(1)%km33, "Error parsing km33")
    end if

    if (mpidir==1) then
        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%a11(1), "Error parsing a11")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%a12(1), "Error parsing a12")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%a13(1), "Error parsing a13")
        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%a22(1), "Error parsing a22")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%a23(1), "Error parsing a23")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%a33(1), "Error parsing a33")

        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%am11(1), "Error parsing am11")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%am12(1), "Error parsing am12")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%am13(1), "Error parsing am13")
        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%am22(1), "Error parsing am22")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%am23(1), "Error parsing am23")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%am33(1), "Error parsing am33")

        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%b11(1), "Error parsing b11")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%b12(1), "Error parsing b12")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%b13(1), "Error parsing b13")
        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%b22(1), "Error parsing b22")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%b23(1), "Error parsing b23")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%b33(1), "Error parsing b33")

        call expect_eq_complx(test_err, (22,0), freqDepMat%Vols(1)%bm11(1), "Error parsing bm11")
        call expect_eq_complx(test_err, (13,0), freqDepMat%Vols(1)%bm12(1), "Error parsing bm12")
        call expect_eq_complx(test_err, (23,0), freqDepMat%Vols(1)%bm13(1), "Error parsing bm13")
        call expect_eq_complx(test_err, (33,0), freqDepMat%Vols(1)%bm22(1), "Error parsing bm22")
        call expect_eq_complx(test_err, (12,0), freqDepMat%Vols(1)%bm23(1), "Error parsing bm23")
        call expect_eq_complx(test_err, (11,0), freqDepMat%Vols(1)%bm33(1), "Error parsing bm33")

        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%eps11, "Error parsing eps11")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%eps12, "Error parsing eps12")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%eps13, "Error parsing eps13")
        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%eps22, "Error parsing eps22")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%eps23, "Error parsing eps23")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%eps33, "Error parsing eps33")

        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%mu11, "Error parsing mu11")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%mu12, "Error parsing mu12")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%mu13, "Error parsing mu13")
        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%mu22, "Error parsing mu22")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%mu23, "Error parsing mu23")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%mu33, "Error parsing mu33")

        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%sigma11, "Error parsing sigma11")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%sigma12, "Error parsing sigma12")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%sigma13, "Error parsing sigma13")
        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%sigma22, "Error parsing sigma22")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%sigma23, "Error parsing sigma23")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%sigma33, "Error parsing sigma33")

        call expect_eq_real(test_err, 22.0_RKIND, freqDepMat%Vols(1)%sigmam11, "Error parsing sigmam11")
        call expect_eq_real(test_err, 13.0_RKIND, freqDepMat%Vols(1)%sigmam12, "Error parsing sigmam12")
        call expect_eq_real(test_err, 23.0_RKIND, freqDepMat%Vols(1)%sigmam13, "Error parsing sigmam13")
        call expect_eq_real(test_err, 33.0_RKIND, freqDepMat%Vols(1)%sigmam22, "Error parsing sigmam22")
        call expect_eq_real(test_err, 12.0_RKIND, freqDepMat%Vols(1)%sigmam23, "Error parsing sigmam23")
        call expect_eq_real(test_err, 11.0_RKIND, freqDepMat%Vols(1)%sigmam33, "Error parsing sigmam33")

        call expect_eq_int(test_err, 22, freqDepMat%Vols(1)%k11, "Error parsing k11")
        call expect_eq_int(test_err, 13, freqDepMat%Vols(1)%k12, "Error parsing k12")
        call expect_eq_int(test_err, 23, freqDepMat%Vols(1)%k13, "Error parsing k13")
        call expect_eq_int(test_err, 33, freqDepMat%Vols(1)%k22, "Error parsing k22")
        call expect_eq_int(test_err, 12, freqDepMat%Vols(1)%k23, "Error parsing k23")
        call expect_eq_int(test_err, 11, freqDepMat%Vols(1)%k33, "Error parsing k33")

        call expect_eq_int(test_err, 22, freqDepMat%Vols(1)%km11, "Error parsing km11")
        call expect_eq_int(test_err, 13, freqDepMat%Vols(1)%km12, "Error parsing km12")
        call expect_eq_int(test_err, 23, freqDepMat%Vols(1)%km13, "Error parsing km13")
        call expect_eq_int(test_err, 33, freqDepMat%Vols(1)%km22, "Error parsing km22")
        call expect_eq_int(test_err, 12, freqDepMat%Vols(1)%km23, "Error parsing km23")
        call expect_eq_int(test_err, 11, freqDepMat%Vols(1)%km33, "Error parsing km33")
    end if
end subroutine verify_fdm_rotation

subroutine cleanup_fdm_test(this)
    type(Parseador), intent(inout) :: this
    deallocate(this%frqDepMats%Vols)

end subroutine cleanup_fdm_test

end module test_rotate_generateFDMs_m