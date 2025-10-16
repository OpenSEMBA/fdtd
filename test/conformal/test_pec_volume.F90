integer function  test_conformal_pec_media() bind(C) result(err)
    use SEMBA_FDTD_mod
    use conformal_mod
    implicit none

    type(semba_fdtd_t) :: semba
   
    err = 0
    call chdir("./testData/cases/conformal/")
    call semba%init("-i conformal_sphere_1mm_rcs_delta.fdtd.json -mapvtk -n 1")
    call semba%launch()
    call semba%end()

    ! check inside/outside
    ! if (semba%media%sggmiEx(1,1,1) /= 0) err = err + 1

end function