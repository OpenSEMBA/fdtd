integer function  test_conformal_pec_media() bind(C) result(err)
    use conformal_mod
    use SEMBA_FDTD_mod
    implicit none

    type(semba_fdtd_t) :: semba
   
    err = 0
    call chdir("./testData/cases/conformal/")
    call semba%init("-i conformal_pec_media.fdtd.json")
    call sleep(1)

end function