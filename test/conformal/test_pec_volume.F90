integer function  test_conformal_pec_media() bind(C) result(err)
    use SEMBA_FDTD_mod
    ! use conformal_mod
    implicit none

    type(semba_fdtd_t) :: semba
    type(solver_t) :: solver

    err = 0
    call chdir("./testData/cases/conformal/")
    call semba%init("-i conformal_sphere_1mm_rcs_delta.fdtd.json -mapvtk -n 1")

    ! check inside/outside and media
    
    !conformal edges
    if (semba%media%sggmiEx(1,2,2) /= 4) err = err + 1
    if (semba%media%sggmiEx(3,2,2) /= 4) err = err + 1
    if (semba%media%sggmiEx(1,2,3) /= 4) err = err + 1
    if (semba%media%sggmiEx(3,2,3) /= 4) err = err + 1
    if (semba%media%sggmiEx(1,3,2) /= 4) err = err + 1
    if (semba%media%sggmiEx(3,3,2) /= 4) err = err + 1
    if (semba%media%sggmiEx(1,3,3) /= 4) err = err + 1
    if (semba%media%sggmiEx(3,3,3) /= 4) err = err + 1

    if (semba%media%sggmiEy(2,1,2) /= 4) err = err + 1
    if (semba%media%sggmiEy(2,3,2) /= 4) err = err + 1
    if (semba%media%sggmiEy(2,1,3) /= 4) err = err + 1
    if (semba%media%sggmiEy(2,3,3) /= 4) err = err + 1
    if (semba%media%sggmiEy(3,1,2) /= 4) err = err + 1
    if (semba%media%sggmiEy(3,3,2) /= 4) err = err + 1
    if (semba%media%sggmiEy(3,1,3) /= 4) err = err + 1
    if (semba%media%sggmiEy(3,3,3) /= 4) err = err + 1

    if (semba%media%sggmiEz(2,2,1) /= 4) err = err + 1
    if (semba%media%sggmiEz(2,2,3) /= 4) err = err + 1
    if (semba%media%sggmiEz(3,2,1) /= 4) err = err + 1
    if (semba%media%sggmiEz(3,2,3) /= 4) err = err + 1
    if (semba%media%sggmiEz(2,3,1) /= 4) err = err + 1
    if (semba%media%sggmiEz(2,3,3) /= 4) err = err + 1
    if (semba%media%sggmiEz(3,3,1) /= 4) err = err + 1
    if (semba%media%sggmiEz(3,3,3) /= 4) err = err + 1

    !pec edges
    if (semba%media%sggmiEx(2,2,2) /= 0) err = err + 1
    if (semba%media%sggmiEy(2,2,2) /= 0) err = err + 1
    if (semba%media%sggmiEz(2,2,2) /= 0) err = err + 1
    if (semba%media%sggmiEx(2,3,2) /= 0) err = err + 1
    if (semba%media%sggmiEy(3,2,2) /= 0) err = err + 1
    if (semba%media%sggmiEz(3,3,2) /= 0) err = err + 1

    ! pec faces
    if (semba%media%sggmiHx(2,2,2) /= 0) err = err + 1
    if (semba%media%sggmiHy(2,2,2) /= 0) err = err + 1
    if (semba%media%sggmiHz(2,2,2) /= 0) err = err + 1
    if (semba%media%sggmiHx(3,2,2) /= 0) err = err + 1
    if (semba%media%sggmiHy(2,3,2) /= 0) err = err + 1
    if (semba%media%sggmiHz(2,2,3) /= 0) err = err + 1
    ! conformal faces
    if (semba%media%sggmiHx(2,1,1) /= 6) err = err + 1
    if (semba%media%sggmiHx(2,2,1) /= 5) err = err + 1
    if (semba%media%sggmiHx(2,3,1) /= 6) err = err + 1
    if (semba%media%sggmiHx(2,1,2) /= 5) err = err + 1
    if (semba%media%sggmiHx(2,3,2) /= 5) err = err + 1
    if (semba%media%sggmiHx(2,1,3) /= 6) err = err + 1
    if (semba%media%sggmiHx(2,2,3) /= 5) err = err + 1
    if (semba%media%sggmiHx(2,3,3) /= 6) err = err + 1

    if (semba%media%sggmiHx(3,1,1) /= 6) err = err + 1
    if (semba%media%sggmiHx(3,2,1) /= 5) err = err + 1
    if (semba%media%sggmiHx(3,3,1) /= 6) err = err + 1
    if (semba%media%sggmiHx(3,1,2) /= 5) err = err + 1
    if (semba%media%sggmiHx(3,3,2) /= 5) err = err + 1
    if (semba%media%sggmiHx(3,1,3) /= 6) err = err + 1
    if (semba%media%sggmiHx(3,2,3) /= 5) err = err + 1
    if (semba%media%sggmiHx(3,3,3) /= 6) err = err + 1

    if (semba%media%sggmiHy(1,2,1) /= 6) err = err + 1
    if (semba%media%sggmiHy(2,2,1) /= 5) err = err + 1
    if (semba%media%sggmiHy(3,2,1) /= 6) err = err + 1
    if (semba%media%sggmiHy(1,2,2) /= 5) err = err + 1
    if (semba%media%sggmiHy(3,2,2) /= 5) err = err + 1
    if (semba%media%sggmiHy(1,2,3) /= 6) err = err + 1
    if (semba%media%sggmiHy(2,2,3) /= 5) err = err + 1
    if (semba%media%sggmiHy(3,2,3) /= 6) err = err + 1

    if (semba%media%sggmiHz(1,1,2) /= 6) err = err + 1
    if (semba%media%sggmiHz(2,1,2) /= 5) err = err + 1
    if (semba%media%sggmiHz(3,1,2) /= 6) err = err + 1
    if (semba%media%sggmiHz(1,2,2) /= 5) err = err + 1
    if (semba%media%sggmiHz(3,2,2) /= 5) err = err + 1
    if (semba%media%sggmiHz(1,3,2) /= 6) err = err + 1
    if (semba%media%sggmiHz(2,3,2) /= 5) err = err + 1
    if (semba%media%sggmiHz(3,3,2) /= 6) err = err + 1

    if (semba%media%sggmiHz(1,1,3) /= 6) err = err + 1
    if (semba%media%sggmiHz(2,1,3) /= 5) err = err + 1
    if (semba%media%sggmiHz(3,1,3) /= 6) err = err + 1
    if (semba%media%sggmiHz(1,2,3) /= 5) err = err + 1
    if (semba%media%sggmiHz(3,2,3) /= 5) err = err + 1
    if (semba%media%sggmiHz(1,3,3) /= 6) err = err + 1
    if (semba%media%sggmiHz(2,3,3) /= 5) err = err + 1
    if (semba%media%sggmiHz(3,3,3) /= 6) err = err + 1

    call semba%launch()
    call semba%end()

    ! check inside/outside
    ! if (semba%media%sggmiEx(1,1,1) /= 0) err = err + 1

end function

integer function  test_conformal_pec_corner() bind(C) result(err)
    use SEMBA_FDTD_mod
    use conformal_mod
    implicit none

    type(semba_fdtd_t) :: semba
   
    err = 0
    call chdir("./testData/cases/conformal/")
    call semba%init("-i conformal_corner.fdtd.json -mapvtk -n 1")
    
    ! check inside/outside and media
    if (semba%media%sggmiEx(2,2,2) /= 4) err = err + 1
    if (semba%media%sggmiEx(2,2,3) /= 4) err = err + 1
    if (semba%media%sggmiEy(2,2,2) /= 4) err = err + 1
    if (semba%media%sggmiEy(2,2,3) /= 4) err = err + 1
    if (semba%media%sggmiEz(2,2,2) /= 0) err = err + 1
    
    if (semba%media%sggmiHx(2,2,2) /= 5) err = err + 1
    if (semba%media%sggmiHy(2,2,2) /= 5) err = err + 1
    if (semba%media%sggmiHz(2,2,2) /= 6) err = err + 1
    if (semba%media%sggmiHz(2,2,3) /= 6) err = err + 1
    
    call semba%launch()
    call semba%end()

end function