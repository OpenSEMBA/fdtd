integer function  test_conformal_pec_media() bind(C) result(err)
    use SEMBA_FDTD_mod
    implicit none

    type(semba_fdtd_t) :: semba
    type(solver_t) :: solver

    err = 0
    call chdir("./testData/cases/conformal/")
    call semba%init("-i conformal_sphere_1mm_rcs_delta.fdtd.json -mapvtk -n 1")

    ! check inside/outside and media
    
    !conformal edges
    if (semba%media%sggmiEx(14,15,15) /= 4) err = err + 1
    if (semba%media%sggmiEx(16,15,15) /= 4) err = err + 1
    if (semba%media%sggmiEx(14,15,16) /= 4) err = err + 1
    if (semba%media%sggmiEx(16,15,16) /= 4) err = err + 1
    if (semba%media%sggmiEx(14,16,15) /= 4) err = err + 1
    if (semba%media%sggmiEx(16,16,15) /= 4) err = err + 1
    if (semba%media%sggmiEx(14,16,16) /= 4) err = err + 1
    if (semba%media%sggmiEx(16,16,16) /= 4) err = err + 1

    if (semba%media%sggmiEy(15,14,15) /= 4) err = err + 1
    if (semba%media%sggmiEy(15,16,15) /= 4) err = err + 1
    if (semba%media%sggmiEy(15,14,16) /= 4) err = err + 1
    if (semba%media%sggmiEy(15,16,16) /= 4) err = err + 1
    if (semba%media%sggmiEy(16,14,15) /= 4) err = err + 1
    if (semba%media%sggmiEy(16,16,15) /= 4) err = err + 1
    if (semba%media%sggmiEy(16,14,16) /= 4) err = err + 1
    if (semba%media%sggmiEy(16,16,16) /= 4) err = err + 1

    if (semba%media%sggmiEz(15,15,14) /= 4) err = err + 1
    if (semba%media%sggmiEz(15,15,16) /= 4) err = err + 1
    if (semba%media%sggmiEz(16,15,14) /= 4) err = err + 1
    if (semba%media%sggmiEz(16,15,16) /= 4) err = err + 1
    if (semba%media%sggmiEz(15,16,14) /= 4) err = err + 1
    if (semba%media%sggmiEz(15,16,16) /= 4) err = err + 1
    if (semba%media%sggmiEz(16,16,14) /= 4) err = err + 1
    if (semba%media%sggmiEz(16,16,16) /= 4) err = err + 1

    !pec edges
    if (semba%media%sggmiEx(15,15,15) /= 0) err = err + 1
    if (semba%media%sggmiEy(15,15,15) /= 0) err = err + 1
    if (semba%media%sggmiEz(15,15,15) /= 0) err = err + 1
    if (semba%media%sggmiEx(15,16,15) /= 0) err = err + 1
    if (semba%media%sggmiEy(16,15,15) /= 0) err = err + 1
    if (semba%media%sggmiEz(16,16,15) /= 0) err = err + 1

    ! pec faces
    if (semba%media%sggmiHx(15,15,15) /= 0) err = err + 1
    if (semba%media%sggmiHy(15,15,15) /= 0) err = err + 1
    if (semba%media%sggmiHz(15,15,15) /= 0) err = err + 1
    if (semba%media%sggmiHx(16,15,15) /= 0) err = err + 1
    if (semba%media%sggmiHy(15,16,15) /= 0) err = err + 1
    if (semba%media%sggmiHz(15,15,16) /= 0) err = err + 1
    ! conformal faces
    if (semba%media%sggmiHx(15,14,14) /= 6) err = err + 1
    if (semba%media%sggmiHx(15,15,14) /= 5) err = err + 1
    if (semba%media%sggmiHx(15,16,14) /= 6) err = err + 1
    if (semba%media%sggmiHx(15,14,15) /= 5) err = err + 1
    if (semba%media%sggmiHx(15,16,15) /= 5) err = err + 1
    if (semba%media%sggmiHx(15,14,16) /= 6) err = err + 1
    if (semba%media%sggmiHx(15,15,16) /= 5) err = err + 1
    if (semba%media%sggmiHx(15,16,16) /= 6) err = err + 1

    if (semba%media%sggmiHx(16,14,14) /= 6) err = err + 1
    if (semba%media%sggmiHx(16,15,14) /= 5) err = err + 1
    if (semba%media%sggmiHx(16,16,14) /= 6) err = err + 1
    if (semba%media%sggmiHx(16,14,15) /= 5) err = err + 1
    if (semba%media%sggmiHx(16,16,15) /= 5) err = err + 1
    if (semba%media%sggmiHx(16,14,16) /= 6) err = err + 1
    if (semba%media%sggmiHx(16,15,16) /= 5) err = err + 1
    if (semba%media%sggmiHx(16,16,16) /= 6) err = err + 1

    if (semba%media%sggmiHy(14,15,14) /= 6) err = err + 1
    if (semba%media%sggmiHy(15,15,14) /= 5) err = err + 1
    if (semba%media%sggmiHy(16,15,14) /= 6) err = err + 1
    if (semba%media%sggmiHy(14,15,15) /= 5) err = err + 1
    if (semba%media%sggmiHy(16,15,15) /= 5) err = err + 1
    if (semba%media%sggmiHy(14,15,16) /= 6) err = err + 1
    if (semba%media%sggmiHy(15,15,16) /= 5) err = err + 1
    if (semba%media%sggmiHy(16,15,16) /= 6) err = err + 1

    if (semba%media%sggmiHz(14,14,15) /= 6) err = err + 1
    if (semba%media%sggmiHz(15,14,15) /= 5) err = err + 1
    if (semba%media%sggmiHz(16,14,15) /= 6) err = err + 1
    if (semba%media%sggmiHz(14,15,15) /= 5) err = err + 1
    if (semba%media%sggmiHz(16,15,15) /= 5) err = err + 1
    if (semba%media%sggmiHz(14,16,15) /= 6) err = err + 1
    if (semba%media%sggmiHz(15,16,15) /= 5) err = err + 1
    if (semba%media%sggmiHz(16,16,15) /= 6) err = err + 1

    if (semba%media%sggmiHz(14,14,16) /= 6) err = err + 1
    if (semba%media%sggmiHz(15,14,16) /= 5) err = err + 1
    if (semba%media%sggmiHz(16,14,16) /= 6) err = err + 1
    if (semba%media%sggmiHz(14,15,16) /= 5) err = err + 1
    if (semba%media%sggmiHz(16,15,16) /= 5) err = err + 1
    if (semba%media%sggmiHz(14,16,16) /= 6) err = err + 1
    if (semba%media%sggmiHz(15,16,16) /= 5) err = err + 1
    if (semba%media%sggmiHz(16,16,16) /= 6) err = err + 1

    call semba%launch()
    call semba%end()

    ! check inside/outside
    ! if (semba%media%sggmiEx(14,14,14) /= 0) err = err + 1

end function

integer function  test_conformal_pec_corner() bind(C) result(err)
    use SEMBA_FDTD_mod
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

    if (semba%l%fatalerror .eqv. .true.) err = err + 1

end function

integer function test_conformal_pec_media_raytracing() bind(C) result(err)
    use SEMBA_FDTD_mod
    implicit none

    type(semba_fdtd_t) :: semba
    type(solver_t) :: solver

    err = 0
    call chdir("./testData/cases/conformal/")
    call semba%init("-i conformal_sphere_30.fdtd.json -mapvtk -n 1")
    call semba%launch()
    call semba%end()

    if (semba%l%fatalerror .eqv. .true.) err = err + 1


end function

