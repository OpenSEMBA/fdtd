module system_testingTools_mod
    use iso_c_binding
    use fdetypes
    implicit none
    
    character(len=*), parameter :: PATH_TO_TEST_DATA = 'testData/'
    character(len=*), parameter :: INPUT_EXAMPLES='input_examples/'

contains
    function initMedia(sgg, sggMiNo, sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz, sggMtag) result(res)
        type(sggfdtdinfo), intent(in) :: sgg
        type(media_matrices_t) :: res
        integer (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in)   ::  &
          sggMiNo(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE), &
          sggMiEx(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), &
          sggMiEy(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), &
          sggMiEz(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), &
          sggMiHx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE), &
          sggMiHy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE), &
          sggMiHz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
        integer (KIND=IKINDMTAG), intent(in)   ::  &
          sggMtag(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)

        res%sggMiNo = sggMiNo
        res%sggMiEx = sggMiEx
        res%sggMiEy = sggMiEy
        res%sggMiEz = sggMiEz
        res%sggMiHx = sggMiHx
        res%sggMiHy = sggMiHy
        res%sggMiHz = sggMiHz
        res%sggMtag = sggMtag

    end function

 end module system_testingTools_mod