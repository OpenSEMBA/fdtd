
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module :  CreateMatrices :  Fills in the media and observation matrices, creates
!                        the pml layers and the PEC borders.
!                        Also creates intermediate media for the boundaries
!                        between different media.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE CreateMatrices
   USE Report
   USE fdetypes
   !
   IMPLICIT NONE
   PRIVATE
   !
   !
   TYPE crosscheck_t
      INTEGER (KIND=4) :: actual, NewActual, NewActual2
      INTEGER (KIND=4), DIMENSION (1:4) :: tent
   END TYPE
   !matriz para controlar lo punietereos indices de cadacomponente
   INTEGER (KIND=4), DIMENSION (6, 3, 2), PARAMETER, PUBLIC :: &
   & in = reshape ( (/ 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, &
   &                   0, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /), (/ 6, 3, 2 /))
   !
   INTEGER (KIND=4), DIMENSION (6, 3, 2), PARAMETER, PUBLIC :: &
   &    on = reshape ( (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
   &                      0, 0,-1, 0, 0, 0,-1,-1, 0,-1, 0,-1, 0,-1, 0, 0, &
   -1,-1,-1, 0 /), (/ 6, 3, 2 /))
   !
   PUBLIC CreatePMLmatrix, Readjust
   PUBLIC CreateVolumeMM, CreateSurfaceMM, CreateLineMM
   PUBLIC CreateSurfaceSlotMM,CreateMagneticSurface
   public CreateConformalPECVolume
   !
    CONTAINS
    
   SUBROUTINE SortInitEndWithIncreasingOrder(p)
        type(XYZlimit_t), intent(inout) :: p
        integer (kind=4) :: aux
        if (p%XI > p%XE) then
            aux = p%XI
            p%XI = p%XE
            p%XE = aux
        endif
        if (p%YI > p%YE) then
            aux = p%YI
            p%YI = p%YE
            p%YE = aux
        endif
        if (p%ZI > p%ZE) then
            aux = p%ZI
            p%ZI = p%ZE
            p%ZE = aux
        endif
   END SUBROUTINE
   
   SUBROUTINE CreateConformalPECVolume (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, BoundingBox,indicemedio)
      character(len=BUFSIZE) :: buff
      TYPE (Shared_t) :: Eshared

      INTEGER (KIND=4) :: NumMedia
      TYPE (MediaData_t), DIMENSION (0:NumMedia) :: med
      !
      TYPE (XYZlimit_t), INTENT(IN) ::  BoundingBox
      integer(kind=4), intent(in) :: indicemedio
      !
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      INTEGER (KIND=IKINDMTAG) numertag
      INTEGER (KIND=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      !
      INTEGER (KIND=4) :: layoutnumber, i, j, k

      ! faces that should be PEC, bc edges are all PEC
      do k = BoundingBox%zi, BoundingBox%ze+1
         do j = BoundingBox%yi, BoundingBox%ye+1
            do i = BoundingBox%xi, BoundingBox%xe+1
               call fillBoundaryFaceIfAllEdgesPEC(i,j,k, FACE_X)
               call fillBoundaryFaceIfAllEdgesPEC(i,j,k, FACE_Y)
               call fillBoundaryFaceIfAllEdgesPEC(i,j,k, FACE_Z)
            end do
         end do
      end do

      ! raytracing x
      do k = BoundingBox%zi, BoundingBox%ze+1
         do j = BoundingBox%yi, BoundingBox%ye+1
            call fillEdgesInsideVolumeX(j,k)
         end do
      end do
      ! raytracing y
      do i = BoundingBox%xi, BoundingBox%xe+1
         do k = BoundingBox%zi, BoundingBox%ze+1
            call fillEdgesInsideVolumeY(i,k)
         end do
      end do
      ! raytracing z
      do j = BoundingBox%yi, BoundingBox%ye+1
         do i = BoundingBox%xi, BoundingBox%xe+1
            call fillEdgesInsideVolumeZ(i,j)
         end do
      end do

      ! faces inside volume, should be PEC
      do k = BoundingBox%zi, BoundingBox%ze+1
         do j = BoundingBox%yi, BoundingBox%ye+1
            do i = BoundingBox%xi, BoundingBox%xe+1
               call fillPECFaceInsideVolume(i,j,k,FACE_X, MMiHx, MMiEy, MMiEz, tags%face%x)
               call fillPECFaceInsideVolume(i,j,k,FACE_Y, MMiHy, MMiEx, MMiEz, tags%face%y)
               call fillPECFaceInsideVolume(i,j,k,FACE_Z, MMiHz, MMiEy, MMiEx, tags%face%z)
            end do
         end do
      end do

   contains
      subroutine fillBoundaryFaceIfAllEdgesPEC(i,j,k,face)
         integer(kind=4), intent(in) :: i, j, k, face
         integer(kind=4) :: m1, m2, m3, m4
         logical :: on_boundary = .false.
         select case(face)
         case(FACE_X)
            m1 = MMiEy(i,j,k)
            m2 = MMiEz(i,j,k)
            m3 = MMiEy(i,j,k+1)
            m4 = MMiEz(i,j+1,k)
         case(FACE_Y)
            m1 = MMiEx(i,j,k)
            m2 = MMiEz(i,j,k)
            m3 = MMiEx(i,j,k+1)
            m4 = MMiEz(i+1,j,k)
         case(FACE_Z)
            m1 = MMiEy(i,j,k)
            m2 = MMiEx(i,j,k)
            m3 = MMiEy(i+1,j,k)
            m4 = MMiEx(i,j+1,k)
         end select
         on_boundary = (med(m1)%Is%PEC) .and.(med(m2)%Is%PEC) .and. &
                       (med(m3)%Is%PEC) .and.(med(m4)%Is%PEC)
         if (on_boundary) then 
            Mtag(i,j,k) = 64*numertag 
            select case(face)
            case(FACE_X)
               MMiHx (i, j, k) = indicemedio
               tags%face%x(i,j,k) = 64*numertag
            case(FACE_Y)
               MMiHy (i, j, k) = indicemedio
               tags%face%y(i,j,k) = 64*numertag
            case(FACE_Z)
               MMiHz (i, j, k) = indicemedio
               tags%face%z(i,j,k) = 64*numertag
            end select
         end if
      end subroutine

      logical function hasCrossedPEC(m1,m2,m3,m4)
         integer (kind=INTEGERSIZEOFMEDIAMATRICES) :: m1,m2,m3,m4
         hasCrossedPEC = (med(m1)%Is%ConformalPEC .or. med(m1)%Is%PEC) .and. &
                         (med(m2)%Is%ConformalPEC .or. med(m2)%Is%PEC) .and. &
                         (med(m3)%Is%ConformalPEC .or. med(m3)%Is%PEC) .and. &
                         (med(m4)%Is%ConformalPEC .or. med(m4)%Is%PEC)
      end function 
      logical function hasCrossedPECOrConformalPEC(m1,m2,m3,m4)
         integer (kind=INTEGERSIZEOFMEDIAMATRICES) :: m1,m2,m3,m4
         hasCrossedPECOrConformalPEC = (med(m1)%Is%ConformalPEC .or. med(m1)%Is%PEC) .or. &
                                       (med(m2)%Is%ConformalPEC .or. med(m2)%Is%PEC) .or. &
                                       (med(m3)%Is%ConformalPEC .or. med(m3)%Is%PEC) .or. &
                                       (med(m4)%Is%ConformalPEC .or. med(m4)%Is%PEC)
      end function 

      subroutine fillPECFaceInsideVolume(i, j, k, face, MH, ME1, ME2, tags)
         integer (kind=4), intent(in) :: i, j, k, face
         integer (kind=INTEGERSIZEOFMEDIAMATRICES), dimension(:,:,:), intent(inout) :: MH
         integer (kind=INTEGERSIZEOFMEDIAMATRICES), dimension(:,:,:), intent(inout) :: ME1
         integer (kind=INTEGERSIZEOFMEDIAMATRICES), dimension(:,:,:), intent(inout) :: ME2
         integer (KIND=IKINDMTAG) , dimension(:,:,:), intent(inout) :: tags
         integer (kind=4) :: m1,m2,m3,m4, m
         logical :: on_boundary
         m1 = ME1 (i, j, k)
         m2 = ME2 (i, j, k)
         m3 = ME1 (i + merge(1,0, face==FACE_Z), j, k + merge(1,0, face /= FACE_Z))
         m4 = ME2 (i + merge(1,0, face==FACE_Y), j + merge(1,0, face /= FACE_Y), k)
         m = MH(i,j,k)
         on_boundary = (med(m1)%Is%PEC .or. med(m1)%is%conformalPEC) .and. &
                       (med(m2)%Is%PEC .or. med(m2)%is%conformalPEC) .and. &
                       (med(m3)%Is%PEC .or. med(m3)%is%conformalPEC) .and. &
                       (med(m4)%Is%PEC .or. med(m4)%is%conformalPEC)

         if (on_boundary .and. .not. (med(m)%Is%PEC .or. med(m)%Is%ConformalPEC)) then 
            MH (i, j, k) = indicemedio
            Mtag(i,j,k)=64*numertag 
            tags(i,j,k) = 64*numertag
         else if (on_boundary .and. (med(m)%Is%PEC .or. med(m)%Is%ConformalPEC)) then 
            Mtag(i,j,k)=64*numertag
            tags(i,j,k) = 64*numertag
         end if

      end subroutine


      subroutine fillEdgesInsideVolumeX(j, k)
         integer(kind=4), intent(in) :: j,k
         integer(kind=4) :: i, ii
         logical :: crossed, inside_volume
         integer(kind=4) :: mE, mEPrev, n_crosses
         integer(kind=4), dimension(:), allocatable :: idx_in, idx_out
         inside_volume = .false.
         crossed = .false.
         n_crosses = countCrossesX(j,k)
         allocate(idx_in(n_crosses/2))
         allocate(idx_out(n_crosses/2))
         idx_in(:) = 0
         idx_out(:) = 0
         do i = BoundingBox%xi, BoundingBox%xe+1
            mE = MMiEx(i,j,k)
            mEPrev = MMiEx(i-1,j,k)
            if (.not. (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC)) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                     crossed = hasCrossedPEC(MMiHx(i,j,k), MMiHx(i,j-1,k), MMiHx(i,j,k-1), MMiHx(i,j-1,k-1))
               else if (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC) then 
                     crossed = hasCrossedPECOrConformalPEC(MMiHx(i,j,k),MMiHx(i,j-1,k), MMiHx(i,j,k-1), MMiHx(i,j-1,k-1))
                     crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHx(i-1,j,k),MMiHx(i-1,j-1,k), MMiHx(i-1,j,k-1), MMiHx(i-1,j-1,k-1))
               end if
            else if (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                     crossed = hasCrossedPECOrConformalPEC(MMiHx(i,j,k),MMiHx(i,j-1,k), MMiHx(i,j,k-1), MMiHx(i,j-1,k-1))
                     crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHx(i+1,j,k),MMiHx(i+1,j-1,k), MMiHx(i+1,j,k-1), MMiHx(i+1,j-1,k-1))
               end if
            end if
            if (crossed) inside_volume = .not. inside_volume
            if (crossed .and. inside_volume) idx_in = i
            if (crossed .and. .not. inside_volume) idx_out = i-1
         end do
         do ii = 1, size(idx_in)
            if (idx_in(ii) /= 0 .and. idx_out(ii) /=0) then 
               do i = idx_in(ii), idx_out(ii)-1
                     if (MMiEx (i, j, k) == 1) then 
                        MMiEx (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag 
                        tags%edge%x(i,j,k) = 64*numertag
                     end if
               end do
            end if
         end do

      end subroutine

      function countCrossesX(j,k) result(res)
         integer(kind=4), intent(in) :: j,k
         integer(kind=4) :: res
         integer(kind=4) :: i, mE, mEPrev
         logical :: crossed = .false.
         res = 0
         do i = BoundingBox%xi, BoundingBox%xe+1
            mE = MMiEx(i,j,k)
            mEPrev = MMiEx(i-1,j,k)
            crossed = .false.
            if (.not. (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC)) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                     crossed = hasCrossedPEC(MMiHx(i,j,k), MMiHx(i,j-1,k), MMiHx(i,j,k-1), MMiHx(i,j-1,k-1))
               else if (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC) then 
                     crossed = hasCrossedPECOrConformalPEC(MMiHx(i,j,k),MMiHx(i,j-1,k), MMiHx(i,j,k-1), MMiHx(i,j-1,k-1))
                     crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHx(i-1,j,k),MMiHx(i-1,j-1,k), MMiHx(i-1,j,k-1), MMiHx(i-1,j-1,k-1))
               end if
            else if (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                     crossed = hasCrossedPECOrConformalPEC(MMiHx(i,j,k),MMiHx(i,j-1,k), MMiHx(i,j,k-1), MMiHx(i,j-1,k-1))
                     crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHx(i+1,j,k),MMiHx(i+1,j-1,k), MMiHx(i+1,j,k-1), MMiHx(i+1,j-1,k-1))
               end if
            end if
            if (crossed) res = res + 1
         end do
         if (res /= 0) then 
            if (modulo(res,2) /= 0) error stop 'uneven number of crosses'
         end if
      end function

      subroutine fillEdgesInsideVolumeY(i, k)
         integer(kind=4), intent(in) :: i,k
         integer(kind=4) :: j, jj
         logical :: crossed, inside_volume
         integer(kind=4) :: mE, mEPrev, n_crosses
         integer(kind=4), dimension(:), allocatable :: idx_in, idx_out
         inside_volume = .false.
         crossed = .false.
         n_crosses = countCrossesY(i,k)
         allocate(idx_in(n_crosses/2))
         allocate(idx_out(n_crosses/2))
         idx_in(:) = 0
         idx_out(:) = 0
         do j = BoundingBox%yi, BoundingBox%ye+1
            ! crossing PEC boundary
            mE = MMiEy(i,j,k)
            mEPrev = MMiEy(i,j-1,k)
            crossed = .false.

            if (.not. (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC)) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPEC(MMiHy(i,j,k), MMiHy(i-1,j,k), MMiHy(i,j,k-1), MMiHy(i-1,j,k-1))
               else if (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHy(i,j,k),MMiHy(i,j,k-1),MMiHy(i-1,j,k),MMiHy(i-1,j,k-1))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHy(i,j-1,k), MMiHy(i,j-1,k-1),MMiHy(i-1,j-1,k),MMiHy(i-1,j-1,k-1))
               end if
            else if (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHy(i,j,k),MMiHy(i,j,k-1),MMiHy(i-1,j,k),MMiHy(i-1,j,k-1))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHy(i,j+1,k), MMiHy(i,j+1,k-1),MMiHy(i-1,j+1,k),MMiHy(i-1,j+1,k-1))
               end if
            end if

            if (crossed) inside_volume = .not. inside_volume
            if (crossed .and. inside_volume) idx_in = j
            if (crossed .and. .not. inside_volume) idx_out = j-1
            
         end do
         do jj = 1, size(idx_in)
            if (idx_in(jj) /= 0 .and. idx_out(jj) /=0) then 
               do j = idx_in(jj), idx_out(jj)-1
                  if (MMiEy (i, j, k) == 1) then 
                     MMiEy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%y(i,j,k) = 64*numertag
                  end if
               end do
            end if
         end do
      end subroutine

      function countCrossesY(i,k) result(res)
         integer(kind=4), intent(in) :: i,k
         integer(kind=4) :: res
         integer(kind=4) :: j, mE, mEPrev
         logical :: crossed = .false.
         res = 0
         do j = BoundingBox%yi, BoundingBox%ye+1
            ! crossing PEC boundary
            mE = MMiEy(i,j,k)
            mEPrev = MMiEy(i,j-1,k)
            crossed = .false.
            if (.not. (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC)) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPEC(MMiHy(i,j,k), MMiHy(i-1,j,k), MMiHy(i,j,k-1), MMiHy(i-1,j,k-1))
               else if (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHy(i,j,k),MMiHy(i,j,k-1),MMiHy(i-1,j,k),MMiHy(i-1,j,k-1))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHy(i,j-1,k), MMiHy(i,j-1,k-1),MMiHy(i-1,j-1,k),MMiHy(i-1,j-1,k-1))
               end if
            else if (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHy(i,j,k),MMiHy(i,j,k-1),MMiHy(i-1,j,k),MMiHy(i-1,j,k-1))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHy(i,j+1,k), MMiHy(i,j+1,k-1),MMiHy(i-1,j+1,k),MMiHy(i-1,j+1,k-1))
               end if
            end if

            if (crossed) res = res + 1
         end do
         if (res /= 0) then 
            if (modulo(res,2) /= 0) then 
               error stop 'uneven number of crosses'
            end if
         end if
      end function


      subroutine fillEdgesInsideVolumeZ(i,j)
         integer(kind=4), intent(in) :: i,j
         integer(kind=4) :: k, kk
         logical :: crossed, inside_volume
         integer(kind=4) :: mE, mEPrev, n_crosses
         integer(kind=4), dimension(:), allocatable :: idx_in, idx_out
         inside_volume = .false.
         crossed = .false.
         n_crosses = countCrossesZ(i,j)
         allocate(idx_in(n_crosses/2))
         allocate(idx_out(n_crosses/2))
         idx_in(:) = 0
         idx_out(:) = 0
         do k = BoundingBox%zi, BoundingBox%ze+1
            ! crossing PEC boundary
            mE = MMiEz(i,j,k)
            mEPrev = MMiEz(i,j,k-1)
            if (.not. (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC)) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPEC(MMiHz(i,j,k), MMiHz(i-1,j,k),MMiHz(i,j-1,k),MMiHz(i-1,j-1,k))
               else if (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHz(i,j,k),MMiHz(i-1,j,k),MMiHz(i,j-1,k),MMiHz(i-1,j-1,k))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHz(i,j,k-1),MMiHz(i-1,j,k-1),MMiHz(i,j-1,k-1),MMiHz(i-1,j-1,k-1))
               end if
            else if (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHz(i,j,k),MMiHz(i-1,j,k),MMiHz(i,j-1,k),MMiHz(i-1,j-1,k))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHz(i,j,k+1),MMiHz(i-1,j,k+1),MMiHz(i,j-1,k+1),MMiHz(i-1,j-1,k+1))
               end if
            end if

            if (crossed) inside_volume = .not. inside_volume
            if (crossed .and. inside_volume) idx_in = k
            if (crossed .and. .not. inside_volume) idx_out = k-1
            
         end do
         do kk = 1, size(idx_in)
            if (idx_in(kk) /= 0 .and. idx_out(kk) /=0) then 
               do k = idx_in(kk), idx_out(kk)-1
                  if (MMiEz (i, j, k) == 1) then 
                     MMiEz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%z(i,j,k) = 64*numertag
                  end if
               end do
            end if
         end do
      end subroutine

      function countCrossesZ(i,j) result(res)
         integer(kind=4), intent(in) :: i,j
         integer(kind=4) :: res
         integer(kind=4) :: k, mE, mEPrev
         logical :: crossed = .false.
         res = 0
         do k = BoundingBox%zi, BoundingBox%ze+1
            ! crossing PEC boundary
            mE = MMiEz(i,j,k)
            mEPrev = MMiEz(i,j,k-1)
            crossed = .false.

            if (.not. (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC)) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPEC(MMiHz(i,j,k), MMiHz(i-1,j,k),MMiHz(i,j-1,k),MMiHz(i-1,j-1,k))
               else if (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHz(i,j,k),MMiHz(i-1,j,k),MMiHz(i,j-1,k),MMiHz(i-1,j-1,k))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHz(i,j,k-1),MMiHz(i-1,j,k-1),MMiHz(i,j-1,k-1),MMiHz(i-1,j-1,k-1))
               end if
            else if (med(mE)%Is%ConformalPEC .or. med(mE)%Is%PEC) then 
               if (.not. (med(mEPrev)%Is%ConformalPEC .or. med(mEPrev)%Is%PEC)) then 
                  crossed = hasCrossedPECOrConformalPEC(MMiHz(i,j,k),MMiHz(i-1,j,k),MMiHz(i,j-1,k),MMiHz(i-1,j-1,k))
                  crossed = crossed .or. hasCrossedPECOrConformalPEC(MMiHz(i,j,k+1),MMiHz(i-1,j,k+1),MMiHz(i,j-1,k+1),MMiHz(i-1,j-1,k+1))
               end if
            end if

            if (crossed) res = res + 1
         end do
         if (res /= 0) then 
            if (modulo(res,2) /= 0) error stop 'uneven number of crosses'
         end if
      end function


   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateVolumeMM :  Sets every field component of a volume voxel to the index of the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each voxel
   !                                        centered at i,j,k (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE CreateVolumeMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, BoundingBox, point, indicemedio)
      character(len=BUFSIZE) :: buff
      TYPE (Shared_t) :: Eshared
      !
      INTEGER (KIND=4) :: NumMedia
      TYPE (MediaData_t), DIMENSION (0:NumMedia) :: med
      INTEGER (KIND=4) :: medio
      !
      TYPE (XYZlimit_t) :: punto, puntoPlus1
      TYPE (XYZlimit_t), INTENT (INOUT) :: point
      TYPE (XYZlimit_t), INTENT(IN) ::  BoundingBox
      !
      INTEGER (KIND=4) :: indicemedio
      !
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      INTEGER (KIND=IKINDMTAG) numertag
      INTEGER (KIND=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      !
      INTEGER (KIND=4) :: layoutnumber, i, j, k
      !
      med(indicemedio)%Is%Volume = .TRUE.
      !
      call SortInitEndWithIncreasingOrder(point)
      !
      punto%XI = Max (point%XI, Min(BoundingBox%XI, BoundingBox%XE))
      punto%YI = Max (point%YI, Min(BoundingBox%YI, BoundingBox%YE))
      punto%ZI = Max (point%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
      !
      punto%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE)-1)
      punto%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE)-1)
      punto%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
      !
      puntoPlus1%XE = Min (point%XE+1, Max(BoundingBox%XI, BoundingBox%XE))
      puntoPlus1%YE = Min (point%YE+1, Max(BoundingBox%YI, BoundingBox%YE))
      puntoPlus1%ZE = Min (point%ZE+1, Max(BoundingBox%ZI, BoundingBox%ZE))
      !!!only for volumes the centroid is assigned  !eliminado 03/07/15
      !!      DO k = punto%ZI, punto%ZE
      !!        DO j = punto%YI, punto%YE
      !!          DO i = punto%XI, punto%XE
      !!            medio = MMcen (i, j, k)
      !!            IF (med(indicemedio)%Priority >= med(medio)%Priority) THEN
      !!              MMcen (i, j, k) = indicemedio
      !!            END IF
      !!          END DO
      !!        END DO
      !!      END DO
      !only take care of the boundaries for interfacing
      DO k = punto%ZI, puntoPlus1%ZE
         DO j = punto%YI, puntoPlus1%YE
            DO i = punto%XI, punto%XE
               medio = MMiEx (i, j, k)
!               IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0); 
                     !ojo no es sumar porque no debe desbordarse. solo hay que poner el bit
                     !solo se pone un tag si ya se habia puesto o si estaba si inicializar (cero) y se shifte 6 bits numerados del 0 (iex) al 5 (ihz) empezando por la derecha (lsb)
                     !lo ponto siempre a .true. y quien llege se lo lleva machacando al que habia. los tags no solucionan el problema de determinar univocamente el medio de una celda
                     !a menos que se definan 6 matrices de tags. esto es un niapa que solo sirve para filtrar celdas incluyendo fallos en celdas compartidas entre medios pero no es fiable para determinar medios en posiciones 161020
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     !no lo detectare en volumenes porque podria llevar tiempos elevados en el preproceso
                     !cuando se actualiza el numero de shared (sept'11)
                     !        OnSurface = (k == punto%ZI).or.(k == puntoPlus1%ZE).or.(j == punto%YI).or.(j == puntoPlus1%YE)
                     !        if (OnSurface) call AddToShared(iEx,i,j,k,indicemedio,medio,Eshared)
                  END IF
!               END IF
            END DO
         END DO
      END DO
      !
      DO k = punto%ZI, puntoPlus1%ZE
         DO j = punto%YI, punto%YE
            DO i = punto%XI, puntoPlus1%XE
               medio = MMiEy (i, j, k)
!               IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     !no lo detectare en volumenes porque podria llevar tiempos elevados en el preproceso
                     !cuando se actualiza el numero de shared (sept'11)
                     !        OnSurface = (k == punto%ZI).or.(k == puntoPlus1%ZE).or.(i == punto%XI).or.(i == puntoPlus1%XE)
                     !        if (OnSurface) call AddToShared(iEy,i,j,k,indicemedio,medio,Eshared)
                  END IF
                  
 !              END IF
            END DO
         END DO
      END DO
      !
      DO k = punto%ZI, punto%ZE
         DO j = punto%YI, puntoPlus1%YE
            DO i = punto%XI, puntoPlus1%XE
               medio = MMiEz (i, j, k)
!               IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     !no lo detectare en volumenes porque podria llevar tiempos elevados en el preproceso
                     !cuando se actualiza el numero de shared (sept'11)
                     !        OnSurface = (i == punto%XI).or.(i == puntoPlus1%XE).or.(j == punto%YI).or.(j == puntoPlus1%YE)
                     !        if (OnSurface) call AddToShared(iEz,i,j,k,indicemedio,medio,Eshared)
                  END IF
!               END IF
            END DO
         END DO
      END DO
      !
      DO k = punto%ZI, punto%ZE
         DO j = punto%YI, punto%YE
            DO i = punto%XI, puntoPlus1%XE
               medio = MMiHx (i, j, k)
!               IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  END IF
!               END IF
            END DO
         END DO
      END DO
      !
      DO k = punto%ZI, punto%ZE
         DO j = punto%YI, puntoPlus1%YE
            DO i = punto%XI, punto%XE
               medio = MMiHy (i, j, k)
!               IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  END IF
 !              END IF
            END DO
         END DO
      END DO
      !
      DO k = punto%ZI, puntoPlus1%ZE
         DO j = punto%YI, punto%YE
            DO i = punto%XI, punto%XE
               medio = MMiHz (i, j, k)
!               IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  END IF
!               END IF
            END DO
         END DO
      END DO
      !
      RETURN
   END SUBROUTINE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateSurfaceMM :  Sets every field component of the lower/back/left surface of a voxel to the index of the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   !          orientacion       : Plane of the surface affected by this medium (iEx,iEy,iEz)
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each voxel centered at i,j,k
   !                                        (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE CreateSurfaceMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  &
   & Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, &
   & Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, &
   & Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, &
   & Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, &
   & Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, &
   & med, NumMedia, Eshared, BoundingBox, point, orientacion, indicemedio)
      character(len=BUFSIZE) :: buff
      INTEGER (KIND=4) :: NumMedia
      TYPE (Shared_t) :: Eshared
      TYPE (MediaData_t), DIMENSION (0:NumMedia) :: med
      !
      TYPE (XYZlimit_t) :: punto, puntoPlus1,puntoBboxplus1
      TYPE (XYZlimit_t), INTENT (INOUT) :: point
      TYPE (XYZlimit_t), INTENT(IN) :: BoundingBox
      !
      INTEGER (KIND=4) :: indicemedio, orientacion
      INTEGER (KIND=4) :: layoutnumber, i, j, k
      INTEGER (KIND=4) :: medio
      !
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI , Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      INTEGER (KIND=IKINDMTAG) numertag
      INTEGER (KIND=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      med(indicemedio)%Is%Surface = .TRUE.

      call SortInitEndWithIncreasingOrder(point)
      !
      punto%XI = Max (point%XI, Min(BoundingBox%XI, BoundingBox%XE))
      punto%YI = Max (point%YI, Min(BoundingBox%YI, BoundingBox%YE))
      punto%ZI = Max (point%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
      !
      punto%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE)-1)
      punto%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE)-1)
      punto%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
      !sgg jun'12 para bug en deteccion medios anisotropos en MPI en flushextrainfo
      puntoBboxplus1%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE))
      puntoBboxplus1%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE))
      puntoBboxplus1%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE))
      !
      puntoPlus1%XE = Min (point%XE+1, Max(BoundingBox%XI, BoundingBox%XE))
      puntoPlus1%YE = Min (point%YE+1, Max(BoundingBox%YI, BoundingBox%YE))
      puntoPlus1%ZE = Min (point%ZE+1, Max(BoundingBox%ZI, BoundingBox%ZE))
      !
      SELECT CASE (Abs(orientacion))
       CASE (iEx)
         !    i=punto%XI
         !    if ((i <= max(BoundingBox%XI,BoundingBox%XE)).and.(i >= min(BoundingBox%XI,BoundingBox%XE))) then
         DO i = punto%XI, puntoBboxplus1%XE
            DO j = punto%YI, punto%YE
               DO k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEy (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEy (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     tags%edge%y(i,j,k) = 64*numertag
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     CALL AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                  END IF
               END DO
            END DO
            DO j = punto%YI, puntoPlus1%YE
               DO k = punto%ZI, punto%ZE
                  medio = MMiEz (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEz (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     tags%edge%z(i,j,k) = 64*numertag
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     CALL AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                  END IF
               END DO
            END DO  
            DO j = punto%YI, punto%YE
               DO k = punto%ZI, punto%ZE
                  medio = MMiHx (i, j, k)
!                  IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     IF (med(indicemedio)%Priority > med(medio)%Priority) then
                         MMiHx (i, j, k) = indicemedio; 
                         Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                         tags%face%x(i,j,k) = 64*numertag
                     endif
!                  END IF
               END DO
            END DO
         END DO
         !    endif
       CASE (iEy)
         !    j=punto%YI
         !    if ((j <= max(BoundingBox%YI,BoundingBox%YE)).and.(j >= min(BoundingBox%YI,BoundingBox%YE))) then
         DO j = punto%YI, puntoBboxplus1%YE
            DO i = punto%XI, puntoPlus1%XE
               DO k = punto%ZI, punto%ZE
                  medio = MMiEz (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEz (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     tags%edge%z(i,j,k) = 64*numertag
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     CALL AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                  END IF
               END DO
            END DO
            DO i = punto%XI, punto%XE
               DO k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEx (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEx (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     tags%edge%x(i,j,k) = 64*numertag
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     CALL AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                  END IF
               END DO
            END DO
            DO i = punto%XI, punto%XE
               DO k = punto%ZI, punto%ZE
                  medio = MMiHy (i, j, k)
!                  IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     IF (med(indicemedio)%Priority > med(medio)%Priority) then
                         MMiHy (i, j, k) = indicemedio; 
                         Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);;
                         tags%face%y(i,j,k) = 64*numertag
                     endif
!                  END IF
               END DO
            END DO
         END DO
         !    endif
       CASE (iEz)
         !    k=punto%ZI
         !    if ((k <= max(BoundingBox%ZI,BoundingBox%ZE)).and.(k >= min(BoundingBox%ZI,BoundingBox%ZE))) then
         DO k = punto%ZI, puntoBboxplus1%ZE
            DO i = punto%XI, punto%XE
               DO j = punto%YI, puntoPlus1%YE
                  medio = MMiEx (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEx (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     tags%edge%x(i,j,k) = 64*numertag
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     CALL AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                  END IF
               END DO
            END DO
            DO i = punto%XI, puntoPlus1%XE
               DO j = punto%YI, punto%YE
                  medio = MMiEy (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiEy (i, j, k) = indicemedio; Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     tags%edge%y(i,j,k) = 64*numertag
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     CALL AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                  END IF
               END DO
            END DO
            DO i = punto%XI, punto%XE
               DO j = punto%YI, punto%YE
                  medio = MMiHz (i, j, k)
!                  IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     IF (med(indicemedio)%Priority > med(medio)%Priority) then
                         MMiHz (i, j, k) = indicemedio; 
                         Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                         tags%face%z(i,j,k) = 64*numertag

                     endif
!                  END IF
               END DO
            END DO
         END DO
         !    endif
      END SELECT
      !
      RETURN
   END SUBROUTINE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateLineMM :  Sets every field component of the inner Y/Y/Z axis of a voxel to the index of the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   !          orientacion       : Axis of the voxel affected by this medium (iEx,iEy,iEz)
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each
   !                                        voxel centered at i,j,k (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE CreateLineMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, BoundingBox, point, orientacion, indicemedio, isathinwire, verbose,numeroasignaciones)
      
      TYPE (Shared_t) :: Eshared
      INTEGER (KIND=4) :: NumMedia
      TYPE (MediaData_t), DIMENSION (0:NumMedia) :: med
      !
      TYPE (XYZlimit_t) :: punto
      TYPE (XYZlimit_t), INTENT (INOUT) :: point
      TYPE (XYZlimit_t), INTENT(IN) :: BoundingBox
      
      INTEGER (KIND=4) :: indicemedio, orientacion,numeroasignaciones
      LOGICAL, INTENT (IN) :: isathinwire, verbose
      INTEGER (KIND=4) :: i, j, k, layoutnumber
      INTEGER (KIND=4) :: medio
      !
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      INTEGER (KIND=IKINDMTAG) numertag
      INTEGER (KIND=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      !
      CHARACTER (LEN=BUFSIZE) :: buff
      med(indicemedio)%Is%Line = .TRUE.
      !
      !
      call SortInitEndWithIncreasingOrder(point)
      !
      punto%XI = Max (point%XI, Min(BoundingBox%XI, BoundingBox%XE))
      punto%YI = Max (point%YI, Min(BoundingBox%YI, BoundingBox%YE))
      punto%ZI = Max (point%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
      !
      punto%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE)-1)
      punto%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE)-1)
      punto%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
      !
      SELECT CASE (Abs(orientacion))
       CASE (iEx)
         !    j=punto%YI
         !    k=punto%ZI
         !    if ((j <= max(BoundingBox%YI,BoundingBox%YE)).and.(j >= min(BoundingBox%YI,BoundingBox%YE)).and. &
         !        (k <= max(BoundingBox%ZI,BoundingBox%ZE)).and.(k >= min(BoundingBox%ZI,BoundingBox%ZE))) then
         DO k = punto%ZI, punto%ZE
            DO j = punto%YI, punto%YE
               DO i = punto%XI, punto%XE
                  medio = MMiEx (i, j, k)
!                  IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        numeroasignaciones=numeroasignaciones+1
                        if (med(indicemedio)%is%lumped) then
                            if (numeroasignaciones==1) then !solo le echa el lumped a 1 segmento !esto es una peticion externa !ojo es agresivo. !solo se pone 1 segmento con la resistencia especificada. me doy cuenta en 040123
                                MMiEx (i, j, k) = indicemedio
                                Mtag(i,j,k)=64*numertag 
                                tags%edge%x(i,j,k) = 64*numertag
                                ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                            else
                                MMiEx (i, j, k) = 0
                                Mtag(i,j,k)=64*numertag 
                                tags%edge%x(i,j,k) = 64*numertag
                                ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);;
                            endif
                        else
                            MMiEx (i, j, k) = indicemedio; 
                            Mtag(i,j,k)=64*numertag 
                            tags%edge%x(i,j,k) = 64*numertag
                            ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                        endif
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        CALL AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                     END IF
               END DO
            END DO
         END DO
         !    endif
       CASE (iEy)
         !    i=punto%XI
         !    k=punto%ZI
         !    if ((i <= max(BoundingBox%XI,BoundingBox%XE)).and.(i >= min(BoundingBox%XI,BoundingBox%XE)).and. &
         !        (k <= max(BoundingBox%ZI,BoundingBox%ZE)).and.(k >= min(BoundingBox%ZI,BoundingBox%ZE))) then
         DO k = punto%ZI, punto%ZE
            DO j = punto%YI, punto%YE
               DO i = punto%XI, punto%XE
                  medio = MMiEy (i, j, k)
!                  IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        numeroasignaciones=numeroasignaciones+1
                        if (med(indicemedio)%is%lumped) then
                            if (numeroasignaciones==1) then !solo le echa el lumped a 1 segmento
                                MMiEy (i, j, k) = indicemedio
                                Mtag(i,j,k)=64*numertag 
                                tags%edge%y(i,j,k) = 64*numertag
                                ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                            else
                                MMiEy (i, j, k)  = 0 
                                Mtag(i,j,k)=64*numertag 
                                tags%edge%y(i,j,k) = 64*numertag
                                ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                            endif
                        else
                            MMiEy (i, j, k) = indicemedio
                            Mtag(i,j,k)=64*numertag
                            tags%edge%y(i,j,k) = 64*numertag
                            ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                        endif
                        
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        CALL AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                     END IF
               END DO
            END DO
         END DO
         !    endif
       CASE (iEz)
         !    i=punto%XI
         !    j=punto%YI
         !    if ((i <= max(BoundingBox%XI,BoundingBox%XE)).and.(i >= min(BoundingBox%XI,BoundingBox%XE)).and. &
         !        (j <= max(BoundingBox%YI,BoundingBox%YE)).and.(j >= min(BoundingBox%YI,BoundingBox%YE))) then
         DO k = punto%ZI, punto%ZE
            DO j = punto%YI, punto%YE
               DO i = punto%XI, punto%XE
                  medio = MMiEz (i, j, k)
!                  IF (medio /= 0) THEN   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        numeroasignaciones=numeroasignaciones+1
                        if (med(indicemedio)%is%lumped) then
                            if (numeroasignaciones==1) then !solo le echa el lumped a 1 segmento
                                MMiEz (i, j, k) = indicemedio
                                Mtag(i,j,k)=64*numertag
                                tags%edge%z(i,j,k) = 64*numertag
                                ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                            else
                                MMiEz (i, j, k) = 0
                                Mtag(i,j,k)=64*numertag
                                tags%edge%z(i,j,k) = 64*numertag
                                ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);;
                            endif
                        else
                            MMiEz (i, j, k) = indicemedio
                            Mtag(i,j,k)=64*numertag
                            tags%edge%z(i,j,k) = 64*numertag
                            ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                        endif
                        
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        CALL AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                     END IF

               END DO
            END DO
         END DO
         !    endif
      END SELECT
      !
      RETURN
   END SUBROUTINE
   !Slot=special case of surface.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateSurfaceSlotMM :  Sets every field component of the lower/back/left surface of a voxel to the index of
   !                                    the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   !          orientacion       : Plane of the surface affected by this medium (iEx,iEy,iEz)
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each voxel centered at i,j,k
   !                                        (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE CreateSurfaceSlotMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, Hshared, BoundingBox, point, orientacion, direccion, indicemedio)
      character(len=BUFSIZE) :: buff
      TYPE (Shared_t) :: Eshared, Hshared
      INTEGER (KIND=4) :: NumMedia
      TYPE (MediaData_t), DIMENSION (0:NumMedia) :: med
      !
      TYPE (XYZlimit_t) :: punto, puntoPlus1,puntoBboxplus1
      TYPE (XYZlimit_t), INTENT (INOUT) :: point
      TYPE (XYZlimit_t), INTENT(IN) :: BoundingBox
      !
      INTEGER (KIND=4) :: indicemedio, orientacion, direccion
      !
      INTEGER (KIND=4) :: layoutnumber, i, j, k, offx, offy, offz
      INTEGER (KIND=4) :: medio
      !
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      INTEGER (KIND=IKINDMTAG) numertag
      INTEGER (KIND=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      med(indicemedio)%Is%Surface = .TRUE.
      !
      call SortInitEndWithIncreasingOrder(point)
      !
      !
      punto%XI = Max (point%XI, Min(BoundingBox%XI, BoundingBox%XE))
      punto%YI = Max (point%YI, Min(BoundingBox%YI, BoundingBox%YE))
      punto%ZI = Max (point%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
      !
      punto%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE)-1)
      punto%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE)-1)
      punto%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
      !sgg jun'12 para bug en deteccion medios anisotropos en MPI en flushextrainfo
      puntoBboxplus1%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE))
      puntoBboxplus1%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE))
      puntoBboxplus1%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE))
      !
      puntoPlus1%XE = Min (point%XE+1, Max(BoundingBox%XI, BoundingBox%XE))
      puntoPlus1%YE = Min (point%YE+1, Max(BoundingBox%YI, BoundingBox%YE))
      puntoPlus1%ZE = Min (point%ZE+1, Max(BoundingBox%ZI, BoundingBox%ZE))
      !
      offx = 0
      offy = 0
      offz = 0
      SELECT CASE (Abs(orientacion))
       CASE (iEx)
         DO i = punto%XI, puntoBboxplus1%XE
            SELECT CASE (direccion)
             CASE (iEz)
               offx = 0
               offy = 0
               offz = 1
               DO j = punto%YI, punto%YE
                  DO k = punto%ZI, puntoPlus1%ZE
                     medio = MMiEy (i, j, k)
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        MMiEy (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%y(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        !CALL AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                     END IF
                  END DO
               END DO
             CASE (iEy)
               offx = 0
               offy = 1
               offz = 0
               DO j = punto%YI, puntoPlus1%YE
                  DO k = punto%ZI, punto%ZE
                     medio = MMiEz (i, j, k)
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        MMiEz (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%z(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        !CALL AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                     END IF
                  END DO
               END DO
            END SELECT
            DO j = Max (punto%YI - offy, Min(BoundingBox%YI, BoundingBox%YE)), &
            &       Min (punto%YE + offy, Max(BoundingBox%YI, BoundingBox%YE)-1)
               DO k = Max (punto%ZI - offz, Min(BoundingBox%ZI, BoundingBox%ZE)),  &
               &       Min (punto%ZE + offz, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
                  medio = MMiHx (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     !CALL AddToShared (iHx, i, j, k, indicemedio, medio, Hshared)

                  END IF
               END DO
            END DO
         END DO
       CASE (iEy)
         DO j = punto%YI, puntoBboxplus1%YE
            SELECT CASE (direccion)
             CASE (iEx)
               offx = 1
               offy = 0
               offz = 0
               DO i = punto%XI, puntoPlus1%XE
                  DO k = punto%ZI, punto%ZE
                     medio = MMiEz (i, j, k)
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        MMiEz (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%z(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        !CALL AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                     END IF
                  END DO
               END DO
             CASE (iEz)
               offx = 0
               offy = 0
               offz = 1
               DO i = punto%XI, punto%XE
                  DO k = punto%ZI, puntoPlus1%ZE
                     medio = MMiEx (i, j, k)
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        MMiEx (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%x(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        !CALL AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                     END IF
                  END DO
               END DO
            END SELECT
            DO i = Max (punto%XI - offx, Min(BoundingBox%XI, BoundingBox%XE)),  &
            &       Min (punto%XE + offx, Max(BoundingBox%XI, BoundingBox%XE)-1)
               DO k = Max (punto%ZI - offz, Min(BoundingBox%ZI, BoundingBox%ZE)),  &
               &       Min (punto%ZE + offz, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
                  medio = MMiHy (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     !CALL AddToShared (iHy, i, j, k, indicemedio, medio, Hshared)
                  END IF
               END DO
            END DO
         END DO
       CASE (iEz)
         DO k = punto%ZI, puntoBboxplus1%ZE
            SELECT CASE (direccion)
             CASE (iEy)
               offx = 0
               offy = 1
               offz = 0
               DO i = punto%XI, punto%XE
                  DO j = punto%YI, puntoPlus1%YE
                     medio = MMiEx (i, j, k)
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        MMiEx (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%x(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        !CALL AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                     END IF
                  END DO
               END DO
             CASE (iEx)
               offx = 1
               offy = 0
               offz = 0
               DO i = punto%XI, puntoPlus1%XE
                  DO j = punto%YI, punto%YE
                     medio = MMiEy (i, j, k)
                     IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                        MMiEy (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%y(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                        !CALL AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                     END IF
                  END DO
               END DO
            END SELECT
            DO i = Max (punto%XI - offx, Min(BoundingBox%XI, BoundingBox%XE)),  &
            &       Min (punto%XE + offx, Max(BoundingBox%XI, BoundingBox%XE)-1)
               DO j = Max (punto%YI - offy, Min(BoundingBox%YI, BoundingBox%YE)),  &
               &       Min (punto%YE + offy, Max(BoundingBox%YI, BoundingBox%YE)-1)
                  medio = MMiHz (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  ELSE IF ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) THEN
                     !CALL AddToShared (iHz, i, j, k, indicemedio, medio, Hshared)
                  END IF
               END DO
            END DO
         END DO
      END SELECT
      !
      RETURN
   END SUBROUTINE
   !!!!special case of magneticsurface (for the multiport padding)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateMagneticSurface :  Sets every field component of the lower/back/left surface of a voxel to the index of the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   !          orientacion       : Plane of the surface affected by this medium (iEx,iEy,iEz)
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each voxel centered at i,j,k
   !                                        (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE CreateMagneticSurface(layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, BoundingBox, point, orientacion, indicemedio)
      character(len=BUFSIZE) :: buff
      INTEGER (KIND=4) :: NumMedia
      TYPE (Shared_t) :: Eshared
      TYPE (MediaData_t), DIMENSION (0:NumMedia) :: med
      !
      TYPE (XYZlimit_t) :: punto, puntoPlus1   !,puntoBboxplus1
      TYPE (XYZlimit_t), INTENT (INOUT) :: point
      TYPE (XYZlimit_t), INTENT(IN) :: BoundingBox
      !
      INTEGER (KIND=4) :: indicemedio, orientacion
      INTEGER (KIND=4) :: layoutnumber, i, j, k
      INTEGER (KIND=4) :: medio
      !
      INTEGER (KIND=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      INTEGER (KIND=IKINDMTAG) numertag
      INTEGER (KIND=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      med(indicemedio)%Is%Surface = .TRUE.
      !
      !
      call SortInitEndWithIncreasingOrder(point)
      !
      punto%XI = Max (point%XI, Min(BoundingBox%XI, BoundingBox%XE))
      punto%YI = Max (point%YI, Min(BoundingBox%YI, BoundingBox%YE))
      punto%ZI = Max (point%ZI, Min(BoundingBox%ZI, BoundingBox%ZE))
      !
      punto%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE)-1)
      punto%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE)-1)
      punto%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
      !!sgg jun'12 para bug en deteccion medios anisotropos en MPI en flushextrainfo
      !      puntoBboxplus1%XE = Min (point%XE, Max(BoundingBox%XI, BoundingBox%XE))
      !      puntoBboxplus1%YE = Min (point%YE, Max(BoundingBox%YI, BoundingBox%YE))
      !      puntoBboxplus1%ZE = Min (point%ZE, Max(BoundingBox%ZI, BoundingBox%ZE))
      ! aqui me da problemas lo comento 20jul 12

      puntoPlus1%XE = Min (point%XE+1, Max(BoundingBox%XI, BoundingBox%XE))
      puntoPlus1%YE = Min (point%YE+1, Max(BoundingBox%YI, BoundingBox%YE))
      puntoPlus1%ZE = Min (point%ZE+1, Max(BoundingBox%ZI, BoundingBox%ZE))
      !
      SELECT CASE (Abs(orientacion))
       CASE (iEx)
         DO i = punto%XI, punto%XE
            DO j = punto%YI, puntoPlus1%YE
               DO k = punto%ZI, punto%ZE
                  medio = MMiHy (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  endif
               END DO
            END DO
            DO j = punto%YI, punto%YE
               DO k = punto%ZI, puntoPlus1%ZE
                  medio = MMiHz (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  endif
               END DO
            END DO
            DO j = punto%YI, puntoPlus1%YE
               DO k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEx (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%edge%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                  endif
               END DO
            END DO
         END DO
       CASE (iEy)
         DO j = punto%YI, punto%YE
            DO i = punto%XI, punto%XE
               DO k = punto%ZI, puntoPlus1%ZE
                  medio = MMiHz (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  endif
               END DO
            END DO
            DO i = punto%XI, puntoPlus1%XE
               DO k = punto%ZI, punto%ZE
                  medio = MMiHx (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  END IF
               END DO
            END DO
            DO i = punto%XI, puntoPlus1%XE
               DO k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEy (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%edge%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                  endif
               END DO
            END DO
         END DO
         !
       CASE (iEz)
         DO k = punto%ZI, punto%ZE
            DO i = punto%XI, puntoPlus1%XE
               DO j = punto%YI, punto%YE
                  medio = MMiHx (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  END IF
               END DO
            END DO
            DO i = punto%XI, punto%XE
               DO j = punto%YI, puntoPlus1%YE
                  medio = MMiHy (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) THEN
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  END IF
               END DO
            END DO
            DO i = punto%XI, puntoPlus1%XE
               DO j = punto%YI, puntoPlus1%YE
                  medio = MMiEz (i, j, k)
                  IF (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%edge%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                  endif
               END DO
            END DO
         END DO

      END SELECT
      !
      RETURN
   END SUBROUTINE


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreatePMLmatrix :  ............
   ! Inputs :   O%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          NumMedia = Previous number of media
   !          Med%Epr,Med%Mur,Med%Sigma,Med%SigmaM = Constitutive parameters
   !          Med%Wire,Med%multiport = Med%Wire and Med%multiport info
   !          Border  = type of borders
   !          PML = PML info
   !          BoundingBox%XE,BoundingBox%YE,BoundingBox%ZE
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = new types of medium taking into account interfaces, PML and PECs
   !          NumMedia   = New number of Media
   !          Med%Epr,Med%Mur,Med%Sigma,Med%SigmaM = New Matrices with average and PML constitutive parameters of new media
   !          Med%Wire,Med%multiport = Same as input but resized accordingly to take into account the increment in NumMedia
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE CreatePMLmatrix (layoutnumber, SIZE, sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, SINPML_fullsize, fullsize, BBox, med, NumMedia, Border, MEDIOEXTRA)
      !
      TYPE (limit_t), DIMENSION (1:6) :: SINPML_fullsize, fullsize
      !Inputs and Outputs
      type (SGGFDTDINFO), intent(INOUT)        :: sgg
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES)   ::  &
      sggMiEx(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE), &
      sggMiEy(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE), &
      sggMiEz(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE), &
      sggMiHx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE), &
      sggMiHy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE), &
      sggMiHz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)
      TYPE (XYZlimit_t) :: BoundingBox, BBox
      TYPE (MediaData_t), POINTER, DIMENSION (:) :: med
      INTEGER (KIND=4), INTENT (INOUT) :: NumMedia
      TYPE (Border_t), INTENT (IN) :: Border
      TYPE (MedioExtra_t), INTENT (IN) :: MEDIOEXTRA
      ! Local stuff
      INTEGER (KIND=4), POINTER, DIMENSION (:) :: tempo
      TYPE (MediaData_t), POINTER, DIMENSION (:) :: NewMed
      INTEGER (KIND=4) :: layoutnumber, SIZE, field, medium, i, j, k, NuevoNumeroMediosConPML
      INTEGER (KIND=4) :: oldNumMedia,oldmed
      INTEGER (KIND=4), DIMENSION (1:6) :: XIPML, XEPML, YIPML, YEPML, ZIPML, ZEPML
      REAL (KIND=RKIND) :: oldepr,oldmur,oldsigma,oldsigmam,newepr,newmur,newsigma,newsigmam
      LOGICAL :: yapuesto

      !Increase mediamatrix with PML regions, and remove the minus sign from the mm matrix (this info is no longer needed)
      !FIRST CLIP THE MATRIX
      !readjust boundingbox for PML correct calculation
      !
      DO field = iEx, iHz
         XIPML (field) = Max (BBox%XI, fullsize(iHx)%XI)
         XEPML (field) = Min (BBox%XE+on(field, icoord, fine), fullsize(iHx)%XE)
         YIPML (field) = Max (BBox%YI, fullsize(iHy)%YI)
         YEPML (field) = Min (BBox%YE+on(field, jcoord, fine), fullsize(iHy)%YE)
         ZIPML (field) = Max (BBox%ZI, fullsize(iHz)%ZI)
         ZEPML (field) = Min (BBox%ZE+on(field, kcoord, fine), fullsize(iHz)%ZE)
      END DO
      !
      BoundingBox%XI = Max (BBox%XI, SINPML_fullsize(iHx)%XI)
      BoundingBox%XE = Min (BBox%XE, SINPML_fullsize(iHx)%XE)
      BoundingBox%YI = Max (BBox%YI, SINPML_fullsize(iHy)%YI)
      BoundingBox%YE = Min (BBox%YE, SINPML_fullsize(iHy)%YE)
      BoundingBox%ZI = Max (BBox%ZI, SINPML_fullsize(iHz)%ZI)
      BoundingBox%ZE = Min (BBox%ZE, SINPML_fullsize(iHz)%ZE)
      ! Build the interior of the PML regions in MediaMatrix
      ! temporarily assing minus sign to PML media
      ! corners are swept twice to assing the correct media (do not remove AbS!!)
      field = iEx
      DO j = YIPML (field), YEPML (field)
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            DO i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiEx (i, j, k) = - Abs (sggmiEx(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Front
            DO i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiEx (i, j, k) = - Abs (sggmiEx(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            DO j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Right
            DO j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            DO k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Up
            DO k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      field = iEy
      DO j = YIPML (field), YEPML (field)
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            DO i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiEy (i, j, k) = - Abs (sggmiEy(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Front
            DO i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiEy (i, j, k) = - Abs (sggmiEy(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            DO j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Right
            DO j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            DO k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Up
            DO k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      field = iEz
      DO j = YIPML (field), YEPML (field)
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            DO i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiEz (i, j, k) = - Abs (sggmiEz(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Front
            DO i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiEz (i, j, k) = - Abs (sggmiEz(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            DO j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Right
            DO j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            DO k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Up
            DO k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      field = iHx
      DO j = YIPML (field), YEPML (field)
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            DO i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiHx (i, j, k) = - Abs (sggmiHx(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Front
            DO i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiHx (i, j, k) = - Abs (sggmiHx(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            DO j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Right
            DO j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            DO k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Up
            DO k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      field = iHy
      DO j = YIPML (field), YEPML (field)
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            DO i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiHy (i, j, k) = - Abs (sggmiHy(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Front
            DO i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiHy (i, j, k) = - Abs (sggmiHy(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            DO j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Right
            DO j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            DO k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Up
            DO k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      field = iHz
      DO j = YIPML (field), YEPML (field)
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            DO i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiHz (i, j, k) = - Abs (sggmiHz(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Front
            DO i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiHz (i, j, k) = - Abs (sggmiHz(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            DO j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Right
            DO j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            END DO
         END DO
      END DO
      !
      DO i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         DO j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            DO k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            END DO
            !!!!!!**Up
            DO k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            END DO
         END DO
      END DO
      !compact the info of PML media
      NuevoNumeroMediosConPML = NumMedia
      ALLOCATE (tempo(0:NumMedia))
      tempo = 0 !temporarily stores the index of the PML medium matching each original media
      field = iEx
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiEx (i, j, k)
               IF (medium < 0) THEN
                  IF (tempo(Abs(medium)) == 0) THEN
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  END IF
               END IF
            END DO
         END DO
      END DO
      field = iEy
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiEy (i, j, k)
               IF (medium < 0) THEN
                  IF (tempo(Abs(medium)) == 0) THEN
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  END IF
               END IF
            END DO
         END DO
      END DO
      field = iEz
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiEz (i, j, k)
               IF (medium < 0) THEN
                  IF (tempo(Abs(medium)) == 0) THEN
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  END IF
               END IF
            END DO
         END DO
      END DO
      field = iHx
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiHx (i, j, k)
               IF (medium < 0) THEN
                  IF (tempo(Abs(medium)) == 0) THEN
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  END IF
               END IF
            END DO
         END DO
      END DO
      field = iHy
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiHy (i, j, k)
               IF (medium < 0) THEN
                  IF (tempo(Abs(medium)) == 0) THEN
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  END IF
               END IF
            END DO
         END DO
      END DO
      field = iHz
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiHz (i, j, k)
               IF (medium < 0) THEN
                  IF (tempo(Abs(medium)) == 0) THEN
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  END IF
               END IF
            END DO
         END DO
      END DO
      !
      ALLOCATE (NewMed(NumMedia+1:NuevoNumeroMediosConPML))
      !Reassing the PML media info with the compact indexes
      field = iEx
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiEx (i, j, k)
               IF (medium < 0) THEN
                  sggmiEx (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               END IF
            END DO
         END DO
      END DO
      field = iEy
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiEy (i, j, k)
               IF (medium < 0) THEN
                  sggmiEy (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               END IF
            END DO
         END DO
      END DO
      field = iEz
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiEz (i, j, k)
               IF (medium < 0) THEN
                  sggmiEz (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               END IF
            END DO
         END DO
      END DO
      field = iHx
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiHx (i, j, k)
               IF (medium < 0) THEN
                  sggmiHx (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               END IF
            END DO
         END DO
      END DO
      field = iHy
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiHy (i, j, k)
               IF (medium < 0) THEN
                  sggmiHy (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               END IF
            END DO
         END DO
      END DO
      field = iHz
      DO k = ZIPML (field), ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               medium = sggmiHz (i, j, k)
               IF (medium < 0) THEN
                  sggmiHz (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               END IF
            END DO
         END DO
      END DO

      !Put PEC and the end if there exists PEC borders in the original problem
      !(PMC are handled with the image technique in the algorithm, no special index is used for PMC)
      !DETRAS Y DELANTE
      field = iEx !!!!!PEC solo en fields donde acabe la red
      !izda y dcha
      IF ((Border%IsLeftPEC)) THEN
         j = YIPML (field)
         DO i = XIPML (field), XEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEx (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsRightPEC)) THEN
         j = YEPML (field)
         DO i = XIPML (field), XEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEx (i, j, k) = 0
            END DO
         END DO
      END IF
      !  !Up y Down
      IF ((Border%IsDownPEC)) THEN
         k = ZIPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               sggmiEx (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsUpPEC)) THEN
         k = ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               sggmiEx (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      field = iEy !!!!!PEC solo en fields donde acabe la red
      !front y back
      IF ((Border%IsBackPEC)) THEN
         i = XIPML (field)
         DO j = YIPML (field), YEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEy (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsFrontPEC)) THEN
         i = XEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEy (i, j, k) = 0
            END DO
         END DO
      END IF
      !  !Up y Down
      IF ((Border%IsDownPEC)) THEN
         k = ZIPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               sggmiEy (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsUpPEC)) THEN
         k = ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               sggmiEy (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      field = iEz !!!!!PEC solo en fields donde acabe la red
      !front y back
      IF ((Border%IsBackPEC)) THEN
         i = XIPML (field)
         DO j = YIPML (field), YEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsFrontPEC)) THEN
         i = XEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            END DO
         END DO
      END IF
      !izda y dcha
      IF ((Border%IsLeftPEC)) THEN
         j = YIPML (field)
         DO i = XIPML (field), XEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsRightPEC)) THEN
         j = YEPML (field)
         DO i = XIPML (field), XEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      field = iHx !!!!!PEC solo en fields donde acabe la red
      !front y back
      IF ((Border%IsBackPEC)) THEN
         i = XIPML (field)
         DO j = YIPML (field), YEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiHx (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsFrontPEC)) THEN
         i = XEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiHx (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      field = iHy !!!!!PEC solo en fields donde acabe la red
      !izda y dcha
      IF ((Border%IsLeftPEC)) THEN
         j = YIPML (field)
         DO i = XIPML (field), XEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiHy (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsRightPEC)) THEN
         j = YEPML (field)
         DO i = XIPML (field), XEPML (field)
            DO k = ZIPML (field), ZEPML (field)
               sggmiHy (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      field = iHz !!!!!PEC solo en fields donde acabe la red
      !  !Up y Down
      IF ((Border%IsDownPEC)) THEN
         k = ZIPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               sggmiHz (i, j, k) = 0
            END DO
         END DO
      END IF
      !
      IF ((Border%IsUpPEC)) THEN
         k = ZEPML (field)
         DO j = YIPML (field), YEPML (field)
            DO i = XIPML (field), XEPML (field)
               sggmiHz (i, j, k) = 0
            END DO
         END DO
      END IF

      !!?!?!?MEDIOEXTRA

      !
      !adjust constitutive parameters, matrices
      oldNumMedia = NumMedia !save it since the next subroutine overwrites it
      CALL Readjust (NumMedia, med, NuevoNumeroMediosConPML)
      sgg%AllocMed=NumMedia
      !copy the new media
      med (1+oldNumMedia:NuevoNumeroMediosConPML) = NewMed (1+oldNumMedia:NuevoNumeroMediosConPML)
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%PML = .TRUE. !all these are PML
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%ThinWire = .FALSE. !put any wire touching the PML to non-wire though treat it with mur
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%SlantedWire = .FALSE. !put any wire touching the PML to non-wire though treat it with mur
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%Needed=.true. !sgg 220817 por defecto lo he puesto en readjust a false
      !
      DEALLOCATE (NewMed, tempo)

      !solo lo creo para las tangenciales electricas
      if (MEDIOEXTRA%exists) then
         !Put MEDIO and the end if there exists PML borders in the original problem
         yapuesto=.false.
         !
         field = iEx
         !izda y dcha
         IF ((Border%IsLeftPML)) THEN
            DO j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsRightPML)) THEN
            DO j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !  !Up y Down
         IF ((Border%IsDownPML)) THEN
            DO k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiEx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsUpPML)) THEN
            DO k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiEx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         field = iEy
         !front y back
         IF ((Border%IsBackPML)) THEN
            DO i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsFrontPML)) THEN
            DO i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !  !Up y Down
         IF ((Border%IsDownPML)) THEN
            DO k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiEy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsUpPML)) THEN
            DO k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiEy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         field = iEz
         !front y back
         IF ((Border%IsBackPML)) THEN
            DO i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsFrontPML)) THEN
            DO i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !izda y dcha
         IF ((Border%IsLeftPML)) THEN
            DO j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsRightPML)) THEN
            DO j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiEz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiEz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         !magneticos

         !
         !
         field = iHx
         !izda y dcha
         IF ((Border%IsLeftPML)) THEN
            DO j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsRightPML)) THEN
            DO j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !  !Up y Down
         IF ((Border%IsDownPML)) THEN
            DO k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiHx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsUpPML)) THEN
            DO k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiHx (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHx (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         field = iHy
         !front y back
         IF ((Border%IsBackPML)) THEN
            DO i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsFrontPML)) THEN
            DO i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !  !Up y Down
         IF ((Border%IsDownPML)) THEN
            DO k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiHy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsUpPML)) THEN
            DO k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO i = XIPML (field), XEPML (field)
                     !
                     oldmed   =sggmiHy (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHy (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         field = iHz
         !front y back
         IF ((Border%IsBackPML)) THEN
            DO i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsFrontPML)) THEN
            DO i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               DO j = YIPML (field), YEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !izda y dcha
         IF ((Border%IsLeftPML)) THEN
            DO j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
         IF ((Border%IsRightPML)) THEN
            DO j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               DO i = XIPML (field), XEPML (field)
                  DO k = ZIPML (field), ZEPML (field)
                     !
                     oldmed   =sggmiHz (i, j, k)
                     if (oldmed /= 0) then
                        oldepr   =sgg%Med(oldmed)%epr
                        oldmur   =sgg%Med(oldmed)%mur
                        oldsigma =sgg%Med(oldmed)%sigma
                        oldsigmam=sgg%Med(oldmed)%sigmam
                        !
                        newepr   =sgg%Med(MEDIOEXTRA%index)%epr
                        newmur   =sgg%Med(MEDIOEXTRA%index)%mur
                        newsigma =sgg%Med(MEDIOEXTRA%index)%sigma
                        newsigmam=sgg%Med(MEDIOEXTRA%index)%sigmam
                        if (yapuesto) then
                           if ((oldmed /= MEDIOEXTRA%index).and.((newepr /= oldepr).or.(newmur /= oldmur).or. &
                           (newsigma /= oldsigma  + MEDIOEXTRA%sigma ).or.(newsigmam /= oldsigmam + MEDIOEXTRA%sigmam))) then
                              CALL STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
                           endif
                        else
                           sgg%Med(MEDIOEXTRA%index)%epr    = oldepr
                           sgg%Med(MEDIOEXTRA%index)%mur    = oldmur
                           sgg%Med(MEDIOEXTRA%index)%sigma  = oldsigma + MEDIOEXTRA%sigma
                           sgg%Med(MEDIOEXTRA%index)%sigmam = oldsigmam + MEDIOEXTRA%sigmam
                        endif
                        !
                        sggmiHz (i, j, k) = MEDIOEXTRA%index
                        yapuesto=.true.
                     endif
                  END DO
               END DO
            END DO
         END IF
         !
      endif !del medioextra%exists
      RETURN
   END SUBROUTINE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  Readjust :  ............
   ! Inputs :   NumMedia = Previous number of media
   !          Med%Epr,Med%Mur,Med%Sigma,Med%SigmaM = Constitutive parameters
   !          Med%Wire,Med%multiport = Med%Wire and Med%multiport info
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = new types of medium taking into account interfaces, PML and PECs
   !          NumMedia   = New number of Media
   !          Med%Epr,Med%Mur,Med%Sigma,Med%SigmaM = New Matrices with average and PML constitutive parameters of new media
   !          Med%Wire,Med%multiport = Same as input but resized accordingly to take into account the increment in NumMedia
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE Readjust (NumMedia, med, NewNumMedia)
      !Inputs and Outputs
      TYPE (MediaData_t), POINTER, DIMENSION (:) :: med
      TYPE (MediaData_t), POINTER, DIMENSION (:) :: dummyMed
      INTEGER (KIND=4), INTENT (INOUT) :: NumMedia
      INTEGER (KIND=4), INTENT (IN) :: NewNumMedia
      INTEGER (KIND=4) :: i
      !
      ALLOCATE (dummyMed(0:NewNumMedia))
      DO i = 0, Min (NumMedia, NewNumMedia)
         dummyMed (i) = med (i)
      END DO
      !
      DEALLOCATE (med)
      ALLOCATE (med(0:NewNumMedia))
      DO i = 0, NewNumMedia
         med (i) = dummyMed (i)
      END DO
      !
      DEALLOCATE (dummyMed)
      DO i = 1 + NumMedia, NewNumMedia
         med(i)%Priority = prior_BV !background
         med(i)%epr=-1.0_RKIND 
         med(i)%MUr=-1.0_RKIND 
         med(i)%SIGMA=-1.0_RKIND 
         med(i)%SIGMAM=-1.0_RKIND 
!
         med(i)%Is%PML = .FALSE.
         med(i)%Is%PEC = .FALSE.
         med(i)%Is%PMC = .FALSE.
         med(i)%Is%ThinWire = .FALSE.
         med(i)%Is%SlantedWire = .FALSE.
         med(i)%Is%EDispersive = .FALSE.
         med(i)%Is%MDispersive = .FALSE.
         Med(i)%Is%EDispersiveAnis= .FALSE.
         Med(i)%Is%MDispersiveAnis= .FALSE.
         med(i)%Is%ThinSlot = .FALSE.
         Med(i)%Is%PMLbody= .FALSE.
         Med(i)%Is%SGBC= .FALSE.
         Med(i)%Is%SGBCDispersive= .FALSE.
         Med(i)%Is%Lumped= .FALSE.
         Med(i)%Is%Lossy= .FALSE.
         med(i)%Is%AnisMultiport = .FALSE.
         med(i)%Is%multiport = .FALSE.
         med(i)%Is%multiportPadding = .FALSE.
         med(i)%Is%Dielectric = .FALSE.
         med(i)%Is%Anisotropic = .FALSE.
         med(i)%Is%Volume = .FALSE.
         med(i)%Is%Line = .FALSE.
         med(i)%Is%Surface = .FALSE.
         med(i)%Is%Needed = .FALSE. !!!.TRUE. !sgg 220817 en principio no es needed. quien llame a readjust debe luego poner needed a true segun sea eso cierto
         med(i)%Is%Interfase = .FALSE.
         Med(i)%Is%already_YEEadvanced_byconformal = .FALSE.
         Med(i)%Is%split_and_useless = .FALSE.
      END DO
      !Update the final number of media
      NumMedia = NewNumMedia
      !
   END SUBROUTINE

   SUBROUTINE AddToShared (campo, i1, j1, k1, Sharedmed, ProPmed,  Shared)
      TYPE (SharedElement_t), POINTER, DIMENSION (:) :: temp
      TYPE (Shared_t), INTENT (INOUT) :: Shared
      INTEGER (KIND=4), INTENT (IN) :: campo, i1, j1, k1, Sharedmed, ProPmed
      INTEGER (KIND=4) :: conta, n
      !
      Shared%conta = Shared%conta + 1
      conta = Shared%conta
      !
      IF (conta > Shared%MaxConta) THEN
         !create space on the fly
         ALLOCATE (temp(1:conta-1))
         DO n = 1, conta - 1
            temp(n)%Sharedmed = Shared%elem(n)%Sharedmed
            temp(n)%ProPmed = Shared%elem(n)%ProPmed
            temp(n)%field = Shared%elem(n)%field
            temp(n)%i = Shared%elem(n)%i
            temp(n)%j = Shared%elem(n)%j
            temp(n)%k = Shared%elem(n)%k
            temp(n)%times = Shared%elem(n)%times
         END DO
         DEALLOCATE (Shared%elem)
         Shared%MaxConta = 2*Shared%MaxConta !!! 040717se atrancaba aqui. Ahora allocateo al doble. Antes era + 10000
         ALLOCATE (Shared%elem(1:Shared%MaxConta))
         DO n = 1, conta - 1
            Shared%elem(n)%Sharedmed = temp(n)%Sharedmed
            Shared%elem(n)%ProPmed = temp(n)%ProPmed
            Shared%elem(n)%field = temp(n)%field
            Shared%elem(n)%i = temp(n)%i
            Shared%elem(n)%j = temp(n)%j
            Shared%elem(n)%k = temp(n)%k
            Shared%elem(n)%times = temp(n)%times
         END DO
         DEALLOCATE (temp)
      END IF
      !
      IF (conta == 1) allocate (Shared%elem(1:Shared%MaxConta))
      Shared%elem(conta)%Sharedmed = Sharedmed
      Shared%elem(conta)%ProPmed = ProPmed
      Shared%elem(conta)%field = campo
      Shared%elem(conta)%i = i1
      Shared%elem(conta)%j = j1
      Shared%elem(conta)%k = k1
      Shared%elem(conta)%times = 2 !it appears two times at least !later updated in preporcess
      RETURN
      !
   END SUBROUTINE
   !
END MODULE
