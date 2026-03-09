
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Module :  CreateMatrices :  Fills in the media and observation matrices, creates
!                        the pml layers and the PEC borders.
!                        Also creates intermediate media for the boundaries
!                        between different media.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module CreateMatrices
   use Report
   use fdetypes
   !
   implicit none
   private
   !
   !
   type crosscheck_t
      integer(kind=4) :: actual, NewActual, NewActual2
      integer(kind=4), dimension(1:4) :: tent
   END type
   !matriz para controlar lo punietereos indices de cadacomponente
   integer(kind=4), dimension(6, 3, 2), parameter, public :: &
   & in = reshape ( (/ 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, &
   &                   0, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1 /), (/ 6, 3, 2 /))
   !
   integer(kind=4), dimension(6, 3, 2), parameter, public :: &
   &    on = reshape ( (/ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
   &                      0, 0,-1, 0, 0, 0,-1,-1, 0,-1, 0,-1, 0,-1, 0, 0, &
   -1,-1,-1, 0 /), (/ 6, 3, 2 /))
   !
   public CreatePMLmatrix, Readjust
   public CreateVolumeMM, CreateSurfaceMM, CreateLineMM
   public CreateSurfaceSlotMM,CreateMagneticSurface
   public CreateConformalPECVolume
   !
    contains
    
   subroutine SortInitEndWithIncreasingOrder(p)
        type(XYZlimit_t), intent(inout) :: p
        integer(kind=4) :: aux
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
   end subroutine
   
   subroutine CreateConformalPECVolume (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, BoundingBox,indicemedio)
      character(len=BUFSIZE) :: buff
      type(Shared_t) :: Eshared

      integer(kind=4) :: NumMedia
      type(MediaData_t), dimension(0:NumMedia) :: med
      !
      type(XYZlimit_t), intent(in) :: BoundingBox
      integer(kind=4), intent(in) :: indicemedio
      !
      integer(kind=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      integer(kind=IKINDMTAG) numertag
      integer(kind=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      !
      integer(kind=4) :: layoutnumber, i, j, k

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
               call fillPECFaceInsideVolume(i,j,k,FACE_X)
               call fillPECFaceInsideVolume(i,j,k,FACE_Y)
               call fillPECFaceInsideVolume(i,j,k,FACE_Z)
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
         integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: m1,m2,m3,m4
         hasCrossedPEC = (med(m1)%Is%ConformalPEC .or. med(m1)%Is%PEC) .and. &
                         (med(m2)%Is%ConformalPEC .or. med(m2)%Is%PEC) .and. &
                         (med(m3)%Is%ConformalPEC .or. med(m3)%Is%PEC) .and. &
                         (med(m4)%Is%ConformalPEC .or. med(m4)%Is%PEC)
      end function 
      logical function hasCrossedPECOrConformalPEC(m1,m2,m3,m4)
         integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: m1,m2,m3,m4
         hasCrossedPECOrConformalPEC = (med(m1)%Is%ConformalPEC .or. med(m1)%Is%PEC) .or. &
                                       (med(m2)%Is%ConformalPEC .or. med(m2)%Is%PEC) .or. &
                                       (med(m3)%Is%ConformalPEC .or. med(m3)%Is%PEC) .or. &
                                       (med(m4)%Is%ConformalPEC .or. med(m4)%Is%PEC)
      end function 

      subroutine fillPECFaceInsideVolume(i, j, k, face)
         integer(kind=4), intent(in) :: i, j, k, face
         integer(kind=4) :: m1,m2,m3,m4, m
         logical :: on_boundary
         select case (face)
         case(FACE_X)
            m1 = MMiEy (i, j, k)
            m2 = MMiEz (i, j, k)
            m3 = MMiEy (i, j, k + 1)
            m4 = MMiEz (i, j + 1, k)
            m = MMiHx(i,j,k)
         case (FACE_Y)
            m1 = MMiEx (i, j, k)
            m2 = MMiEz (i, j, k)
            m3 = MMiEx (i, j, k + 1)
            m4 = MMiEz (i + 1, j, k)
            m =  MMiHy(i,j,k)
         case (FACE_Z)
            m1 = MMiEy (i, j, k)
            m2 = MMiEx (i, j, k)
            m3 = MMiEy (i + 1, j, k)
            m4 = MMiEx (i, j + 1, k)
            m =  MMiHz(i,j,k)
         end select

         on_boundary = (med(m1)%Is%PEC .or. med(m1)%is%conformalPEC) .and. &
                       (med(m2)%Is%PEC .or. med(m2)%is%conformalPEC) .and. &
                       (med(m3)%Is%PEC .or. med(m3)%is%conformalPEC) .and. &
                       (med(m4)%Is%PEC .or. med(m4)%is%conformalPEC)

         if (on_boundary .and. .not. (med(m)%Is%PEC .or. med(m)%Is%ConformalPEC)) then 
            Mtag(i,j,k)=64*numertag 
            select case (face)
            case(FACE_X)
               MMiHx (i, j, k) = indicemedio
               tags%face%x(i,j,k) = 64*numertag
            case(FACE_Y)
               MMiHy (i, j, k) = indicemedio
               tags%face%y(i,j,k) = 64*numertag
            case(FACE_z)
               MMiHz (i, j, k) = indicemedio
               tags%face%z(i,j,k) = 64*numertag
            end select
         else if (on_boundary .and. (med(m)%Is%PEC .or. med(m)%Is%ConformalPEC)) then 
            Mtag(i,j,k)=64*numertag
            select case (face)
            case(FACE_X)
               tags%face%x(i,j,k) = 64*numertag
            case(FACE_Y)
               tags%face%y(i,j,k) = 64*numertag
            case(FACE_z)
               tags%face%z(i,j,k) = 64*numertag
            end select
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
   subroutine CreateVolumeMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz, Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, BoundingBox, point, indicemedio)
      character(len=BUFSIZE) :: buff
      type(Shared_t) :: Eshared
      !
      integer(kind=4) :: NumMedia
      type(MediaData_t), dimension(0:NumMedia) :: med
      integer(kind=4) :: medio
      !
      type(XYZlimit_t) :: punto, puntoPlus1
      type(XYZlimit_t), intent(inout) :: point
      type(XYZlimit_t), intent(in) :: BoundingBox
      !
      integer(kind=4) :: indicemedio
      !
      integer(kind=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      integer(kind=IKINDMTAG) numertag
      integer(kind=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      !
      integer(kind=4) :: layoutnumber, i, j, k
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
      !!      do k = punto%ZI, punto%ZE
      !!        do j = punto%YI, punto%YE
      !!          do i = punto%XI, punto%XE
      !!            medio = MMcen (i, j, k)
      !!            if (med(indicemedio)%Priority >= med(medio)%Priority) then
      !!              MMcen (i, j, k) = indicemedio
      !!            end if
      !!          end do
      !!        end do
      !!      end do
      !only take care of the boundaries for interfacing
      do k = punto%ZI, puntoPlus1%ZE
         do j = punto%YI, puntoPlus1%YE
            do i = punto%XI, punto%XE
               medio = MMiEx (i, j, k)
!               if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0); 
                     !ojo no es sumar porque no debe desbordarse. solo hay que poner el bit
                     !solo se pone un tag si ya se habia puesto o si estaba si inicializar (cero) y se shifte 6 bits numerados del 0 (iex) al 5 (ihz) empezando por la derecha (lsb)
                     !lo ponto siempre a .true. y quien llege se lo lleva machacando al que habia. los tags no solucionan el problema de determinar univocamente el medio de una celda
                     !a menos que se definan 6 matrices de tags. esto es un niapa que solo sirve para filtrar celdas incluyendo fallos en celdas compartidas entre medios pero no es fiable para determinar medios en posiciones 161020
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     !no lo detectare en volumenes porque podria llevar tiempos elevados en el preproceso
                     !cuando se actualiza el numero de shared (sept'11)
                     !        OnSurface = (k == punto%ZI).or.(k == puntoPlus1%ZE).or.(j == punto%YI).or.(j == puntoPlus1%YE)
                     !        if (OnSurface) call AddToShared(iEx,i,j,k,indicemedio,medio,Eshared)
                  end if
!               end if
            end do
         end do
      end do
      !
      do k = punto%ZI, puntoPlus1%ZE
         do j = punto%YI, punto%YE
            do i = punto%XI, puntoPlus1%XE
               medio = MMiEy (i, j, k)
!               if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     !no lo detectare en volumenes porque podria llevar tiempos elevados en el preproceso
                     !cuando se actualiza el numero de shared (sept'11)
                     !        OnSurface = (k == punto%ZI).or.(k == puntoPlus1%ZE).or.(i == punto%XI).or.(i == puntoPlus1%XE)
                     !        if (OnSurface) call AddToShared(iEy,i,j,k,indicemedio,medio,Eshared)
                  end if
                  
 !              end if
            end do
         end do
      end do
      !
      do k = punto%ZI, punto%ZE
         do j = punto%YI, puntoPlus1%YE
            do i = punto%XI, puntoPlus1%XE
               medio = MMiEz (i, j, k)
!               if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%edge%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     !no lo detectare en volumenes porque podria llevar tiempos elevados en el preproceso
                     !cuando se actualiza el numero de shared (sept'11)
                     !        OnSurface = (i == punto%XI).or.(i == puntoPlus1%XE).or.(j == punto%YI).or.(j == puntoPlus1%YE)
                     !        if (OnSurface) call AddToShared(iEz,i,j,k,indicemedio,medio,Eshared)
                  end if
!               end if
            end do
         end do
      end do
      !
      do k = punto%ZI, punto%ZE
         do j = punto%YI, punto%YE
            do i = punto%XI, puntoPlus1%XE
               medio = MMiHx (i, j, k)
!               if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  end if
!               end if
            end do
         end do
      end do
      !
      do k = punto%ZI, punto%ZE
         do j = punto%YI, puntoPlus1%YE
            do i = punto%XI, punto%XE
               medio = MMiHy (i, j, k)
!               if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  end if
 !              end if
            end do
         end do
      end do
      !
      do k = punto%ZI, puntoPlus1%ZE
         do j = punto%YI, punto%YE
            do i = punto%XI, punto%XE
               medio = MMiHz (i, j, k)
!               if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag 
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  end if
!               end if
            end do
         end do
      end do
      !
      return
   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateSurfaceMM :  Sets every field component of the lower/back/left surface of a voxel to the index of the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   !          orientacion       : Plane of the surface affected by this medium (iEx,iEy,iEz)
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each voxel centered at i,j,k
   !                                        (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine CreateSurfaceMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  &
   & Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, &
   & Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, &
   & Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, &
   & Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, &
   & Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, &
   & med, NumMedia, Eshared, BoundingBox, point, orientacion, indicemedio)
      character(len=BUFSIZE) :: buff
      integer(kind=4) :: NumMedia
      type(Shared_t) :: Eshared
      type(MediaData_t), dimension(0:NumMedia) :: med
      !
      type(XYZlimit_t) :: punto, puntoPlus1,puntoBboxplus1
      type(XYZlimit_t), intent(inout) :: point
      type(XYZlimit_t), intent(in) :: BoundingBox
      !
      integer(kind=4) :: indicemedio, orientacion
      integer(kind=4) :: layoutnumber, i, j, k
      integer(kind=4) :: medio
      !
      integer(kind=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI , Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      integer(kind=IKINDMTAG) numertag
      integer(kind=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
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
         do i = punto%XI, puntoBboxplus1%XE
            do j = punto%YI, punto%YE
               do k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEy (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEy (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     tags%edge%y(i,j,k) = 64*numertag
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     call AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                  end if
               end do
            end do
            do j = punto%YI, puntoPlus1%YE
               do k = punto%ZI, punto%ZE
                  medio = MMiEz (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEz (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     tags%edge%z(i,j,k) = 64*numertag
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     call AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                  end if
               end do
            end do  
            do j = punto%YI, punto%YE
               do k = punto%ZI, punto%ZE
                  medio = MMiHx (i, j, k)
!                  if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                         MMiHx (i, j, k) = indicemedio; 
                         Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                         tags%face%x(i,j,k) = 64*numertag
                     endif
!                  end if
               end do
            end do
         end do
         !    endif
       CASE (iEy)
         !    j=punto%YI
         !    if ((j <= max(BoundingBox%YI,BoundingBox%YE)).and.(j >= min(BoundingBox%YI,BoundingBox%YE))) then
         do j = punto%YI, puntoBboxplus1%YE
            do i = punto%XI, puntoPlus1%XE
               do k = punto%ZI, punto%ZE
                  medio = MMiEz (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEz (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     tags%edge%z(i,j,k) = 64*numertag
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     call AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                  end if
               end do
            end do
            do i = punto%XI, punto%XE
               do k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEx (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEx (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     tags%edge%x(i,j,k) = 64*numertag
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     call AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                  end if
               end do
            end do
            do i = punto%XI, punto%XE
               do k = punto%ZI, punto%ZE
                  medio = MMiHy (i, j, k)
!                  if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                         MMiHy (i, j, k) = indicemedio; 
                         Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);;
                         tags%face%y(i,j,k) = 64*numertag
                     endif
!                  end if
               end do
            end do
         end do
         !    endif
       CASE (iEz)
         !    k=punto%ZI
         !    if ((k <= max(BoundingBox%ZI,BoundingBox%ZE)).and.(k >= min(BoundingBox%ZI,BoundingBox%ZE))) then
         do k = punto%ZI, puntoBboxplus1%ZE
            do i = punto%XI, punto%XE
               do j = punto%YI, puntoPlus1%YE
                  medio = MMiEx (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEx (i, j, k) = indicemedio; 
                     Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     tags%edge%x(i,j,k) = 64*numertag
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     call AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                  end if
               end do
            end do
            do i = punto%XI, puntoPlus1%XE
               do j = punto%YI, punto%YE
                  medio = MMiEy (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEy (i, j, k) = indicemedio; Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     tags%edge%y(i,j,k) = 64*numertag
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     call AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                  end if
               end do
            end do
            do i = punto%XI, punto%XE
               do j = punto%YI, punto%YE
                  medio = MMiHz (i, j, k)
!                  if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                         MMiHz (i, j, k) = indicemedio; 
                         Mtag(i,j,k)=64*numertag ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                         tags%face%z(i,j,k) = 64*numertag

                     endif
!                  end if
               end do
            end do
         end do
         !    endif
      end select
      !
      return
   end subroutine
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Routine :  CreateLineMM :  Sets every field component of the inner Y/Y/Z axis of a voxel to the index of the medium
   ! Inputs :   M(field)%Mediamatrix(i,j,k)  : type of medium at each i,j,k, for each field
   !          punto%XI,punto%XE,punto%YI,punto%YE,punto%ZI,punto%ZE : initial and end coordinates of the voxel
   !          indicemedio       : index of the voxel medium
   !          orientacion       : Axis of the voxel affected by this medium (iEx,iEy,iEz)
   ! Outputs :  M(field)%Mediamatrix(i,j,k) = type of medium indicemedio set for all the fields at each
   !                                        voxel centered at i,j,k (usual convention)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine CreateLineMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, BoundingBox, point, orientacion, indicemedio, isathinwire, verbose,numeroasignaciones)
      
      type(Shared_t) :: Eshared
      integer(kind=4) :: NumMedia
      type(MediaData_t), dimension(0:NumMedia) :: med
      !
      type(XYZlimit_t) :: punto
      type(XYZlimit_t), intent(inout) :: point
      type(XYZlimit_t), intent(in) :: BoundingBox
      
      integer(kind=4) :: indicemedio, orientacion,numeroasignaciones
      LOGICAL, intent(in) :: isathinwire, verbose
      integer(kind=4) :: i, j, k, layoutnumber
      integer(kind=4) :: medio
      !
      integer(kind=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      integer(kind=IKINDMTAG) numertag
      integer(kind=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      !
      character(len=BUFSIZE) :: buff
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
         do k = punto%ZI, punto%ZE
            do j = punto%YI, punto%YE
               do i = punto%XI, punto%XE
                  medio = MMiEx (i, j, k)
!                  if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
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
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        call AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                     end if
               end do
            end do
         end do
         !    endif
       CASE (iEy)
         !    i=punto%XI
         !    k=punto%ZI
         !    if ((i <= max(BoundingBox%XI,BoundingBox%XE)).and.(i >= min(BoundingBox%XI,BoundingBox%XE)).and. &
         !        (k <= max(BoundingBox%ZI,BoundingBox%ZE)).and.(k >= min(BoundingBox%ZI,BoundingBox%ZE))) then
         do k = punto%ZI, punto%ZE
            do j = punto%YI, punto%YE
               do i = punto%XI, punto%XE
                  medio = MMiEy (i, j, k)
!                  if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
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
                        
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        call AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                     end if
               end do
            end do
         end do
         !    endif
       CASE (iEz)
         !    i=punto%XI
         !    j=punto%YI
         !    if ((i <= max(BoundingBox%XI,BoundingBox%XE)).and.(i >= min(BoundingBox%XI,BoundingBox%XE)).and. &
         !        (j <= max(BoundingBox%YI,BoundingBox%YE)).and.(j >= min(BoundingBox%YI,BoundingBox%YE))) then
         do k = punto%ZI, punto%ZE
            do j = punto%YI, punto%YE
               do i = punto%XI, punto%XE
                  medio = MMiEz (i, j, k)
!                  if (medio /= 0) then   !ojo esto estaba antes de 031016 y daba maxima prioridad al medio 0 PEC. Ahora puedo tener medios con mas prioridad!!! !?!? cambio agresivo 031016!!!
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
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
                        
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        call AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                     end if

               end do
            end do
         end do
         !    endif
      end select
      !
      return
   end subroutine
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
   subroutine CreateSurfaceSlotMM (layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, Hshared, BoundingBox, point, orientacion, direccion, indicemedio)
      character(len=BUFSIZE) :: buff
      type(Shared_t) :: Eshared, Hshared
      integer(kind=4) :: NumMedia
      type(MediaData_t), dimension(0:NumMedia) :: med
      !
      type(XYZlimit_t) :: punto, puntoPlus1,puntoBboxplus1
      type(XYZlimit_t), intent(inout) :: point
      type(XYZlimit_t), intent(in) :: BoundingBox
      !
      integer(kind=4) :: indicemedio, orientacion, direccion
      !
      integer(kind=4) :: layoutnumber, i, j, k, offx, offy, offz
      integer(kind=4) :: medio
      !
      integer(kind=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      integer(kind=IKINDMTAG) numertag
      integer(kind=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
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
         do i = punto%XI, puntoBboxplus1%XE
            SELECT CASE (direccion)
             CASE (iEz)
               offx = 0
               offy = 0
               offz = 1
               do j = punto%YI, punto%YE
                  do k = punto%ZI, puntoPlus1%ZE
                     medio = MMiEy (i, j, k)
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                        MMiEy (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%y(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        !call AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                     end if
                  end do
               end do
             CASE (iEy)
               offx = 0
               offy = 1
               offz = 0
               do j = punto%YI, puntoPlus1%YE
                  do k = punto%ZI, punto%ZE
                     medio = MMiEz (i, j, k)
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                        MMiEz (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%z(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        !call AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                     end if
                  end do
               end do
            end select
            do j = Max (punto%YI - offy, Min(BoundingBox%YI, BoundingBox%YE)), &
            &       Min (punto%YE + offy, Max(BoundingBox%YI, BoundingBox%YE)-1)
               do k = Max (punto%ZI - offz, Min(BoundingBox%ZI, BoundingBox%ZE)),  &
               &       Min (punto%ZE + offz, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
                  medio = MMiHx (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     !call AddToShared (iHx, i, j, k, indicemedio, medio, Hshared)

                  end if
               end do
            end do
         end do
       CASE (iEy)
         do j = punto%YI, puntoBboxplus1%YE
            SELECT CASE (direccion)
             CASE (iEx)
               offx = 1
               offy = 0
               offz = 0
               do i = punto%XI, puntoPlus1%XE
                  do k = punto%ZI, punto%ZE
                     medio = MMiEz (i, j, k)
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                        MMiEz (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%z(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        !call AddToShared (iEz, i, j, k, indicemedio, medio, Eshared)
                     end if
                  end do
               end do
             CASE (iEz)
               offx = 0
               offy = 0
               offz = 1
               do i = punto%XI, punto%XE
                  do k = punto%ZI, puntoPlus1%ZE
                     medio = MMiEx (i, j, k)
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                        MMiEx (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%x(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        !call AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                     end if
                  end do
               end do
            end select
            do i = Max (punto%XI - offx, Min(BoundingBox%XI, BoundingBox%XE)),  &
            &       Min (punto%XE + offx, Max(BoundingBox%XI, BoundingBox%XE)-1)
               do k = Max (punto%ZI - offz, Min(BoundingBox%ZI, BoundingBox%ZE)),  &
               &       Min (punto%ZE + offz, Max(BoundingBox%ZI, BoundingBox%ZE)-1)
                  medio = MMiHy (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     !call AddToShared (iHy, i, j, k, indicemedio, medio, Hshared)
                  end if
               end do
            end do
         end do
       CASE (iEz)
         do k = punto%ZI, puntoBboxplus1%ZE
            SELECT CASE (direccion)
             CASE (iEy)
               offx = 0
               offy = 1
               offz = 0
               do i = punto%XI, punto%XE
                  do j = punto%YI, puntoPlus1%YE
                     medio = MMiEx (i, j, k)
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                        MMiEx (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%x(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        !call AddToShared (iEx, i, j, k, indicemedio, medio, Eshared)
                     end if
                  end do
               end do
             CASE (iEx)
               offx = 1
               offy = 0
               offz = 0
               do i = punto%XI, puntoPlus1%XE
                  do j = punto%YI, punto%YE
                     medio = MMiEy (i, j, k)
                     if (med(indicemedio)%Priority > med(medio)%Priority) then
                        MMiEy (i, j, k) = indicemedio
                        Mtag(i,j,k)=64*numertag
                        tags%edge%y(i,j,k) = 64*numertag
                        ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                     ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                        !call AddToShared (iEy, i, j, k, indicemedio, medio, Eshared)
                     end if
                  end do
               end do
            end select
            do i = Max (punto%XI - offx, Min(BoundingBox%XI, BoundingBox%XE)),  &
            &       Min (punto%XE + offx, Max(BoundingBox%XI, BoundingBox%XE)-1)
               do j = Max (punto%YI - offy, Min(BoundingBox%YI, BoundingBox%YE)),  &
               &       Min (punto%YE + offy, Max(BoundingBox%YI, BoundingBox%YE)-1)
                  medio = MMiHz (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  ELSE if ((med(indicemedio)%Priority == med(medio)%Priority) .AND. (medio /= indicemedio)) then
                     !call AddToShared (iHz, i, j, k, indicemedio, medio, Hshared)
                  end if
               end do
            end do
         end do
      end select
      !
      return
   end subroutine
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
   subroutine CreateMagneticSurface(layoutnumber, Mtag, tags, numertag, MMiEx, MMiEy, MMiEz, MMiHx, &
   & MMiHy, MMiHz,  Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, &
   & Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, &
   & Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, &
   & Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, &
   & Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, med, &
   & NumMedia, Eshared, BoundingBox, point, orientacion, indicemedio)
      character(len=BUFSIZE) :: buff
      integer(kind=4) :: NumMedia
      type(Shared_t) :: Eshared
      type(MediaData_t), dimension(0:NumMedia) :: med
      !
      type(XYZlimit_t) :: punto, puntoPlus1   !,puntoBboxplus1
      type(XYZlimit_t), intent(inout) :: point
      type(XYZlimit_t), intent(in) :: BoundingBox
      !
      integer(kind=4) :: indicemedio, orientacion
      integer(kind=4) :: layoutnumber, i, j, k
      integer(kind=4) :: medio
      !
      integer(kind=4) :: Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, &
      & Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, &
      & Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, &
      & Alloc_iHx_ZE, Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, &
      & Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE
      !
      type(taglist_t) :: tags
      integer(kind=IKINDMTAG) numertag
      integer(kind=IKINDMTAG ) :: Mtag  (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEx (Alloc_iEx_XI:Alloc_iEx_XE, Alloc_iEx_YI:Alloc_iEx_YE, Alloc_iEx_ZI:Alloc_iEx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEy (Alloc_iEy_XI:Alloc_iEy_XE, Alloc_iEy_YI:Alloc_iEy_YE, Alloc_iEy_ZI:Alloc_iEy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiEz (Alloc_iEz_XI:Alloc_iEz_XE, Alloc_iEz_YI:Alloc_iEz_YE, Alloc_iEz_ZI:Alloc_iEz_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHx (Alloc_iHx_XI:Alloc_iHx_XE, Alloc_iHx_YI:Alloc_iHx_YE, Alloc_iHx_ZI:Alloc_iHx_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHy (Alloc_iHy_XI:Alloc_iHy_XE, Alloc_iHy_YI:Alloc_iHy_YE, Alloc_iHy_ZI:Alloc_iHy_ZE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: MMiHz (Alloc_iHz_XI:Alloc_iHz_XE, Alloc_iHz_YI:Alloc_iHz_YE, Alloc_iHz_ZI:Alloc_iHz_ZE)
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
         do i = punto%XI, punto%XE
            do j = punto%YI, puntoPlus1%YE
               do k = punto%ZI, punto%ZE
                  medio = MMiHy (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  endif
               end do
            end do
            do j = punto%YI, punto%YE
               do k = punto%ZI, puntoPlus1%ZE
                  medio = MMiHz (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  endif
               end do
            end do
            do j = punto%YI, puntoPlus1%YE
               do k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEx (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%edge%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,0);
                  endif
               end do
            end do
         end do
       CASE (iEy)
         do j = punto%YI, punto%YE
            do i = punto%XI, punto%XE
               do k = punto%ZI, puntoPlus1%ZE
                  medio = MMiHz (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,5);
                  endif
               end do
            end do
            do i = punto%XI, puntoPlus1%XE
               do k = punto%ZI, punto%ZE
                  medio = MMiHx (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  end if
               end do
            end do
            do i = punto%XI, puntoPlus1%XE
               do k = punto%ZI, puntoPlus1%ZE
                  medio = MMiEy (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%edge%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,1);
                  endif
               end do
            end do
         end do
         !
       CASE (iEz)
         do k = punto%ZI, punto%ZE
            do i = punto%XI, puntoPlus1%XE
               do j = punto%YI, punto%YE
                  medio = MMiHx (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHx (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%x(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,3);
                  end if
               end do
            end do
            do i = punto%XI, punto%XE
               do j = punto%YI, puntoPlus1%YE
                  medio = MMiHy (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiHy (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%face%y(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,4);
                  end if
               end do
            end do
            do i = punto%XI, puntoPlus1%XE
               do j = punto%YI, puntoPlus1%YE
                  medio = MMiEz (i, j, k)
                  if (med(indicemedio)%Priority > med(medio)%Priority) then
                     MMiEz (i, j, k) = indicemedio
                     Mtag(i,j,k)=64*numertag
                     tags%edge%z(i,j,k) = 64*numertag
                     ! if (.true..or.(Mtag(i,j,k)==0).or.(int(Mtag(i,j,k)/64) == numertag)) Mtag(i,j,k) = IBSET(64*numertag,2);
                  endif
               end do
            end do
         end do

      end select
      !
      return
   end subroutine


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
   subroutine CreatePMLmatrix (layoutnumber, SIZE, sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, SINPML_fullsize, fullsize, BBox, med, NumMedia, Border, MEDIOEXTRA)
      !
      type(limit_t), dimension(1:6) :: SINPML_fullsize, fullsize
      !Inputs and Outputs
      type(SGGFDTDINFO_t), intent(INOUT) :: sgg
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: &
      sggMiEx(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE), &
      sggMiEy(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE), &
      sggMiEz(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE), &
      sggMiHx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE), &
      sggMiHy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE), &
      sggMiHz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)
      type(XYZlimit_t) :: BoundingBox, BBox
      type(MediaData_t), pointer, dimension(:) :: med
      integer(kind=4), intent(inout) :: NumMedia
      type(Border_t), intent(in) :: Border
      type(MedioExtra_t), intent(in) :: MEDIOEXTRA
      ! Local stuff
      integer(kind=4), pointer, dimension(:) :: tempo
      type(MediaData_t), pointer, dimension(:) :: NewMed
      integer(kind=4) :: layoutnumber, SIZE, field, medium, i, j, k, NuevoNumeroMediosConPML
      integer(kind=4) :: oldNumMedia,oldmed
      integer(kind=4), dimension(1:6) :: XIPML, XEPML, YIPML, YEPML, ZIPML, ZEPML
      real(kind=RKIND) :: oldepr,oldmur,oldsigma,oldsigmam,newepr,newmur,newsigma,newsigmam
      LOGICAL :: yapuesto

      !Increase mediamatrix with PML regions, and remove the minus sign from the mm matrix (this info is no longer needed)
      !FIRST CLIP THE MATRIX
      !readjust boundingbox for PML correct calculation
      !
      do field = iEx, iHz
         XIPML (field) = Max (BBox%XI, fullsize(iHx)%XI)
         XEPML (field) = Min (BBox%XE+on(field, icoord, fine), fullsize(iHx)%XE)
         YIPML (field) = Max (BBox%YI, fullsize(iHy)%YI)
         YEPML (field) = Min (BBox%YE+on(field, jcoord, fine), fullsize(iHy)%YE)
         ZIPML (field) = Max (BBox%ZI, fullsize(iHz)%ZI)
         ZEPML (field) = Min (BBox%ZE+on(field, kcoord, fine), fullsize(iHz)%ZE)
      end do
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
      do j = YIPML (field), YEPML (field)
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            do i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiEx (i, j, k) = - Abs (sggmiEx(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Front
            do i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiEx (i, j, k) = - Abs (sggmiEx(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            do j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Right
            do j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            do k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Up
            do k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiEx (i, j, k) = - Abs (sggmiEx(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      field = iEy
      do j = YIPML (field), YEPML (field)
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            do i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiEy (i, j, k) = - Abs (sggmiEy(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Front
            do i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiEy (i, j, k) = - Abs (sggmiEy(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            do j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Right
            do j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            do k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Up
            do k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiEy (i, j, k) = - Abs (sggmiEy(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      field = iEz
      do j = YIPML (field), YEPML (field)
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            do i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiEz (i, j, k) = - Abs (sggmiEz(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Front
            do i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiEz (i, j, k) = - Abs (sggmiEz(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            do j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Right
            do j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            do k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Up
            do k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiEz (i, j, k) = - Abs (sggmiEz(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      field = iHx
      do j = YIPML (field), YEPML (field)
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            do i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiHx (i, j, k) = - Abs (sggmiHx(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Front
            do i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiHx (i, j, k) = - Abs (sggmiHx(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            do j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Right
            do j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            do k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Up
            do k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiHx (i, j, k) = - Abs (sggmiHx(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      field = iHy
      do j = YIPML (field), YEPML (field)
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            do i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiHy (i, j, k) = - Abs (sggmiHy(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Front
            do i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiHy (i, j, k) = - Abs (sggmiHy(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            do j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Right
            do j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            do k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Up
            do k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiHy (i, j, k) = - Abs (sggmiHy(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      field = iHz
      do j = YIPML (field), YEPML (field)
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Back
            do i = XIPML (field), BoundingBox%XI + in (field, icoord, fine)
               sggmiHz (i, j, k) = - Abs (sggmiHz(BoundingBox%XI+in(field, icoord, fine)+1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Front
            do i = BoundingBox%XE + in (field, icoord, comi), XEPML (field)
               sggmiHz (i, j, k) = - Abs (sggmiHz(BoundingBox%XE+in(field, icoord, comi)-1, j, k))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do k = ZIPML (field), ZEPML (field)
            !!!!!!**Left
            do j = YIPML (field), BoundingBox%YI + in (field, jcoord, fine)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, BoundingBox%YI+in(field, jcoord, fine)+1, k))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Right
            do j = BoundingBox%YE + in (field, jcoord, comi), YEPML (field)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, BoundingBox%YE+in(field, jcoord, comi)-1, k))
            end do
         end do
      end do
      !
      do i = XIPML (field), XEPML (field)!barre ahora el total, para incluir aristas y corners
         do j = YIPML (field), YEPML (field)!barre ahora el total, para incluir aristas y corners
            !!!!!!**Down
            do k = ZIPML (field), BoundingBox%ZI + in (field, kcoord, fine)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, j, BoundingBox%ZI+in(field, kcoord, fine)+1))
               !para notar medio PML se le cambia el signo al medio
            end do
            !!!!!!**Up
            do k = BoundingBox%ZE + in (field, kcoord, comi), ZEPML (field)
               sggmiHz (i, j, k) = - Abs (sggmiHz(i, j, BoundingBox%ZE+in(field, kcoord, comi)-1))
               !para notar medio PML se le cambia el signo al medio
            end do
         end do
      end do
      !compact the info of PML media
      NuevoNumeroMediosConPML = NumMedia
     allocate(tempo(0:NumMedia))
      tempo = 0 !temporarily stores the index of the PML medium matching each original media
      field = iEx
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiEx (i, j, k)
               if (medium < 0) then
                  if (tempo(Abs(medium)) == 0) then
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  end if
               end if
            end do
         end do
      end do
      field = iEy
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiEy (i, j, k)
               if (medium < 0) then
                  if (tempo(Abs(medium)) == 0) then
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  end if
               end if
            end do
         end do
      end do
      field = iEz
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiEz (i, j, k)
               if (medium < 0) then
                  if (tempo(Abs(medium)) == 0) then
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  end if
               end if
            end do
         end do
      end do
      field = iHx
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiHx (i, j, k)
               if (medium < 0) then
                  if (tempo(Abs(medium)) == 0) then
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  end if
               end if
            end do
         end do
      end do
      field = iHy
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiHy (i, j, k)
               if (medium < 0) then
                  if (tempo(Abs(medium)) == 0) then
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  end if
               end if
            end do
         end do
      end do
      field = iHz
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiHz (i, j, k)
               if (medium < 0) then
                  if (tempo(Abs(medium)) == 0) then
                     NuevoNumeroMediosConPML = NuevoNumeroMediosConPML + 1
                     tempo (Abs(medium)) = NuevoNumeroMediosConPML
                  end if
               end if
            end do
         end do
      end do
      !
     allocate(NewMed(NumMedia+1:NuevoNumeroMediosConPML))
      !Reassing the PML media info with the compact indexes
      field = iEx
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiEx (i, j, k)
               if (medium < 0) then
                  sggmiEx (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               end if
            end do
         end do
      end do
      field = iEy
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiEy (i, j, k)
               if (medium < 0) then
                  sggmiEy (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               end if
            end do
         end do
      end do
      field = iEz
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiEz (i, j, k)
               if (medium < 0) then
                  sggmiEz (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               end if
            end do
         end do
      end do
      field = iHx
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiHx (i, j, k)
               if (medium < 0) then
                  sggmiHx (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               end if
            end do
         end do
      end do
      field = iHy
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiHy (i, j, k)
               if (medium < 0) then
                  sggmiHy (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               end if
            end do
         end do
      end do
      field = iHz
      do k = ZIPML (field), ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               medium = sggmiHz (i, j, k)
               if (medium < 0) then
                  sggmiHz (i, j, k) = tempo (Abs(medium))
                  NewMed (tempo(Abs(medium))) = med (Abs(medium))
               end if
            end do
         end do
      end do

      !Put PEC and the end if there exists PEC borders in the original problem
      !(PMC are handled with the image technique in the algorithm, no special index is used for PMC)
      !DETRAS Y DELANTE
      field = iEx !!!!!PEC solo en fields donde acabe la red
      !izda y dcha
      if ((Border%IsLeftPEC)) then
         j = YIPML (field)
         do i = XIPML (field), XEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEx (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsRightPEC)) then
         j = YEPML (field)
         do i = XIPML (field), XEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEx (i, j, k) = 0
            end do
         end do
      end if
      !  !Up y Down
      if ((Border%IsDownPEC)) then
         k = ZIPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               sggmiEx (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsUpPEC)) then
         k = ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               sggmiEx (i, j, k) = 0
            end do
         end do
      end if
      !
      field = iEy !!!!!PEC solo en fields donde acabe la red
      !front y back
      if ((Border%IsBackPEC)) then
         i = XIPML (field)
         do j = YIPML (field), YEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEy (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsFrontPEC)) then
         i = XEPML (field)
         do j = YIPML (field), YEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEy (i, j, k) = 0
            end do
         end do
      end if
      !  !Up y Down
      if ((Border%IsDownPEC)) then
         k = ZIPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               sggmiEy (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsUpPEC)) then
         k = ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               sggmiEy (i, j, k) = 0
            end do
         end do
      end if
      !
      field = iEz !!!!!PEC solo en fields donde acabe la red
      !front y back
      if ((Border%IsBackPEC)) then
         i = XIPML (field)
         do j = YIPML (field), YEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsFrontPEC)) then
         i = XEPML (field)
         do j = YIPML (field), YEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            end do
         end do
      end if
      !izda y dcha
      if ((Border%IsLeftPEC)) then
         j = YIPML (field)
         do i = XIPML (field), XEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsRightPEC)) then
         j = YEPML (field)
         do i = XIPML (field), XEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiEz (i, j, k) = 0
            end do
         end do
      end if
      !
      field = iHx !!!!!PEC solo en fields donde acabe la red
      !front y back
      if ((Border%IsBackPEC)) then
         i = XIPML (field)
         do j = YIPML (field), YEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiHx (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsFrontPEC)) then
         i = XEPML (field)
         do j = YIPML (field), YEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiHx (i, j, k) = 0
            end do
         end do
      end if
      !
      field = iHy !!!!!PEC solo en fields donde acabe la red
      !izda y dcha
      if ((Border%IsLeftPEC)) then
         j = YIPML (field)
         do i = XIPML (field), XEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiHy (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsRightPEC)) then
         j = YEPML (field)
         do i = XIPML (field), XEPML (field)
            do k = ZIPML (field), ZEPML (field)
               sggmiHy (i, j, k) = 0
            end do
         end do
      end if
      !
      field = iHz !!!!!PEC solo en fields donde acabe la red
      !  !Up y Down
      if ((Border%IsDownPEC)) then
         k = ZIPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               sggmiHz (i, j, k) = 0
            end do
         end do
      end if
      !
      if ((Border%IsUpPEC)) then
         k = ZEPML (field)
         do j = YIPML (field), YEPML (field)
            do i = XIPML (field), XEPML (field)
               sggmiHz (i, j, k) = 0
            end do
         end do
      end if

      !!?!?!?MEDIOEXTRA

      !
      !adjust constitutive parameters, matrices
      oldNumMedia = NumMedia !save it since the next subroutine overwrites it
      call Readjust (NumMedia, med, NuevoNumeroMediosConPML)
      sgg%AllocMed=NumMedia
      !copy the new media
      med (1+oldNumMedia:NuevoNumeroMediosConPML) = NewMed (1+oldNumMedia:NuevoNumeroMediosConPML)
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%PML = .TRUE. !all these are PML
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%ThinWire = .FALSE. !put any wire touching the PML to non-wire though treat it with mur
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%SlantedWire = .FALSE. !put any wire touching the PML to non-wire though treat it with mur
      med(1+oldNumMedia:NuevoNumeroMediosConPML)%Is%Needed=.true. !sgg 220817 por defecto lo he puesto en readjust a false
      !
      deallocate(NewMed, tempo)

      !solo lo creo para las tangenciales electricas
      if (MEDIOEXTRA%exists) then
         !Put MEDIO and the end if there exists PML borders in the original problem
         yapuesto=.false.
         !
         field = iEx
         !izda y dcha
         if ((Border%IsLeftPML)) then
            do j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsRightPML)) then
            do j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !  !Up y Down
         if ((Border%IsDownPML)) then
            do k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsUpPML)) then
            do k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         field = iEy
         !front y back
         if ((Border%IsBackPML)) then
            do i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsFrontPML)) then
            do i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !  !Up y Down
         if ((Border%IsDownPML)) then
            do k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsUpPML)) then
            do k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         field = iEz
         !front y back
         if ((Border%IsBackPML)) then
            do i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsFrontPML)) then
            do i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !izda y dcha
         if ((Border%IsLeftPML)) then
            do j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsRightPML)) then
            do j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         !magneticos

         !
         !
         field = iHx
         !izda y dcha
         if ((Border%IsLeftPML)) then
            do j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsRightPML)) then
            do j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !  !Up y Down
         if ((Border%IsDownPML)) then
            do k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsUpPML)) then
            do k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         field = iHy
         !front y back
         if ((Border%IsBackPML)) then
            do i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsFrontPML)) then
            do i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !  !Up y Down
         if ((Border%IsDownPML)) then
            do k = ZIPML (field),ZIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsUpPML)) then
            do k = ZEPML (field)- MEDIOEXTRA%size,ZEPML (field)
               do j = YIPML (field), YEPML (field)
                  do i = XIPML (field), XEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         field = iHz
         !front y back
         if ((Border%IsBackPML)) then
            do i = XIPML (field),XIPML (field)+ MEDIOEXTRA%size
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsFrontPML)) then
            do i = XEPML (field)- MEDIOEXTRA%size,XEPML (field)
               do j = YIPML (field), YEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !izda y dcha
         if ((Border%IsLeftPML)) then
            do j = YIPML (field),YIPML (field)+ MEDIOEXTRA%size
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
         if ((Border%IsRightPML)) then
            do j = YEPML (field)- MEDIOEXTRA%size,YEPML (field)
               do i = XIPML (field), XEPML (field)
                  do k = ZIPML (field), ZEPML (field)
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
                              call STOPONERROR (layoutnumber,size,'Multilayer corrected PML unsupported. Relaunch without -pmlcorr')
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
                  end do
               end do
            end do
         end if
         !
      endif !del medioextra%exists
      return
   end subroutine
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
   subroutine Readjust (NumMedia, med, NewNumMedia)
      !Inputs and Outputs
      type(MediaData_t), pointer, dimension(:) :: med
      type(MediaData_t), pointer, dimension(:) :: dummyMed
      integer(kind=4), intent(inout) :: NumMedia
      integer(kind=4), intent(in) :: NewNumMedia
      integer(kind=4) :: i
      !
     allocate(dummyMed(0:NewNumMedia))
      do i = 0, Min (NumMedia, NewNumMedia)
         dummyMed (i) = med (i)
      end do
      !
      deallocate(med)
     allocate(med(0:NewNumMedia))
      do i = 0, NewNumMedia
         med (i) = dummyMed (i)
      end do
      !
      deallocate(dummyMed)
      do i = 1 + NumMedia, NewNumMedia
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
         med(i)%Is%Multiwire = .FALSE.
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
      end do
      !Update the final number of media
      NumMedia = NewNumMedia
      !
   end subroutine

   subroutine AddToShared (campo, i1, j1, k1, Sharedmed, ProPmed,  Shared)
      type(SharedElement_t), pointer, dimension(:) :: temp
      type(Shared_t), intent(inout) :: Shared
      integer(kind=4), intent(in) :: campo, i1, j1, k1, Sharedmed, ProPmed
      integer(kind=4) :: conta, n
      !
      Shared%conta = Shared%conta + 1
      conta = Shared%conta
      !
      if (conta > Shared%MaxConta) then
         !create space on the fly
        allocate(temp(1:conta-1))
         do n = 1, conta - 1
            temp(n)%Sharedmed = Shared%elem(n)%Sharedmed
            temp(n)%ProPmed = Shared%elem(n)%ProPmed
            temp(n)%field = Shared%elem(n)%field
            temp(n)%i = Shared%elem(n)%i
            temp(n)%j = Shared%elem(n)%j
            temp(n)%k = Shared%elem(n)%k
            temp(n)%times = Shared%elem(n)%times
         end do
         deallocate(Shared%elem)
         Shared%MaxConta = 2*Shared%MaxConta !!! 040717se atrancaba aqui. Ahora allocateo al doble. Antes era + 10000
        allocate(Shared%elem(1:Shared%MaxConta))
         do n = 1, conta - 1
            Shared%elem(n)%Sharedmed = temp(n)%Sharedmed
            Shared%elem(n)%ProPmed = temp(n)%ProPmed
            Shared%elem(n)%field = temp(n)%field
            Shared%elem(n)%i = temp(n)%i
            Shared%elem(n)%j = temp(n)%j
            Shared%elem(n)%k = temp(n)%k
            Shared%elem(n)%times = temp(n)%times
         end do
         deallocate(temp)
      end if
      !
      if (conta == 1) allocate (Shared%elem(1:Shared%MaxConta))
      Shared%elem(conta)%Sharedmed = Sharedmed
      Shared%elem(conta)%ProPmed = ProPmed
      Shared%elem(conta)%field = campo
      Shared%elem(conta)%i = i1
      Shared%elem(conta)%j = j1
      Shared%elem(conta)%k = k1
      Shared%elem(conta)%times = 2 !it appears two times at least !later updated in preporcess
      return
      !
   end subroutine
   !
end module
