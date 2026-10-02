!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module for the Maloney-Smith thin material sheet subcell model.
!
! Reference:
!   J. G. Maloney and G. S. Smith, "The efficient modeling of thin material
!   sheets in the finite-difference time-domain (FDTD) method,"
!   IEEE Trans. Antennas Propag., vol. 40, no. 3, pp. 323-330, Mar. 1992.
!
! Materials of type maloneySheet are assigned to oriented cell faces. The sheet
! is centered on the face and is split evenly between the two adjacent cells.
! For every special cell (the cells adjacent to the sheet) one additional
! normal electric field unknown Eintern is advanced with the sheet constitutive
! parameters (eq. 8 of the reference). The ordinary normal field Efield is
! advanced by the main FDTD update (eq. 7), and the tangential magnetic field
! updates are corrected by the weighted difference between Eintern and Efield
! with the fraction of the cell occupied by the sheet (eqs. 15-16).
! The tangential electric fields use an effective medium over the face
! (eqs. 10-13), created during preprocessing.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module maloney_thin_sheet_m

   use Report_m

   use FDETYPES_m
   use SGBC_m, only: g1g2

   implicit none
   private


   public :: maloneySheetNode_t
   public :: maloneySheets_t
   public :: InitMaloneySheets
   public :: AdvanceMaloneySheetE
   public :: AdvanceMaloneySheetH
   public :: DestroyMaloneySheets
   public :: StoreFieldsMaloneySheets
   public :: GetMaloneySheets


   ! One special cell of the thin sheet model.
   type :: maloneySheetNode_t
      integer(kind=4) :: orientation = 0
      real(kind=RKIND), pointer :: Efield => null ()     ! ordinary normal field
      real(kind=RKIND) :: Eintern = 0.0_RKIND            ! in-sheet normal field
      real(kind=RKIND) :: beta = 0.0_RKIND               ! sheet fraction in the cell
      real(kind=RKIND) :: g1 = 1.0_RKIND                 ! in-sheet update coefficients
      real(kind=RKIND) :: g2 = 0.0_RKIND
      real(kind=RKIND), pointer :: H1 => null ()
      real(kind=RKIND), pointer :: H2 => null ()
      real(kind=RKIND), pointer :: H3 => null ()
      real(kind=RKIND), pointer :: H4 => null ()
      real(kind=RKIND) :: c1 = 0.0_RKIND                 ! magnetic field corrections
      real(kind=RKIND) :: c2 = 0.0_RKIND
      real(kind=RKIND) :: c3 = 0.0_RKIND
      real(kind=RKIND) :: c4 = 0.0_RKIND
      real(kind=RKIND) :: deltaA = 0.0_RKIND             ! stencil spacings for Eintern
      real(kind=RKIND) :: deltaB = 0.0_RKIND
   end type maloneySheetNode_t


   type :: maloneySheets_t
      integer(kind=4) :: numNodes = 0
      type(maloneySheetNode_t), allocatable :: nodes(:)
   end type maloneySheets_t


   type(maloneySheets_t), save, target :: maloneySheets


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initializes the list of special cells from the media flagged as maloneySheet.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMaloneySheets(sgg, media, Ex, Ey, Ez, Hx, Hy, Hz, &
                                IDxe, IDye, IDze, IDxh, IDyh, IDzh, &
                                g, eps00, ThereAreMaloneySheets, resume)
      type(media_matrices_t), intent(in) :: media
      type(constants_t), intent(in) :: g
      real(kind=RKIND), intent(in) :: eps00
      logical, intent(out) :: ThereAreMaloneySheets
      logical, intent(in) :: resume
      type(SGGFDTDINFO_t), intent(inout) :: sgg
      real(kind=RKIND), intent(in), target :: &
         Ex(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE, sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE, sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), &
         Ey(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE, sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE, sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), &
         Ez(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE, sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE, sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), &
         Hx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE, sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE, sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE), &
         Hy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE, sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE, sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE), &
         Hz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE, sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE, sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      real(kind=RKIND), intent(in) :: &
         IDxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
         IDye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
         IDze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
         IDxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
         IDyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
         IDzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)

      integer(kind=4) :: jmed, conta, i1, j1, k1, cidx, l, maxNodes
      integer(kind=4) :: orientacion
      real(kind=RKIND) :: eps, sigma, g1s, g2s, delta, gm2v
      real(kind=RKIND) :: thk, epr
      character(len=BUFSIZE) :: buff
      logical, allocatable, dimension(:,:,:) :: addedEx, addedEy, addedEz

      ThereAreMaloneySheets = .false.
      do jmed = 1, sgg%NumMedia
         if (sgg%Med(jmed)%Is%MaloneySheet) then
            ThereAreMaloneySheets = .true.
         end if
      end do
      if (.not. ThereAreMaloneySheets) then
         return
      end if

      ! Upper bound: two special cells per tangential E node on a sheet face.
      maxNodes = 0
      do k1 = sgg%SINPMLSweep(iEx)%ZI, sgg%SINPMLSweep(iEx)%ZE
         do j1 = sgg%SINPMLSweep(iEx)%YI, sgg%SINPMLSweep(iEx)%YE
            do i1 = sgg%SINPMLSweep(iEx)%XI, sgg%SINPMLSweep(iEx)%XE
               jmed = media%sggMiEx(i1,j1,k1)
               if (sgg%Med(jmed)%Is%MaloneySheet) maxNodes = maxNodes + 2
            end do
         end do
      end do
      do k1 = sgg%SINPMLSweep(iEy)%ZI, sgg%SINPMLSweep(iEy)%ZE
         do j1 = sgg%SINPMLSweep(iEy)%YI, sgg%SINPMLSweep(iEy)%YE
            do i1 = sgg%SINPMLSweep(iEy)%XI, sgg%SINPMLSweep(iEy)%XE
               jmed = media%sggMiEy(i1,j1,k1)
               if (sgg%Med(jmed)%Is%MaloneySheet) maxNodes = maxNodes + 2
            end do
         end do
      end do
      do k1 = sgg%SINPMLSweep(iEz)%ZI, sgg%SINPMLSweep(iEz)%ZE
         do j1 = sgg%SINPMLSweep(iEz)%YI, sgg%SINPMLSweep(iEz)%YE
            do i1 = sgg%SINPMLSweep(iEz)%XI, sgg%SINPMLSweep(iEz)%XE
               jmed = media%sggMiEz(i1,j1,k1)
               if (sgg%Med(jmed)%Is%MaloneySheet) maxNodes = maxNodes + 2
            end do
         end do
      end do

      maloneySheets%numNodes = 0
      if (maxNodes == 0) then
         ThereAreMaloneySheets = .false.
         return
      end if
      allocate (maloneySheets%nodes(1 : maxNodes))

      allocate (addedEx(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE, sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE, &
                        sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), source=.false.)
      allocate (addedEy(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE, sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE, &
                        sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), source=.false.)
      allocate (addedEz(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE, sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE, &
                        sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), source=.false.)

      ! Tangential E nodes lying on a sheet face. Each one yields the two
      ! special cells adjacent to the face in the normal direction.
      do k1 = sgg%SINPMLSweep(iEx)%ZI, sgg%SINPMLSweep(iEx)%ZE
         do j1 = sgg%SINPMLSweep(iEx)%YI, sgg%SINPMLSweep(iEx)%YE
            do i1 = sgg%SINPMLSweep(iEx)%XI, sgg%SINPMLSweep(iEx)%XE
               jmed = media%sggMiEx(i1,j1,k1)
               if (.not.sgg%Med(jmed)%Is%MaloneySheet) cycle
               orientacion = abs(sgg%Med(jmed)%Multiport(1)%Multiportdir)
               if (orientacion == iEx) cycle
               if (orientacion == iEy) then
                  call addSheetCell(orientacion, i1, j1,   k1, jmed)
                  call addSheetCell(orientacion, i1, j1-1, k1, jmed)
               else if (orientacion == iEz) then
                  call addSheetCell(orientacion, i1, j1, k1,   jmed)
                  call addSheetCell(orientacion, i1, j1, k1-1, jmed)
               end if
            end do
         end do
      end do
      do k1 = sgg%SINPMLSweep(iEy)%ZI, sgg%SINPMLSweep(iEy)%ZE
         do j1 = sgg%SINPMLSweep(iEy)%YI, sgg%SINPMLSweep(iEy)%YE
            do i1 = sgg%SINPMLSweep(iEy)%XI, sgg%SINPMLSweep(iEy)%XE
               jmed = media%sggMiEy(i1,j1,k1)
               if (.not.sgg%Med(jmed)%Is%MaloneySheet) cycle
               orientacion = abs(sgg%Med(jmed)%Multiport(1)%Multiportdir)
               if (orientacion == iEy) cycle
               if (orientacion == iEx) then
                  call addSheetCell(orientacion, i1,   j1, k1, jmed)
                  call addSheetCell(orientacion, i1-1, j1, k1, jmed)
               else if (orientacion == iEz) then
                  call addSheetCell(orientacion, i1, j1, k1,   jmed)
                  call addSheetCell(orientacion, i1, j1, k1-1, jmed)
               end if
            end do
         end do
      end do
      do k1 = sgg%SINPMLSweep(iEz)%ZI, sgg%SINPMLSweep(iEz)%ZE
         do j1 = sgg%SINPMLSweep(iEz)%YI, sgg%SINPMLSweep(iEz)%YE
            do i1 = sgg%SINPMLSweep(iEz)%XI, sgg%SINPMLSweep(iEz)%XE
               jmed = media%sggMiEz(i1,j1,k1)
               if (.not.sgg%Med(jmed)%Is%MaloneySheet) cycle
               orientacion = abs(sgg%Med(jmed)%Multiport(1)%Multiportdir)
               if (orientacion == iEz) cycle
               if (orientacion == iEx) then
                  call addSheetCell(orientacion, i1,   j1, k1, jmed)
                  call addSheetCell(orientacion, i1-1, j1, k1, jmed)
               else if (orientacion == iEy) then
                  call addSheetCell(orientacion, i1, j1,   k1, jmed)
                  call addSheetCell(orientacion, i1, j1-1, k1, jmed)
               end if
            end do
         end do
      end do

      deallocate (addedEx, addedEy, addedEz)

      if (.not.resume) then
         do conta = 1, maloneySheets%numNodes
            maloneySheets%nodes(conta)%Eintern = 0.0_RKIND
         end do
      else
         do conta = 1, maloneySheets%numNodes
            read (14) maloneySheets%nodes(conta)%Eintern
         end do
      end if

      write (buff, *) ' Maximum maloneySheet nodes= ', maloneySheets%numNodes
      call WarnErrReport(buff)

      return

   contains

      ! Adds the special cell (i,j,k), where the index in the direction of
      ! orientacion is a cell index and the other two are node indices.
      subroutine addSheetCell(or, i, j, k, jmed)
         integer(kind=4), intent(in) :: or, i, j, k, jmed
         integer(kind=4) :: ic, jc, kc
         real(kind=RKIND) :: gm2H1, gm2H2, gm2H3, gm2H4

         ic = i
         jc = j
         kc = k
         select case (or)
          case (iEx)
            if (.not. owned(iEx, ic, j, k)) return
            if (addedEx(ic,j,k)) return
            addedEx(ic,j,k) = .true.
          case (iEy)
            if (.not. owned(iEy, i, jc, k)) return
            if (addedEy(i,jc,k)) return
            addedEy(i,jc,k) = .true.
          case (iEz)
            if (.not. owned(iEz, i, j, kc)) return
            if (addedEz(i,j,kc)) return
            addedEz(i,j,kc) = .true.
         end select

         thk = sgg%Med(jmed)%Multiport(1)%width(1)
         epr = sgg%Med(jmed)%Multiport(1)%epr(1)
         sigma = sgg%Med(jmed)%Multiport(1)%sigma(1)

         maloneySheets%numNodes = maloneySheets%numNodes + 1
         l = maloneySheets%numNodes
         maloneySheets%nodes(l)%orientation = or

         select case (or)
          case (iEx)
            maloneySheets%nodes(l)%Efield => Ex(ic,j,k)
            delta = sgg%DX(ic)
            maloneySheets%nodes(l)%deltaA = IDyh(j)
            maloneySheets%nodes(l)%deltaB = IDzh(k)
            maloneySheets%nodes(l)%H1 => Hz(ic,j,k)
            maloneySheets%nodes(l)%H2 => Hz(ic,j-1,k)
            maloneySheets%nodes(l)%H3 => Hy(ic,j,k)
            maloneySheets%nodes(l)%H4 => Hy(ic,j,k-1)
            gm2H1 = g%gm2(media%sggMiHz(ic,j,k))
            gm2H2 = g%gm2(media%sggMiHz(ic,j-1,k))
            gm2H3 = g%gm2(media%sggMiHy(ic,j,k))
            gm2H4 = g%gm2(media%sggMiHy(ic,j,k-1))
            maloneySheets%nodes(l)%c1 = -gm2H1 * IDye(j)
            maloneySheets%nodes(l)%c2 = +gm2H2 * IDye(j-1)
            maloneySheets%nodes(l)%c3 = +gm2H3 * IDze(k)
            maloneySheets%nodes(l)%c4 = -gm2H4 * IDze(k-1)
          case (iEy)
            maloneySheets%nodes(l)%Efield => Ey(i,jc,k)
            delta = sgg%DY(jc)
            maloneySheets%nodes(l)%deltaA = IDzh(k)
            maloneySheets%nodes(l)%deltaB = IDxh(i)
            maloneySheets%nodes(l)%H1 => Hx(i,jc,k)
            maloneySheets%nodes(l)%H2 => Hx(i,jc,k-1)
            maloneySheets%nodes(l)%H3 => Hz(i,jc,k)
            maloneySheets%nodes(l)%H4 => Hz(i-1,jc,k)
            gm2H1 = g%gm2(media%sggMiHx(i,jc,k))
            gm2H2 = g%gm2(media%sggMiHx(i,jc,k-1))
            gm2H3 = g%gm2(media%sggMiHz(i,jc,k))
            gm2H4 = g%gm2(media%sggMiHz(i-1,jc,k))
            maloneySheets%nodes(l)%c1 = -gm2H1 * IDze(k)
            maloneySheets%nodes(l)%c2 = +gm2H2 * IDze(k-1)
            maloneySheets%nodes(l)%c3 = +gm2H3 * IDxe(i)
            maloneySheets%nodes(l)%c4 = -gm2H4 * IDxe(i-1)
          case (iEz)
            maloneySheets%nodes(l)%Efield => Ez(i,j,kc)
            delta = sgg%DZ(kc)
            maloneySheets%nodes(l)%deltaA = IDxh(i)
            maloneySheets%nodes(l)%deltaB = IDyh(j)
            maloneySheets%nodes(l)%H1 => Hy(i,j,kc)
            maloneySheets%nodes(l)%H2 => Hy(i-1,j,kc)
            maloneySheets%nodes(l)%H3 => Hx(i,j,kc)
            maloneySheets%nodes(l)%H4 => Hx(i,j-1,kc)
            gm2H1 = g%gm2(media%sggMiHy(i,j,kc))
            gm2H2 = g%gm2(media%sggMiHy(i-1,j,kc))
            gm2H3 = g%gm2(media%sggMiHx(i,j,kc))
            gm2H4 = g%gm2(media%sggMiHx(i,j-1,kc))
            maloneySheets%nodes(l)%c1 = -gm2H1 * IDxe(i)
            maloneySheets%nodes(l)%c2 = +gm2H2 * IDxe(i-1)
            maloneySheets%nodes(l)%c3 = +gm2H3 * IDye(j)
            maloneySheets%nodes(l)%c4 = -gm2H4 * IDye(j-1)
         end select

         maloneySheets%nodes(l)%beta = 0.5_RKIND * thk / delta
         if (maloneySheets%nodes(l)%beta > 1.0_RKIND) then
            write (buff, *) 'ERROR: maloneySheet is thicker than twice the cell: ', thk, delta
            call StopOnError(0, 0, buff)
         end if

         eps = epr * eps00
         call g1g2(sgg%dt, eps, sigma, g1s, g2s)
         maloneySheets%nodes(l)%g1 = g1s
         maloneySheets%nodes(l)%g2 = g2s

      end subroutine addSheetCell

      logical function owned(or, i, j, k) result(res)
         integer(kind=4), intent(in) :: or, i, j, k
         res = (i >= sgg%SINPMLSweep(or)%XI).and.(i <= sgg%SINPMLSweep(or)%XE).and. &
               (j >= sgg%SINPMLSweep(or)%YI).and.(j <= sgg%SINPMLSweep(or)%YE).and. &
               (k >= sgg%SINPMLSweep(or)%ZI).and.(k <= sgg%SINPMLSweep(or)%ZE)
      end function owned

   end subroutine InitMaloneySheets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Advances the in-sheet normal electric field of each special cell.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine AdvanceMaloneySheetE()
      integer(kind=4) :: conta
      real(kind=RKIND) :: curl

      if (maloneySheets%numNodes == 0) return
#ifdef CompileWithOpenMP
!$OMP  PARALLEL do DEFAULT(SHARED) private (conta,curl) schedule(guided)
#endif
      do conta = 1, maloneySheets%numNodes
         associate (node => maloneySheets%nodes(conta))
            curl = (node%H1 - node%H2) * node%deltaA - (node%H3 - node%H4) * node%deltaB
            node%Eintern = node%g1 * node%Eintern + node%g2 * curl
         end associate
      end do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
   end subroutine AdvanceMaloneySheetE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Corrects the tangential magnetic fields with the weighted difference between
! the in-sheet and the ordinary normal electric fields.
! This loop is not parallelized because several nodes can share H components.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine AdvanceMaloneySheetH()
      integer(kind=4) :: conta
      real(kind=RKIND) :: dE

      if (maloneySheets%numNodes == 0) return
      do conta = 1, maloneySheets%numNodes
         associate (node => maloneySheets%nodes(conta))
            dE = node%beta * (node%Eintern - node%Efield)
            node%H1 = node%H1 + node%c1 * dE
            node%H2 = node%H2 + node%c2 * dE
            node%H3 = node%H3 + node%c3 * dE
            node%H4 = node%H4 + node%c4 * dE
         end associate
      end do
   end subroutine AdvanceMaloneySheetH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine StoreFieldsMaloneySheets()
      integer(kind=4) :: conta

      do conta = 1, maloneySheets%numNodes
         write (14, err=634) maloneySheets%nodes(conta)%Eintern
      end do
      goto 635
634   call print11(0, SEPARADOR//separador//separador)
      call print11(0, 'maloneySheet: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0, SEPARADOR//separador//separador)
635   return
   end subroutine StoreFieldsMaloneySheets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyMaloneySheets(sgg)
      type(SGGFDTDINFO_t), intent(in) :: sgg

      if (allocated(maloneySheets%nodes)) deallocate (maloneySheets%nodes)
      maloneySheets%numNodes = 0
   end subroutine DestroyMaloneySheets

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function GetMaloneySheets() result(r)
      type(maloneySheets_t), pointer :: r

      r => maloneySheets
   end function GetMaloneySheets

end module maloney_thin_sheet_m
