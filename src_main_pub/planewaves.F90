module ilumina_m
   use FDETYPES_m
   use Report_m

   implicit none
   private

   real(kind=RKIND), allocatable, dimension(:,:,:) :: fpw
   real(kind=RKIND), allocatable, dimension(:,:) :: distanciaInicial,pxpw,pypw,pzpw,INCERT
   real(kind=RKIND), allocatable, dimension(:,:) :: evol
   real(kind=RKIND), allocatable, dimension(:) :: deltaevol
   integer(kind=4), allocatable, dimension(:) :: numus


   type ehxyz_t
      integer(kind=4) :: Ex=-15,Ey=-15,Ez=-15,Hx=-15,Hy=-15,Hz=-15
   end type
   type tfidaa_t
      type(ehxyz_t) :: com,fin,backDir,frontDir,leftDir,rightDir,downDir,arr
   end type
   type ijk_t
      type(tfidaa_t) :: i,j,k
   end type

   !!! global variables
   real(kind=RKIND) :: cluz,zvac
   real(kind=RKIND) :: eps0,mu0

   !!! local variables
   type(coorsxyzP_t) , save :: gridPoint
   type(ijk_t), allocatable, dimension(:)       , save  :: TrFr,IzDe,AbAr
   logical  , allocatable, dimension(:)        , save  :: IluminaTr,IluminaFr,IluminaIz,IluminaDe,IluminaAr,IluminaAb
   public Incid,AdvancePlaneWaveE,AdvancePlaneWaveH,InitPlaneWave,DestroyIlumina,storeplanewaves,calc_planewaveconstants,corrigeondaplanaH



contains
   subroutine InitPlaneWave(sgg,media,layoutnumber,num_procs,SINPML_Fullsize,ThereArePlaneWaveBoxes,resume,eps00,mu00)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(media_matrices_t), intent(in) :: media
      integer(kind=4), intent(in) :: layoutnumber,num_procs
      type(limit_t), dimension(1:6), intent(in) :: SINPML_fullsize
      integer j,k,field,i,jjj,maxnumus,maxmodes,kkk
      real(kind=RKIND) :: modulus,Xd0,Yd0,Zd0,diagonalcaja
      logical, intent(out) :: ThereArePlaneWaveBoxes
      logical  :: abortar, resume
      character(len=BUFSIZE) :: buff
      real(kind=RKIND), intent(in) :: eps00,mu00
      eps0=eps00; mu0=mu00; !hack to turn the step variables into globals
      cluz=1.0_RKIND/sqrt(eps0*mu0) !incid will need it
      zvac=sqrt(mu0/eps0) !the variables below need it

      do field=IEX,IHZ
         allocate (gridPoint%PhysCoor(field)%x(sgg%Sweep(field)%XI-1 : sgg%Sweep(field)%XE+1), &
         gridPoint%PhysCoor(field)%y(sgg%Sweep(field)%YI-1 : sgg%Sweep(field)%YE+1), &
         gridPoint%PhysCoor(field)%z(sgg%Sweep(field)%ZI-1 : sgg%Sweep(field)%ZE+1))
      end do

      field=IEX
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE
         gridPoint%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=sgg%Liney(j)
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=IEY
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE
         gridPoint%PhysCoor(field)%y(j)=(sgg%Liney(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=IEZ
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=sgg%Liney(j)
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE
         gridPoint%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=IHX
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE+1
         gridPoint%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE
         gridPoint%PhysCoor(field)%y(j)=(sgg%Liney(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE
         gridPoint%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=IHY
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE
         gridPoint%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE+1
         gridPoint%PhysCoor(field)%y(j)=sgg%Liney(j)
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE
         gridPoint%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=IHZ
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE
         gridPoint%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE
         gridPoint%PhysCoor(field)%y(j)=(sgg%Liney(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE+1
         gridPoint%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      !
      if (sgg%NumPlaneWaves >= 1) then
          TherearePlaneWaveBoxes=.true.
          continue
      else
          TherearePlaneWaveBoxes=.false.
          return
      end if
      if (ThereArePlaneWaveBoxes) then
         thereareplanewaveboxes=.false. !reset it because the MPI slice may not have any
         allocate (TrFr(1:sgg%numplanewaves), &
                   IzDe(1:sgg%numplanewaves), &
                   AbAr(1:sgg%numplanewaves), &
                   IluminaTr(1:sgg%numplanewaves), &
                   IluminaFr(1:sgg%numplanewaves), &
                   IluminaIz(1:sgg%numplanewaves), &
                   IluminaDe(1:sgg%numplanewaves), &
                   IluminaAr(1:sgg%numplanewaves), &
                   IluminaAb(1:sgg%numplanewaves), &
                   numus(1:sgg%numplanewaves), &
                   deltaevol(1:sgg%numplanewaves))
         do jjj=1,sgg%NumPlaneWaves
             numus(jjj)=sgg%PlaneWave(jjj)%sourceFile%NumSamples
             !by OLD's request I abort if there is nothing to illuminate
             abortar= &
             (sgg%PlaneWave(jjj)%esqx1 <=  SINPML_fullsize(IHX)%XI).and. &
             (sgg%PlaneWave(jjj)%esqx2 >=  SINPML_fullsize(IHX)%XE).and. &
             (sgg%PlaneWave(jjj)%esqy1 <=  SINPML_fullsize(IHY)%YI).and. &
             (sgg%PlaneWave(jjj)%esqy2 >=  SINPML_fullsize(IHY)%YE).and. &
             (sgg%PlaneWave(jjj)%esqz1 <=  SINPML_fullsize(IHZ)%ZI).and. &
             (sgg%PlaneWave(jjj)%esqz2 >=  SINPML_fullsize(IHZ)%ZE)
             if (abortar) then
                write (buff,'(a)') 'At least one of TF/SF planes must be 1 cell inside the simulation region. Aborting'
                call stoponerror(layoutnumber,num_procs,buff)
             end if
             !!!!!!!

             IluminaTr(jjj)=.false.
             IluminaFr(jjj)=.false.
             IluminaIz(jjj)=.false.
             IluminaDe(jjj)=.false.
             IluminaAr(jjj)=.false.
             IluminaAb(jjj)=.false.
             if ((sgg%PlaneWave(jjj)%esqx1 >= sgg%SINPMLSweep(IHX)%XI).and.(sgg%PlaneWave(jjj)%esqx1 <= sgg%SINPMLSweep(IHX)%XE)) &
             IluminaTr(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqx2 <= sgg%SINPMLSweep(IHX)%XE).and.(sgg%PlaneWave(jjj)%esqx2 >= sgg%SINPMLSweep(IHX)%XI)) &
             IluminaFr(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqy1 >= sgg%SINPMLSweep(IHY)%YI).and.(sgg%PlaneWave(jjj)%esqy1 <= sgg%SINPMLSweep(IHY)%YE)) &
             IluminaIz(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqy2 <= sgg%SINPMLSweep(IHY)%YE).and.(sgg%PlaneWave(jjj)%esqy2 >= sgg%SINPMLSweep(IHY)%YI)) &
             IluminaDe(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqz1 >= sgg%SINPMLSweep(IHZ)%ZI).and.(sgg%PlaneWave(jjj)%esqz1 <= sgg%SINPMLSweep(IHZ)%ZE)) &
             IluminaAb(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqz2 <= sgg%SINPMLSweep(IHZ)%ZE).and.(sgg%PlaneWave(jjj)%esqz2 >= sgg%SINPMLSweep(IHZ)%ZI)) &
             IluminaAr(jjj)=.true.
             !
             !find the coordinate limits of the Huygens Box for each component
             TrFr(jjj)%I%backDir%Ez=Max(sgg%SINPMLSweep(IEZ)%XI, sgg%PlaneWave(jjj)%esqx1)
             TrFr(jjj)%I%frontDir%Ez=Min(sgg%SINPMLSweep(IEZ)%XE, sgg%PlaneWave(jjj)%esqx2)
             TrFr(jjj)%J%com%Ez=Max(sgg%SINPMLSweep(IEZ)%YI, sgg%PlaneWave(jjj)%esqy1)
             TrFr(jjj)%J%fin%Ez=Min(sgg%SINPMLSweep(IEZ)%YE, sgg%PlaneWave(jjj)%esqy2)
             TrFr(jjj)%K%com%Ez=Max(sgg%SINPMLSweep(IEZ)%ZI, sgg%PlaneWave(jjj)%esqz1)
             TrFr(jjj)%K%fin%Ez=MIn(sgg%SINPMLSweep(IEZ)%ZE, sgg%PlaneWave(jjj)%esqz2-1)
             !
             TrFr(jjj)%I%backDir%Ey=Max(sgg%SINPMLSweep(IEY)%XI, sgg%PlaneWave(jjj)%esqx1)
             TrFr(jjj)%I%frontDir%Ey=Min(sgg%SINPMLSweep(IEY)%XE, sgg%PlaneWave(jjj)%esqx2)
             TrFr(jjj)%J%com%Ey=Max(sgg%SINPMLSweep(IEY)%YI, sgg%PlaneWave(jjj)%esqy1)
             TrFr(jjj)%J%fin%Ey=Min(sgg%SINPMLSweep(IEY)%YE, sgg%PlaneWave(jjj)%esqy2-1)
             TrFr(jjj)%K%com%Ey=Max(sgg%SINPMLSweep(IEY)%ZI ,sgg%PlaneWave(jjj)%esqz1)
             TrFr(jjj)%K%fin%Ey=MIn(sgg%SINPMLSweep(IEY)%ZE ,sgg%PlaneWave(jjj)%esqz2)
             !
             TrFr(jjj)%I%backDir%Hy= TrFr(jjj)%I%backDir%Ez - 1
             TrFr(jjj)%I%frontDir%Hy= TrFr(jjj)%I%frontDir%Ez
             TrFr(jjj)%J%com%Hy= TrFr(jjj)%J%com%Ez
             TrFr(jjj)%J%fin%Hy= TrFr(jjj)%J%fin%Ez
             TrFr(jjj)%K%com%Hy= TrFr(jjj)%K%com%Ez
             TrFr(jjj)%K%fin%Hy= TrFr(jjj)%K%fin%Ez
             !
             TrFr(jjj)%I%backDir%Hz= TrFr(jjj)%I%backDir%Ey - 1
             TrFr(jjj)%I%frontDir%Hz= TrFr(jjj)%I%frontDir%Ey
             TrFr(jjj)%J%com%Hz= TrFr(jjj)%J%com%Ey
             TrFr(jjj)%J%fin%Hz= TrFr(jjj)%J%fin%Ey
             TrFr(jjj)%K%com%Hz= TrFr(jjj)%K%com%Ey
             TrFr(jjj)%K%fin%Hz= TrFr(jjj)%K%fin%Ey
             !
             !
             IzDe(jjj)%J%leftDir%Ex=Max(sgg%SINPMLSweep(IEX)%yI, sgg%PlaneWave(jjj)%esqy1)
             IzDe(jjj)%J%rightDir%Ex=Min(sgg%SINPMLSweep(IEX)%yE, sgg%PlaneWave(jjj)%esqy2)
             IzDe(jjj)%I%com%Ex=Max(sgg%SINPMLSweep(IEX)%xI, sgg%PlaneWave(jjj)%esqx1)
             IzDe(jjj)%I%fin%Ex=Min(sgg%SINPMLSweep(IEX)%xE, sgg%PlaneWave(jjj)%esqx2-1)
             IzDe(jjj)%K%com%Ex=Max(sgg%SINPMLSweep(IEX)%ZI ,sgg%PlaneWave(jjj)%esqz1)
             IzDe(jjj)%K%fin%Ex=MIn(sgg%SINPMLSweep(IEX)%ZE ,sgg%PlaneWave(jjj)%esqz2)
             !
             IzDe(jjj)%J%leftDir%Ez=Max(sgg%SINPMLSweep(IEZ)%yI, sgg%PlaneWave(jjj)%esqy1)
             IzDe(jjj)%J%rightDir%Ez=Min(sgg%SINPMLSweep(IEZ)%yE, sgg%PlaneWave(jjj)%esqy2)
             IzDe(jjj)%I%com%Ez=Max(sgg%SINPMLSweep(IEZ)%xI, sgg%PlaneWave(jjj)%esqx1)
             IzDe(jjj)%I%fin%Ez=Min(sgg%SINPMLSweep(IEZ)%xE, sgg%PlaneWave(jjj)%esqx2)
             IzDe(jjj)%K%com%Ez=Max(sgg%SINPMLSweep(IEZ)%ZI ,sgg%PlaneWave(jjj)%esqz1)
             IzDe(jjj)%K%fin%Ez=MIn(sgg%SINPMLSweep(IEZ)%ZE ,sgg%PlaneWave(jjj)%esqz2-1)
             !
             IzDe(jjj)%J%leftDir%Hz= IzDe(jjj)%J%leftDir%Ex - 1
             IzDe(jjj)%J%rightDir%Hz= IzDe(jjj)%J%rightDir%Ex
             IzDe(jjj)%I%com%Hz= IzDe(jjj)%I%com%Ex
             IzDe(jjj)%I%fin%Hz= IzDe(jjj)%I%fin%Ex
             IzDe(jjj)%K%com%Hz= IzDe(jjj)%K%com%Ex
             IzDe(jjj)%K%fin%Hz= IzDe(jjj)%K%fin%Ex
             !
             IzDe(jjj)%J%leftDir%Hx= IzDe(jjj)%J%leftDir%Ez - 1
             IzDe(jjj)%J%rightDir%Hx= IzDe(jjj)%J%rightDir%Ez
             IzDe(jjj)%I%com%Hx= IzDe(jjj)%I%com%Ez
             IzDe(jjj)%I%fin%Hx= IzDe(jjj)%I%fin%Ez
             IzDe(jjj)%K%com%Hx= IzDe(jjj)%K%com%Ez
             IzDe(jjj)%K%fin%Hx= IzDe(jjj)%K%fin%Ez
             !
             !
             AbAr(jjj)%K%downDir%Ey=Max(sgg%SINPMLSweep(IEY)%ZI, sgg%PlaneWave(jjj)%esqz1)
             AbAr(jjj)%K%arr%Ey=Min(sgg%SINPMLSweep(IEY)%ZE, sgg%PlaneWave(jjj)%esqz2)
             AbAr(jjj)%I%com%Ey=Max(sgg%SINPMLSweep(IEY)%XI, sgg%PlaneWave(jjj)%esqx1)
             AbAr(jjj)%I%fin%Ey=Min(sgg%SINPMLSweep(IEY)%XE, sgg%PlaneWave(jjj)%esqx2)
             AbAr(jjj)%J%com%Ey=Max(sgg%SINPMLSweep(IEY)%YI, sgg%PlaneWave(jjj)%esqy1)
             AbAr(jjj)%J%fin%Ey=Min(sgg%SINPMLSweep(IEY)%YE, sgg%PlaneWave(jjj)%esqy2-1)
             !
             AbAr(jjj)%K%downDir%Ex=Max(sgg%SINPMLSweep(IEX)%ZI, sgg%PlaneWave(jjj)%esqz1)
             AbAr(jjj)%K%arr%Ex=Min(sgg%SINPMLSweep(IEX)%ZE, sgg%PlaneWave(jjj)%esqz2)
             AbAr(jjj)%I%com%Ex=Max(sgg%SINPMLSweep(IEX)%XI, sgg%PlaneWave(jjj)%esqx1)
             AbAr(jjj)%I%fin%Ex=Min(sgg%SINPMLSweep(IEX)%XE, sgg%PlaneWave(jjj)%esqx2-1)
             AbAr(jjj)%J%com%Ex=Max(sgg%SINPMLSweep(IEX)%YI, sgg%PlaneWave(jjj)%esqy1)
             AbAr(jjj)%J%fin%Ex=Min(sgg%SINPMLSweep(IEX)%YE, sgg%PlaneWave(jjj)%esqy2)
             !
             AbAr(jjj)%K%downDir%Hx= AbAr(jjj)%K%downDir%Ey - 1
             AbAr(jjj)%K%arr%Hx= AbAr(jjj)%K%arr%Ey
             AbAr(jjj)%I%com%Hx= AbAr(jjj)%I%com%Ey
             AbAr(jjj)%I%fin%Hx= AbAr(jjj)%I%fin%Ey
             AbAr(jjj)%J%com%Hx= AbAr(jjj)%J%com%Ey
             AbAr(jjj)%J%fin%Hx= AbAr(jjj)%J%fin%Ey
             !
             AbAr(jjj)%K%downDir%Hy= AbAr(jjj)%K%downDir%Ex -1
             AbAr(jjj)%K%arr%Hy= AbAr(jjj)%K%arr%Ex
             AbAr(jjj)%I%com%Hy= AbAr(jjj)%I%com%Ex
             AbAr(jjj)%I%fin%Hy= AbAr(jjj)%I%fin%Ex
             AbAr(jjj)%J%com%Hy= AbAr(jjj)%J%com%Ex
             AbAr(jjj)%J%fin%Hy= AbAr(jjj)%J%fin%Ex
             thereareplanewaveboxes=thereareplanewaveboxes.or.IluminaTr(jjj).or.IluminaFr(jjj).or.IluminaIz(jjj).or. &
             IluminaDe(jjj).or.IluminaAr(jjj).or.IluminaAb(jjj)

         end do !sweep j planewaves
      end if  !ThereArePlaneWaveBoxes

       maxnumus=maxval(numus)
       allocate (evol(1:sgg%numplanewaves,0 : maxnumus))
       do jjj=1,sgg%numplanewaves
           do k=0,numus(jjj)
              evol(jjj,k)=sgg%PlaneWave(jjj)%sourceFile%Samples(k)
           end do
           deltaevol(jjj)=sgg%PlaneWave(jjj)%sourceFile%deltaSamples
           if (deltaevol(jjj) > sgg%dt) then
              write (buff,'(a,e12.2e3)')  'WARNING: '//trim(adjustl(sgg%PlaneWave(jjj)%sourceFile%Name))// &
              ' undersampled by a factor ',deltaevol(jjj)/sgg%dt
              call WarnErrReport(buff)
           end if
       end do
!!
       maxmodes=maxval(sgg%PlaneWave(1:sgg%numplanewaves)%nummodes)
       allocate(pxpw(1:sgg%numplanewaves,maxmodes), &
                pypw(1:sgg%numplanewaves,maxmodes), &
                pzpw(1:sgg%numplanewaves,maxmodes), &
                fpw(1:sgg%numplanewaves,1:6,maxmodes), &
                INCERT(1:sgg%numplanewaves,maxmodes), &
                distanciaInicial(1:sgg%numplanewaves,maxmodes))
       do jjj=1,sgg%numplanewaves
         do kkk=1,sgg%PlaneWave(jjj)%nummodes
             if (.not.resume) then
                 pxpw(jjj,kkk)=sgg%PlaneWave(jjj)%px(kkk)
                 pypw(jjj,kkk)=sgg%PlaneWave(jjj)%py(kkk)
                 pzpw(jjj,kkk)=sgg%PlaneWave(jjj)%pz(kkk)
                 fpw(jjj,1,kkk)=sgg%PlaneWave(jjj)%ex(kkk)
                 fpw(jjj,2,kkk)=sgg%PlaneWave(jjj)%ey(kkk)
                 fpw(jjj,3,kkk)=sgg%PlaneWave(jjj)%ez(kkk)
!
                 modulus=sqrt(pxpw(jjj,kkk)**2+pypw(jjj,kkk)**2+pzpw(jjj,kkk)**2.0_RKIND)
                 pxpw(jjj,kkk)=pxpw(jjj,kkk)/modulus
                 pypw(jjj,kkk)=pypw(jjj,kkk)/modulus
                 pzpw(jjj,kkk)=pzpw(jjj,kkk)/modulus  
                 INCERT(jjj,kkk)=sgg%PlaneWave(jjj)%incert(kkk)
             else
                 if (sgg%PlaneWave(jjj)%isRC) then
                     read(14) pxpw(jjj,kkk),pypw(jjj,kkk),pzpw(jjj,kkk),fpw(jjj,1,kkk),fpw(jjj,2,kkk),fpw(jjj,3,kkk),INCERT(jjj,kkk)
                 else !initialize it as usual
                     pxpw(jjj,kkk)=sgg%PlaneWave(jjj)%px(kkk)
                     pypw(jjj,kkk)=sgg%PlaneWave(jjj)%py(kkk)
                     pzpw(jjj,kkk)=sgg%PlaneWave(jjj)%pz(kkk)
                     fpw(jjj,1,kkk)=sgg%PlaneWave(jjj)%ex(kkk)
                     fpw(jjj,2,kkk)=sgg%PlaneWave(jjj)%ey(kkk)
                     fpw(jjj,3,kkk)=sgg%PlaneWave(jjj)%ez(kkk)
    !
                     modulus=sqrt(pxpw(jjj,kkk)**2+pypw(jjj,kkk)**2+pzpw(jjj,kkk)**2.0_RKIND)
                     pxpw(jjj,kkk)=pxpw(jjj,kkk)/modulus
                     pypw(jjj,kkk)=pypw(jjj,kkk)/modulus
                     pzpw(jjj,kkk)=pzpw(jjj,kkk)/modulus  
                     INCERT(jjj,kkk)=sgg%PlaneWave(jjj)%incert(kkk)
                 end if
             end if
         end do
       end do
       do jjj=1,sgg%numplanewaves
         !Find the angles and amplitudes
         do kkk=1,sgg%PlaneWave(jjj)%nummodes
!!!! moved to where it is needed for permit scaling 081118
!!             fpw(jjj,4,kkk)=(pypw(jjj,kkk)*fpw(jjj,3,kkk)-pzpw(jjj,kkk)*fpw(jjj,2,kkk))/zvac
!!             fpw(jjj,5,kkk)=(pzpw(jjj,kkk)*fpw(jjj,1,kkk)-pxpw(jjj,kkk)*fpw(jjj,3,kkk))/zvac
!!             fpw(jjj,6,kkk)=(pxpw(jjj,kkk)*fpw(jjj,2,kkk)-pypw(jjj,kkk)*fpw(jjj,1,kkk))/zvac
             !
             !Find the null-phase corner depending on the angle of propagation
             if ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(IHX)%XI))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(IHY)%YI))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(IHZ)%ZI))
             else if ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(IHX)%XI))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(IHY)%YI))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(IHZ)%ZE))
             else if ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(IHX)%XI))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(IHY)%YE))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(IHZ)%ZI))
             else if ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(IHX)%XE))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(IHY)%YI))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(IHZ)%ZI))
             else if ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(IHX)%XI))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(IHY)%YE))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(IHZ)%ZE))
             else if ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(IHX)%XE))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(IHY)%YE))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(IHZ)%ZI))
             else if ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(IHX)%XE))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(IHY)%YI))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(IHZ)%ZE))
             else if ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(IHX)%XE))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(IHY)%YE))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(IHZ)%ZE))
             else
                call stoponerror(layoutnumber,num_procs,'buggy xo,yo,z0')
             end if
             diagonalcaja=sqrt( (sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(IHX)%XI)) - sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(IHX)%XE)))**2.0_RKIND  + &
                                (sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(IHY)%YI)) - sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(IHY)%YE)))**2.0_RKIND  + &
                                (sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(IHZ)%ZI)) - sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(IHZ)%ZE)))**2.0_RKIND  ) 
             distanciaInicial(jjj,kkk)=((XD0*pxpw(jjj,kkk)+YD0*pypw(jjj,kkk)+ZD0*pzpw(jjj,kkk)))-INCERT(jjj,kkk)*diagonalcaja !I THINK I HAVE TO SUBTRACT IT SO THAT THE UNCERTAINTY ONLY ENLARGES THE BOX (IT DELAYS THE SIGNAL)
                                                                           !!!! I confirm on 150419 that the uncertainty must be subtracted, after doubting about the sign, because then t-d/c, d=n.r-(n.r0-incert)>0 iff incert>0, since n.r-n.r0>0 always

             
         end do !of kkk
      end do !of maxmodes

      !check if materials are crossed by the box
      do jjj=1, sgg%numplanewaves
          if(IluminaTr(jjj)) then
             !Ez Back
             i = TrFr(jjj)%I%backDir%Ez !Back
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   if (media%sggMiEz(i, j, k) /=1) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Ez ',i,j,k
                      if (((media%sggMiEz(i,j,k) ==0).or.(sgg%med(media%sggMiEz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Ey Back
             i = TrFr(jjj)%I%backDir%Ey
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   if (media%sggMiEy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Ey ',i,j,k
                      if (((media%sggMiEy(i,j,k) ==0).or.(sgg%med(media%sggMiEy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !--->
          if(IluminaFr(jjj)) then
             !Ez  Front
             i = TrFr(jjj)%I%frontDir%Ez !Front
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   if (media%sggMiEz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Ez ',i,j,k
                      if (((media%sggMiEz(i,j,k) ==0).or.(sgg%med(media%sggMiEz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Ey  Front
             i = TrFr(jjj)%I%frontDir%Ey !Front
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   if (media%sggMiEy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Ey ',i,j,k
                      if (((media%sggMiEy(i,j,k) ==0).or.(sgg%med(media%sggMiEy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !--->
          if(IluminaIz(jjj)) then
             !Ex Left
             j = IzDe(jjj)%J%leftDir%Ex  !Left
             do k = IzDe(jjj)%K%com%Ex, IzDe(jjj)%K%fin%Ex
                do i = IzDe(jjj)%I%com%Ex, IzDe(jjj)%I%fin%Ex
                   if (media%sggMiEx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Ex ',i,j,k
                      if (((media%sggMiEx(i,j,k) ==0).or.(sgg%med(media%sggMiEx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Ez Left
             j = IzDe(jjj)%J%leftDir%Ez  !Left
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   if (media%sggMiEz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Ez ',i,j,k
                      if (((media%sggMiEz(i,j,k) ==0).or.(sgg%med(media%sggMiEz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !--->
          if(IluminaDe(jjj)) then
             !Ez  Right
             j = IzDe(jjj)%J%rightDir%Ez !Right
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   if (media%sggMiEz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Ez ',i,j,k
                      if (((media%sggMiEz(i,j,k) ==0).or.(sgg%med(media%sggMiEz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Ex  Right
             j = IzDe(jjj)%J%rightDir%Ex !Right
             do k = IzDe(jjj)%K%com%Ex,IzDe(jjj)%K%fin%Ex
                do i=IzDe(jjj)%I%com%Ex,IzDe(jjj)%I%fin%Ex
                   if (media%sggMiEx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Ex ',i,j,k
                      if (((media%sggMiEx(i,j,k) ==0).or.(sgg%med(media%sggMiEx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !--->
          if(IluminaAb(jjj)) then
             !Ex  Down
             k = AbAr(jjj)%K%downDir%Ex  !Down
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                do i=AbAr(jjj)%I%com%Ex,AbAr(jjj)%I%fin%Ex
                   if (media%sggMiEx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Ex ',i,j,k
                      if (((media%sggMiEx(i,j,k) ==0).or.(sgg%med(media%sggMiEx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Ey Down
             k = AbAr(jjj)%K%downDir%Ey  !Down
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   if (media%sggMiEy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Ey ',i,j,k
                      if (((media%sggMiEy(i,j,k) ==0).or.(sgg%med(media%sggMiEy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !--->
          if(IluminaAr(jjj)) then
             !Ex Up
             k = AbAr(jjj)%K%arr%Ex  !Up
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                do i = AbAr(jjj)%I%com%Ex, AbAr(jjj)%I%fin%Ex
                   if (media%sggMiEx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Ex ',i,j,k
                      if (((media%sggMiEx(i,j,k) ==0).or.(sgg%med(media%sggMiEx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Ey Up
             k = AbAr(jjj)%K%arr%Ey
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   if (media%sggMiEy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Ey ',i,j,k
                      if (((media%sggMiEy(i,j,k) ==0).or.(sgg%med(media%sggMiEy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !!!
          if(IluminaTr(jjj)) then
             !Hz Back
             i = TrFr(jjj)%I%backDir%Hz  !Back
             do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                   if (media%sggMiHz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Hz ',i,j,k
                      if (((media%sggMiHz(i,j,k) ==0).or.(sgg%med(media%sggMiHz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Hy Back
             i = TrFr(jjj)%I%backDir%Hy  !Back
             do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                   if (media%sggMiHy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Hy ',i,j,k
                      if (((media%sggMiHy(i,j,k) ==0).or.(sgg%med(media%sggMiHy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          if(IluminaFr(jjj)) then
             !Hz  Front
             i = TrFr(jjj)%I%frontDir%Hz !Front
             do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                   if (media%sggMiHz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Hz ',i,j,k
                      if (((media%sggMiHz(i,j,k) ==0).or.(sgg%med(media%sggMiHz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Hy  Front
             i = TrFr(jjj)%I%frontDir%Hy !Front
             do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                   if (media%sggMiHy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Hy ',i,j,k
                      if (((media%sggMiHy(i,j,k) ==0).or.(sgg%med(media%sggMiHy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do

          end if
          if(IluminaIz(jjj)) then
             !Hx Left
             j = IzDe(jjj)%J%leftDir%Hx  !Left
             do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                   if (media%sggMiHx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Hx ',i,j,k
                      if (((media%sggMiHx(i,j,k) ==0).or.(sgg%med(media%sggMiHx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Hz Left
             j = IzDe(jjj)%J%leftDir%Hz  !Left
             do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                   if (media%sggMiHz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Hz ',i,j,k
                      if (((media%sggMiHz(i,j,k) ==0).or.(sgg%med(media%sggMiHz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          if(IluminaDe(jjj)) then
             !Hx  Right
             j = IzDe(jjj)%J%rightDir%Hx !Right
             do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                   if (media%sggMiHx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Hx ',i,j,k
                      if (((media%sggMiHx(i,j,k) ==0).or.(sgg%med(media%sggMiHx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Hz  Right
             j = IzDe(jjj)%J%rightDir%Hz !Right
             do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                   if (media%sggMiHz(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Hz ',i,j,k
                      if (((media%sggMiHz(i,j,k) ==0).or.(sgg%med(media%sggMiHz(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          if(IluminaAb(jjj)) then
             !Hx  Down
             k = AbAr(jjj)%K%downDir%Hx  !Down
             do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                   if (media%sggMiHx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Hx ',i,j,k
                      if (((media%sggMiHx(i,j,k) ==0).or.(sgg%med(media%sggMiHx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Hy  Down
             k = AbAr(jjj)%K%downDir%Hy  !Down
             do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                do i=AbAr(jjj)%I%com%Hy,AbAr(jjj)%I%fin%Hy
                   if (media%sggMiHy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Hy ',i,j,k
                      if (((media%sggMiHy(i,j,k) ==0).or.(sgg%med(media%sggMiHy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
          !--->
          if(IluminaAr(jjj)) then
             !Hx Up
             k = AbAr(jjj)%K%arr%Hx  !Up
             do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                   if (media%sggMiHx(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Hx ',i,j,k
                      if (((media%sggMiHx(i,j,k) ==0).or.(sgg%med(media%sggMiHx(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
             !Hy Up
             k=AbAr(jjj)%K%arr%Hy  !Up
             do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                do i = AbAr(jjj)%I%com%Hy, AbAr(jjj)%I%fin%Hy
                   if (media%sggMiHy(i,j,k) /=1) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Hy ',i,j,k
                      if (((media%sggMiHy(i,j,k) ==0).or.(sgg%med(media%sggMiHy(i,j,k))%is%PEC)).and. .not. &
                      ((i == sgg%SINPMLSweep(IHX)%XI).or.(j == sgg%SINPMLSweep(IHY)%YI).or.(k == sgg%SINPMLSweep(IHZ)%ZI).or. &
                      (i == sgg%SINPMLSweep(IHX)%XE).or.(j == sgg%SINPMLSweep(IHY)%YE).or.(k == sgg%SINPMLSweep(IHZ)%ZE))) &
                      call stoponerror(layoutnumber,num_procs,buff)
                   end if
                end do
             end do
          end if
      end do !of j numplanewaves

!!!!
      call calc_planewaveconstants(sgg,eps0,mu0)
!!!
      return
   end subroutine InitPlaneWave




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Calculate the incident field at a given time/space point
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function Incid(sgg,jjj, nfield,time,i,j,k,still_planewave_time,calledfromobservation)    result(EHI)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      logical :: still_planewave_time,calledfromobservation
      integer(kind=4) i,j,k,nfield,jjj,kkk,jdum
      real(kind=RKIND) :: EHI
      real(kind=RKIND) :: time,d,xf,yf,zf
      !
      xf=gridPoint%PhysCoor(nfield)%x(i)
      yf=gridPoint%PhysCoor(nfield)%y(j)
      zf=gridPoint%PhysCoor(nfield)%z(k)
      ehi=0.0_RKIND

      if (calledfromobservation) then     
#ifdef CompileWithOpenMP
!$xMP   PARALLEL do DEFAULT(SHARED) private (d,kkk,jjj) REDUCTION(+:EhI)
#endif
            do jdum=1, sgg%numplanewaves !150419 observation must sum the planewaves; it has been moved here from the call
              do kkk=1,sgg%PlaneWave(jdum)%nummodes
                 d=(xf*pxpw(jdum,kkk)+yf*pypw(jdum,kkk)+zf*pzpw(jdum,kkk))-distanciaInicial(jdum,kkk)
                 EhI=EhI + fpw(jdum,nfield,kkk)*evolucion(jdum,time,d,still_planewave_time)
              !!!!!!!!!!!!!!!!!!!Ehi=Ehi*exp(-0.2*((-20 + i)**2.0_RKIND + (-20 + j)**2.0_RKIND ))
              end do
            end do
#ifdef CompileWithOpenMP
!$xMP   END PARALLEL DO
#endif
      else !if observation does not call it, jjj is already specified
#ifdef CompileWithOpenMP
!$xMP   PARALLEL do DEFAULT(SHARED) private (d,kkk,) REDUCTION(+:EhI)
#endif
              do kkk=1,sgg%PlaneWave(jjj)%nummodes
                 d=(xf*pxpw(jjj,kkk)+yf*pypw(jjj,kkk)+zf*pzpw(jjj,kkk))-distanciaInicial(jjj,kkk)
                 EhI=EhI + fpw(jjj,nfield,kkk)*evolucion(jjj,time,d,still_planewave_time)
              !!!!!!!!!!!!!!!!!!!Ehi=Ehi*exp(-0.2*((-20 + i)**2.0_RKIND + (-20 + j)**2.0_RKIND ))
              end do
#ifdef CompileWithOpenMP
!$xMP   END PARALLEL DO
#endif      
      end if
      return

   contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Evolution function to interpolate from the input file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real(kind=RKIND) function evolucion(jjj,t,d,still_planewave_time)
         real(kind=RKIND) t,d
         integer(kind=8) :: nprev
         integer(kind=4) :: jjj
         logical  :: still_planewave_time
!         if (d<=0.0_RKIND) then
!             print *,layr,' buggy error in d planewaves.evolucion. ' !beware because it slows things down. remove when we are sure about RC
!         end if

         evolucion=0.0_RKIND
         nprev=int((t-d/cluz)/deltaevol(jjj))
         if ((nprev+1 <= numus(jjj))) then 
           still_planewave_time=.true. !there may still be activity
           if (nprev > 0) then
            !first order interpolation
               evolucion=(evol(jjj,nprev+1)-evol(jjj,nprev))/deltaevol(jjj)*((t-d/cluz)-nprev*deltaevol(jjj))+evol(jjj,nprev) !linear interpolation
            !second order !no advantages over first order
            !  if (nprev+2 > numus(jjj)) then
            !      evolucion=0.0_RKIND !it is assumed that the input file contains an excitation that vanishes afterwards
            !  else
            !      evolucion=evol(jjj,nprev+2) * ( ((t-d/cluz)-nprev    *deltaevol(jjj)) * ((t-d/cluz)-(nprev+1)*deltaevol(jjj)) ) /(2.0_RKIND * deltaevol(jjj)**2.0_RKIND ) - &
            !                evol(jjj,nprev+1) * ( ((t-d/cluz)-nprev    *deltaevol(jjj)) * ((t-d/cluz)-(nprev+2)*deltaevol(jjj)) ) /(   deltaevol(jjj)**2.0_RKIND ) + &
            !                evol(jjj,nprev  ) * ( ((t-d/cluz)-(nprev+2)*deltaevol(jjj)) * ((t-d/cluz)-(nprev+1)*deltaevol(jjj)) ) /(2.0_RKIND * deltaevol(jjj)**2.0_RKIND )
            !  end if
           end if
         end if
         return
      end function evolucion
      !
   end function incid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyIlumina(sgg)
      type(SGGFDTDINFO_t), intent(inout) :: sgg
      integer(kind=4) :: field

      do field=IEX,IHZ
         if (associated(gridPoint%PhysCoor(field)%x)) deallocate(gridPoint%PhysCoor(field)%x)
         if (associated(gridPoint%PhysCoor(field)%y)) deallocate(gridPoint%PhysCoor(field)%y)
         if (associated(gridPoint%PhysCoor(field)%z)) deallocate(gridPoint%PhysCoor(field)%z)
      end do

      if (sgg%numplanewaves >=1) then
       deallocate(TrFr, IzDe,AbAr, IluminaTr, IluminaFr, IluminaIz,IluminaDe, IluminaAr,IluminaAb, pxpw, pypw, pzpw,   fpw, INCERT, numus,deltaevol,distanciainicial)
      end if
      if (allocated(evol)) deallocate(evol)
      if (associated(sgg%PlaneWave)) deallocate(sgg%PlaneWave)
   end subroutine DestroyIlumina



   !**************************************************************************************************
   subroutine AdvancePlaneWaveE(sgg, timeinstant, b, g2, Idxh, Idyh, Idzh, Ex, Ey, Ez,still_planewave_time)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      logical :: still_planewave_time
      logical :: called_fromobservation
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent(in) :: timeinstant
      !!!
      type(bounds_t), intent(in) :: b
      !--->
      real(kind = RKIND), dimension(0 :  sgg%NumMedia), intent(in) :: g2
      !--->
      real(kind = RKIND), dimension(0 :  b%dxh%NX-1), intent(in) :: Idxh
      real(kind = RKIND), dimension(0 :  b%dyh%NY-1), intent(in) :: Idyh
      real(kind = RKIND), dimension(0 :  b%dzh%NZ-1), intent(in) :: Idzh
      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent(inout) :: Ex
      real(kind = RKIND), dimension(0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent(inout) :: Ey
      real(kind = RKIND), dimension(0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent(inout) :: Ez
      !---------------------------> local variables <-----------------------------------------------
      real(kind = RKIND) :: timei, G2_1, Id,incidente
      integer  :: i, j, k, i_m, j_m, k_m,jjj
      character(len=BUFSIZE) :: dubuf
      !---------------------------> begins AdvancePlaneWaveE <---------------------------------------
!!!!

      !!!
      still_planewave_time=.false. !by default there will be no more plane wave activity, unless it goes through some non-trivial incid
      called_fromobservation=.false. !210419 
      
      timei = sgg%time(timeinstant)
      !!!! deprecated in pscale and the +3 of the sync with ORIGINAL is broken forever 110219 
      !!! timei = (timeinstant +3) * sgg%dt !ORIGINAL sync
      
      G2_1 = G2(1)
      !--->
      do jjj=1, sgg%numplanewaves
          if(IluminaTr(jjj)) then
             !Ez Back
             i = TrFr(jjj)%I%backDir%Ez !Back
             i_m = i - b%Ez%XI
             Id = Idxh(i_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   j_m = j - b%Ez%YI
                   !--->
                   incidente = Incid(sgg,jjj, IHY, timei, i-1, j, k,still_planewave_time,called_fromobservation)
                   Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) - G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey Back
             i = TrFr(jjj)%I%backDir%Ey  !Back
             i_m = i - b%Ey%XI
             Id = Idxh(i_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                k_m = k - b%Ey%ZI
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   j_m = j - b%Ey%YI
                   !--->
                   incidente = Incid(sgg,jjj, IHZ, timei, i-1, j, k,still_planewave_time,called_fromobservation)
                   Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) + G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          end if
          !--->
          if(IluminaFr(jjj)) then
             !Ez  Front
             i = TrFr(jjj)%I%frontDir%Ez !Front
             i_m = i - b%Ez%XI
             Id = Idxh(i_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   j_m = j - b%Ez%YI
                   !--->
                   incidente = Incid(sgg,jjj, IHY, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) + G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey  Front
             i = TrFr(jjj)%I%frontDir%Ey !Front
             i_m = i - b%Ey%XI
             Id = Idxh(i_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                k_m = k - b%Ey%ZI
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   j_m = j - b%Ey%YI
                   !--->
                   incidente = Incid(sgg,jjj, IHZ, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) - G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          end if
          !--->
          if(IluminaIz(jjj)) then
             !Ex Left
             j = IzDe(jjj)%J%leftDir%Ex  !Left
             j_m = j - b%Ex%YI
             Id = Idyh(j_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ex, IzDe(jjj)%K%fin%Ex
                k_m = k - b%Ex%ZI
                do i = IzDe(jjj)%I%com%Ex, IzDe(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHZ, timei, i, j-1, k,still_planewave_time,called_fromobservation)
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) - G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ez Left
             j = IzDe(jjj)%J%leftDir%Ez  !Left
             j_m = j - b%Ez%YI
             Id = Idyh(j_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   i_m = i - b%Ez%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHX, timei, i, j-1, k,still_planewave_time,called_fromobservation)
                   Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) + G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          end if
          !--->
          if(IluminaDe(jjj)) then
             !Ez  Right
             j = IzDe(jjj)%J%rightDir%Ez !Right
             j_m = j - b%Ez%YI
             Id = Idyh(j_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   i_m = i - b%Ez%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHX, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ez(i_m, j_m, k_m) = Ez(i_m, j_m, k_m) - G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ex  Right
             j = IzDe(jjj)%J%rightDir%Ex !Right
             j_m = j - b%Ex%YI
             Id = Idyh(j_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ex,IzDe(jjj)%K%fin%Ex
                k_m = k - b%Ex%ZI
                do i=IzDe(jjj)%I%com%Ex,IzDe(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHZ, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) + G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          end if
          !--->
          if(IluminaAb(jjj)) then
             !Ex  Down
             k = AbAr(jjj)%K%downDir%Ex  !Down
             k_m = k - b%Ex%ZI
             Id = Idzh(k_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                j_m = j - b%Ex%YI
                do i=AbAr(jjj)%I%com%Ex,AbAr(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHY, timei, i, j, k-1,still_planewave_time,called_fromobservation)
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) + G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey Down
             k = AbAr(jjj)%K%downDir%Ey  !Down
             k_m = k - b%Ey%ZI
             Id = Idzh(k_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                j_m = j - b%Ey%YI
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   i_m = i - b%Ey%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHX, timei, i, j, k-1,still_planewave_time,called_fromobservation)
                   Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) - G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          end if
          !--->
          if(IluminaAr(jjj)) then
             !Ex Up
             k = AbAr(jjj)%K%arr%Ex  !Up
             k_m = k - b%Ex%ZI
             Id = Idzh(k_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                j_m = j - b%Ex%YI
                do i = AbAr(jjj)%I%com%Ex, AbAr(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHY, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ex(i_m, j_m, k_m) = Ex(i_m, j_m, k_m) - G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey Up
             k = AbAr(jjj)%K%arr%Ey  !Up
             k_m = k - b%Ey%ZI
             Id = Idzh(k_m)
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                j_m = j - b%Ey%YI
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   i_m = i - b%Ey%XI
                   !--->
                   incidente = Incid(sgg,jjj, IHX, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ey(i_m, j_m, k_m) = Ey(i_m, j_m, k_m) + G2_1 * incidente * Id
                end do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          end if
      end do
      !---------------------------> ends AdvancePlaneWaveE <-----------------------------------------
      return
   end subroutine AdvancePlaneWaveE
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the H-field at n+0.5_RKIND
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**************************************************************************************************
   subroutine AdvancePlaneWaveH(sgg, timeinstant,  b, gm2, Idxe, Idye, Idze, Hx, Hy, Hz,still_planewave_time)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      logical :: still_planewave_time
      logical :: called_fromobservation
      
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent(in) :: timeinstant
      !!!
      type(bounds_t), intent(in) :: b
      !--->
      real(kind = RKIND), dimension(0 :  sgg%NumMedia), intent(in) :: gm2
      !--->
      real(kind = RKIND), dimension(0 :  b%dxe%NX-1), intent(in) :: Idxe
      real(kind = RKIND), dimension(0 :  b%dye%NY-1), intent(in) :: Idye
      real(kind = RKIND), dimension(0 :  b%dze%NZ-1), intent(in) :: Idze
      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(inout) :: Hx
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(inout) :: Hy
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(inout) :: Hz
      !---------------------------> local variables <-----------------------------------------------
      real(kind = RKIND) :: timei, Gm2_1, Id,incidente
      integer(kind=4) :: i, j, k, i_m, j_m, k_m,jjj
      character(len=BUFSIZE) :: dubuf
      !---------------------------> begins AdvancePlaneWaveH <---------------------------------------
      still_planewave_time=.false. !by default there will be no more plane wave activity, unless it goes through some non-trivial incid
      called_fromobservation=.false. !210419 
      !!!
      !!!
      
      timei = sgg%time(timeinstant) + 0.5_RKIND  * sgg%dt
      !!!! deprecated in pscale and the +3 of the sync with ORIGINAL is broken forever 110219 
      !!! timei = ( timeinstant + 0.5_RKIND  +3.0_RKIND) * sgg%dt  !ORIGINAL sync
      Gm2_1 = Gm2(1)
      !--->
     do jjj=1, sgg%numplanewaves
              if(IluminaTr(jjj)) then
                 !Hz Back
                 i = TrFr(jjj)%I%backDir%Hz  !Back
                 i_m = i - b%Hz%XI
                 Id = Idxe(i_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       !--->
                       incidente = Incid(sgg,jjj, IEY, timei, i+1, j, k,still_planewave_time,called_fromobservation)
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Back
                 i = TrFr(jjj)%I%backDir%Hy  !Back
                 i_m = i - b%Hy%XI
                 Id = Idxe(i_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       !--->
                       incidente = Incid(sgg,jjj,  IEZ, timei, i+1, j, k,still_planewave_time,called_fromobservation)
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaFr(jjj)) then
                 !Hz  Front
                 i = TrFr(jjj)%I%frontDir%Hz !Front
                 i_m = i - b%Hz%XI
                 Id = Idxe(i_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       !--->
                       incidente = Incid(sgg,jjj,  IEY, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Front
                 i = TrFr(jjj)%I%frontDir%Hy !Front
                 i_m = i - b%Hy%XI
                 Id = Idxe(i_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       !--->
                       incidente = Incid(sgg,jjj,  IEZ, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaIz(jjj)) then
                 !Hx Left
                 j = IzDe(jjj)%J%leftDir%Hx  !Left
                 j_m = j - b%Hx%YI
                 Id = Idye(j_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEZ, timei, i, j+1, k,still_planewave_time,called_fromobservation)
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz Left
                 j = IzDe(jjj)%J%leftDir%Hz  !Left
                 j_m = j - b%Hz%YI
                 Id = Idye(j_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEX, timei, i, j+1, k,still_planewave_time,called_fromobservation)
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaDe(jjj)) then
                 !Hx  Right
                 j = IzDe(jjj)%J%rightDir%Hx !Right
                 j_m = j - b%Hx%YI
                 Id = Idye(j_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEZ, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz  Right
                 j = IzDe(jjj)%J%rightDir%Hz !Right
                 j_m = j - b%Hz%YI
                 Id = Idye(j_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEX, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hz(i_m, j_m, k_m)=Hz(i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaAb(jjj)) then
                 !Hx  Down
                 k = AbAr(jjj)%K%downDir%Hx  !Down
                 k_m = k - b%Hx%ZI
                 Id = Idze(k_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEY, timei, i, j, k+1,still_planewave_time,called_fromobservation)
                       Hx(i_m, j_m, k_m)=Hx(i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Down
                 k = AbAr(jjj)%K%downDir%Hy  !Down
                 k_m = k - b%Hy%ZI
                 Id = Idze(k_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i=AbAr(jjj)%I%com%Hy,AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEX, timei, i, j, k+1,still_planewave_time,called_fromobservation)
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaAr(jjj)) then
                 !Hx Up
                 k = AbAr(jjj)%K%arr%Hx  !Up
                 k_m = k - b%Hx%ZI
                 Id = Idze(k_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEY, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Up
                 k=AbAr(jjj)%K%arr%Hy  !Up
                 k_m = k - b%Hy%ZI
                 Id = Idze(k_m)
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i = AbAr(jjj)%I%com%Hy, AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       !--->
                       incidente = Incid(sgg,jjj,  IEX, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
      end do 
      !---------------------------> ends AdvancePlaneWaveH <-----------------------------------------
      return
   end subroutine AdvancePlaneWaveH

    subroutine storeplanewaves(sgg)
      type(SGGFDTDINFO_t), intent(in) :: sgg
       integer(kind=4) :: jjj,kkk
       do jjj=1,sgg%numplanewaves
         do kkk=1,sgg%PlaneWave(jjj)%nummodes
            if (sgg%PlaneWave(jjj)%isRC) then
                 write(14,err=634) pxpw(jjj,kkk),pypw(jjj,kkk),pzpw(jjj,kkk),fpw(jjj,1,kkk),fpw(jjj,2,kkk),fpw(jjj,3,kkk),sgg%PlaneWave(jjj)%incert(kkk)
            end if
         end do
       end do
      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'PLANEWAVES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
    end subroutine storeplanewaves

    subroutine calc_planewaveconstants(sgg,eps00,mu00)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      real(kind = RKIND), intent(in) :: eps00,mu00
      integer :: jjj,kkk
      eps0=eps00; mu0=mu00; !hack to turn the step variables into globals
      cluz=1.0_RKIND/sqrt(eps0*mu0) !incid will need it
      zvac=sqrt(mu0/eps0) !the variables below need it
!!!!

      do jjj=1, sgg%numplanewaves
        do kkk=1,sgg%PlaneWave(jjj)%nummodes
             fpw(jjj,4,kkk)=(pypw(jjj,kkk)*fpw(jjj,3,kkk)-pzpw(jjj,kkk)*fpw(jjj,2,kkk))/zvac
             fpw(jjj,5,kkk)=(pzpw(jjj,kkk)*fpw(jjj,1,kkk)-pxpw(jjj,kkk)*fpw(jjj,3,kkk))/zvac
             fpw(jjj,6,kkk)=(pxpw(jjj,kkk)*fpw(jjj,2,kkk)-pypw(jjj,kkk)*fpw(jjj,1,kkk))/zvac
        end do
      end do
    end subroutine  calc_planewaveconstants

    
    subroutine corrigeondaplanaH(sgg,b,Hx,Hy,Hz,Hxvac, Hyvac, Hzvac)
      !!!
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(bounds_t), intent(in) :: b
      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(inout) :: Hx
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(inout) :: Hy
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(inout) :: Hz
      !---------------------------> local variables <-----------------------------------------------
      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(inout) :: Hxvac
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(inout) :: Hyvac
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(inout) :: Hzvac
      !---------------------------> local variables <-----------------------------------------------
      integer(kind=4) :: i, j, k, i_m, j_m, k_m,jjj

      do jjj=1, sgg%numplanewaves
              if(IluminaTr(jjj)) then
                 !Hz Back
                 i = TrFr(jjj)%I%backDir%Hz  !Back
                 i_m = i - b%Hz%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       !--->
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Hzvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Back
                 i = TrFr(jjj)%I%backDir%Hy  !Back
                 i_m = i - b%Hy%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Hyvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaFr(jjj)) then
                 !Hz  Front
                 i = TrFr(jjj)%I%frontDir%Hz !Front
                 i_m = i - b%Hz%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Hzvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Front
                 i = TrFr(jjj)%I%frontDir%Hy !Front
                 i_m = i - b%Hy%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Hyvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaIz(jjj)) then
                 !Hx Left
                 j = IzDe(jjj)%J%leftDir%Hx  !Left
                 j_m = j - b%Hx%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Hxvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz Left
                 j = IzDe(jjj)%J%leftDir%Hz  !Left
                 j_m = j - b%Hz%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       Hz(i_m, j_m, k_m) = Hz(i_m, j_m, k_m) - Hzvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaDe(jjj)) then
                 !Hx  Right
                 j = IzDe(jjj)%J%rightDir%Hx !Right
                 j_m = j - b%Hx%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Hxvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz  Right
                 j = IzDe(jjj)%J%rightDir%Hz !Right
                 j_m = j - b%Hz%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       Hz(i_m, j_m, k_m)=Hz(i_m, j_m, k_m) - Hzvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaAb(jjj)) then
                 !Hx  Down
                 k = AbAr(jjj)%K%downDir%Hx  !Down
                 k_m = k - b%Hx%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx(i_m, j_m, k_m)=Hx(i_m, j_m, k_m) - Hxvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Down
                 k = AbAr(jjj)%K%downDir%Hy  !Down
                 k_m = k - b%Hy%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i=AbAr(jjj)%I%com%Hy,AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Hyvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
              !--->
              if(IluminaAr(jjj)) then
                 !Hx Up
                 k = AbAr(jjj)%K%arr%Hx  !Up
                 k_m = k - b%Hx%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx(i_m, j_m, k_m) = Hx(i_m, j_m, k_m) - Hxvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Up
                 k=AbAr(jjj)%K%arr%Hy  !Up
                 k_m = k - b%Hy%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i = AbAr(jjj)%I%com%Hy, AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       Hy(i_m, j_m, k_m) = Hy(i_m, j_m, k_m) - Hyvac(i_m, j_m, k_m)
                    end do
                 end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              end if
      end do     
      
      
      
    return
    end subroutine corrigeondaplanaH
    

end module ilumina_m
