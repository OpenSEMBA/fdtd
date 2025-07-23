
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Plane Wave Feeding modules
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module ilumina
   use fdetypes
   USE REPORT

   IMPLICIT NONE
   private

   REAL (KIND=RKIND), allocatable, dimension (:,:,:)    :: fpw
   REAL (KIND=RKIND), allocatable, dimension(:,:)       :: distanciaInicial,pxpw,pypw,pzpw,INCERT
   REAL (KIND=RKIND), allocatable, dimension ( :,: )  ::  evol
   REAL (KIND=RKIND), allocatable, dimension ( : )  ::  deltaevol
   integer (kind=4), allocatable, dimension(:)        ::  numus


   type ehxyz
      integer (kind=4)  ::  Ex=-15,Ey=-15,Ez=-15,Hx=-15,Hy=-15,Hz=-15
   end type
   type tfidaa
      type (ehxyz)  ::   com,fin,tra,fro,izq,der,aba,arr
   end type
   type ijk
      type (tfidaa)  ::  i,j,k
   end type

!!!variables globales del modulo
      REAL (KIND=RKIND)           ::  cluz,zvac
      REAL (KIND=RKIND)           ::  eps0,mu0
!!!
   !!!!local variables
   type (coorsxyzP) , save ::   Punto
   type (ijk), allocatable, dimension(:)       , SAVE  ::  TrFr,IzDe,AbAr
   logical  , allocatable, dimension(:)        , save  ::  IluminaTr,IluminaFr,IluminaIz,IluminaDe,IluminaAr,IluminaAb
   public Incid,AdvancePlaneWaveE,AdvancePlaneWaveH,InitPlaneWave,DestroyIlumina,storeplanewaves,calc_planewaveconstants,corrigeondaplanaH



contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes Plane Wave data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitPlaneWave(sgg,media,layoutnumber,size,SINPML_Fullsize,ThereArePlaneWaveBoxes,resume,eps00,mu00)
   ! subroutine InitPlaneWave(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,layoutnumber,size,SINPML_Fullsize,ThereArePlaneWaveBoxes,resume,eps00,mu00)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      type(media_matrices_t), intent(in) :: media
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES)  ::  &
      sggMiEx(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), &
      sggMiEy(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), &
      sggMiEz(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), &
      sggMiHx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE), &
      sggMiHy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE), &
      sggMiHz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      !!!
      integer (kind=4), intent(in) :: layoutnumber,size
      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      integer j,k,field,i,jjj,maxnumus,maxmodes,kkk
      REAL (KIND=RKIND) :: modulus,Xd0,Yd0,Zd0,diagonalcaja
      logical, intent(out)  ::  ThereArePlaneWaveBoxes
      logical  ::  abortar, resume
      character(len=BUFSIZE) :: buff
      REAL (KIND=RKIND), intent(in)   :: eps00,mu00
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      cluz=1.0_RKIND/sqrt(eps0*mu0) !lo necesitara incid
      zvac=sqrt(mu0/eps0) !lo necesitan las variables de mas abajo

      sggMiEx = media%sggMiEx
      sggMiEy = media%sggMiEy
      sggMiEz = media%sggMiEz
      sggMiHx = media%sggMiHx
      sggMiHy = media%sggMiHy
      sggMiHz = media%sggMiHz

!!!
      do field=iEx,iHz
         allocate (Punto%PhysCoor(field)%x(sgg%Sweep(field)%XI-1 : sgg%Sweep(field)%XE+1), &
         Punto%PhysCoor(field)%y(sgg%Sweep(field)%YI-1 : sgg%Sweep(field)%YE+1), &
         Punto%PhysCoor(field)%z(sgg%Sweep(field)%ZI-1 : sgg%Sweep(field)%ZE+1))
      end do

      field=iEx
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE
         Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=sgg%Liney(j)
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=iEy
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE
         Punto%PhysCoor(field)%y(j)=(sgg%Liney(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=iEz
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=sgg%Liney(j)
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE
         Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHx
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE+1
         Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE
         Punto%PhysCoor(field)%y(j)=(sgg%Liney(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE
         Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHy
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE
         Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE+1
         Punto%PhysCoor(field)%y(j)=sgg%Liney(j)
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE
         Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHz
      do i=sgg%Sweep(field)%XI-1,sgg%Sweep(field)%XE
         Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=sgg%Sweep(field)%YI-1,sgg%Sweep(field)%YE
         Punto%PhysCoor(field)%y(j)=(sgg%Liney(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=sgg%Sweep(field)%ZI-1,sgg%Sweep(field)%ZE+1
         Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      !
      if (sgg%NumPlaneWaves >= 1) then
          TherearePlaneWaveBoxes=.true.
          continue
      else
          TherearePlaneWaveBoxes=.false.
          return
      endif
      If (ThereArePlaneWaveBoxes) then
         thereareplanewaveboxes=.false. !resetealo porque puede que el slice MPI no tenga
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
             numus(jjj)=sgg%PlaneWave(jjj)%Fichero%NumSamples
             !a peticion de OLD aborto si no hay nada que iluminar
             abortar= &
             (sgg%PlaneWave(jjj)%esqx1 <=  SINPML_fullsize(iHx)%XI).and. &
             (sgg%PlaneWave(jjj)%esqx2 >=  SINPML_fullsize(iHx)%XE).and. &
             (sgg%PlaneWave(jjj)%esqy1 <=  SINPML_fullsize(iHy)%YI).and. &
             (sgg%PlaneWave(jjj)%esqy2 >=  SINPML_fullsize(iHy)%YE).and. &
             (sgg%PlaneWave(jjj)%esqz1 <=  SINPML_fullsize(iHz)%ZI).and. &
             (sgg%PlaneWave(jjj)%esqz2 >=  SINPML_fullsize(iHz)%ZE)
             if (abortar) then
                write (buff,'(a)') 'At least one of TF/SF planes must be 1 cell inside the simulation region. Aborting'
                call stoponerror(layoutnumber,size,buff)
             endif
             !!!!!!!

             IluminaTr(jjj)=.false.
             IluminaFr(jjj)=.false.
             IluminaIz(jjj)=.false.
             IluminaDe(jjj)=.false.
             IluminaAr(jjj)=.false.
             IluminaAb(jjj)=.false.
             if ((sgg%PlaneWave(jjj)%esqx1 >= sgg%SINPMLSweep(iHx)%XI).and.(sgg%PlaneWave(jjj)%esqx1 <= sgg%SINPMLSweep(IHx)%XE)) &
             IluminaTr(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqx2 <= sgg%SINPMLSweep(IHx)%XE).and.(sgg%PlaneWave(jjj)%esqx2 >= sgg%SINPMLSweep(iHx)%XI)) &
             IluminaFr(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqy1 >= sgg%SINPMLSweep(iHy)%YI).and.(sgg%PlaneWave(jjj)%esqy1 <= sgg%SINPMLSweep(IHy)%YE)) &
             IluminaIz(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqy2 <= sgg%SINPMLSweep(IHy)%YE).and.(sgg%PlaneWave(jjj)%esqy2 >= sgg%SINPMLSweep(iHy)%YI)) &
             IluminaDe(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqz1 >= sgg%SINPMLSweep(iHz)%ZI).and.(sgg%PlaneWave(jjj)%esqz1 <= sgg%SINPMLSweep(IHz)%ZE)) &
             IluminaAb(jjj)=.true.
             if ((sgg%PlaneWave(jjj)%esqz2 <= sgg%SINPMLSweep(IHz)%ZE).and.(sgg%PlaneWave(jjj)%esqz2 >= sgg%SINPMLSweep(iHz)%ZI)) &
             IluminaAr(jjj)=.true.
             !
             !find the coordinate limits of the Huygens Box for each component
             TrFr(jjj)%I%tra%Ez=Max( sgg%SINPMLSweep(iEz)%XI, sgg%PlaneWave(jjj)%esqx1    )
             TrFr(jjj)%I%fro%Ez=Min( sgg%SINPMLSweep(iEz)%XE, sgg%PlaneWave(jjj)%esqx2    )
             TrFr(jjj)%J%com%Ez=Max( sgg%SINPMLSweep(iEz)%YI, sgg%PlaneWave(jjj)%esqy1    )
             TrFr(jjj)%J%fin%Ez=Min( sgg%SINPMLSweep(iEz)%YE, sgg%PlaneWave(jjj)%esqy2    )
             TrFr(jjj)%K%com%Ez=Max( sgg%SINPMLSweep(iEz)%ZI, sgg%PlaneWave(jjj)%esqz1    )
             TrFr(jjj)%K%fin%Ez=MIn( sgg%SINPMLSweep(iEz)%ZE, sgg%PlaneWave(jjj)%esqz2-1  )
             !
             TrFr(jjj)%I%tra%Ey=Max( sgg%SINPMLSweep(iEy)%XI, sgg%PlaneWave(jjj)%esqx1    )
             TrFr(jjj)%I%fro%Ey=Min( sgg%SINPMLSweep(iEy)%XE, sgg%PlaneWave(jjj)%esqx2    )
             TrFr(jjj)%J%com%Ey=Max( sgg%SINPMLSweep(iEy)%YI, sgg%PlaneWave(jjj)%esqy1    )
             TrFr(jjj)%J%fin%Ey=Min( sgg%SINPMLSweep(iEy)%YE, sgg%PlaneWave(jjj)%esqy2-1  )
             TrFr(jjj)%K%com%Ey=Max( sgg%SINPMLSweep(iEy)%ZI ,sgg%PlaneWave(jjj)%esqz1    )
             TrFr(jjj)%K%fin%Ey=MIn( sgg%SINPMLSweep(iEy)%ZE ,sgg%PlaneWave(jjj)%esqz2    )
             !
             TrFr(jjj)%I%tra%Hy= TrFr(jjj)%I%tra%Ez - 1
             TrFr(jjj)%I%fro%Hy= TrFr(jjj)%I%fro%Ez
             TrFr(jjj)%J%com%Hy= TrFr(jjj)%J%com%Ez
             TrFr(jjj)%J%fin%Hy= TrFr(jjj)%J%fin%Ez
             TrFr(jjj)%K%com%Hy= TrFr(jjj)%K%com%Ez
             TrFr(jjj)%K%fin%Hy= TrFr(jjj)%K%fin%Ez
             !
             TrFr(jjj)%I%tra%Hz= TrFr(jjj)%I%tra%Ey - 1
             TrFr(jjj)%I%fro%Hz= TrFr(jjj)%I%fro%Ey
             TrFr(jjj)%J%com%Hz= TrFr(jjj)%J%com%Ey
             TrFr(jjj)%J%fin%Hz= TrFr(jjj)%J%fin%Ey
             TrFr(jjj)%K%com%Hz= TrFr(jjj)%K%com%Ey
             TrFr(jjj)%K%fin%Hz= TrFr(jjj)%K%fin%Ey
             !
             !
             IzDe(jjj)%J%izq%Ex=Max( sgg%SINPMLSweep(iEx)%yI, sgg%PlaneWave(jjj)%esqy1     )
             IzDe(jjj)%J%der%Ex=Min( sgg%SINPMLSweep(iEx)%yE, sgg%PlaneWave(jjj)%esqy2     )
             IzDe(jjj)%I%com%Ex=Max( sgg%SINPMLSweep(iEx)%xI, sgg%PlaneWave(jjj)%esqx1     )
             IzDe(jjj)%I%fin%Ex=Min( sgg%SINPMLSweep(iEx)%xE, sgg%PlaneWave(jjj)%esqx2-1   )
             IzDe(jjj)%K%com%Ex=Max( sgg%SINPMLSweep(iEx)%ZI ,sgg%PlaneWave(jjj)%esqz1     )
             IzDe(jjj)%K%fin%Ex=MIn( sgg%SINPMLSweep(iEx)%ZE ,sgg%PlaneWave(jjj)%esqz2     )
             !
             IzDe(jjj)%J%izq%Ez=Max( sgg%SINPMLSweep(iEz)%yI, sgg%PlaneWave(jjj)%esqy1     )
             IzDe(jjj)%J%der%Ez=Min( sgg%SINPMLSweep(iEz)%yE, sgg%PlaneWave(jjj)%esqy2     )
             IzDe(jjj)%I%com%Ez=Max( sgg%SINPMLSweep(iEz)%xI, sgg%PlaneWave(jjj)%esqx1     )
             IzDe(jjj)%I%fin%Ez=Min( sgg%SINPMLSweep(iEz)%xE, sgg%PlaneWave(jjj)%esqx2     )
             IzDe(jjj)%K%com%Ez=Max( sgg%SINPMLSweep(iEz)%ZI ,sgg%PlaneWave(jjj)%esqz1     )
             IzDe(jjj)%K%fin%Ez=MIn( sgg%SINPMLSweep(iEz)%ZE ,sgg%PlaneWave(jjj)%esqz2-1   )
             !
             IzDe(jjj)%J%izq%Hz= IzDe(jjj)%J%izq%Ex - 1
             IzDe(jjj)%J%der%Hz= IzDe(jjj)%J%der%Ex
             IzDe(jjj)%I%com%Hz= IzDe(jjj)%I%com%Ex
             IzDe(jjj)%I%fin%Hz= IzDe(jjj)%I%fin%Ex
             IzDe(jjj)%K%com%Hz= IzDe(jjj)%K%com%Ex
             IzDe(jjj)%K%fin%Hz= IzDe(jjj)%K%fin%Ex
             !
             IzDe(jjj)%J%izq%Hx= IzDe(jjj)%J%izq%Ez - 1
             IzDe(jjj)%J%der%Hx= IzDe(jjj)%J%der%Ez
             IzDe(jjj)%I%com%Hx= IzDe(jjj)%I%com%Ez
             IzDe(jjj)%I%fin%Hx= IzDe(jjj)%I%fin%Ez
             IzDe(jjj)%K%com%Hx= IzDe(jjj)%K%com%Ez
             IzDe(jjj)%K%fin%Hx= IzDe(jjj)%K%fin%Ez
             !
             !
             AbAr(jjj)%K%aba%Ey=Max( sgg%SINPMLSweep(iEy)%ZI, sgg%PlaneWave(jjj)%esqz1     )
             AbAr(jjj)%K%arr%Ey=Min( sgg%SINPMLSweep(iEy)%ZE, sgg%PlaneWave(jjj)%esqz2     )
             AbAr(jjj)%I%com%Ey=Max( sgg%SINPMLSweep(iEy)%XI, sgg%PlaneWave(jjj)%esqx1     )
             AbAr(jjj)%I%fin%Ey=Min( sgg%SINPMLSweep(iEy)%XE, sgg%PlaneWave(jjj)%esqx2     )
             AbAr(jjj)%J%com%Ey=Max( sgg%SINPMLSweep(iEy)%YI, sgg%PlaneWave(jjj)%esqy1     )
             AbAr(jjj)%J%fin%Ey=Min( sgg%SINPMLSweep(iEy)%YE, sgg%PlaneWave(jjj)%esqy2-1   )
             !
             AbAr(jjj)%K%aba%Ex=Max( sgg%SINPMLSweep(iEx)%ZI, sgg%PlaneWave(jjj)%esqz1     )
             AbAr(jjj)%K%arr%Ex=Min( sgg%SINPMLSweep(iEx)%ZE, sgg%PlaneWave(jjj)%esqz2     )
             AbAr(jjj)%I%com%Ex=Max( sgg%SINPMLSweep(iEx)%XI, sgg%PlaneWave(jjj)%esqx1     )
             AbAr(jjj)%I%fin%Ex=Min( sgg%SINPMLSweep(iEx)%XE, sgg%PlaneWave(jjj)%esqx2-1   )
             AbAr(jjj)%J%com%Ex=Max( sgg%SINPMLSweep(iEx)%YI, sgg%PlaneWave(jjj)%esqy1     )
             AbAr(jjj)%J%fin%Ex=Min( sgg%SINPMLSweep(iEx)%YE, sgg%PlaneWave(jjj)%esqy2     )
             !
             AbAr(jjj)%K%aba%Hx= AbAr(jjj)%K%aba%Ey - 1
             AbAr(jjj)%K%arr%Hx= AbAr(jjj)%K%arr%Ey
             AbAr(jjj)%I%com%Hx= AbAr(jjj)%I%com%Ey
             AbAr(jjj)%I%fin%Hx= AbAr(jjj)%I%fin%Ey
             AbAr(jjj)%J%com%Hx= AbAr(jjj)%J%com%Ey
             AbAr(jjj)%J%fin%Hx= AbAr(jjj)%J%fin%Ey
             !
             AbAr(jjj)%K%aba%Hy= AbAr(jjj)%K%aba%Ex -1
             AbAr(jjj)%K%arr%Hy= AbAr(jjj)%K%arr%Ex
             AbAr(jjj)%I%com%Hy= AbAr(jjj)%I%com%Ex
             AbAr(jjj)%I%fin%Hy= AbAr(jjj)%I%fin%Ex
             AbAr(jjj)%J%com%Hy= AbAr(jjj)%J%com%Ex
             AbAr(jjj)%J%fin%Hy= AbAr(jjj)%J%fin%Ex
             thereareplanewaveboxes=thereareplanewaveboxes.or.IluminaTr(jjj).or.IluminaFr(jjj).or.IluminaIz(jjj).or. &
             IluminaDe(jjj).or.IluminaAr(jjj).or.IluminaAb(jjj)

         enddo !barrido j planewaves
      endif  !ThereArePlaneWaveBoxes

       maxnumus=maxval(numus)
       allocate (evol(1:sgg%numplanewaves,0 : maxnumus))
       do jjj=1,sgg%numplanewaves
           do k=0,numus(jjj)
              evol(jjj,k)=sgg%PlaneWave(jjj)%fichero%Samples(k)
           end do
           deltaevol(jjj)=sgg%PlaneWave(jjj)%fichero%deltaSamples
           if (deltaevol(jjj) > sgg%dt) then
              write (buff,'(a,e12.2e3)')  'WARNING: '//trim(adjustl(sgg%PlaneWave(jjj)%Fichero%Name))// &
              ' undersampled by a factor ',deltaevol(jjj)/sgg%dt
              call WarnErrReport(buff)
           endif
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
                 modulus=sqrt(pxpw(jjj,kkk)**2+pypw(jjj,kkk)**2+pzpw(jjj,kkk)**2.0_RKIND )
                 pxpw(jjj,kkk)=pxpw(jjj,kkk)/modulus
                 pypw(jjj,kkk)=pypw(jjj,kkk)/modulus
                 pzpw(jjj,kkk)=pzpw(jjj,kkk)/modulus  
                 INCERT(jjj,kkk)=sgg%PlaneWave(jjj)%incert(kkk)
             else
                 if (sgg%PlaneWave(jjj)%isRC) then
                     read(14) pxpw(jjj,kkk),pypw(jjj,kkk),pzpw(jjj,kkk),fpw(jjj,1,kkk),fpw(jjj,2,kkk),fpw(jjj,3,kkk),INCERT(jjj,kkk)
                 else !inicializalo como siempre
                     pxpw(jjj,kkk)=sgg%PlaneWave(jjj)%px(kkk)
                     pypw(jjj,kkk)=sgg%PlaneWave(jjj)%py(kkk)
                     pzpw(jjj,kkk)=sgg%PlaneWave(jjj)%pz(kkk)
                     fpw(jjj,1,kkk)=sgg%PlaneWave(jjj)%ex(kkk)
                     fpw(jjj,2,kkk)=sgg%PlaneWave(jjj)%ey(kkk)
                     fpw(jjj,3,kkk)=sgg%PlaneWave(jjj)%ez(kkk)
    !
                     modulus=sqrt(pxpw(jjj,kkk)**2+pypw(jjj,kkk)**2+pzpw(jjj,kkk)**2.0_RKIND )
                     pxpw(jjj,kkk)=pxpw(jjj,kkk)/modulus
                     pypw(jjj,kkk)=pypw(jjj,kkk)/modulus
                     pzpw(jjj,kkk)=pzpw(jjj,kkk)/modulus  
                     INCERT(jjj,kkk)=sgg%PlaneWave(jjj)%incert(kkk)
                 endif
             endif
         end do
       end do
       do jjj=1,sgg%numplanewaves
         !Find the angles and amplitudes
         do kkk=1,sgg%PlaneWave(jjj)%nummodes
!!!! movido a donde se precisa para permit scaling 081118
!!             fpw(jjj,4,kkk)=(pypw(jjj,kkk)*fpw(jjj,3,kkk)-pzpw(jjj,kkk)*fpw(jjj,2,kkk))/zvac
!!             fpw(jjj,5,kkk)=(pzpw(jjj,kkk)*fpw(jjj,1,kkk)-pxpw(jjj,kkk)*fpw(jjj,3,kkk))/zvac
!!             fpw(jjj,6,kkk)=(pxpw(jjj,kkk)*fpw(jjj,2,kkk)-pypw(jjj,kkk)*fpw(jjj,1,kkk))/zvac
             !
             !Find the null-phase corner depending on the angle of propagation
             IF ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(iHx)%XI))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(iHy)%YI))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(iHz)%ZI))
             elseIF ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(iHx)%XI))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(iHy)%YI))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(iHz)%ZE))
             elseIF ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(iHx)%XI))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(iHy)%YE))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(iHz)%ZI))
             elseIF ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(iHx)%XE))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(iHy)%YI))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(iHz)%ZI))
             elseIF ((pxpw(jjj,kkk) >= 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(iHx)%XI))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(iHy)%YE))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(iHz)%ZE))
             elseIF ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) >= 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(iHx)%XE))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(iHy)%YE))
                ZD0=sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(iHz)%ZI))
             elseIF ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) >= 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(iHx)%XE))
                YD0=sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(iHy)%YI))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(iHz)%ZE))
             elseIF ((pxpw(jjj,kkk) < 0).and.(pypw(jjj,kkk) < 0).and.(pzpw(jjj,kkk) < 0)) then
                XD0=sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(iHx)%XE))
                YD0=sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(iHy)%YE))
                ZD0=sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(iHz)%ZE))
             else
                call stoponerror(layoutnumber,size,'buggy xo,yo,z0')
             endif
             diagonalcaja=sqrt( (sgg%Linex(max(sgg%PlaneWave(jjj)%esqx1-1,SINPML_fullsize(iHx)%XI)) - sgg%Linex(min(sgg%PlaneWave(jjj)%esqx2+1,SINPML_fullsize(iHx)%XE)))**2.0_RKIND  + &
                                (sgg%Liney(max(sgg%PlaneWave(jjj)%esqy1-1,SINPML_fullsize(iHy)%YI)) - sgg%Liney(min(sgg%PlaneWave(jjj)%esqy2+1,SINPML_fullsize(iHy)%YE)))**2.0_RKIND  + &
                                (sgg%Linez(max(sgg%PlaneWave(jjj)%esqz1-1,SINPML_fullsize(iHz)%ZI)) - sgg%Linez(min(sgg%PlaneWave(jjj)%esqz2+1,SINPML_fullsize(iHz)%ZE)))**2.0_RKIND  ) 
             distanciaInicial(jjj,kkk)=((XD0*pxpw(jjj,kkk)+YD0*pypw(jjj,kkk)+ZD0*pzpw(jjj,kkk)))-INCERT(jjj,kkk)*diagonalcaja !CREO QUE LA TENGO QUE RESTAR PARA QUE LA INCERTIDUMBRE SOLO AGRANDE LA CAJA (RETRASE LA SEniAL)
                                                                           !!!! corroboro a 150419 que hay que restar la incertidumbre, despues de dudar sobre el signo, pq luego t-d/c, d=n.r-(n.r0-incert)>0 sii incert>0, ya que n.r-n.r0>0 siempre

             
         end do !del kkk
      end do !del maxmodes

      !check if materials are crossed by the box
      do jjj=1, sgg%numplanewaves
          If( IluminaTr(jjj)) then
             !Ez Back
             i = TrFr(jjj)%I%tra%Ez !Back
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   if (sggMiEz( i, j, k) /=1 ) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Ez ',i,j,k
                      if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Ey Back
             i = TrFr(jjj)%I%tra%Ey
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   if (sggMiEy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Ey ',i,j,k
                      if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                End do
             end do
          endif
          !--->
          If( IluminaFr(jjj)) then
             !Ez  Front
             i = TrFr(jjj)%I%fro%Ez !Front
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   if (sggMiEz(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Ez ',i,j,k
                      if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Ey  Front
             i = TrFr(jjj)%I%fro%Ey !Front
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   if (sggMiEy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Ey ',i,j,k
                      if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          !--->
          If( IluminaIz(jjj)) then
             !Ex Left
             j = IzDe(jjj)%J%izq%Ex  !Left
             do k = IzDe(jjj)%K%com%Ex, IzDe(jjj)%K%fin%Ex
                do i = IzDe(jjj)%I%com%Ex, IzDe(jjj)%I%fin%Ex
                   if (sggMiEx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Ex ',i,j,k
                      if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Ez Left
             j = IzDe(jjj)%J%izq%Ez  !Left
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   if (sggMiEz(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Ez ',i,j,k
                      if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          !--->
          If( IluminaDe(jjj)) then
             !Ez  Right
             j = IzDe(jjj)%J%der%Ez !Right
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   if (sggMiEz(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Ez ',i,j,k
                      if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Ex  Right
             j = IzDe(jjj)%J%der%Ex !Right
             do k = IzDe(jjj)%K%com%Ex,IzDe(jjj)%K%fin%Ex
                do i=IzDe(jjj)%I%com%Ex,IzDe(jjj)%I%fin%Ex
                   if (sggMiEx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Ex ',i,j,k
                      if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                End do
             end do
          endif
          !--->
          If( IluminaAb(jjj)) then
             !Ex  Down
             k = AbAr(jjj)%K%aba%Ex  !Down
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                Do i=AbAr(jjj)%I%com%Ex,AbAr(jjj)%I%fin%Ex
                   if (sggMiEx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Ex ',i,j,k
                      if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Ey Down
             k = AbAr(jjj)%K%aba%Ey  !Down
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   if (sggMiEy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Ey ',i,j,k
                      if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          !--->
          If( IluminaAr(jjj)) then
             !Ex Up
             k = AbAr(jjj)%K%arr%Ex  !Up
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                do i = AbAr(jjj)%I%com%Ex, AbAr(jjj)%I%fin%Ex
                   if (sggMiEx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Ex ',i,j,k
                      if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Ey Up
             k = AbAr(jjj)%K%arr%Ey
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   if (sggMiEy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Ey ',i,j,k
                      if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          !!!
          if( IluminaTr(jjj)) then
             !Hz Back
             i = TrFr(jjj)%I%tra%Hz  !Back
             do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                   if (sggMiHz(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Hz ',i,j,k
                      if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Hy Back
             i = TrFr(jjj)%I%tra%Hy  !Back
             do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                   if (sggMiHy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Back TF/SF region intersects a material at Hy ',i,j,k
                      if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          if( IluminaFr(jjj)) then
             !Hz  Front
             i = TrFr(jjj)%I%fro%Hz !Front
             do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                   if (sggMiHz(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Hz ',i,j,k
                      if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Hy  Front
             i = TrFr(jjj)%I%fro%Hy !Front
             do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                   if (sggMiHy( i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Front TF/SF region intersects a material at Hy ',i,j,k
                      if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo

          endif
          if( IluminaIz(jjj)) then
             !Hx Left
             j = IzDe(jjj)%J%izq%Hx  !Left
             do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                   if (sggMiHx( i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Hx ',i,j,k
                      if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Hz Left
             j = IzDe(jjj)%J%izq%Hz  !Left
             do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                   if (sggMiHz( i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Left TF/SF region intersects a material at Hz ',i,j,k
                      if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          if( IluminaDe(jjj)) then
             !Hx  Right
             j = IzDe(jjj)%J%der%Hx !Right
             do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                   if (sggMiHx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Hx ',i,j,k
                      if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Hz  Right
             j = IzDe(jjj)%J%der%Hz !Right
             do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                   if (sggMiHz(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Right TF/SF region intersects a material at Hz ',i,j,k
                      if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          if( IluminaAb(jjj)) then
             !Hx  Down
             k = AbAr(jjj)%K%aba%Hx  !Down
             do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                   if (sggMiHx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Hx ',i,j,k
                      if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Hy  Down
             k = AbAr(jjj)%K%aba%Hy  !Down
             do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                do i=AbAr(jjj)%I%com%Hy,AbAr(jjj)%I%fin%Hy
                   if (sggMiHy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Down TF/SF region intersects a material at Hy ',i,j,k
                      if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
          !--->
          if( IluminaAr(jjj)) then
             !Hx Up
             k = AbAr(jjj)%K%arr%Hx  !Up
             do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                   if (sggMiHx(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Hx ',i,j,k
                      if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
             !Hy Up
             k=AbAr(jjj)%K%arr%Hy  !Up
             do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                do i = AbAr(jjj)%I%com%Hy, AbAr(jjj)%I%fin%Hy
                   if (sggMiHy(i,j,k) /=1 ) then
                      write (buff,'(a,3i7)') 'Up TF/SF region intersects a material at Hy ',i,j,k
                      if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                      ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                      (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                      call stoponerror(layoutnumber,size,buff)
                   endif
                enddo
             enddo
          endif
      end do !del j numplanewaves

!!!!
      call calc_planewaveconstants(sgg,eps0,mu0)
!!!
      RETURN
   end subroutine InitPlaneWave




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Calculate the incident field at a given time/space point
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   function Incid(sgg,jjj, nfield,time,i,j,k,still_planewave_time,calledfromobservation)    RESULT(EHI)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      logical :: still_planewave_time,calledfromobservation
      integer (KIND=4) i,j,k,nfield,jjj,kkk,jdum
      REAL (KIND=RKIND)  ::  EHI
      REAL (KIND=RKIND)  ::   time,d,xf,yf,zf
      !
      xf=Punto%PhysCoor(nfield)%x(i)
      yf=Punto%PhysCoor(nfield)%y(j)
      zf=Punto%PhysCoor(nfield)%z(k)
      ehi=0.0_RKIND

      if (calledfromobservation) then     
#ifdef CompileWithOpenMP
!$xMP   PARALLEL DO DEFAULT(SHARED) private (d,kkk,jjj) REDUCTION(+:EhI)
#endif
            do jdum=1, sgg%numplanewaves !150419 observation debe sumar las planewaves se ha movido aqui desde la llamada
              do kkk=1,sgg%PlaneWave(jdum)%nummodes
                 d=(xf*pxpw(jdum,kkk)+yf*pypw(jdum,kkk)+zf*pzpw(jdum,kkk))-distanciaInicial(jdum,kkk)
                 EhI=EhI + fpw(jdum,nfield,kkk)*evolucion(jdum,time,d,still_planewave_time)
              !!!!!!!!!!!!!!!!!!!Ehi=Ehi*exp(-0.2*((-20 + i)**2.0_RKIND + (-20 + j)**2.0_RKIND ))
              end do
            end do
#ifdef CompileWithOpenMP
!$xMP   END PARALLEL DO
#endif
      else !si no lo llama observation el jjj ya viene especificado
#ifdef CompileWithOpenMP
!$xMP   PARALLEL DO DEFAULT(SHARED) private (d,kkk,) REDUCTION(+:EhI)
#endif
              do kkk=1,sgg%PlaneWave(jjj)%nummodes
                 d=(xf*pxpw(jjj,kkk)+yf*pypw(jjj,kkk)+zf*pzpw(jjj,kkk))-distanciaInicial(jjj,kkk)
                 EhI=EhI + fpw(jjj,nfield,kkk)*evolucion(jjj,time,d,still_planewave_time)
              !!!!!!!!!!!!!!!!!!!Ehi=Ehi*exp(-0.2*((-20 + i)**2.0_RKIND + (-20 + j)**2.0_RKIND ))
              end do
#ifdef CompileWithOpenMP
!$xMP   END PARALLEL DO
#endif      
      endif
      return

   contains

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Evolution function to interpolate from the input file
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL (KIND=RKIND) function evolucion(jjj,t,d,still_planewave_time)
         REAL (KIND=RKIND) t,d
         integer (kind=8)  ::  nprev
         integer (kind=4)  ::  jjj
         logical  ::  still_planewave_time
!         if (d<=0.0_RKIND) then
!             print *,layr,' buggy error in d planewaves.evolucion. ' !ojo porque ralentiza. quitar cuando estemos seguros de RC
!         endif

         evolucion=0.0_RKIND
         nprev=int((t-d/cluz)/deltaevol(jjj))
         if ((nprev+1 <= numus(jjj))) then 
           still_planewave_time=.true. !todavia puede haber actividad
           if (nprev > 0) then
            !first order interpolation
               evolucion=(evol(jjj,nprev+1)-evol(jjj,nprev))/deltaevol(jjj)*((t-d/cluz)-nprev*deltaevol(jjj))+evol(jjj,nprev) !interpolacion lineal
            !second order !no advantages over first order
            !  if (nprev+2 > numus(jjj)) then
            !      evolucion=0.0_RKIND !se asume que el fichero de entrada contiene una excitacion que se anula despues
            !  else
            !      evolucion=evol(jjj,nprev+2) * ( ((t-d/cluz)-nprev    *deltaevol(jjj)) * ((t-d/cluz)-(nprev+1)*deltaevol(jjj)) ) /(2.0_RKIND * deltaevol(jjj)**2.0_RKIND ) - &
            !                evol(jjj,nprev+1) * ( ((t-d/cluz)-nprev    *deltaevol(jjj)) * ((t-d/cluz)-(nprev+2)*deltaevol(jjj)) ) /(   deltaevol(jjj)**2.0_RKIND ) + &
            !                evol(jjj,nprev  ) * ( ((t-d/cluz)-(nprev+2)*deltaevol(jjj)) * ((t-d/cluz)-(nprev+1)*deltaevol(jjj)) ) /(2.0_RKIND * deltaevol(jjj)**2.0_RKIND )
            !  endif
           endif
         endif
         return
      end function evolucion
      !
   end function incid

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyIlumina(sgg)
      type (SGGFDTDINFO), intent(INOUT)         ::  sgg
      integer (kind=4)  ::  field

      do field=iEx,iHz
         if (associated(Punto%PhysCoor(field)%x)) deallocate (Punto%PhysCoor(field)%x)
         if (associated(Punto%PhysCoor(field)%y)) deallocate (Punto%PhysCoor(field)%y)
         if (associated(Punto%PhysCoor(field)%z)) deallocate (Punto%PhysCoor(field)%z)
      end do

      if (sgg%numplanewaves >=1) then
       deallocate (TrFr, IzDe,AbAr, IluminaTr, IluminaFr, IluminaIz,IluminaDe, IluminaAr,IluminaAb, pxpw, pypw, pzpw,   fpw, INCERT, numus,deltaevol,distanciainicial)
      endif
      if (allocated(evol)) deallocate (evol)
      if (associated(sgg%PlaneWave)) deallocate (sgg%PlaneWave)
   end subroutine DestroyIlumina



   !**************************************************************************************************
   subroutine AdvancePlaneWaveE( sgg, timeinstant, b, g2, Idxh, Idyh, Idzh, Ex, Ey, Ez,still_planewave_time)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      logical :: still_planewave_time
      logical :: called_fromobservation
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent( IN)  ::   timeinstant
      !!!
      type( bounds_t), intent( IN)  ::  b
      !--->
      real (kind = RKIND), dimension( 0 :  sgg%NumMedia), intent( IN)  ::  g2
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN)  ::  Idxh
      real (kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN)  ::  Idyh
      real (kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN)  ::  Idzh
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT)  ::  Ez
      !---------------------------> variables locales <-----------------------------------------------
      real (kind = RKIND)  ::  timei, G2_1, Id,incidente
      integer  ::  i, j, k, i_m, j_m, k_m,jjj
      character (LEN=BUFSIZE)     ::  dubuf
      !---------------------------> empieza AdvancePlaneWaveE <---------------------------------------
!!!!

      !!!
      still_planewave_time=.false. !por defecto no va a haber mas actividad de onda plana, a menos que pase por algun incid no trivial
      called_fromobservation=.false. !210419 
      
      timei = sgg%tiempo(timeinstant)
      !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219 
      !!! timei = (timeinstant +3) * sgg%dt !ORIGINAL sync
      
      G2_1 = G2( 1)
      !--->
      do jjj=1, sgg%numplanewaves
          If( IluminaTr(jjj)) then
             !Ez Back
             i = TrFr(jjj)%I%tra%Ez !Back
             i_m = i - b%Ez%XI
             Id = Idxh( i_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   j_m = j - b%Ez%YI
                   !--->
                   incidente = Incid(sgg,jjj, iHy, timei, i-1, j, k,still_planewave_time,called_fromobservation)
                   Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) - G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey Back
             i = TrFr(jjj)%I%tra%Ey  !Back
             i_m = i - b%Ey%XI
             Id = Idxh( i_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                k_m = k - b%Ey%ZI
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   j_m = j - b%Ey%YI
                   !--->
                   incidente = Incid(sgg,jjj, iHz, timei, i-1, j, k,still_planewave_time,called_fromobservation)
                   Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) + G2_1 * incidente * Id
                End do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          endif
          !--->
          If( IluminaFr(jjj)) then
             !Ez  Front
             i = TrFr(jjj)%I%fro%Ez !Front
             i_m = i - b%Ez%XI
             Id = Idxh( i_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ez, TrFr(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do j = TrFr(jjj)%J%com%Ez, TrFr(jjj)%J%fin%Ez
                   j_m = j - b%Ez%YI
                   !--->
                   incidente = Incid(sgg,jjj, iHy, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) + G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey  Front
             i = TrFr(jjj)%I%fro%Ey !Front
             i_m = i - b%Ey%XI
             Id = Idxh( i_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
             do k = TrFr(jjj)%K%com%Ey, TrFr(jjj)%K%fin%Ey
                k_m = k - b%Ey%ZI
                do j = TrFr(jjj)%J%com%Ey, TrFr(jjj)%J%fin%Ey
                   j_m = j - b%Ey%YI
                   !--->
                   incidente = Incid(sgg,jjj, iHz, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) - G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          endif
          !--->
          If( IluminaIz(jjj)) then
             !Ex Left
             j = IzDe(jjj)%J%izq%Ex  !Left
             j_m = j - b%Ex%YI
             Id = Idyh( j_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ex, IzDe(jjj)%K%fin%Ex
                k_m = k - b%Ex%ZI
                do i = IzDe(jjj)%I%com%Ex, IzDe(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHz, timei, i, j-1, k,still_planewave_time,called_fromobservation)
                   Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) - G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ez Left
             j = IzDe(jjj)%J%izq%Ez  !Left
             j_m = j - b%Ez%YI
             Id = Idyh( j_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   i_m = i - b%Ez%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHx, timei, i, j-1, k,still_planewave_time,called_fromobservation)
                   Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) + G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          endif
          !--->
          If( IluminaDe(jjj)) then
             !Ez  Right
             j = IzDe(jjj)%J%der%Ez !Right
             j_m = j - b%Ez%YI
             Id = Idyh( j_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ez, IzDe(jjj)%K%fin%Ez
                k_m = k - b%Ez%ZI
                do i = IzDe(jjj)%I%com%Ez, IzDe(jjj)%I%fin%Ez
                   i_m = i - b%Ez%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHx, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) - G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ex  Right
             j = IzDe(jjj)%J%der%Ex !Right
             j_m = j - b%Ex%YI
             Id = Idyh( j_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
             do k = IzDe(jjj)%K%com%Ex,IzDe(jjj)%K%fin%Ex
                k_m = k - b%Ex%ZI
                do i=IzDe(jjj)%I%com%Ex,IzDe(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHz, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) + G2_1 * incidente * Id
                End do
             end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          endif
          !--->
          If( IluminaAb(jjj)) then
             !Ex  Down
             k = AbAr(jjj)%K%aba%Ex  !Down
             k_m = k - b%Ex%ZI
             Id = Idzh( k_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                j_m = j - b%Ex%YI
                Do i=AbAr(jjj)%I%com%Ex,AbAr(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHy, timei, i, j, k-1,still_planewave_time,called_fromobservation)
                   Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) + G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey Down
             k = AbAr(jjj)%K%aba%Ey  !Down
             k_m = k - b%Ey%ZI
             Id = Idzh( k_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                j_m = j - b%Ey%YI
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   i_m = i - b%Ey%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHx, timei, i, j, k-1,still_planewave_time,called_fromobservation)
                   Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) - G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          endif
          !--->
          If( IluminaAr(jjj)) then
             !Ex Up
             k = AbAr(jjj)%K%arr%Ex  !Up
             k_m = k - b%Ex%ZI
             Id = Idzh( k_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ex, AbAr(jjj)%J%fin%Ex
                j_m = j - b%Ex%YI
                do i = AbAr(jjj)%I%com%Ex, AbAr(jjj)%I%fin%Ex
                   i_m = i - b%Ex%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHy, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) - G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
             !Ey Up
             k = AbAr(jjj)%K%arr%Ey  !Up
             k_m = k - b%Ey%ZI
             Id = Idzh( k_m )
             !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
             do j = AbAr(jjj)%J%com%Ey, AbAr(jjj)%J%fin%Ey
                j_m = j - b%Ey%YI
                do i = AbAr(jjj)%I%com%Ey, AbAr(jjj)%I%fin%Ey
                   i_m = i - b%Ey%XI
                   !--->
                   incidente = Incid(sgg,jjj, iHx, timei, i, j, k,still_planewave_time,called_fromobservation)
                   Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) + G2_1 * incidente * Id
                enddo
             enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
          endif
      end do
      !---------------------------> acaba AdvancePlaneWaveE <-----------------------------------------
      return
   endsubroutine AdvancePlaneWaveE
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the H-field at n+0.5_RKIND
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**************************************************************************************************
   subroutine AdvancePlaneWaveH(sgg, timeinstant,  b, gm2, Idxe, Idye, Idze, Hx, Hy, Hz,still_planewave_time)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      logical :: still_planewave_time
      logical :: called_fromobservation
      
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent( IN)  ::   timeinstant
      !!!
      type( bounds_t), intent( IN)  ::  b
      !--->
      real (kind = RKIND), dimension( 0 :  sgg%NumMedia), intent( IN)  ::  gm2
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxe%NX-1), intent( IN)  ::  Idxe
      real (kind = RKIND), dimension( 0 :  b%dye%NY-1), intent( IN)  ::  Idye
      real (kind = RKIND), dimension( 0 :  b%dze%NZ-1), intent( IN)  ::  Idze
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      real (kind = RKIND)  ::  timei, Gm2_1, Id,incidente
      integer (kind=4)  ::  i, j, k, i_m, j_m, k_m,jjj
      character (LEN=BUFSIZE)     ::  dubuf
      !---------------------------> empieza AdvancePlaneWaveH <---------------------------------------
      still_planewave_time=.false. !por defecto no va a haber mas actividad de onda plana, a menos que pase por algun incid no trivial
      called_fromobservation=.false. !210419 
      !!!
      !!!
      
      timei = sgg%tiempo(timeinstant) + 0.5_RKIND  * sgg%dt
      !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219 
      !!! timei = ( timeinstant + 0.5_RKIND  +3.0_RKIND) * sgg%dt  !ORIGINAL sync
      Gm2_1 = Gm2(1)
      !--->
     do jjj=1, sgg%numplanewaves
              if( IluminaTr(jjj)) then
                 !Hz Back
                 i = TrFr(jjj)%I%tra%Hz  !Back
                 i_m = i - b%Hz%XI
                 Id = Idxe( i_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       !--->
                       incidente = Incid(sgg,jjj, iEy, timei, i+1, j, k,still_planewave_time,called_fromobservation)
                       Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Back
                 i = TrFr(jjj)%I%tra%Hy  !Back
                 i_m = i - b%Hy%XI
                 Id = Idxe( i_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       !--->
                       incidente = Incid(sgg,jjj,  iEz, timei, i+1, j, k,still_planewave_time,called_fromobservation)
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaFr(jjj)) then
                 !Hz  Front
                 i = TrFr(jjj)%I%fro%Hz !Front
                 i_m = i - b%Hz%XI
                 Id = Idxe( i_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       !--->
                       incidente = Incid(sgg,jjj,  iEy, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Front
                 i = TrFr(jjj)%I%fro%Hy !Front
                 i_m = i - b%Hy%XI
                 Id = Idxe( i_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       !--->
                       incidente = Incid(sgg,jjj,  iEz, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaIz(jjj)) then
                 !Hx Left
                 j = IzDe(jjj)%J%izq%Hx  !Left
                 j_m = j - b%Hx%YI
                 Id = Idye( j_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEz, timei, i, j+1, k,still_planewave_time,called_fromobservation)
                       Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz Left
                 j = IzDe(jjj)%J%izq%Hz  !Left
                 j_m = j - b%Hz%YI
                 Id = Idye( j_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEx, timei, i, j+1, k,still_planewave_time,called_fromobservation)
                       Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaDe(jjj)) then
                 !Hx  Right
                 j = IzDe(jjj)%J%der%Hx !Right
                 j_m = j - b%Hx%YI
                 Id = Idye( j_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEz, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz  Right
                 j = IzDe(jjj)%J%der%Hz !Right
                 j_m = j - b%Hz%YI
                 Id = Idye( j_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEx, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hz( i_m, j_m, k_m)=Hz( i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaAb(jjj)) then
                 !Hx  Down
                 k = AbAr(jjj)%K%aba%Hx  !Down
                 k_m = k - b%Hx%ZI
                 Id = Idze( k_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEy, timei, i, j, k+1,still_planewave_time,called_fromobservation)
                       Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Down
                 k = AbAr(jjj)%K%aba%Hy  !Down
                 k_m = k - b%Hy%ZI
                 Id = Idze( k_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i=AbAr(jjj)%I%com%Hy,AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEx, timei, i, j, k+1,still_planewave_time,called_fromobservation)
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaAr(jjj)) then
                 !Hx Up
                 k = AbAr(jjj)%K%arr%Hx  !Up
                 k_m = k - b%Hx%ZI
                 Id = Idze( k_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEy, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Up
                 k=AbAr(jjj)%K%arr%Hy  !Up
                 k_m = k - b%Hy%ZI
                 Id = Idze( k_m )
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (incidente,i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i = AbAr(jjj)%I%com%Hy, AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       !--->
                       incidente = Incid(sgg,jjj,  iEx, timei, i, j, k,still_planewave_time,called_fromobservation)
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - Gm2_1 * incidente * Id
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
      end do 
      !---------------------------> acaba AdvancePlaneWaveH <-----------------------------------------
      return
   endsubroutine AdvancePlaneWaveH

    subroutine storeplanewaves(sgg)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
       integer (kind=4) :: jjj,kkk
       do jjj=1,sgg%numplanewaves
         do kkk=1,sgg%PlaneWave(jjj)%nummodes
            if (sgg%PlaneWave(jjj)%isRC) then
                 write(14,err=634) pxpw(jjj,kkk),pypw(jjj,kkk),pzpw(jjj,kkk),fpw(jjj,1,kkk),fpw(jjj,2,kkk),fpw(jjj,3,kkk),sgg%PlaneWave(jjj)%incert(kkk)
            endif
         end do
       end do
      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'PLANEWAVES: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
    end subroutine storeplanewaves

    subroutine calc_planewaveconstants(sgg,eps00,mu00)
      type (SGGFDTDINFO), intent(IN)   ::  sgg
      real (kind = RKIND), intent(in)  ::  eps00,mu00
      integer :: jjj,kkk
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      cluz=1.0_RKIND/sqrt(eps0*mu0) !lo necesitara incid
      zvac=sqrt(mu0/eps0) !lo necesitan las variables de mas abajo
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
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      type( bounds_t), intent( IN)  ::  b
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hxvac
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hyvac
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hzvac
      !---------------------------> variables locales <-----------------------------------------------
      integer (kind=4)  ::  i, j, k, i_m, j_m, k_m,jjj

      do jjj=1, sgg%numplanewaves
              if( IluminaTr(jjj)) then
                 !Hz Back
                 i = TrFr(jjj)%I%tra%Hz  !Back
                 i_m = i - b%Hz%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       !--->
                       Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - Hzvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Back
                 i = TrFr(jjj)%I%tra%Hy  !Back
                 i_m = i - b%Hy%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - Hyvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaFr(jjj)) then
                 !Hz  Front
                 i = TrFr(jjj)%I%fro%Hz !Front
                 i_m = i - b%Hz%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hz, TrFr(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do j = TrFr(jjj)%J%com%Hz, TrFr(jjj)%J%fin%Hz
                       j_m = j - b%Hz%YI
                       Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - Hzvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Front
                 i = TrFr(jjj)%I%fro%Hy !Front
                 i_m = i - b%Hy%XI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m)
#endif
                 do k = TrFr(jjj)%K%com%Hy, TrFr(jjj)%K%fin%Hy
                    k_m = k - b%Hy%ZI
                    do j = TrFr(jjj)%J%com%Hy, TrFr(jjj)%J%fin%Hy
                       j_m = j - b%Hy%YI
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - Hyvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaIz(jjj)) then
                 !Hx Left
                 j = IzDe(jjj)%J%izq%Hx  !Left
                 j_m = j - b%Hx%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) - Hxvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz Left
                 j = IzDe(jjj)%J%izq%Hz  !Left
                 j_m = j - b%Hz%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - Hzvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaDe(jjj)) then
                 !Hx  Right
                 j = IzDe(jjj)%J%der%Hx !Right
                 j_m = j - b%Hx%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hx, IzDe(jjj)%K%fin%Hx
                    k_m = k - b%Hx%ZI
                    do i = IzDe(jjj)%I%com%Hx, IzDe(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) - Hxvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hz  Right
                 j = IzDe(jjj)%J%der%Hz !Right
                 j_m = j - b%Hz%YI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (k,i,k_m,i_m)
#endif
                 do k = IzDe(jjj)%K%com%Hz, IzDe(jjj)%K%fin%Hz
                    k_m = k - b%Hz%ZI
                    do i = IzDe(jjj)%I%com%Hz, IzDe(jjj)%I%fin%Hz
                       i_m = i - b%Hz%XI
                       Hz( i_m, j_m, k_m)=Hz( i_m, j_m, k_m) - Hzvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaAb(jjj)) then
                 !Hx  Down
                 k = AbAr(jjj)%K%aba%Hx  !Down
                 k_m = k - b%Hx%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m) - Hxvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy  Down
                 k = AbAr(jjj)%K%aba%Hy  !Down
                 k_m = k - b%Hy%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i=AbAr(jjj)%I%com%Hy,AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - Hyvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
              !--->
              if( IluminaAr(jjj)) then
                 !Hx Up
                 k = AbAr(jjj)%K%arr%Hx  !Up
                 k_m = k - b%Hx%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hx, AbAr(jjj)%J%fin%Hx
                    j_m = j - b%Hx%YI
                    do i = AbAr(jjj)%I%com%Hx, AbAr(jjj)%I%fin%Hx
                       i_m = i - b%Hx%XI
                       Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) - Hxvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
                 !Hy Up
                 k=AbAr(jjj)%K%arr%Hy  !Up
                 k_m = k - b%Hy%ZI
                 !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m)
#endif
                 do j = AbAr(jjj)%J%com%Hy, AbAr(jjj)%J%fin%Hy
                    j_m = j - b%Hy%YI
                    do i = AbAr(jjj)%I%com%Hy, AbAr(jjj)%I%fin%Hy
                       i_m = i - b%Hy%XI
                       Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - Hyvac( i_m, j_m, k_m)
                    enddo
                 enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
              endif
      end do     
      
      
      
    return
    end subroutine corrigeondaplanaH
    

END MODULE ILUMINA
