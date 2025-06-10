
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Borders :  PML  handling
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BORDERS_CPML
   use fdetypes
   use report
   implicit none
   private
   !
   !
   REAL (KIND=RKIND), parameter  ::  StaticFrequency=1.0e14_RKIND
   ! Limits of the PML region
   type XYZlimit_tvar
      integer (kind=4), dimension(1:6)  ::  XI,XE,YI,YE,ZI,ZE
   end type XYZlimit_tvar
   type (XYZlimit_tvar), dimension(1:6)  ::    PMLc

   type LR
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Psi_Exy,Psi_Ezy,Psi_Hxy,Psi_Hzy
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Psi_Exyvac,Psi_Ezyvac,Psi_Hxyvac,Psi_Hzyvac
   end type
   type DU
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Psi_Eyz,Psi_Exz,Psi_Hyz,Psi_Hxz
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Psi_Eyzvac,Psi_Exzvac,Psi_Hyzvac,Psi_Hxzvac
   end type
   type BF
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Psi_Ezx,Psi_Eyx,Psi_Hzx,Psi_Hyx
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Psi_Ezxvac,Psi_Eyxvac,Psi_Hzxvac,Psi_Hyxvac
   end type



   !LOCAL VARIABLES
   type (LR), dimension(left : right) , save ::  regLR
   type (DU), dimension(down : up)    , save ::  regDU
   type (BF), dimension(back : front) , save ::  regBF
   REAL (KIND=RKIND) , pointer, dimension ( : , : ) , SAVE  ::   sig_max
   REAL (KIND=RKIND) , pointer, dimension ( : , : ) , SAVE  ::   aPar_max ,kPar_max
   REAL (KIND=RKIND) , pointer, dimension ( : ) , SAVE ::   P_ce_x ,P_ce_y ,P_ce_z ,P_be_x ,P_be_y ,P_be_z,&
   P_cm_x ,P_cm_y,P_cm_z ,P_bm_x ,P_bm_y ,P_bm_z
   REAL (KIND=RKIND) , pointer, dimension ( : ), SAVE  ::  ce_x ,ce_y ,ce_z ,cm_x ,cm_y ,cm_z , &
   Ice_x ,Ice_y ,Ice_z ,Icm_x ,Icm_y ,Icm_z

!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  zvac
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   REAL (KIND=RKIND), save     ::  alphamaxpar,alphaOrden,kappamaxpar
!cpml stretching maximum parameters !!alphamaxpar=StaticFrequency*2*pi*Eps0
   type (limit_t), dimension(1:6), save  ::  SINPML_fullsize
   REAL (KIND=RKIND) , dimension (:)   ,  allocatable, save    :: dxe, dye,dze,dxh,dyh,dzh
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   public  ::  InitCPMLBorders, AdvanceelectricCPML,AdvanceMagneticCPML,StoreFieldsCPMLBorders,DestroyCPMLBorders,AdvanceelectricCPML_freespace,AdvanceMagneticCPML_freespace
   public  ::  calc_cpmlconstants
   !!!public  ::  FreeSpace_AdvanceMagneticCPML,calc_cpmlconstants
contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes PML data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitCPMLBorders(sgg,layoutnumber,temp_SINPML_Fullsize,ThereArePMLBorders,resume, &
   temp_dxe,temp_dye,temp_dze,temp_dxh,temp_dyh,temp_dzh,Idxe,Idye,Idze,Idxh,Idyh,Idzh,temp_alphamaxpar,temp_alphaOrden,temp_kappamaxpar,eps00,mu00)
      REAL (KIND=RKIND)           ::  eps00,mu00
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      REAL (KIND=RKIND) , dimension (:)   ,  intent(in)    ::  &
      temp_dxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      temp_dye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      temp_dze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      temp_dxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      temp_dyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      temp_dzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)
      REAL (KIND=RKIND) , dimension (:)   , intent(inout)      ::  &
      Idxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      Idye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      Idze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)
       REAL (KIND=RKIND),  intent(in)   ::  temp_alphamaxpar,temp_alphaOrden,temp_kappamaxpar
       type (limit_t), dimension(1:6),  intent(in)  ::  temp_SINPML_fullsize
      !!!
      !
      logical  ::  ThereArePMLBorders,resume
      integer (kind=4)  ::  i,j,k,region,field,layoutnumber
!      character(len=BUFSIZE) :: buff
!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      SINPML_fullsize =  temp_SINPML_fullsize
      alphamaxpar = temp_alphamaxpar
      alphaOrden  = temp_alphaOrden
      kappamaxpar = temp_kappamaxpar
      allocate (dxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE))
      allocate (dye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE))
      allocate (dze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE))
      allocate (dxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE))
      allocate (dyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE))
      allocate (dzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE))
      dxe = temp_dxe
      dye = temp_dye
      dze = temp_dze
      dxh = temp_dxh
      dyh = temp_dyh
      dzh = temp_dzh
!
      !!!
      !character (len=BUFSIZE) :: donde
      !integer (KIND=4) :: layoutnumber

      !
      ThereArePMLBorders=.false.
      if (sgg%Border%IsBackPML.or.sgg%Border%IsFrontPML.or.sgg%Border%IsLeftPML.or.sgg%Border%IsRightPML.or. &
      sgg%Border%IsUpPML.or.sgg%Border%IsDownPML) ThereArePMLBorders=.true.
      IF (.not.(ThereArePMLBorders)) return


      !Find the limits of each of the 6 padding PML regions for each field component


      do field=iEx,iHz
         !
         PMLc(field)%XI(Down)  =                                sgg%Sweep(field)%XI
         PMLc(field)%XE(Down)  =                                sgg%Sweep(field)%XE
         PMLc(field)%YI(Down)  =                                sgg%Sweep(field)%YI
         PMLc(field)%YE(Down)  =                                sgg%Sweep(field)%YE
         PMLc(field)%ZI(Down)  =                                sgg%Sweep(field)%ZI
         PMLc(field)%ZE(Down)  =min(SINPML_Fullsize(field)%ZI-1,sgg%Sweep(field)%ZE)
         !
         PMLc(field)%XI(Up)    =                                sgg%Sweep(field)%XI
         PMLc(field)%XE(Up)    =                                sgg%Sweep(field)%XE
         PMLc(field)%YI(Up)    =                                sgg%Sweep(field)%YI
         PMLc(field)%YE(Up)    =                                sgg%Sweep(field)%YE
         PMLc(field)%ZI(Up)    =max(SINPML_Fullsize(field)%ZE+1,sgg%Sweep(field)%ZI)
         PMLc(field)%ZE(Up)    =                                sgg%Sweep(field)%ZE
         !
         PMLc(field)%XI(Left)  =                                sgg%Sweep(field)%XI
         PMLc(field)%XE(Left)  =                                sgg%Sweep(field)%XE
         PMLc(field)%YI(Left)  =                                sgg%Sweep(field)%YI
         PMLc(field)%YE(Left)  =min(SINPML_Fullsize(field)%YI-1,sgg%Sweep(field)%YE)
         PMLc(field)%ZI(Left)  =                                sgg%Sweep(field)%ZI
         PMLc(field)%ZE(Left)  =                                sgg%Sweep(field)%ZE
         !
         PMLc(field)%XI(Right) =                                sgg%Sweep(field)%XI
         PMLc(field)%XE(Right) =                                sgg%Sweep(field)%XE
         PMLc(field)%YI(Right) =max(SINPML_Fullsize(field)%YE+1,sgg%Sweep(field)%YI)
         PMLc(field)%YE(Right) =                                sgg%Sweep(field)%YE
         PMLc(field)%ZI(Right) =                                sgg%Sweep(field)%ZI
         PMLc(field)%ZE(Right) =                                sgg%Sweep(field)%ZE
         !
         PMLc(field)%XI(Back)  =                                sgg%Sweep(field)%XI
         PMLc(field)%XE(Back)  =min(SINPML_Fullsize(field)%XI-1,sgg%Sweep(field)%XE)
         PMLc(field)%YI(Back)  =                                sgg%Sweep(field)%YI
         PMLc(field)%YE(Back)  =                                sgg%Sweep(field)%YE
         PMLc(field)%ZI(Back)  =                                sgg%Sweep(field)%ZI
         PMLc(field)%ZE(Back)  =                                sgg%Sweep(field)%ZE
         !
         PMLc(field)%XI(Front) =max(SINPML_Fullsize(field)%XE+1,sgg%Sweep(field)%XI)
         PMLc(field)%XE(Front) =                                sgg%Sweep(field)%XE
         PMLc(field)%YI(Front) =                                sgg%Sweep(field)%YI
         PMLc(field)%YE(Front) =                                sgg%Sweep(field)%YE
         PMLc(field)%ZI(Front) =                                sgg%Sweep(field)%ZI
         PMLc(field)%ZE(Front) =                                sgg%Sweep(field)%ZE
         !
      end do


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! PML stuff
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ALLOCATE (sig_max(1 : 3,1 : 2),aPar_max(1 : 3,1 : 2),kPar_max(1 : 3,1 : 2))
      ALLOCATE (P_ce_x(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE ), &
      P_ce_y(sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE ), &
      P_ce_z(sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE ), &
      P_be_x(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE ), &
      P_be_y(sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE ), &
      P_be_z(sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE ), &
      P_cm_x(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE ), &
      P_cm_y(sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE ), &
      P_cm_z(sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE ), &
      P_bm_x(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE ), &
      P_bm_y(sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE ), &
      P_bm_z(sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE ))
      allocate ( ce_x(sgg%alloc(iHx)%XI  : sgg%alloc(iHx)%XE ), &
      ce_y(sgg%alloc(iHy)%YI  : sgg%alloc(iHy)%YE ), &
      ce_z(sgg%alloc(iHz)%ZI  : sgg%alloc(iHz)%ZE ),  &
      cm_x(sgg%alloc(iHx)%XI  : sgg%alloc(iHx)%XE ), &
      cm_y(sgg%alloc(iHy)%YI  : sgg%alloc(iHy)%YE ), &
      cm_z(sgg%alloc(iHz)%ZI  : sgg%alloc(iHz)%ZE ))
      allocate (Ice_x(sgg%alloc(iHx)%XI  : sgg%alloc(iHx)%XE ),&
      Ice_y(sgg%alloc(iHy)%YI  : sgg%alloc(iHy)%YE ),&
      Ice_z(sgg%alloc(iHz)%ZI  : sgg%alloc(iHz)%ZE ), &
      Icm_x(sgg%alloc(iHx)%XI  : sgg%alloc(iHx)%XE ),&
      Icm_y(sgg%alloc(iHy)%YI  : sgg%alloc(iHy)%YE ),&
      Icm_z(sgg%alloc(iHz)%ZI  : sgg%alloc(iHz)%ZE ))
      ce_x=0;  ce_y=0;  ce_z=0;  cm_x=0;  cm_y=0;  cm_z=0;
      Ice_x=0; Ice_y=0; Ice_z=0; Icm_x=0; Icm_y=0; Icm_z=0;

      !depth information matrices
      do i=sgg%ALLOC(iHx)%XI,sgg%ALLOC(iHx)%XE
         if (i <= SINPML_Fullsize(iHx)%XI .and. sgg%PML%NumLayers(1,1) /= 0) then
            ce_x (i)=                        1.0_RKIND * (SINPML_Fullsize(iHx)%XI-i) /(sgg%PML%NumLayers(1,1))
            Ice_x(i)=1.0_RKIND * (sgg%PML%NumLayers(1,1)-(SINPML_Fullsize(iHx)%XI-i))/(sgg%PML%NumLayers(1,1))
         elseif (i >= SINPML_Fullsize(iHx)%XE .and. sgg%PML%NumLayers(1,2) /= 0) then
            ce_x(i)=                         1.0_RKIND * (i-SINPML_Fullsize(iHx)%XE) /(sgg%PML%NumLayers(1,2))
            Ice_x(i)=1.0_RKIND * (sgg%PML%NumLayers(1,2)-(i-SINPML_Fullsize(iHx)%XE))/(sgg%PML%NumLayers(1,2))
         else
            ce_x(i) =0.0_RKIND
            Ice_x(i)=0.0_RKIND
         endif
      end do
      do i=sgg%ALLOC(iHx)%XI,sgg%ALLOC(iHx)%XE
         if (i <= SINPML_Fullsize(iHx)%XI-1 .and. sgg%PML%NumLayers(1,1) /= 0) then
            cm_x(i)=                         1.0_RKIND * (SINPML_Fullsize(iHx)%XI-(i+0.5_RKIND )) /(sgg%PML%NumLayers(1,1))
            Icm_x(i)=1.0_RKIND * (sgg%PML%NumLayers(1,1)-(SINPML_Fullsize(iHx)%XI-(i+0.5_RKIND )))/(sgg%PML%NumLayers(1,1))
         elseif (i >= SINPML_Fullsize(iHx)%XE .and. sgg%PML%NumLayers(1,2) /= 0) then  !!cuidado pues los h empiezan antes
            cm_x(i)=                         1.0_RKIND * (i-SINPML_Fullsize(iHx)%XE+0.5_RKIND ) /(sgg%PML%NumLayers(1,2))
            Icm_x(i)=1.0_RKIND * (sgg%PML%NumLayers(1,2)-(i-SINPML_Fullsize(iHx)%XE+0.5_RKIND ))/(sgg%PML%NumLayers(1,2))
         else
            cm_x(i)=0.0_RKIND
            Icm_x(i)=0.0_RKIND
         endif
      end do
      do j=sgg%ALLOC(iHy)%YI,sgg%ALLOC(iHy)%YE
         if (j <= SINPML_Fullsize(iHy)%YI .and. sgg%PML%NumLayers(2,1) /= 0) then
            ce_y(j)=                         1.0_RKIND * (SINPML_Fullsize(iHy)%YI-j) /(sgg%PML%NumLayers(2,1))
            Ice_y(j)=1.0_RKIND * (sgg%PML%NumLayers(2,1)-(SINPML_Fullsize(iHy)%YI-j))/(sgg%PML%NumLayers(2,1))
         elseif (j >= SINPML_Fullsize(iHy)%YE .and. sgg%PML%NumLayers(2,2) /= 0) then
            ce_y(j)=                         1.0_RKIND * (j-SINPML_Fullsize(iHy)%YE) /(sgg%PML%NumLayers(2,2))
            Ice_y(j)=1.0_RKIND * (sgg%PML%NumLayers(2,2)-(j-SINPML_Fullsize(iHy)%YE))/(sgg%PML%NumLayers(2,2))
         else
            ce_y(j)=0.0_RKIND
            Ice_y(j)=0.0_RKIND
         endif
      end do
      do j=sgg%ALLOC(iHy)%YI,sgg%ALLOC(iHy)%YE
         if (j <= SINPML_Fullsize(iHy)%YI-1 .and. sgg%PML%NumLayers(2,1) /= 0) then
            cm_y(j)=                         1.0_RKIND * (SINPML_Fullsize(iHy)%YI-(j+0.5_RKIND )) /(sgg%PML%NumLayers(2,1))
            Icm_y(j)=1.0_RKIND * (sgg%PML%NumLayers(2,1)-(SINPML_Fullsize(iHy)%YI-(j+0.5_RKIND )))/(sgg%PML%NumLayers(2,1))
         elseif (j >= SINPML_Fullsize(iHy)%YE .and. sgg%PML%NumLayers(2,2) /= 0) then
            cm_y(j)=                         1.0_RKIND * (j-SINPML_Fullsize(iHy)%YE+0.5_RKIND ) /(sgg%PML%NumLayers(2,2))
            Icm_y(j)=1.0_RKIND * (sgg%PML%NumLayers(2,2)-(j-SINPML_Fullsize(iHy)%YE+0.5_RKIND ))/(sgg%PML%NumLayers(2,2))
         else
            cm_y(j)=0.0_RKIND
            Icm_y(j)=0.0_RKIND
         endif
      end do
      do k=sgg%ALLOC(iHz)%ZI,sgg%ALLOC(iHz)%ZE
         if (k <= SINPML_Fullsize(iHz)%ZI .and. sgg%PML%NumLayers(3,1) /= 0) then
            ce_z(k)=                         1.0_RKIND * (SINPML_Fullsize(iHz)%ZI-k)/(sgg%PML%NumLayers(3,1))
            Ice_z(k)=1.0_RKIND * (sgg%PML%NumLayers(3,1)-(SINPML_Fullsize(iHz)%ZI-k))/(sgg%PML%NumLayers(3,1))
         elseif (k >= SINPML_Fullsize(iHz)%ZE .and. sgg%PML%NumLayers(3,2) /= 0) then
            ce_z(k)=                         1.0_RKIND * (k-SINPML_Fullsize(iHz)%ZE) /(sgg%PML%NumLayers(3,2))
            Ice_z(k)=1.0_RKIND * (sgg%PML%NumLayers(3,2)-(k-SINPML_Fullsize(iHz)%ZE))/(sgg%PML%NumLayers(3,2))
         else
            ce_z(k)=0.0_RKIND
            Ice_z(k)=0.0_RKIND
         endif
      end do
      do k=sgg%ALLOC(iHz)%ZI,sgg%ALLOC(iHz)%ZE
         if (k <= SINPML_Fullsize(iHz)%ZI-1 .and. sgg%PML%NumLayers(3,1) /= 0) then
            cm_z(k)=                         1.0_RKIND * (SINPML_Fullsize(iHz)%ZI-(k+0.5_RKIND )) /(sgg%PML%NumLayers(3,1))
            Icm_z(k)=1.0_RKIND * (sgg%PML%NumLayers(3,1)-(SINPML_Fullsize(iHz)%ZI-(k+0.5_RKIND )))/(sgg%PML%NumLayers(3,1))
         elseif (k >= SINPML_Fullsize(iHz)%ZE .and. sgg%PML%NumLayers(3,2) /= 0) then
            cm_z(k)=                         1.0_RKIND * (k-SINPML_Fullsize(iHz)%ZE+0.5_RKIND ) /(sgg%PML%NumLayers(3,2))
            Icm_z(k)=1.0_RKIND * (sgg%PML%NumLayers(3,2)-(k-SINPML_Fullsize(iHz)%ZE+0.5_RKIND ))/(sgg%PML%NumLayers(3,2))
         else
            cm_z(k)=0.0_RKIND
            Icm_z(k)=0.0_RKIND
         endif
      end do



      call calc_cpmlconstants(sgg,Idxe,Idye,Idze,Idxh,Idyh,Idzh,eps0,mu0)

      !Fake coms and ends IN CASE OF NO pml SO THAT NEVER ENTER THE DO FOR THESE CASES
      IF (.not.(sgg%Border%IsDownPML)) PMLc(1:6)%ZI(down)=PMLc(1:6)%ZE(down)+100
      IF (.not.(sgg%Border%IsUpPML))   PMLc(1:6)%ZI(up)  =PMLc(1:6)%ZE(up)  +100
      !
      IF (.not.(sgg%Border%IsLeftPML))  PMLc(1:6)%ZI(left) =PMLc(1:6)%ZE(left) +100
      IF (.not.(sgg%Border%IsRightPML)) PMLc(1:6)%ZI(right)=PMLc(1:6)%ZE(right)+100
      !
      IF (.not.(sgg%Border%IsFrontPML)) PMLc(1:6)%ZI(front)=PMLc(1:6)%ZE(front)+100
      IF (.not.(sgg%Border%IsBackPML))  PMLc(1:6)%ZI(back) =PMLc(1:6)%ZE(back) +100





      !PML Field component matrix allocation
      do REGION=left,right
         allocate (regLR(region)%Psi_Exy(PMLc(iEx)%XI(region) : PMLc(iEx)%XE(region), &
         PMLc(iEx)%YI(region) : PMLc(iEx)%YE(region), &
         PMLc(iEx)%ZI(region) : PMLc(iEx)%ZE(region)),&
         regLR(region)%Psi_Ezy(PMLc(iEz)%XI(region) : PMLc(iEz)%XE(region), &
         PMLc(iEz)%YI(region) : PMLc(iEz)%YE(region), &
         PMLc(iEz)%ZI(region) : PMLc(iEz)%ZE(region)),&
         regLR(region)%Psi_Hxy(PMLc(iHx)%XI(region) : PMLc(iHx)%XE(region), &
         PMLc(iHx)%YI(region) : PMLc(iHx)%YE(region), &
         PMLc(iHx)%ZI(region) : PMLc(iHx)%ZE(region)),&
         regLR(region)%Psi_Hzy(PMLc(iHz)%XI(region) : PMLc(iHz)%XE(region), &
         PMLc(iHz)%YI(region) : PMLc(iHz)%YE(region), &
         PMLc(iHz)%ZI(region) : PMLc(iHz)%ZE(region)))
         
         if (.not.resume) then
            regLR(REGION)%Psi_Exy=0.0_RKIND  ;  regLR(REGION)%Psi_Ezy=0.0_RKIND ;  regLR(REGION)%Psi_Hxy=0.0_RKIND ; regLR(REGION)%Psi_Hzy=0.0_RKIND ;
         else
            Do k=PMLc(iEx)%ZI(region),PMLc(iEx)%ZE(region)
               Do j=PMLc(iEx)%YI(region),PMLc(iEx)%YE(region)
                  READ (14) (regLR(region)%Psi_Exy(i,j,k),i=PMLc(iEx)%XI(region),PMLc(iEx)%XE(region))
               End do
            End do
            Do k=PMLc(iEz)%ZI(region),PMLc(iEz)%ZE(region)
               Do j=PMLc(iEz)%YI(region),PMLc(iEz)%YE(region)
                  READ (14) ( regLR(region)%Psi_Ezy(i,j,k),i=PMLc(iEz)%XI(region),PMLc(iEz)%XE(region))
               End do
            End do
            Do k=PMLc(iHx)%ZI(region),PMLc(iHx)%ZE(region)
               Do j=PMLc(iHx)%YI(region),PMLc(iHx)%YE(region)
                  READ (14) (regLR(region)%Psi_Hxy(i,j,k),i=PMLc(iHx)%XI(region),PMLc(iHx)%XE(region))
               End do
            End do
            Do k=PMLc(iHz)%ZI(region),PMLc(iHz)%ZE(region)
               Do j=PMLc(iHz)%YI(region),PMLc(iHz)%YE(region)
                  READ (14) ( regLR(region)%Psi_Hzy(i,j,k),i=PMLc(iHz)%XI(region),PMLc(iHz)%XE(region))
               End do
            End do
         endif
      end do
      do REGION=down,up
         allocate (regDU(region)%Psi_Eyz(PMLc(iEy)%XI(region) : PMLc(iEy)%XE(region), &
         PMLc(iEy)%YI(region) : PMLc(iEy)%YE(region), &
         PMLc(iEy)%ZI(region) : PMLc(iEy)%ZE(region)),&
         regDU(region)%Psi_Exz(PMLc(iEx)%XI(region) : PMLc(iEx)%XE(region), &
         PMLc(iEx)%YI(region) : PMLc(iEx)%YE(region), &
         PMLc(iEx)%ZI(region) : PMLc(iEx)%ZE(region)),&
         regDU(region)%Psi_Hyz(PMLc(iHy)%XI(region) : PMLc(iHy)%XE(region), &
         PMLc(iHy)%YI(region) : PMLc(iHy)%YE(region), &
         PMLc(iHy)%ZI(region) : PMLc(iHy)%ZE(region)),&
         regDU(region)%Psi_Hxz(PMLc(iHx)%XI(region) : PMLc(iHx)%XE(region), &
         PMLc(iHx)%YI(region) : PMLc(iHx)%YE(region), &
         PMLc(iHx)%ZI(region) : PMLc(iHx)%ZE(region)))
         if (.not.resume) then
            regDU(REGION)%Psi_Eyz=0.0_RKIND  ;  regDU(REGION)%Psi_Exz=0.0_RKIND ;  regDU(REGION)%Psi_Hyz=0.0_RKIND ; regDU(REGION)%Psi_Hxz=0.0_RKIND ;
         else
            Do k=PMLc(iEy)%ZI(region),PMLc(iEy)%ZE(region)
               Do j=PMLc(iEy)%YI(region),PMLc(iEy)%YE(region)
                  READ (14) (regDU(region)%Psi_Eyz(i,j,k),i=PMLc(iEy)%XI(region),PMLc(iEy)%XE(region))
               End do
            End do
            Do k=PMLc(iEx)%ZI(region),PMLc(iEx)%ZE(region)
               Do j=PMLc(iEx)%YI(region),PMLc(iEx)%YE(region)
                  READ (14) ( regDU(region)%Psi_Exz(i,j,k),i=PMLc(iEx)%XI(region),PMLc(iEx)%XE(region))
               End do
            End do
            Do k=PMLc(iHy)%ZI(region),PMLc(iHy)%ZE(region)
               Do j=PMLc(iHy)%YI(region),PMLc(iHy)%YE(region)
                  READ (14) (regDU(region)%Psi_Hyz(i,j,k),i=PMLc(iHy)%XI(region),PMLc(iHy)%XE(region))
               End do
            End do
            Do k=PMLc(iHx)%ZI(region),PMLc(iHx)%ZE(region)
               Do j=PMLc(iHx)%YI(region),PMLc(iHx)%YE(region)
                  READ (14) ( regDU(region)%Psi_Hxz(i,j,k),i=PMLc(iHx)%XI(region),PMLc(iHx)%XE(region))
               End do
            End do
         endif
      end do
      do REGION=back,front
         allocate (regBF(region)%Psi_Ezx(PMLc(iEz)%XI(region) : PMLc(iEz)%XE(region), &
         PMLc(iEz)%YI(region) : PMLc(iEz)%YE(region), &
         PMLc(iEz)%ZI(region) : PMLc(iEz)%ZE(region)),&
         regBF(region)%Psi_Eyx(PMLc(iEy)%XI(region) : PMLc(iEy)%XE(region), &
         PMLc(iEy)%YI(region) : PMLc(iEy)%YE(region), &
         PMLc(iEy)%ZI(region) : PMLc(iEy)%ZE(region)),&
         regBF(region)%Psi_Hzx(PMLc(iHz)%XI(region) : PMLc(iHz)%XE(region), &
         PMLc(iHz)%YI(region) : PMLc(iHz)%YE(region), &
         PMLc(iHz)%ZI(region) : PMLc(iHz)%ZE(region)),&
         regBF(region)%Psi_Hyx(PMLc(iHy)%XI(region) : PMLc(iHy)%XE(region), &
         PMLc(iHy)%YI(region) : PMLc(iHy)%YE(region), &
         PMLc(iHy)%ZI(region) : PMLc(iHy)%ZE(region)))
         if (.not.resume) then
            regBF(REGION)%Psi_Ezx=0.0_RKIND  ;  regBF(REGION)%Psi_Eyx=0.0_RKIND ;  regBF(REGION)%Psi_Hzx=0.0_RKIND ; regBF(REGION)%Psi_Hyx=0.0_RKIND ;
         else
            Do k=PMLc(iEz)%ZI(region),PMLc(iEz)%ZE(region)
               Do j=PMLc(iEz)%YI(region),PMLc(iEz)%YE(region)
                  READ (14) (regBF(region)%Psi_Ezx(i,j,k),i=PMLc(iEz)%XI(region),PMLc(iEz)%XE(region))
               End do
            End do
            Do k=PMLc(iEy)%ZI(region),PMLc(iEy)%ZE(region)
               Do j=PMLc(iEy)%YI(region),PMLc(iEy)%YE(region)
                  READ (14) ( regBF(region)%Psi_Eyx(i,j,k),i=PMLc(iEy)%XI(region),PMLc(iEy)%XE(region))
               End do
            End do
            Do k=PMLc(iHz)%ZI(region),PMLc(iHz)%ZE(region)
               Do j=PMLc(iHz)%YI(region),PMLc(iHz)%YE(region)
                  READ (14) (regBF(region)%Psi_Hzx(i,j,k),i=PMLc(iHz)%XI(region),PMLc(iHz)%XE(region))
               End do
            End do
            Do k=PMLc(iHy)%ZI(region),PMLc(iHy)%ZE(region)
               Do j=PMLc(iHy)%YI(region),PMLc(iHy)%YE(region)
                  READ (14) ( regBF(region)%Psi_Hyx(i,j,k),i=PMLc(iHy)%XI(region),PMLc(iHy)%XE(region))
               End do
            End do
         endif
      end do

      return
   end subroutine InitCPMLBorders




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Flush the PML data to disk for resuming purposes
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine StoreFieldsCPMLBorders

      integer (kind=4)  ::  region,i,j,k


      do REGION=left,right
         Do k=PMLc(iEx)%ZI(region),PMLc(iEx)%ZE(region)
            Do j=PMLc(iEx)%YI(region),PMLc(iEx)%YE(region)
               write(14,err=634) (regLR(region)%Psi_Exy(i,j,k),i=PMLc(iEx)%XI(region),PMLc(iEx)%XE(region))
            End do
         End do
         Do k=PMLc(iEz)%ZI(region),PMLc(iEz)%ZE(region)
            Do j=PMLc(iEz)%YI(region),PMLc(iEz)%YE(region)
               write(14,err=634) ( regLR(region)%Psi_Ezy(i,j,k),i=PMLc(iEz)%XI(region),PMLc(iEz)%XE(region))
            End do
         End do
         Do k=PMLc(iHx)%ZI(region),PMLc(iHx)%ZE(region)
            Do j=PMLc(iHx)%YI(region),PMLc(iHx)%YE(region)
               write(14,err=634) (regLR(region)%Psi_Hxy(i,j,k),i=PMLc(iHx)%XI(region),PMLc(iHx)%XE(region))
            End do
         End do
         Do k=PMLc(iHz)%ZI(region),PMLc(iHz)%ZE(region)
            Do j=PMLc(iHz)%YI(region),PMLc(iHz)%YE(region)
               write(14,err=634) ( regLR(region)%Psi_Hzy(i,j,k),i=PMLc(iHz)%XI(region),PMLc(iHz)%XE(region))
            End do
         End do
      end do
      do REGION=down,up
         Do k=PMLc(iEy)%ZI(region),PMLc(iEy)%ZE(region)
            Do j=PMLc(iEy)%YI(region),PMLc(iEy)%YE(region)
               write(14,err=634) (regDU(region)%Psi_Eyz(i,j,k),i=PMLc(iEy)%XI(region),PMLc(iEy)%XE(region))
            End do
         End do
         Do k=PMLc(iEx)%ZI(region),PMLc(iEx)%ZE(region)
            Do j=PMLc(iEx)%YI(region),PMLc(iEx)%YE(region)
               write(14,err=634) ( regDU(region)%Psi_Exz(i,j,k),i=PMLc(iEx)%XI(region),PMLc(iEx)%XE(region))
            End do
         End do
         Do k=PMLc(iHy)%ZI(region),PMLc(iHy)%ZE(region)
            Do j=PMLc(iHy)%YI(region),PMLc(iHy)%YE(region)
               write(14,err=634) (regDU(region)%Psi_Hyz(i,j,k),i=PMLc(iHy)%XI(region),PMLc(iHy)%XE(region))
            End do
         End do
         Do k=PMLc(iHx)%ZI(region),PMLc(iHx)%ZE(region)
            Do j=PMLc(iHx)%YI(region),PMLc(iHx)%YE(region)
               write(14,err=634) ( regDU(region)%Psi_Hxz(i,j,k),i=PMLc(iHx)%XI(region),PMLc(iHx)%XE(region))
            End do
         End do
      end do
      do REGION=back,front
         Do k=PMLc(iEz)%ZI(region),PMLc(iEz)%ZE(region)
            Do j=PMLc(iEz)%YI(region),PMLc(iEz)%YE(region)
               write(14,err=634) (regBF(region)%Psi_Ezx(i,j,k),i=PMLc(iEz)%XI(region),PMLc(iEz)%XE(region))
            End do
         End do
         Do k=PMLc(iEy)%ZI(region),PMLc(iEy)%ZE(region)
            Do j=PMLc(iEy)%YI(region),PMLc(iEy)%YE(region)
               write(14,err=634) ( regBF(region)%Psi_Eyx(i,j,k),i=PMLc(iEy)%XI(region),PMLc(iEy)%XE(region))
            End do
         End do
         Do k=PMLc(iHz)%ZI(region),PMLc(iHz)%ZE(region)
            Do j=PMLc(iHz)%YI(region),PMLc(iHz)%YE(region)
               write(14,err=634) (regBF(region)%Psi_Hzx(i,j,k),i=PMLc(iHz)%XI(region),PMLc(iHz)%XE(region))
            End do
         End do
         Do k=PMLc(iHy)%ZI(region),PMLc(iHy)%ZE(region)
            Do j=PMLc(iHy)%YI(region),PMLc(iHy)%YE(region)
               write(14,err=634) ( regBF(region)%Psi_Hyx(i,j,k),i=PMLc(iHy)%XI(region),PMLc(iHy)%XE(region))
            End do
         End do
      end do

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'BORDERSPML: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsCPMLBorders


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyCPMLBorders
      integer (kind=4)  ::  region
      if (associated(sig_max)) deALLOCATE (sig_max,aPar_max,kPar_max)
      if (associated(P_ce_x )) deallocate (P_ce_x,P_ce_y,P_ce_z,P_be_x,P_be_y,P_be_z,P_cm_x,P_cm_y,P_cm_z,P_bm_x,P_bm_y,P_bm_z)
      if (associated(ce_x   )) deallocate ( ce_x, ce_y, ce_z, cm_x, cm_y, cm_z,Ice_x,Ice_y,Ice_z,Icm_x,Icm_y,Icm_z)

      do REGION=left,right
         if (associated(regLR(region)%Psi_Exy)) deallocate (regLR(region)%Psi_Exy,regLR(region)%Psi_Ezy,regLR(region)%Psi_Hxy,regLR(region)%Psi_Hzy)
      end do
      do REGION=down,up
         if (associated(regDU(region)%Psi_Eyz)) deallocate (regDU(region)%Psi_Eyz,regDU(region)%Psi_Exz,regDU(region)%Psi_Hyz,regDU(region)%Psi_Hxz)
      end do
      do REGION=back,front
         if (associated(regBF(region)%Psi_Ezx)) deallocate (regBF(region)%Psi_Ezx,regBF(region)%Psi_Eyx,regBF(region)%Psi_Hzx,regBF(region)%Psi_Hyx)
      end do

      return
   end subroutine DestroyCPMLBorders

   !**************************************************************************************************
   subroutine AdvanceelectricCPML( NumMedia, b, sggMiEx, sggMiEy, sggMiEz, g2, Ex, Ey, Ez, Hx, Hy, Hz)
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent( IN)  ::  NumMedia
      type( bounds_t), intent( IN)  ::  b
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN)  ::  sggMiEx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent( IN)  ::  sggMiEy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent( IN)  ::  sggMiEz
      !--->
      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  g2
      !--->
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hz
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT)  ::  Ez
      !---------------------------> variables locales <-----------------------------------------------
      integer (kind=4)  ::  REGION, i, j, k, medio, i_m, j_m, k_m
      !---------------------------> empieza AdvanceelectricCPML <-------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               !--->
               medio = sggMiEx( i_m , j_m , k_m )
               regLR( REGION)%Psi_Exy( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Exy( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) + G2( medio) * regLR( REGION)%Psi_Exy( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = sggMiEz( i_m , j_m , k_m )
               regLR( REGION)%Psi_Ezy( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Ezy( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) - G2( medio) * regLR( REGION)%Psi_Ezy( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION =  right
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               !--->
               medio = sggMiEx( i_m , j_m , k_m )
               regLR( REGION)%Psi_Exy( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Exy( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) + G2( medio) * regLR( REGION)%Psi_Exy( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = sggMiEz( i_m , j_m , k_m )
               regLR( REGION)%Psi_Ezy( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Ezy( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) - G2( medio) * regLR( REGION)%Psi_Ezy( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION), PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION),PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = sggMiEy( i_m , j_m , k_m )
               regDU( REGION)%Psi_Eyz( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Eyz( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) + G2( medio) * regDU( REGION)%Psi_Eyz( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               medio = sggMiEx( i_m , j_m , k_m )
               regDU( REGION)%Psi_Exz( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Exz( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) - G2( medio) * regDU( REGION)%Psi_Exz( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION), PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION),PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = sggMiEy( i_m , j_m , k_m )
               regDU( REGION)%Psi_Eyz( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Eyz( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) + G2( medio) * regDU( REGION)%Psi_Eyz( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               medio = sggMiEx( i_m , j_m , k_m )
               regDU( REGION)%Psi_Exz( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Exz( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) - G2( medio) * regDU( REGION)%Psi_Exz( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif




      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = sggMiEz( i_m , j_m , k_m )
               regBF( REGION)%Psi_Ezx( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Ezx( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) + G2( medio) * regBF( REGION)%Psi_Ezx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION) ,PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION), PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = sggMiEy( i_m , j_m , k_m )
               regBF( REGION)%Psi_Eyx( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Eyx( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) - G2( medio) * regBF( REGION)%Psi_Eyx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = sggMiEz( i_m , j_m , k_m )
               regBF( REGION)%Psi_Ezx( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Ezx( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) + G2( medio) * regBF( REGION)%Psi_Ezx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION) ,PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION), PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = sggMiEy( i_m , j_m , k_m )
               regBF( REGION)%Psi_Eyx( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Eyx( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) - G2( medio) * regBF( REGION)%Psi_Eyx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !---------------------------> acaba AdvanceelectricCPML <---------------------------------------
      return
   endsubroutine AdvanceelectricCPML
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Advances the magnetic field in the PML
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine AdvanceMagneticCPML( NumMedia, b, sggMiHx, sggMiHy, sggMiHz, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent( IN)  ::  NumMedia
      type( bounds_t), intent( IN)  ::  b
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN)  ::  sggMiHx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN)  ::  sggMiHy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN)  ::  sggMiHz
      !--->
      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  gm2
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer (kind=4)  ::  REGION, i, j, k, medio, i_m, j_m, k_m
      !---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
      !Hetic Fields PML Zone
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
               medio = sggMiHx( i_m , j_m , k_m )
               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2(medio)*regLR(REGION)%Psi_Hxy(i,j,k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
               medio = sggMiHz( i_m , j_m , k_m )
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2( medio) * regLR( REGION)%Psi_Hzy( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REGION = right
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
               medio = sggMiHx( i_m , j_m , k_m )
               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2(medio)*regLR(REGION)%Psi_Hxy(i,j,k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
               medio = sggMiHz( i_m , j_m , k_m )
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2( medio) * regLR( REGION)%Psi_Hzy( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
               medio = sggMiHy( i_m , j_m , k_m )
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2( medio) * regDU( REGION)%Psi_Hyz( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
               medio = sggMiHx( i_m , j_m , k_m )
               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2(medio) * regDU( REGION)%Psi_Hxz( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
               medio = sggMiHy( i_m , j_m , k_m )
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2( medio) * regDU( REGION)%Psi_Hyz( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
               medio = sggMiHx( i_m , j_m , k_m )
               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2(medio) * regDU( REGION)%Psi_Hxz( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION=back
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
               medio = sggMiHz( i_m , j_m , k_m )
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2( medio) * regBF( REGION)%Psi_Hzx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
               medio = sggMiHy( i_m , j_m , k_m )
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2( medio) * regBF( REGION)%Psi_Hyx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION=front
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
               medio = sggMiHz( i_m , j_m , k_m )
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2( medio) * regBF( REGION)%Psi_Hzx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
               medio = sggMiHy( i_m , j_m , k_m )
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2( medio) * regBF( REGION)%Psi_Hyx( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
      return
   endsubroutine AdvanceMagneTicCPML
!!!
!!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   !!! Advances the magnetic field in the PML
!!!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!   subroutine FreeSpace_AdvanceMagneticCPML( NumMedia, b, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
!!!      !---------------------------> inputs <----------------------------------------------------------
!!!      integer, intent( IN)  ::  NumMedia
!!!      type( bounds_t), intent( IN)  ::  b
!!!      !--->
!!!      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  gm2
!!!      !--->
!!!      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
!!!      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
!!!      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
!!!      !---------------------------> inputs/outputs <--------------------------------------------------
!!!      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
!!!      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
!!!      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
!!!      !---------------------------> variables locales <-----------------------------------------------
!!!      integer (kind=4)  ::  REGION, i, j, k, i_m, j_m, k_m
!!!      real (kind = RKIND)  ::  GM2_1
!!!      GM2_1=GM2(1)
!!!      !---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
!!!      !Hetic Fields PML Zone
!!!      !
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      REGION = left
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
!!!         k_m = k - b%Hx%ZI
!!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
!!!            j_m = j - b%Hx%YI
!!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
!!!               i_m = i - b%Hx%XI
!!!               !--->
!!!               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
!!!               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
!!!               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2_1*regLR(REGION)%Psi_Hxy(i,j,k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
!!!         k_m = k - b%Hz%ZI
!!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
!!!            j_m = j - b%Hz%YI
!!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
!!!               i_m = i - b%Hz%XI
!!!               !--->
!!!               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
!!!               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
!!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2_1 * regLR( REGION)%Psi_Hzy( i, j, k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!!      REGION = right
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
!!!         k_m = k - b%Hx%ZI
!!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
!!!            j_m = j - b%Hx%YI
!!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
!!!               i_m = i - b%Hx%XI
!!!               !--->
!!!               regLR( REGION)%Psi_Hxy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxy( i, j, k) +  &
!!!               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
!!!               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2_1*regLR(REGION)%Psi_Hxy(i,j,k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
!!!         k_m = k - b%Hz%ZI
!!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
!!!            j_m = j - b%Hz%YI
!!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
!!!               i_m = i - b%Hz%XI
!!!               !--->
!!!               regLR( REGION)%Psi_Hzy( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzy( i, j, k) +  &
!!!               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
!!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2_1 * regLR( REGION)%Psi_Hzy( i, j, k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      REGION = down
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
!!!         k_m = k - b%Hy%ZI
!!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
!!!            j_m = j - b%Hy%YI
!!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
!!!               i_m = i - b%Hy%XI
!!!               !--->
!!!               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
!!!               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
!!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2_1 * regDU( REGION)%Psi_Hyz( i, j, k)
!!!            enddo !bucle i
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
!!!         k_m = k - b%Hx%ZI
!!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
!!!            j_m = j - b%Hx%YI
!!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
!!!               i_m = i - b%Hx%XI
!!!               !--->
!!!               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
!!!               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
!!!               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2_1 * regDU( REGION)%Psi_Hxz( i, j, k)
!!!            enddo !bucle i
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      REGION = up
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
!!!         k_m = k - b%Hy%ZI
!!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
!!!            j_m = j - b%Hy%YI
!!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
!!!               i_m = i - b%Hy%XI
!!!               !--->
!!!               regDU( REGION)%Psi_Hyz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyz( i, j, k) +  &
!!!               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
!!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2_1 * regDU( REGION)%Psi_Hyz( i, j, k)
!!!            enddo !bucle i
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
!!!         k_m = k - b%Hx%ZI
!!!         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
!!!            j_m = j - b%Hx%YI
!!!            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
!!!               i_m = i - b%Hx%XI
!!!               !--->
!!!               regDU( REGION)%Psi_Hxz( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxz( i, j, k) +  &
!!!               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
!!!               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2_1 * regDU( REGION)%Psi_Hxz( i, j, k)
!!!            enddo !bucle i
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      REGION=back
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
!!!         k_m = k - b%Hz%ZI
!!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
!!!            j_m = j - b%Hz%YI
!!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
!!!               i_m = i - b%Hz%XI
!!!               !--->
!!!               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
!!!               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
!!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2_1 * regBF( REGION)%Psi_Hzx( i, j, k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
!!!         k_m = k - b%Hy%ZI
!!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
!!!            j_m = j - b%Hy%YI
!!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
!!!               i_m = i - b%Hy%XI
!!!               !--->
!!!               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
!!!               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
!!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2_1 * regBF( REGION)%Psi_Hyx( i, j, k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!      REGION=front
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
!!!         k_m = k - b%Hz%ZI
!!!         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
!!!            j_m = j - b%Hz%YI
!!!            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
!!!               i_m = i - b%Hz%XI
!!!               !--->
!!!               regBF( REGION)%Psi_Hzx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzx( i, j, k) +  &
!!!               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
!!!               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2_1 * regBF( REGION)%Psi_Hzx( i, j, k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!#ifdef CompileWithOpenMP
!!!!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
!!!#endif
!!!      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
!!!         k_m = k - b%Hy%ZI
!!!         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
!!!            j_m = j - b%Hy%YI
!!!            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
!!!               i_m = i - b%Hy%XI
!!!               !--->
!!!               regBF( region)%Psi_Hyx( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyx( i, j, k) +  &
!!!               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
!!!               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2_1 * regBF( REGION)%Psi_Hyx( i, j, k)
!!!            enddo
!!!         enddo
!!!      enddo
!!!#ifdef CompileWithOpenMP
!!!!$OMP END PARALLEL DO
!!!#endif
!!!
!!!
!!!      !---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
!!!      return
!!!   endsubroutine FreeSpace_AdvanceMagneTicCPML


subroutine calc_cpmlconstants(sgg, Idxe,Idye,Idze,Idxh,Idyh,Idzh,eps00,mu00)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      REAL (KIND=RKIND), intent (in)  ::  eps00,mu00
      REAL (KIND=RKIND) , dimension (:)   , intent(inout)      ::  &
      Idxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      Idye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      Idze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)
      integer :: i,j,k,o,p
      REAL (KIND=RKIND) :: del,sigmae,kpare,apare,sigmam,kparm,aparm
      character(len=BUFSIZE) :: buff
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)

      del=1.0 !!!una simple inicializacion para que gfortran no se queje
      !Find the maximum conductivity for each direcion o=1,2,3 and for the starting and ending layer p=1,2

      sig_max=0.0_RKIND; aPar_max=0.0_RKIND; kPar_max=0.0_RKIND;
      do o=1,3
         do p=1,2
            if ((o == 1).and.(p == 1)) del=dxe(sgg%ALLOC(iEx)%XI)
            if ((o == 1).and.(p == 2)) del=dxe(sgg%ALLOC(iEx)%XE)
            if ((o == 2).and.(p == 1)) del=dye(sgg%ALLOC(iEy)%YI)
            if ((o == 2).and.(p == 2)) del=dye(sgg%ALLOC(iEy)%YE)
            if ((o == 3).and.(p == 1)) del=dze(sgg%ALLOC(iEz)%ZI)
            if ((o == 3).and.(p == 2)) del=dze(sgg%ALLOC(iEz)%ZE)
            if (sgg%PML%NumLayers(o,p) /= 0) then
               if ((sgg%PML%NumLayers(o,p) == 10).or.(sgg%PML%NumLayers(o,p) == 5)) then
                  !gedney sigma optimo no tiene en cuenta el refle (taflove 3 ed, pag 294)
                  !para 5 celdas es exp(-8) y para 16 celdas es exp(-16)
                  sig_max(o,p)=0.8*(sgg%PML%orden(o,p)+1)/(sqrt(Mu0/eps0)*del) !cambio tonto no afecta a nada 260919
               else
                  sig_max(o,p)=-((log(sgg%PML%CoeffReflPML(o,p))*(sgg%PML%orden(o,p)+1))/ &
                  (2*sqrt(Mu0/eps0)*sgg%PML%NumLayers(o,p)*del))
                  !ojo LO SIGUIENTE  estaba maaaaallllllllll porque NO ES MATERIAL INDEPENDENT
                  !                          (2*sqrt(Mu0/eps0)*sqrt(sgg%Med(jmed)%Epr*sgg%Med(jmed)%Mur)*sgg%PML%NumLayers(o,p)*del))
                  !los multilayer petan- !!! !Viene de Gedney, pero esta maaaalllll!! !corregido 20marzo 2011
               endif
            else
               sig_max(o,p)=1.0e29_RKIND
            endif
         end do
      end do

      !readjust relatively the alphamaxpar to the maximum conductivity
      do o=1,3
         do p=1,2
            if ((o == 1).and.(p == 1)) del=dxe(sgg%ALLOC(iEx)%XI)
            if ((o == 1).and.(p == 2)) del=dxe(sgg%ALLOC(iEx)%XE)
            if ((o == 2).and.(p == 1)) del=dye(sgg%ALLOC(iEy)%YI)
            if ((o == 2).and.(p == 2)) del=dye(sgg%ALLOC(iEy)%YE)
            if ((o == 3).and.(p == 1)) del=dze(sgg%ALLOC(iEz)%ZI)
            if ((o == 3).and.(p == 2)) del=dze(sgg%ALLOC(iEz)%ZE)
            if (sgg%PML%NumLayers(o,p) /= 0) then
               aPar_max(o,p)=alphamaxpar * sig_max(o,p)
               kPar_max(o,p)=kappamaxpar
            else
               aPar_max(o,p)=0.0_RKIND
               kPar_max(o,p)=1.0_RKIND
            endif
         end do
      end do

      !

      P_ce_x=0.0_RKIND; P_ce_y=0.0_RKIND; P_ce_z=0.0_RKIND; P_cm_x=0.0_RKIND; P_cm_y=0.0_RKIND; P_cm_z=0.0_RKIND;
      P_be_x=1.0_RKIND; P_be_y=1.0_RKIND; P_be_z=1.0_RKIND; P_bm_x=1.0_RKIND; P_bm_y=1.0_RKIND; P_bm_z=1.0_RKIND;
      !Calculate the coefficients (in CPML they are the same for every possible medium OJOOOOOOOOOOOOOOOOO)

      !Default !2011 already done in main
      !Idxh        ( : )=    1.0_RKIND / dxh( : )
      !Idyh        ( : )=    1.0_RKIND / dyh( : )
      !Idzh        ( : )=    1.0_RKIND / dzh( : )
      !Idxe        ( : )=    1.0_RKIND / dxe( : )
      !Idye        ( : )=    1.0_RKIND / dye( : )
      !Idze        ( : )=    1.0_RKIND / dze( : )
      ! Calculate
      do i=sgg%ALLOC(iEx)%XI,sgg%ALLOC(iEx)%XE
         if (i <= SINPML_Fullsize(iHx)%XI-1) then !Back
            if ((sgg%PML%orden(1,1) == 0)) then
               Sigmae=    Sig_max(1,1)
               kPare=1.0_RKIND+(kPar_max(1,1)-1)
            else
               Sigmae=    Sig_max(1,1)   * ce_x(i)**sgg%PML%orden(1,1)
               kPare=1.0_RKIND+(kPar_max(1,1)-1)* ce_x(i)**sgg%PML%orden(1,1)
            endif
            aPare=    aPar_max(1,1)*Ice_x(i)**alphaOrden !!  **sgg%PML%orden(1,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente !gedney lo escala linealmente
            P_be_x(i)=exp(-(sigmae/kPare+aPare)*sgg%dt/Eps0)
            P_ce_x(i)=(sigmae*(P_be_x(i)-1.0_RKIND)/(sigmae+kPare*aPare)/kpare)/dxh(i)
            Idxh(i)=1.0_RKIND / (kPare*dxh(i))
         elseif (i >= SINPML_Fullsize(iHx)%XE +1) then  !Front
            if ((sgg%PML%orden(1,2) == 0)) then
               Sigmae=    Sig_max(1,2)
               kPare=1.0_RKIND+(kPar_max(1,2)-1)
            else
               Sigmae=    Sig_max(1,2)   * ce_x(i)**sgg%PML%orden(1,2)
               kPare=1.0_RKIND+(kPar_max(1,2)-1)* ce_x(i)**sgg%PML%orden(1,2)
            endif
            aPare=    aPar_max(1,2)*Ice_x(i)**alphaOrden !!  **sgg%PML%orden(1,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_be_x(i)=exp(-(sigmae/kPare+aPare)*sgg%dt/Eps0)
            P_ce_x(i)=(sigmae*(P_be_x(i)-1.0_RKIND)/(sigmae+kPare*aPare)/kpare)/dxh(i)
            Idxh(i)=1.0_RKIND / (kPare*dxh(i))
         endif
      end do
      do j=sgg%ALLOC(iEy)%YI,sgg%ALLOC(iEy)%YE
         if (j <= SINPML_Fullsize(iHy)%YI-1) then !Left
            if ((sgg%PML%orden(2,1) == 0)) then
               Sigmae=    Sig_max(2,1)
               kPare=1.0_RKIND+(kPar_max(2,1)-1)
            else
               Sigmae=    Sig_max(2,1)   * ce_y(j)**sgg%PML%orden(2,1)
               kPare=1.0_RKIND+(kPar_max(2,1)-1)* ce_y(j)**sgg%PML%orden(2,1)
            endif
            aPare=    aPar_max(2,1)*Ice_y(j)**alphaOrden !!  **sgg%PML%orden(2,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_be_y     (j)=exp(-(sigmae/kPare+aPare)*sgg%dt/Eps0)
            P_ce_y     (j)=(sigmae*(P_be_y(j)-1.0_RKIND)/(sigmae+kPare*aPare)/kpare)/dyh(j)
            IdYh(j) =1.0_RKIND / (kPare*dyh(j))
         elseif (j >= SINPML_Fullsize(iHy)%YE +1) then  !Right
            if ((sgg%PML%orden(2,2) == 0)) then
               Sigmae=    Sig_max(2,2)
               kPare=1.0_RKIND+(kPar_max(2,2)-1)
            else
               Sigmae=    Sig_max(2,2)   * ce_y(j)**sgg%PML%orden(2,2)
               kPare=1.0_RKIND+(kPar_max(2,2)-1)* ce_y(j)**sgg%PML%orden(2,2)
            endif
            aPare=    aPar_max(2,2)*Ice_y(j)**alphaOrden !!  **sgg%PML%orden(2,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_be_y     (j)=exp(-(sigmae/kPare+aPare)*sgg%dt/Eps0)
            P_ce_y     (j)=(sigmae*(P_be_y(j)-1.0_RKIND)/(sigmae+kPare*aPare)/kpare)/dyh(j)
            IdYh(j) =1.0_RKIND / (kPare*dyh(j))
         endif
      end do
      do k=sgg%ALLOC(iEz)%ZI,sgg%ALLOC(iEz)%ZE
         if (k <= SINPML_Fullsize(iHz)%ZI-1) then !Down
            if ((sgg%PML%orden(3,1) == 0)) then
               sigmae=    Sig_max(3,1)
               kPare=1.0_RKIND+(kPar_max(3,1)-1)
            else
               sigmae=    Sig_max(3,1)   * ce_z(k)**sgg%PML%orden(3,1)
               kPare=1.0_RKIND+(kPar_max(3,1)-1)* ce_z(k)**sgg%PML%orden(3,1)
            endif
            aPare=    aPar_max(3,1)*Ice_z(k)**alphaOrden !!  **sgg%PML%orden(3,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_be_z     (k)=exp(-(sigmae/kPare+aPare)*sgg%dt/Eps0)
            P_ce_z     (k)=(sigmae*(P_be_z(k)-1.0_RKIND)/(sigmae+kPare*aPare)/kpare)/dzh(k)
            Idzh(k)=1.0_RKIND / (kPare*dzh(k))
         elseif (k >= SINPML_Fullsize(iHz)%ZE +1) then  !Up
            if ((sgg%PML%orden(3,2) == 0)) then
               sigmae=    Sig_max(3,2)
               kPare=1.0_RKIND+(kPar_max(3,2)-1)
            else
               sigmae=    Sig_max(3,2)   * ce_z(k)**sgg%PML%orden(3,2)
               kPare=1.0_RKIND+(kPar_max(3,2)-1)* ce_z(k)**sgg%PML%orden(3,2)
            endif
            aPare=    aPar_max(3,2)*Ice_z(k)**alphaOrden !!  **sgg%PML%orden(3,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_be_z     (k)=exp(-(sigmae/kPare+aPare)*sgg%dt/Eps0)
            P_ce_z     (k)=(sigmae*(P_be_z(k)-1.0_RKIND)/(sigmae+kPare*aPare)/kpare)/dzh(k)
            Idzh(k)=1.0_RKIND / (kPare*dzh(k))
         endif
      end do
      !magnetic
      do i=sgg%ALLOC(iHx)%XI,sgg%ALLOC(iHx)%XE
         if (i <= SINPML_Fullsize(iHx)%XI-1) then !back
            if ((sgg%PML%orden(1,1) == 0)) then
               Sigmam=    Sig_max(1,1)
               kParm=1.0_RKIND+(kPar_max(1,1)-1)
            else
               Sigmam=    Sig_max(1,1)   * cm_x(i)**sgg%PML%orden(1,1)
               kParm=1.0_RKIND+(kPar_max(1,1)-1)* cm_x(i)**sgg%PML%orden(1,1)
            endif
            aParm=    aPar_max(1,1)*Icm_x(i)**alphaOrden !!  **sgg%PML%orden(1,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_x(i)=exp(-(sigmam/kParm+aParm)*sgg%dt/Eps0)
            P_cm_x(i)=(sigmam*(P_bm_x(i)-1.0_RKIND)/(sigmam+kParm*aParm)/kparm)/dxe(i)
            Idxe(i)=1.0_RKIND / (kParm*dxe(i))
            !
            WRITE (buff,'(a,i4,a,5e9.2e2)') 'back(',i,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            (sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            !IF ((sgg%Border%IsBackPML).and.(i>sgg%ALLOC(iHx)%XI)) CALL print11 (layoutnumber, buff)
         elseif (i >= SINPML_Fullsize(iHx)%XE )  then !front
            if ((sgg%PML%orden(1,2) == 0)) then
               Sigmam=    Sig_max(1,2)
               kParm=1.0_RKIND+(kPar_max(1,2)-1)
            else
               Sigmam=    Sig_max(1,2)   * cm_x(i)**sgg%PML%orden(1,2)
               kParm=1.0_RKIND+(kPar_max(1,2)-1)* cm_x(i)**sgg%PML%orden(1,2)
            endif
            aParm=    aPar_max(1,2)*Icm_x(i)**alphaOrden !!  **sgg%PML%orden(1,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_x(i)=exp(-(sigmam/kParm+aParm)*sgg%dt/Eps0)
            P_cm_x(i)=(sigmam*(P_bm_x(i)-1.0_RKIND)/(sigmam+kParm*aParm)/kparm)/dxe(i)
            Idxe(i)=1.0_RKIND / (kParm*dxe(i))
            !
            !WRITE (buff,'(a,i4,a,5e9.2e2)') 'front(',i,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            !IF ((sgg%Border%IsFrontPML).and.(i<sgg%ALLOC(iHx)%XE-1))  CALL print11 (layoutnumber, buff)
         endif
      end do
      do j=sgg%ALLOC(iHy)%YI,sgg%ALLOC(iHy)%YE
         if (j <= SINPML_Fullsize(iHy)%YI-1) then !Left
            if ((sgg%PML%orden(2,1) == 0)) then
               Sigmam=    Sig_max(2,1)
               kParm=1.0_RKIND+(kPar_max(2,1)-1)
            else
               Sigmam=    Sig_max(2,1)   * cm_y(j)**sgg%PML%orden(2,1)
               kParm=1.0_RKIND+(kPar_max(2,1)-1)* cm_y(j)**sgg%PML%orden(2,1)
            endif
            aParm=    aPar_max(2,1)*Icm_y(j)**alphaOrden !!  **sgg%PML%orden(2,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_y(j)=exp(-(sigmam/kParm+aParm)*sgg%dt/Eps0)
            P_cm_y     (j)=(sigmam*(P_bm_y(j)-1.0_RKIND)/(sigmam+kParm*aParm)/kparm)/dye(j)
            Idye(j)=1.0_RKIND / (kParm*dye(j))
            !
            !WRITE (buff,'(a,i4,a,5e9.2e2)') 'left(',j,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            !IF ((sgg%Border%IsLeftPML).and.(j>sgg%ALLOC(iHy)%YI)) CALL print11 (layoutnumber, buff)
         elseif (j >= SINPML_Fullsize(iHy)%YE ) then  !Right
            if ((sgg%PML%orden(2,2) == 0)) then
               Sigmam=    Sig_max(2,2)
               kParm=1.0_RKIND+(kPar_max(2,2)-1)
            else
               Sigmam=    Sig_max(2,2)   * cm_y(j)**sgg%PML%orden(2,2)
               kParm=1.0_RKIND+(kPar_max(2,2)-1)* cm_y(j)**sgg%PML%orden(2,2)
            endif
            aParm=    aPar_max(2,2)*Icm_y(j)**alphaOrden !!  **sgg%PML%orden(2,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_y(j)=exp(-(sigmam/kParm+aParm)*sgg%dt/Eps0)
            P_cm_y     (j)=(sigmam*(P_bm_y(j)-1.0_RKIND)/(sigmam+kParm*aParm)/kparm)/dye(j)
            Idye(j)=1.0_RKIND / (kParm*dye(j))
            !
            !WRITE (buff,'(a,i4,a,5e9.2e2)') 'right(',j,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            !IF ((sgg%Border%IsRightPML).and.(j<sgg%ALLOC(iHy)%YE-1)) CALL print11 (layoutnumber, buff)
         endif
      end do
      do k=sgg%ALLOC(iHz)%ZI,sgg%ALLOC(iHz)%ZE
         if (k <= SINPML_Fullsize(iHz)%ZI-1) then !Down
            if ((sgg%PML%orden(3,1) == 0)) then
               Sigmam=    Sig_max(3,1)
               kParm=1.0_RKIND+(kPar_max(3,1)-1)
            else
               Sigmam=    Sig_max(3,1)   * cm_z(k)**sgg%PML%orden(3,1)
               kParm=1.0_RKIND+(kPar_max(3,1)-1)* cm_z(k)**sgg%PML%orden(3,1)
            endif
            aParm=    aPar_max(3,1)*Icm_z(k)**alphaOrden !!  **sgg%PML%orden(3,1) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_z(k)=exp(-(sigmam/kParm+aParm)*sgg%dt/Eps0)
            P_cm_z     (k)=(sigmam*(P_bm_z(k)-1.0_RKIND)/(sigmam+kParm*aParm)/kparm)/dze(k)
            Idze(k)=1.0_RKIND / (kParm*dze(k))
            !
            !WRITE (buff,'(a,i4,a,5e9.2e2)') 'down(',k,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            !IF ((sgg%Border%IsDownPML).and.(k>sgg%ALLOC(iHz)%ZI)) CALL print11 (layoutnumber, buff)
         elseif (k >= SINPML_Fullsize(iHz)%ZE ) then !Up
            if ((sgg%PML%orden(3,2) == 0)) then
               Sigmam=    Sig_max(3,2)
               kParm=1.0_RKIND+(kPar_max(3,2)-1)
            else
               Sigmam=    Sig_max(3,2)   * cm_z(k)**sgg%PML%orden(3,2)
               kParm=1.0_RKIND+(kPar_max(3,2)-1)* cm_z(k)**sgg%PML%orden(3,2)
            endif
            aParm=    aPar_max(3,2)*Icm_z(k)**alphaOrden !!  **sgg%PML%orden(3,2) !!**1.0_RKIND !perfil lineal propuesto por Gedney originalmente
            P_bm_z(k)=exp(-(sigmam/kParm+aParm)*sgg%dt/Eps0)
            P_cm_z     (k)=(sigmam*(P_bm_z(k)-1.0_RKIND)/(sigmam+kParm*aParm)/kparm)/dze(k)
            Idze(k)=1.0_RKIND / (kParm*dze(k))
            !
            !WRITE (buff,'(a,i4,a,5e9.2e2)') 'up   (',k,'+d/2), A,S,FcS,FcA,refleLF=',aparm,sigmam,sigmam/(2.0_RKIND * pi*eps0),aparm/(2.0_RKIND * pi*eps0), &
            !(sqrt(kparm+sigmam/(aParm+1d-15))-1.0_RKIND)/(sqrt(kparm+sigmam/(aParm+1d-15))+1.0_RKIND)
            !IF ((sgg%Border%IsUpPML).and.(k<sgg%ALLOC(iHz)%ZE-1)) CALL print11 (layoutnumber, buff)
         endif
      end do


end subroutine calc_cpmlconstants

!!!!!!!!!!!!!!!


   !**************************************************************************************************
   subroutine AdvanceelectricCPML_freespace( NumMedia, b, sggMiEx, sggMiEy, sggMiEz, g2, Ex, Ey, Ez, Hx, Hy, Hz)
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent( IN)  ::  NumMedia
      type( bounds_t), intent( IN)  ::  b
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN)  ::  sggMiEx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent( IN)  ::  sggMiEy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent( IN)  ::  sggMiEz
      !--->
      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  g2
      !--->
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hz
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT)  ::  Ez
      !---------------------------> variables locales <-----------------------------------------------
      integer (kind=4)  ::  REGION, i, j, k, medio, i_m, j_m, k_m
      !---------------------------> empieza AdvanceelectricCPML <-------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               !--->
               medio = 1
               regLR( REGION)%Psi_Exyvac( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Exyvac( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) + G2( medio) * regLR( REGION)%Psi_Exyvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = 1
               regLR( REGION)%Psi_Ezyvac( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Ezyvac( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) - G2( medio) * regLR( REGION)%Psi_Ezyvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION =  right
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               !--->
               medio = 1
               regLR( REGION)%Psi_Exyvac( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Exyvac( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) + G2( medio) * regLR( REGION)%Psi_Exyvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = 1
               regLR( REGION)%Psi_Ezyvac( i, j, k) = P_be_y( j) * regLR( REGION)%Psi_Ezyvac( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m-1, k_m)) * P_ce_y( j)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) - G2( medio) * regLR( REGION)%Psi_Ezyvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION), PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION),PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = 1
               regDU( REGION)%Psi_Eyzvac( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Eyzvac( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) + G2( medio) * regDU( REGION)%Psi_Eyzvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               medio = 1
               regDU( REGION)%Psi_Exzvac( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Exzvac( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) - G2( medio) * regDU( REGION)%Psi_Exzvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION), PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION),PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = 1
               regDU( REGION)%Psi_Eyzvac( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Eyzvac( i, j, k) +  &
               (Hx( i_m, j_m, k_m) - Hx( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) + G2( medio) * regDU( REGION)%Psi_Eyzvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEx)%ZI( REGION), PMLc(iEx)%ZE( REGION)
         k_m = k - b%Ex%ZI
         do j = PMLc(iEx)%YI( REGION), PMLc(iEx)%YE( REGION)
            j_m = j - b%Ex%YI
            do i = PMLc(iEx)%XI( REGION), PMLc(iEx)%XE( REGION)
               i_m = i - b%Ex%XI
               medio = 1
               regDU( REGION)%Psi_Exzvac( i, j, k) = P_be_z( k) * regDU( REGION)%Psi_Exzvac( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m, j_m, k_m-1)) * P_ce_z( k)
               Ex( i_m, j_m, k_m) = Ex( i_m, j_m, k_m) - G2( medio) * regDU( REGION)%Psi_Exzvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif




      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = back
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = 1
               regBF( REGION)%Psi_Ezxvac( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Ezxvac( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) + G2( medio) * regBF( REGION)%Psi_Ezxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION) ,PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION), PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = 1
               regBF( REGION)%Psi_Eyxvac( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Eyxvac( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) - G2( medio) * regBF( REGION)%Psi_Eyxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = front
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEz)%ZI( REGION), PMLc(iEz)%ZE( REGION)
         k_m = k - b%Ez%ZI
         do j = PMLc(iEz)%YI( REGION), PMLc(iEz)%YE( REGION)
            j_m = j - b%Ez%YI
            do i = PMLc(iEz)%XI( REGION), PMLc(iEz)%XE( REGION)
               i_m = i - b%Ez%XI
               medio = 1
               regBF( REGION)%Psi_Ezxvac( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Ezxvac( i, j, k) +  &
               (Hy( i_m, j_m, k_m) - Hy( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ez( i_m, j_m, k_m) = Ez( i_m, j_m, k_m) + G2( medio) * regBF( REGION)%Psi_Ezxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iEy)%ZI( REGION), PMLc(iEy)%ZE( REGION)
         k_m = k - b%Ey%ZI
         do j = PMLc(iEy)%YI( REGION) ,PMLc(iEy)%YE( REGION)
            j_m = j - b%Ey%YI
            do i = PMLc(iEy)%XI( REGION), PMLc(iEy)%XE( REGION)
               i_m = i - b%Ey%XI
               medio = 1
               regBF( REGION)%Psi_Eyxvac( i, j, k) = P_be_x( i) * regBF( REGION)%Psi_Eyxvac( i, j, k) +  &
               (Hz( i_m, j_m, k_m) - Hz( i_m-1, j_m, k_m)) * P_ce_x( i)
               Ey( i_m, j_m, k_m) = Ey( i_m, j_m, k_m) - G2( medio) * regBF( REGION)%Psi_Eyxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !---------------------------> acaba AdvanceelectricCPML <---------------------------------------
      return
   endsubroutine AdvanceelectricCPML_freespace
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Advances the magnetic field in the PML
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine AdvanceMagneticCPML_freespace( NumMedia, b, sggMiHx, sggMiHy, sggMiHz, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
      !---------------------------> inputs <----------------------------------------------------------
      integer, intent( IN)  ::  NumMedia
      type( bounds_t), intent( IN)  ::  b
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN)  ::  sggMiHx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN)  ::  sggMiHy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN)  ::  sggMiHz
      !--->
      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  gm2
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer (kind=4)  ::  REGION, i, j, k, medio, i_m, j_m, k_m
      !---------------------------> empieza AdvanceMagneTicCPML <-------------------------------------
      !Hetic Fields PML Zone
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = left
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regLR( REGION)%Psi_Hxyvac( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxyvac( i, j, k) +  &
               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
               medio = 1
               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2(medio)*regLR(REGION)%Psi_Hxyvac(i,j,k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regLR( REGION)%Psi_Hzyvac( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzyvac( i, j, k) +  &
               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
               medio = 1
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2( medio) * regLR( REGION)%Psi_Hzyvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      REGION = right
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regLR( REGION)%Psi_Hxyvac( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hxyvac( i, j, k) +  &
               (Ez( i_m, j_m+1, k_m) - Ez( i_m, j_m, k_m)) * P_cm_y( j)
               medio = 1
               Hx( i_m, j_m, k_m)=Hx( i_m, j_m, k_m)-GM2(medio)*regLR(REGION)%Psi_Hxyvac(i,j,k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regLR( REGION)%Psi_Hzyvac( i, j, k) = P_bm_y( j) * regLR( REGION)%Psi_Hzyvac( i, j, k) +  &
               (Ex( i_m, j_m+1, k_m) - Ex( i_m, j_m, k_m)) * P_cm_y( j)
               medio = 1
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) + GM2( medio) * regLR( REGION)%Psi_Hzyvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = down
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regDU( REGION)%Psi_Hyzvac( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyzvac( i, j, k) +  &
               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
               medio = 1
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2( medio) * regDU( REGION)%Psi_Hyzvac( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regDU( REGION)%Psi_Hxzvac( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxzvac( i, j, k) +  &
               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
               medio = 1
               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2(medio) * regDU( REGION)%Psi_Hxzvac( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION = up
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regDU( REGION)%Psi_Hyzvac( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hyzvac( i, j, k) +  &
               (Ex( i_m, j_m, k_m+1) - Ex( i_m, j_m, k_m)) * P_cm_z( k)
               medio = 1
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) - GM2( medio) * regDU( REGION)%Psi_Hyzvac( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHx)%ZI( REGION), PMLc(iHx)%ZE( REGION)
         k_m = k - b%Hx%ZI
         do j = PMLc(iHx)%YI( REGION), PMLc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
            do i = PMLc(iHx)%XI( REGION), PMLc(iHx)%XE( REGION)
               i_m = i - b%Hx%XI
               !--->
               regDU( REGION)%Psi_Hxzvac( i, j, k) = P_bm_z( k) * regDU( REGION)%Psi_Hxzvac( i, j, k) +  &
               (Ey( i_m, j_m, k_m+1) - Ey( i_m, j_m, k_m)) * P_cm_z( k)
               medio = 1
               Hx( i_m, j_m, k_m) = Hx( i_m, j_m, k_m) + GM2(medio) * regDU( REGION)%Psi_Hxzvac( i, j, k)
            enddo !bucle i
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION=back
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regBF( REGION)%Psi_Hzxvac( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzxvac( i, j, k) +  &
               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
               medio = 1
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2( medio) * regBF( REGION)%Psi_Hzxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regBF( region)%Psi_Hyxvac( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyxvac( i, j, k) +  &
               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
               medio = 1
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2( medio) * regBF( REGION)%Psi_Hyxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REGION=front
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHz)%ZI( REGION), PMLc(iHz)%ZE( REGION)
         k_m = k - b%Hz%ZI
         do j = PMLc(iHz)%YI( REGION), PMLc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
            do i = PMLc(iHz)%XI( REGION), PMLc(iHz)%XE( REGION)
               i_m = i - b%Hz%XI
               !--->
               regBF( REGION)%Psi_Hzxvac( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hzxvac( i, j, k) +  &
               (Ey( i_m+1, j_m, k_m) - Ey( i_m, j_m, k_m)) * P_cm_x( i)
               medio = 1
               Hz( i_m, j_m, k_m) = Hz( i_m, j_m, k_m) - GM2( medio) * regBF( REGION)%Psi_Hzxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m,medio)
#endif
      do k = PMLc(iHy)%ZI( REGION), PMLc(iHy)%ZE( REGION)
         k_m = k - b%Hy%ZI
         do j = PMLc(iHy)%YI( REGION), PMLc(iHy)%YE( REGION)
            j_m = j - b%Hy%YI
            do i = PMLc(iHy)%XI( REGION), PMLc(iHy)%XE( REGION)
               i_m = i - b%Hy%XI
               !--->
               regBF( region)%Psi_Hyxvac( i, j, k) = P_bm_x( i) * regBF( REGION)%Psi_Hyxvac( i, j, k) +  &
               (Ez( i_m+1, j_m, k_m) - Ez( i_m, j_m, k_m)) * P_cm_x( i)
               medio = 1
               Hy( i_m, j_m, k_m) = Hy( i_m, j_m, k_m) + GM2( medio) * regBF( REGION)%Psi_Hyxvac( i, j, k)
            enddo
         enddo
      enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif


      !---------------------------> acaba AdvanceMagneTicCPML <---------------------------------------
      return
   endsubroutine AdvanceMagneTicCPML_freespace


end Module Borders_CPML


