!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MIT License
! 
! Copyright (c) 2023 University of Granada
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Current on Huygens Box
!  Taken from planewaves
!  Creation date Date :  Oct, 10, 2012
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module farfield_m

#ifdef CompileWithNF2FF

   use fdetypes
   USE REPORT

   IMPLICIT NONE
   private


   type ehxyz
      integer (kind=4)  ::  Ex=-15,Ey=-15,Ez=-15,Hx=-15,Hy=-15,Hz=-15
   end type
   type tfidaa
      type (ehxyz)  ::   com,fin,tra,fro,izq,der,aba,arr
   end type
   type ijk
      type (tfidaa)  ::  i,j,k
   end type

   type co_t
      real (kind=RKIND) x_Mx,y_Mx,z_Mx
      real (kind=RKIND) x_My,y_My,z_My
      real (kind=RKIND) x_Mz,y_Mz,z_Mz
      real (kind=RKIND) x_Jx,y_Jx,z_Jx
      real (kind=RKIND) x_Jy,y_Jy,z_Jy
      real (kind=RKIND) x_Jz,y_Jz,z_Jz
   end type

   type farfield_t
      type (ijk)  ::  TrFr,IzDe,AbAr
      logical  ::  farfieldTr,farfieldFr,farfieldIz,farfieldDe,farfieldAr,farfieldAb

      logical  ::  farfieldTr_ClonePEC_Front,farfieldTr_ClonePEC_Left,farfieldTr_ClonePEC_Right,farfieldTr_ClonePEC_Up,farfieldTr_ClonePEC_Down
      logical  ::  farfieldTr_ClonePMC_Front,farfieldTr_ClonePMC_Left,farfieldTr_ClonePMC_Right,farfieldTr_ClonePMC_Up,farfieldTr_ClonePMC_Down

      logical  ::  farfieldFr_ClonePEC_Back,farfieldFr_ClonePEC_Left,farfieldFr_ClonePEC_Right,farfieldFr_ClonePEC_Up,farfieldFr_ClonePEC_Down
      logical  ::  farfieldFr_ClonePMC_Back,farfieldFr_ClonePMC_Left,farfieldFr_ClonePMC_Right,farfieldFr_ClonePMC_Up,farfieldFr_ClonePMC_Down

      logical  ::  farfieldIz_ClonePEC_Back,farfieldIz_ClonePEC_Front,farfieldIz_ClonePEC_Right,farfieldIz_ClonePEC_Up,farfieldIz_ClonePEC_Down
      logical  ::  farfieldIz_ClonePMC_Back,farfieldIz_ClonePMC_Front,farfieldIz_ClonePMC_Right,farfieldIz_ClonePMC_Up,farfieldIz_ClonePMC_Down

      logical  ::  farfieldDe_ClonePEC_Back,farfieldDe_ClonePEC_Front,farfieldDe_ClonePEC_Left,farfieldDe_ClonePEC_Up,farfieldDe_ClonePEC_Down
      logical  ::  farfieldDe_ClonePMC_Back,farfieldDe_ClonePMC_Front,farfieldDe_ClonePMC_Left,farfieldDe_ClonePMC_Up,farfieldDe_ClonePMC_Down

      logical  ::  farfieldAr_ClonePEC_Back,farfieldAr_ClonePEC_Front,farfieldAr_ClonePEC_Left,farfieldAr_ClonePEC_Right,farfieldAr_ClonePEC_Down
      logical  ::  farfieldAr_ClonePMC_Back,farfieldAr_ClonePMC_Front,farfieldAr_ClonePMC_Left,farfieldAr_ClonePMC_Right,farfieldAr_ClonePMC_Down

      logical  ::  farfieldAb_ClonePEC_Back,farfieldAb_ClonePEC_Front,farfieldAb_ClonePEC_Left,farfieldAb_ClonePEC_Right,farfieldAb_ClonePEC_Up
      logical  ::  farfieldAb_ClonePMC_Back,farfieldAb_ClonePMC_Front,farfieldAb_ClonePMC_Left,farfieldAb_ClonePMC_Right,farfieldAb_ClonePMC_Up

      complex( kind = CKIND), dimension( :,:,:), allocatable ::  ExIz,ExDe,ExAb,ExAr,EyFr,EyTr,EyAb,EyAr,EzIz,EzDe,EzFr,EzTr
      complex( kind = CKIND), dimension( :,:,:), allocatable ::  HxIz,HxDe,HxAb,HxAr,HyFr,HyTr,HyAb,HyAr,HzIz,HzDe,HzFr,HzTr
      complex( kind = CKIND), dimension( :,:,:), allocatable ::  HxIz2,HxDe2,HxAb2,HxAr2,HyFr2,HyTr2,HyAb2,HyAr2,HzIz2,HzDe2,HzFr2,HzTr2 !to compute the scheneider geometric mean
      complex( kind = CKIND), dimension( :), allocatable  :: expIwdt,auxExp_E,auxExp_H,dftEntrada
      integer (KIND=4) :: NumFreqs,esqx1,esqx2,esqy1,esqy2,esqz1,esqz2, Ndecim
      type (coorsxyzP)  ::  Punto
      real (kind=Rkind) :: InitialFreq,FinalFreq,FreqStep,dtDecim
      REAL (KIND=RKIND)  ::  thetaStart,thetaStop,thetaStep
      REAL (KIND=RKIND)  ::  phiStart,phiStop,phiStep
      character (LEN=BUFSIZE)  ::   FileNormalize
      integer (KIND=4)    :: unitfarfield
      character (LEN=BUFSIZE)   ::  filefarfield
      REAL (KIND=RKIND)  :: XDobleAncho,YDobleAncho,ZDobleAncho
      REAL (KIND=RKIND)  :: XOffsetPlus,YOffsetPlus,ZOffsetPlus
      REAL (KIND=RKIND)  :: XOffsetMinus,YOffsetMinus,ZOffsetMinus
#ifdef CompileWithMPI
      integer (kind=4)  :: MPISubComm,MPIRoot
#endif
   end type
!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  cluz,zvac
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!
   !
   public UpdateFarField,InitFarField,Destroyfarfield,FlushFarfield,StoreFarfields
   public farfield_t
   !
   type (farfield_t), save, target :: FF

contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes Plane Wave data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitFarField(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,layoutnumber,size, &
   b,resume, unitfarfield,filefarfield, &
   esqx1,esqx2,esqy1,esqy2,esqz1,esqz2, &
   InitialFreq,FinalFreq,FreqStep,phiStart ,phiStop,phiStep,thetaStart,thetaStop ,thetaStep,FileNormalize, &
   SINPML_fullsize,facesNF2FF,NF2FFDecim &
#ifdef CompileWithMPI
   ,MPISubComm,MpiRoot &
#endif
   ,eps00,mu00)
      REAL (KIND=RKIND)           ::  eps00,mu00

      type (nf2ff_t) :: facesNF2FF
      LOGICAL :: NF2FFDecim
      type (limit_t), dimension(1:6), intent(in) :: SINPML_fullsize
      real (kind=Rkind) :: InitialFreq,FinalFreq,FreqStep
      REAL (KIND=RKIND)  ::  thetaStart,thetaStop,thetaStep
      REAL (KIND=RKIND)  ::  phiStart,phiStop,phiStep
      character (LEN=BUFSIZE)  ::   FileNormalize
      integer (KIND=4)    :: unitfarfield
      character (LEN=BUFSIZE)   ::  filefarfield
      integer (kind=4)  :: MPISubComm,MPIRoot

      !---------------------------> inputs <----------------------------------------------------------
      type( bounds_t), intent( IN)  ::  b
      !
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      logical , intent(in)   ::  resume
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in)   ::  &
      sggMiEx(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), &
      sggMiEy(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), &
      sggMiEz(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), &
      sggMiHx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE), &
      sggMiHy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE), &
      sggMiHz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      REAL (KIND=RKIND)  ::tiempo1,tiempo2,field1,field2,dtevol
      integer j,k,field,i,layoutnumber,size,ii,esqx1,esqx2,esqy1,esqy2,esqz1,esqz2,pozi
      character(len=BUFSIZE) :: buFF
      logical :: errnofile,error

!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      cluz=1.0_RKIND/sqrt(eps0*mu0)
      zvac=sqrt(mu0/eps0)
!

      do field=iEx,iHz
        FF%Punto%PhysCoor(field)%x => null()
        FF%Punto%PhysCoor(field)%y => null()
        FF%Punto%PhysCoor(field)%z => null()
      end do
      !!!
      !store absolute limits to later correct edge contributions
      FF%esqx1=max(esqx1,SINPML_fullsize(iHx)%XI)
      FF%esqx2=min(esqx2,SINPML_fullsize(iHx)%XE)
      FF%esqy1=max(esqy1,SINPML_fullsize(iHy)%YI)
      FF%esqy2=min(esqy2,SINPML_fullsize(iHy)%YE)
      FF%esqz1=max(esqz1,SINPML_fullsize(iHz)%ZI)
      FF%esqz2=min(esqz2,SINPML_fullsize(iHz)%ZE)
      !!!!!!!!
      FF%unitfarfield =    unitfarfield
      FF%filefarfield =    filefarfield
      FF%InitialFreq  =    InitialFreq
      FF%FinalFreq    =    FinalFreq
      FF%FreqStep     =    FreqStep
      FF%phiStart     =    phiStart
      FF%phiStop      =    phiStop
      FF%phiStep      =    phiStep
      FF%thetaStart   =    thetaStart
      FF%thetaStop    =    thetaStop
      FF%thetaStep    =    thetaStep
      FF%FileNormalize   = FileNormalize
#ifdef CompileWithMPI
      FF%MPISubComm  =  MPISubComm
      FF%MPIRoot  =  MPIRoot
#endif
      !!!decimacion 30/01/15
      if (NF2FFDecim) then
         FF%NDecim = max(int((0.5_RKIND/FinalFreq)/sgg%dt) - 1,1) ! el -1 es para curarme en salud
      else
         FF%NDecim = 1
      endif
      FF%dtDecim=FF%NDecim * sgg%dt

      !!!

      do field=iEx,iHz
         allocate (FF%Punto%PhysCoor(field)%x(SINPML_fullsize(field)%XI-1 : SINPML_fullsize(field)%XE+1), &
         FF%Punto%PhysCoor(field)%y(SINPML_fullsize(field)%YI-1 : SINPML_fullsize(field)%YE+1), &
         FF%Punto%PhysCoor(field)%z(SINPML_fullsize(field)%ZI-1 : SINPML_fullsize(field)%ZE+1))
      end do

      field=iEx
      do i=SINPML_fullsize(field)%XI-1,SINPML_fullsize(field)%XE+1
         FF%Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=SINPML_fullsize(field)%YI-1,SINPML_fullsize(field)%YE+1
         FF%Punto%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=SINPML_fullsize(field)%ZI-1,SINPML_fullsize(field)%ZE+1
         FF%Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=iEy
      do i=SINPML_fullsize(field)%XI-1,SINPML_fullsize(field)%XE+1
         FF%Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=SINPML_fullsize(field)%YI-1,SINPML_fullsize(field)%YE+1
         FF%Punto%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=SINPML_fullsize(field)%ZI-1,SINPML_fullsize(field)%ZE+1
         FF%Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      field=iEz
      do i=SINPML_fullsize(field)%XI-1,SINPML_fullsize(field)%XE+1
         FF%Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=SINPML_fullsize(field)%YI-1,SINPML_fullsize(field)%YE+1
         FF%Punto%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=SINPML_fullsize(field)%ZI-1,SINPML_fullsize(field)%ZE+1
         FF%Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHx
      do i=SINPML_fullsize(field)%XI-1,SINPML_fullsize(field)%XE+1
         FF%Punto%PhysCoor(field)%x(i)=sgg%LineX(i)
      end do
      do j=SINPML_fullsize(field)%YI-1,SINPML_fullsize(field)%YE+1
         FF%Punto%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=SINPML_fullsize(field)%ZI-1,SINPML_fullsize(field)%ZE+1
         FF%Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHy
      do i=SINPML_fullsize(field)%XI-1,SINPML_fullsize(field)%XE+1
         FF%Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=SINPML_fullsize(field)%YI-1,SINPML_fullsize(field)%YE+1
         FF%Punto%PhysCoor(field)%y(j)=sgg%LineY(j)
      end do
      do k=SINPML_fullsize(field)%ZI-1,SINPML_fullsize(field)%ZE+1
         FF%Punto%PhysCoor(field)%z(k)=(sgg%LineZ(k)+sgg%LineZ(k+1))*0.5_RKIND
      end do
      field=iHz
      do i=SINPML_fullsize(field)%XI-1,SINPML_fullsize(field)%XE+1
         FF%Punto%PhysCoor(field)%x(i)=(sgg%LineX(i)+sgg%LineX(i+1))*0.5_RKIND
      end do
      do j=SINPML_fullsize(field)%YI-1,SINPML_fullsize(field)%YE+1
         FF%Punto%PhysCoor(field)%y(j)=(sgg%LineY(j)+sgg%LineY(j+1))*0.5_RKIND
      end do
      do k=SINPML_fullsize(field)%ZI-1,SINPML_fullsize(field)%ZE+1
         FF%Punto%PhysCoor(field)%z(k)=sgg%LineZ(k)
      end do
      !

      FF%farfieldTr_clonePEC_Front=.false.
      FF%farfieldTr_clonePMC_Front=.false.
      FF%farfieldFr_clonePEC_Back=.false.
      FF%farfieldFr_clonePMC_Back=.false.
      FF%farfieldIz_clonePEC_Right=.false.
      FF%farfieldIz_clonePMC_Right=.false.
      FF%farfieldDe_clonePEC_Left=.false.
      FF%farfieldDe_clonePMC_Left=.false.
      FF%farfieldAb_clonePEC_Up=.false.
      FF%farfieldAb_clonePMC_Up=.false.
      FF%farfieldAr_clonePEC_Down=.false.
      FF%farfieldAr_clonePMC_Down=.false.
      !
      FF%farfieldTr=.false.
      FF%farfieldFr=.false.
      FF%farfieldIz=.false.
      FF%farfieldDe=.false.
      FF%farfieldAr=.false.
      FF%farfieldAb=.false.
      !MPI NO DUPLICAR CALCULOS !revisar cuando se haga lo de las geometrias
      if ((FF%esqx1 > sgg%SINPMLSweep(iHx)%XI).and.(FF%esqx1 <= sgg%SINPMLSweep(IHx)%XE)) &
      FF%farfieldTr=.true.
      if ((FF%esqx2 < sgg%SINPMLSweep(IHx)%XE).and.(FF%esqx2 >= sgg%SINPMLSweep(iHx)%XI)) &
      FF%farfieldFr=.true.
      if ((FF%esqy1 > sgg%SINPMLSweep(iHy)%YI).and.(FF%esqy1 <= sgg%SINPMLSweep(IHy)%YE)) &
      FF%farfieldIz=.true.
      if ((FF%esqy2 < sgg%SINPMLSweep(IHy)%YE).and.(FF%esqy2 >= sgg%SINPMLSweep(iHy)%YI)) &
      FF%farfieldDe=.true.
      if ((FF%esqz1 > sgg%SINPMLSweep(iHz)%ZI).and.(FF%esqz1 <= sgg%SINPMLSweep(IHz)%ZE)) &
      FF%farfieldAb=.true.
      if ((FF%esqz2 < sgg%SINPMLSweep(IHz)%ZE).and.(FF%esqz2 >= sgg%SINPMLSweep(iHz)%ZI)) &
      FF%farfieldAr=.true.
      !

      FF%XDobleAncho= 2*( FF%Punto%PhysCoor(iHx)%x(FF%esqx2)-FF%Punto%PhysCoor(iHx)%x(FF%esqx1) )
      FF%YDobleAncho= 2*( FF%Punto%PhysCoor(iHy)%y(FF%esqy2)-FF%Punto%PhysCoor(iHy)%y(FF%esqy1) )
      FF%ZDobleAncho= 2*( FF%Punto%PhysCoor(iHz)%z(FF%esqz2)-FF%Punto%PhysCoor(iHz)%z(FF%esqz1) )
      FF%XOffsetMinus=2*( FF%Punto%PhysCoor(iHx)%x(FF%esqx1) )
      FF%YOffsetMinus=2*( FF%Punto%PhysCoor(iHy)%y(FF%esqy1) )
      FF%ZOffsetMinus=2*( FF%Punto%PhysCoor(iHz)%z(FF%esqz1) )
      FF%XOffsetPlus= 2*( FF%Punto%PhysCoor(iHx)%x(FF%esqx2))
      FF%YOffsetPlus= 2*( FF%Punto%PhysCoor(iHy)%y(FF%esqy2))
      FF%ZOffsetPlus= 2*( FF%Punto%PhysCoor(iHz)%z(FF%esqz2))

      !manejo de simetrias PEC y PMC
      if (FF%esqx1 <= SINPML_fullsize(IHx)%XI) then
         FF%esqx1  = SINPML_fullsize(IHx)%XI
         if (FF%farfieldTr) then
            error=.TRUE.
         else
            if ((FF%farfieldFr).and.(sgg%Border%IsBackPEC)) FF%farfieldFr_clonePEC_Back=.true.
            if ((FF%farfieldFr).and.(sgg%Border%IsBackPMC)) FF%farfieldFr_clonePMC_Back=.true.
            !
            if ((FF%farfieldIz).and.(sgg%Border%IsBackPEC)) FF%farfieldIz_clonePEC_Back=.true.
            if ((FF%farfieldIz).and.(sgg%Border%IsBackPMC)) FF%farfieldIz_clonePMC_Back=.true.
            !
            if ((FF%farfieldDe).and.(sgg%Border%IsBackPEC)) FF%farfieldDe_clonePEC_Back=.true.
            if ((FF%farfieldDe).and.(sgg%Border%IsBackPMC)) FF%farfieldDe_clonePMC_Back=.true.
            !
            if ((FF%farfieldAr).and.(sgg%Border%IsBackPEC)) FF%farfieldAr_clonePEC_Back=.true.
            if ((FF%farfieldAr).and.(sgg%Border%IsBackPMC)) FF%farfieldAr_clonePMC_Back=.true.
            !
            if ((FF%farfieldAb).and.(sgg%Border%IsBackPEC)) FF%farfieldAb_clonePEC_Back=.true.
            if ((FF%farfieldAb).and.(sgg%Border%IsBackPMC)) FF%farfieldAb_clonePMC_Back=.true.
         endif
      endif
      if (FF%esqx2 >= SINPML_fullsize(IHx)%XE) then
         FF%esqx2  = SINPML_fullsize(IHx)%XE
         if (FF%farfieldFr) then
            error=.TRUE.
         else
            if ((FF%farfieldTr).and.(sgg%Border%IsFrontPEC)) FF%farfieldTr_clonePEC_Front=.true.
            if ((FF%farfieldTr).and.(sgg%Border%IsFrontPMC)) FF%farfieldTr_clonePMC_Front=.true.

            if ((FF%farfieldIz).and.(sgg%Border%IsFrontPEC)) FF%farfieldIz_clonePEC_Front=.true.
            if ((FF%farfieldIz).and.(sgg%Border%IsFrontPMC)) FF%farfieldIz_clonePMC_Front=.true.
            !
            if ((FF%farfieldDe).and.(sgg%Border%IsFrontPEC)) FF%farfieldDe_clonePEC_Front=.true.
            if ((FF%farfieldDe).and.(sgg%Border%IsFrontPMC)) FF%farfieldDe_clonePMC_Front=.true.
            !
            if ((FF%farfieldAr).and.(sgg%Border%IsFrontPEC)) FF%farfieldAr_clonePEC_Front=.true.
            if ((FF%farfieldAr).and.(sgg%Border%IsFrontPMC)) FF%farfieldAr_clonePMC_Front=.true.
            !
            if ((FF%farfieldAb).and.(sgg%Border%IsFrontPEC)) FF%farfieldAb_clonePEC_Front=.true.
            if ((FF%farfieldAb).and.(sgg%Border%IsFrontPMC)) FF%farfieldAb_clonePMC_Front=.true.
         endif
      endif
      if (FF%esqy1 <= SINPML_fullsize(iHy)%YI) then
         FF%esqy1  = SINPML_fullsize(iHy)%YI
         if (FF%farfieldIz) then
            error=.TRUE.
         else
            if ((FF%farfieldFr).and.(sgg%Border%IsLeftPEC)) FF%farfieldFr_clonePEC_Left=.true.
            if ((FF%farfieldFr).and.(sgg%Border%IsLeftPMC)) FF%farfieldFr_clonePMC_Left=.true.

            if ((FF%farfieldTr).and.(sgg%Border%IsLeftPEC)) FF%farfieldTr_clonePEC_Left=.true.
            if ((FF%farfieldTr).and.(sgg%Border%IsLeftPMC)) FF%farfieldTr_clonePMC_Left=.true.
            !
            if ((FF%farfieldDe).and.(sgg%Border%IsLeftPEC)) FF%farfieldDe_clonePEC_Left=.true.
            if ((FF%farfieldDe).and.(sgg%Border%IsLeftPMC)) FF%farfieldDe_clonePMC_Left=.true.
            !
            if ((FF%farfieldAr).and.(sgg%Border%IsLeftPEC)) FF%farfieldAr_clonePEC_Left=.true.
            if ((FF%farfieldAr).and.(sgg%Border%IsLeftPMC)) FF%farfieldAr_clonePMC_Left=.true.
            !
            if ((FF%farfieldAb).and.(sgg%Border%IsLeftPEC)) FF%farfieldAb_clonePEC_Left=.true.
            if ((FF%farfieldAb).and.(sgg%Border%IsLeftPMC)) FF%farfieldAb_clonePMC_Left=.true.
         endif
      endif
      if (FF%esqy2 >= SINPML_fullsize(IHy)%YE) then
         FF%esqy2  = SINPML_fullsize(IHy)%YE
         if (FF%farfieldDe) then
            error=.TRUE.
         else
            if ((FF%farfieldFr).and.(sgg%Border%IsRightPEC)) FF%farfieldFr_clonePEC_Right=.true.
            if ((FF%farfieldFr).and.(sgg%Border%IsRightPMC)) FF%farfieldFr_clonePMC_Right=.true.

            if ((FF%farfieldTr).and.(sgg%Border%IsRightPEC)) FF%farfieldTr_clonePEC_Right=.true.
            if ((FF%farfieldTr).and.(sgg%Border%IsRightPMC)) FF%farfieldTr_clonePMC_Right=.true.
            !
            if ((FF%farfieldIz).and.(sgg%Border%IsRightPEC)) FF%farfieldIz_clonePEC_Right=.true.
            if ((FF%farfieldIz).and.(sgg%Border%IsRightPMC)) FF%farfieldIz_clonePMC_Right=.true.
            !
            if ((FF%farfieldAr).and.(sgg%Border%IsRightPEC)) FF%farfieldAr_clonePEC_Right=.true.
            if ((FF%farfieldAr).and.(sgg%Border%IsRightPMC)) FF%farfieldAr_clonePMC_Right=.true.
            !
            if ((FF%farfieldAb).and.(sgg%Border%IsRightPEC)) FF%farfieldAb_clonePEC_Right=.true.
            if ((FF%farfieldAb).and.(sgg%Border%IsRightPMC)) FF%farfieldAb_clonePMC_Right=.true.
         endif
      endif
      if (FF%esqz1 <= SINPML_fullsize(iHz)%ZI) then
         FF%esqz1  = SINPML_fullsize(iHz)%ZI
         if (FF%farfieldAb) then
            error=.TRUE.
         else
            if ((FF%farfieldFr).and.(sgg%Border%IsDownPEC)) FF%farfieldFr_clonePEC_Down=.true.
            if ((FF%farfieldFr).and.(sgg%Border%IsDownPMC)) FF%farfieldFr_clonePMC_Down=.true.

            if ((FF%farfieldTr).and.(sgg%Border%IsDownPEC)) FF%farfieldTr_clonePEC_Down=.true.
            if ((FF%farfieldTr).and.(sgg%Border%IsDownPMC)) FF%farfieldTr_clonePMC_Down=.true.
            !
            if ((FF%farfieldIz).and.(sgg%Border%IsDownPEC)) FF%farfieldIz_clonePEC_Down=.true.
            if ((FF%farfieldIz).and.(sgg%Border%IsDownPMC)) FF%farfieldIz_clonePMC_Down=.true.
            !
            if ((FF%farfieldDe).and.(sgg%Border%IsDownPEC)) FF%farfieldDe_clonePEC_Down=.true.
            if ((FF%farfieldDe).and.(sgg%Border%IsDownPMC)) FF%farfieldDe_clonePMC_Down=.true.
            !
            if ((FF%farfieldAr).and.(sgg%Border%IsDownPEC)) FF%farfieldAr_clonePEC_Down=.true.
            if ((FF%farfieldAr).and.(sgg%Border%IsDownPMC)) FF%farfieldAr_clonePMC_Down=.true.

         endif
      endif
      if (FF%esqz2 >= SINPML_fullsize(IHz)%ZE) then
         FF%esqz2  = SINPML_fullsize(IHz)%ZE
         if (FF%farfieldAr) then
            error=.TRUE.
         else
            if ((FF%farfieldFr).and.(sgg%Border%IsUpPEC)) FF%farfieldFr_clonePEC_Up=.true.
            if ((FF%farfieldFr).and.(sgg%Border%IsUpPMC)) FF%farfieldFr_clonePMC_Up=.true.

            if ((FF%farfieldTr).and.(sgg%Border%IsUpPEC)) FF%farfieldTr_clonePEC_Up=.true.
            if ((FF%farfieldTr).and.(sgg%Border%IsUpPMC)) FF%farfieldTr_clonePMC_Up=.true.
            !
            if ((FF%farfieldIz).and.(sgg%Border%IsUpPEC)) FF%farfieldIz_clonePEC_Up=.true.
            if ((FF%farfieldIz).and.(sgg%Border%IsUpPMC)) FF%farfieldIz_clonePMC_Up=.true.
            !
            if ((FF%farfieldDe).and.(sgg%Border%IsUpPEC)) FF%farfieldDe_clonePEC_Up=.true.
            if ((FF%farfieldDe).and.(sgg%Border%IsUpPMC)) FF%farfieldDe_clonePMC_Up=.true.
            !
            if ((FF%farfieldAb).and.(sgg%Border%IsUpPEC)) FF%farfieldAb_clonePEC_Up=.true.
            if ((FF%farfieldAb).and.(sgg%Border%IsUpPMC)) FF%farfieldAb_clonePMC_Up=.true.
         endif
      endif
      !
      FF%farfieldFr_clonePEC_Back =((.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePMC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePMC_Back)) ) .and. FF%farfieldFr_clonePEC_Back
      FF%farfieldFr_clonePMC_Back =((.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePMC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePMC_Back)) ) .and. FF%farfieldFr_clonePMC_Back
      FF%farfieldTr_clonePEC_Front=((.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePMC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePMC_Back)) ) .and. FF%farfieldTr_clonePEC_Front
      FF%farfieldTr_clonePMC_Front=((.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePEC_Front.and.FF%farfieldFr_clonePMC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePEC_Back)).or. &
      (.not.(FF%farfieldTr_clonePMC_Front.and.FF%farfieldFr_clonePMC_Back)) ) .and. FF%farfieldTr_clonePMC_Front
      !
      FF%farfieldDe_clonePEC_Left =((.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePMC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePMC_Left)) ) .and. FF%farfieldDe_clonePEC_Left
      FF%farfieldDe_clonePMC_Left =((.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePMC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePMC_Left)) ) .and. FF%farfieldDe_clonePMC_Left
      FF%farfieldIz_clonePEC_Right=((.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePMC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePMC_Left)) ) .and. FF%farfieldIz_clonePEC_Right
      FF%farfieldIz_clonePMC_Right=((.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePEC_Right.and.FF%farfieldDe_clonePMC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePEC_Left)).or. &
      (.not.(FF%farfieldIz_clonePMC_Right.and.FF%farfieldDe_clonePMC_Left)) ) .and. FF%farfieldIz_clonePMC_Right
      !
      FF%farfieldAb_clonePEC_Up  = ((.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePMC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePMC_Down)) ) .and. FF%farfieldAb_clonePEC_Up
      FF%farfieldAb_clonePMC_Up  = ((.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePMC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePMC_Down)) ) .and. FF%farfieldAb_clonePMC_Up
      FF%farfieldAr_clonePEC_Down= ((.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePMC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePMC_Down)) ) .and. FF%farfieldAr_clonePEC_Down
      FF%farfieldAr_clonePMC_Down= ((.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePEC_Up.and.FF%farfieldAr_clonePMC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePEC_Down)).or. &
      (.not.(FF%farfieldAb_clonePMC_Up.and.FF%farfieldAr_clonePMC_Down)) ) .and. FF%farfieldAr_clonePMC_Down


      !
      error=FF%farfieldFr_clonePEC_Back.OR.FF%farfieldFr_clonePMC_Back.OR. &
      FF%farfieldTr_clonePEC_Front.OR.FF%farfieldTr_clonePMC_Front.OR. &
      FF%farfieldDe_clonePEC_Left.OR.FF%farfieldDe_clonePMC_Left.OR. &
      FF%farfieldIz_clonePEC_Right.OR.FF%farfieldIz_clonePMC_Right.OR. &
      FF%farfieldAb_clonePEC_Up.OR.FF%farfieldAb_clonePMC_Up.OR. &
      FF%farfieldAr_clonePEC_Down.OR.FF%farfieldAr_clonePMC_Down
      if (error) then
         buff='NF2FF STILL UNSUPPORTED'
         call stoponerror(layoutnumber,size,Buff)
      endif
      !!!fin casuistica de clonados periodicos



      !find the coordinate limits of the Huygens Box for each component
      FF%TrFr%I%tra%Ez=Max( sgg%SINPMLSweep(iEz)%XI,       FF%esqx1     )
      FF%TrFr%I%fro%Ez=Min( sgg%SINPMLSweep(iEz)%XE,       FF%esqx2     )
      FF%TrFr%J%com%Ez=Max( sgg%SINPMLSweep(iEz)%YI,       FF%esqy1     )
      FF%TrFr%J%fin%Ez=Min( sgg%SINPMLSweep(iEz)%YE,       FF%esqy2     )
      FF%TrFr%K%com%Ez=Max( sgg%SINPMLSweep(iEz)%ZI,       FF%esqz1     )
      FF%TrFr%K%fin%Ez=MIn( sgg%SINPMLSweep(iEz)%ZE,       FF%esqz2-1   )
      !
      FF%TrFr%I%tra%Ey=Max( sgg%SINPMLSweep(iEy)%XI,       FF%esqx1     )
      FF%TrFr%I%fro%Ey=Min( sgg%SINPMLSweep(iEy)%XE,       FF%esqx2     )
      FF%TrFr%J%com%Ey=Max( sgg%SINPMLSweep(iEy)%YI,       FF%esqy1     )
      FF%TrFr%J%fin%Ey=Min( sgg%SINPMLSweep(iEy)%YE,       FF%esqy2-1   )
      FF%TrFr%K%com%Ey=Max( sgg%SINPMLSweep(iEy)%ZI,       FF%esqz1     )
      FF%TrFr%K%fin%Ey=MIn( sgg%SINPMLSweep(iEy)%ZE-01,    FF%esqz2     ) !MPI NO DUPLICAR CALCULOS
      !
      FF%TrFr%I%tra%Hy= FF%TrFr%I%tra%Ez - 1
      FF%TrFr%I%fro%Hy= FF%TrFr%I%fro%Ez
      FF%TrFr%J%com%Hy= FF%TrFr%J%com%Ez
      FF%TrFr%J%fin%Hy= FF%TrFr%J%fin%Ez
      FF%TrFr%K%com%Hy= FF%TrFr%K%com%Ez
      FF%TrFr%K%fin%Hy= FF%TrFr%K%fin%Ez
      !
      FF%TrFr%I%tra%Hz= FF%TrFr%I%tra%Ey -1
      FF%TrFr%I%fro%Hz= FF%TrFr%I%fro%Ey
      FF%TrFr%J%com%Hz= FF%TrFr%J%com%Ey
      FF%TrFr%J%fin%Hz= FF%TrFr%J%fin%Ey
      FF%TrFr%K%com%Hz= FF%TrFr%K%com%Ey
      FF%TrFr%K%fin%Hz= FF%TrFr%K%fin%Ey
      !
      !
      FF%IzDe%J%izq%Ex=Max( sgg%SINPMLSweep(iEx)%yI,       FF%esqy1     )
      FF%IzDe%J%der%Ex=Min( sgg%SINPMLSweep(iEx)%yE,       FF%esqy2     )
      FF%IzDe%I%com%Ex=Max( sgg%SINPMLSweep(iEx)%xI,       FF%esqx1     )
      FF%IzDe%I%fin%Ex=Min( sgg%SINPMLSweep(iEx)%xE,       FF%esqx2-1   )
      FF%IzDe%K%com%Ex=Max( sgg%SINPMLSweep(iEx)%ZI,       FF%esqz1     )
      FF%IzDe%K%fin%Ex=MIn( sgg%SINPMLSweep(iEx)%ZE-01,    FF%esqz2     ) !MPI NO DUPLICAR CALCULOS
      !
      FF%IzDe%J%izq%Ez=Max( sgg%SINPMLSweep(iEz)%yI,    FF%esqy1     )
      FF%IzDe%J%der%Ez=Min( sgg%SINPMLSweep(iEz)%yE,    FF%esqy2     )
      FF%IzDe%I%com%Ez=Max( sgg%SINPMLSweep(iEz)%xI,    FF%esqx1     )
      FF%IzDe%I%fin%Ez=Min( sgg%SINPMLSweep(iEz)%xE,    FF%esqx2     )
      FF%IzDe%K%com%Ez=Max( sgg%SINPMLSweep(iEz)%ZI,    FF%esqz1     )
      FF%IzDe%K%fin%Ez=MIn( sgg%SINPMLSweep(iEz)%ZE,    FF%esqz2-1   )
      !
      FF%IzDe%J%izq%Hz= FF%IzDe%J%izq%Ex - 1
      FF%IzDe%J%der%Hz= FF%IzDe%J%der%Ex
      FF%IzDe%I%com%Hz= FF%IzDe%I%com%Ex
      FF%IzDe%I%fin%Hz= FF%IzDe%I%fin%Ex
      FF%IzDe%K%com%Hz= FF%IzDe%K%com%Ex
      FF%IzDe%K%fin%Hz= FF%IzDe%K%fin%Ex
      !
      FF%IzDe%J%izq%Hx= FF%IzDe%J%izq%Ez - 1
      FF%IzDe%J%der%Hx= FF%IzDe%J%der%Ez
      FF%IzDe%I%com%Hx= FF%IzDe%I%com%Ez
      FF%IzDe%I%fin%Hx= FF%IzDe%I%fin%Ez
      FF%IzDe%K%com%Hx= FF%IzDe%K%com%Ez
      FF%IzDe%K%fin%Hx= FF%IzDe%K%fin%Ez
      !
      !
      FF%AbAr%K%aba%Ey=Max( sgg%SINPMLSweep(iEy)%ZI,    FF%esqz1     )
      FF%AbAr%K%arr%Ey=Min( sgg%SINPMLSweep(iEy)%ZE,    FF%esqz2     )
      FF%AbAr%I%com%Ey=Max( sgg%SINPMLSweep(iEy)%XI,    FF%esqx1     )
      FF%AbAr%I%fin%Ey=Min( sgg%SINPMLSweep(iEy)%XE,    FF%esqx2     )
      FF%AbAr%J%com%Ey=Max( sgg%SINPMLSweep(iEy)%YI,    FF%esqy1     )
      FF%AbAr%J%fin%Ey=Min( sgg%SINPMLSweep(iEy)%YE,    FF%esqy2-1   )
      !
      FF%AbAr%K%aba%Ex=Max( sgg%SINPMLSweep(iEx)%ZI,    FF%esqz1     )
      FF%AbAr%K%arr%Ex=Min( sgg%SINPMLSweep(iEx)%ZE,    FF%esqz2     )
      FF%AbAr%I%com%Ex=Max( sgg%SINPMLSweep(iEx)%XI,    FF%esqx1     )
      FF%AbAr%I%fin%Ex=Min( sgg%SINPMLSweep(iEx)%XE,    FF%esqx2-1   )
      FF%AbAr%J%com%Ex=Max( sgg%SINPMLSweep(iEx)%YI,    FF%esqy1     )
      FF%AbAr%J%fin%Ex=Min( sgg%SINPMLSweep(iEx)%YE,    FF%esqy2     )
      !
      FF%AbAr%K%aba%Hx= FF%AbAr%K%aba%Ey - 1
      FF%AbAr%K%arr%Hx= FF%AbAr%K%arr%Ey
      FF%AbAr%I%com%Hx= FF%AbAr%I%com%Ey
      FF%AbAr%I%fin%Hx= FF%AbAr%I%fin%Ey
      FF%AbAr%J%com%Hx= FF%AbAr%J%com%Ey
      FF%AbAr%J%fin%Hx= FF%AbAr%J%fin%Ey
      !
      FF%AbAr%K%aba%Hy= FF%AbAr%K%aba%Ex - 1
      FF%AbAr%K%arr%Hy= FF%AbAr%K%arr%Ex
      FF%AbAr%I%com%Hy= FF%AbAr%I%com%Ex
      FF%AbAr%I%fin%Hy= FF%AbAr%I%fin%Ex
      FF%AbAr%J%com%Hy= FF%AbAr%J%com%Ex
      FF%AbAr%J%fin%Hy= FF%AbAr%J%fin%Ex


      !check if materials are crossed by the box
      If( FF%farfieldTr) then
         !Ez Back
         i = FF%TrFr%I%tra%Ez !Back
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               if (sggMiEz( i, j, k) /=1 ) then
                  write (buff,'(a,3i7)') 'Back NF/FF region intersects a material'//' at Ez',i,j,k
                  if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Ey Back
         i = FF%TrFr%I%tra%Ey
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               if (sggMiEy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Back NF/FF region intersects a material'//' at Ey',i,j,k
                  if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            End do
         end do
      endif
      !--->
      If( FF%farfieldFr) then
         !Ez  Front
         i = FF%TrFr%I%fro%Ez !Front
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               if (sggMiEz(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Front NF/FF region intersects a material'//' at Ez',i,j,k
                  if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Ey  Front
         i = FF%TrFr%I%fro%Ey !Front
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               if (sggMiEy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Front NF/FF region intersects a material'//' at Ey',i,j,k
                  if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldIz) then
         !Ex Left
         j = FF%IzDe%J%izq%Ex  !Left
         do k = FF%IzDe%K%com%Ex, FF%IzDe%K%fin%Ex
            do i = FF%IzDe%I%com%Ex, FF%IzDe%I%fin%Ex
               if (sggMiEx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Left NF/FF region intersects a material'//' at Ex',i,j,k
                  if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Ez Left
         j = FF%IzDe%J%izq%Ez  !Left
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               if (sggMiEz(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Left NF/FF region intersects a material'//' at Ez',i,j,k
                  if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldDe) then
         !Ez  Right
         j = FF%IzDe%J%der%Ez !Right
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               if (sggMiEz(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Right NF/FF region intersects a material'//' at Ez',i,j,k
                  if (((sggMiEz(i,j,k) ==0).or.(sgg%med(sggMiEz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Ex  Right
         j = FF%IzDe%J%der%Ex !Right
         do k = FF%IzDe%K%com%Ex,FF%IzDe%K%fin%Ex
            do i=FF%IzDe%I%com%Ex,FF%IzDe%I%fin%Ex
               if (sggMiEx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Right NF/FF region intersects a material'//' at Ex',i,j,k
                  if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            End do
         end do
      endif
      !--->
      If( FF%farfieldAb) then
         !Ex  Down
         k = FF%AbAr%K%aba%Ex  !Down
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            Do i=FF%AbAr%I%com%Ex,FF%AbAr%I%fin%Ex
               if (sggMiEx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Down NF/FF region intersects a material'//' at Ex',i,j,k
                  if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Ey Down
         k = FF%AbAr%K%aba%Ey  !Down
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               if (sggMiEy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Down NF/FF region intersects a material'//' at Ey',i,j,k
                  if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldAr) then
         !Ex Up
         k = FF%AbAr%K%arr%Ex  !Up
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            do i = FF%AbAr%I%com%Ex, FF%AbAr%I%fin%Ex
               if (sggMiEx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Up NF/FF region intersects a material'//' at Ex',i,j,k
                  if (((sggMiEx(i,j,k) ==0).or.(sgg%med(sggMiEx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Ey Up
         k = FF%AbAr%K%arr%Ey
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               if (sggMiEy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Up NF/FF region intersects a material'//' at Ey',i,j,k
                  if (((sggMiEy(i,j,k) ==0).or.(sgg%med(sggMiEy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      !!!
      if( FF%farfieldTr) then
         !Hz Back
         i = FF%TrFr%I%tra%Hz  !Back
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               if (sggMiHz(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Back NF/FF region intersects a material'//' at Hz',i,j,k
                  if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Hy Back
         i = FF%TrFr%I%tra%Hy  !Back
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               if (sggMiHy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Back NF/FF region intersects a material'//' at Hy',i,j,k
                  if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      if( FF%farfieldFr) then
         !Hz  Front
         i = FF%TrFr%I%fro%Hz !Front
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               if (sggMiHz(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Front NF/FF region intersects a material'//' at Hz',i,j,k
                  if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Hy  Front
         i = FF%TrFr%I%fro%Hy !Front
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               if (sggMiHy( i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Front NF/FF region intersects a material'//' at Hy',i,j,k
                  if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo

      endif
      if( FF%farfieldIz) then
         !Hx Left
         j = FF%IzDe%J%izq%Hx  !Left
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               if (sggMiHx( i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Left NF/FF region intersects a material'//' at Hx',i,j,k
                  if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Hz Left
         j = FF%IzDe%J%izq%Hz  !Left
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               if (sggMiHz( i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Left NF/FF region intersects a material'//' at Hz',i,j,k
                  if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      if( FF%farfieldDe) then
         !Hx  Right
         j = FF%IzDe%J%der%Hx !Right
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               if (sggMiHx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Right NF/FF region intersects a material'//' at Hx',i,j,k
                  if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Hz  Right
         j = FF%IzDe%J%der%Hz !Right
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               if (sggMiHz(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Right NF/FF region intersects a material'//' at Hz',i,j,k
                  if (((sggMiHz(i,j,k) ==0).or.(sgg%med(sggMiHz(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      if( FF%farfieldAb) then
         !Hx  Down
         k = FF%AbAr%K%aba%Hx  !Down
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               if (sggMiHx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Down NF/FF region intersects a material'//' at Hx',i,j,k
                  if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Hy  Down
         k = FF%AbAr%K%aba%Hy  !Down
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            do i=FF%AbAr%I%com%Hy,FF%AbAr%I%fin%Hy
               if (sggMiHy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Down NF/FF region intersects a material'//' at Hy',i,j,k
                  if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldAr) then
         !Hx Up
         k = FF%AbAr%K%arr%Hx  !Up
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               if (sggMiHx(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Up NF/FF region intersects a material'//' at Hx',i,j,k
                  if (((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
         !Hy Up
         k=FF%AbAr%K%arr%Hy  !Up
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            do i = FF%AbAr%I%com%Hy, FF%AbAr%I%fin%Hy
               if (sggMiHy(i,j,k) /=1 ) then
                  write (buff,'(a,3i7)') 'Up NF/FF region intersects a material'//' at Hy',i,j,k
                  if (((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec)).and. .not. &
                  ((i == sgg%SINPMLSweep(iHx)%XI).or.(j == sgg%SINPMLSweep(iHy)%YI).or.(k == sgg%SINPMLSweep(iHz)%ZI).or. &
                  (i == sgg%SINPMLSweep(iHx)%XE).or.(j == sgg%SINPMLSweep(iHy)%YE).or.(k == sgg%SINPMLSweep(iHz)%ZE))) &
                  call stoponerror(layoutnumber,size,Buff)
               endif
            enddo
         enddo
      endif

      !allocatea espacio

      if (FF%FreqStep/=0) then
         FF%NumFreqs=int(abs(FF%InitialFreq-FF%FinalFreq)/FF%FreqStep)+1
      else
         FF%NumFreqs=1 !default
      endif


      if ((FF%NumFreqs < 0)) then
         write(Buff,*) 'Freq. range for NF/FF invalid',FF%NumFreqs,FF%InitialFreq,FF%FinalFreq,FF%FreqStep
         call stoponerror(layoutnumber,size,Buff)
      endif
      if ((FF%NumFreqs > 100000)) then
         Buff='Too many NF/FF frequencies requested (>100000)'
         call stoponerror(layoutnumber,size,Buff)
      endif

      ALLOCATE (FF%expIwdt(1:FF%NumFreqs),FF%auxExp_E(1:FF%NumFreqs),FF%auxExp_H(1:FF%NumFreqs),FF%dftEntrada(1:FF%NumFreqs))
      !
      if( FF%farfieldIz) then
         allocate (FF%ExIz( 0 :  b%Ex%NX-1, 0 :  b%Ex%NZ-1, 1:FF%NumFreqs))
         allocate (FF%EzIz( 0 :  b%Ez%NX-1, 0 :  b%Ez%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HxIz( 0 :  b%Hx%NX-1, 0 :  b%Hx%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzIz( 0 :  b%Hz%NX-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HxIz2( 0 :  b%Hx%NX-1, 0 :  b%Hx%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzIz2( 0 :  b%Hz%NX-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
      endif
      if( FF%farfieldDe) then
         allocate (FF%ExDe( 0 :  b%Ex%NX-1, 0 :  b%Ex%NZ-1, 1:FF%NumFreqs))
         allocate (FF%EzDe( 0 :  b%Ez%NX-1, 0 :  b%Ez%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HxDe( 0 :  b%Hx%NX-1, 0 :  b%Hx%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzDe( 0 :  b%Hz%NX-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HxDe2( 0 :  b%Hx%NX-1, 0 :  b%Hx%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzDe2( 0 :  b%Hz%NX-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
      endif
      if( FF%farfieldAb) then
         allocate (FF%ExAb( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 1:FF%NumFreqs))
         allocate (FF%EyAb( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 1:FF%NumFreqs))
         allocate (FF%HxAb( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 1:FF%NumFreqs))
         allocate (FF%HyAb( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 1:FF%NumFreqs))
         allocate (FF%HxAb2( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 1:FF%NumFreqs))
         allocate (FF%HyAb2( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 1:FF%NumFreqs))
      endif
      if( FF%farfieldAr) then
         allocate (FF%ExAr( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 1:FF%NumFreqs))
         allocate (FF%EyAr( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 1:FF%NumFreqs))
         allocate (FF%HxAr( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 1:FF%NumFreqs))
         allocate (FF%HyAr( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 1:FF%NumFreqs))
         allocate (FF%HxAr2( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 1:FF%NumFreqs))
         allocate (FF%HyAr2( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 1:FF%NumFreqs))
      endif
      if( FF%farfieldFr) then
         allocate (FF%EyFr( 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1, 1:FF%NumFreqs))
         allocate (FF%EzFr( 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HyFr( 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzFr( 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HyFr2( 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzFr2( 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
      endif
      if( FF%farfieldTr) then
         allocate (FF%EyTr( 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1, 1:FF%NumFreqs))
         allocate (FF%EzTr( 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HyTr( 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzTr( 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HyTr2( 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1, 1:FF%NumFreqs))
         allocate (FF%HzTr2( 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1, 1:FF%NumFreqs))
      endif


      !Read the time evolution file
      !
      errnofile = .FALSE.
      INQUIRE (FILE=trim(adjustl(FF%FileNormalize)), EXIST=errnofile)
      IF ( .NOT. errnofile) THEN
         Buff=trim(adjustl(FF%FileNormalize))//' DOES NOT EXIST'
         CALL STOPONERROR (layoutnumber,size,Buff)
      END IF
      OPEN (15, FILE=trim(adjustl(FF%FileNormalize)))
      READ (15,*) tiempo1, field1
      READ (15,*) tiempo2, field2
      CLOSE (15)
      dtevol = tiempo2 - tiempo1 !!!ojo tocar para permit scaling pq. no estan sampleadas uniformemente 06118
      FF%dftEntrada=0.0_RKIND



      pozi=index(FF%filefarfield,'_log_')
      if (pozi/=0) then
         FF%InitialFreq=log10(FF%InitialFreq)
         FF%FinalFreq=log10(FF%FinalFreq)
         FF%FreqStep=abs(FF%InitialFreq-FF%FinalFreq)/(FF%NumFreqs)
      endif

      if (pozi == 0) then
         !vector con exponenciales con el sgg%dt del fichero de entrada
         do ii=1,FF%NumFreqs
            FF%expIwdt(ii)=     Exp(mcpi2*(FF%InitialFreq + (ii-1) *FF%FreqStep) *dtevol     )
            !iniciales
            FF%auxExp_E(ii)=   dtevol * (1.0E0_RKIND, 0.0E0_RKIND)   !hay que multiplicar por sgg%dt !bug
         END DO
      else !logaritmico
         do ii=1,FF%NumFreqs
            FF%expIwdt(ii)=     Exp(mcpi2*(10.0_RKIND **(FF%InitialFreq + (ii-1) *FF%FreqStep)) *dtevol     )
            !iniciales
            FF%auxExp_E(ii)=   dtevol * (1.0E0_RKIND, 0.0E0_RKIND)   !hay que multiplicar por sgg%dt !bug 
         END DO
      endif


      !read the normalization file and find its DFT


      OPEN (15, FILE=trim(adjustl(FF%FileNormalize)))
      READ (15,*) tiempo1, field1
      DO
         READ (15,*, end=98) tiempo1, field1
         do ii=1,FF%NumFreqs
            FF%dftEntrada(ii) = FF%dftEntrada(ii) + field1 * FF%auxExp_E(ii)
         END DO
         !solo los samples despues de 1 actualizan el valor
         !ver rutina dtft en postprocesws
         do ii=1,FF%NumFreqs
            FF%auxExp_E(ii)=FF%auxExp_E(ii) * FF%expIwdt(ii)
         end do
      END DO
98    CONTINUE
      CLOSE (15)

      !!
      !machaca con el vector con exponenciales para el sgg%dt de la simulacion
      if (pozi == 0) then
         do ii=1,FF%NumFreqs
            FF%expIwdt(ii)= Exp(-(0.0_RKIND,1.0_RKIND)*2.0_RKIND * pi*(FF%InitialFreq + (ii-1) *FF%FreqStep) *FF%dtDecim     ) !en vez uso el dtdecimado 30/01/15
         END DO
      else !logaritmico
         do ii=1,FF%NumFreqs
            FF%expIwdt(ii)= Exp(-(0.0_RKIND,1.0_RKIND)*2.0_RKIND * pi*(10.0_RKIND **(FF%InitialFreq + (ii-1) *FF%FreqStep)) *FF%dtDecim  ) !en vez uso el decimado 30/01/15
         END DO
      endif

      !!
      !machaca con el vector con exponenciales para el sgg%dt de la simulacion
      if (pozi == 0) then
         do ii=1,FF%NumFreqs
            FF%auxExp_E(ii)=   FF%dtDecim   * (1.0E0_RKIND, 0.0E0_RKIND)   !hay que multiplicar por FF%dtDecim   
            FF%auxExp_H(ii)=   FF%dtDecim   * (1.0E0_RKIND, 0.0E0_RKIND) * &
            Exp(-(0.0_RKIND,1.0_RKIND)*2.0_RKIND * pi*(FF%InitialFreq + (ii-1) *FF%FreqStep) *(sgg%dt * 0.5_RKIND) )
         END DO
      else !logaritmico
         do ii=1,FF%NumFreqs
            FF%auxExp_E(ii)=   FF%dtDecim   * (1.0E0_RKIND, 0.0E0_RKIND)   !hay que multiplicar por FF%dtDecim 
            FF%auxExp_H(ii)=   FF%dtDecim   * (1.0E0_RKIND, 0.0E0_RKIND) * &
            Exp(-(0.0_RKIND,1.0_RKIND)*2.0_RKIND * pi*(10.0_RKIND **(FF%InitialFreq + (ii-1) *FF%FreqStep)) *(sgg%dt * 0.5_RKIND) )
         END DO
      endif

      if (.not.resume) then
         !
         if( FF%farfieldIz) then
            FF%ExIz=0.0_RKIND
            FF%EzIz=0.0_RKIND
            FF%HxIz=0.0_RKIND
            FF%HzIz=0.0_RKIND
            FF%HxIz2=0.0_RKIND
            FF%HzIz2=0.0_RKIND
         endif
         if( FF%farfieldDe) then
            FF%ExDe=0.0_RKIND
            FF%EzDe=0.0_RKIND
            FF%HxDe=0.0_RKIND
            FF%HzDe=0.0_RKIND
            FF%HxDe2=0.0_RKIND
            FF%HzDe2=0.0_RKIND
         endif
         if( FF%farfieldAb) then
            FF%ExAb=0.0_RKIND
            FF%EyAb=0.0_RKIND
            FF%HxAb=0.0_RKIND
            FF%HyAb=0.0_RKIND
            FF%HxAb2=0.0_RKIND
            FF%HyAb2=0.0_RKIND
         endif
         if( FF%farfieldAr) then
            FF%ExAr=0.0_RKIND
            FF%EyAr=0.0_RKIND
            FF%HxAr=0.0_RKIND
            FF%HyAr=0.0_RKIND
            FF%HxAr2=0.0_RKIND
            FF%HyAr2=0.0_RKIND
         endif
         if( FF%farfieldFr) then
            FF%EyFr=0.0_RKIND
            FF%EzFr=0.0_RKIND
            FF%HyFr=0.0_RKIND
            FF%HzFr=0.0_RKIND
            FF%HyFr2=0.0_RKIND
            FF%HzFr2=0.0_RKIND
         endif
         if( FF%farfieldTr) then
            FF%EyTr=0.0_RKIND
            FF%EzTr=0.0_RKIND
            FF%HyTr=0.0_RKIND
            FF%HzTr=0.0_RKIND
            FF%HyTr2=0.0_RKIND
            FF%HzTr2=0.0_RKIND
         endif
      else
         call readFarfield(b)
      endif


      RETURN
   end subroutine InitFarField



   !**************************************************************************************************
   subroutine UpdateFarField(ntime, b, Ex, Ey, Ez,Hx,Hy,Hz)
      type( bounds_t), intent( IN)  ::  b
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer  ::  i, j, k, i_m, j_m, k_m,ii, ntime
      !---------------------------> empieza UpdateFarField <---------------------------------------


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (mod(ntime-1,FF%Ndecim) /=0) return !decimacion 30/01%15
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!
      !electricos
      !!!!!!!!!!!!!!!!!!!!
      If( FF%farfieldTr) then
         !Ez Back
         i = FF%TrFr%I%tra%Ez !Back
         i_m = i - b%Ez%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            k_m = k - b%Ez%ZI
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               j_m = j - b%Ez%YI
               do ii=1,FF%NumFreqs
                  FF%EzTr( j_m, k_m,ii) = FF%EzTr(j_m, k_m,ii) + FF%auxExp_E(ii) * Ez( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Ey Back
         i = FF%TrFr%I%tra%Ey  !Back
         i_m = i - b%Ey%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            k_m = k - b%Ey%ZI
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               j_m = j - b%Ey%YI
               do ii=1,FF%NumFreqs
                  FF%EyTr( j_m, k_m,ii) = FF%EyTr(j_m, k_m,ii) + FF%auxExp_E(ii) *  Ey( i_m, j_m, k_m)

               end do
            End do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !--->
      If( FF%farfieldFr) then
         !Ez  Front
         i = FF%TrFr%I%fro%Ez !Front
         i_m = i - b%Ez%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            k_m = k - b%Ez%ZI
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               j_m = j - b%Ez%YI
               do ii=1,FF%NumFreqs
                  FF%EzFr(j_m, k_m,ii) = FF%EzFr( j_m, k_m,ii) + FF%auxExp_E(ii) *  Ez( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Ey  Front
         i = FF%TrFr%I%fro%Ey !Front
         i_m = i - b%Ey%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            k_m = k - b%Ey%ZI
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               j_m = j - b%Ey%YI
               do ii=1,FF%NumFreqs
                  FF%EyFr( j_m, k_m,ii) = FF%EyFr(j_m, k_m,ii) + FF%auxExp_E(ii) *  Ey( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !--->
      If( FF%farfieldIz) then
         !Ex Left
         j = FF%IzDe%J%izq%Ex  !Left
         j_m = j - b%Ex%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Ex, FF%IzDe%K%fin%Ex
            k_m = k - b%Ex%ZI
            do i = FF%IzDe%I%com%Ex, FF%IzDe%I%fin%Ex
               i_m = i - b%Ex%XI
               do ii=1,FF%NumFreqs
                  FF%ExIz( i_m, k_m,ii) = FF%ExIz( i_m, k_m,ii) + FF%auxExp_E(ii) *  Ex( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Ez Left
         j = FF%IzDe%J%izq%Ez  !Left
         j_m = j - b%Ez%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            k_m = k - b%Ez%ZI
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               i_m = i - b%Ez%XI
               do ii=1,FF%NumFreqs
                  FF%EzIz( i_m, k_m,ii) = FF%EzIz( i_m, k_m,ii) + FF%auxExp_E(ii) *  Ez( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !--->
      If( FF%farfieldDe) then
         !Ez  Right
         j = FF%IzDe%J%der%Ez !Right
         j_m = j - b%Ez%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            k_m = k - b%Ez%ZI
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               i_m = i - b%Ez%XI
               do ii=1,FF%NumFreqs
                  FF%EzDe( i_m, k_m,ii) = FF%EzDe( i_m, k_m,ii) + FF%auxExp_E(ii) *  Ez( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Ex  Right
         j = FF%IzDe%J%der%Ex !Right
         j_m = j - b%Ex%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Ex,FF%IzDe%K%fin%Ex
            k_m = k - b%Ex%ZI
            do i=FF%IzDe%I%com%Ex,FF%IzDe%I%fin%Ex
               i_m = i - b%Ex%XI
               do ii=1,FF%NumFreqs
                  FF%ExDe( i_m, k_m,ii) = FF%ExDe( i_m, k_m,ii) + FF%auxExp_E(ii) *  Ex( i_m, j_m, k_m)

               end do
            End do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !--->
      If( FF%farfieldAb) then
         !Ex  Down
         k = FF%AbAr%K%aba%Ex  !Down
         k_m = k - b%Ex%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            j_m = j - b%Ex%YI
            Do i=FF%AbAr%I%com%Ex,FF%AbAr%I%fin%Ex
               i_m = i - b%Ex%XI
               do ii=1,FF%NumFreqs
                  FF%ExAb( i_m, j_m,ii) = FF%ExAb( i_m, j_m,ii) + FF%auxExp_E(ii) *  Ex( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Ey Down
         k = FF%AbAr%K%aba%Ey  !Down
         k_m = k - b%Ey%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            j_m = j - b%Ey%YI
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               i_m = i - b%Ey%XI
               do ii=1,FF%NumFreqs
                  FF%EyAb( i_m, j_m,ii) = FF%EyAb( i_m, j_m,ii) + FF%auxExp_E(ii) *  Ey( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !--->
      If( FF%farfieldAr) then
         !Ex Up
         k = FF%AbAr%K%arr%Ex  !Up
         k_m = k - b%Ex%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            j_m = j - b%Ex%YI
            do i = FF%AbAr%I%com%Ex, FF%AbAr%I%fin%Ex
               i_m = i - b%Ex%XI
               do ii=1,FF%NumFreqs
                  FF%ExAr( i_m, j_m,ii) = FF%ExAr( i_m, j_m,ii) + FF%auxExp_E(ii) *  Ex( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Ey Up
         k = FF%AbAr%K%arr%Ey  !Up
         k_m = k - b%Ey%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            j_m = j - b%Ey%YI
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               i_m = i - b%Ey%XI
               do ii=1,FF%NumFreqs
                  FF%EyAr( i_m, j_m,ii) = FF%EyAr( i_m, j_m,ii) + FF%auxExp_E(ii) *  Ey( i_m, j_m, k_m)

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !---------------------------> acaba UpdateFarFieldE <-----------------------------------------

      !!!!!!!!!!!!!!!!!!!!!!!!!!
      !magnetic field
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if( FF%farfieldTr) then
         !Hz Back
         i = FF%TrFr%I%tra%Hz  !Back
         i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            k_m = k - b%Hz%ZI
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               j_m = j - b%Hz%YI
               do ii=1,FF%NumFreqs
                  FF%HzTr(  j_m, k_m,ii) = FF%HzTr(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  ( Hz( i_m, j_m, k_m)  )
                  FF%HzTr2(  j_m, k_m,ii) = FF%HzTr2(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  ( Hz( i_m+1, j_m, k_m) )

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Hy Back
         i = FF%TrFr%I%tra%Hy  !Back
         i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            k_m = k - b%Hy%ZI
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               j_m = j - b%Hy%YI
               do ii=1,FF%NumFreqs
                  FF%HyTr(  j_m, k_m,ii) = FF%HyTr(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  (  Hy( i_m, j_m, k_m) )
                  FF%HyTr2(  j_m, k_m,ii) = FF%HyTr2(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  (  Hy( i_m+1, j_m, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !--->
      if( FF%farfieldFr) then
         !Hz  Front
         i = FF%TrFr%I%fro%Hz !Front
         i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            k_m = k - b%Hz%ZI
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               j_m = j - b%Hz%YI
               do ii=1,FF%NumFreqs
                  FF%HzFr(  j_m, k_m,ii) = FF%HzFr(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hz( i_m, j_m, k_m) )
                  FF%HzFr2(  j_m, k_m,ii) = FF%HzFr2(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hz( i_m-1, j_m, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Hy  Front
         i = FF%TrFr%I%fro%Hy !Front
         i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,j,k,j_m,k_m)
#endif
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            k_m = k - b%Hy%ZI
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               j_m = j - b%Hy%YI
               do ii=1,FF%NumFreqs
                  FF%HyFr(  j_m, k_m,ii) =  FF%HyFr(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hy( i_m, j_m, k_m) )
                  FF%HyFr2(  j_m, k_m,ii) =  FF%HyFr2(  j_m, k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hy( i_m - 1 , j_m, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !--->
      if( FF%farfieldIz) then
         !Hx Left
         j = FF%IzDe%J%izq%Hx  !Left
         j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            k_m = k - b%Hx%ZI
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               i_m = i - b%Hx%XI
               do ii=1,FF%NumFreqs
                  FF%HxIz( i_m,  k_m,ii) = FF%HxIz( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m, k_m) )
                  FF%HxIz2( i_m,  k_m,ii) = FF%HxIz2( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m +1, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Hz Left
         j = FF%IzDe%J%izq%Hz  !Left
         j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            k_m = k - b%Hz%ZI
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               i_m = i - b%Hz%XI
               do ii=1,FF%NumFreqs
                  FF%HzIz( i_m,  k_m,ii) = FF%HzIz( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hz( i_m, j_m, k_m))
                  FF%HzIz2( i_m,  k_m,ii) = FF%HzIz2( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hz( i_m, j_m +1, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !--->
      if( FF%farfieldDe) then
         !Hx  Right
         j = FF%IzDe%J%der%Hx !Right
         j_m = j - b%Hx%YI
         !--->
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            k_m = k - b%Hx%ZI
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               i_m = i - b%Hx%XI
               do ii=1,FF%NumFreqs
                  FF%HxDe( i_m,  k_m,ii) = FF%HxDe( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m, k_m))
                  FF%HxDe2( i_m,  k_m,ii) = FF%HxDe2( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m-1, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Hz  Right
         j = FF%IzDe%J%der%Hz !Right
         j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,k,i,k_m,i_m)
#endif
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            k_m = k - b%Hz%ZI
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               i_m = i - b%Hz%XI
               do ii=1,FF%NumFreqs
                  FF%HzDe( i_m,  k_m,ii) = FF%HzDe( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hz( i_m, j_m, k_m))
                  FF%HzDe2( i_m,  k_m,ii) = FF%HzDe2( i_m,  k_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hz( i_m, j_m-1, k_m))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!
      !--->
      if( FF%farfieldAb) then
         !Hx  Down
         k = FF%AbAr%K%aba%Hx  !Down
         k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            j_m = j - b%Hx%YI
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               i_m = i - b%Hx%XI
               do ii=1,FF%NumFreqs
                  FF%HxAb( i_m, j_m,ii) = FF%HxAb( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m, k_m))
                  FF%HxAb2( i_m, j_m,ii) = FF%HxAb2( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m, k_m+1))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Hy  Down
         k = FF%AbAr%K%aba%Hy  !Down
         k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            j_m = j - b%Hy%YI
            do i=FF%AbAr%I%com%Hy,FF%AbAr%I%fin%Hy
               i_m = i - b%Hy%XI
               do ii=1,FF%NumFreqs
                  FF%HyAb( i_m, j_m,ii) = FF%HyAb( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hy( i_m, j_m, k_m))
                  FF%HyAb2( i_m, j_m,ii) = FF%HyAb2( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hy( i_m, j_m, k_m+1))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !--->
      if( FF%farfieldAr) then
         !Hx Up
         k = FF%AbAr%K%arr%Hx  !Up
         k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            j_m = j - b%Hx%YI
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               i_m = i - b%Hx%XI
               do ii=1,FF%NumFreqs
                  FF%HxAr( i_m, j_m,ii) = FF%HxAr( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m, k_m))
                  FF%HxAr2( i_m, j_m,ii) = FF%HxAr2( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hx( i_m, j_m, k_m-1))

               end do
            end do
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         !Hy Up
         k=FF%AbAr%K%arr%Hy  !Up
         k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (ii,i,j,i_m,j_m)
#endif
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            j_m = j - b%Hy%YI
            do i = FF%AbAr%I%com%Hy, FF%AbAr%I%fin%Hy
               i_m = i - b%Hy%XI
               do ii=1,FF%NumFreqs
                  FF%HyAr( i_m, j_m,ii) =  FF%HyAr( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hy( i_m, j_m, k_m))
                  FF%HyAr2( i_m, j_m,ii) =  FF%HyAr2( i_m, j_m,ii) + FF%auxExp_H(ii) *  &
                  (   Hy( i_m, j_m, k_m-1))

               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif


      !solo los samples despues de 1 actualizan el valor
      !ver rutina dtft en postprocesws
      do ii=1,FF%NumFreqs
         FF%auxExp_E(ii)=FF%auxExp_E(ii) * FF%expIwdt(ii)
         FF%auxExp_H(ii)=FF%auxExp_H(ii) * FF%expIwdt(ii)
      end do

      !---------------------------> acaba UpdateFarFieldH <-----------------------------------------
      return

   endsubroutine UpdateFarField
   !
   !
   !
   subroutine StoreFarFields(B)
      type( bounds_t), intent( IN)  ::  b

      integer (kind=4)  ::  i, j, k, i_m, j_m, k_m,ii

      !---------------------------> empieza  <---------------------------------------
      do ii=1,FF%NumFreqs
         write (14,err=634) FF%auxExp_E(ii)
         write (14,err=634) FF%auxExp_H(ii)
      end do
      !
      If( FF%farfieldTr) then
         !Ez Back
         i = FF%TrFr%I%tra%Ez !Back
         i_m = i - b%Ez%XI
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            k_m = k - b%Ez%ZI
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               j_m = j - b%Ez%YI
               write (14,err=634) ( FF%EzTr(j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey Back
         i = FF%TrFr%I%tra%Ey  !Back
         i_m = i - b%Ey%XI
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            k_m = k - b%Ey%ZI
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               j_m = j - b%Ey%YI
               write (14,err=634) ( FF%EyTr(j_m, k_m,ii), ii=1,FF%NumFreqs)
            End do
         end do
      endif
      !--->
      If( FF%farfieldFr) then
         !Ez  Front
         i = FF%TrFr%I%fro%Ez !Front
         i_m = i - b%Ez%XI
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            k_m = k - b%Ez%ZI
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               j_m = j - b%Ez%YI
               write (14,err=634) ( FF%EzFr(j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey  Front
         i = FF%TrFr%I%fro%Ey !Front
         i_m = i - b%Ey%XI
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            k_m = k - b%Ey%ZI
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               j_m = j - b%Ey%YI
               write (14,err=634) ( FF%EyFr(j_m, k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldIz) then
         !Ex Left
         j = FF%IzDe%J%izq%Ex  !Left
         j_m = j - b%Ex%YI
         do k = FF%IzDe%K%com%Ex, FF%IzDe%K%fin%Ex
            k_m = k - b%Ex%ZI
            do i = FF%IzDe%I%com%Ex, FF%IzDe%I%fin%Ex
               i_m = i - b%Ex%XI
               write (14,err=634) ( FF%ExIz( i_m, k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ez Left
         j = FF%IzDe%J%izq%Ez  !Left
         j_m = j - b%Ez%YI
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            k_m = k - b%Ez%ZI
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               i_m = i - b%Ez%XI
               write (14,err=634) ( FF%EzIz( i_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldDe) then
         !Ez  Right
         j = FF%IzDe%J%der%Ez !Right
         j_m = j - b%Ez%YI
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            k_m = k - b%Ez%ZI
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               i_m = i - b%Ez%XI
               write (14,err=634) ( FF%EzDe( i_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ex  Right
         j = FF%IzDe%J%der%Ex !Right
         j_m = j - b%Ex%YI
         do k = FF%IzDe%K%com%Ex,FF%IzDe%K%fin%Ex
            k_m = k - b%Ex%ZI
            do i=FF%IzDe%I%com%Ex,FF%IzDe%I%fin%Ex
               i_m = i - b%Ex%XI
               write (14,err=634) ( FF%ExDe( i_m, k_m,ii) , ii=1,FF%NumFreqs)
            End do
         end do
      endif
      !--->
      If( FF%farfieldAb) then
         !Ex  Down
         k = FF%AbAr%K%aba%Ex  !Down
         k_m = k - b%Ex%ZI
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            j_m = j - b%Ex%YI
            Do i=FF%AbAr%I%com%Ex,FF%AbAr%I%fin%Ex
               i_m = i - b%Ex%XI
               write (14,err=634) ( FF%ExAb( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey Down
         k = FF%AbAr%K%aba%Ey  !Down
         k_m = k - b%Ey%ZI
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            j_m = j - b%Ey%YI
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               i_m = i - b%Ey%XI
               write (14,err=634) ( FF%EyAb( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldAr) then
         !Ex Up
         k = FF%AbAr%K%arr%Ex  !Up
         k_m = k - b%Ex%ZI
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            j_m = j - b%Ex%YI
            do i = FF%AbAr%I%com%Ex, FF%AbAr%I%fin%Ex
               i_m = i - b%Ex%XI
               write (14,err=634) ( FF%ExAr( i_m, j_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey Up
         k = FF%AbAr%K%arr%Ey  !Up
         k_m = k - b%Ey%ZI
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            j_m = j - b%Ey%YI
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               i_m = i - b%Ey%XI
               write (14,err=634) ( FF%EyAr( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !magneticos
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !--->
      if( FF%farfieldTr) then
         !Hz Back
         i = FF%TrFr%I%tra%Hz  !Back
         i_m = i - b%Hz%XI
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            k_m = k - b%Hz%ZI
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               j_m = j - b%Hz%YI
               write (14,err=634) ( FF%HzTr(  j_m, k_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HzTr2(  j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hy Back
         i = FF%TrFr%I%tra%Hy  !Back
         i_m = i - b%Hy%XI
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            k_m = k - b%Hy%ZI
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               j_m = j - b%Hy%YI
               write (14,err=634) ( FF%HyTr(  j_m, k_m,ii) , ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HyTr2(  j_m, k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldFr) then
         !Hz  Front
         i = FF%TrFr%I%fro%Hz !Front
         i_m = i - b%Hz%XI
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            k_m = k - b%Hz%ZI
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               j_m = j - b%Hz%YI
               write (14,err=634) ( FF%HzFr(  j_m, k_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HzFr2(  j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hy  Front
         i = FF%TrFr%I%fro%Hy !Front
         i_m = i - b%Hy%XI
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            k_m = k - b%Hy%ZI
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               j_m = j - b%Hy%YI
               write (14,err=634) ( FF%HyFr(  j_m, k_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HyFr2(  j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldIz) then
         !Hx Left
         j = FF%IzDe%J%izq%Hx  !Left
         j_m = j - b%Hx%YI
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            k_m = k - b%Hx%ZI
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               i_m = i - b%Hx%XI
               write (14,err=634) ( FF%HxIz( i_m,  k_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HxIz2( i_m,  k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hz Left
         j = FF%IzDe%J%izq%Hz  !Left
         j_m = j - b%Hz%YI
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            k_m = k - b%Hz%ZI
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               i_m = i - b%Hz%XI
               write (14,err=634) ( FF%HzIz( i_m,  k_m,ii) , ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HzIz2( i_m,  k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldDe) then
         !Hx  Right
         j = FF%IzDe%J%der%Hx !Right
         j_m = j - b%Hx%YI
         !--->
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            k_m = k - b%Hx%ZI
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               i_m = i - b%Hx%XI
               write (14,err=634) ( FF%HxDe( i_m,  k_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HxDe2( i_m,  k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hz  Right
         j = FF%IzDe%J%der%Hz !Right
         j_m = j - b%Hz%YI
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            k_m = k - b%Hz%ZI
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               i_m = i - b%Hz%XI
               write (14,err=634) ( FF%HzDe( i_m,  k_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HzDe2( i_m,  k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldAb) then
         !Hx  Down
         k = FF%AbAr%K%aba%Hx  !Down
         k_m = k - b%Hx%ZI
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            j_m = j - b%Hx%YI
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               i_m = i - b%Hx%XI
               write (14,err=634) ( FF%HxAb( i_m, j_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HxAb2( i_m, j_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hy  Down
         k = FF%AbAr%K%aba%Hy  !Down
         k_m = k - b%Hy%ZI
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            j_m = j - b%Hy%YI
            do i=FF%AbAr%I%com%Hy,FF%AbAr%I%fin%Hy
               i_m = i - b%Hy%XI
               write (14,err=634) ( FF%HyAb( i_m, j_m,ii) , ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HyAb2( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldAr) then
         !Hx Up
         k = FF%AbAr%K%arr%Hx  !Up
         k_m = k - b%Hx%ZI
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            j_m = j - b%Hx%YI
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               i_m = i - b%Hx%XI
               write (14,err=634) ( FF%HxAr( i_m, j_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) ( FF%HxAr2( i_m, j_m,ii), ii=1,FF%NumFreqs)
            end do
         enddo
         !Hy Up
         k=FF%AbAr%K%arr%Hy  !Up
         k_m = k - b%Hy%ZI
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            j_m = j - b%Hy%YI
            do i = FF%AbAr%I%com%Hy, FF%AbAr%I%fin%Hy
               i_m = i - b%Hy%XI
               write (14,err=634) (FF%HyAr( i_m, j_m,ii), ii=1,FF%NumFreqs)
               write (14,err=634) (FF%HyAr2( i_m, j_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif


      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'FARFIELD: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   endsubroutine
   !

   subroutine ReadFarfield(b)
      type( bounds_t), intent( IN)  ::  b


      integer ::  i, j, k, i_m, j_m, k_m,ii

      !---------------------------> empieza  <---------------------------------------
      do ii=1,FF%NumFreqs
         READ (14) FF%auxExp_E(ii)
         READ (14) FF%auxExp_H(ii)
      end do
      If( FF%farfieldTr) then
         !Ez Back
         i = FF%TrFr%I%tra%Ez !Back
         i_m = i - b%Ez%XI
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            k_m = k - b%Ez%ZI
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               j_m = j - b%Ez%YI
               READ (14) ( FF%EzTr(j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey Back
         i = FF%TrFr%I%tra%Ey  !Back
         i_m = i - b%Ey%XI
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            k_m = k - b%Ey%ZI
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               j_m = j - b%Ey%YI
               READ (14) ( FF%EyTr(j_m, k_m,ii), ii=1,FF%NumFreqs)
            End do
         end do
      endif
      !--->
      If( FF%farfieldFr) then
         !Ez  Front
         i = FF%TrFr%I%fro%Ez !Front
         i_m = i - b%Ez%XI
         do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ez
            k_m = k - b%Ez%ZI
            do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
               j_m = j - b%Ez%YI
               READ (14) ( FF%EzFr(j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey  Front
         i = FF%TrFr%I%fro%Ey !Front
         i_m = i - b%Ey%XI
         do k = FF%TrFr%K%com%Ey, FF%TrFr%K%fin%Ey
            k_m = k - b%Ey%ZI
            do j = FF%TrFr%J%com%Ey, FF%TrFr%J%fin%Ey
               j_m = j - b%Ey%YI
               READ (14) ( FF%EyFr(j_m, k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldIz) then
         !Ex Left
         j = FF%IzDe%J%izq%Ex  !Left
         j_m = j - b%Ex%YI
         do k = FF%IzDe%K%com%Ex, FF%IzDe%K%fin%Ex
            k_m = k - b%Ex%ZI
            do i = FF%IzDe%I%com%Ex, FF%IzDe%I%fin%Ex
               i_m = i - b%Ex%XI
               READ (14) ( FF%ExIz( i_m, k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ez Left
         j = FF%IzDe%J%izq%Ez  !Left
         j_m = j - b%Ez%YI
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            k_m = k - b%Ez%ZI
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               i_m = i - b%Ez%XI
               READ (14) ( FF%EzIz( i_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldDe) then
         !Ez  Right
         j = FF%IzDe%J%der%Ez !Right
         j_m = j - b%Ez%YI
         do k = FF%IzDe%K%com%Ez, FF%IzDe%K%fin%Ez
            k_m = k - b%Ez%ZI
            do i = FF%IzDe%I%com%Ez, FF%IzDe%I%fin%Ez
               i_m = i - b%Ez%XI
               READ (14) ( FF%EzDe( i_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ex  Right
         j = FF%IzDe%J%der%Ex !Right
         j_m = j - b%Ex%YI
         do k = FF%IzDe%K%com%Ex,FF%IzDe%K%fin%Ex
            k_m = k - b%Ex%ZI
            do i=FF%IzDe%I%com%Ex,FF%IzDe%I%fin%Ex
               i_m = i - b%Ex%XI
               READ (14) ( FF%ExDe( i_m, k_m,ii) , ii=1,FF%NumFreqs)
            End do
         end do
      endif
      !--->
      If( FF%farfieldAb) then
         !Ex  Down
         k = FF%AbAr%K%aba%Ex  !Down
         k_m = k - b%Ex%ZI
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            j_m = j - b%Ex%YI
            Do i=FF%AbAr%I%com%Ex,FF%AbAr%I%fin%Ex
               i_m = i - b%Ex%XI
               READ (14) ( FF%ExAb( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey Down
         k = FF%AbAr%K%aba%Ey  !Down
         k_m = k - b%Ey%ZI
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            j_m = j - b%Ey%YI
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               i_m = i - b%Ey%XI
               READ (14) ( FF%EyAb( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      If( FF%farfieldAr) then
         !Ex Up
         k = FF%AbAr%K%arr%Ex  !Up
         k_m = k - b%Ex%ZI
         do j = FF%AbAr%J%com%Ex, FF%AbAr%J%fin%Ex
            j_m = j - b%Ex%YI
            do i = FF%AbAr%I%com%Ex, FF%AbAr%I%fin%Ex
               i_m = i - b%Ex%XI
               READ (14) ( FF%ExAr( i_m, j_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Ey Up
         k = FF%AbAr%K%arr%Ey  !Up
         k_m = k - b%Ey%ZI
         do j = FF%AbAr%J%com%Ey, FF%AbAr%J%fin%Ey
            j_m = j - b%Ey%YI
            do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
               i_m = i - b%Ey%XI
               READ (14) ( FF%EyAr( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !magneticos
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !--->
      if( FF%farfieldTr) then
         !Hz Back
         i = FF%TrFr%I%tra%Hz  !Back
         i_m = i - b%Hz%XI
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            k_m = k - b%Hz%ZI
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               j_m = j - b%Hz%YI
               READ (14) ( FF%HzTr(  j_m, k_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HzTr2(  j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hy Back
         i = FF%TrFr%I%tra%Hy  !Back
         i_m = i - b%Hy%XI
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            k_m = k - b%Hy%ZI
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               j_m = j - b%Hy%YI
               READ (14) ( FF%HyTr(  j_m, k_m,ii) , ii=1,FF%NumFreqs)
               READ (14) ( FF%HyTr2(  j_m, k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldFr) then
         !Hz  Front
         i = FF%TrFr%I%fro%Hz !Front
         i_m = i - b%Hz%XI
         do k = FF%TrFr%K%com%Hz, FF%TrFr%K%fin%Hz
            k_m = k - b%Hz%ZI
            do j = FF%TrFr%J%com%Hz, FF%TrFr%J%fin%Hz
               j_m = j - b%Hz%YI
               READ (14) ( FF%HzFr(  j_m, k_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HzFr2(  j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hy  Front
         i = FF%TrFr%I%fro%Hy !Front
         i_m = i - b%Hy%XI
         do k = FF%TrFr%K%com%Hy, FF%TrFr%K%fin%Hy
            k_m = k - b%Hy%ZI
            do j = FF%TrFr%J%com%Hy, FF%TrFr%J%fin%Hy
               j_m = j - b%Hy%YI
               READ (14) ( FF%HyFr(  j_m, k_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HyFr2(  j_m, k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldIz) then
         !Hx Left
         j = FF%IzDe%J%izq%Hx  !Left
         j_m = j - b%Hx%YI
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            k_m = k - b%Hx%ZI
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               i_m = i - b%Hx%XI
               READ (14) ( FF%HxIz( i_m,  k_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HxIz2( i_m,  k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hz Left
         j = FF%IzDe%J%izq%Hz  !Left
         j_m = j - b%Hz%YI
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            k_m = k - b%Hz%ZI
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               i_m = i - b%Hz%XI
               READ (14) ( FF%HzIz( i_m,  k_m,ii) , ii=1,FF%NumFreqs)
               READ (14) ( FF%HzIz2( i_m,  k_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldDe) then
         !Hx  Right
         j = FF%IzDe%J%der%Hx !Right
         j_m = j - b%Hx%YI
         !--->
         do k = FF%IzDe%K%com%Hx, FF%IzDe%K%fin%Hx
            k_m = k - b%Hx%ZI
            do i = FF%IzDe%I%com%Hx, FF%IzDe%I%fin%Hx
               i_m = i - b%Hx%XI
               READ (14) ( FF%HxDe( i_m,  k_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HxDe2( i_m,  k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hz  Right
         j = FF%IzDe%J%der%Hz !Right
         j_m = j - b%Hz%YI
         do k = FF%IzDe%K%com%Hz, FF%IzDe%K%fin%Hz
            k_m = k - b%Hz%ZI
            do i = FF%IzDe%I%com%Hz, FF%IzDe%I%fin%Hz
               i_m = i - b%Hz%XI
               READ (14) ( FF%HzDe( i_m,  k_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HzDe2( i_m,  k_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldAb) then
         !Hx  Down
         k = FF%AbAr%K%aba%Hx  !Down
         k_m = k - b%Hx%ZI
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            j_m = j - b%Hx%YI
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               i_m = i - b%Hx%XI
               READ (14) ( FF%HxAb( i_m, j_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HxAb2( i_m, j_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
         !Hy  Down
         k = FF%AbAr%K%aba%Hy  !Down
         k_m = k - b%Hy%ZI
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            j_m = j - b%Hy%YI
            do i=FF%AbAr%I%com%Hy,FF%AbAr%I%fin%Hy
               i_m = i - b%Hy%XI
               READ (14) ( FF%HyAb( i_m, j_m,ii) , ii=1,FF%NumFreqs)
               READ (14) ( FF%HyAb2( i_m, j_m,ii) , ii=1,FF%NumFreqs)
            enddo
         enddo
      endif
      !--->
      if( FF%farfieldAr) then
         !Hx Up
         k = FF%AbAr%K%arr%Hx  !Up
         k_m = k - b%Hx%ZI
         do j = FF%AbAr%J%com%Hx, FF%AbAr%J%fin%Hx
            j_m = j - b%Hx%YI
            do i = FF%AbAr%I%com%Hx, FF%AbAr%I%fin%Hx
               i_m = i - b%Hx%XI
               READ (14) ( FF%HxAr( i_m, j_m,ii), ii=1,FF%NumFreqs)
               READ (14) ( FF%HxAr2( i_m, j_m,ii), ii=1,FF%NumFreqs)
            end do
         enddo
         !Hy Up
         k=FF%AbAr%K%arr%Hy  !Up
         k_m = k - b%Hy%ZI
         do j = FF%AbAr%J%com%Hy, FF%AbAr%J%fin%Hy
            j_m = j - b%Hy%YI
            do i = FF%AbAr%I%com%Hy, FF%AbAr%I%fin%Hy
               i_m = i - b%Hy%XI
               READ (14) (FF%HyAr( i_m, j_m,ii), ii=1,FF%NumFreqs)
               READ (14) (FF%HyAr2( i_m, j_m,ii), ii=1,FF%NumFreqs)
            enddo
         enddo
      endif

      return
   endsubroutine ReadFarField

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Calculate the physical point coordinates of index i,j,k
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!!Subroutine PuntoPhys(nfield,i,j,k,xf,yf,zf)
   !!!!integer (KIND=4) i,j,k,nfield
   !!!!REAL (KIND=RKIND)  ::   xf,yf,zf
   !!!!!
   !!!!xf=FF%Punto%PhysCoor(nfield)%x(i)
   !!!!yf=FF%Punto%PhysCoor(nfield)%y(j)
   !!!!zf=FF%Punto%PhysCoor(nfield)%z(k)
   !!!!
   !!!!return
   !!!!!
   !!!!end Subroutine PuntoPhys
   !embebido en el codigo para evitar errores con OMP 29/01/15
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine Destroyfarfield
      integer (kind=4)  ::  field

      do field=iEx,iHz
         if (associated(FF%Punto%PhysCoor(field)%x)) deallocate (FF%Punto%PhysCoor(field)%x)
         if (associated(FF%Punto%PhysCoor(field)%y)) deallocate (FF%Punto%PhysCoor(field)%y)
         if (associated(FF%Punto%PhysCoor(field)%z)) deallocate (FF%Punto%PhysCoor(field)%z)
      end do
      !
      if (allocated(FF%expIwdt)) deALLOCATE (FF%expIwdt)
      !
      if( FF%farfieldIz) then
         if (allocated(FF%ExIz)) deallocate (FF%ExIz,FF%EzIz,FF%HxIz,FF%HzIz,FF%HxIz2,FF%HzIz2)
      endif
      if( FF%farfieldDe) then
         if (allocated(FF%ExDe)) deallocate (FF%ExDe,FF%EzDe,FF%HxDe,FF%HzDe,FF%HxDe2,FF%HzDe2)
      endif
      if( FF%farfieldAb) then
         if (allocated(FF%ExAb)) deallocate (FF%ExAb,FF%EyAb,FF%HxAb,FF%HyAb,FF%HxAb2,FF%HyAb2)
      endif
      if( FF%farfieldAr) then
         if (allocated(FF%ExAr)) deallocate (FF%ExAr,FF%EyAr,FF%HxAr,FF%HyAr,FF%HxAr2,FF%HyAr2)
      endif
      if( FF%farfieldFr) then
         if (allocated(FF%EyFr)) deallocate (FF%EyFr,FF%EzFr,FF%HyFr,FF%HzFr,FF%HyFr2,FF%HzFr2)
      endif
      if( FF%farfieldTr) then
         if (allocated(FF%EyTr)) deallocate (FF%EyTr,FF%EzTr,FF%HyTr,FF%HzTr,FF%HyTr2,FF%HzTr2)
      endif

   end subroutine Destroyfarfield


   !!!!!!!!!!!!!!!!!!!!!!!calculo del far field propiamente dicho y volcado en fichero

   subroutine FlushFarfield(layoutnumber,size, b, dxe, dye, dze, dxh, dyh, dzh,facesNF2FF,rinstant)
      type (co_t) :: co,new_co
      type (nf2ff_t) :: facesNF2FF
      type( bounds_t), intent( IN)      ::  b
      real (kind = RKIND), dimension( 0 :  b%dxe%NX-1), intent( IN)  ::  dxe
      real (kind = RKIND), dimension( 0 :  b%dye%NY-1), intent( IN)  ::  dye
      real (kind = RKIND), dimension( 0 :  b%dze%NZ-1), intent( IN)  ::  dze
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN)  ::  dxh
      real (kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN)  ::  dyh
      real (kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN)  ::  dzh

      integer (kind=4), intent(in) :: layoutnumber,size

      integer (kind=4)  ::  ntheta,nphi,ithe,iphi,ii,i,j,k,i_m,j_m,k_m,pasadas
      real (kind = RKIND) :: theta,phi,sintheta_sinphi,sintheta_cosphi, &
      costheta,cosphi,costheta_cosphi,costheta_sinphi,sintheta,sinphi,&
      freq, NORMAL, SIGNO,  dummy,newdummy1,newdummy2,RCS(1:2)
      real (kind = RKIND_tiempo) :: rinstant
      integer (kind=4)  ::  ierr,pozi,donde
      complex( kind = CKIND) :: L_theta,L_phi,N_theta,N_phi,Etheta(1:2),Ephi(1:2),Mx,My,Mz,Jx,Jy,Jz,comun
      complex( kind = CKIND) :: new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz
      complex( kind = CKIND) :: L_theta_final,L_phi_final,N_theta_final,N_phi_final
      character (LEN=BUFSIZE)     ::  dubuf
      logical :: GOahead
      complex( kind = CKIND), dimension(:,:,:), pointer :: EcampoX,EcampoY,EcampoZ,HcampoX,HcampoY,HcampoZ,Hcampo2X,Hcampo2Y,Hcampo2Z
      !---------------------------> empieza  <---------------------------------------
      integer (kind=4)  :: number,mt,mf,contaje,m,n
      real (kind = RKIND), allocatable :: Phimatrix(:,:),sizephi(:),Thetavector(:)
      real (kind = RKIND) :: thetaini,thetafin,phiini,phifin,solid,ddd,aaa,dfi,dth
      integer :: my_iostat
!
      character(LEN=BUFSIZE) :: chninstant
      write(chninstant,fmt) rinstant

      write(dubuf,'(a)')  ' NF2FF: INIT ' ! mpimaster=',LAYOUTNUMBER
      if (layoutnumber == 0) call print11(layoutnumber,dubuf,.TRUE.)

#ifdef CompileWithMPI
      if (FF%MPIRoot == layoutnumber) then
#endif
      open (FF%unitfarfield,file=trim(adjustl(FF%filefarfield)))
      CLOSE (FF%unitfarfield, STATUS='delete')
      
      
      my_iostat=0
9138  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9138,'.',layoutnumber,trim(adjustl(FF%filefarfield))
      open (FF%unitfarfield,file=trim(adjustl(FF%filefarfield)),form='formatted',position='append',err=9138,iostat=my_iostat,status='new',action='write') !lista de todos los .h5bin   
      write(FF%unitfarfield,'(a)') ' f_at_'//trim(adjustl(chninstant))//'   Theta    Phi    Etheta_mod    Etheta_phase    Ephi_mod    Ephi_phase    RCS(ARIT) RCS(GEOM)'

#ifdef CompileWithMPI
   endif
#endif
      !calculo y escritura
      ntheta=int(abs(FF%ThetaStop-FF%ThetaStart)/FF%ThetaStep)+1
      nphi=  int(abs(FF%PhiStop-FF%PhiStart)/FF%PhiStep)+1

!
!!!!algoritmo de reparto por angulo solido de puntos 080317 !basado en el script matlab esfera.m
      number=ntheta*nphi;
      thetaini=FF%ThetaStart
      thetafin=FF%ThetaStop;
      if (thetafin==thetaini) thetafin=thetaini+0.01_RKIND
      phiini=FF%PhiStart;
      phifin=FF%PhiStop;
      if (phifin==phiini) phifin=phiini+0.01_RKIND
      solid=abs((phifin-phiini)*(cos(thetaini)-cos(thetafin)))
      aaa=solid/number;
      ddd=sqrt(aaa);
      mt=nint((thetafin-thetaini)/ddd);
      if (mt<=0) mt=1
      dth=(thetafin-thetaini)/mt;
      dfi=aaa/dth;
      allocate (Thetavector(1:mt+1))
      allocate (sizephi(1:mt+1))
      do contaje=0,1
         if (contaje==1) allocate (Phimatrix(1:mt+1,1:int(maxval(sizephi)))) !tiro espacio pero no hay otro remedio
         do m=1,mt+1
              Thetavector(m)=thetaini+(thetafin-thetaini)*((m-1))/mt;
              mf=nint((phifin-phiini)*sin(Thetavector(m))/dfi);
              if (mf<=0) mf=1
              sizephi(m)=mf+1;
              if (contaje==1) then
                  do n=1,mf+1
                        if ((n==1).and.(m==1)) then
                             Phimatrix(m,n)=phiini
                        else
                            Phimatrix(m,n)=phiini+(phifin-phiini)*(n-1)/mf
                        endif
                  end do
              endif
         end do
     end do
!!!!! fin cambios 080317

      pozi=index(FF%filefarfield,'_log_')


      do ii=1,FF%NumFreqs
         freq=FF%InitialFreq + (ii-1)*FF%FreqStep
         if (pozi /=0) then !logaritmico
            freq=10.0_RKIND **freq
         endif
!
         write(dubuf,'(a,i9,a,i9,a,e19.9e3)')  ' NF2FF: Start processing freq (',ii,'/',FF%NumFreqs,')= ',freq
         call print11(layoutnumber,dubuf,.TRUE.)
!
         comun=(0.0_RKIND,1.0_RKIND)*2.0_RKIND * pi*freq/cluz
!el barrido isoesferico no me funciona bien todavia 080317
!!!         do ithe=1,mt
!!!          theta=Thetavector(ithe)
!!!          do iphi=1,sizephi(ithe)
!!!             phi=Phimatrix(ithe,iphi)
         !do ithe=0,ntheta-1
         !   theta=FF%ThetaStart+ithe*FF%ThetaStep
         !   do iphi=0,nphi-1
         !      phi=FF%phiStart+iphi*FF%PhiStep
           theta=FF%ThetaStart-FF%ThetaStep
           do while (theta<FF%ThetaStop)
             theta=min(theta+FF%ThetaStep,FF%ThetaStop)
             phi=FF%phiStart-FF%PhiStep
             do while (phi<FF%phiStop)
               phi=min(phi+FF%PhiStep,FF%PhiStop)
!
               do pasadas=2,1,-1 !PRIMERO aritmetica (LOS CAMPOS SON ARITMETICOS) Y OTRA GEOMETRICA (LAS RCS SON AMBAS)
                  !
                  sintheta_sinphi= Sin(theta)*sin(phi)
                  sintheta_cosphi= Sin(theta)*cos(phi)
                  costheta=    Cos(theta)
                  cosphi  =    Cos(phi)
                  costheta_cosphi= Cos(theta)*cos(phi)
                  costheta_sinphi= cos(theta)*sin(phi)
                  sintheta=   Sin(theta)
                  sinphi  =   Sin(phi)
                  !
                  !para curarme en salud de los problemas de reduccion 29/01/15
                  L_theta_final=0.0_RKIND
                  L_phi_final=0.0_RKIND
                  N_theta_final=0.0_RKIND
                  N_phi_final=0.0_RKIND
                  !
                  !!!!!!!!!!!!!
                  ! Back front


                  L_theta=0.0_RKIND ; L_phi=0.0_RKIND ; N_theta=0.0_RKIND ; N_phi=0.0_RKIND ;
                  do donde=1,2
                     co%x_Mx=0;co%y_Mx=0;co%z_Mx=0;
                     co%x_My=0;co%y_My=0;co%z_My=0;
                     co%x_Mz=0;co%y_Mz=0;co%z_Mz=0;
                     co%x_Jx=0;co%y_Jx=0;co%z_Jx=0;
                     co%x_Jy=0;co%y_Jy=0;co%z_Jy=0;
                     co%x_Jz=0;co%y_Jz=0;co%z_Jz=0;
                     if (donde==1) then
                        i = FF%TrFr%I%tra%Ez !Back !el del Ey coincide. Lo hago asi para no picar tanto codigo!!!
                        normal=-1.0_RKIND
                        GOahead = ( FF%farfieldTr .and. facesNF2FF%Tr)
                        EcampoZ =>  FF%EzTr
                        EcampoY =>  FF%EyTr
                        HcampoZ =>  FF%HzTr
                        HcampoY =>  FF%HyTr
                        Hcampo2Z => FF%HzTr2
                        Hcampo2Y => FF%HyTr2
                     else
                        i = FF%TrFr%I%fro%Ez !Front !el del Ey coincide. Lo hago asi para no picar tanto codigo!!!
                        GOahead = ( FF%farfieldFr .and. facesNF2FF%Fr)
                        normal=+1.0_RKIND
                        EcampoZ  => FF%EzFr
                        EcampoY  => FF%EyFr
                        HcampoZ  => FF%HzFr
                        HcampoY  => FF%HyFr
                        Hcampo2Z => FF%HzFr2
                        Hcampo2Y => FF%HyFr2
                     endif
                     if (GOahead) then
#ifdef CompileWithOpenMP
!$OMP                 PARALLEL DO DEFAULT(SHARED) private (co,new_co,Mx,My,Mz,Jx,Jy,Jz,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,j,k,j_m,k_m) REDUCTION(+:L_theta) REDUCTION(+:L_phi) REDUCTION(+:N_theta) REDUCTION(+:N_phi)
#endif
                        do k = FF%TrFr%K%com%Ez, FF%TrFr%K%fin%Ey
                           k_m = k - b%Ez%ZI
                           do j = FF%TrFr%J%com%Ez, FF%TrFr%J%fin%Ez
                              j_m = j - b%Ez%YI
                              Mx=0.0_RKIND; My=0.0_RKIND; Mz=0.0_RKIND; Jx=0.0_RKIND; Jy=0.0_RKIND; Jz=0.0_RKIND
                              if (k.le.FF%TrFr%K%fin%Ez) My = - EcampoZ( j_m, k_m,ii) *dyh(j_m)*dze(k_m)*NORMAL !los finales si varian
                              if (j.le.FF%TrFr%J%fin%Ey) Mz = + EcampoY( j_m, k_m,ii) *dye(j_m)*dzh(k_m)*NORMAL !los finales si varian
                              if (k.le.FF%TrFr%K%fin%Hz) Jy = + (Average(pasadas, HcampoZ( j_m, k_m,ii) , Hcampo2Z( j_m, k_m,ii))) *dye(j_m)*dzh(k_m)*NORMAL
                              if (j.le.FF%TrFr%J%fin%Hy) Jz = - (Average(pasadas, HcampoY( j_m, k_m,ii) , Hcampo2Y( j_m, k_m,ii))) *dyh(j_m)*dze(k_m)*NORMAL
                              co%x_My=FF%Punto%PhysCoor(iEz)%x(i); co%y_My=FF%Punto%PhysCoor(iEz)%y(j); co%z_My=FF%Punto%PhysCoor(iEz)%z(k)
                              co%x_Mz=FF%Punto%PhysCoor(iEy)%x(i); co%y_Mz=FF%Punto%PhysCoor(iEy)%y(j); co%z_Mz=FF%Punto%PhysCoor(iEy)%z(k)
                              co%x_Jy=co%x_Mz;                     co%y_Jy=co%y_Mz;                     co%z_Jy=co%z_Mz;
                              co%x_Jz=co%x_My;                     co%y_Jz=co%y_My;                     co%z_Jz=co%z_My;
                              call update_LN(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi)
                              !! simetrias
                              !!!la trfr hay que llamarla para cada caso
                              new_My = My
                              new_Mz = Mz
                              new_Jy = Jy
                              new_Jz = Jz
                              new_co = co
                              call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              !
                              If( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN) then
                                 new_My = + My
                                 new_Mz = - Mz
                                 new_Jy = - Jy
                                 new_Jz = + Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%z_Jy= new_co%z_Mz;
                                 new_co%z_Jz= new_co%z_My;
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldTr_ClonePMC_DOWN .or. FF%farfieldFr_ClonePMC_DOWN) then
                                 new_My = - My
                                 new_Mz = + Mz
                                 new_Jy = + Jy
                                 new_Jz = - Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%z_Jy= new_co%z_Mz;
                                 new_co%z_Jz= new_co%z_My;
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP) then
                                 new_My = + My
                                 new_Mz = - Mz
                                 new_Jy = - Jy
                                 new_Jz = + Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%z_Jy= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP) then
                                 new_My = - My
                                 new_Mz = + Mz
                                 new_Jy = + Jy
                                 new_Jz = - Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%z_Jy= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT) then
                                 !!!!!!!!!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
                                 new_My = - My
                                 new_Mz = + Mz
                                 new_Jy = + Jy
                                 new_Jz = - Jz
                                 new_co%y_My=     -co%y_My +FF%YOffsetMinus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetMinus
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT) then
                                 new_My = + My
                                 new_Mz = - Mz
                                 new_Jy = - Jy
                                 new_Jz = + Jz
                                 new_co%y_My=     -co%y_My +FF%YOffsetMinus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetMinus
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT) then
                                 new_My = - My
                                 new_Mz = + Mz
                                 new_Jy = + Jy
                                 new_Jz = - Jz
                                 new_co%y_My=     -co%y_My +FF%YOffsetPlus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetPlus
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT) then
                                 new_My = + My
                                 new_Mz = - Mz
                                 new_Jy = - Jy
                                 new_Jz = + Jz
                                 new_co%y_My=     -co%y_My +FF%YOffsetPlus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetPlus
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              !!!!!!!!!CASOS MIXTOS esquinas
                              If ((( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN).and.( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT)).or. &
                              (( FF%farfieldTr_ClonePMC_DOWN.or.FF%farfieldFr_ClonePMC_DOWN).and.( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT)).or. &
                              (( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN).and.( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT)).or. &
                              (( FF%farfieldTr_ClonePMC_DOWN.or.FF%farfieldFr_ClonePMC_DOWN).and.( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT)) )  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN).and.( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT)).or. &
                                 (( FF%farfieldTr_ClonePMC_DOWN.or.FF%farfieldFr_ClonePMC_DOWN).and.( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT)) )  signo=-1.0_RKIND
                                 new_My = signo * My
                                 new_Mz = signo * Mz
                                 new_Jy = signo * Jy
                                 new_Jz = signo * Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%y_My=     -co%y_My +FF%YOffsetMinus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetMinus
                                 new_co%z_Jy= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_My
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              If ((( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN).and.( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT)).or. &
                              (( FF%farfieldTr_ClonePMC_DOWN.or.FF%farfieldFr_ClonePMC_DOWN).and.( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT)).or. &
                              (( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN).and.( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT)).or. &
                              (( FF%farfieldTr_ClonePMC_DOWN.or.FF%farfieldFr_ClonePMC_DOWN).and.( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT)) ) THEN
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldTr_ClonePEC_DOWN.or.FF%farfieldFr_ClonePEC_DOWN).and.( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT)).or. &
                                 (( FF%farfieldTr_ClonePMC_DOWN.or.FF%farfieldFr_ClonePMC_DOWN).and.( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT)) ) signo=-1.0_RKIND
                                 new_My = signo * My
                                 new_Mz = signo * Mz
                                 new_Jy = signo * Jy
                                 new_Jz = signo * Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%y_My=     -co%y_My +FF%YOffsetPlus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetPlus
                                 new_co%z_Jy= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_My
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF

                              If ((( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP).and.( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT)).or. &
                              (( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP).and.( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT)).or. &
                              (( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP).and.( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT)).or. &
                              (( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP).and.( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT)) ) THEN
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP).and.( FF%farfieldTr_ClonePEC_LEFT.or.FF%farfieldFr_ClonePEC_LEFT)).or. &
                                 (( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP).and.( FF%farfieldTr_ClonePMC_LEFT.or.FF%farfieldFr_ClonePMC_LEFT)) )  signo=-1.0_RKIND
                                 new_My = signo * My
                                 new_Mz = signo * Mz
                                 new_Jy = signo * Jy
                                 new_Jz = signo * Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%y_My=     -co%y_My +FF%YOffsetMinus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetMinus
                                 new_co%z_Jy= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_My
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              If ((( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP).and.( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT)).or. &
                              (( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP).and.( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT)).or. &
                              (( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP).and.( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT)).or. &
                              (( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP).and.( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT))) THEN
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldTr_ClonePEC_UP.or.FF%farfieldFr_ClonePEC_UP).and.( FF%farfieldTr_ClonePEC_RIGHT.or.FF%farfieldFr_ClonePEC_RIGHT)).or. &
                                 (( FF%farfieldTr_ClonePMC_UP.or.FF%farfieldFr_ClonePMC_UP).and.( FF%farfieldTr_ClonePMC_RIGHT.or.FF%farfieldFr_ClonePMC_RIGHT)))  signo=-1.0_RKIND
                                 new_My = signo * My
                                 new_Mz = signo * Mz
                                 new_Jy = signo * Jy
                                 new_Jz = signo * Jz
                                 new_co%z_My=     -co%z_My +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%y_My=     -co%y_My +FF%YOffsetPlus
                                 new_co%y_Mz=     -co%y_Mz +FF%YOffsetPlus
                                 new_co%z_Jy= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_My
                                 new_co%y_Jy= new_co%y_Mz
                                 new_co%y_Jz= new_co%y_My
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneTrFr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !!! fin simetrias
                           enddo
                        enddo
#ifdef CompileWithOpenMP
!$OMP                 END PARALLEL DO
!$OMP BARRIER
#endif
                        L_theta_final = L_theta_final + L_theta ; L_phi_final   = L_phi_final   + L_phi
                        N_theta_final = N_theta_final + N_theta ; N_phi_final   = N_phi_final   + N_phi
                        L_theta=0.0_RKIND ; L_phi=0.0_RKIND ; N_theta=0.0_RKIND ; N_phi=0.0_RKIND ;
                     endif !del goAhead
                  end do !del donde
                  !--->
                  !--->
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  ! Left Right
                  L_theta=0.0_RKIND ; L_phi=0.0_RKIND ; N_theta=0.0_RKIND ; N_phi=0.0_RKIND ;
                  do donde=1,2
                     co%x_Mx=0;co%y_Mx=0;co%z_Mx=0;
                     co%x_My=0;co%y_My=0;co%z_My=0;
                     co%x_Mz=0;co%y_Mz=0;co%z_Mz=0;
                     co%x_Jx=0;co%y_Jx=0;co%z_Jx=0;
                     co%x_Jy=0;co%y_Jy=0;co%z_Jy=0;
                     co%x_Jz=0;co%y_Jz=0;co%z_Jz=0;
                     if (donde==1) then
                        j = FF%IzDe%J%izq%Ex
                        normal=-1.0_RKIND
                        GOahead = ( FF%farfieldIz .and. facesNF2FF%Iz)
                        EcampoZ  => FF%EzIz
                        EcampoX  => FF%ExIz
                        HcampoZ  => FF%HzIz
                        HcampoX  => FF%HxIz
                        Hcampo2Z => FF%HzIz2
                        Hcampo2X => FF%HxIz2
                     else
                        j = FF%IzDe%J%der%Ex
                        normal=+1.0_RKIND
                        GOahead = ( FF%farfieldDe .and. facesNF2FF%De)
                        EcampoZ  => FF%EzDe
                        EcampoX  => FF%ExDe
                        HcampoZ  => FF%HzDe
                        HcampoX  => FF%HxDe
                        Hcampo2Z => FF%HzDe2
                        Hcampo2X => FF%HxDe2
                     endif
                     if (GOahead) then
#ifdef CompileWithOpenMP
!$OMP                 PARALLEL DO DEFAULT(SHARED) private (co,new_co,Mx,My,Mz,Jx,Jy,Jz,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,k,i,k_m,i_m) REDUCTION(+:L_theta) REDUCTION(+:L_phi) REDUCTION(+:N_theta) REDUCTION(+:N_phi)
#endif
                        do k = FF%IzDe%K%com%Ex, FF%IzDe%K%fin%Ex
                           k_m = k - b%Ex%ZI
                           do i = FF%IzDe%I%com%Ex, FF%IzDe%I%fin%Ez
                              i_m = i - b%Ex%XI
                              Mx=0.0_RKIND; My=0.0_RKIND; Mz=0.0_RKIND; Jx=0.0_RKIND; Jy=0.0_RKIND; Jz=0.0_RKIND
                              if (k.le.FF%IzDe%K%fin%Ez)  Mx = + EcampoZ( i_m, k_m,ii) *dxh(i_m)*dze(k_m)*NORMAL
                              if (i.le.FF%IzDe%I%fin%Ex)  Mz = - EcampoX( i_m, k_m,ii) *dxe(i_m)*dzh(k_m)*NORMAL
                              if (k.le.FF%IzDe%K%fin%Hz)  Jx = - (Average(pasadas, HcampoZ( i_m, k_m,ii) , Hcampo2Z( i_m, k_m,ii))) *dxe(i_m)*dzh(k_m)*NORMAL
                              if (i.le.FF%IzDe%I%fin%Hx)  Jz = + (Average(pasadas, HcampoX( i_m, k_m,ii) , Hcampo2X( i_m, k_m,ii))) *dxh(i_m)*dze(k_m)*NORMAL
                              co%x_Mx=FF%Punto%PhysCoor(iEz)%x(i); co%y_Mx=FF%Punto%PhysCoor(iEz)%y(j); co%z_Mx=FF%Punto%PhysCoor(iEz)%z(k)
                              co%x_Mz=FF%Punto%PhysCoor(iEx)%x(i); co%y_Mz=FF%Punto%PhysCoor(iEx)%y(j); co%z_Mz=FF%Punto%PhysCoor(iEx)%z(k)
                              co%x_Jz=co%x_Mx;                     co%y_Jz=co%y_Mx;                     co%z_Jz=co%z_Mx;
                              co%x_Jx=co%x_Mz;                     co%y_Jx=co%y_Mz;                     co%z_Jx=co%z_Mz;
                              call update_LN(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi)
                              !! simetrias
                              !!!la IzDe hay que llamarla para cada caso
                              new_Mx = Mx
                              new_Mz = Mz
                              new_Jx = Jx
                              new_Jz = Jz
                              new_co = co
                              call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              !
                              If( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN) then
                                 new_Mx = + Mx
                                 new_Mz = - Mz
                                 new_Jx = - Jx
                                 new_Jz = + Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%z_Jx= new_co%z_Mz;
                                 new_co%z_Jz= new_co%z_Mx;
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldIz_ClonePMC_DOWN .or. FF%farfieldDe_ClonePMC_DOWN) then
                                 new_Mx = - Mx
                                 new_Mz = + Mz
                                 new_Jx = + Jx
                                 new_Jz = - Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%z_Jx= new_co%z_Mz;
                                 new_co%z_Jz= new_co%z_Mx;
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePEC_UP) then
                                 new_Mx = + Mx
                                 new_Mz = - Mz
                                 new_Jx = - Jx
                                 new_Jz = + Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%z_Jx= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePMC_UP) then
                                 new_Mx = - Mx
                                 new_Mz = + Mz
                                 new_Jx = + Jx
                                 new_Jz = - Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%z_Jx= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK) then
                                 new_Mx = - Mx
                                 new_Mz = + Mz
                                 new_Jx = + Jx
                                 new_Jz = - Jz
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetMinus
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK) then
                                 new_Mx = + Mx
                                 new_Mz = - Mz
                                 new_Jx = - Jx
                                 new_Jz = + Jz
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetMinus
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT) then
                                 new_Mx = - Mx
                                 new_Mz = + Mz
                                 new_Jx = + Jx
                                 new_Jz = - Jz
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetPlus
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT) then
                                 new_Mx = + Mx
                                 new_Mz = - Mz
                                 new_Jx = - Jx
                                 new_Jz = + Jz
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetPlus
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              !!!!!!!!!CASOS MIXTOS esquinas
                              If ((( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN).and.( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK)).or. &
                              (( FF%farfieldIz_ClonePMC_DOWN.or.FF%farfieldDe_ClonePMC_DOWN).and.( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK)).or. &
                              (( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN).and.( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK)).or. &
                              (( FF%farfieldIz_ClonePMC_DOWN.or.FF%farfieldDe_ClonePMC_DOWN).and.( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK)) )  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN).and.( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK)).or. &
                                 (( FF%farfieldIz_ClonePMC_DOWN.or.FF%farfieldDe_ClonePMC_DOWN).and.( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK)) )  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_Mz = signo * Mz
                                 new_Jx = signo * Jx
                                 new_Jz = signo * Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetMinus
                                 new_co%z_Jx= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_Mx
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              If ((( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN).and.( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT)).or. &
                              (( FF%farfieldIz_ClonePMC_DOWN.or.FF%farfieldDe_ClonePMC_DOWN).and.( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT)).or. &
                              (( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN).and.( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT)).or. &
                              (( FF%farfieldIz_ClonePMC_DOWN.or.FF%farfieldDe_ClonePMC_DOWN).and.( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT)) ) then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldIz_ClonePEC_DOWN.or.FF%farfieldDe_ClonePEC_DOWN).and.( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT)).or. &
                                 (( FF%farfieldIz_ClonePMC_DOWN.or.FF%farfieldDe_ClonePMC_DOWN).and.( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT)) )   sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_Mz = signo * Mz
                                 new_Jx = signo * Jx
                                 new_Jz = signo * Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetMinus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetMinus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetPlus
                                 new_co%z_Jx= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_Mx
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF

                              If ((( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePEC_UP).and.( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK)).or. &
                              (( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePMC_UP).and.( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK)).or. &
                              (( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePEC_UP).and.( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK)).or. &
                              (( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePMC_UP).and.( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK)) )  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePEC_UP).and.( FF%farfieldIz_ClonePEC_BACK.or.FF%farfieldDe_ClonePEC_BACK)).or. &
                                 (( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePMC_UP).and.( FF%farfieldIz_ClonePMC_BACK.or.FF%farfieldDe_ClonePMC_BACK)) )  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_Mz = signo * Mz
                                 new_Jx = signo * Jx
                                 new_Jz = signo * Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetMinus
                                 new_co%z_Jx= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_Mx
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              If ((( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePEC_UP).and.( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT)).or. &
                              (( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePMC_UP).and.( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT)).or. &
                              (( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePMC_UP).and.( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT)).or. &
                              (( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePEC_UP).and.( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT)))  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldIz_ClonePEC_UP.or.FF%farfieldDe_ClonePEC_UP).and.( FF%farfieldIz_ClonePEC_FRONT.or.FF%farfieldDe_ClonePEC_FRONT)).or. &
                                 (( FF%farfieldIz_ClonePMC_UP.or.FF%farfieldDe_ClonePMC_UP).and.( FF%farfieldIz_ClonePMC_FRONT.or.FF%farfieldDe_ClonePMC_FRONT)))  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_Mz = signo * Mz
                                 new_Jx = signo * Jx
                                 new_Jz = signo * Jz
                                 new_co%z_Mx=     -co%z_Mx +FF%ZOffsetPlus
                                 new_co%z_Mz=     -co%z_Mz +FF%ZOffsetPlus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_Mz=     -co%x_Mz +FF%xOffsetPlus
                                 new_co%z_Jx= new_co%z_Mz
                                 new_co%z_Jz= new_co%z_Mx
                                 new_co%x_Jx= new_co%x_Mz
                                 new_co%x_Jz= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneIzDe(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !!! fin simetrias
                           enddo
                        enddo
#ifdef CompileWithOpenMP
!$OMP                 END PARALLEL DO
!$OMP BARRIER
#endif
                        L_theta_final = L_theta_final + L_theta ; L_phi_final   = L_phi_final   + L_phi
                        N_theta_final = N_theta_final + N_theta ; N_phi_final   = N_phi_final   + N_phi
                        L_theta=0.0_RKIND ; L_phi=0.0_RKIND ; N_theta=0.0_RKIND ; N_phi=0.0_RKIND ;
                     endif !del goAhead
                  end do !del donde

                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                  !--->
                  !Down Up
                  L_theta=0.0_RKIND ; L_phi=0.0_RKIND ; N_theta=0.0_RKIND ; N_phi=0.0_RKIND ;
                  do donde=1,2
                     co%x_Mx=0;co%y_Mx=0;co%z_Mx=0;
                     co%x_My=0;co%y_My=0;co%z_My=0;
                     co%x_Mz=0;co%y_Mz=0;co%z_Mz=0;
                     co%x_Jx=0;co%y_Jx=0;co%z_Jx=0;
                     co%x_Jy=0;co%y_Jy=0;co%z_Jy=0;
                     co%x_Jz=0;co%y_Jz=0;co%z_Jz=0;
                     if (donde==1) then
                        k = FF%AbAr%K%aba%Ey
                        normal=-1.0_RKIND
                        GOahead = ( FF%farfieldAb .and. facesNF2FF%Ab)
                        EcampoY  => FF%EyAb
                        EcampoX  => FF%ExAb
                        HcampoY  => FF%HyAb
                        HcampoX  => FF%HxAb
                        Hcampo2Y => FF%HyAb2
                        Hcampo2X => FF%HxAb2
                     else
                        k = FF%AbAr%K%arr%Ey
                        normal=+1.0_RKIND
                        GOahead = ( FF%farfieldAr .and. facesNF2FF%Ar)
                        EcampoY  => FF%EyAr
                        EcampoX  => FF%ExAr
                        HcampoY  => FF%HyAr
                        HcampoX  => FF%HxAr
                        Hcampo2Y => FF%HyAr2
                        Hcampo2X => FF%HxAr2
                     endif
                     if (GOahead) then
#ifdef CompileWithOpenMP
!$OMP                 PARALLEL DO DEFAULT(SHARED) private (co,new_co,Mx,My,Mz,Jx,Jy,Jz,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,i,j,i_m,j_m) REDUCTION(+:L_theta) REDUCTION(+:L_phi) REDUCTION(+:N_theta) REDUCTION(+:N_phi)
#endif
                        do j = FF%AbAr%J%com%Ey,FF%AbAr%J%fin%Ex
                           j_m = j - b%Ey%YI
                           do i = FF%AbAr%I%com%Ey, FF%AbAr%I%fin%Ey
                              i_m = i - b%Ey%XI
                              Mx=0.0_RKIND; My=0.0_RKIND; Mz=0.0_RKIND; Jx=0.0_RKIND; Jy=0.0_RKIND; Jz=0.0_RKIND
                              if (j.le.FF%AbAr%J%fin%Ey)  Mx = - EcampoY( i_m, j_m,ii) *dxh(i_m)*dye(j_m)*NORMAL
                              if (i.le.FF%AbAr%I%fin%Ex)  My = + EcampoX( i_m, j_m,ii) *dxe(i_m)*dyh(j_m)*NORMAL
                              if (j.le.FF%AbAr%J%fin%Hy)  Jx = + (Average(pasadas, HcampoY( i_m, j_m,ii) , Hcampo2Y( i_m, j_m,ii))) *dxe(i_m)*dyh(j_m)*NORMAL
                              if (i.le.FF%AbAr%I%fin%Hx)  Jy = - (Average(pasadas, HcampoX( i_m, j_m,ii) , Hcampo2X( i_m, j_m,ii))) *dxh(i_m)*dye(j_m)*NORMAL
                              co%x_Mx=FF%Punto%PhysCoor(iEy)%x(i); co%y_Mx=FF%Punto%PhysCoor(iEy)%y(j); co%z_Mx=FF%Punto%PhysCoor(iEy)%z(k)
                              co%x_My=FF%Punto%PhysCoor(iEx)%x(i); co%y_My=FF%Punto%PhysCoor(iEx)%y(j); co%z_My=FF%Punto%PhysCoor(iEx)%z(k)
                              co%x_Jx=co%x_My;                     co%y_Jx=co%y_My;                     co%z_Jx=co%z_My;
                              co%x_Jy=co%x_Mx;                     co%y_Jy=co%y_Mx;                     co%z_Jy=co%z_Mx;
                              call update_LN(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi)!! simetrias
                              !!!!!!simetrias
                              !!!la AbAr hay que llamarla para cada caso
                              new_Mx = Mx
                              new_My = My
                              new_Jx = Jx
                              new_Jy = Jy
                              new_co = co
                              call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              !
                              If( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT) then
                                 !!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?!?
                                 new_Mx = + Mx
                                 new_My = - My
                                 new_Jx = - Jx
                                 new_Jy = + Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetMinus
                                 new_co%y_My=     -co%y_My +FF%yOffsetMinus
                                 new_co%y_Jx= new_co%y_My;
                                 new_co%y_Jy= new_co%y_Mx;
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldAb_ClonePMC_LEFT .or. FF%farfieldAr_ClonePMC_LEFT) then
                                 new_Mx = - Mx
                                 new_My = + My
                                 new_Jx = + Jx
                                 new_Jy = - Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetMinus
                                 new_co%y_My=     -co%y_My +FF%yOffsetMinus
                                 new_co%y_Jx= new_co%y_My;
                                 new_co%y_Jy= new_co%y_Mx;
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT) then
                                 new_Mx = + Mx
                                 new_My = - My
                                 new_Jx = - Jx
                                 new_Jy = + Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetPlus
                                 new_co%y_My=     -co%y_My +FF%yOffsetPlus
                                 new_co%y_Jx= new_co%y_My
                                 new_co%y_Jy= new_co%y_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT) then
                                 new_Mx = - Mx
                                 new_My = + My
                                 new_Jx = + Jx
                                 new_Jy = - Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetPlus
                                 new_co%y_My=     -co%y_My +FF%yOffsetPlus
                                 new_co%y_Jx= new_co%y_My
                                 new_co%y_Jy= new_co%y_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK) then
                                 new_Mx = - Mx
                                 new_My = + My
                                 new_Jx = + Jx
                                 new_Jy = - Jy
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_My=     -co%x_My +FF%xOffsetMinus
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK) then
                                 new_Mx = + Mx
                                 new_My = - My
                                 new_Jx = - Jx
                                 new_Jy = + Jy
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_My=     -co%x_My +FF%xOffsetMinus
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !
                              If( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT) then
                                 new_Mx = - Mx
                                 new_My = + My
                                 new_Jx = + Jx
                                 new_Jy = - Jy
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_My=     -co%x_My +FF%xOffsetPlus
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF
                              IF  ( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT) then
                                 new_Mx = + Mx
                                 new_My = - My
                                 new_Jx = - Jx
                                 new_Jy = + Jy
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_My=     -co%x_My +FF%xOffsetPlus
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              !!!!!!!!!CASOS MIXTOS esquinas
                              If ((( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT).and.( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK)).or. &
                              (( FF%farfieldAb_ClonePMC_LEFT.or.FF%farfieldAr_ClonePMC_LEFT).and.( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK)).or. &
                              (( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT).and.( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK)).or. &
                              (( FF%farfieldAb_ClonePMC_LEFT.or.FF%farfieldAr_ClonePMC_LEFT).and.( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK)) )  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT).and.( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK)).or. &
                                 (( FF%farfieldAb_ClonePMC_LEFT.or.FF%farfieldAr_ClonePMC_LEFT).and.( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK)) )  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_My = signo * My
                                 new_Jx = signo * Jx
                                 new_Jy = signo * Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetMinus
                                 new_co%y_My=     -co%y_My +FF%yOffsetMinus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_My=     -co%x_My +FF%xOffsetMinus
                                 new_co%y_Jx= new_co%y_My
                                 new_co%y_Jy= new_co%y_Mx
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              If ((( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT).and.( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT)).or. &
                              (( FF%farfieldAb_ClonePMC_LEFT.or.FF%farfieldAr_ClonePMC_LEFT).and.( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT)).or. &
                              (( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT).and.( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT)).or. &
                              (( FF%farfieldAb_ClonePMC_LEFT.or.FF%farfieldAr_ClonePMC_LEFT).and.( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT)) ) then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldAb_ClonePEC_LEFT.or.FF%farfieldAr_ClonePEC_LEFT).and.( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT)).or. &
                                 (( FF%farfieldAb_ClonePMC_LEFT.or.FF%farfieldAr_ClonePMC_LEFT).and.( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT)) )  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_My = signo * My
                                 new_Jx = signo * Jx
                                 new_Jy = signo * Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetMinus
                                 new_co%y_My=     -co%y_My +FF%yOffsetMinus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_My=     -co%x_My +FF%xOffsetPlus
                                 new_co%y_Jx= new_co%y_My
                                 new_co%y_Jy= new_co%y_Mx
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              ENDIF

                              If ((( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT).and.( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK)).or. &
                              (( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT).and.( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK)).or. &
                              (( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT).and.( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK)).or. &
                              (( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT).and.( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK)) )  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT).and.( FF%farfieldAb_ClonePEC_BACK.or.FF%farfieldAr_ClonePEC_BACK)).or. &
                                 (( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT).and.( FF%farfieldAb_ClonePMC_BACK.or.FF%farfieldAr_ClonePMC_BACK)) )  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_My = signo * My
                                 new_Jx = signo * Jx
                                 new_Jy = signo * Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetPlus
                                 new_co%y_My=     -co%y_My +FF%yOffsetPlus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetMinus
                                 new_co%x_My=     -co%x_My +FF%xOffsetMinus
                                 new_co%y_Jx= new_co%y_My
                                 new_co%y_Jy= new_co%y_Mx
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif

                              If ((( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT).and.( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT)).or. &
                              (( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT).and.( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT)).or. &
                              (( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT).and.( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT)).or. &
                              (( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT).and.( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT)))  then
                                 sigNo=+1.0_RKIND
                                 If ((( FF%farfieldAb_ClonePEC_RIGHT.or.FF%farfieldAr_ClonePEC_RIGHT).and.( FF%farfieldAb_ClonePEC_FRONT.or.FF%farfieldAr_ClonePEC_FRONT)).or. &
                                 (( FF%farfieldAb_ClonePMC_RIGHT.or.FF%farfieldAr_ClonePMC_RIGHT).and.( FF%farfieldAb_ClonePMC_FRONT.or.FF%farfieldAr_ClonePMC_FRONT)))  sigNo=-1.0_RKIND
                                 new_Mx = signo * Mx
                                 new_My = signo * My
                                 new_Jx = signo * Jx
                                 new_Jy = signo * Jy
                                 new_co%y_Mx=     -co%y_Mx +FF%yOffsetPlus
                                 new_co%y_My=     -co%y_My +FF%yOffsetPlus
                                 new_co%x_Mx=     -co%x_Mx +FF%xOffsetPlus
                                 new_co%x_My=     -co%x_My +FF%xOffsetPlus
                                 new_co%y_Jx= new_co%y_My
                                 new_co%y_Jy= new_co%y_Mx
                                 new_co%x_Jx= new_co%x_My
                                 new_co%x_Jy= new_co%x_Mx
                                 call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
                                 call cloneAbAr(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
                              endif
                              !!! fin simetrias
                           enddo
                        enddo
#ifdef CompileWithOpenMP
!$OMP                 END PARALLEL DO
!$OMP BARRIER
#endif
                        L_theta_final = L_theta_final + L_theta ; L_phi_final   = L_phi_final   + L_phi
                        N_theta_final = N_theta_final + N_theta ; N_phi_final   = N_phi_final   + N_phi
                        L_theta=0.0_RKIND ; L_phi=0.0_RKIND ; N_theta=0.0_RKIND ; N_phi=0.0_RKIND ;
                     endif !del goAhead
                  end do !del donde


                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef CompileWithMPI
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  dummy=real(L_theta_final)
                  call MPI_AllReduce(dummy, newdummy1, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  dummy=AIMAG(L_theta_final)
                  call MPI_AllReduce(dummy, newdummy2, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)

                  L_theta_final= newdummy1+(0.0_RKIND,1.0_RKIND)*newdummy2
                  !
                  dummy=real(L_phi_final)
                  call MPI_AllReduce(dummy, newdummy1, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  dummy=AIMAG(L_phi_final)
                  call MPI_AllReduce(dummy, newdummy2, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  !
                  L_phi_final= newdummy1+(0.0_RKIND,1.0_RKIND)*newdummy2

                  dummy=real(N_theta_final)
                  call MPI_AllReduce(dummy, newdummy1, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  dummy=AIMAG(N_theta_final)
                  call MPI_AllReduce(dummy, newdummy2, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)

                  N_theta_final= newdummy1+(0.0_RKIND,1.0_RKIND)*newdummy2
                  !
                  dummy=real(N_phi_final)
                  call MPI_AllReduce(dummy, newdummy1, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  dummy=AIMAG(N_phi_final)
                  call MPI_AllReduce(dummy, newdummy2, 1_4, REALSIZE, MPI_SUM, FF%MPISubComm, ierr)
                  call MPI_Barrier(FF%MPISubComm,ierr)
                  !
                  N_phi_final= newdummy1+(0.0_RKIND,1.0_RKIND)*newdummy2
                  call MPI_Barrier(FF%MPISubComm,ierr)
#endif
                  Etheta(pasadas) = -(0,1.0_RKIND)*freq/(2.0_RKIND * cluz)*(L_phi_final + zvac * N_theta_final) !/FF%dftEntrada(ii) !no normalizar para calcular bien potencia
                  Ephi(pasadas)   =  (0,1.0_RKIND)*freq/(2.0_RKIND * cluz)*(L_theta_final - zvac * N_phi_final) !/FF%dftEntrada(ii) !no normalizar para calcular bien potencia
                  RCS(pasadas)    =  (2.0_RKIND * pi*freq/cluz)**2.0_RKIND  / (4.0_RKIND * pi*Abs(FF%dftEntrada(ii))**2.0_RKIND ) * &
                  (abs(L_phi_final + zvac * N_theta_final)**2.0_RKIND + abs(L_theta_final - zvac * N_phi_final)**2.0_RKIND )


#ifdef CompileWithMPI
                  if (FF%MPIRoot == layoutnumber)  then
#endif
                  if (pasadas==1) write(FF%unitfarfield,fmt) freq,theta,phi,&
                  abs(Etheta(2)),ATAN2( AIMAG( Etheta(2)) , REAL( Etheta(2) ) ), & !!! PASADAS=2=GEOMETRICA,, PASADAS=1=ARITMETICA
                  abs(Ephi(2)) , ATAN2( AIMAG( Ephi(2)  ) , REAL( Ephi(2)   ) ), RCS(1),RCS(2)

#ifdef CompileWithMPI
               endif
#endif
               end do !de pasadas
            end do ! del bucle en phi
         end do !del bucle en theta
!
         write(dubuf,'(a,e19.9e3)')  ' NF2FF: End processing freq= ',freq
         ! call print11(layoutnumber,dubuf,.TRUE.)
!
      end do !del de frecuencias
      !fin calculo
      if (allocated(Thetavector)) deallocate (Thetavector)
      if (allocated(sizephi)) deallocate (sizephi)
      if (allocated(Phimatrix)) deallocate (Phimatrix)

      CLOSE (FF%unitfarfield)



      write(dubuf,'(a)')  ' NF2FF: END '
      if (layoutnumber == 0) call print11(layoutnumber,dubuf,.TRUE.)

   end subroutine


   subroutine update_LN(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi)

      type (co_t) :: co
      complex (KIND=CKIND) :: L_theta,L_phi,Mx,My,Mz,comunMx,comunMy,comunMz,comunJx,comunJy,comunJz,comun
      complex (KIND=CKIND) :: N_theta,N_phi,Jx,Jy,Jz
      real (KIND=RKIND) :: sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi
      !!puede que me falte un signo en la exponencial aunque no afecta a nada (creo que es -exp[-j k r] taflove 3rd 372 nf2ff) 02/03/15
      comunMx=exp(comun*(co%x_Mx*sintheta_cosphi + co%y_Mx*sintheta_sinphi + co%z_Mx*costheta))
      comunMy=exp(comun*(co%x_My*sintheta_cosphi + co%y_My*sintheta_sinphi + co%z_My*costheta))
      comunMz=exp(comun*(co%x_Mz*sintheta_cosphi + co%y_Mz*sintheta_sinphi + co%z_Mz*costheta))
      comunJx=exp(comun*(co%x_Jx*sintheta_cosphi + co%y_Jx*sintheta_sinphi + co%z_Jx*costheta))
      comunJy=exp(comun*(co%x_Jy*sintheta_cosphi + co%y_Jy*sintheta_sinphi + co%z_Jy*costheta))
      comunJz=exp(comun*(co%x_Jz*sintheta_cosphi + co%y_Jz*sintheta_sinphi + co%z_Jz*costheta))
      !
      L_theta=L_theta + (  Mx * costheta_cosphi * comunMx +   My * costheta_sinphi * comunMy -  Mz  * sintheta * comunMz)
      L_phi = L_phi   + ( -Mx * sinphi * comunMx   +   My * cosphi * comunMy)
      !
      N_theta=N_theta + (  Jx * costheta_cosphi * comunJx +   Jy * costheta_sinphi * comunJy -  Jz  *sintheta * comunJz)
      N_phi = N_phi   + ( -Jx * sinphi * comunJx   +   Jy * cosphi * comunJy)
   end subroutine



   subroutine cloneTrFr(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
      type (co_t) :: co,new_co
      complex( kind = CKIND) :: new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz
      complex (KIND=CKIND) :: L_theta,L_phi,Mx,My,Mz,comun
      complex (KIND=CKIND) :: N_theta,N_phi,Jx,Jy,Jz
      real (KIND=RKIND) :: sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,NORMAL

      new_co=co; new_Mx=Mx; new_My=My; new_Mz=Mz; new_Jx=Jx; new_Jy=Jy; new_Jz=Jz;
      If( FF%farfieldTr_ClonePEC_Front.or.FF%farfieldFr_ClonePEC_Back) then
         new_My = + My  !solo en este caso cambian las normales
         new_Mz = + Mz
         new_Jy = - Jy
         new_Jz = - Jz
         new_co%x_My=     co%x_My + FF%XDobleAncho*NORMAL !cambio de signo resto o sumo distancia
         new_co%x_Mz=     co%x_Mz + FF%XDobleAncho*NORMAL
         new_co%x_Jy= new_co%x_Mz
         new_co%x_Jz= new_co%x_My
         call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
      ENDIF
      IF ( FF%farfieldTr_ClonePMC_Front.or.FF%farfieldFr_ClonePMC_Back) then
         new_My = - My
         new_Mz = - Mz
         new_Jy = + Jy
         new_Jz = + Jz
         new_co%x_My=     co%x_My + FF%XDobleAncho*NORMAL !cambio de signo resto o sumo distancia
         new_co%x_Mz=     co%x_Mz + FF%XDobleAncho*NORMAL
         new_co%x_Jy= new_co%x_Mz
         new_co%x_Jz= new_co%x_My
         call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
      endif
      return
   end subroutine



   subroutine cloneIzDe(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
      type (co_t) :: co,new_co
      complex( kind = CKIND) :: new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz
      complex (KIND=CKIND) :: L_theta,L_phi,Mx,My,Mz,comun
      complex (KIND=CKIND) :: N_theta,N_phi,Jx,Jy,Jz
      real (KIND=RKIND) :: sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,NORMAL

      new_co=co; new_Mx=Mx; new_My=My; new_Mz=Mz; new_Jx=Jx; new_Jy=Jy; new_Jz=Jz;
      If( FF%farfieldIz_ClonePEC_Right.or.FF%farfieldDe_ClonePEC_Left) then
         new_Mx = + Mx !solo en este caso cambian las normales
         new_Mz = + Mz
         new_Jx = - Jx
         new_Jz = - Jz
         new_co%y_Mx=     co%y_Mx + FF%YDobleAncho*NORMAL !cambio de signo resto o sumo distancia
         new_co%y_Mz=     co%y_Mz + FF%YDobleAncho*NORMAL
         new_co%y_Jx= new_co%y_Mz
         new_co%y_Jz= new_co%y_Mx
         call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
      ENDIF
      IF ( FF%farfieldIz_ClonePMC_Right.or.FF%farfieldDe_ClonePMC_Left) then
         new_Mx = - Mx !solo en este caso cambian las normales
         new_Mz = - Mz
         new_Jx = + Jx
         new_Jz = + Jz
         new_co%y_Mx=     co%y_Mx + FF%YDobleAncho*NORMAL !cambio de signo resto o sumo distancia
         new_co%y_Mz=     co%y_Mz + FF%YDobleAncho*NORMAL
         new_co%y_Jx= new_co%y_Mz
         new_co%y_Jz= new_co%y_Mx
         call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
      endif
      return
   end subroutine


   subroutine cloneAbAr(comun,co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,Mx,My,Mz,Jx,Jy,Jz,L_theta,L_phi,N_theta,N_phi,NORMAL)
      type (co_t) :: co,new_co
      complex( kind = CKIND) :: new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz
      complex (KIND=CKIND) :: L_theta,L_phi,Mx,My,Mz,comun
      complex (KIND=CKIND) :: N_theta,N_phi,Jx,Jy,Jz
      real (KIND=RKIND) :: sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,NORMAL

      new_co=co; new_Mx=Mx; new_My=My; new_Mz=Mz; new_Jx=Jx; new_Jy=Jy; new_Jz=Jz;
      If( FF%farfieldAb_ClonePEC_UP.or.FF%farfieldAr_ClonePEC_DOWN) then
         new_Mx = + Mx !solo en este caso cambian las normales
         new_My = + My
         new_Jx = - Jx
         new_Jy = - Jy
         new_co%Z_Mx=     co%Z_Mx + FF%ZDobleAncho*NORMAL !cambio de signo resto o sumo distancia
         new_co%Z_My=     co%Z_My + FF%ZDobleAncho*NORMAL
         new_co%Z_Jx= new_co%Z_My
         new_co%Z_Jy= new_co%Z_Mx
         call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
      ENDIF
      IF ( FF%farfieldAb_ClonePMC_UP.or.FF%farfieldAr_ClonePMC_DOWN) then
         new_Mx = - Mx !solo en este caso cambian las normales
         new_My = - My
         new_Jx = + Jx
         new_Jy = + Jy
         new_co%Z_Mx=     co%Z_Mx + FF%ZDobleAncho*NORMAL !cambio de signo resto o sumo distancia
         new_co%Z_My=     co%Z_My + FF%ZDobleAncho*NORMAL
         new_co%Z_Jx= new_co%Z_My
         new_co%Z_Jy= new_co%Z_Mx
         call update_LN(comun,new_co,sintheta_cosphi,sintheta_sinphi,costheta,costheta_cosphi,costheta_sinphi,sintheta,sinphi,cosphi,new_Mx,new_My,new_Mz,new_Jx,new_Jy,new_Jz,L_theta,L_phi,N_theta,N_phi)
      endif
      return
   end subroutine


   !!!!!!!COMPUTE THE GEOMETRIC AVERAGE OF TWO COMPLEX NUMBERS
   function average (pasadas,z1,z2) result (z)
      complex (KIND=CKIND) :: Z1,Z2,Z
      real  (KIND=RKIND) :: phi1,phi2,nphi1,nphi2
      integer(4) :: pasadas
      Z=(0.0_RKIND,0.0_RKIND)
      if (pasadas ==2 ) then !geometrica
         phi1=ATAN2(AIMAG(Z1),REAL(Z1))
         phi2=ATAN2(AIMAG(Z2),REAL(Z2))


         !TRAMPA CHINA PARA EVITAR LOS BRANCH CUT EN 0 Y PI
         nphi1 =phi1
         nphi2 =phi2
         if ((phi1 < -pi/2.0_RKIND).AND.(PHI2 > PI/2.0_RKIND)) nphi1 =phi1 -2.0_RKIND * pi
         if ((phi2 < -pi/2.0_RKIND).AND.(PHI1 > PI/2.0_RKIND)) nphi2 =phi2 -2.0_RKIND * pi

         Z=SQRT(ABS(Z1*Z2)) * EXP((0.0_RKIND,1.0_RKIND)*(nPHI1+nPHI2)/2.0_RKIND)
      elseif (pasadas==1) then !aritmetica
         Z=(z1+z2)/2.0_RKIND
      endif

      !!if (abs(z)/abs(zarit) > 1e1) write(3555,'(e18.6e3,a,2e18.6e3,a,2e18.6e3)') abs(z)/abs(zarit),'---> ',abs(z),abs(zarit),' -- ',atan2(AIMAG(z),real(z)),atan2(AIMAG(zarit),real(zarit))

      RETURN
   END FUNCTION

#endif
END MODULE farfield_m
