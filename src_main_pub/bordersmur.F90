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
!  Borders :  MUR  handling
!  Creation date Date :  January, 8, 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BORDERS_MUR
   use fdetypes
   use Report
   implicit none
   private
   !
   !
   ! Limits of the MUR region
   type XYZlimit_tvar
      integer (kind=4), dimension(1:6)  ::  XI,XE,YI,YE,ZI,ZE
   end type XYZlimit_tvar
   type (XYZlimit_tvar), dimension(4:6)  ::    MURc



   type LR
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Past_Hx,Past_Hz,PastPast_Hx,PastPast_Hz
   end type
   type DU
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Past_Hy,Past_Hx,PastPast_Hy,PastPast_Hx
   end type
   type BF
      REAL (KIND=RKIND) , pointer, dimension ( : , : , : )  ::  Past_Hz,Past_Hy,PastPast_Hz,PastPast_Hy
   end type

   !!!LOCAL VARIABLES
   type (LR), dimension(left : right) , save ::  regLR
   type (DU), dimension(down : up)    , save ::  regDU
   type (BF), dimension(back : front) , save ::  regBF


   real (kind = RKIND), dimension(  :), allocatable, SAVE ::  back_CAB1, back_CAB3, back_cab4, &
   front_CAB1,front_CAB3,front_cab4, &
   left_CAB1, left_CAB3, left_cab4, &
   right_CAB1,right_CAB3,right_cab4, &
   down_CAB1, down_CAB3, down_cab4, &
   up_CAB1,   up_CAB3,   up_cab4
!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  cluz
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!
   !
   public  ::  InitMURBorders, AdvanceMagneticMUR,StoreFieldsMURBorders,DestroyMURBorders,calc_murconstants


contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes MUR data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMURBorders(sgg,ThereAreMURBorders,resume,Idxh,Idyh,Idzh,eps00,mu00)
      REAL (KIND=RKIND)           ::  eps00,mu00

      type (SGGFDTDINFO), intent(IN)         ::  sgg

      REAL (KIND=RKIND) , dimension (:)   , intent(in)      ::  &
      Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)
      !!!
      !
      logical  ::  ThereAreMURBorders,resume
      integer (kind=4)  ::  i,j,k,region,field,i1

      !character (len=BUFSIZE) :: donde
      !integer (KIND=4) :: layoutnumber
!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      cluz=1.0_RKIND/sqrt(eps0*mu0)
!


      !
      ThereAreMURBorders=.false.
      if (sgg%Border%IsBackMUR.or.sgg%Border%IsFrontMUR.or.sgg%Border%IsLeftMUR.or.sgg%Border%IsRightMUR.or. &
      sgg%Border%IsUpMUR.or.sgg%Border%IsDownMUR) ThereAreMURBorders=.true.
      IF (.not.(ThereAreMURBorders)) return

      allocate( back_CAB1( 0 :  sgg%NumMedia), back_CAB3( 0 :  sgg%NumMedia), back_cab4( 0 :  sgg%NumMedia), &
      front_CAB1( 0 :  sgg%NumMedia),front_CAB3( 0 :  sgg%NumMedia),front_cab4( 0 :  sgg%NumMedia), &
      left_CAB1( 0 :  sgg%NumMedia), left_CAB3( 0 :  sgg%NumMedia), left_cab4( 0 :  sgg%NumMedia), &
      right_CAB1( 0 :  sgg%NumMedia),right_CAB3( 0 :  sgg%NumMedia),right_cab4( 0 :  sgg%NumMedia), &
      down_CAB1( 0 :  sgg%NumMedia), down_CAB3( 0 :  sgg%NumMedia), down_cab4( 0 :  sgg%NumMedia), &
      up_CAB1( 0 :  sgg%NumMedia),   up_CAB3( 0 :  sgg%NumMedia),   up_cab4( 0 :  sgg%NumMedia) )
      !Find the limits of each of the 6 padding MUR regions for each field component


      do field=iHx,iHz
         !
         MURc(field)%XI(Down)  =                                sgg%Sweep(field)%XI
         MURc(field)%XE(Down)  =                                sgg%Sweep(field)%XE
         MURc(field)%YI(Down)  =                                sgg%Sweep(field)%YI
         MURc(field)%YE(Down)  =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(Down)  =                                sgg%Sweep(field)%ZI-1
         MURc(field)%ZE(Down)  = MURc(field)%ZI(Down) + 1
         !
         MURc(field)%XI(Up)    =                                sgg%Sweep(field)%XI
         MURc(field)%XE(Up)    =                                sgg%Sweep(field)%XE
         MURc(field)%YI(Up)    =                                sgg%Sweep(field)%YI
         MURc(field)%YE(Up)    =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(Up)    =                                sgg%Sweep(field)%ZE
         MURc(field)%ZE(Up)    = MURc(field)%ZI(Up) + 1
         !
         MURc(field)%XI(Left)  =                                sgg%Sweep(field)%XI
         MURc(field)%XE(Left)  =                                sgg%Sweep(field)%XE
         MURc(field)%YI(Left)  =                                sgg%Sweep(field)%YI-1
         MURc(field)%YE(Left)  = MURc(field)%YI(Left) + 1
         MURc(field)%ZI(Left)  =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(Left)  =                                sgg%Sweep(field)%ZE
         !
         MURc(field)%XI(Right) =                                sgg%Sweep(field)%XI
         MURc(field)%XE(Right) =                                sgg%Sweep(field)%XE
         MURc(field)%YI(Right) =                                sgg%Sweep(field)%YE
         MURc(field)%YE(Right) = MURc(field)%YI(Right) + 1
         MURc(field)%ZI(Right) =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(Right) =                                sgg%Sweep(field)%ZE
         !
         MURc(field)%XI(Back)  =                                sgg%Sweep(field)%XI-1
         MURc(field)%XE(Back)  = MURc(field)%XI(Back) + 1
         MURc(field)%YI(Back)  =                                sgg%Sweep(field)%YI
         MURc(field)%YE(Back)  =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(Back)  =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(Back)  =                                sgg%Sweep(field)%ZE
         !
         MURc(field)%XI(Front) =                                sgg%Sweep(field)%XE
         MURc(field)%XE(Front) = MURc(field)%XI(Front) + 1
         MURc(field)%YI(Front) =                                sgg%Sweep(field)%YI
         MURc(field)%YE(Front) =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(Front) =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(Front) =                                sgg%Sweep(field)%ZE
         !
      end do

      !Fake coms and ends IN CASE OF NO MUR SO THAT NEVER ENTER THE DO FOR THESE CASES
      IF (.not.(sgg%Border%IsDownMUR)) MURc(4:6)%ZI(down)=MURc(4:6)%ZE(down)+100
      IF (.not.(sgg%Border%IsUpMUR))   MURc(4:6)%ZI(up)  =MURc(4:6)%ZE(up)  +100
      !
      IF (.not.(sgg%Border%IsLeftMUR))  MURc(4:6)%ZI(left) =MURc(4:6)%ZE(left) +100
      IF (.not.(sgg%Border%IsRightMUR)) MURc(4:6)%ZI(right)=MURc(4:6)%ZE(right)+100
      !
      IF (.not.(sgg%Border%IsFrontMUR)) MURc(4:6)%ZI(front)=MURc(4:6)%ZE(front)+100
      IF (.not.(sgg%Border%IsBackMUR))  MURc(4:6)%ZI(back) =MURc(4:6)%ZE(back) +100

      !MUR Field component matrix allocation
      do REGION =left,right
         allocate (regLR(region)%Past_Hx(MURc(iHx)%XI(region) : MURc(iHx)%XE(region), &
         MURc(iHx)%YI(region) : MURc(iHx)%YE(region), &
         MURc(iHx)%ZI(region) : MURc(iHx)%ZE(region)),&
         regLR(region)%Past_Hz(MURc(iHz)%XI(region) : MURc(iHz)%XE(region), &
         MURc(iHz)%YI(region) : MURc(iHz)%YE(region), &
         MURc(iHz)%ZI(region) : MURc(iHz)%ZE(region)))
         if (.not.resume) then
            regLR(REGION)%Past_Hx=0.0_RKIND ; regLR(REGION)%Past_Hz=0.0_RKIND ;
         else
            Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
               Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
                  READ (14) (regLR(region)%Past_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
               End do
            End do
            Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
               Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
                  READ (14) ( regLR(region)%Past_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
               End do
            End do
         endif
      end do
      do REGION =down,up
         allocate (regDU(region)%Past_Hy(MURc(iHy)%XI(region) : MURc(iHy)%XE(region), &
         MURc(iHy)%YI(region) : MURc(iHy)%YE(region), &
         MURc(iHy)%ZI(region) : MURc(iHy)%ZE(region)),&
         regDU(region)%Past_Hx(MURc(iHx)%XI(region) : MURc(iHx)%XE(region), &
         MURc(iHx)%YI(region) : MURc(iHx)%YE(region), &
         MURc(iHx)%ZI(region) : MURc(iHx)%ZE(region)))
         if (.not.resume) then
            regDU(REGION)%Past_Hy=0.0_RKIND ; regDU(REGION)%Past_Hx=0.0_RKIND ;
         else
            Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
               Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
                  READ (14) (regDU(region)%Past_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
               End do
            End do
            Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
               Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
                  READ (14) ( regDU(region)%Past_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
               End do
            End do
         endif
      end do
      do REGION =back,front
         allocate (regBF(region)%Past_Hz(MURc(iHz)%XI(region) : MURc(iHz)%XE(region), &
         MURc(iHz)%YI(region) : MURc(iHz)%YE(region), &
         MURc(iHz)%ZI(region) : MURc(iHz)%ZE(region)),&
         regBF(region)%Past_Hy(MURc(iHy)%XI(region) : MURc(iHy)%XE(region), &
         MURc(iHy)%YI(region) : MURc(iHy)%YE(region), &
         MURc(iHy)%ZI(region) : MURc(iHy)%ZE(region)))
         if (.not.resume) then
            regBF(REGION)%Past_Hz=0.0_RKIND ; regBF(REGION)%Past_Hy=0.0_RKIND ;
         else
            Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
               Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
                  READ (14) (regBF(region)%Past_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
               End do
            End do
            Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
               Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
                  READ (14) ( regBF(region)%Past_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
               End do
            End do
         endif
      end do

      !past past

      !MUR Field component matrix allocation
      do REGION =left,right
         allocate (regLR(region)%PastPast_Hx(MURc(iHx)%XI(region) : MURc(iHx)%XE(region), &
         MURc(iHx)%YI(region) : MURc(iHx)%YE(region), &
         MURc(iHx)%ZI(region) : MURc(iHx)%ZE(region)),&
         regLR(region)%PastPast_Hz(MURc(iHz)%XI(region) : MURc(iHz)%XE(region), &
         MURc(iHz)%YI(region) : MURc(iHz)%YE(region), &
         MURc(iHz)%ZI(region) : MURc(iHz)%ZE(region)))
         if (.not.resume) then
            regLR(REGION)%PastPast_Hx=0.0_RKIND ; regLR(REGION)%PastPast_Hz=0.0_RKIND ;
         else
            Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
               Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
                  READ (14) (regLR(region)%PastPast_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
               End do
            End do
            Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
               Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
                  READ (14) ( regLR(region)%PastPast_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
               End do
            End do
         endif
      end do
      do REGION =down,up
         allocate (regDU(region)%PastPast_Hy(MURc(iHy)%XI(region) : MURc(iHy)%XE(region), &
         MURc(iHy)%YI(region) : MURc(iHy)%YE(region), &
         MURc(iHy)%ZI(region) : MURc(iHy)%ZE(region)),&
         regDU(region)%PastPast_Hx(MURc(iHx)%XI(region) : MURc(iHx)%XE(region), &
         MURc(iHx)%YI(region) : MURc(iHx)%YE(region), &
         MURc(iHx)%ZI(region) : MURc(iHx)%ZE(region)))
         if (.not.resume) then
            regDU(REGION)%PastPast_Hy=0.0_RKIND ; regDU(REGION)%PastPast_Hx=0.0_RKIND ;
         else
            Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
               Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
                  READ (14) (regDU(region)%PastPast_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
               End do
            End do
            Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
               Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
                  READ (14) ( regDU(region)%PastPast_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
               End do
            End do
         endif
      end do
      do REGION =back,front
         allocate (regBF(region)%PastPast_Hz(MURc(iHz)%XI(region) : MURc(iHz)%XE(region), &
         MURc(iHz)%YI(region) : MURc(iHz)%YE(region), &
         MURc(iHz)%ZI(region) : MURc(iHz)%ZE(region)),&
         regBF(region)%PastPast_Hy(MURc(iHy)%XI(region) : MURc(iHy)%XE(region), &
         MURc(iHy)%YI(region) : MURc(iHy)%YE(region), &
         MURc(iHy)%ZI(region) : MURc(iHy)%ZE(region)))
         if (.not.resume) then
            regBF(REGION)%PastPast_Hz=0.0_RKIND ; regBF(REGION)%PastPast_Hy=0.0_RKIND ;
         else
            Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
               Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
                  READ (14) (regBF(region)%PastPast_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
               End do
            End do
            Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
               Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
                  READ (14) ( regBF(region)%PastPast_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
               End do
            End do
         endif
      end do

!!incializa constantes
      call calc_murconstants(sgg,Idxh,Idyh,Idzh,eps0,mu0)
 
      return
   end subroutine InitMURBorders

   subroutine calc_murconstants(sgg,Idxh,Idyh,Idzh,eps00,mu00)
        type (SGGFDTDINFO), intent(IN)         ::  sgg
        REAL (KIND=RKIND)           ::  eps00,mu00
        integer (kind=4)  ::  i,j,k,region,field,i1
        REAL (KIND=RKIND)   ::  cnum
        REAL (KIND=RKIND) , dimension (:)   , intent(in)      ::  &
        Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
        Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
        Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE)
!
        eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
        cluz=1.0_RKIND/sqrt(eps0*mu0)
!

        do i1=0,sgg%NumMedia
            !SE CREAN MAS DE LA CUENTA PERO LUEGO SE UTILIZAN SOLO LAS QUE SE NECESITEN
            cnum=(1.0_RKIND/Idxh(sgg%ALLOC(iEx)%XI))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            back_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            back_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            back_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
            !
            cnum=(1.0_RKIND/Idxh(sgg%ALLOC(iEx)%XE))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            front_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            front_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            front_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
            !!
            cnum=(1.0_RKIND/Idyh(sgg%ALLOC(iEy)%YI))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            left_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            left_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            left_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
            !
            cnum=(1.0_RKIND/Idyh(sgg%ALLOC(iEy)%YE))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            right_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            right_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            right_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
            !!
            cnum=(1.0_RKIND/Idzh(sgg%ALLOC(iEz)%ZI))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            down_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            down_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            down_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
            !
            cnum=(1.0_RKIND/Idzh(sgg%ALLOC(iEz)%ZE))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            up_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            up_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            up_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
        end do
        return
   end subroutine calc_murconstants


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Flush the MUR data to disk for resuming purposes
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine StoreFieldsMURBorders

      integer (kind=4)  ::  region,i,j,k


      do REGION =left,right
         Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
            Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
               write(14,err=634) (regLR(region)%Past_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
            End do
         End do
         Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
            Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
               write(14,err=634) ( regLR(region)%Past_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
            End do
         End do
      end do
      do REGION =down,up
         Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
            Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
               write(14,err=634) (regDU(region)%Past_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
            End do
         End do
         Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
            Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
               write(14,err=634) ( regDU(region)%Past_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
            End do
         End do
      end do
      do REGION =back,front
         Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
            Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
               write(14,err=634) (regBF(region)%Past_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
            End do
         End do
         Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
            Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
               write(14,err=634) ( regBF(region)%Past_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
            End do
         End do
      end do


      do REGION =left,right
         Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
            Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
               write(14,err=634) (regLR(region)%PastPast_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
            End do
         End do
         Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
            Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
               write(14,err=634) ( regLR(region)%PastPast_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
            End do
         End do
      end do
      do REGION =down,up
         Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
            Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
               write(14,err=634) (regDU(region)%PastPast_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
            End do
         End do
         Do k=MURc(iHx)%ZI(region),MURc(iHx)%ZE(region)
            Do j=MURc(iHx)%YI(region),MURc(iHx)%YE(region)
               write(14,err=634) ( regDU(region)%PastPast_Hx(i,j,k),i=MURc(iHx)%XI(region),MURc(iHx)%XE(region))
            End do
         End do
      end do
      do REGION =back,front
         Do k=MURc(iHz)%ZI(region),MURc(iHz)%ZE(region)
            Do j=MURc(iHz)%YI(region),MURc(iHz)%YE(region)
               write(14,err=634) (regBF(region)%PastPast_Hz(i,j,k),i=MURc(iHz)%XI(region),MURc(iHz)%XE(region))
            End do
         End do
         Do k=MURc(iHy)%ZI(region),MURc(iHy)%ZE(region)
            Do j=MURc(iHy)%YI(region),MURc(iHy)%YE(region)
               write(14,err=634) ( regBF(region)%PastPast_Hy(i,j,k),i=MURc(iHy)%XI(region),MURc(iHy)%XE(region))
            End do
         End do
      end do

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'BORDERSMUR: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsMURBorders


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyMURBorders
      integer (kind=4)  ::  region

      do REGION =left,right
         if (associated(regLR(region)%Past_Hx)) deallocate (regLR(region)%Past_Hx,regLR(region)%Past_Hz)
      end do
      do REGION =down,up
         if (associated(regDU(region)%Past_Hy)) deallocate (regDU(region)%Past_Hy,regDU(region)%Past_Hx)
      end do
      do REGION =back,front
         if (associated(regBF(region)%Past_Hz)) deallocate (regBF(region)%Past_Hz,regBF(region)%Past_Hy)
      end do


      do REGION =left,right
         if (associated(regLR(region)%PastPast_Hx)) deallocate (regLR(region)%PastPast_Hx,regLR(region)%PastPast_Hz)
      end do
      do REGION =down,up
         if (associated(regDU(region)%PastPast_Hy)) deallocate (regDU(region)%PastPast_Hy,regDU(region)%PastPast_Hx)
      end do
      do REGION =back,front
         if (associated(regBF(region)%PastPast_Hz)) deallocate (regBF(region)%PastPast_Hz,regBF(region)%PastPast_Hy)
      end do


      if (allocated(back_CAB1)) &
      deallocate(back_CAB1 ,back_CAB3 ,back_cab4 , &
      front_CAB1,front_CAB3,front_cab4, &
      left_CAB1 ,left_CAB3 ,left_cab4 , &
      right_CAB1,right_CAB3,right_cab4, &
      down_CAB1 ,down_CAB3 ,down_cab4 , &
      up_CAB1   ,up_CAB3   ,up_cab4    )

      return
   end subroutine DestroyMURBorders

   !**************************************************************************************************

   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Advances the magnetic field in the MUR
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine AdvanceMagneticMUR(b, sgg,sggMiHx, sggMiHy, sggMiHz, Hx, Hy, Hz,mur_second)
      !---------------------------> inputs <----------------------------------------------------------
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      type( bounds_t), intent( IN)  ::  b
      logical :: mur_second
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN)  ::  sggMiHx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN)  ::  sggMiHy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN)  ::  sggMiHz
      !--->
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer (kind=4)  ::  REGION, i, j, k, medio, i_m, j_m, k_m
      !---------------------------> empieza AdvanceMagneTicMUR <-------------------------------------

      !Hetic Fields MUR Zone
      !primero hay que updatear los edges porque las caras los utilizan



      if (mur_second) then


         call stoponerror(0,0,'ERROR: MUR SECOND not correctly implemented')
         !!!!!!?!?!?!?!?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!?!?!?!?!?!!!!!!!!!!!!!!!!!!!Edges!!!!!!!!!!!!!!!!!!!!!Mur primer orden
         !!!!!!?!?!?!?!?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsLeftMUR) then
            REGION = left
            j = MURc(iHx)%YI( REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHx)%ZI( REGION), MURc(iHx)%ZE( REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m + 1, k_m    )
                  Hx( i_m, j_m, k_m)=                                           + regLR( REGION)%Past_Hx(i    ,j + 1,k)          &
                  +  left_CAB1(medio)*(                    Hx(i_m  ,j_m + 1,k_m) - regLR( REGION)%Past_Hx(i    ,j    ,k))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(iHz)%YI( REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(iHz)%XI( REGION), MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz( i_m    , j_m + 1, k_m    )
                  Hz( i_m, j_m, k_m) =                                             + regLR( REGION)%Past_Hz(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                   Hz(i_m  ,j_m + 1,k_m) - regLR( REGION)%Past_Hz(i    ,j    ,k))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsRightMUR) then
            REGION = right
            j = MURc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHx)%ZI( REGION), MURc(iHx)%ZE( REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m - 1, k_m    )
                  Hx( i_m, j_m, k_m)=                                            + regLR( REGION)%Past_Hx(i    ,j - 1,k)          &
                  + right_CAB1(medio)*(                   Hx(i_m  ,j_m - 1,k_m) - regLR( REGION)%Past_Hx(i    ,j    ,k))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(iHz)%XI( REGION), MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz( i_m    , j_m - 1, k_m    )
                  Hz( i_m, j_m, k_m) =                                              + regLR( REGION)%Past_Hz(i    ,j - 1,k    )      &
                  + right_CAB1(medio)*(                   Hz(i_m  ,j_m - 1,k_m) - regLR( REGION)%Past_Hz(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsDownMUR) then
            REGION = down
            k = MURc(iHy)%ZI( REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION), MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy( i_m    , j_m    , k_m + 1)
                  Hy( i_m, j_m, k_m) =                                             + regDU( REGION)%Past_Hy(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hy(i_m  ,j_m,k_m + 1) - regDU( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(iHx)%ZI( REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHx)%YI( REGION), MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m    , k_m + 1)
                  Hx( i_m, j_m, k_m) =                                             + regDU( REGION)%Past_Hx(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hx(i_m  ,j_m,k_m + 1) - regDU( REGION)%Past_Hx(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsUpMUR) then
            REGION = up
            k = MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION), MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy( i_m    , j_m    , k_m - 1)
                  Hy( i_m, j_m, k_m) =                                            + regDU( REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hy(i_m  ,j_m    ,k_m - 1) - regDU( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHx)%YI( REGION), MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m    , k_m - 1)
                  Hx( i_m, j_m, k_m) =                                               + regDU( REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                   Hx(i_m  ,j_m  ,k_m - 1)   - regDU( REGION)%Past_Hx(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsBackMUR) then
            REGION =back
            i = MURc(iHz)%XI( REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(iHz)%YI( REGION), MURc(iHz)%YE( REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz( i_m + 1, j_m    , k_m    )
                  Hz( i_m, j_m, k_m) =                                              + regBF( REGION)%Past_Hz(i + 1,j    ,k    )      &
                  + back_CAB1(medio)*(                     Hz(i_m + 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hz(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(iHy)%XI( REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHy)%ZI( REGION), MURc(iHy)%ZE( REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
                  j_m = j - b%Hy%YI
                  !--->orig
                  medio = sggMiHy( i_m + 1, j_m    , k_m    )
                  Hy( i_m, j_m, k_m) =                                              + regBF( REGION)%Past_Hy(i + 1,j    ,k    )      &
                  + back_CAB1(medio)*(                   Hy(i_m + 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsFrontMUR) then
            REGION =front
            i = MURc(iHz)%XE( REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(iHz)%YI( REGION), MURc(iHz)%YE( REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz( i_m - 1, j_m    , k_m    )
                  Hz( i_m, j_m, k_m) =                                               + regBF( REGION)%Past_Hz(i - 1,j    ,k    )      &
                  + front_CAB1(medio)*(                   Hz(i_m - 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hz(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(iHy)%XE( REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHy)%ZI( REGION), MURc(iHy)%ZE( REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
                  j_m = j - b%Hy%YI
                  !--->
                  medio = sggMiHy( i_m - 1, j_m    , k_m    )
                  Hy( i_m, j_m, k_m) =                                               + regBF( REGION)%Past_Hy(i - 1,j    ,k    )      &
                  + front_CAB1(medio)*(                   Hy(i_m - 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !

         !!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!Faces!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (sgg%Border%IsLeftMUR) then
            REGION = left
            j = MURc(iHx)%YI( REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHx)%ZI( REGION) + 1, MURc(iHx)%ZE( REGION) - 1
               k_m = k - b%Hx%ZI
               do i = MURc(iHx)%XI( REGION) + 1, MURc(iHx)%XE( REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m + 1, k_m    )
                  Hx( i_m, j_m, k_m)=                                           - regLR( REGION)%PastPast_Hx(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                     Hx(i_m  ,j_m + 1,k_m) + regLR( REGION)%PastPast_Hx(i    ,j    ,k))         &
                  + left_CAB4(medio)*( regLR( REGION)%Past_Hx(i    ,j    ,k    ) +     regLR( REGION)%Past_Hx(i    ,j + 1,k    ))     &
                  + left_CAB3(medio)*( regLR( REGION)%Past_Hx(i + 1,j    ,k    ) +     regLR( REGION)%Past_Hx(i - 1,j    ,k    )      &
                  +                    regLR( REGION)%Past_Hx(i + 1,j + 1,k    ) +     regLR( REGION)%Past_Hx(i - 1,j + 1,k    )      &
                  +                    regLR( REGION)%Past_Hx(i    ,j    ,k +1 ) +     regLR( REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  +                    regLR( REGION)%Past_Hx(i    ,j + 1,k +1 ) +     regLR( REGION)%Past_Hx(i    ,j + 1,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(iHz)%YI( REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION) + 1, MURc(iHz)%ZE( REGION) - 1
               k_m = k - b%Hz%ZI
               do i = MURc(iHz)%XI( REGION) + 1, MURc(iHz)%XE( REGION) - 1
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz( i_m    , j_m + 1, k_m    )
                  Hz( i_m, j_m, k_m) =                                             - regLR( REGION)%PastPast_Hz(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                     Hz(i_m  ,j_m + 1,k_m) + regLR( REGION)%PastPast_Hz(i    ,j    ,k))         &
                  + left_CAB4(medio)*( regLR( REGION)%Past_Hz(i    ,j    ,k    ) +     regLR( REGION)%Past_Hz(i    ,j + 1,k    ))     &
                  + left_CAB3(medio)*( regLR( REGION)%Past_Hz(i + 1,j    ,k    ) +     regLR( REGION)%Past_Hz(i - 1,j    ,k    )      &
                  +                    regLR( REGION)%Past_Hz(i + 1,j + 1,k    ) +     regLR( REGION)%Past_Hz(i - 1,j + 1,k    )      &
                  +                    regLR( REGION)%Past_Hz(i    ,j    ,k +1 ) +     regLR( REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                    regLR( REGION)%Past_Hz(i    ,j + 1,k +1 ) +     regLR( REGION)%Past_Hz(i    ,j + 1,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsRightMUR) then
            REGION = right
            j = MURc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHx)%ZI( REGION) + 1, MURc(iHx)%ZE( REGION) - 1
               k_m = k - b%Hx%ZI
               do i = MURc(iHx)%XI( REGION) + 1, MURc(iHx)%XE( REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m - 1, k_m    )
                  Hx( i_m, j_m- 1, k_m)=                                            - regLR( REGION)%PastPast_Hx(i    ,j - 1,k)          &
                  + right_CAB1(medio)*(                     Hx(i_m  ,j_m - 1,k_m) + regLR( REGION)%PastPast_Hx(i    ,j    ,k))         &
                  + right_CAB4(medio)*( regLR( REGION)%Past_Hx(i    ,j    ,k    ) +     regLR( REGION)%Past_Hx(i    ,j - 1,k    ))     &
                  + right_CAB3(medio)*( regLR( REGION)%Past_Hx(i + 1,j    ,k    ) +     regLR( REGION)%Past_Hx(i - 1,j    ,k    )      &
                  +                     regLR( REGION)%Past_Hx(i + 1,j - 1,k    ) +     regLR( REGION)%Past_Hx(i - 1,j - 1,k    )      &
                  +                     regLR( REGION)%Past_Hx(i    ,j    ,k +1 ) +     regLR( REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  +                     regLR( REGION)%Past_Hx(i    ,j - 1,k +1 ) +     regLR( REGION)%Past_Hx(i    ,j - 1,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION) + 1, MURc(iHz)%ZE( REGION) - 1
               k_m = k - b%Hz%ZI
               do i = MURc(iHz)%XI( REGION) + 1, MURc(iHz)%XE( REGION) - 1
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz( i_m    , j_m - 1, k_m    )
                  Hz( i_m, j_m, k_m) =                                              - regLR( REGION)%PastPast_Hz(i    ,j - 1,k    )      &
                  + right_CAB1(medio)*(                     Hz(i_m  ,j_m - 1,k_m) + regLR( REGION)%PastPast_Hz(i    ,j    ,k    ))     &
                  + right_CAB4(medio)*( regLR( REGION)%Past_Hz(i    ,j    ,k    ) +     regLR( REGION)%Past_Hz(i    ,j - 1,k    ))     &
                  + right_CAB3(medio)*( regLR( REGION)%Past_Hz(i + 1,j    ,k    ) +     regLR( REGION)%Past_Hz(i - 1,j    ,k    )      &
                  +                     regLR( REGION)%Past_Hz(i + 1,j - 1,k    ) +     regLR( REGION)%Past_Hz(i - 1,j - 1,k    )      &
                  +                     regLR( REGION)%Past_Hz(i    ,j    ,k +1 ) +     regLR( REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                     regLR( REGION)%Past_Hz(i    ,j - 1,k +1 ) +     regLR( REGION)%Past_Hz(i    ,j - 1,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsDownMUR) then
            REGION = down
            k = MURc(iHy)%ZI( REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHy)%YI( REGION) + 1, MURc(iHy)%YE( REGION) - 1
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION) + 1, MURc(iHy)%XE( REGION) - 1
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy( i_m    , j_m    , k_m + 1)
                  Hy( i_m, j_m, k_m) =                                             - regDU( REGION)%PastPast_Hy(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                     Hy(i_m  ,j_m,k_m + 1) + regDU( REGION)%PastPast_Hy(i    ,j    ,k    ))     &
                  + down_CAB4(medio)*( regDU( REGION)%Past_Hy(i    ,j    ,k    ) +     regDU( REGION)%Past_Hy(i    ,j    ,k + 1))     &
                  + down_CAB3(medio)*( regDU( REGION)%Past_Hy(i + 1,j    ,k    ) +     regDU( REGION)%Past_Hy(i - 1,j    ,k    )      &
                  +                    regDU( REGION)%Past_Hy(i + 1,j    ,k + 1) +     regDU( REGION)%Past_Hy(i - 1,j    ,k + 1)      &
                  +                    regDU( REGION)%Past_Hy(i    ,j +1 ,k    ) +     regDU( REGION)%Past_Hy(i    ,j - 1,k    )      &
                  +                    regDU( REGION)%Past_Hy(i    ,j +1 ,k + 1) +     regDU( REGION)%Past_Hy(i    ,j - 1,k + 1))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(iHx)%ZI( REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHx)%YI( REGION) + 1, MURc(iHx)%YE( REGION) - 1
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION) + 1, MURc(iHx)%XE( REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m    , k_m + 1)
                  Hx( i_m, j_m, k_m) =                                             - regDU( REGION)%PastPast_Hx(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                     Hx(i_m  ,j_m,k_m + 1) + regDU( REGION)%PastPast_Hx(i    ,j    ,k    ))     &
                  + down_CAB4(medio)*( regDU( REGION)%Past_Hx(i    ,j    ,k    ) +     regDU( REGION)%Past_Hx(i    ,j    ,k + 1))     &
                  + down_CAB3(medio)*( regDU( REGION)%Past_Hx(i + 1,j    ,k    ) +     regDU( REGION)%Past_Hx(i - 1,j    ,k    )      &
                  +                    regDU( REGION)%Past_Hx(i + 1,j    ,k + 1) +     regDU( REGION)%Past_Hx(i - 1,j    ,k + 1)      &
                  +                    regDU( REGION)%Past_Hx(i    ,j +1 ,k    ) +     regDU( REGION)%Past_Hx(i    ,j - 1,k    )      &
                  +                    regDU( REGION)%Past_Hx(i    ,j +1 ,k + 1) +     regDU( REGION)%Past_Hx(i    ,j - 1,k + 1))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsUpMUR) then
            REGION = up
            k = MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHy)%YI( REGION) + 1, MURc(iHy)%YE( REGION) - 1
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION) + 1, MURc(iHy)%XE( REGION) - 1
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy( i_m    , j_m    , k_m - 1)
                  Hy( i_m, j_m, k_m) =                                            - regDU( REGION)%PastPast_Hy(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hy(i_m  ,j_m    ,k_m - 1) + regDU( REGION)%PastPast_Hy(i    ,j    ,k    ))     &
                  + up_CAB4(medio)*( regDU( REGION)%Past_Hy(i    ,j     ,k    ) +     regDU( REGION)%Past_Hy(i    ,j    ,k - 1))     &
                  + up_CAB3(medio)*( regDU( REGION)%Past_Hy(i + 1,j     ,k    ) +     regDU( REGION)%Past_Hy(i - 1,j    ,k    )      &
                  +                  regDU( REGION)%Past_Hy(i + 1,j     ,k - 1) +     regDU( REGION)%Past_Hy(i - 1,j    ,k - 1)      &
                  +                  regDU( REGION)%Past_Hy(i    ,j +1  ,k    ) +     regDU( REGION)%Past_Hy(i    ,j - 1,k    )      &
                  +                  regDU( REGION)%Past_Hy(i    ,j +1  ,k - 1) +     regDU( REGION)%Past_Hy(i    ,j - 1,k - 1))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHx)%YI( REGION) + 1, MURc(iHx)%YE( REGION) - 1
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION) + 1, MURc(iHx)%XE( REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m    , k_m - 1)
                  Hx( i_m, j_m, k_m) =                                               - regDU( REGION)%PastPast_Hx(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hx(i_m  ,j_m  ,k_m - 1)   + regDU( REGION)%PastPast_Hx(i    ,j    ,k    ))     &
                  + up_CAB4(medio)*( regDU( REGION)%Past_Hx(i    ,j    ,k      )   +     regDU( REGION)%Past_Hx(i    ,j    ,k - 1))     &
                  + up_CAB3(medio)*( regDU( REGION)%Past_Hx(i + 1,j    ,k      )   +     regDU( REGION)%Past_Hx(i - 1,j    ,k    )      &
                  +                  regDU( REGION)%Past_Hx(i + 1,j    ,k - 1  )   +     regDU( REGION)%Past_Hx(i - 1,j    ,k - 1)      &
                  +                  regDU( REGION)%Past_Hx(i    ,j +1 ,k      )   +     regDU( REGION)%Past_Hx(i    ,j - 1,k    )      &
                  +                  regDU( REGION)%Past_Hx(i    ,j +1 ,k - 1  )   +     regDU( REGION)%Past_Hx(i    ,j - 1,k - 1))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsBackMUR) then
            REGION =back
            i = MURc(iHz)%XI( REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION) + 1, MURc(iHz)%ZE( REGION) - 1
               k_m = k - b%Hz%ZI
               do j = MURc(iHz)%YI( REGION) + 1, MURc(iHz)%YE( REGION) - 1
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz( i_m + 1, j_m    , k_m    )
                  Hz( i_m, j_m, k_m) =                                              - regBF( REGION)%PastPast_Hz(i + 1,j    ,k    )      &
                  + back_CAB1(medio)*(                     Hz(i_m + 1,j_m  ,k_m    ) + regBF( REGION)%PastPast_Hz(i    ,j    ,k    ))     &
                  + back_CAB4(medio)*( regBF( REGION)%Past_Hz(i      ,j    ,k   ) +     regBF( REGION)%Past_Hz(i + 1,j    ,k    ))     &
                  + back_CAB3(medio)*( regBF( REGION)%Past_Hz(i      ,j + 1,k   ) +     regBF( REGION)%Past_Hz(i    ,j - 1,k    )      &
                  +                    regBF( REGION)%Past_Hz(i + 1  ,j + 1,k   ) +     regBF( REGION)%Past_Hz(i + 1,j - 1,k    )      &
                  +                    regBF( REGION)%Past_Hz(i      ,j    ,k +1) +     regBF( REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                    regBF( REGION)%Past_Hz(i + 1  ,j    ,k +1) +     regBF( REGION)%Past_Hz(i + 1,j    ,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(iHy)%XI( REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHy)%ZI( REGION) + 1, MURc(iHy)%ZE( REGION) - 1
               k_m = k - b%Hy%ZI
               do j = MURc(iHy)%YI( REGION) + 1, MURc(iHy)%YE( REGION) - 1
                  j_m = j - b%Hy%YI
                  !--->orig
                  medio = sggMiHy( i_m + 1, j_m    , k_m    )
                  Hy( i_m, j_m, k_m) =                                              - regBF( REGION)%PastPast_Hy(i + 1,j    ,k    )      &
                  + back_CAB1(medio)*(                     Hy(i_m + 1,j_m  ,k_m    ) + regBF( REGION)%PastPast_Hy(i    ,j    ,k    ))     &
                  + back_CAB4(medio)*( regBF( REGION)%Past_Hy(i      ,j    ,k   ) +     regBF( REGION)%Past_Hy(i + 1,j    ,k    ))     &
                  + back_CAB3(medio)*( regBF( REGION)%Past_Hy(i      ,j + 1,k   ) +     regBF( REGION)%Past_Hy(i    ,j - 1,k    )      &
                  +                    regBF( REGION)%Past_Hy(i + 1  ,j + 1,k   ) +     regBF( REGION)%Past_Hy(i + 1,j - 1,k    )      &
                  +                    regBF( REGION)%Past_Hy(i      ,j    ,k +1) +     regBF( REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  +                    regBF( REGION)%Past_Hy(i + 1  ,j    ,k +1) +     regBF( REGION)%Past_Hy(i + 1,j    ,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsFrontMUR) then
            REGION =front
            i = MURc(iHz)%XE( REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION) + 1, MURc(iHz)%ZE( REGION) - 1
               k_m = k - b%Hz%ZI
               do j = MURc(iHz)%YI( REGION) + 1, MURc(iHz)%YE( REGION) - 1
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz( i_m - 1, j_m    , k_m    )
                  Hz( i_m, j_m, k_m) =                                               - regBF( REGION)%PastPast_Hz(i - 1,j    ,k    )      &
                  + front_CAB1(medio)*(                     Hz(i_m - 1,j_m  ,k_m    ) + regBF( REGION)%PastPast_Hz(i    ,j    ,k    ))     &
                  + front_CAB4(medio)*( regBF( REGION)%Past_Hz(i      ,j    ,k   ) +     regBF( REGION)%Past_Hz(i - 1,j    ,k    ))     &
                  + front_CAB3(medio)*( regBF( REGION)%Past_Hz(i      ,j + 1,k   ) +     regBF( REGION)%Past_Hz(i    ,j - 1,k    )      &
                  +                     regBF( REGION)%Past_Hz(i - 1  ,j + 1,k   ) +     regBF( REGION)%Past_Hz(i - 1,j - 1,k    )      &
                  +                     regBF( REGION)%Past_Hz(i      ,j    ,k +1) +     regBF( REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                     regBF( REGION)%Past_Hz(i - 1  ,j    ,k +1) +     regBF( REGION)%Past_Hz(i - 1,j    ,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(iHy)%XE( REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHy)%ZI( REGION) + 1, MURc(iHy)%ZE( REGION) - 1
               k_m = k - b%Hy%ZI
               do j = MURc(iHy)%YI( REGION) + 1, MURc(iHy)%YE( REGION) - 1
                  j_m = j - b%Hy%YI
                  !--->
                  medio = sggMiHy( i_m - 1, j_m    , k_m    )
                  Hy( i_m, j_m, k_m) =                                               - regBF( REGION)%PastPast_Hy(i - 1,j    ,k    )      &
                  + front_CAB1(medio)*(                     Hy(i_m - 1,j_m  ,k_m    ) + regBF( REGION)%PastPast_Hy(i    ,j    ,k    ))     &
                  + front_CAB4(medio)*( regBF( REGION)%Past_Hy(i      ,j    ,k   ) +     regBF( REGION)%Past_Hy(i - 1,j    ,k    ))     &
                  + front_CAB3(medio)*( regBF( REGION)%Past_Hy(i      ,j + 1,k   ) +     regBF( REGION)%Past_Hy(i    ,j - 1,k    )      &
                  +                     regBF( REGION)%Past_Hy(i - 1  ,j + 1,k   ) +     regBF( REGION)%Past_Hy(i - 1,j - 1,k    )      &
                  +                     regBF( REGION)%Past_Hy(i      ,j    ,k +1) +     regBF( REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  +                     regBF( REGION)%Past_Hy(i - 1  ,j    ,k +1) +     regBF( REGION)%Past_Hy(i - 1,j    ,k - 1))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!FIRST ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!FIRST ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!FIRST ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!FIRST ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!FIRST ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!FIRST ORDER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      else !first order mur
         IF (sgg%Border%IsLeftMUR) then
            REGION = left
            j = MURc(iHx)%YI( REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHx)%ZI( REGION), MURc(iHx)%ZE( REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m + 1, k_m    )
                  Hx( i_m, j_m, k_m)=                                           + regLR( REGION)%Past_Hx(i    ,j + 1,k)          &
                  +  left_CAB1(medio)*(                    Hx(i_m  ,j_m + 1,k_m) - regLR( REGION)%Past_Hx(i    ,j    ,k))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(iHz)%YI( REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(iHz)%XI( REGION), MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz( i_m    , j_m + 1, k_m    )
                  Hz( i_m, j_m, k_m) =                                             + regLR( REGION)%Past_Hz(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                   Hz(i_m  ,j_m + 1,k_m) - regLR( REGION)%Past_Hz(i    ,j    ,k))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsRightMUR) then
            REGION = right
            j = MURc(iHx)%YE( REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHx)%ZI( REGION), MURc(iHx)%ZE( REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m - 1, k_m    )
                  Hx( i_m, j_m, k_m)=                                            + regLR( REGION)%Past_Hx(i    ,j - 1,k)          &
                  + right_CAB1(medio)*(                   Hx(i_m  ,j_m - 1,k_m) - regLR( REGION)%Past_Hx(i    ,j    ,k))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(iHz)%YE( REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(iHz)%XI( REGION), MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz( i_m    , j_m - 1, k_m    )
                  Hz( i_m, j_m, k_m) =                                              + regLR( REGION)%Past_Hz(i    ,j - 1,k    )      &
                  + right_CAB1(medio)*(                   Hz(i_m  ,j_m - 1,k_m) - regLR( REGION)%Past_Hz(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsDownMUR) then
            REGION = down
            k = MURc(iHy)%ZI( REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION), MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy( i_m    , j_m    , k_m + 1)
                  Hy( i_m, j_m, k_m) =                                             + regDU( REGION)%Past_Hy(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hy(i_m  ,j_m,k_m + 1) - regDU( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(iHx)%ZI( REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHx)%YI( REGION), MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m    , k_m + 1)
                  Hx( i_m, j_m, k_m) =                                             + regDU( REGION)%Past_Hx(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hx(i_m  ,j_m,k_m + 1) - regDU( REGION)%Past_Hx(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsUpMUR) then
            REGION = up
            k = MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION), MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy( i_m    , j_m    , k_m - 1)
                  Hy( i_m, j_m, k_m) =                                            + regDU( REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hy(i_m  ,j_m    ,k_m - 1) - regDU( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(iHx)%YI( REGION), MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx( i_m    , j_m    , k_m - 1)
                  Hx( i_m, j_m, k_m) =                                               + regDU( REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                   Hx(i_m  ,j_m  ,k_m - 1)   - regDU( REGION)%Past_Hx(i    ,j    ,k    ))
               enddo !bucle i
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsBackMUR) then
            REGION =back
            i = MURc(iHz)%XI( REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(iHz)%YI( REGION), MURc(iHz)%YE( REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz( i_m + 1, j_m    , k_m    )
                  Hz( i_m, j_m, k_m) =                                              + regBF( REGION)%Past_Hz(i + 1,j    ,k    )      &
                  + back_CAB1(medio)*(                     Hz(i_m + 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hz(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(iHy)%XI( REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHy)%ZI( REGION), MURc(iHy)%ZE( REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
                  j_m = j - b%Hy%YI
                  !--->orig
                  medio = sggMiHy( i_m + 1, j_m    , k_m    )
                  Hy( i_m, j_m, k_m) =                                              + regBF( REGION)%Past_Hy(i + 1,j    ,k    )      &
                  + back_CAB1(medio)*(                   Hy(i_m + 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         IF (sgg%Border%IsFrontMUR) then
            REGION =front
            i = MURc(iHz)%XE( REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(iHz)%YI( REGION), MURc(iHz)%YE( REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz( i_m - 1, j_m    , k_m    )
                  Hz( i_m, j_m, k_m) =                                               + regBF( REGION)%Past_Hz(i - 1,j    ,k    )      &
                  + front_CAB1(medio)*(                   Hz(i_m - 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hz(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(iHy)%XE( REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(iHy)%ZI( REGION), MURc(iHy)%ZE( REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
                  j_m = j - b%Hy%YI
                  !--->
                  medio = sggMiHy( i_m - 1, j_m    , k_m    )
                  Hy( i_m, j_m, k_m) =                                               + regBF( REGION)%Past_Hy(i - 1,j    ,k    )      &
                  + front_CAB1(medio)*(                   Hy(i_m - 1,j_m  ,k_m    ) - regBF( REGION)%Past_Hy(i    ,j    ,k    ))
               enddo
            enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         endif
         !
      endif !del if mur_second_order


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !guardar los past y pastpast
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!Total!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (sgg%Border%IsLeftMUR) then
         REGION = left
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHx)%ZI( REGION) , MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(iHx)%YI( REGION)  , MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION) , MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regLR( REGION)%PastPast_Hx(i,j,k) = regLR( REGION)%Past_Hx(i  ,j   ,k   )
                  regLR( REGION)%Past_Hx    (i,j,k) =                     Hx(i_m, j_m, k_m)
               enddo
            end do
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHz)%ZI( REGION) , MURc(iHz)%ZE( REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(iHz)%YI( REGION) , MURc(iHz)%YE( REGION)
               j_m = j - b%Hz%YI
               do i = MURc(iHz)%XI( REGION) , MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regLR( REGION)%PastPast_Hz(i,j,k) = regLR( REGION)%Past_Hz(i  ,j   ,k   )
                  regLR( REGION)%Past_Hz    (i,j,k) =                     Hz( i_m, j_m, k_m)
               enddo
            end do
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (sgg%Border%IsRightMUR) then
         REGION = right
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHx)%ZI( REGION), MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(iHx)%YI( REGION) , MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regLR( REGION)%PastPast_Hx(i,j,k) = regLR( REGION)%Past_Hx(i  ,j   ,k   )
                  regLR( REGION)%Past_Hx    (i,j,k) =                     Hx( i_m, j_m, k_m)
               enddo
            end do
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(iHz)%YI( REGION) , MURc(iHz)%YE( REGION)
               j_m = j - b%Hz%YI
               do i = MURc(iHz)%XI( REGION), MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regLR( REGION)%PastPast_Hz(i,j,k) = regLR( REGION)%Past_Hz(i  ,j   ,k   )
                  regLR( REGION)%Past_Hz    (i,j,k) =                     Hz( i_m, j_m, k_m)
               enddo
            end do
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (sgg%Border%IsDownMUR) then
         REGION = down
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHy)%ZI( REGION) , MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(iHy)%YI( REGION) , MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION) , MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  regDU( REGION)%PastPast_Hy(i,j,k) = regDU( REGION)%Past_Hy(i  ,j   ,k   )
                  regDU( REGION)%Past_Hy    (i,j,k) =                     Hy( i_m, j_m, k_m)
               enddo !bucle i
            enddo
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHx)%ZI( REGION) , MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(iHx)%YI( REGION), MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION) , MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regDU( REGION)%PastPast_Hx(i,j,k) = regDU( REGION)%Past_Hx(i  ,j   ,k   )
                  regDU( REGION)%Past_Hx    (i,j,k) =                     Hx( i_m, j_m, k_m)
               enddo !bucle i
            enddo
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (sgg%Border%IsUpMUR) then
         REGION = up
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHy)%ZI( REGION) , MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(iHy)%YI( REGION) , MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION), MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  regDU( REGION)%PastPast_Hy(i,j,k) = regDU( REGION)%Past_Hy(i  ,j   ,k   )
                  regDU( REGION)%Past_Hy    (i,j,k) =                     Hy( i_m, j_m, k_m)
               enddo !bucle i
            enddo
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHx)%ZI( REGION) , MURc(iHx)%ZE( REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(iHx)%YI( REGION), MURc(iHx)%YE( REGION)
               j_m = j - b%Hx%YI
               do i = MURc(iHx)%XI( REGION), MURc(iHx)%XE( REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regDU( REGION)%PastPast_Hx(i,j,k) = regDU( REGION)%Past_Hx(i  ,j   ,k   )
                  regDU( REGION)%Past_Hx    (i,j,k) =                     Hx( i_m, j_m, k_m)
               enddo !bucle i
            enddo
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (sgg%Border%IsBackMUR) then
         REGION =back
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(iHz)%YI( REGION), MURc(iHz)%YE( REGION)
               j_m = j - b%Hz%YI
               do i = MURc(iHz)%XI( REGION) , MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regBF( REGION)%PastPast_Hz(i,j,k) = regBF( REGION)%Past_Hz(i  ,j   ,k   )
                  regBF( REGION)%Past_Hz    (i,j,k) =                     Hz( i_m, j_m, k_m)
               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHy)%ZI( REGION), MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION) , MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->orig
                  regBF( REGION)%PastPast_Hy(i,j,k) = regBF( REGION)%Past_Hy(i  ,j   ,k   )
                  regBF( REGION)%Past_Hy    (i,j,k) =                     Hy( i_m, j_m, k_m)
               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (sgg%Border%IsFrontMUR) then
         REGION =front
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHz)%ZI( REGION), MURc(iHz)%ZE( REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(iHz)%YI( REGION), MURc(iHz)%YE( REGION)
               j_m = j - b%Hz%YI
               do i = MURc(iHz)%XI( REGION) , MURc(iHz)%XE( REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regBF( REGION)%PastPast_Hz(i,j,k) = regBF( REGION)%Past_Hz(i  ,j   ,k   )
                  regBF( REGION)%Past_Hz    (i,j,k) =                     Hz( i_m, j_m, k_m)
               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(iHy)%ZI( REGION) , MURc(iHy)%ZE( REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(iHy)%YI( REGION), MURc(iHy)%YE( REGION)
               j_m = j - b%Hy%YI
               do i = MURc(iHy)%XI( REGION) , MURc(iHy)%XE( REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  regBF( REGION)%PastPast_Hy(i,j,k) = regBF( REGION)%Past_Hy(i  ,j   ,k   )
                  regBF( REGION)%Past_Hy    (i,j,k) =                     Hy( i_m, j_m, k_m)
               end do
            enddo
         enddo
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      endif

      !---------------------------> acaba AdvanceMagneTicMUR <---------------------------------------
      return
   endsubroutine AdvanceMagneTicMUR


end Module Borders_Mur


