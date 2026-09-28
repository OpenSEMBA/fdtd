
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Borders :  MUR  handling
!  Creation date Date :  January, 8, 2013
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BORDERS_MUR_m
   use FDETYPES_m
   use Report_m
   implicit none
   private
   !
   !
   ! Limits of the MUR region
   type xyzlimit_var_t
      integer(kind=4), dimension(1:6) :: XI,XE,YI,YE,ZI,ZE
   end type xyzlimit_var_t
   type(xyzlimit_var_t), dimension(4:6) :: MURc



   type LR_t
      real(kind=RKIND) , pointer, dimension(: , : , :) :: Past_Hx,Past_Hz,PastPast_Hx,PastPast_Hz
   end type
   type DU_t
      real(kind=RKIND) , pointer, dimension(: , : , :) :: Past_Hy,Past_Hx,PastPast_Hy,PastPast_Hx
   end type
   type BF_t
      real(kind=RKIND) , pointer, dimension(: , : , :) :: Past_Hz,Past_Hy,PastPast_Hz,PastPast_Hy
   end type

   !!!LOCAL VARIABLES
   type(LR_t), dimension(LEFT : RIGHT) , save :: regLR
   type(DU_t), dimension(DOWN : UP)    , save :: regDU
   type(BF_t), dimension(BACK : front) , save :: regBF


   real(kind = RKIND), dimension(:), allocatable, save :: back_CAB1, back_CAB3, back_cab4, &
   front_CAB1,front_CAB3,front_cab4, &
   left_CAB1, left_CAB3, left_cab4, &
   right_CAB1,right_CAB3,right_cab4, &
   down_CAB1, down_CAB3, down_cab4, &
   up_CAB1,   up_CAB3,   up_cab4
!!!variables globales del modulo
   real(kind=RKIND), save           :: cluz
   real(kind=RKIND), save           :: eps0,mu0
!!!
   !
   public  :: InitMURBorders, AdvanceMagneticMUR,StoreFieldsMURBorders,DestroyMURBorders,calc_murconstants


contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes MUR data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMURBorders(sgg,ThereAreMURBorders,resume,Idxh,Idyh,Idzh,eps00,mu00)
      real(kind=RKIND) :: eps00,mu00

      type(SGGFDTDINFO_t), intent(in) :: sgg

      real(kind=RKIND) , dimension(:)   , intent(in) :: &
      Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      Idzh(sgg%ALLOC(IEZ)%ZI : sgg%ALLOC(IEZ)%ZE)
      !!!
      !
      logical  :: ThereAreMURBorders,resume
      integer(kind=4) :: i,j,k,region,field,i1

      !character(len=BUFSIZE) :: donde
      !integer(kind=4) :: layoutnumber
!
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      cluz=1.0_RKIND/sqrt(eps0*mu0)
!


      !
      ThereAreMURBorders=.false.
      if (sgg%Border%IsBackMUR.or.sgg%Border%IsFrontMUR.or.sgg%Border%IsLeftMUR.or.sgg%Border%IsRightMUR.or. &
      sgg%Border%IsUpMUR.or.sgg%Border%IsDownMUR) ThereAreMURBorders=.true.
      if (.not.(ThereAreMURBorders)) return

      allocate( back_CAB1(0 :  sgg%NumMedia), back_CAB3(0 :  sgg%NumMedia), back_cab4(0 :  sgg%NumMedia), &
      front_CAB1(0 :  sgg%NumMedia),front_CAB3(0 :  sgg%NumMedia),front_cab4(0 :  sgg%NumMedia), &
      left_CAB1(0 :  sgg%NumMedia), left_CAB3(0 :  sgg%NumMedia), left_cab4(0 :  sgg%NumMedia), &
      right_CAB1(0 :  sgg%NumMedia),right_CAB3(0 :  sgg%NumMedia),right_cab4(0 :  sgg%NumMedia), &
      down_CAB1(0 :  sgg%NumMedia), down_CAB3(0 :  sgg%NumMedia), down_cab4(0 :  sgg%NumMedia), &
      up_CAB1(0 :  sgg%NumMedia),   up_CAB3(0 :  sgg%NumMedia),   up_cab4(0 :  sgg%NumMedia) )
      !Find the limits of each of the 6 padding MUR regions for each field component


      do field=IHX,IHZ
         !
         MURc(field)%XI(DOWN)  =                                sgg%Sweep(field)%XI
         MURc(field)%XE(DOWN)  =                                sgg%Sweep(field)%XE
         MURc(field)%YI(DOWN)  =                                sgg%Sweep(field)%YI
         MURc(field)%YE(DOWN)  =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(DOWN)  =                                sgg%Sweep(field)%ZI-1
         MURc(field)%ZE(DOWN)  = MURc(field)%ZI(DOWN) + 1
         !
         MURc(field)%XI(UP)    =                                sgg%Sweep(field)%XI
         MURc(field)%XE(UP)    =                                sgg%Sweep(field)%XE
         MURc(field)%YI(UP)    =                                sgg%Sweep(field)%YI
         MURc(field)%YE(UP)    =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(UP)    =                                sgg%Sweep(field)%ZE
         MURc(field)%ZE(UP)    = MURc(field)%ZI(UP) + 1
         !
         MURc(field)%XI(LEFT)  =                                sgg%Sweep(field)%XI
         MURc(field)%XE(LEFT)  =                                sgg%Sweep(field)%XE
         MURc(field)%YI(LEFT)  =                                sgg%Sweep(field)%YI-1
         MURc(field)%YE(LEFT)  = MURc(field)%YI(LEFT) + 1
         MURc(field)%ZI(LEFT)  =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(LEFT)  =                                sgg%Sweep(field)%ZE
         !
         MURc(field)%XI(RIGHT) =                                sgg%Sweep(field)%XI
         MURc(field)%XE(RIGHT) =                                sgg%Sweep(field)%XE
         MURc(field)%YI(RIGHT) =                                sgg%Sweep(field)%YE
         MURc(field)%YE(RIGHT) = MURc(field)%YI(RIGHT) + 1
         MURc(field)%ZI(RIGHT) =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(RIGHT) =                                sgg%Sweep(field)%ZE
         !
         MURc(field)%XI(BACK)  =                                sgg%Sweep(field)%XI-1
         MURc(field)%XE(BACK)  = MURc(field)%XI(BACK) + 1
         MURc(field)%YI(BACK)  =                                sgg%Sweep(field)%YI
         MURc(field)%YE(BACK)  =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(BACK)  =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(BACK)  =                                sgg%Sweep(field)%ZE
         !
         MURc(field)%XI(Front) =                                sgg%Sweep(field)%XE
         MURc(field)%XE(Front) = MURc(field)%XI(Front) + 1
         MURc(field)%YI(Front) =                                sgg%Sweep(field)%YI
         MURc(field)%YE(Front) =                                sgg%Sweep(field)%YE
         MURc(field)%ZI(Front) =                                sgg%Sweep(field)%ZI
         MURc(field)%ZE(Front) =                                sgg%Sweep(field)%ZE
         !
      end do

      !Fake coms and ends IN CASE OF NO MUR SO THAT NEVER ENTER THE do FOR THESE CASES
      if (.not.(sgg%Border%IsDownMUR)) MURc(4:6)%ZI(DOWN)=MURc(4:6)%ZE(DOWN)+100
      if (.not.(sgg%Border%IsUpMUR))   MURc(4:6)%ZI(UP)  =MURc(4:6)%ZE(UP)  +100
      !
      if (.not.(sgg%Border%IsLeftMUR))  MURc(4:6)%ZI(LEFT) =MURc(4:6)%ZE(LEFT) +100
      if (.not.(sgg%Border%IsRightMUR)) MURc(4:6)%ZI(RIGHT)=MURc(4:6)%ZE(RIGHT)+100
      !
      if (.not.(sgg%Border%IsFrontMUR)) MURc(4:6)%ZI(front)=MURc(4:6)%ZE(front)+100
      if (.not.(sgg%Border%IsBackMUR))  MURc(4:6)%ZI(BACK) =MURc(4:6)%ZE(BACK) +100

      !MUR Field component matrix allocation
      do REGION =LEFT,RIGHT
         allocate (regLR(region)%Past_Hx(MURc(IHX)%XI(region) : MURc(IHX)%XE(region), &
         MURc(IHX)%YI(region) : MURc(IHX)%YE(region), &
         MURc(IHX)%ZI(region) : MURc(IHX)%ZE(region)),&
         regLR(region)%Past_Hz(MURc(IHZ)%XI(region) : MURc(IHZ)%XE(region), &
         MURc(IHZ)%YI(region) : MURc(IHZ)%YE(region), &
         MURc(IHZ)%ZI(region) : MURc(IHZ)%ZE(region)))
         if (.not.resume) then
            regLR(REGION)%Past_Hx=0.0_RKIND ; regLR(REGION)%Past_Hz=0.0_RKIND ;
         else
            do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
               do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
                  read (14) (regLR(region)%Past_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
               end do
            end do
            do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
               do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
                  read (14) (regLR(region)%Past_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
               end do
            end do
         end if
      end do
      do REGION =DOWN,UP
         allocate (regDU(region)%Past_Hy(MURc(IHY)%XI(region) : MURc(IHY)%XE(region), &
         MURc(IHY)%YI(region) : MURc(IHY)%YE(region), &
         MURc(IHY)%ZI(region) : MURc(IHY)%ZE(region)),&
         regDU(region)%Past_Hx(MURc(IHX)%XI(region) : MURc(IHX)%XE(region), &
         MURc(IHX)%YI(region) : MURc(IHX)%YE(region), &
         MURc(IHX)%ZI(region) : MURc(IHX)%ZE(region)))
         if (.not.resume) then
            regDU(REGION)%Past_Hy=0.0_RKIND ; regDU(REGION)%Past_Hx=0.0_RKIND ;
         else
            do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
               do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
                  read (14) (regDU(region)%Past_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
               end do
            end do
            do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
               do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
                  read (14) (regDU(region)%Past_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
               end do
            end do
         end if
      end do
      do REGION =BACK,front
         allocate (regBF(region)%Past_Hz(MURc(IHZ)%XI(region) : MURc(IHZ)%XE(region), &
         MURc(IHZ)%YI(region) : MURc(IHZ)%YE(region), &
         MURc(IHZ)%ZI(region) : MURc(IHZ)%ZE(region)),&
         regBF(region)%Past_Hy(MURc(IHY)%XI(region) : MURc(IHY)%XE(region), &
         MURc(IHY)%YI(region) : MURc(IHY)%YE(region), &
         MURc(IHY)%ZI(region) : MURc(IHY)%ZE(region)))
         if (.not.resume) then
            regBF(REGION)%Past_Hz=0.0_RKIND ; regBF(REGION)%Past_Hy=0.0_RKIND ;
         else
            do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
               do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
                  read (14) (regBF(region)%Past_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
               end do
            end do
            do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
               do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
                  read (14) (regBF(region)%Past_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
               end do
            end do
         end if
      end do

      !past past

      !MUR Field component matrix allocation
      do REGION =LEFT,RIGHT
         allocate (regLR(region)%PastPast_Hx(MURc(IHX)%XI(region) : MURc(IHX)%XE(region), &
         MURc(IHX)%YI(region) : MURc(IHX)%YE(region), &
         MURc(IHX)%ZI(region) : MURc(IHX)%ZE(region)),&
         regLR(region)%PastPast_Hz(MURc(IHZ)%XI(region) : MURc(IHZ)%XE(region), &
         MURc(IHZ)%YI(region) : MURc(IHZ)%YE(region), &
         MURc(IHZ)%ZI(region) : MURc(IHZ)%ZE(region)))
         if (.not.resume) then
            regLR(REGION)%PastPast_Hx=0.0_RKIND ; regLR(REGION)%PastPast_Hz=0.0_RKIND ;
         else
            do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
               do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
                  read (14) (regLR(region)%PastPast_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
               end do
            end do
            do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
               do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
                  read (14) (regLR(region)%PastPast_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
               end do
            end do
         end if
      end do
      do REGION =DOWN,UP
         allocate (regDU(region)%PastPast_Hy(MURc(IHY)%XI(region) : MURc(IHY)%XE(region), &
         MURc(IHY)%YI(region) : MURc(IHY)%YE(region), &
         MURc(IHY)%ZI(region) : MURc(IHY)%ZE(region)),&
         regDU(region)%PastPast_Hx(MURc(IHX)%XI(region) : MURc(IHX)%XE(region), &
         MURc(IHX)%YI(region) : MURc(IHX)%YE(region), &
         MURc(IHX)%ZI(region) : MURc(IHX)%ZE(region)))
         if (.not.resume) then
            regDU(REGION)%PastPast_Hy=0.0_RKIND ; regDU(REGION)%PastPast_Hx=0.0_RKIND ;
         else
            do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
               do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
                  read (14) (regDU(region)%PastPast_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
               end do
            end do
            do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
               do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
                  read (14) (regDU(region)%PastPast_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
               end do
            end do
         end if
      end do
      do REGION =BACK,front
         allocate (regBF(region)%PastPast_Hz(MURc(IHZ)%XI(region) : MURc(IHZ)%XE(region), &
         MURc(IHZ)%YI(region) : MURc(IHZ)%YE(region), &
         MURc(IHZ)%ZI(region) : MURc(IHZ)%ZE(region)),&
         regBF(region)%PastPast_Hy(MURc(IHY)%XI(region) : MURc(IHY)%XE(region), &
         MURc(IHY)%YI(region) : MURc(IHY)%YE(region), &
         MURc(IHY)%ZI(region) : MURc(IHY)%ZE(region)))
         if (.not.resume) then
            regBF(REGION)%PastPast_Hz=0.0_RKIND ; regBF(REGION)%PastPast_Hy=0.0_RKIND ;
         else
            do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
               do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
                  read (14) (regBF(region)%PastPast_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
               end do
            end do
            do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
               do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
                  read (14) (regBF(region)%PastPast_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
               end do
            end do
         end if
      end do

!!incializa constantes
      call calc_murconstants(sgg,Idxh,Idyh,Idzh,eps0,mu0)
 
      return
   end subroutine InitMURBorders

   subroutine calc_murconstants(sgg,Idxh,Idyh,Idzh,eps00,mu00)
        type(SGGFDTDINFO_t), intent(in) :: sgg
        real(kind=RKIND) :: eps00,mu00
        integer(kind=4) :: i,j,k,region,field,i1
        real(kind=RKIND) :: cnum
        real(kind=RKIND) , dimension(:)   , intent(in) :: &
        Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
        Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
        Idzh(sgg%ALLOC(IEZ)%ZI : sgg%ALLOC(IEZ)%ZE)
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
            cnum=(1.0_RKIND/Idzh(sgg%ALLOC(IEZ)%ZI))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
            down_CAB1(i1) = (1.0_RKIND-CNUM)/(1.0_RKIND+CNUM)
            down_CAB3(i1) = 1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))
            down_cab4(i1) = (2.0_RKIND * CNUM/(1.0_RKIND+CNUM)-4.0_RKIND * (1.0_RKIND / (2.0_RKIND * CNUM*(1.0_RKIND+CNUM))))
            !
            cnum=(1.0_RKIND/Idzh(sgg%ALLOC(IEZ)%ZE))/(sgg%dt * cluz/sqrt(sgg%Med(i1)%Epr * sgg%Med(i1)%Mur))
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

      integer(kind=4) :: region,i,j,k


      do REGION =LEFT,RIGHT
         do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
            do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
               write(14,err=634) (regLR(region)%Past_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
            end do
         end do
         do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
            do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
               write(14,err=634) (regLR(region)%Past_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
            end do
         end do
      end do
      do REGION =DOWN,UP
         do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
            do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
               write(14,err=634) (regDU(region)%Past_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
            end do
         end do
         do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
            do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
               write(14,err=634) (regDU(region)%Past_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
            end do
         end do
      end do
      do REGION =BACK,front
         do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
            do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
               write(14,err=634) (regBF(region)%Past_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
            end do
         end do
         do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
            do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
               write(14,err=634) (regBF(region)%Past_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
            end do
         end do
      end do


      do REGION =LEFT,RIGHT
         do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
            do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
               write(14,err=634) (regLR(region)%PastPast_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
            end do
         end do
         do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
            do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
               write(14,err=634) (regLR(region)%PastPast_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
            end do
         end do
      end do
      do REGION =DOWN,UP
         do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
            do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
               write(14,err=634) (regDU(region)%PastPast_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
            end do
         end do
         do k=MURc(IHX)%ZI(region),MURc(IHX)%ZE(region)
            do j=MURc(IHX)%YI(region),MURc(IHX)%YE(region)
               write(14,err=634) (regDU(region)%PastPast_Hx(i,j,k),i=MURc(IHX)%XI(region),MURc(IHX)%XE(region))
            end do
         end do
      end do
      do REGION =BACK,front
         do k=MURc(IHZ)%ZI(region),MURc(IHZ)%ZE(region)
            do j=MURc(IHZ)%YI(region),MURc(IHZ)%YE(region)
               write(14,err=634) (regBF(region)%PastPast_Hz(i,j,k),i=MURc(IHZ)%XI(region),MURc(IHZ)%XE(region))
            end do
         end do
         do k=MURc(IHY)%ZI(region),MURc(IHY)%ZE(region)
            do j=MURc(IHY)%YI(region),MURc(IHY)%YE(region)
               write(14,err=634) (regBF(region)%PastPast_Hy(i,j,k),i=MURc(IHY)%XI(region),MURc(IHY)%XE(region))
            end do
         end do
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
      integer(kind=4) :: region

      do REGION =LEFT,RIGHT
         if (associated(regLR(region)%Past_Hx)) deallocate(regLR(region)%Past_Hx,regLR(region)%Past_Hz)
      end do
      do REGION =DOWN,UP
         if (associated(regDU(region)%Past_Hy)) deallocate(regDU(region)%Past_Hy,regDU(region)%Past_Hx)
      end do
      do REGION =BACK,front
         if (associated(regBF(region)%Past_Hz)) deallocate(regBF(region)%Past_Hz,regBF(region)%Past_Hy)
      end do


      do REGION =LEFT,RIGHT
         if (associated(regLR(region)%PastPast_Hx)) deallocate(regLR(region)%PastPast_Hx,regLR(region)%PastPast_Hz)
      end do
      do REGION =DOWN,UP
         if (associated(regDU(region)%PastPast_Hy)) deallocate(regDU(region)%PastPast_Hy,regDU(region)%PastPast_Hx)
      end do
      do REGION =BACK,front
         if (associated(regBF(region)%PastPast_Hz)) deallocate(regBF(region)%PastPast_Hz,regBF(region)%PastPast_Hy)
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
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(bounds_t), intent(in) :: b
      logical :: mur_second
      !--->
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent(in) :: sggMiHx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent(in) :: sggMiHy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent(in) :: sggMiHz
      !--->
      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(inout) :: Hx
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(inout) :: Hy
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(inout) :: Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer(kind=4) :: REGION, i, j, k, medio, i_m, j_m, k_m
      !---------------------------> empieza AdvanceMagneTicMUR <-------------------------------------

      !Hetic Fields MUR Zone
      !primero hay que updatear los edges porque las caras los utilizan



      if (mur_second) then


         call stoponerror(0,0,'ERROR: MUR SECOND not correctly implemented')
         !!!!!!?!?!?!?!?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!?!?!?!?!?!!!!!!!!!!!!!!!!!!!Edges!!!!!!!!!!!!!!!!!!!!!Mur primer orden
         !!!!!!?!?!?!?!?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsLeftMUR) then
            REGION = LEFT
            j = MURc(IHX)%YI(REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHX)%ZI(REGION), MURc(IHX)%ZE(REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m + 1, k_m)
                  Hx(i_m, j_m, k_m)=                                           + regLR(REGION)%Past_Hx(i    ,j + 1,k)          &
                  +  left_CAB1(medio)*(                    Hx(i_m  ,j_m + 1,k_m) - regLR(REGION)%Past_Hx(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(IHZ)%YI(REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(IHZ)%XI(REGION), MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz(i_m    , j_m + 1, k_m)
                  Hz(i_m, j_m, k_m) =                                             + regLR(REGION)%Past_Hz(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                   Hz(i_m  ,j_m + 1,k_m) - regLR(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsRightMUR) then
            REGION = RIGHT
            j = MURc(IHX)%YE(REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHX)%ZI(REGION), MURc(IHX)%ZE(REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m - 1, k_m)
                  Hx(i_m, j_m, k_m)=                                            + regLR(REGION)%Past_Hx(i    ,j - 1,k)          &
                  + right_CAB1(medio)*(                   Hx(i_m  ,j_m - 1,k_m) - regLR(REGION)%Past_Hx(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(IHZ)%YE(REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(IHZ)%XI(REGION), MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz(i_m    , j_m - 1, k_m)
                  Hz(i_m, j_m, k_m) =                                              + regLR(REGION)%Past_Hz(i    ,j - 1,k)      &
                  + right_CAB1(medio)*(                   Hz(i_m  ,j_m - 1,k_m) - regLR(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsDownMUR) then
            REGION = DOWN
            k = MURc(IHY)%ZI(REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION), MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy(i_m    , j_m    , k_m + 1)
                  Hy(i_m, j_m, k_m) =                                             + regDU(REGION)%Past_Hy(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hy(i_m  ,j_m,k_m + 1) - regDU(REGION)%Past_Hy(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(IHX)%ZI(REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHX)%YI(REGION), MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m    , k_m + 1)
                  Hx(i_m, j_m, k_m) =                                             + regDU(REGION)%Past_Hx(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hx(i_m  ,j_m,k_m + 1) - regDU(REGION)%Past_Hx(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsUpMUR) then
            REGION = UP
            k = MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION), MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy(i_m    , j_m    , k_m - 1)
                  Hy(i_m, j_m, k_m) =                                            + regDU(REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hy(i_m  ,j_m    ,k_m - 1) - regDU(REGION)%Past_Hy(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHX)%YI(REGION), MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m    , k_m - 1)
                  Hx(i_m, j_m, k_m) =                                               + regDU(REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                   Hx(i_m  ,j_m  ,k_m - 1)   - regDU(REGION)%Past_Hx(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsBackMUR) then
            REGION =BACK
            i = MURc(IHZ)%XI(REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(IHZ)%YI(REGION), MURc(IHZ)%YE(REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz(i_m + 1, j_m    , k_m)
                  Hz(i_m, j_m, k_m) =                                              + regBF(REGION)%Past_Hz(i + 1,j    ,k)      &
                  + back_CAB1(medio)*(                     Hz(i_m + 1,j_m  ,k_m) - regBF(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(IHY)%XI(REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHY)%ZI(REGION), MURc(IHY)%ZE(REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
                  j_m = j - b%Hy%YI
                  !--->orig
                  medio = sggMiHy(i_m + 1, j_m    , k_m)
                  Hy(i_m, j_m, k_m) =                                              + regBF(REGION)%Past_Hy(i + 1,j    ,k)      &
                  + back_CAB1(medio)*(                   Hy(i_m + 1,j_m  ,k_m) - regBF(REGION)%Past_Hy(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsFrontMUR) then
            REGION =front
            i = MURc(IHZ)%XE(REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(IHZ)%YI(REGION), MURc(IHZ)%YE(REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz(i_m - 1, j_m    , k_m)
                  Hz(i_m, j_m, k_m) =                                               + regBF(REGION)%Past_Hz(i - 1,j    ,k)      &
                  + front_CAB1(medio)*(                   Hz(i_m - 1,j_m  ,k_m) - regBF(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(IHY)%XE(REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHY)%ZI(REGION), MURc(IHY)%ZE(REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
                  j_m = j - b%Hy%YI
                  !--->
                  medio = sggMiHy(i_m - 1, j_m    , k_m)
                  Hy(i_m, j_m, k_m) =                                               + regBF(REGION)%Past_Hy(i - 1,j    ,k)      &
                  + front_CAB1(medio)*(                   Hy(i_m - 1,j_m  ,k_m) - regBF(REGION)%Past_Hy(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
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
         if (sgg%Border%IsLeftMUR) then
            REGION = LEFT
            j = MURc(IHX)%YI(REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHX)%ZI(REGION) + 1, MURc(IHX)%ZE(REGION) - 1
               k_m = k - b%Hx%ZI
               do i = MURc(IHX)%XI(REGION) + 1, MURc(IHX)%XE(REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m + 1, k_m)
                  Hx(i_m, j_m, k_m)=                                           - regLR(REGION)%PastPast_Hx(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                     Hx(i_m  ,j_m + 1,k_m) + regLR(REGION)%PastPast_Hx(i    ,j    ,k))         &
                  + left_CAB4(medio)*( regLR(REGION)%Past_Hx(i    ,j    ,k) +     regLR(REGION)%Past_Hx(i    ,j + 1,k))     &
                  + left_CAB3(medio)*( regLR(REGION)%Past_Hx(i + 1,j    ,k) +     regLR(REGION)%Past_Hx(i - 1,j    ,k)      &
                  +                    regLR(REGION)%Past_Hx(i + 1,j + 1,k) +     regLR(REGION)%Past_Hx(i - 1,j + 1,k)      &
                  +                    regLR(REGION)%Past_Hx(i    ,j    ,k +1) +     regLR(REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  +                    regLR(REGION)%Past_Hx(i    ,j + 1,k +1) +     regLR(REGION)%Past_Hx(i    ,j + 1,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(IHZ)%YI(REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION) + 1, MURc(IHZ)%ZE(REGION) - 1
               k_m = k - b%Hz%ZI
               do i = MURc(IHZ)%XI(REGION) + 1, MURc(IHZ)%XE(REGION) - 1
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz(i_m    , j_m + 1, k_m)
                  Hz(i_m, j_m, k_m) =                                             - regLR(REGION)%PastPast_Hz(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                     Hz(i_m  ,j_m + 1,k_m) + regLR(REGION)%PastPast_Hz(i    ,j    ,k))         &
                  + left_CAB4(medio)*( regLR(REGION)%Past_Hz(i    ,j    ,k) +     regLR(REGION)%Past_Hz(i    ,j + 1,k))     &
                  + left_CAB3(medio)*( regLR(REGION)%Past_Hz(i + 1,j    ,k) +     regLR(REGION)%Past_Hz(i - 1,j    ,k)      &
                  +                    regLR(REGION)%Past_Hz(i + 1,j + 1,k) +     regLR(REGION)%Past_Hz(i - 1,j + 1,k)      &
                  +                    regLR(REGION)%Past_Hz(i    ,j    ,k +1) +     regLR(REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                    regLR(REGION)%Past_Hz(i    ,j + 1,k +1) +     regLR(REGION)%Past_Hz(i    ,j + 1,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsRightMUR) then
            REGION = RIGHT
            j = MURc(IHX)%YE(REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHX)%ZI(REGION) + 1, MURc(IHX)%ZE(REGION) - 1
               k_m = k - b%Hx%ZI
               do i = MURc(IHX)%XI(REGION) + 1, MURc(IHX)%XE(REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m - 1, k_m)
                  Hx(i_m, j_m- 1, k_m)=                                            - regLR(REGION)%PastPast_Hx(i    ,j - 1,k)          &
                  + right_CAB1(medio)*(                     Hx(i_m  ,j_m - 1,k_m) + regLR(REGION)%PastPast_Hx(i    ,j    ,k))         &
                  + right_CAB4(medio)*( regLR(REGION)%Past_Hx(i    ,j    ,k) +     regLR(REGION)%Past_Hx(i    ,j - 1,k))     &
                  + right_CAB3(medio)*( regLR(REGION)%Past_Hx(i + 1,j    ,k) +     regLR(REGION)%Past_Hx(i - 1,j    ,k)      &
                  +                     regLR(REGION)%Past_Hx(i + 1,j - 1,k) +     regLR(REGION)%Past_Hx(i - 1,j - 1,k)      &
                  +                     regLR(REGION)%Past_Hx(i    ,j    ,k +1) +     regLR(REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  +                     regLR(REGION)%Past_Hx(i    ,j - 1,k +1) +     regLR(REGION)%Past_Hx(i    ,j - 1,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(IHZ)%YE(REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION) + 1, MURc(IHZ)%ZE(REGION) - 1
               k_m = k - b%Hz%ZI
               do i = MURc(IHZ)%XI(REGION) + 1, MURc(IHZ)%XE(REGION) - 1
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz(i_m    , j_m - 1, k_m)
                  Hz(i_m, j_m, k_m) =                                              - regLR(REGION)%PastPast_Hz(i    ,j - 1,k)      &
                  + right_CAB1(medio)*(                     Hz(i_m  ,j_m - 1,k_m) + regLR(REGION)%PastPast_Hz(i    ,j    ,k))     &
                  + right_CAB4(medio)*( regLR(REGION)%Past_Hz(i    ,j    ,k) +     regLR(REGION)%Past_Hz(i    ,j - 1,k))     &
                  + right_CAB3(medio)*( regLR(REGION)%Past_Hz(i + 1,j    ,k) +     regLR(REGION)%Past_Hz(i - 1,j    ,k)      &
                  +                     regLR(REGION)%Past_Hz(i + 1,j - 1,k) +     regLR(REGION)%Past_Hz(i - 1,j - 1,k)      &
                  +                     regLR(REGION)%Past_Hz(i    ,j    ,k +1) +     regLR(REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                     regLR(REGION)%Past_Hz(i    ,j - 1,k +1) +     regLR(REGION)%Past_Hz(i    ,j - 1,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsDownMUR) then
            REGION = DOWN
            k = MURc(IHY)%ZI(REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHY)%YI(REGION) + 1, MURc(IHY)%YE(REGION) - 1
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION) + 1, MURc(IHY)%XE(REGION) - 1
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy(i_m    , j_m    , k_m + 1)
                  Hy(i_m, j_m, k_m) =                                             - regDU(REGION)%PastPast_Hy(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                     Hy(i_m  ,j_m,k_m + 1) + regDU(REGION)%PastPast_Hy(i    ,j    ,k))     &
                  + down_CAB4(medio)*( regDU(REGION)%Past_Hy(i    ,j    ,k) +     regDU(REGION)%Past_Hy(i    ,j    ,k + 1))     &
                  + down_CAB3(medio)*( regDU(REGION)%Past_Hy(i + 1,j    ,k) +     regDU(REGION)%Past_Hy(i - 1,j    ,k)      &
                  +                    regDU(REGION)%Past_Hy(i + 1,j    ,k + 1) +     regDU(REGION)%Past_Hy(i - 1,j    ,k + 1)      &
                  +                    regDU(REGION)%Past_Hy(i    ,j +1 ,k) +     regDU(REGION)%Past_Hy(i    ,j - 1,k)      &
                  +                    regDU(REGION)%Past_Hy(i    ,j +1 ,k + 1) +     regDU(REGION)%Past_Hy(i    ,j - 1,k + 1))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(IHX)%ZI(REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHX)%YI(REGION) + 1, MURc(IHX)%YE(REGION) - 1
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION) + 1, MURc(IHX)%XE(REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m    , k_m + 1)
                  Hx(i_m, j_m, k_m) =                                             - regDU(REGION)%PastPast_Hx(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                     Hx(i_m  ,j_m,k_m + 1) + regDU(REGION)%PastPast_Hx(i    ,j    ,k))     &
                  + down_CAB4(medio)*( regDU(REGION)%Past_Hx(i    ,j    ,k) +     regDU(REGION)%Past_Hx(i    ,j    ,k + 1))     &
                  + down_CAB3(medio)*( regDU(REGION)%Past_Hx(i + 1,j    ,k) +     regDU(REGION)%Past_Hx(i - 1,j    ,k)      &
                  +                    regDU(REGION)%Past_Hx(i + 1,j    ,k + 1) +     regDU(REGION)%Past_Hx(i - 1,j    ,k + 1)      &
                  +                    regDU(REGION)%Past_Hx(i    ,j +1 ,k) +     regDU(REGION)%Past_Hx(i    ,j - 1,k)      &
                  +                    regDU(REGION)%Past_Hx(i    ,j +1 ,k + 1) +     regDU(REGION)%Past_Hx(i    ,j - 1,k + 1))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsUpMUR) then
            REGION = UP
            k = MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHY)%YI(REGION) + 1, MURc(IHY)%YE(REGION) - 1
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION) + 1, MURc(IHY)%XE(REGION) - 1
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy(i_m    , j_m    , k_m - 1)
                  Hy(i_m, j_m, k_m) =                                            - regDU(REGION)%PastPast_Hy(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hy(i_m  ,j_m    ,k_m - 1) + regDU(REGION)%PastPast_Hy(i    ,j    ,k))     &
                  + up_CAB4(medio)*( regDU(REGION)%Past_Hy(i    ,j     ,k) +     regDU(REGION)%Past_Hy(i    ,j    ,k - 1))     &
                  + up_CAB3(medio)*( regDU(REGION)%Past_Hy(i + 1,j     ,k) +     regDU(REGION)%Past_Hy(i - 1,j    ,k)      &
                  +                  regDU(REGION)%Past_Hy(i + 1,j     ,k - 1) +     regDU(REGION)%Past_Hy(i - 1,j    ,k - 1)      &
                  +                  regDU(REGION)%Past_Hy(i    ,j +1  ,k) +     regDU(REGION)%Past_Hy(i    ,j - 1,k)      &
                  +                  regDU(REGION)%Past_Hy(i    ,j +1  ,k - 1) +     regDU(REGION)%Past_Hy(i    ,j - 1,k - 1))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHX)%YI(REGION) + 1, MURc(IHX)%YE(REGION) - 1
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION) + 1, MURc(IHX)%XE(REGION) - 1
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m    , k_m - 1)
                  Hx(i_m, j_m, k_m) =                                               - regDU(REGION)%PastPast_Hx(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hx(i_m  ,j_m  ,k_m - 1)   + regDU(REGION)%PastPast_Hx(i    ,j    ,k))     &
                  + up_CAB4(medio)*( regDU(REGION)%Past_Hx(i    ,j    ,k)   +     regDU(REGION)%Past_Hx(i    ,j    ,k - 1))     &
                  + up_CAB3(medio)*( regDU(REGION)%Past_Hx(i + 1,j    ,k)   +     regDU(REGION)%Past_Hx(i - 1,j    ,k)      &
                  +                  regDU(REGION)%Past_Hx(i + 1,j    ,k - 1)   +     regDU(REGION)%Past_Hx(i - 1,j    ,k - 1)      &
                  +                  regDU(REGION)%Past_Hx(i    ,j +1 ,k)   +     regDU(REGION)%Past_Hx(i    ,j - 1,k)      &
                  +                  regDU(REGION)%Past_Hx(i    ,j +1 ,k - 1)   +     regDU(REGION)%Past_Hx(i    ,j - 1,k - 1))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsBackMUR) then
            REGION =BACK
            i = MURc(IHZ)%XI(REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION) + 1, MURc(IHZ)%ZE(REGION) - 1
               k_m = k - b%Hz%ZI
               do j = MURc(IHZ)%YI(REGION) + 1, MURc(IHZ)%YE(REGION) - 1
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz(i_m + 1, j_m    , k_m)
                  Hz(i_m, j_m, k_m) =                                              - regBF(REGION)%PastPast_Hz(i + 1,j    ,k)      &
                  + back_CAB1(medio)*(                     Hz(i_m + 1,j_m  ,k_m) + regBF(REGION)%PastPast_Hz(i    ,j    ,k))     &
                  + back_CAB4(medio)*( regBF(REGION)%Past_Hz(i      ,j    ,k) +     regBF(REGION)%Past_Hz(i + 1,j    ,k))     &
                  + back_CAB3(medio)*( regBF(REGION)%Past_Hz(i      ,j + 1,k) +     regBF(REGION)%Past_Hz(i    ,j - 1,k)      &
                  +                    regBF(REGION)%Past_Hz(i + 1  ,j + 1,k) +     regBF(REGION)%Past_Hz(i + 1,j - 1,k)      &
                  +                    regBF(REGION)%Past_Hz(i      ,j    ,k +1) +     regBF(REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                    regBF(REGION)%Past_Hz(i + 1  ,j    ,k +1) +     regBF(REGION)%Past_Hz(i + 1,j    ,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(IHY)%XI(REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHY)%ZI(REGION) + 1, MURc(IHY)%ZE(REGION) - 1
               k_m = k - b%Hy%ZI
               do j = MURc(IHY)%YI(REGION) + 1, MURc(IHY)%YE(REGION) - 1
                  j_m = j - b%Hy%YI
                  !--->orig
                  medio = sggMiHy(i_m + 1, j_m    , k_m)
                  Hy(i_m, j_m, k_m) =                                              - regBF(REGION)%PastPast_Hy(i + 1,j    ,k)      &
                  + back_CAB1(medio)*(                     Hy(i_m + 1,j_m  ,k_m) + regBF(REGION)%PastPast_Hy(i    ,j    ,k))     &
                  + back_CAB4(medio)*( regBF(REGION)%Past_Hy(i      ,j    ,k) +     regBF(REGION)%Past_Hy(i + 1,j    ,k))     &
                  + back_CAB3(medio)*( regBF(REGION)%Past_Hy(i      ,j + 1,k) +     regBF(REGION)%Past_Hy(i    ,j - 1,k)      &
                  +                    regBF(REGION)%Past_Hy(i + 1  ,j + 1,k) +     regBF(REGION)%Past_Hy(i + 1,j - 1,k)      &
                  +                    regBF(REGION)%Past_Hy(i      ,j    ,k +1) +     regBF(REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  +                    regBF(REGION)%Past_Hy(i + 1  ,j    ,k +1) +     regBF(REGION)%Past_Hy(i + 1,j    ,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsFrontMUR) then
            REGION =front
            i = MURc(IHZ)%XE(REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION) + 1, MURc(IHZ)%ZE(REGION) - 1
               k_m = k - b%Hz%ZI
               do j = MURc(IHZ)%YI(REGION) + 1, MURc(IHZ)%YE(REGION) - 1
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz(i_m - 1, j_m    , k_m)
                  Hz(i_m, j_m, k_m) =                                               - regBF(REGION)%PastPast_Hz(i - 1,j    ,k)      &
                  + front_CAB1(medio)*(                     Hz(i_m - 1,j_m  ,k_m) + regBF(REGION)%PastPast_Hz(i    ,j    ,k))     &
                  + front_CAB4(medio)*( regBF(REGION)%Past_Hz(i      ,j    ,k) +     regBF(REGION)%Past_Hz(i - 1,j    ,k))     &
                  + front_CAB3(medio)*( regBF(REGION)%Past_Hz(i      ,j + 1,k) +     regBF(REGION)%Past_Hz(i    ,j - 1,k)      &
                  +                     regBF(REGION)%Past_Hz(i - 1  ,j + 1,k) +     regBF(REGION)%Past_Hz(i - 1,j - 1,k)      &
                  +                     regBF(REGION)%Past_Hz(i      ,j    ,k +1) +     regBF(REGION)%Past_Hz(i    ,j    ,k - 1)      &
                  +                     regBF(REGION)%Past_Hz(i - 1  ,j    ,k +1) +     regBF(REGION)%Past_Hz(i - 1,j    ,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(IHY)%XE(REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHY)%ZI(REGION) + 1, MURc(IHY)%ZE(REGION) - 1
               k_m = k - b%Hy%ZI
               do j = MURc(IHY)%YI(REGION) + 1, MURc(IHY)%YE(REGION) - 1
                  j_m = j - b%Hy%YI
                  !--->
                  medio = sggMiHy(i_m - 1, j_m    , k_m)
                  Hy(i_m, j_m, k_m) =                                               - regBF(REGION)%PastPast_Hy(i - 1,j    ,k)      &
                  + front_CAB1(medio)*(                     Hy(i_m - 1,j_m  ,k_m) + regBF(REGION)%PastPast_Hy(i    ,j    ,k))     &
                  + front_CAB4(medio)*( regBF(REGION)%Past_Hy(i      ,j    ,k) +     regBF(REGION)%Past_Hy(i - 1,j    ,k))     &
                  + front_CAB3(medio)*( regBF(REGION)%Past_Hy(i      ,j + 1,k) +     regBF(REGION)%Past_Hy(i    ,j - 1,k)      &
                  +                     regBF(REGION)%Past_Hy(i - 1  ,j + 1,k) +     regBF(REGION)%Past_Hy(i - 1,j - 1,k)      &
                  +                     regBF(REGION)%Past_Hy(i      ,j    ,k +1) +     regBF(REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  +                     regBF(REGION)%Past_Hy(i - 1  ,j    ,k +1) +     regBF(REGION)%Past_Hy(i - 1,j    ,k - 1))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
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
         if (sgg%Border%IsLeftMUR) then
            REGION = LEFT
            j = MURc(IHX)%YI(REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHX)%ZI(REGION), MURc(IHX)%ZE(REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m + 1, k_m)
                  Hx(i_m, j_m, k_m)=                                           + regLR(REGION)%Past_Hx(i    ,j + 1,k)          &
                  +  left_CAB1(medio)*(                    Hx(i_m  ,j_m + 1,k_m) - regLR(REGION)%Past_Hx(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(IHZ)%YI(REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(IHZ)%XI(REGION), MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz(i_m    , j_m + 1, k_m)
                  Hz(i_m, j_m, k_m) =                                             + regLR(REGION)%Past_Hz(i    ,j + 1,k)          &
                  + left_CAB1(medio)*(                   Hz(i_m  ,j_m + 1,k_m) - regLR(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsRightMUR) then
            REGION = RIGHT
            j = MURc(IHX)%YE(REGION)
            j_m = j - b%Hx%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHX)%ZI(REGION), MURc(IHX)%ZE(REGION)
               k_m = k - b%Hx%ZI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m - 1, k_m)
                  Hx(i_m, j_m, k_m)=                                            + regLR(REGION)%Past_Hx(i    ,j - 1,k)          &
                  + right_CAB1(medio)*(                   Hx(i_m  ,j_m - 1,k_m) - regLR(REGION)%Past_Hx(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            j = MURc(IHZ)%YE(REGION)
            j_m = j - b%Hz%YI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,k,i_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do i = MURc(IHZ)%XI(REGION), MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  medio = sggMiHz(i_m    , j_m - 1, k_m)
                  Hz(i_m, j_m, k_m) =                                              + regLR(REGION)%Past_Hz(i    ,j - 1,k)      &
                  + right_CAB1(medio)*(                   Hz(i_m  ,j_m - 1,k_m) - regLR(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsDownMUR) then
            REGION = DOWN
            k = MURc(IHY)%ZI(REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION), MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy(i_m    , j_m    , k_m + 1)
                  Hy(i_m, j_m, k_m) =                                             + regDU(REGION)%Past_Hy(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hy(i_m  ,j_m,k_m + 1) - regDU(REGION)%Past_Hy(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(IHX)%ZI(REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHX)%YI(REGION), MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m    , k_m + 1)
                  Hx(i_m, j_m, k_m) =                                             + regDU(REGION)%Past_Hx(i    ,j    ,k + 1)      &
                  + down_CAB1(medio)*(                   Hx(i_m  ,j_m,k_m + 1) - regDU(REGION)%Past_Hx(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsUpMUR) then
            REGION = UP
            k = MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION), MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  medio = sggMiHy(i_m    , j_m    , k_m - 1)
                  Hy(i_m, j_m, k_m) =                                            + regDU(REGION)%Past_Hy(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                     Hy(i_m  ,j_m    ,k_m - 1) - regDU(REGION)%Past_Hy(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            k = MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,i_m,j_m,medio)
#endif
            do j = MURc(IHX)%YI(REGION), MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  medio = sggMiHx(i_m    , j_m    , k_m - 1)
                  Hx(i_m, j_m, k_m) =                                               + regDU(REGION)%Past_Hx(i    ,j    ,k - 1)      &
                  + up_CAB1(medio)*(                   Hx(i_m  ,j_m  ,k_m - 1)   - regDU(REGION)%Past_Hx(i    ,j    ,k))
               end do !bucle i
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsBackMUR) then
            REGION =BACK
            i = MURc(IHZ)%XI(REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(IHZ)%YI(REGION), MURc(IHZ)%YE(REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz(i_m + 1, j_m    , k_m)
                  Hz(i_m, j_m, k_m) =                                              + regBF(REGION)%Past_Hz(i + 1,j    ,k)      &
                  + back_CAB1(medio)*(                     Hz(i_m + 1,j_m  ,k_m) - regBF(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(IHY)%XI(REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHY)%ZI(REGION), MURc(IHY)%ZE(REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
                  j_m = j - b%Hy%YI
                  !--->orig
                  medio = sggMiHy(i_m + 1, j_m    , k_m)
                  Hy(i_m, j_m, k_m) =                                              + regBF(REGION)%Past_Hy(i + 1,j    ,k)      &
                  + back_CAB1(medio)*(                   Hy(i_m + 1,j_m  ,k_m) - regBF(REGION)%Past_Hy(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         if (sgg%Border%IsFrontMUR) then
            REGION =front
            i = MURc(IHZ)%XE(REGION)
            i_m = i - b%Hz%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
               k_m = k - b%Hz%ZI
               do j = MURc(IHZ)%YI(REGION), MURc(IHZ)%YE(REGION)
                  j_m = j - b%Hz%YI
                  !--->
                  medio = sggMiHz(i_m - 1, j_m    , k_m)
                  Hz(i_m, j_m, k_m) =                                               + regBF(REGION)%Past_Hz(i - 1,j    ,k)      &
                  + front_CAB1(medio)*(                   Hz(i_m - 1,j_m  ,k_m) - regBF(REGION)%Past_Hz(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
            i = MURc(IHY)%XE(REGION)
            i_m = i - b%Hy%XI
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (j,k,j_m,k_m,medio)
#endif
            do k = MURc(IHY)%ZI(REGION), MURc(IHY)%ZE(REGION)
               k_m = k - b%Hy%ZI
               do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
                  j_m = j - b%Hy%YI
                  !--->
                  medio = sggMiHy(i_m - 1, j_m    , k_m)
                  Hy(i_m, j_m, k_m) =                                               + regBF(REGION)%Past_Hy(i - 1,j    ,k)      &
                  + front_CAB1(medio)*(                   Hy(i_m - 1,j_m  ,k_m) - regBF(REGION)%Past_Hy(i    ,j    ,k))
               end do
            end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
         end if
         !
      end if !del if mur_second_order


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !guardar los past y pastpast
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!Total!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (sgg%Border%IsLeftMUR) then
         REGION = LEFT
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHX)%ZI(REGION) , MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(IHX)%YI(REGION)  , MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION) , MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regLR(REGION)%PastPast_Hx(i,j,k) = regLR(REGION)%Past_Hx(i  ,j   ,k)
                  regLR(REGION)%Past_Hx    (i,j,k) =                     Hx(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHZ)%ZI(REGION) , MURc(IHZ)%ZE(REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(IHZ)%YI(REGION) , MURc(IHZ)%YE(REGION)
               j_m = j - b%Hz%YI
               do i = MURc(IHZ)%XI(REGION) , MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regLR(REGION)%PastPast_Hz(i,j,k) = regLR(REGION)%Past_Hz(i  ,j   ,k)
                  regLR(REGION)%Past_Hz    (i,j,k) =                     Hz(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (sgg%Border%IsRightMUR) then
         REGION = RIGHT
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHX)%ZI(REGION), MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(IHX)%YI(REGION) , MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regLR(REGION)%PastPast_Hx(i,j,k) = regLR(REGION)%Past_Hx(i  ,j   ,k)
                  regLR(REGION)%Past_Hx    (i,j,k) =                     Hx(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(IHZ)%YI(REGION) , MURc(IHZ)%YE(REGION)
               j_m = j - b%Hz%YI
               do i = MURc(IHZ)%XI(REGION), MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regLR(REGION)%PastPast_Hz(i,j,k) = regLR(REGION)%Past_Hz(i  ,j   ,k)
                  regLR(REGION)%Past_Hz    (i,j,k) =                     Hz(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (sgg%Border%IsDownMUR) then
         REGION = DOWN
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHY)%ZI(REGION) , MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(IHY)%YI(REGION) , MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION) , MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  regDU(REGION)%PastPast_Hy(i,j,k) = regDU(REGION)%Past_Hy(i  ,j   ,k)
                  regDU(REGION)%Past_Hy    (i,j,k) =                     Hy(i_m, j_m, k_m)
               end do !bucle i
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHX)%ZI(REGION) , MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(IHX)%YI(REGION), MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION) , MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regDU(REGION)%PastPast_Hx(i,j,k) = regDU(REGION)%Past_Hx(i  ,j   ,k)
                  regDU(REGION)%Past_Hx    (i,j,k) =                     Hx(i_m, j_m, k_m)
               end do !bucle i
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (sgg%Border%IsUpMUR) then
         REGION = UP
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHY)%ZI(REGION) , MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(IHY)%YI(REGION) , MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION), MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  regDU(REGION)%PastPast_Hy(i,j,k) = regDU(REGION)%Past_Hy(i  ,j   ,k)
                  regDU(REGION)%Past_Hy    (i,j,k) =                     Hy(i_m, j_m, k_m)
               end do !bucle i
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHX)%ZI(REGION) , MURc(IHX)%ZE(REGION)
            k_m = k - b%Hx%ZI
            do j = MURc(IHX)%YI(REGION), MURc(IHX)%YE(REGION)
               j_m = j - b%Hx%YI
               do i = MURc(IHX)%XI(REGION), MURc(IHX)%XE(REGION)
                  i_m = i - b%Hx%XI
                  !--->
                  regDU(REGION)%PastPast_Hx(i,j,k) = regDU(REGION)%Past_Hx(i  ,j   ,k)
                  regDU(REGION)%Past_Hx    (i,j,k) =                     Hx(i_m, j_m, k_m)
               end do !bucle i
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (sgg%Border%IsBackMUR) then
         REGION =BACK
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(IHZ)%YI(REGION), MURc(IHZ)%YE(REGION)
               j_m = j - b%Hz%YI
               do i = MURc(IHZ)%XI(REGION) , MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regBF(REGION)%PastPast_Hz(i,j,k) = regBF(REGION)%Past_Hz(i  ,j   ,k)
                  regBF(REGION)%Past_Hz    (i,j,k) =                     Hz(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHY)%ZI(REGION), MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION) , MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->orig
                  regBF(REGION)%PastPast_Hy(i,j,k) = regBF(REGION)%Past_Hy(i  ,j   ,k)
                  regBF(REGION)%Past_Hy    (i,j,k) =                     Hy(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (sgg%Border%IsFrontMUR) then
         REGION =front
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHZ)%ZI(REGION), MURc(IHZ)%ZE(REGION)
            k_m = k - b%Hz%ZI
            do j = MURc(IHZ)%YI(REGION), MURc(IHZ)%YE(REGION)
               j_m = j - b%Hz%YI
               do i = MURc(IHZ)%XI(REGION) , MURc(IHZ)%XE(REGION)
                  i_m = i - b%Hz%XI
                  !--->
                  regBF(REGION)%PastPast_Hz(i,j,k) = regBF(REGION)%Past_Hz(i  ,j   ,k)
                  regBF(REGION)%Past_Hz    (i,j,k) =                     Hz(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
#ifdef CompileWithOpenMP
!$OMP PARALLEL do DEFAULT(SHARED) private (i,j,k,i_m,j_m,k_m)
#endif
         do k = MURc(IHY)%ZI(REGION) , MURc(IHY)%ZE(REGION)
            k_m = k - b%Hy%ZI
            do j = MURc(IHY)%YI(REGION), MURc(IHY)%YE(REGION)
               j_m = j - b%Hy%YI
               do i = MURc(IHY)%XI(REGION) , MURc(IHY)%XE(REGION)
                  i_m = i - b%Hy%XI
                  !--->
                  regBF(REGION)%PastPast_Hy(i,j,k) = regBF(REGION)%Past_Hy(i  ,j   ,k)
                  regBF(REGION)%Past_Hy    (i,j,k) =                     Hy(i_m, j_m, k_m)
               end do
            end do
         end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif
      end if

      !---------------------------> acaba AdvanceMagneTicMUR <---------------------------------------
      return
   end subroutine AdvanceMagneTicMUR


end module BORDERS_MUR_m


