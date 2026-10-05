
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Mdispersives !note the conjugate poles MUST APPEAR explicitly 20JUNE'12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Beware in MDutton's model BOTH pair OF Complex conjugate poles/residues
! in input from .nfde MUST APPEAR (this is why the factor /2 in the algorithm part
! for instance a 1 real-pole and 2 cComplex-pole material would require in nfde
! 5 poles (not 3) !UNTESTED SGG JUN'12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Mdispersives_m

   !same electric and magnetic switch

   use FDETYPES_m
   use Report_m
   implicit none
   private

   !structures needed by the Mdispersive



   type field_t
      integer(kind=4) :: i,j,k
      integer(kind=4) :: WhatField
      real(kind=RKIND), pointer                 :: FieldPresent !points to the background field
      real(kind=RKIND)                          :: FieldPrevious
      complex(kind=CKIND), pointer, dimension(:) :: Current
   end type

   type Mdispersive_t
      integer(kind=4) :: indexmed,numnodesHx,numnodesHy,numnodesHz,numpolres11
      complex(kind=CKIND), pointer, dimension(:) :: Beta,Kappa,GM3
      type(field_t), pointer, dimension(:) :: NodesHx,NodesHy,NodesHz
   end type Mdispersive_t


   type  Mdispersive2_t
      integer(kind=4) :: NumMdispersives
      type(Mdispersive_t), pointer, dimension(:) :: Medium
   end type

   !!!LOCAL VARIABLES
   type(Mdispersive2_t) , save :: MDutton


   public AdvanceMdispersiveH,InitMdispersives,StoreFieldsMdispersives,DestroyMdispersives

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMdispersives(sgg,media,ThereAreMdispersives,resume,GM1,GM2,Hx,Hy,Hz)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(media_matrices_t), intent(in) :: media
      real(kind=RKIND)     , intent(inout) :: &
      GM1(0 : sgg%NumMedia),GM2(0 : sgg%NumMedia)
      real(kind=RKIND)   , intent(inout), target      :: &
      Hx(sgg%Alloc(IHX)%XI : sgg%Alloc(IHX)%XE,sgg%Alloc(IHX)%YI : sgg%Alloc(IHX)%YE,sgg%Alloc(IHX)%ZI : sgg%Alloc(IHX)%ZE),&
      Hy(sgg%Alloc(IHY)%XI : sgg%Alloc(IHY)%XE,sgg%Alloc(IHY)%YI : sgg%Alloc(IHY)%YE,sgg%Alloc(IHY)%ZI : sgg%Alloc(IHY)%ZE),&
      Hz(sgg%Alloc(IHZ)%XI : sgg%Alloc(IHZ)%XE,sgg%Alloc(IHZ)%YI : sgg%Alloc(IHZ)%YE,sgg%Alloc(IHZ)%ZI : sgg%Alloc(IHZ)%ZE)

      !passing Mdispersive, etc. should be deprecated because there is direct access to sgg%Med%Dispersiv
      logical, intent(out) :: ThereAreMdispersives
      logical, intent(in) :: resume
      integer(kind=4) :: jmed,j1,conta,k1,i1,tempindex
      real(kind=RKIND) :: tempo
      integer(kind=4) :: numpolres
      MDutton%Medium => null()

      ThereAreMdispersives=.FALSE.
      conta=0
      do jmed=1,sgg%NumMedia
         if ((sgg%Med(jmed)%Is%Mdispersive).and.(.not.sgg%Med(jmed)%Is%MdispersiveANIS)) then
            conta=conta+1
         end if
      end do


      MDutton%NumMdispersives=conta
      allocate (MDutton%Medium(1 : MDutton%NumMdispersives))
      conta=0
      do jmed=1,sgg%NumMedia
         if ((sgg%Med(jmed)%Is%Mdispersive).and.(.not.sgg%Med(jmed)%Is%MdispersiveANIS)) then
            conta=conta+1
            MDutton%Medium(conta)%indexmed=jmed !correspondence with the main medium
            MDutton%Medium(conta)%numpolres11=sgg%Med(jmed)%Mdispersive(1)%numpolres11
            allocate (MDutton%Medium(conta)%Beta (1 : sgg%Med(jmed)%Mdispersive(1)%numpolres11),&
            MDutton%Medium(conta)%Kappa(1 : sgg%Med(jmed)%Mdispersive(1)%numpolres11), &
            MDutton%Medium(conta)%GM3  (1 : sgg%Med(jmed)%Mdispersive(1)%numpolres11))
            MDutton%Medium(conta)%Beta (1 : sgg%Med(jmed)%Mdispersive(1)%numpolres11)=0.0_RKIND
            MDutton%Medium(conta)%Kappa(1 : sgg%Med(jmed)%Mdispersive(1)%numpolres11)=0.0_RKIND
            MDutton%Medium(conta)%GM3  (1 : sgg%Med(jmed)%Mdispersive(1)%numpolres11)=0.0_RKIND
            do i1=1,sgg%Med(jmed)%Mdispersive(1)%numpolres11
               MDutton%Medium(conta)%Kappa(i1) =(1.0_RKIND+sgg%Med(jmed)%Mdispersive(1)%a11(i1)*sgg%dt/2.0_RKIND)/&
               (1.0_RKIND-sgg%Med(jmed)%Mdispersive(1)%a11(i1)*sgg%dt/2.0_RKIND)
               MDutton%Medium(conta)%Beta(i1)=  (   sgg%Med(jmed)%Mdispersive(1)%c11(i1)*sgg%dt) /&
               (1.0_RKIND-sgg%Med(jmed)%Mdispersive(1)%a11(i1)*sgg%dt/2.0_RKIND)
            end do
         end if
      end do

      !calculate the coefficients
      do jmed=1,MDutton%NumMdispersives
         tempindex=MDutton%Medium(jmed)%indexmed
         numpolres=sgg%Med(tempindex)%Mdispersive(1)%numpolres11
         tempo=0.0_RKIND
         do i1=1,NumPolRes
            tempo=tempo+real(MDutton%Medium(jmed)%Beta(i1))
         end do
         GM1(tempindex)=        (2.0_RKIND * sgg%Med(tempindex)%Mdispersive(1)%mu11+tempo-sgg%Med(tempindex)%Mdispersive(1)%Sigmam11*sgg%dt)/ &
         (2.0_RKIND * sgg%Med(tempindex)%Mdispersive(1)%mu11+tempo+sgg%Med(tempindex)%Mdispersive(1)%Sigmam11*sgg%dt)
         GM2(tempindex)=(2.0_RKIND * sgg%dt)/(2.0_RKIND * sgg%Med(tempindex)%Mdispersive(1)%mu11+tempo+sgg%Med(tempindex)%Mdispersive(1)%Sigmam11*sgg%dt)
         do i1=1,NumPolRes
            MDutton%Medium(jmed)%GM3(i1)=GM2(tempindex)/2.0_RKIND * (1.0_RKIND+MDutton%Medium(jmed)%Kappa(i1))
         end do
      end do

      do jmed=1,MDutton%NumMdispersives
         tempindex=MDutton%Medium(jmed)%indexmed
         !!!Hx
         conta=0
         do k1=sgg%Sweep(IHX)%ZI,sgg%Sweep(IHX)%ZE
            do j1=sgg%Sweep(IHX)%YI,sgg%Sweep(IHX)%YE
               do i1=sgg%Sweep(IHX)%XI,sgg%Sweep(IHX)%XE
                  if ((media%sggMiHx(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do

         ThereAreMdispersives=ThereAreMdispersives.or.(conta /=0)
         MDutton%Medium(jmed)%NumNodesHx=conta
         allocate (MDutton%Medium(jmed)%NodesHx(1 : conta))
         do i1=1, conta
            allocate (MDutton%Medium(jmed)%NodesHx(i1)%Current(1 : sgg%Med(tempindex)%Mdispersive(1)%numpolres11))
         end do
         conta=0
         do k1=sgg%Sweep(IHX)%ZI,sgg%Sweep(IHX)%ZE
            do j1=sgg%Sweep(IHX)%YI,sgg%Sweep(IHX)%YE
               do i1=sgg%Sweep(IHX)%XI,sgg%Sweep(IHX)%XE
                  if ((media%sggMiHx(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     MDutton%Medium(jmed)%NodesHx(conta)%i=i1
                     MDutton%Medium(jmed)%NodesHx(conta)%j=j1
                     MDutton%Medium(jmed)%NodesHx(conta)%k=k1
                     MDutton%Medium(jmed)%NodesHx(conta)%WhatField=IHX
                     MDutton%Medium(jmed)%NodesHx(conta)%FieldPresent=>Hx(i1,j1,k1)
                  end if
               end do
            end do
         end do
         !!!Hy
         conta=0
         do k1=sgg%Sweep(IHY)%ZI,sgg%Sweep(IHY)%ZE
            do j1=sgg%Sweep(IHY)%YI,sgg%Sweep(IHY)%YE
               do i1=sgg%Sweep(IHY)%XI,sgg%Sweep(IHY)%XE
                  if ((media%sggMiHy(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do

         ThereAreMdispersives=ThereAreMdispersives.or.(conta /=0)
         MDutton%Medium(jmed)%NumNodesHy=conta
         allocate (MDutton%Medium(jmed)%NodesHy(1 : conta))
         do i1=1, conta
            allocate (MDutton%Medium(jmed)%NodesHy(i1)%Current(1 : sgg%Med(tempindex)%Mdispersive(1)%numpolres11))
         end do
         conta=0
         do k1=sgg%Sweep(IHY)%ZI,sgg%Sweep(IHY)%ZE
            do j1=sgg%Sweep(IHY)%YI,sgg%Sweep(IHY)%YE
               do i1=sgg%Sweep(IHY)%XI,sgg%Sweep(IHY)%XE
                  if ((media%sggMiHy(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     MDutton%Medium(jmed)%NodesHy(conta)%i=i1
                     MDutton%Medium(jmed)%NodesHy(conta)%j=j1
                     MDutton%Medium(jmed)%NodesHy(conta)%k=k1
                     MDutton%Medium(jmed)%NodesHy(conta)%WhatField=IHY
                     MDutton%Medium(jmed)%NodesHy(conta)%FieldPresent=>Hy(i1,j1,k1)
                  end if
               end do
            end do
         end do
         !!!Hz
         conta=0
         do k1=sgg%Sweep(IHZ)%ZI,sgg%Sweep(IHZ)%ZE
            do j1=sgg%Sweep(IHZ)%YI,sgg%Sweep(IHZ)%YE
               do i1=sgg%Sweep(IHZ)%XI,sgg%Sweep(IHZ)%XE
                  if ((media%sggMiHz(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do


         ThereAreMdispersives=ThereAreMdispersives.or.(conta /=0)
         MDutton%Medium(jmed)%NumNodesHz=conta
         allocate (MDutton%Medium(jmed)%NodesHz(1 : conta))
         do i1=1, conta
            allocate (MDutton%Medium(jmed)%NodesHz(i1)%Current(1 : sgg%Med(tempindex)%Mdispersive(1)%numpolres11))
         end do
         conta=0
         do k1=sgg%Sweep(IHZ)%ZI,sgg%Sweep(IHZ)%ZE
            do j1=sgg%Sweep(IHZ)%YI,sgg%Sweep(IHZ)%YE
               do i1=sgg%Sweep(IHZ)%XI,sgg%Sweep(IHZ)%XE
                  if ((media%sggMiHz(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     MDutton%Medium(jmed)%NodesHz(conta)%i=i1
                     MDutton%Medium(jmed)%NodesHz(conta)%j=j1
                     MDutton%Medium(jmed)%NodesHz(conta)%k=k1
                     MDutton%Medium(jmed)%NodesHz(conta)%WhatField=IHZ
                     MDutton%Medium(jmed)%NodesHz(conta)%FieldPresent=>Hz(i1,j1,k1)
                  end if
               end do
            end do
         end do
      end do

      !resume or start
      do jmed=1,MDutton%NumMdispersives
         numpolres=sgg%Med(MDutton%Medium(jmed)%indexmed)%Mdispersive(1)%numpolres11
         if (.not.resume) then
            !Hx,Jx
            do i1=1,MDutton%Medium(jmed)%NumNodesHx
               MDutton%Medium(jmed)%NodesHx(i1)%fieldPrevious=0.0_RKIND
               do k1=1,NumPolRes
                  MDutton%Medium(jmed)%NodesHx(i1)%current(k1)=0.0_RKIND
               end do
            end do
            !Hy,Jy
            do i1=1,MDutton%Medium(jmed)%NumNodesHy
               MDutton%Medium(jmed)%NodesHy(i1)%fieldPrevious=0.0_RKIND
               do k1=1,NumPolRes
                  MDutton%Medium(jmed)%NodesHy(i1)%current(k1)=0.0_RKIND
               end do
            end do

            !Hz,Jz
            do i1=1,MDutton%Medium(jmed)%NumNodesHz
               MDutton%Medium(jmed)%NodesHz(i1)%fieldPrevious=0.0_RKIND
               do k1=1,NumPolRes
                  MDutton%Medium(jmed)%NodesHz(i1)%current(k1)=0.0_RKIND
               end do
            end do
         else
            !Hx,Jx
            do i1=1,MDutton%Medium(jmed)%NumNodesHx
               read (14) MDutton%Medium(jmed)%NodesHx(i1)%fieldPrevious
               do k1=1,NumPolRes
                  read (14) MDutton%Medium(jmed)%NodesHx(i1)%current(k1)
               end do
            end do
            !Hy,Jy
            do i1=1,MDutton%Medium(jmed)%NumNodesHy
               read (14) MDutton%Medium(jmed)%NodesHy(i1)%fieldPrevious
               do k1=1,NumPolRes
                  read (14) MDutton%Medium(jmed)%NodesHy(i1)%current(k1)
               end do
            end do

            !Hz,Jz
            do i1=1,MDutton%Medium(jmed)%NumNodesHz
               read (14) MDutton%Medium(jmed)%NodesHz(i1)%fieldPrevious
               do k1=1,NumPolRes
                  read (14) MDutton%Medium(jmed)%NodesHz(i1)%current(k1)
               end do
            end do
         end if
      end do
      return
   end subroutine InitMdispersives

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! subroutine to advance the E field in the Mdispersive (no need to advance the magnetic field)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine AdvanceMdispersiveH(sgg)
      type(SGGFDTDINFO_t), intent(in)              :: sgg              ! Simulation data.
      !!!

      integer(kind=4) :: jmed,i1,k1,numpolres

      type(field_t), pointer  :: tempnode
      !!!

      do jmed=1,MDutton%NumMdispersives
         numpolres=MDutton%Medium(jmed)%numpolres11
         !Hx,Jx
         do i1=1,MDutton%Medium(jmed)%NumNodesHx
            tempnode=>MDutton%Medium(jmed)%NodesHx(i1)
            do k1=1,NumPolRes
               tempnode%fieldPresent=tempnode%FieldPresent-real(MDutton%Medium(jmed)%GM3(k1)*tempnode%current(k1))
            end do
            do k1=1,NumPolRes
               tempnode%current(k1)=MDutton%Medium(jmed)%Kappa(k1)  *tempnode%current(k1) + &
               MDutton%Medium(jmed)%Beta(k1)/sgg%dt*(tempnode%fieldPresent-tempnode%fieldPrevious)
            end do
            tempnode%fieldPrevious=tempnode%fieldPresent
            !stores previous field (careful, it is not a pointer but an assignment of values)
            !before the background algorithm starts computing it again
         end do
         !Hy,Jy
         do i1=1,MDutton%Medium(jmed)%NumNodesHy
            tempnode=>MDutton%Medium(jmed)%NodesHy(i1)
            do k1=1,NumPolRes
               tempnode%FieldPresent=tempnode%FieldPresent-real(MDutton%Medium(jmed)%GM3(k1)*tempnode%current(k1))
            end do
            do k1=1,NumPolRes
               tempnode%current(k1)=MDutton%Medium(jmed)%Kappa(k1)  *tempnode%current(k1)+ &
               MDutton%Medium(jmed)%Beta(k1)/sgg%dt*(tempnode%fieldPresent-tempnode%fieldPrevious)
            end do
            tempnode%fieldPrevious=tempnode%fieldPresent
         end do

         !Hz,Jz
         do i1=1,MDutton%Medium(jmed)%NumNodesHz
            tempnode=>MDutton%Medium(jmed)%NodesHz(i1)
            do k1=1,NumPolRes
               tempnode%FieldPresent=tempnode%FieldPresent-real(MDutton%Medium(jmed)%GM3(k1)*tempnode%current(k1))
            end do
            do k1=1,NumPolRes
               tempnode%current(k1)=MDutton%Medium(jmed)%Kappa(k1)   *tempnode%current(k1)+ &
               MDutton%Medium(jmed)%Beta(k1)/sgg%dt*(tempnode%fieldPresent-tempnode%fieldPrevious)
            end do
            tempnode%fieldPrevious=tempnode%fieldPresent
         end do



      end do

   end subroutine AdvanceMdispersiveH


   subroutine StoreFieldsMdispersives

      integer(kind=4) :: jmed,numpolres,i1,k1


      do jmed=1,MDutton%NumMdispersives
         numpolres=MDutton%Medium(jmed)%numpolres11
         !Hx,Jx
         do i1=1,MDutton%Medium(jmed)%NumNodesHx
            write(14,err=634) MDutton%Medium(jmed)%NodesHx(i1)%fieldPrevious
            do k1=1,NumPolRes
               write(14,err=634) MDutton%Medium(jmed)%NodesHx(i1)%current(k1)
            end do
         end do
         !Hy,Jy
         do i1=1,MDutton%Medium(jmed)%NumNodesHy
            write(14,err=634) MDutton%Medium(jmed)%NodesHy(i1)%fieldPrevious
            do k1=1,NumPolRes
               write(14,err=634) MDutton%Medium(jmed)%NodesHy(i1)%current(k1)
            end do
         end do

         !Hz,Jz
         do i1=1,MDutton%Medium(jmed)%NumNodesHz
            write(14,err=634) MDutton%Medium(jmed)%NodesHz(i1)%fieldPrevious
            do k1=1,NumPolRes
               write(14,err=634) MDutton%Medium(jmed)%NodesHz(i1)%current(k1)
            end do
         end do
      end do

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'MAGNETICDISPERSIVE: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsMdispersives

   subroutine DestroyMdispersives(sgg)
      type(SGGFDTDINFO_t), intent(inout) :: sgg

      integer(kind=4) :: jmed,i1,i


      !free up memory
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%Mdispersive).and.(.not.sgg%Med(i)%Is%PML).and.(.not.sgg%Med(i)%Is%MdispersiveANIS)) then
            deallocate(sgg%Med(i)%Mdispersive(1)%c11,sgg%Med(i)%Mdispersive(1)%a11)
         end if
      end do
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%Mdispersive).and.(.not.sgg%Med(i)%Is%PML).and.(.not.sgg%Med(i)%Is%MdispersiveANIS)) then
            deallocate(sgg%Med(i)%Mdispersive)
         end if
      end do

      do jmed=1,MDutton%NumMdispersives
         deallocate(MDutton%Medium(jmed)%Beta,MDutton%Medium(jmed)%Kappa,MDutton%Medium(jmed)%GM3)

         do i1=1,MDutton%Medium(jmed)%NumNodesHx
            deallocate(MDutton%Medium(jmed)%NodesHx(i1)%Current)
         end do
         deallocate(MDutton%Medium(jmed)%NodesHx)
         !
         do i1=1,MDutton%Medium(jmed)%NumNodesHy
            deallocate(MDutton%Medium(jmed)%NodesHy(i1)%Current)
         end do
         deallocate(MDutton%Medium(jmed)%NodesHy)
         !
         do i1=1,MDutton%Medium(jmed)%NumNodesHz
            deallocate(MDutton%Medium(jmed)%NodesHz(i1)%Current)
         end do
         deallocate(MDutton%Medium(jmed)%NodesHz)
      end do


      if (associated(MDutton%Medium))  deallocate(MDutton%Medium)

   end subroutine

end module Mdispersives_m
