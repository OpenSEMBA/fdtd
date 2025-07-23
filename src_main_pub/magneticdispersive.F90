
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Mdispersives !ojo los polos conjugados DEBEN APARECER Hxplicitamente 20JUNE'12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Beware in MDutton's model BOTH pair OF Complex conjugate poles/residues
! in input from .nfde MUST APPEAR (this is why the factor /2 in the algorithm part
! for instance a 1 real-pole and 2 cComplex-pole material would require in nfde
! 5 poles (not 3) !UNTESTED SGG JUN'12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Mdispersives

   !mismo switch electrico y magnetico

   use fdetypes
   USE REPORT
   implicit none
   private

   !structures needed by the Mdispersive



   type field_t
      integer (kind=4)   ::  i,j,k
      integer (kind=4)   ::  WhatField
      REAL (KIND=RKIND), pointer                 ::  FieldPresent !apunta al campo del background
      REAL (KIND=RKIND)                          ::  FieldPrevious
      Complex (Kind=CKIND), pointer, dimension ( : )   ::  Current
   end type

   TYPE Mdispersive_t
      integer (kind=4)  ::  indexmed,numnodesHx,numnodesHy,numnodesHz,numpolres11
      Complex (Kind=CKIND), pointer, dimension ( : )      ::  Beta,Kappa,GM3
      type (field_t), pointer, dimension ( : )      ::   NodesHx,NodesHy,NodesHz
   END TYPE Mdispersive_t


   type  Mdispersive_t2
      integer (kind=4)  ::  NumMdispersives
      type (Mdispersive_t), pointer, dimension( : )  ::  Medium
   end type

   !!!LOCAL VARIABLES
   type (Mdispersive_t2) , save ::  MDutton


   public AdvanceMdispersiveH,InitMdispersives,StoreFieldsMdispersives,DestroyMdispersives

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMdispersives(sgg,media,ThereAreMdispersives,resume,GM1,GM2,Hx,Hy,Hz)
   ! subroutine InitMdispersives(sgg,sggmiHx,sggmiHy,sggmiHz,ThereAreMdispersives,resume,GM1,GM2,Hx,Hy,Hz)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      type(media_matrices_t), intent(in) :: media
      !!!
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) ::  &
      sggMiHx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE), &
      sggMiHy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE), &
      sggMiHz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)
      REAL (KIND=RKIND)     , intent(inout)      ::  &
      GM1(0 : sgg%NumMedia),GM2(0 : sgg%NumMedia)
      REAL (KIND=RKIND)   , intent(inout), target      :: &
      Hx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE),&
      Hy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE),&
      Hz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)

      !habria que ir deprecando lo de pasar el Mdispersive, etc.. porque hay acceso directo a sgg%Med%Dispersiv
      logical, INTENT(OUT)  ::  ThereAreMdispersives
      logical, INTENT(in)  ::  resume
      integer (kind=4)  ::  jmed,j1,conta,k1,i1,tempindex
      REAL (KIND=RKIND)   ::  tempo
      integer (kind=4)  ::  numpolres
      MDutton%Medium => null()
      !!!
      sggMiHx = media%sggMiHx
      sggMiHy = media%sggMiHy
      sggMiHz = media%sggMiHz

      ThereAreMdispersives=.FALSE.
      conta=0
      do jmed=1,sgg%NumMedia
         if ((sgg%Med(jmed)%Is%Mdispersive).and.(.not.sgg%Med(jmed)%Is%MdispersiveANIS)) then
            conta=conta+1
         endif
      end do


      MDutton%NumMdispersives=conta
      allocate (MDutton%Medium(1 : MDutton%NumMdispersives))
      conta=0
      do jmed=1,sgg%NumMedia
         if ((sgg%Med(jmed)%Is%Mdispersive).and.(.not.sgg%Med(jmed)%Is%MdispersiveANIS)) then
            conta=conta+1
            MDutton%Medium(conta)%indexmed=jmed !correspondencia con el medio principal
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
         endif
      end do

      !calculate the coefficients
      do jmed=1,MDutton%NumMdispersives
         tempindex=MDutton%Medium(jmed)%indexmed
         numpolres=sgg%Med(tempindex)%Mdispersive(1)%numpolres11
         tempo=0.0_RKIND
         Do i1=1,NumPolRes
            tempo=tempo+REAL (MDutton%Medium(jmed)%Beta(i1))
         end do
         GM1(tempindex)=        (2.0_RKIND * sgg%Med(tempindex)%Mdispersive(1)%mu11+tempo-sgg%Med(tempindex)%Mdispersive(1)%Sigmam11*sgg%dt)/ &
         (2.0_RKIND * sgg%Med(tempindex)%Mdispersive(1)%mu11+tempo+sgg%Med(tempindex)%Mdispersive(1)%Sigmam11*sgg%dt)
         GM2(tempindex)=(2.0_RKIND * sgg%dt)/(2.0_RKIND * sgg%Med(tempindex)%Mdispersive(1)%mu11+tempo+sgg%Med(tempindex)%Mdispersive(1)%Sigmam11*sgg%dt)
         Do i1=1,NumPolRes
            MDutton%Medium(jmed)%GM3(i1)=GM2(tempindex)/2.0_RKIND * (1.0_RKIND+MDutton%Medium(jmed)%Kappa(i1))
         end do
      end do

      do jmed=1,MDutton%NumMdispersives
         tempindex=MDutton%Medium(jmed)%indexmed
         !!!Hx
         conta=0
         Do k1=sgg%Sweep(iHx)%ZI,sgg%Sweep(iHx)%ZE
            Do j1=sgg%Sweep(iHx)%YI,sgg%Sweep(iHx)%YE
               Do i1=sgg%Sweep(iHx)%XI,sgg%Sweep(iHx)%XE
                  if ((sggMiHx(i1,j1,k1)) == tempindex)  conta=conta+1
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
         Do k1=sgg%Sweep(iHx)%ZI,sgg%Sweep(iHx)%ZE
            Do j1=sgg%Sweep(iHx)%YI,sgg%Sweep(iHx)%YE
               Do i1=sgg%Sweep(iHx)%XI,sgg%Sweep(iHx)%XE
                  if ((sggMiHx(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     MDutton%Medium(jmed)%NodesHx(conta)%i=i1
                     MDutton%Medium(jmed)%NodesHx(conta)%j=j1
                     MDutton%Medium(jmed)%NodesHx(conta)%k=k1
                     MDutton%Medium(jmed)%NodesHx(conta)%WhatField=iHx
                     MDutton%Medium(jmed)%NodesHx(conta)%FieldPresent=>Hx(i1,j1,k1)
                  endif
               end do
            end do
         end do
         !!!Hy
         conta=0
         Do k1=sgg%Sweep(iHy)%ZI,sgg%Sweep(iHy)%ZE
            Do j1=sgg%Sweep(iHy)%YI,sgg%Sweep(iHy)%YE
               Do i1=sgg%Sweep(iHy)%XI,sgg%Sweep(iHy)%XE
                  if ((sggMiHy(i1,j1,k1)) == tempindex)  conta=conta+1
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
         Do k1=sgg%Sweep(iHy)%ZI,sgg%Sweep(iHy)%ZE
            Do j1=sgg%Sweep(iHy)%YI,sgg%Sweep(iHy)%YE
               Do i1=sgg%Sweep(iHy)%XI,sgg%Sweep(iHy)%XE
                  if ((sggMiHy(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     MDutton%Medium(jmed)%NodesHy(conta)%i=i1
                     MDutton%Medium(jmed)%NodesHy(conta)%j=j1
                     MDutton%Medium(jmed)%NodesHy(conta)%k=k1
                     MDutton%Medium(jmed)%NodesHy(conta)%WhatField=iHy
                     MDutton%Medium(jmed)%NodesHy(conta)%FieldPresent=>Hy(i1,j1,k1)
                  endif
               end do
            end do
         end do
         !!!Hz
         conta=0
         Do k1=sgg%Sweep(iHz)%ZI,sgg%Sweep(iHz)%ZE
            Do j1=sgg%Sweep(iHz)%YI,sgg%Sweep(iHz)%YE
               Do i1=sgg%Sweep(iHz)%XI,sgg%Sweep(iHz)%XE
                  if ((sggMiHz(i1,j1,k1)) == tempindex)  conta=conta+1
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
         Do k1=sgg%Sweep(iHz)%ZI,sgg%Sweep(iHz)%ZE
            Do j1=sgg%Sweep(iHz)%YI,sgg%Sweep(iHz)%YE
               Do i1=sgg%Sweep(iHz)%XI,sgg%Sweep(iHz)%XE
                  if ((sggMiHz(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     MDutton%Medium(jmed)%NodesHz(conta)%i=i1
                     MDutton%Medium(jmed)%NodesHz(conta)%j=j1
                     MDutton%Medium(jmed)%NodesHz(conta)%k=k1
                     MDutton%Medium(jmed)%NodesHz(conta)%WhatField=iHz
                     MDutton%Medium(jmed)%NodesHz(conta)%FieldPresent=>Hz(i1,j1,k1)
                  endif
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
               Do k1=1,NumPolRes
                  MDutton%Medium(jmed)%NodesHx(i1)%current(k1)=0.0_RKIND
               enddo
            end do
            !Hy,Jy
            do i1=1,MDutton%Medium(jmed)%NumNodesHy
               MDutton%Medium(jmed)%NodesHy(i1)%fieldPrevious=0.0_RKIND
               Do k1=1,NumPolRes
                  MDutton%Medium(jmed)%NodesHy(i1)%current(k1)=0.0_RKIND
               enddo
            end do

            !Hz,Jz
            do i1=1,MDutton%Medium(jmed)%NumNodesHz
               MDutton%Medium(jmed)%NodesHz(i1)%fieldPrevious=0.0_RKIND
               Do k1=1,NumPolRes
                  MDutton%Medium(jmed)%NodesHz(i1)%current(k1)=0.0_RKIND
               enddo
            end do
         else
            !Hx,Jx
            do i1=1,MDutton%Medium(jmed)%NumNodesHx
               READ (14) MDutton%Medium(jmed)%NodesHx(i1)%fieldPrevious
               Do k1=1,NumPolRes
                  READ (14) MDutton%Medium(jmed)%NodesHx(i1)%current(k1)
               enddo
            end do
            !Hy,Jy
            do i1=1,MDutton%Medium(jmed)%NumNodesHy
               READ (14) MDutton%Medium(jmed)%NodesHy(i1)%fieldPrevious
               Do k1=1,NumPolRes
                  READ (14) MDutton%Medium(jmed)%NodesHy(i1)%current(k1)
               enddo
            end do

            !Hz,Jz
            do i1=1,MDutton%Medium(jmed)%NumNodesHz
               READ (14) MDutton%Medium(jmed)%NodesHz(i1)%fieldPrevious
               Do k1=1,NumPolRes
                  READ (14) MDutton%Medium(jmed)%NodesHz(i1)%current(k1)
               enddo
            end do
         endif
      end do
      return
   end subroutine InitMdispersives

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to advance the E field in the Mdispersive (no need to advance the magnetic field)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine AdvanceMdispersiveH(sgg)
      type (SGGFDTDINFO), intent(IN)              :: sgg              ! Simulation data.
      !!!

      integer (kind=4)  ::  jmed,i1,k1,numpolres

      type (field_t), pointer  ::  tempnode
      !!!

      do jmed=1,MDutton%NumMdispersives
         numpolres=MDutton%Medium(jmed)%numpolres11
         !Hx,Jx
         do i1=1,MDutton%Medium(jmed)%NumNodesHx
            tempnode=>MDutton%Medium(jmed)%NodesHx(i1)
            Do k1=1,NumPolRes
               tempnode%fieldPresent=tempnode%FieldPresent-REAL (MDutton%Medium(jmed)%GM3(k1)*tempnode%current(k1))
            enddo
            Do k1=1,NumPolRes
               tempnode%current(k1)=MDutton%Medium(jmed)%Kappa(k1)  *tempnode%current(k1) + &
               MDutton%Medium(jmed)%Beta(k1)/sgg%dt*(tempnode%fieldPresent-tempnode%fieldPrevious)
            enddo
            tempnode%fieldPrevious=tempnode%fieldPresent
            !stores previous field (cuidado no es un apuntamiento sino una igualdad de valores)
            !antes de que re-empiHze a calcularlo el algoritmo del background
         end do
         !Hy,Jy
         do i1=1,MDutton%Medium(jmed)%NumNodesHy
            tempnode=>MDutton%Medium(jmed)%NodesHy(i1)
            Do k1=1,NumPolRes
               tempnode%FieldPresent=tempnode%FieldPresent-REAL (MDutton%Medium(jmed)%GM3(k1)*tempnode%current(k1))
            enddo
            Do k1=1,NumPolRes
               tempnode%current(k1)=MDutton%Medium(jmed)%Kappa(k1)  *tempnode%current(k1)+ &
               MDutton%Medium(jmed)%Beta(k1)/sgg%dt*(tempnode%fieldPresent-tempnode%fieldPrevious)
            enddo
            tempnode%fieldPrevious=tempnode%fieldPresent
         end do

         !Hz,Jz
         do i1=1,MDutton%Medium(jmed)%NumNodesHz
            tempnode=>MDutton%Medium(jmed)%NodesHz(i1)
            Do k1=1,NumPolRes
               tempnode%FieldPresent=tempnode%FieldPresent-REAL (MDutton%Medium(jmed)%GM3(k1)*tempnode%current(k1))
            enddo
            Do k1=1,NumPolRes
               tempnode%current(k1)=MDutton%Medium(jmed)%Kappa(k1)   *tempnode%current(k1)+ &
               MDutton%Medium(jmed)%Beta(k1)/sgg%dt*(tempnode%fieldPresent-tempnode%fieldPrevious)
            enddo
            tempnode%fieldPrevious=tempnode%fieldPresent
         end do



      end do

   end subroutine AdvanceMdispersiveH


   subroutine StoreFieldsMdispersives

      integer (kind=4)  ::  jmed,numpolres,i1,k1


      do jmed=1,MDutton%NumMdispersives
         numpolres=MDutton%Medium(jmed)%numpolres11
         !Hx,Jx
         do i1=1,MDutton%Medium(jmed)%NumNodesHx
            write(14,err=634) MDutton%Medium(jmed)%NodesHx(i1)%fieldPrevious
            Do k1=1,NumPolRes
               write(14,err=634) MDutton%Medium(jmed)%NodesHx(i1)%current(k1)
            enddo
         end do
         !Hy,Jy
         do i1=1,MDutton%Medium(jmed)%NumNodesHy
            write(14,err=634) MDutton%Medium(jmed)%NodesHy(i1)%fieldPrevious
            Do k1=1,NumPolRes
               write(14,err=634) MDutton%Medium(jmed)%NodesHy(i1)%current(k1)
            enddo
         end do

         !Hz,Jz
         do i1=1,MDutton%Medium(jmed)%NumNodesHz
            write(14,err=634) MDutton%Medium(jmed)%NodesHz(i1)%fieldPrevious
            Do k1=1,NumPolRes
               write(14,err=634) MDutton%Medium(jmed)%NodesHz(i1)%current(k1)
            enddo
         end do
      end do

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'MAGNETICDISPERSIVE: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsMdispersives

   subroutine DestroyMdispersives(sgg)
      type (SGGFDTDINFO), intent(INout)         ::  sgg

      integer (kind=4)  ::  jmed,i1,i


      !free up memory
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%Mdispersive).and.(.not.sgg%Med(i)%Is%PML).and.(.not.sgg%Med(i)%Is%MdispersiveANIS)) then
            deallocate (sgg%Med(i)%Mdispersive(1)%c11,sgg%Med(i)%Mdispersive(1)%a11)
         endif
      end do
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%Mdispersive).and.(.not.sgg%Med(i)%Is%PML).and.(.not.sgg%Med(i)%Is%MdispersiveANIS)) then
            deallocate (sgg%Med(i)%Mdispersive)
         endif
      end do

      do jmed=1,MDutton%NumMdispersives
         deallocate (MDutton%Medium(jmed)%Beta,MDutton%Medium(jmed)%Kappa,MDutton%Medium(jmed)%GM3)

         do i1=1,MDutton%Medium(jmed)%NumNodesHx
            deallocate (MDutton%Medium(jmed)%NodesHx(i1)%Current)
         end do
         deallocate (MDutton%Medium(jmed)%NodesHx)
         !
         do i1=1,MDutton%Medium(jmed)%NumNodesHy
            deallocate (MDutton%Medium(jmed)%NodesHy(i1)%Current)
         end do
         deallocate (MDutton%Medium(jmed)%NodesHy)
         !
         do i1=1,MDutton%Medium(jmed)%NumNodesHz
            deallocate (MDutton%Medium(jmed)%NodesHz(i1)%Current)
         end do
         deallocate (MDutton%Medium(jmed)%NodesHz)
      end do


      if (associated(MDutton%Medium))  deallocate (MDutton%Medium)

   end subroutine

end module Mdispersives
