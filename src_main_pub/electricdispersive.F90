

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module EDispersives !ojo los polos conjugados DEBEN APARECER explicitamente 20JUNE'12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Beware in Dutton's model BOTH pair OF complex conjugate poles/residues
! in input from .nfde MUST APPEAR (this is why the factor /2 in the algorithm part
! for instance a 1 real-pole and 2 ccomplex-pole material would require in nfde
! 5 poles (not 3) !UNTESTED SGG JUN'12
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module EDispersives


   use fdetypes
   USE REPORT
   implicit none
   private

   !structures needed by the EDispersive



   type field_t
      integer (kind=4)   ::  i,j,k
      integer (kind=4)   ::  WhatField

      REAL (KIND=RKIND), pointer                 ::  FieldPresent !apunta al campo del background
      REAL (KIND=RKIND)                          ::  FieldPrevious
      Complex (Kind=CKIND), pointer, dimension ( : )   ::  Current
   end type

   TYPE EDispersive_t
      integer (kind=4)  ::  indexmed,numnodesEx,numnodesEy,numnodesEz,numpolres11
      Complex (Kind=CKIND), pointer, dimension ( : )      ::  Beta,Kappa,G3
      type (field_t), pointer, dimension ( : )      ::   NodesEx,NodesEy,NodesEz
   END TYPE EDispersive_t


   type  EDispersive_t2
      integer (kind=4)  ::  NumEDispersives
      type (EDispersive_t), pointer, dimension( : )  ::  Medium
   end type
   type (EDispersive_t2) , save , target ::  Dutton


   public AdvanceEDispersiveE,InitEDispersives,StoreFieldsEDispersives,DestroyEDispersives

contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitEDispersives(sgg,media,ThereAreEDispersives,resume,g1,g2,ex,ey,ez)
   ! subroutine InitEDispersives(sgg,sggmiex,sggmiey,sggmiez,ThereAreEDispersives,resume,g1,g2,ex,ey,ez)
      type (SGGFDTDINFO), intent(IN)     ::  sgg
      type(media_matrices_t), intent(in) :: media
      !!!
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) ::  &
      sggMiEx(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE), &
      sggMiEy(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE), &
      sggMiEz(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE)
      REAL (KIND=RKIND)     , intent(inout)      ::  &
      G1(0 : sgg%NumMedia),G2(0 : sgg%NumMedia)
      REAL (KIND=RKIND)   , intent(inout), target      :: &
      Ex(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE),&
      Ey(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE),&
      Ez(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE)

      !habria que ir deprecando lo de pasar el EDispersive, etc.. porque hay acceso directo a sgg%Med%Dispersiv
      logical, INTENT(OUT)  ::  ThereAreEDispersives
      logical, INTENT(in)  ::  resume
      integer (kind=4)  ::  jmed,j1,conta,k1,i1,tempindex
      REAL (KIND=RKIND)   ::  tempo
      integer (kind=4)  ::  numpolres

      sggMiEx = media%sggMiEx
      sggMiEy = media%sggMiEy
      sggMiEz = media%sggMiEz

      ThereAreEDispersives=.FALSE.
      conta=0
      do jmed=1,sgg%NumMedia
         if ((sgg%Med(jmed)%Is%EDispersive).and.(.not.sgg%Med(jmed)%Is%EDispersiveAnis)) then
            conta=conta+1
         endif
      end do


      Dutton%NumEDispersives=conta
      allocate (Dutton%Medium(1 : Dutton%NumEDispersives))
      conta=0
      do jmed=1,sgg%NumMedia
         if ((sgg%Med(jmed)%Is%EDispersive).and.(.not.sgg%Med(jmed)%Is%EDispersiveAnis)) then
            conta=conta+1
            Dutton%Medium(conta)%indexmed=jmed !correspondencia con el medio principal
            Dutton%Medium(conta)%numpolres11=sgg%Med(jmed)%EDispersive(1)%numpolres11
            allocate (Dutton%Medium(conta)%Beta (1 : sgg%Med(jmed)%EDispersive(1)%numpolres11),&
            Dutton%Medium(conta)%Kappa(1 : sgg%Med(jmed)%EDispersive(1)%numpolres11), &
            Dutton%Medium(conta)%G3   (1 : sgg%Med(jmed)%EDispersive(1)%numpolres11))
            Dutton%Medium(conta)%Beta (1 : sgg%Med(jmed)%EDispersive(1)%numpolres11)=0.0_RKIND
            Dutton%Medium(conta)%Kappa(1 : sgg%Med(jmed)%EDispersive(1)%numpolres11)=0.0_RKIND
            Dutton%Medium(conta)%G3   (1 : sgg%Med(jmed)%EDispersive(1)%numpolres11)=0.0_RKIND
            do i1=1,sgg%Med(jmed)%EDispersive(1)%numpolres11
               Dutton%Medium(conta)%Kappa(i1) =(1.0_RKIND+sgg%Med(jmed)%EDispersive(1)%a11(i1)*sgg%dt/2.0_RKIND)/&
               (1.0_RKIND-sgg%Med(jmed)%EDispersive(1)%a11(i1)*sgg%dt/2.0_RKIND)
               Dutton%Medium(conta)%Beta(i1)=  (   sgg%Med(jmed)%EDispersive(1)%C11(i1)*sgg%dt) /&
               (1.0_RKIND-sgg%Med(jmed)%EDispersive(1)%a11(i1)*sgg%dt/2.0_RKIND)
            end do
         endif
      end do

      !calculate the coefficients
      do jmed=1,Dutton%NumEDispersives
         tempindex=Dutton%Medium(jmed)%indexmed
         numpolres=sgg%Med(tempindex)%EDispersive(1)%numpolres11
         tempo=0.0_RKIND
         Do i1=1,NumPolRes
            tempo=tempo+REAL (Dutton%Medium(jmed)%Beta(i1))
         end do
         G1(tempindex)=        (2.0_RKIND * sgg%Med(tempindex)%Edispersive(1)%eps11+tempo-sgg%Med(tempindex)%Edispersive(1)%Sigma11*sgg%dt)/ &
         (2.0_RKIND * sgg%Med(tempindex)%Edispersive(1)%eps11+tempo+sgg%Med(tempindex)%Edispersive(1)%Sigma11*sgg%dt)
         G2(tempindex)= 2.0_RKIND * sgg%dt/ (2.0_RKIND * sgg%Med(tempindex)%Edispersive(1)%eps11+tempo+sgg%Med(tempindex)%Edispersive(1)%Sigma11*sgg%dt)
         Do i1=1,NumPolRes
            Dutton%Medium(jmed)%G3(i1)=G2(tempindex)/2.0_RKIND * (1.0_RKIND+Dutton%Medium(jmed)%Kappa(i1))
         end do
      end do

      do jmed=1,Dutton%NumEDispersives
         tempindex=Dutton%Medium(jmed)%indexmed
         !!!Ex
         conta=0
         Do k1=sgg%Sweep(iEx)%ZI,sgg%Sweep(iEx)%ZE
            Do j1=sgg%Sweep(iEx)%YI,sgg%Sweep(iEx)%YE
               Do i1=sgg%Sweep(iEx)%XI,sgg%Sweep(iEx)%XE
                  if ((sggMiEx(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do

         ThereAreEdispersives=ThereAreEdispersives.or.(conta /=0)
         Dutton%Medium(jmed)%NumNodesEx=conta
         allocate (Dutton%Medium(jmed)%NodesEx(1 : conta))
         do i1=1, conta
            allocate (Dutton%Medium(jmed)%NodesEx(i1)%Current(1 : sgg%Med(tempindex)%EDispersive(1)%numpolres11))
         end do
         conta=0
         Do k1=sgg%Sweep(iEx)%ZI,sgg%Sweep(iEx)%ZE
            Do j1=sgg%Sweep(iEx)%YI,sgg%Sweep(iEx)%YE
               Do i1=sgg%Sweep(iEx)%XI,sgg%Sweep(iEx)%XE
                  if ((sggMiEx(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     Dutton%Medium(jmed)%NodesEx(conta)%i=i1
                     Dutton%Medium(jmed)%NodesEx(conta)%j=j1
                     Dutton%Medium(jmed)%NodesEx(conta)%k=k1
                     Dutton%Medium(jmed)%NodesEx(conta)%WhatField=iEx
                     Dutton%Medium(jmed)%NodesEx(conta)%FieldPresent=>Ex(i1,j1,k1)
                  endif
               end do
            end do
         end do
         !!!Ey
         conta=0
         Do k1=sgg%Sweep(iEy)%ZI,sgg%Sweep(iEy)%ZE
            Do j1=sgg%Sweep(iEy)%YI,sgg%Sweep(iEy)%YE
               Do i1=sgg%Sweep(iEy)%XI,sgg%Sweep(iEy)%XE
                  if ((sggMiEy(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do

         ThereAreEdispersives=ThereAreEdispersives.or.(conta /=0)
         Dutton%Medium(jmed)%NumNodesEy=conta
         allocate (Dutton%Medium(jmed)%NodesEy(1 : conta))
         do i1=1, conta
            allocate (Dutton%Medium(jmed)%NodesEy(i1)%Current(1 : sgg%Med(tempindex)%EDispersive(1)%numpolres11))
         end do
         conta=0
         Do k1=sgg%Sweep(iEy)%ZI,sgg%Sweep(iEy)%ZE
            Do j1=sgg%Sweep(iEy)%YI,sgg%Sweep(iEy)%YE
               Do i1=sgg%Sweep(iEy)%XI,sgg%Sweep(iEy)%XE
                  if ((sggMiEy(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     Dutton%Medium(jmed)%NodesEy(conta)%i=i1
                     Dutton%Medium(jmed)%NodesEy(conta)%j=j1
                     Dutton%Medium(jmed)%NodesEy(conta)%k=k1
                     Dutton%Medium(jmed)%NodesEy(conta)%WhatField=iEy
                     Dutton%Medium(jmed)%NodesEy(conta)%FieldPresent=>Ey(i1,j1,k1)
                  endif
               end do
            end do
         end do
         !!!Ez
         conta=0
         Do k1=sgg%Sweep(iEz)%ZI,sgg%Sweep(iEz)%ZE
            Do j1=sgg%Sweep(iEz)%YI,sgg%Sweep(iEz)%YE
               Do i1=sgg%Sweep(iEz)%XI,sgg%Sweep(iEz)%XE
                  if ((sggMiEz(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do


         ThereAreEdispersives=ThereAreEdispersives.or.(conta /=0)
         Dutton%Medium(jmed)%NumNodesEz=conta
         allocate (Dutton%Medium(jmed)%NodesEz(1 : conta))
         do i1=1, conta
            allocate (Dutton%Medium(jmed)%NodesEz(i1)%Current(1 : sgg%Med(tempindex)%EDispersive(1)%numpolres11))
         end do
         conta=0
         Do k1=sgg%Sweep(iEz)%ZI,sgg%Sweep(iEz)%ZE
            Do j1=sgg%Sweep(iEz)%YI,sgg%Sweep(iEz)%YE
               Do i1=sgg%Sweep(iEz)%XI,sgg%Sweep(iEz)%XE
                  if ((sggMiEz(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     Dutton%Medium(jmed)%NodesEz(conta)%i=i1
                     Dutton%Medium(jmed)%NodesEz(conta)%j=j1
                     Dutton%Medium(jmed)%NodesEz(conta)%k=k1
                     Dutton%Medium(jmed)%NodesEz(conta)%WhatField=iEz
                     Dutton%Medium(jmed)%NodesEz(conta)%FieldPresent=>Ez(i1,j1,k1)
                  endif
               end do
            end do
         end do
      end do

      !resume or start
      do jmed=1,Dutton%NumEDispersives
         numpolres=sgg%Med(Dutton%Medium(jmed)%indexmed)%EDispersive(1)%numpolres11
         if (.not.resume) then
            !Ex,Jx
            do i1=1,Dutton%Medium(jmed)%NumNodesEx
               Dutton%Medium(jmed)%NodesEx(i1)%fieldPrevious=0.0_RKIND
               Do k1=1,NumPolRes
                  Dutton%Medium(jmed)%NodesEx(i1)%current(k1)=0.0_RKIND
               enddo
            end do
            !Ey,Jy
            do i1=1,Dutton%Medium(jmed)%NumNodesEy
               Dutton%Medium(jmed)%NodesEy(i1)%fieldPrevious=0.0_RKIND
               Do k1=1,NumPolRes
                  Dutton%Medium(jmed)%NodesEy(i1)%current(k1)=0.0_RKIND
               enddo
            end do

            !Ez,Jz
            do i1=1,Dutton%Medium(jmed)%NumNodesEz
               Dutton%Medium(jmed)%NodesEz(i1)%fieldPrevious=0.0_RKIND
               Do k1=1,NumPolRes
                  Dutton%Medium(jmed)%NodesEz(i1)%current(k1)=0.0_RKIND
               enddo
            end do
         else
            !Ex,Jx
            do i1=1,Dutton%Medium(jmed)%NumNodesEx
               READ (14) Dutton%Medium(jmed)%NodesEx(i1)%fieldPrevious
               Do k1=1,NumPolRes
                  READ (14) Dutton%Medium(jmed)%NodesEx(i1)%current(k1)
               enddo
            end do
            !Ey,Jy
            do i1=1,Dutton%Medium(jmed)%NumNodesEy
               READ (14) Dutton%Medium(jmed)%NodesEy(i1)%fieldPrevious
               Do k1=1,NumPolRes
                  READ (14) Dutton%Medium(jmed)%NodesEy(i1)%current(k1)
               enddo
            end do

            !Ez,Jz
            do i1=1,Dutton%Medium(jmed)%NumNodesEz
               READ (14) Dutton%Medium(jmed)%NodesEz(i1)%fieldPrevious
               Do k1=1,NumPolRes
                  READ (14) Dutton%Medium(jmed)%NodesEz(i1)%current(k1)
               enddo
            end do
         endif
      end do
      return
   end subroutine InitEDispersives

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to advance the E field in the EDispersive (no need to advance the magnetic field)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine AdvanceEDispersiveE(sgg)

      type (SGGFDTDINFO), intent(IN)              :: sgg              ! Simulation data.
      !!!

      integer (kind=4)  ::  jmed,i1,k1,numpolres

      type (field_t), pointer  ::  tempnode
      !!!

      do jmed=1,Dutton%NumEDispersives
         numpolres=Dutton%Medium(jmed)%numpolres11
         !Ex,Jx
         do i1=1,Dutton%Medium(jmed)%NumNodesEx
            tempnode=>Dutton%Medium(jmed)%NodesEx(i1)
            Do k1=1,NumPolRes
               tempnode%fieldPresent=tempnode%FieldPresent-REAL (Dutton%Medium(jmed)%G3(k1)*tempnode%current(k1) )
            enddo
            Do k1=1,NumPolRes
               tempnode%current(k1)=Dutton%Medium(jmed)%Kappa(k1)  *tempnode%current(k1) + &
               Dutton%Medium(jmed)%Beta(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious) /sgg%dt
            enddo
            tempnode%fieldPrevious=tempnode%fieldPresent
            !stores previous field (cuidado no es un apuntamiento sino una igualdad de valores)
            !antes de que re-empieze a calcularlo el algoritmo del background
         end do
         !Ey,Jy
         do i1=1,Dutton%Medium(jmed)%NumNodesEy
            tempnode=>Dutton%Medium(jmed)%NodesEy(i1)
            Do k1=1,NumPolRes
               tempnode%FieldPresent=tempnode%FieldPresent-REAL (Dutton%Medium(jmed)%G3(k1)*tempnode%current(k1))
            enddo


            Do k1=1,NumPolRes
               tempnode%current(k1)=Dutton%Medium(jmed)%Kappa(k1)  *tempnode%current(k1)+ &
               Dutton%Medium(jmed)%Beta(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious) /sgg%dt
            enddo
            tempnode%fieldPrevious=tempnode%fieldPresent
         end do

         !Ez,Jz
         do i1=1,Dutton%Medium(jmed)%NumNodesEz
            tempnode=>Dutton%Medium(jmed)%NodesEz(i1)
            Do k1=1,NumPolRes
               tempnode%FieldPresent=tempnode%FieldPresent-REAL (Dutton%Medium(jmed)%G3(k1)*tempnode%current(k1))
            enddo
            Do k1=1,NumPolRes
               tempnode%current(k1)=Dutton%Medium(jmed)%Kappa(k1)   *tempnode%current(k1)+ &
               Dutton%Medium(jmed)%Beta(k1)*(tempnode%fieldPresent-tempnode%fieldPrevious)  /sgg%dt
            enddo
            tempnode%fieldPrevious=tempnode%fieldPresent
         end do



      end do

   end subroutine AdvanceEDispersiveE


   subroutine StoreFieldsEDispersives

      integer (kind=4)  ::  jmed,numpolres,i1,k1


      do jmed=1,Dutton%NumEDispersives
         numpolres=Dutton%Medium(jmed)%numpolres11
         !Ex,Jx
         do i1=1,Dutton%Medium(jmed)%NumNodesEx
            write(14,err=634) Dutton%Medium(jmed)%NodesEx(i1)%fieldPrevious
            Do k1=1,NumPolRes
               write(14,err=634) Dutton%Medium(jmed)%NodesEx(i1)%current(k1)
            enddo
         end do
         !Ey,Jy
         do i1=1,Dutton%Medium(jmed)%NumNodesEy
            write(14,err=634) Dutton%Medium(jmed)%NodesEy(i1)%fieldPrevious
            Do k1=1,NumPolRes
               write(14,err=634) Dutton%Medium(jmed)%NodesEy(i1)%current(k1)
            enddo
         end do

         !Ez,Jz
         do i1=1,Dutton%Medium(jmed)%NumNodesEz
            write(14,err=634) Dutton%Medium(jmed)%NodesEz(i1)%fieldPrevious
            Do k1=1,NumPolRes
               write(14,err=634) Dutton%Medium(jmed)%NodesEz(i1)%current(k1)
            enddo
         end do
      end do
    
      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'ELECTRICDISPERSIVE: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine StoreFieldsEDispersives

   subroutine DestroyEDispersives(sgg)
      type (SGGFDTDINFO), intent(INOUT)         ::  sgg

      integer (kind=4)  ::  jmed,i1,i


      !free up memory
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%EDispersive).and.(.not.sgg%Med(i)%Is%PML).and.(.not.sgg%Med(i)%Is%EDispersiveAnis))  then
            deallocate (sgg%Med(i)%EDispersive(1)%C11,sgg%Med(i)%EDispersive(1)%a11)
         endif
      end do
      do i=1,sgg%NumMedia
         if ((sgg%Med(i)%Is%EDispersive).and.(.not.sgg%Med(i)%Is%PML).and.(.not.sgg%Med(i)%Is%EDispersiveAnis))  then
            deallocate (sgg%Med(i)%EDispersive)
         endif
      end do

      do jmed=1,Dutton%NumEDispersives
         deallocate (Dutton%Medium(jmed)%Beta,Dutton%Medium(jmed)%Kappa,Dutton%Medium(jmed)%G3)

         do i1=1,Dutton%Medium(jmed)%NumNodesEx
            deallocate (Dutton%Medium(jmed)%NodesEx(i1)%Current)
         end do
         deallocate (Dutton%Medium(jmed)%NodesEx)
         !
         do i1=1,Dutton%Medium(jmed)%NumNodesEy
            deallocate (Dutton%Medium(jmed)%NodesEy(i1)%Current)
         end do
         deallocate (Dutton%Medium(jmed)%NodesEy)
         !
         do i1=1,Dutton%Medium(jmed)%NumNodesEz
            deallocate (Dutton%Medium(jmed)%NodesEz(i1)%Current)
         end do
         deallocate (Dutton%Medium(jmed)%NodesEz)
      end do


      if (associated(Dutton%Medium)) deallocate (Dutton%Medium)

   end subroutine

end module EDispersives
