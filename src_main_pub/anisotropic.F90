
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Module Anisotropic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Beware in Anis's model only one of the pair complex conjugate poles/residues
! in input from .nfde (this is why the factor /2 in the parser part
! for instance a 1 real-pole and 2 ccomplex-pole material would require in nfde
! ONLY 2 poles (not 3) and provide the real one and any one of the couple of conjugates)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Anisotropic
   use fdetypes
   implicit none
   private

   !structures needed by the Anisotropic
   type, public::  Coeff_t
      real (kind=rkind) :: eexx,eexy,eexz,eeyx,eeyy,eeyz,eezx,eezy,eezz
      real (kind=rkind) :: ehxx,ehxy,ehxz,ehyx,ehyy,ehyz,ehzx,ehzy,ehzz
      real (kind=rkind) :: hexx,hexy,hexz,heyx,heyy,heyz,hezx,hezy,hezz
      real (kind=rkind) :: hhxx,hhxy,hhxz,hhyx,hhyy,hhyz,hhzx,hhzy,hhzz
   end type
   type, public::  LocalSharedElement_t
      integer (kind=4) :: times
      integer (kind=4), dimension(:), pointer ::  SharedMed
      !
      type (coeff_t) :: coeff
   end type

   TYPE, public :: Anisotropicinfo_t
      integer (kind=4)  ::  indexmed,numnodesEx,numnodesEy,numnodesEz
      integer (kind=4)  ::           numnodesHx,numnodesHy,numnodesHz
      integer (kind=4), dimension(:), pointer :: Ex_i,Ey_i,Ez_i,Hx_i,Hy_i,Hz_i
      integer (kind=4), dimension(:), pointer :: Ex_j,Ey_j,Ez_j,Hx_j,Hy_j,Hz_j
      integer (kind=4), dimension(:), pointer :: Ex_k,Ey_k,Ez_k,Hx_k,Hy_k,Hz_k
      real (kind=rkind), dimension(:), pointer :: Ex_value,Ey_value,Ez_value,Hx_value,Hy_value,Hz_value
      type (LocalSharedElement_t), dimension(:), pointer :: Ex_Shared,Ey_Shared,Ez_Shared
      type (LocalSharedElement_t), dimension(:), pointer :: Hx_Shared,Hy_Shared,Hz_Shared
      !
      type (coeff_t) :: coeff
      logical :: IsOnlyThinSlot
!
      REAL (KIND=RKIND),  DIMENSION(3,3)  ::  sigma,epr,mur,sigmaM  
   END TYPE Anisotropicinfo_t


   type, public ::  AnisotropicMed_t
      integer (kind=4)  ::  NumMed
      type (Anisotropicinfo_t), pointer, dimension( : )  ::  info
   end type
   type (AnisotropicMed_t),save, target :: AniMed

!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  eps0,mu0,cluz,zvac
!!!
   public AdvanceAnisotropicE,AdvanceAnisotropich,InitAnisotropic,DestroyAnisotropic,calc_anisotropicconstants


contains

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to initialize the parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitAnisotropic(sgg,media,ThereAreAnisotropic,ThereAreThinSlot,eps00,mu00)
   ! subroutine InitAnisotropic(sgg,sggmiex,sggmiey,sggmiez,sggMiHx ,sggMiHy ,sggMiHz,ThereAreAnisotropic,ThereAreThinSlot,eps00,mu00)
      REAL (KIND=RKIND)           ::  eps00,mu00
      type (SGGFDTDINFO), intent(IN) , target      ::  sgg
      type(media_matrices_t), intent(in) :: media
      !!!

      type (Anisotropic_t), pointer :: DummyAnisProp,dummyAnisShared

      logical, INTENT(OUT)  ::  ThereAreAnisotropic,ThereAreThinSlot
      integer (kind=4)  ::  jmed,j1,conta,k1,i1,tempindex

      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales

      !!!

      ThereAreAnisotropic=.FALSE.
      ThereAreThinSlot=.FALSE.
      conta=0
      do jmed=1,sgg%NumMedia
         if (sgg%Med(jmed)%Is%Anisotropic) then
            conta=conta+1
         endif
      end do


      AniMed%NumMed=conta
      allocate (AniMed%Info(1 : AniMed%NumMed))
      conta=0
      do jmed=1,sgg%NumMedia
         if (sgg%Med(jmed)%Is%Anisotropic) then
            conta=conta+1
            AniMed%info(conta)%indexmed=jmed !correspondencia con el medio principal
            if (sgg%Med(jmed)%Is%ThinSlot) then
               AniMed%info(conta)%IsOnlyThinSlot=.true.
            else
               AniMed%info(conta)%IsOnlyThinSlot=.false.
            endif
         endif
      end do

      do jmed=1,AniMed%NumMed
         AniMed%info(jmed)%sigma  = sgg%med(AniMed%Info(jmed)%indexmed)%Anisotropic(1)%sigma
         AniMed%info(jmed)%sigmam = sgg%med(AniMed%Info(jmed)%indexmed)%Anisotropic(1)%sigmam
         AniMed%info(jmed)%mur    = sgg%med(AniMed%Info(jmed)%indexmed)%Anisotropic(1)%mur
         AniMed%info(jmed)%epr    = sgg%med(AniMed%Info(jmed)%indexmed)%Anisotropic(1)%epr
      end do

      do jmed=1,AniMed%NumMed
         tempindex=AniMed%info(jmed)%indexmed
         !!!Ex
         conta=0
         Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
            Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
               Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
                  if ((media%sggMiEx(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do
         thereareanisotropic=thereareanisotropic.or.(conta /=0)
         thereareThinSlot=thereareanisotropic.and.AniMed%info(jmed)%IsOnlyThinSlot
         AniMed%info(jmed)%NumNodesEx=conta
         allocate (AniMed%info(jmed)%Ex_i(1 : conta))
         allocate (AniMed%info(jmed)%Ex_j(1 : conta))
         allocate (AniMed%info(jmed)%Ex_k(1 : conta))
         allocate (AniMed%info(jmed)%Ex_value(1 : conta))
         allocate (AniMed%info(jmed)%Ex_Shared(1 : conta))
         AniMed%info(jmed)%Ex_value=0.0_RKIND
         conta=0
         Do k1=sgg%SINPMLSweep(iEx)%ZI,sgg%SINPMLSweep(iEx)%ZE
            Do j1=sgg%SINPMLSweep(iEx)%YI,sgg%SINPMLSweep(iEx)%YE
               Do i1=sgg%SINPMLSweep(iEx)%XI,sgg%SINPMLSweep(iEx)%XE
                  if ((media%sggMiEx(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     AniMed%info(jmed)%Ex_Shared(conta)%times = 1
                     AniMed%info(jmed)%Ex_i(conta)=i1
                     AniMed%info(jmed)%Ex_j(conta)=j1
                     AniMed%info(jmed)%Ex_k(conta)=k1
                  endif
               end do
            end do
         end do
         !!!Ey
         conta=0
         Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
            Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
               Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
                  if ((media%sggMiEy(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do
         thereareanisotropic=thereareanisotropic.or.(conta /=0)
         thereareThinSlot=thereareanisotropic.and.AniMed%info(jmed)%IsOnlyThinSlot
         AniMed%info(jmed)%NumNodesEy=conta
         allocate (AniMed%info(jmed)%Ey_i(1 : conta))
         allocate (AniMed%info(jmed)%Ey_j(1 : conta))
         allocate (AniMed%info(jmed)%Ey_k(1 : conta))
         allocate (AniMed%info(jmed)%Ey_value(1 : conta))
         allocate (AniMed%info(jmed)%Ey_Shared(1 : conta))
         AniMed%info(jmed)%Ey_value=0.0_RKIND
         conta=0
         Do k1=sgg%SINPMLSweep(iEy)%ZI,sgg%SINPMLSweep(iEy)%ZE
            Do j1=sgg%SINPMLSweep(iEy)%YI,sgg%SINPMLSweep(iEy)%YE
               Do i1=sgg%SINPMLSweep(iEy)%XI,sgg%SINPMLSweep(iEy)%XE
                  if ((media%sggMiEy(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     AniMed%info(jmed)%Ey_Shared(conta)%times = 1
                     AniMed%info(jmed)%Ey_i(conta)=i1
                     AniMed%info(jmed)%Ey_j(conta)=j1
                     AniMed%info(jmed)%Ey_k(conta)=k1
                  endif
               end do
            end do
         end do
         !!!Ez
         conta=0
         Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
            Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
               Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
                  if ((media%sggMiEz(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do
         thereareanisotropic=thereareanisotropic.or.(conta /=0)
         thereareThinSlot=thereareanisotropic.and.AniMed%info(jmed)%IsOnlyThinSlot
         AniMed%info(jmed)%NumNodesEz=conta
         allocate (AniMed%info(jmed)%Ez_i(1 : conta))
         allocate (AniMed%info(jmed)%Ez_j(1 : conta))
         allocate (AniMed%info(jmed)%Ez_k(1 : conta))
         allocate (AniMed%info(jmed)%Ez_value(1 : conta))
         allocate (AniMed%info(jmed)%Ez_Shared(1 : conta))
         AniMed%info(jmed)%Ez_value=0.0_RKIND
         conta=0
         Do k1=sgg%SINPMLSweep(iEz)%ZI,sgg%SINPMLSweep(iEz)%ZE
            Do j1=sgg%SINPMLSweep(iEz)%YI,sgg%SINPMLSweep(iEz)%YE
               Do i1=sgg%SINPMLSweep(iEz)%XI,sgg%SINPMLSweep(iEz)%XE
                  if ((media%sggMiEz(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     AniMed%info(jmed)%Ez_Shared(conta)%times = 1
                     AniMed%info(jmed)%Ez_i(conta)=i1
                     AniMed%info(jmed)%Ez_j(conta)=j1
                     AniMed%info(jmed)%Ez_k(conta)=k1
                  endif
               end do
            end do
         end do
         !magneticos
         !!!Hx
         conta=0
         Do k1=sgg%SINPMLSweep(iHx)%ZI,sgg%SINPMLSweep(iHx)%ZE
            Do j1=sgg%SINPMLSweep(iHx)%YI,sgg%SINPMLSweep(iHx)%YE
               Do i1=sgg%SINPMLSweep(iHx)%XI,sgg%SINPMLSweep(iHx)%XE
                  if ((media%sggMiHx(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do
         thereareanisotropic=thereareanisotropic.or.(conta /=0)
         thereareThinSlot=thereareanisotropic.and.AniMed%info(jmed)%IsOnlyThinSlot
         AniMed%info(jmed)%NumNodesHx=conta
         allocate (AniMed%info(jmed)%Hx_i(1 : conta))
         allocate (AniMed%info(jmed)%Hx_j(1 : conta))
         allocate (AniMed%info(jmed)%Hx_k(1 : conta))
         allocate (AniMed%info(jmed)%Hx_value(1 : conta))
         allocate (AniMed%info(jmed)%Hx_Shared(1 : conta))
         AniMed%info(jmed)%Hx_value=0.0_RKIND
         conta=0
         Do k1=sgg%SINPMLSweep(iHx)%ZI,sgg%SINPMLSweep(iHx)%ZE
            Do j1=sgg%SINPMLSweep(iHx)%YI,sgg%SINPMLSweep(iHx)%YE
               Do i1=sgg%SINPMLSweep(iHx)%XI,sgg%SINPMLSweep(iHx)%XE
                  if ((media%sggMiHx(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     AniMed%info(jmed)%Hx_Shared(conta)%times = 1
                     AniMed%info(jmed)%Hx_i(conta)=i1
                     AniMed%info(jmed)%Hx_j(conta)=j1
                     AniMed%info(jmed)%Hx_k(conta)=k1
                  endif
               end do
            end do
         end do
         !!!Hy
         conta=0
         Do k1=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
            Do j1=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
               Do i1=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
                  if ((media%sggMiHy(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do
         thereareanisotropic=thereareanisotropic.or.(conta /=0)
         thereareThinSlot=thereareanisotropic.and.AniMed%info(jmed)%IsOnlyThinSlot
         AniMed%info(jmed)%NumNodesHy=conta
         allocate (AniMed%info(jmed)%Hy_i(1 : conta))
         allocate (AniMed%info(jmed)%Hy_j(1 : conta))
         allocate (AniMed%info(jmed)%Hy_k(1 : conta))
         allocate (AniMed%info(jmed)%Hy_value(1 : conta))
         allocate (AniMed%info(jmed)%Hy_Shared(1 : conta))
         AniMed%info(jmed)%Hy_value=0.0_RKIND
         conta=0
         Do k1=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
            Do j1=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
               Do i1=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
                  if ((media%sggMiHy(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     AniMed%info(jmed)%Hy_Shared(conta)%times = 1
                     AniMed%info(jmed)%Hy_i(conta)=i1
                     AniMed%info(jmed)%Hy_j(conta)=j1
                     AniMed%info(jmed)%Hy_k(conta)=k1
                  endif
               end do
            end do
         end do
         !!!Hz
         conta=0
         Do k1=sgg%SINPMLSweep(iHz)%ZI,sgg%SINPMLSweep(iHz)%ZE
            Do j1=sgg%SINPMLSweep(iHz)%YI,sgg%SINPMLSweep(iHz)%YE
               Do i1=sgg%SINPMLSweep(iHz)%XI,sgg%SINPMLSweep(iHz)%XE
                  if ((media%sggMiHz(i1,j1,k1)) == tempindex)  conta=conta+1
               end do
            end do
         end do
         thereareanisotropic=thereareanisotropic.or.(conta /=0)
         thereareThinSlot=thereareanisotropic.and.AniMed%info(jmed)%IsOnlyThinSlot
         AniMed%info(jmed)%NumNodesHz=conta
         allocate (AniMed%info(jmed)%Hz_i(1 : conta))
         allocate (AniMed%info(jmed)%Hz_j(1 : conta))
         allocate (AniMed%info(jmed)%Hz_k(1 : conta))
         allocate (AniMed%info(jmed)%Hz_value(1 : conta))
         allocate (AniMed%info(jmed)%Hz_Shared(1 : conta))
         AniMed%info(jmed)%Hz_value=0.0_RKIND
         conta=0
         Do k1=sgg%SINPMLSweep(iHz)%ZI,sgg%SINPMLSweep(iHz)%ZE
            Do j1=sgg%SINPMLSweep(iHz)%YI,sgg%SINPMLSweep(iHz)%YE
               Do i1=sgg%SINPMLSweep(iHz)%XI,sgg%SINPMLSweep(iHz)%XE
                  if ((media%sggMiHz(i1,j1,k1))==tempindex)  then
                     conta=conta+1
                     AniMed%info(jmed)%Hz_Shared(conta)%times = 1
                     AniMed%info(jmed)%Hz_i(conta)=i1
                     AniMed%info(jmed)%Hz_j(conta)=j1
                     AniMed%info(jmed)%Hz_k(conta)=k1
                  endif
               end do
            end do
         end do
      end do


      !actualiza el numero de shared

      do j1=1,sgg%Eshared%conta
         BuscaEx: do jmed=1,AniMed%NumMed !uno que encuentra ya contiene toda la informacion de times
            do i1=1,AniMed%info(jmed)%NumNodesEx
               if ((sgg%Eshared%elem(j1)%i == AniMed%info(jmed)%Ex_i(i1)).and. &
               (sgg%Eshared%elem(j1)%j == AniMed%info(jmed)%Ex_j(i1)).and. &
               (sgg%Eshared%elem(j1)%k == AniMed%info(jmed)%Ex_k(i1)).and. &
               (sgg%Eshared%elem(j1)%Field == iEx)) then
                  AniMed%info(jmed)%Ex_Shared(i1)%times = sgg%Eshared%elem(j1)%Times
                  if (sgg%Eshared%elem(j1)%Times > 1 ) then
                     allocate (AniMed%info(jmed)%Ex_Shared(i1)%SharedMed(1:sgg%Eshared%elem(j1)%Times)) !OJO, EN EL ORIGINAL LO ALLOCATEABA HASTA -1
                     !LO CAMBIO EL 26/10 POR ERROR EN LA NACELLE DEL SS5, PERO NO ESTOY SEGURO QUE HAYA UN BUG CON ESTO DE LOS SHARED !!?
                  endif
                  !PRINT *,'---> eX',sgg%Eshared%elem(j1)%i,sgg%Eshared%elem(j1)%J,sgg%Eshared%elem(j1)%k,sgg%Eshared%elem(j1)%times
                  exit buscaEx
               endif
            end do
         end do BuscaEx
         !
         BuscaEy: do jmed=1,AniMed%NumMed
            do i1=1,AniMed%info(jmed)%NumNodesEy
               if ((sgg%Eshared%elem(j1)%i == AniMed%info(jmed)%Ey_i(i1)).and. &
               (sgg%Eshared%elem(j1)%j == AniMed%info(jmed)%Ey_j(i1)).and. &
               (sgg%Eshared%elem(j1)%k == AniMed%info(jmed)%Ey_k(i1)).and. &
               (sgg%Eshared%elem(j1)%Field == iEy)) then
                  AniMed%info(jmed)%Ey_Shared(i1)%times = sgg%Eshared%elem(j1)%Times
                  if (sgg%Eshared%elem(j1)%Times > 1 ) then
                     allocate (AniMed%info(jmed)%Ey_Shared(i1)%SharedMed(1:sgg%Eshared%elem(j1)%Times))
                  endif
                  !PRINT *,'---> ey',sgg%Eshared%elem(j1)%i,sgg%Eshared%elem(j1)%J,sgg%Eshared%elem(j1)%k,sgg%Eshared%elem(j1)%times
                  exit buscaEy
               endif
            end do
         end do BuscaEy
         !
         BuscaEz: do jmed=1,AniMed%NumMed
            do i1=1,AniMed%info(jmed)%NumNodesEz
               if ((sgg%Eshared%elem(j1)%i == AniMed%info(jmed)%Ez_i(i1)).and. &
               (sgg%Eshared%elem(j1)%j == AniMed%info(jmed)%Ez_j(i1)).and. &
               (sgg%Eshared%elem(j1)%k == AniMed%info(jmed)%Ez_k(i1)).and. &
               (sgg%Eshared%elem(j1)%Field == iEz)) then
                  AniMed%info(jmed)%Ez_Shared(i1)%times = sgg%Eshared%elem(j1)%Times
                  if (sgg%Eshared%elem(j1)%Times > 1 ) then
                     allocate (AniMed%info(jmed)%Ez_Shared(i1)%SharedMed(1:sgg%Eshared%elem(j1)%Times))
                  endif
                  !PRINT *,'---> ez',sgg%Eshared%elem(j1)%i,sgg%Eshared%elem(j1)%J,sgg%Eshared%elem(j1)%k,sgg%Eshared%elem(j1)%times
                  exit buscaez
               endif
            end do
         end do BuscaEz
         !
      End do
      do j1=1,sgg%Hshared%conta
         BuscaHx: do jmed=1,AniMed%NumMed
            do i1=1,AniMed%info(jmed)%NumNodesHx
               if ((sgg%Hshared%elem(j1)%i == AniMed%info(jmed)%Hx_i(i1)).and. &
               (sgg%Hshared%elem(j1)%j == AniMed%info(jmed)%Hx_j(i1)).and. &
               (sgg%Hshared%elem(j1)%k == AniMed%info(jmed)%Hx_k(i1)).and. &
               (sgg%Hshared%elem(j1)%Field == iHx)) then
                  AniMed%info(jmed)%Hx_Shared(i1)%times = sgg%Hshared%elem(j1)%Times
                  if (sgg%Hshared%elem(j1)%Times > 1 ) then
                     allocate (AniMed%info(jmed)%Hx_Shared(i1)%SharedMed(1:sgg%Hshared%elem(j1)%Times))
                  endif
                  !PRINT *,'---> Hx',sgg%Hshared%elem(j1)%i,sgg%Hshared%elem(j1)%J,sgg%Hshared%elem(j1)%k,sgg%Hshared%elem(j1)%times
                  Exit buscaHx
               endif
            end do
         end do BuscaHx
         !
         BuscaHy: do jmed=1,AniMed%NumMed
            do i1=1,AniMed%info(jmed)%NumNodesHy
               if ((sgg%Hshared%elem(j1)%i == AniMed%info(jmed)%Hy_i(i1)).and. &
               (sgg%Hshared%elem(j1)%j == AniMed%info(jmed)%Hy_j(i1)).and. &
               (sgg%Hshared%elem(j1)%k == AniMed%info(jmed)%Hy_k(i1)).and. &
               (sgg%Hshared%elem(j1)%Field == iHy)) then
                  AniMed%info(jmed)%Hy_Shared(i1)%times = sgg%Hshared%elem(j1)%Times
                  if (sgg%Hshared%elem(j1)%Times > 1 ) then
                     allocate (AniMed%info(jmed)%Hy_Shared(i1)%SharedMed(1:sgg%Hshared%elem(j1)%Times))
                  endif
                  !PRINT *,'---> Hy',sgg%Hshared%elem(j1)%i,sgg%Hshared%elem(j1)%J,sgg%Hshared%elem(j1)%k,sgg%Hshared%elem(j1)%times
                  Exit buscaHy
               endif
            end do
         end do BuscaHy
         !
         BuscaHz: do jmed=1,AniMed%NumMed
            do i1=1,AniMed%info(jmed)%NumNodesHz
               if ((sgg%Hshared%elem(j1)%i == AniMed%info(jmed)%Hz_i(i1)).and. &
               (sgg%Hshared%elem(j1)%j == AniMed%info(jmed)%Hz_j(i1)).and. &
               (sgg%Hshared%elem(j1)%k == AniMed%info(jmed)%Hz_k(i1)).and. &
               (sgg%Hshared%elem(j1)%Field == iHz)) then
                  AniMed%info(jmed)%Hz_Shared(i1)%times = sgg%Hshared%elem(j1)%Times
                  if (sgg%Hshared%elem(j1)%Times > 1 ) then
                     allocate (AniMed%info(jmed)%Hz_Shared(i1)%SharedMed(1:sgg%Hshared%elem(j1)%Times))
                  endif
                  !PRINT *,'---> Hz',sgg%Hshared%elem(j1)%i,sgg%Hshared%elem(j1)%J,sgg%Hshared%elem(j1)%k,sgg%Hshared%elem(j1)%times
                  Exit buscaHz
               endif
            end do
         end do BuscaHz
         !
      end do

      !STORE THE INDEXES OF THE SHARED MEDIA


      do jmed=1,AniMed%NumMed !barrelos ahora todos
         do i1=1,AniMed%info(jmed)%NumNodesEx
            conta=0
            do j1=1,sgg%Eshared%conta
               if ((sgg%Eshared%elem(j1)%i == AniMed%info(jmed)%Ex_i(i1)).and. &
               (sgg%Eshared%elem(j1)%j == AniMed%info(jmed)%Ex_j(i1)).and. &
               (sgg%Eshared%elem(j1)%k == AniMed%info(jmed)%Ex_k(i1)).and. &
               (sgg%Eshared%elem(j1)%Field == iEx)) then
                  conta=conta+1
                  AniMed%info(jmed)%Ex_Shared(i1)%SharedMed(conta) =  sgg%Eshared%elem(j1)%SharedMed
               endif
            end do
         end do
      end do
      !
      do jmed=1,AniMed%NumMed
         do i1=1,AniMed%info(jmed)%NumNodesEy
            conta=0
            do j1=1,sgg%Eshared%conta
               if ((sgg%Eshared%elem(j1)%i == AniMed%info(jmed)%Ey_i(i1)).and. &
               (sgg%Eshared%elem(j1)%j == AniMed%info(jmed)%Ey_j(i1)).and. &
               (sgg%Eshared%elem(j1)%k == AniMed%info(jmed)%Ey_k(i1)).and. &
               (sgg%Eshared%elem(j1)%Field == iEy)) then
                  conta=conta+1
                  AniMed%info(jmed)%Ey_Shared(i1)%SharedMed(conta) =  sgg%Eshared%elem(j1)%SharedMed
                  continue
               endif
            end do
         end do
      end do
      !
      do jmed=1,AniMed%NumMed
         do i1=1,AniMed%info(jmed)%NumNodesEz
            conta=0
            do j1=1,sgg%Eshared%conta
               if ((sgg%Eshared%elem(j1)%i == AniMed%info(jmed)%Ez_i(i1)).and. &
               (sgg%Eshared%elem(j1)%j == AniMed%info(jmed)%Ez_j(i1)).and. &
               (sgg%Eshared%elem(j1)%k == AniMed%info(jmed)%Ez_k(i1)).and. &
               (sgg%Eshared%elem(j1)%Field == iEz)) then
                  conta=conta+1
                  AniMed%info(jmed)%Ez_Shared(i1)%SharedMed(conta) =  sgg%Eshared%elem(j1)%SharedMed
               endif
            end do
         end do
      end do
      !
      do jmed=1,AniMed%NumMed
         do i1=1,AniMed%info(jmed)%NumNodesHx
            conta = 0
            do j1=1,sgg%Hshared%conta
               if ((sgg%Hshared%elem(j1)%i == AniMed%info(jmed)%Hx_i(i1)).and. &
               (sgg%Hshared%elem(j1)%j == AniMed%info(jmed)%Hx_j(i1)).and. &
               (sgg%Hshared%elem(j1)%k == AniMed%info(jmed)%Hx_k(i1)).and. &
               (sgg%Hshared%elem(j1)%Field == iHx)) then
                  conta=conta+1
                  AniMed%info(jmed)%Hx_Shared(i1)%SharedMed(conta) =  sgg%Hshared%elem(j1)%SharedMed
                  continue
               endif
            end do
         end do
      end do
      !
      do jmed=1,AniMed%NumMed
         do i1=1,AniMed%info(jmed)%NumNodesHy
            conta = 0
            do j1=1,sgg%Hshared%conta
               if ((sgg%Hshared%elem(j1)%i == AniMed%info(jmed)%Hy_i(i1)).and. &
               (sgg%Hshared%elem(j1)%j == AniMed%info(jmed)%Hy_j(i1)).and. &
               (sgg%Hshared%elem(j1)%k == AniMed%info(jmed)%Hy_k(i1)).and. &
               (sgg%Hshared%elem(j1)%Field == iHy)) then
                  conta=conta+1
                  AniMed%info(jmed)%Hy_Shared(i1)%SharedMed(conta) =  sgg%Hshared%elem(j1)%SharedMed
               endif
            end do
         end do
      end do
      do jmed=1,AniMed%NumMed
         do i1=1,AniMed%info(jmed)%NumNodesHz
            conta = 0
            do j1=1,sgg%Hshared%conta
               if ((sgg%Hshared%elem(j1)%i == AniMed%info(jmed)%Hz_i(i1)).and. &
               (sgg%Hshared%elem(j1)%j == AniMed%info(jmed)%Hz_j(i1)).and. &
               (sgg%Hshared%elem(j1)%k == AniMed%info(jmed)%Hz_k(i1)).and. &
               (sgg%Hshared%elem(j1)%Field == iHz)) then
                  conta=conta+1
                  AniMed%info(jmed)%Hz_Shared(i1)%SharedMed(conta) =  sgg%Hshared%elem(j1)%SharedMed
               endif
            end do
         end do
      end do

      !crea las matrices



      do jmed=1,AniMed%NumMed !barrelos ahora todos
         dummyAnisProp => sgg%med( AniMed%Info(jmed)%indexmed )%Anisotropic(1)
         do i1=1,AniMed%info(jmed)%NumNodesEx
            conta=AniMed%info(jmed)%Ex_Shared(i1)%times
            if (conta > 1) then
               AniMed%info(jmed)%sigma  = dummyAnisProp%sigma  /conta
               AniMed%info(jmed)%sigmam = dummyAnisProp%sigmam /conta
               AniMed%info(jmed)%mur    = dummyAnisProp%mur    /conta
               AniMed%info(jmed)%epr    = dummyAnisProp%epr    /conta
               do k1=1,conta-1
                  dummyAnisShared => sgg%med( AniMed%info(jmed)%Ex_Shared(i1)%SharedMed(k1))%Anisotropic(1)
                  AniMed%info(jmed)%sigma  = AniMed%info(jmed)%sigma  + dummyAnisShared%sigma  /conta
                  AniMed%info(jmed)%sigmam = AniMed%info(jmed)%sigmam + dummyAnisShared%sigmam /conta
                  AniMed%info(jmed)%mur    = AniMed%info(jmed)%mur    + dummyAnisShared%mur    /conta
                  AniMed%info(jmed)%epr    = AniMed%info(jmed)%epr    + dummyAnisShared%epr    /conta
               end do
            endif
         end do


         do i1=1,AniMed%info(jmed)%NumNodesEy
            conta=AniMed%info(jmed)%Ey_Shared(i1)%times
            if (conta > 1) then
               AniMed%info(jmed)%sigma  = dummyAnisProp%sigma  /conta
               AniMed%info(jmed)%sigmam = dummyAnisProp%sigmam /conta
               AniMed%info(jmed)%mur    = dummyAnisProp%mur    /conta
               AniMed%info(jmed)%epr    = dummyAnisProp%epr    /conta
               do k1=1,conta-1
                  dummyAnisShared => sgg%med( AniMed%info(jmed)%Ey_Shared(i1)%SharedMed(k1))%Anisotropic(1)
                  AniMed%info(jmed)%sigma  = AniMed%info(jmed)%sigma  + dummyAnisShared%sigma  /conta
                  AniMed%info(jmed)%sigmam = AniMed%info(jmed)%sigmam + dummyAnisShared%sigmam /conta
                  AniMed%info(jmed)%mur    = AniMed%info(jmed)%mur    + dummyAnisShared%mur    /conta
                  AniMed%info(jmed)%epr    = AniMed%info(jmed)%epr    + dummyAnisShared%epr    /conta
               end do
            endif
         end do



         do i1=1,AniMed%info(jmed)%NumNodesEz
            conta=AniMed%info(jmed)%Ez_Shared(i1)%times
            if (conta > 1) then
               AniMed%info(jmed)%sigma  = dummyAnisProp%sigma  /conta
               AniMed%info(jmed)%sigmam = dummyAnisProp%sigmam /conta
               AniMed%info(jmed)%mur    = dummyAnisProp%mur    /conta
               AniMed%info(jmed)%epr    = dummyAnisProp%epr    /conta
               do k1=1,conta-1
                  dummyAnisShared => sgg%med( AniMed%info(jmed)%Ez_Shared(i1)%SharedMed(k1))%Anisotropic(1)
                  AniMed%info(jmed)%sigma  = AniMed%info(jmed)%sigma  + dummyAnisShared%sigma  /conta
                  AniMed%info(jmed)%sigmam = AniMed%info(jmed)%sigmam + dummyAnisShared%sigmam /conta
                  AniMed%info(jmed)%mur    = AniMed%info(jmed)%mur    + dummyAnisShared%mur    /conta
                  AniMed%info(jmed)%epr    = AniMed%info(jmed)%epr    + dummyAnisShared%epr    /conta
               end do
            endif
         end do



         do i1=1,AniMed%info(jmed)%NumNodesHx
            conta=AniMed%info(jmed)%Hx_Shared(i1)%times
            if (conta > 1) then
               AniMed%info(jmed)%sigma  = dummyAnisProp%sigma  /conta
               AniMed%info(jmed)%sigmam = dummyAnisProp%sigmam /conta
               AniMed%info(jmed)%mur    = dummyAnisProp%mur    /conta
               AniMed%info(jmed)%epr    = dummyAnisProp%epr    /conta
               do k1=1,conta-1
                  dummyAnisShared => sgg%med( AniMed%info(jmed)%Hx_Shared(i1)%SharedMed(k1))%Anisotropic(1)
                  AniMed%info(jmed)%sigma  = AniMed%info(jmed)%sigma  + dummyAnisShared%sigma  /conta
                  AniMed%info(jmed)%sigmam = AniMed%info(jmed)%sigmam + dummyAnisShared%sigmam /conta
                  AniMed%info(jmed)%mur    = AniMed%info(jmed)%mur    + dummyAnisShared%mur    /conta
                  AniMed%info(jmed)%epr    = AniMed%info(jmed)%epr    + dummyAnisShared%epr    /conta
               end do
            endif
         end do



         do i1=1,AniMed%info(jmed)%NumNodesHy
            conta=AniMed%info(jmed)%Hy_Shared(i1)%times
            if (conta > 1) then
               AniMed%info(jmed)%sigma  = dummyAnisProp%sigma  /conta
               AniMed%info(jmed)%sigmam = dummyAnisProp%sigmam /conta
               AniMed%info(jmed)%mur    = dummyAnisProp%mur    /conta
               AniMed%info(jmed)%epr    = dummyAnisProp%epr    /conta
               do k1=1,conta-1
                  dummyAnisShared => sgg%med( AniMed%info(jmed)%Hy_Shared(i1)%SharedMed(k1))%Anisotropic(1)
                  AniMed%info(jmed)%sigma  = AniMed%info(jmed)%sigma  + dummyAnisShared%sigma  /conta
                  AniMed%info(jmed)%sigmam = AniMed%info(jmed)%sigmam + dummyAnisShared%sigmam /conta
                  AniMed%info(jmed)%mur    = AniMed%info(jmed)%mur    + dummyAnisShared%mur    /conta
                  AniMed%info(jmed)%epr    = AniMed%info(jmed)%epr    + dummyAnisShared%epr    /conta
               end do
            endif
         end do



         do i1=1,AniMed%info(jmed)%NumNodesHz
            conta=AniMed%info(jmed)%Hz_Shared(i1)%times
            if (conta > 1) then
               AniMed%info(jmed)%sigma  = dummyAnisProp%sigma  /conta
               AniMed%info(jmed)%sigmam = dummyAnisProp%sigmam /conta
               AniMed%info(jmed)%mur    = dummyAnisProp%mur    /conta
               AniMed%info(jmed)%epr    = dummyAnisProp%epr    /conta
               do k1=1,conta-1
                  dummyAnisShared => sgg%med( AniMed%info(jmed)%Hz_Shared(i1)%SharedMed(k1))%Anisotropic(1)
                  AniMed%info(jmed)%sigma  = AniMed%info(jmed)%sigma  + dummyAnisShared%sigma  /conta
                  AniMed%info(jmed)%sigmam = AniMed%info(jmed)%sigmam + dummyAnisShared%sigmam /conta
                  AniMed%info(jmed)%mur    = AniMed%info(jmed)%mur    + dummyAnisShared%mur    /conta
                  AniMed%info(jmed)%epr    = AniMed%info(jmed)%epr    + dummyAnisShared%epr    /conta
               end do
               continue
            endif
         end do
      end do

      call calc_anisotropicconstants(sgg,eps0,mu0)

      return

   end subroutine InitAnisotropic

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to advance the E field in the Anisotropic (no need to advance the magnetic field)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine AdvanceAnisotropicE(sggalloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)

      type (XYZlimit_t), dimension (1:6), intent(in)                      ::  sggAlloc

      REAL (KIND=RKIND)   , intent(inout)      :: &
      Ex(sggalloc(iEx)%XI : sggalloc(iEx)%XE,sggalloc(iEx)%YI : sggalloc(iEx)%YE,sggalloc(iEx)%ZI : sggalloc(iEx)%ZE),&
      Ey(sggalloc(iEy)%XI : sggalloc(iEy)%XE,sggalloc(iEy)%YI : sggalloc(iEy)%YE,sggalloc(iEy)%ZI : sggalloc(iEy)%ZE),&
      Ez(sggalloc(iEz)%XI : sggalloc(iEz)%XE,sggalloc(iEz)%YI : sggalloc(iEz)%YE,sggalloc(iEz)%ZI : sggalloc(iEz)%ZE),&
      Hx(sggalloc(iHx)%XI : sggalloc(iHx)%XE,sggalloc(iHx)%YI : sggalloc(iHx)%YE,sggalloc(iHx)%ZI : sggalloc(iHx)%ZE),&
      Hy(sggalloc(iHy)%XI : sggalloc(iHy)%XE,sggalloc(iHy)%YI : sggalloc(iHy)%YE,sggalloc(iHy)%ZI : sggalloc(iHy)%ZE),&
      Hz(sggalloc(iHz)%XI : sggalloc(iHz)%XE,sggalloc(iHz)%YI : sggalloc(iHz)%YE,sggalloc(iHz)%ZI : sggalloc(iHz)%ZE)
      REAL (KIND=RKIND) , dimension (:)   , intent(in)      ::  &
      Idxe(sggalloc(iHx)%XI : sggalloc(iHx)%XE), &
      Idye(sggalloc(iHy)%YI : sggalloc(iHy)%YE), &
      Idze(sggalloc(iHz)%ZI : sggalloc(iHz)%ZE), &
      Idxh(sggalloc(iEx)%XI : sggalloc(iEx)%XE), &
      Idyh(sggalloc(iEy)%YI : sggalloc(iEy)%YE), &
      Idzh(sggalloc(iEz)%ZI : sggalloc(iEz)%ZE)
      !



      type (Coeff_t), pointer :: coeff
      type (Anisotropicinfo_t), pointer :: dummy
      integer (kind=4)  ::  jmed,i1,i,j,k

      do jmed=1,AniMed%NumMed
         dummy => AniMed%info(jmed)
         do i1=1,dummy%NumNodesEx
            if (Dummy%Ex_Shared(i1)%times == 1) then
               Coeff => AniMed%info(jmed)%coeff
            else
               Coeff => AniMed%info(jmed)%Ex_Shared(i1)%Coeff
            endif

            i=dummy%Ex_i(i1)
            j=dummy%Ex_j(i1)
            k=dummy%Ex_k(i1)
            !
            if  (dummy%IsOnlyThinSlot) then
               dummy%Ex_value(i1)= ex(i,j,k) - Coeff%ehxx*(hz(i,-1 + j,k) - hz(i,j,k))*Idyh(j) + &
               Coeff%ehxx*(hy(i,j,-1 + k) - hy(i,j,k))*Idzh(k)
            else
               dummy%Ex_value(i1)=                                                                                  &
               (4. * Coeff%eexx * ex(i,j,k) + Coeff%eexy * ey(i,-1 + j,k) + Coeff%eexy * ey(i,j,k) +                     &
               Coeff%eexy * ey(1 + i,-1 + j,k) + Coeff%eexy * ey(1 + i,j,k) + Coeff%eexz * ez(i,j,-1 + k) +        &
               Coeff%eexz * ez(i,j,k) + Coeff%eexz * ez(1 + i,j,-1 + k) + Coeff%eexz * ez(1 + i,j,k))/4. -         &
               ((Coeff%ehxz * hy(-1 + i,j,-1 + k) + Coeff%ehxz * hy(-1 + i,j,k) - Coeff%ehxz * hy(1 + i,j,-1 + k) -   &
               Coeff%ehxz * hy(1 + i,j,k) - Coeff%ehxy * hz(-1 + i,-1 + j,k) - Coeff%ehxy * hz(-1 + i,j,k) +     &
               Coeff%ehxy * hz(1 + i,-1 + j,k) + Coeff%ehxy * hz(1 + i,j,k)) * Idxe(i))/4. +                    &
               ((Coeff%ehxz * hx(i,-1 + j,-1 + k) + Coeff%ehxz * hx(i,-1 + j,k) - Coeff%ehxz * hx(i,j,-1 + k) -       &
               Coeff%ehxz * hx(i,j,k) + Coeff%ehxz * hx(1 + i,-1 + j,-1 + k) +                                &
               Coeff%ehxz * hx(1 + i,-1 + j,k) - Coeff%ehxz * hx(1 + i,j,-1 + k) - Coeff%ehxz * hx(1 + i,j,k) -  &
               4. * Coeff%ehxx * hz(i,-1 + j,k) + 4. * Coeff%ehxx * hz(i,j,k)) * Idyh(j))/4. -                        &
               ((Coeff%ehxy * hx(i,-1 + j,-1 + k) - Coeff%ehxy * hx(i,-1 + j,k) + Coeff%ehxy * hx(i,j,-1 + k) -       &
               Coeff%ehxy * hx(i,j,k) + Coeff%ehxy * hx(1 + i,-1 + j,-1 + k) -                                &
               Coeff%ehxy * hx(1 + i,-1 + j,k) + Coeff%ehxy * hx(1 + i,j,-1 + k) - Coeff%ehxy * hx(1 + i,j,k) -  &
               4. * Coeff%ehxx * hy(i,j,-1 + k) + 4. * Coeff%ehxx * hy(i,j,k)) * Idzh(k))/4.0_RKIND
            endif
         end do
         do i1=1,dummy%NumNodesEy
            if (Dummy%Ey_Shared(i1)%times == 1) then
               Coeff => AniMed%info(jmed)%coeff
            else
               Coeff => AniMed%info(jmed)%Ey_Shared(i1)%Coeff
            endif
            i=dummy%Ey_i(i1)
            j=dummy%Ey_j(i1)
            k=dummy%Ey_k(i1)
            !!
            if  (dummy%IsOnlyThinSlot) then
               dummy%Ey_value(i1)=ey(i,j,k) + Coeff%ehyy*(hz(-1 + i,j,k) - hz(i,j,k))*Idxh(i) -  &
               Coeff%ehyy*(hx(i,j,-1 + k) - hx(i,j,k))*Idzh(k)
            else
               dummy%Ey_value(i1)=                                                                                            &
               (Coeff%eeyx * ex(-1 + i,j,k) + Coeff%eeyx * ex(-1 + i,1 + j,k) + Coeff%eeyx * ex(i,j,k) +                      &
               Coeff%eeyx * ex(i,1 + j,k) + 4. * Coeff%eeyy * ey(i,j,k) + Coeff%eeyz * ez(i,j,-1 + k) +                       &
               Coeff%eeyz * ez(i,j,k) + Coeff%eeyz * ez(i,1 + j,-1 + k) + Coeff%eeyz * ez(i,1 + j,k))/4. -                 &
               ((Coeff%ehyz * hy(-1 + i,j,-1 + k) + Coeff%ehyz * hy(-1 + i,j,k) +                                          &
               Coeff%ehyz * hy(-1 + i,1 + j,-1 + k) + Coeff%ehyz * hy(-1 + i,1 + j,k) -                               &
               Coeff%ehyz * hy(i,j,-1 + k) - Coeff%ehyz * hy(i,j,k) - Coeff%ehyz * hy(i,1 + j,-1 + k) -                  &
               Coeff%ehyz * hy(i,1 + j,k) - 4. * Coeff%ehyy * hz(-1 + i,j,k) + 4. * Coeff%ehyy * hz(i,j,k)) * Idxh(i))/       &
               4. + ((Coeff%ehyz * hx(i,-1 + j,-1 + k) + Coeff%ehyz * hx(i,-1 + j,k) -                                    &
               Coeff%ehyz * hx(i,1 + j,-1 + k) - Coeff%ehyz * hx(i,1 + j,k) -                                         &
               Coeff%ehyx * hz(-1 + i,-1 + j,k) + Coeff%ehyx * hz(-1 + i,1 + j,k) -                                   &
               Coeff%ehyx * hz(i,-1 + j,k) + Coeff%ehyx * hz(i,1 + j,k)) * Idye(j))/4. -                                &
               ((4. * Coeff%ehyy * hx(i,j,-1 + k) - 4. * Coeff%ehyy * hx(i,j,k) - Coeff%ehyx * hy(-1 + i,j,-1 + k) +   &
               Coeff%ehyx * hy(-1 + i,j,k) - Coeff%ehyx * hy(-1 + i,1 + j,-1 + k) +                                   &
               Coeff%ehyx * hy(-1 + i,1 + j,k) - Coeff%ehyx * hy(i,j,-1 + k) + Coeff%ehyx * hy(i,j,k) -                  &
               Coeff%ehyx * hy(i,1 + j,-1 + k) + Coeff%ehyx * hy(i,1 + j,k)) * Idzh(k))/4.0_RKIND
            endif
         end do
         do i1=1,dummy%NumNodesEz
            if (Dummy%Ez_Shared(i1)%times == 1) then
               Coeff => AniMed%info(jmed)%coeff
            else
               Coeff => AniMed%info(jmed)%Ez_Shared(i1)%Coeff
            endif
            i=dummy%Ez_i(i1)
            j=dummy%Ez_j(i1)
            k=dummy%Ez_k(i1)
            !
            if  (dummy%IsOnlyThinSlot) then
               dummy%Ez_value(i1)= ez(i,j,k) - Coeff%ehzz*(hy(-1 + i,j,k) - hy(i,j,k))*Idxh(i) +  &
               Coeff%ehzz*(hx(i,-1 + j,k) - hx(i,j,k))*Idyh(j)
            else
               dummy%Ez_value(i1)=                                                                                        &
               (Coeff%eezx * ex(-1 + i,j,k) + Coeff%eezx * ex(-1 + i,j,1 + k) + Coeff%eezx * ex(i,j,k) +                    &
               Coeff%eezx * ex(i,j,1 + k) + Coeff%eezy * ey(i,-1 + j,k) + Coeff%eezy * ey(i,-1 + j,1 + k) +              &
               Coeff%eezy * ey(i,j,k) + Coeff%eezy * ey(i,j,1 + k) + 4. * Coeff%eezz * ez(i,j,k))/4. -                      &
               ((4. * Coeff%ehzz * hy(-1 + i,j,k) - 4. * Coeff%ehzz * hy(i,j,k) - Coeff%ehzy * hz(-1 + i,-1 + j,k) -              &
               Coeff%ehzy * hz(-1 + i,-1 + j,1 + k) - Coeff%ehzy * hz(-1 + i,j,k) -                                 &
               Coeff%ehzy * hz(-1 + i,j,1 + k) + Coeff%ehzy * hz(i,-1 + j,k) +                                      &
               Coeff%ehzy * hz(i,-1 + j,1 + k) + Coeff%ehzy * hz(i,j,k) + Coeff%ehzy * hz(i,j,1 + k)) * Idxh(i))/        &
               4. + ((4. * Coeff%ehzz * hx(i,-1 + j,k) - 4. * Coeff%ehzz * hx(i,j,k) - Coeff%ehzx * hz(-1 + i,-1 + j,k) -        &
               Coeff%ehzx * hz(-1 + i,-1 + j,1 + k) + Coeff%ehzx * hz(-1 + i,j,k) +                                 &
               Coeff%ehzx * hz(-1 + i,j,1 + k) - Coeff%ehzx * hz(i,-1 + j,k) -                                      &
               Coeff%ehzx * hz(i,-1 + j,1 + k) + Coeff%ehzx * hz(i,j,k) + Coeff%ehzx * hz(i,j,1 + k)) * Idyh(j))/        &
               4. - ((Coeff%ehzy * hx(i,-1 + j,-1 + k) - Coeff%ehzy * hx(i,-1 + j,1 + k) +                              &
               Coeff%ehzy * hx(i,j,-1 + k) - Coeff%ehzy * hx(i,j,1 + k) - Coeff%ehzx * hy(-1 + i,j,-1 + k) +           &
               Coeff%ehzx * hy(-1 + i,j,1 + k) - Coeff%ehzx * hy(i,j,-1 + k) + Coeff%ehzx * hy(i,j,1 + k)) *             &
               Idze(k))/4.0_RKIND
            endif
         end do
      end do

      !store it

      !
      do jmed=1,AniMed%NumMed
         dummy => AniMed%info(jmed)
         do i1=1,dummy%NumNodesEx
            i=dummy%Ex_i(i1)
            j=dummy%Ex_j(i1)
            k=dummy%Ex_k(i1)
            !
            ex(i,j,k)=dummy%Ex_value(i1)
         end do
         do i1=1,dummy%NumNodesEy
            i=dummy%Ey_i(i1)
            j=dummy%Ey_j(i1)
            k=dummy%Ey_k(i1)
            !!
            ey(i,j,k)=dummy%Ey_value(i1)
         end do
         do i1=1,dummy%NumNodesEz
            i=dummy%Ez_i(i1)
            j=dummy%Ez_j(i1)
            k=dummy%Ez_k(i1)
            !
            ez(i,j,k)=dummy%Ez_value(i1)
         end do
      end do

      RETURN

   end subroutine AdvanceAnisotropicE


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Subroutine to advance the E field in the Anisotropic (no need to advance the magnetic field)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine AdvanceAnisotropicH(sggalloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)

      type (XYZlimit_t), dimension (1:6), intent(in)                      ::  sggAlloc


      REAL (KIND=RKIND)   , intent(inout)      :: &
      Ex(sggAlloc(iEx)%XI : sggAlloc(iEx)%XE,sggAlloc(iEx)%YI : sggAlloc(iEx)%YE,sggAlloc(iEx)%ZI : sggAlloc(iEx)%ZE),&
      Ey(sggAlloc(iEy)%XI : sggAlloc(iEy)%XE,sggAlloc(iEy)%YI : sggAlloc(iEy)%YE,sggAlloc(iEy)%ZI : sggAlloc(iEy)%ZE),&
      Ez(sggAlloc(iEz)%XI : sggAlloc(iEz)%XE,sggAlloc(iEz)%YI : sggAlloc(iEz)%YE,sggAlloc(iEz)%ZI : sggAlloc(iEz)%ZE),&
      Hx(sggAlloc(iHx)%XI : sggAlloc(iHx)%XE,sggAlloc(iHx)%YI : sggAlloc(iHx)%YE,sggAlloc(iHx)%ZI : sggAlloc(iHx)%ZE),&
      Hy(sggAlloc(iHy)%XI : sggAlloc(iHy)%XE,sggAlloc(iHy)%YI : sggAlloc(iHy)%YE,sggAlloc(iHy)%ZI : sggAlloc(iHy)%ZE),&
      Hz(sggAlloc(iHz)%XI : sggAlloc(iHz)%XE,sggAlloc(iHz)%YI : sggAlloc(iHz)%YE,sggAlloc(iHz)%ZI : sggAlloc(iHz)%ZE)
      REAL (KIND=RKIND) , dimension (:)   , intent(in)      ::  &
      Idxe(sggALLOC(iHx)%XI : sggALLOC(iHx)%XE), &
      Idye(sggALLOC(iHy)%YI : sggALLOC(iHy)%YE), &
      Idze(sggALLOC(iHz)%ZI : sggALLOC(iHz)%ZE), &
      Idxh(sggALLOC(iEx)%XI : sggALLOC(iEx)%XE), &
      Idyh(sggALLOC(iEy)%YI : sggALLOC(iEy)%YE), &
      Idzh(sggALLOC(iEz)%ZI : sggALLOC(iEz)%ZE)
      !

      type (Anisotropicinfo_t), pointer :: dummy
      type (coeff_t), pointer :: coeff
      integer (kind=4)  ::  jmed,i1,i,j,k


      do jmed=1,AniMed%NumMed
         dummy => AniMed%info(jmed)
         !
         do i1=1,dummy%NumNodesHx
            if (Dummy%Hx_Shared(i1)%times == 1) then
               Coeff => AniMed%info(jmed)%coeff
            else
               Coeff => AniMed%info(jmed)%Hx_Shared(i1)%Coeff
            endif
            i=dummy%Hx_i(i1)
            j=dummy%Hx_j(i1)
            k=dummy%Hx_k(i1)
            !
            if  (dummy%IsOnlyThinSlot) then
               dummy%hx_value(i1) = hx(i,j,k) + Coeff%hexx*(ez(i,j,k) - ez(i,1 + j,k))*Idye(j) - &
               Coeff%hexx*(ey(i,j,k) - ey(i,j,1 + k))*Idze(k)
            else
               dummy%hx_value(i1) =                                                                                         &
               (4. * Coeff%hhxx * hx(i,j,k) + Coeff%hhxy * hy(-1 + i,j,k) + Coeff%hhxy * hy(-1 + i,1 + j,k) +               &
               Coeff%hhxy * hy(i,j,k) + Coeff%hhxy * hy(i,1 + j,k) + Coeff%hhxz * hz(-1 + i,j,k) +                       &
               Coeff%hhxz * hz(-1 + i,j,1 + k) + Coeff%hhxz * hz(i,j,k) + Coeff%hhxz * hz(i,j,1 + k))/4. +               &
               ((Coeff%hexz * ey(-1 + i,j,k) + Coeff%hexz * ey(-1 + i,j,1 + k) - Coeff%hexz * ey(1 + i,j,k) -               &
               Coeff%hexz * ey(1 + i,j,1 + k) - Coeff%hexy * ez(-1 + i,j,k) - Coeff%hexy * ez(-1 + i,1 + j,k) +        &
               Coeff%hexy * ez(1 + i,j,k) + Coeff%hexy * ez(1 + i,1 + j,k)) * Idxh(i))/4. -                            &
               ((Coeff%hexz * ex(-1 + i,j,k) + Coeff%hexz * ex(-1 + i,j,1 + k) - Coeff%hexz * ex(-1 + i,1 + j,k) -          &
               Coeff%hexz * ex(-1 + i,1 + j,1 + k) + Coeff%hexz * ex(i,j,k) + Coeff%hexz * ex(i,j,1 + k) -             &
               Coeff%hexz * ex(i,1 + j,k) - Coeff%hexz * ex(i,1 + j,1 + k) - 4. * Coeff%hexx * ez(i,j,k) +             &
               4. * Coeff%hexx * ez(i,1 + j,k)) * Idye(j))/4. +                                                        &
               ((Coeff%hexy * ex(-1 + i,j,k) - Coeff%hexy * ex(-1 + i,j,1 + k) + Coeff%hexy * ex(-1 + i,1 + j,k) -          &
               Coeff%hexy * ex(-1 + i,1 + j,1 + k) + Coeff%hexy * ex(i,j,k) - Coeff%hexy * ex(i,j,1 + k) +             &
               Coeff%hexy * ex(i,1 + j,k) - Coeff%hexy * ex(i,1 + j,1 + k) - 4. * Coeff%hexx * ey(i,j,k) +             &
               4. * Coeff%hexx * ey(i,j,1 + k)) * Idze(k))/4.0_RKIND
            endif
         end do
         do i1=1,dummy%NumNodesHy
            if (Dummy%Hy_Shared(i1)%times == 1) then
               Coeff => AniMed%info(jmed)%coeff
            else
               Coeff => AniMed%info(jmed)%Hy_Shared(i1)%Coeff
            endif
            i=dummy%Hy_i(i1)
            j=dummy%Hy_j(i1)
            k=dummy%Hy_k(i1)
            !
            if  (dummy%IsOnlyThinSlot) then
               dummy%hy_value(i1) =  hy(i,j,k) - Coeff%heyy*(ez(i,j,k) - ez(1 + i,j,k))*Idxe(i) + &
               Coeff%heyy*(ex(i,j,k) - ex(i,j,1 + k))*Idze(k)
            else
               dummy%hy_value(i1) =                                                                                            &
               (Coeff%hhyx * hx(i,-1 + j,k) + Coeff%hhyx * hx(i,j,k) + Coeff%hhyx * hx(1 + i,-1 + j,k) +                       &
               Coeff%hhyx * hx(1 + i,j,k) + 4. * Coeff%hhyy * hy(i,j,k) + Coeff%hhyz * hz(i,-1 + j,k) +                     &
               Coeff%hhyz * hz(i,-1 + j,1 + k) + Coeff%hhyz * hz(i,j,k) + Coeff%hhyz * hz(i,j,1 + k))/4. +                  &
               ((Coeff%heyz * ey(i,-1 + j,k) + Coeff%heyz * ey(i,-1 + j,1 + k) + Coeff%heyz * ey(i,j,k) +                      &
               Coeff%heyz * ey(i,j,1 + k) - Coeff%heyz * ey(1 + i,-1 + j,k) -                                             &
               Coeff%heyz * ey(1 + i,-1 + j,1 + k) - Coeff%heyz * ey(1 + i,j,k) -                                         &
               Coeff%heyz * ey(1 + i,j,1 + k) - 4. * Coeff%heyy * ez(i,j,k) + 4. * Coeff%heyy * ez(1 + i,j,k)) * Idxe(i)  &
               )/4. - ((Coeff%heyz * ex(i,-1 + j,k) + Coeff%heyz * ex(i,-1 + j,1 + k) -                                     &
               Coeff%heyz * ex(i,1 + j,k) - Coeff%heyz * ex(i,1 + j,1 + k) - Coeff%heyx * ez(i,-1 + j,k) +                &
               Coeff%heyx * ez(i,1 + j,k) - Coeff%heyx * ez(1 + i,-1 + j,k) + Coeff%heyx * ez(1 + i,1 + j,k)) *           &
               Idyh(j))/4. + ((4. * Coeff%heyy * ex(i,j,k) - 4. * Coeff%heyy * ex(i,j,1 + k) -                              &
               Coeff%heyx * ey(i,-1 + j,k) + Coeff%heyx * ey(i,-1 + j,1 + k) - Coeff%heyx * ey(i,j,k) +                   &
               Coeff%heyx * ey(i,j,1 + k) - Coeff%heyx * ey(1 + i,-1 + j,k) +                                             &
               Coeff%heyx * ey(1 + i,-1 + j,1 + k) - Coeff%heyx * ey(1 + i,j,k) + Coeff%heyx * ey(1 + i,j,1 + k)          &
               ) * Idze(k))/4.0_RKIND
            endif
         end do
         do i1=1,dummy%NumNodesHz
            if (Dummy%Hz_Shared(i1)%times == 1) then
               Coeff => AniMed%info(jmed)%coeff
            else
               Coeff => AniMed%info(jmed)%Hz_Shared(i1)%Coeff
            endif
            i=dummy%Hz_i(i1)
            j=dummy%Hz_j(i1)
            k=dummy%Hz_k(i1)
            !
            if  (dummy%IsOnlyThinSlot) then
               dummy%hz_value(i1) = hz(i,j,k) + Coeff%hezz*(ey(i,j,k) - ey(1 + i,j,k))*Idxe(i) - &
               Coeff%hezz*(ex(i,j,k) - ex(i,1 + j,k))*Idye(j)
            else
               dummy%hz_value(i1) =                                                                                        &
               (Coeff%hhzx * hx(i,j,-1 + k) + Coeff%hhzx * hx(i,j,k) + Coeff%hhzx * hx(1 + i,j,-1 + k) +                    &
               Coeff%hhzx * hx(1 + i,j,k) + Coeff%hhzy * hy(i,j,-1 + k) + Coeff%hhzy * hy(i,j,k) +                      &
               Coeff%hhzy * hy(i,1 + j,-1 + k) + Coeff%hhzy * hy(i,1 + j,k) + 4. * Coeff%hhzz * hz(i,j,k))/4. +            &
               ((4. * Coeff%hezz * ey(i,j,k) - 4. * Coeff%hezz * ey(1 + i,j,k) - Coeff%hezy * ez(i,j,-1 + k) -                   &
               Coeff%hezy * ez(i,j,k) - Coeff%hezy * ez(i,1 + j,-1 + k) - Coeff%hezy * ez(i,1 + j,k) +                &
               Coeff%hezy * ez(1 + i,j,-1 + k) + Coeff%hezy * ez(1 + i,j,k) +                                      &
               Coeff%hezy * ez(1 + i,1 + j,-1 + k) + Coeff%hezy * ez(1 + i,1 + j,k)) * Idxe(i))/4. -                 &
               ((4. * Coeff%hezz * ex(i,j,k) - 4. * Coeff%hezz * ex(i,1 + j,k) - Coeff%hezx * ez(i,j,-1 + k) -                   &
               Coeff%hezx * ez(i,j,k) + Coeff%hezx * ez(i,1 + j,-1 + k) + Coeff%hezx * ez(i,1 + j,k) -                &
               Coeff%hezx * ez(1 + i,j,-1 + k) - Coeff%hezx * ez(1 + i,j,k) +                                      &
               Coeff%hezx * ez(1 + i,1 + j,-1 + k) + Coeff%hezx * ez(1 + i,1 + j,k)) * Idye(j))/4. +                 &
               ((Coeff%hezy * ex(i,j,-1 + k) - Coeff%hezy * ex(i,j,1 + k) + Coeff%hezy * ex(i,1 + j,-1 + k) -              &
               Coeff%hezy * ex(i,1 + j,1 + k) - Coeff%hezx * ey(i,j,-1 + k) + Coeff%hezx * ey(i,j,1 + k) -            &
               Coeff%hezx * ey(1 + i,j,-1 + k) + Coeff%hezx * ey(1 + i,j,1 + k)) * Idzh(k))/4.0_RKIND
            endif
         end do
      end do


      !store it

      !
      do jmed=1,AniMed%NumMed
         dummy => AniMed%info(jmed)
         do i1=1,dummy%NumNodesHx
            i=dummy%Hx_i(i1)
            j=dummy%Hx_j(i1)
            k=dummy%Hx_k(i1)
            !
            Hx(i,j,k)=dummy%Hx_value(i1)
         end do
         do i1=1,dummy%NumNodesHy
            i=dummy%Hy_i(i1)
            j=dummy%Hy_j(i1)
            k=dummy%Hy_k(i1)
            !!
            Hy(i,j,k)=dummy%Hy_value(i1)
         end do
         do i1=1,dummy%NumNodesHz
            i=dummy%Hz_i(i1)
            j=dummy%Hz_j(i1)
            k=dummy%Hz_k(i1)
            !
            Hz(i,j,k)=dummy%Hz_value(i1)
         end do
      end do

      RETURN

   end subroutine AdvanceAnisotropicH


   subroutine DestroyAnisotropic(sgg)
      type (SGGFDTDINFO), intent(INOUT)         ::  sgg
      !
      integer (kind=4)  ::  jmed,i
      !free up memory
      do jmed=1,AniMed%NumMed
         deallocate (AniMed%info(jmed)%Ex_i,AniMed%info(jmed)%Ex_j,AniMed%info(jmed)%Ex_k)
         deallocate (AniMed%info(jmed)%Ey_i,AniMed%info(jmed)%Ey_j,AniMed%info(jmed)%Ey_k)
         deallocate (AniMed%info(jmed)%Ez_i,AniMed%info(jmed)%Ez_j,AniMed%info(jmed)%Ez_k)
         deallocate (AniMed%info(jmed)%Hx_i,AniMed%info(jmed)%Hx_j,AniMed%info(jmed)%Hx_k)
         deallocate (AniMed%info(jmed)%Hy_i,AniMed%info(jmed)%Hy_j,AniMed%info(jmed)%Hy_k)
         deallocate (AniMed%info(jmed)%Hz_i,AniMed%info(jmed)%Hz_j,AniMed%info(jmed)%Hz_k)
         deallocate (AniMed%info(jmed)%Ex_value,AniMed%info(jmed)%Ex_Shared)
         deallocate (AniMed%info(jmed)%Ey_value,AniMed%info(jmed)%Ey_Shared)
         deallocate (AniMed%info(jmed)%Ez_value,AniMed%info(jmed)%Ez_Shared)
         deallocate (AniMed%info(jmed)%Hx_value,AniMed%info(jmed)%Hx_Shared)
         deallocate (AniMed%info(jmed)%Hy_value,AniMed%info(jmed)%Hy_Shared)
         deallocate (AniMed%info(jmed)%Hz_value,AniMed%info(jmed)%Hz_Shared)
      end do
      do i=1,sgg%NumMedia
         if (sgg%Med(i)%Is%Anisotropic)  then
            deallocate (sgg%Med(i)%Anisotropic)
         endif
      end do
      deallocate (AniMed%info)

   end subroutine


   !funcion para publicar el valor de Med

   function GetMed() result(r)
      type(AnisotropicMed_t), pointer  ::  r

      r=>AniMed

      return
   end function

   !found by the mathematica notebook
   Subroutine calc_anisotropicconstants(sgg,eps00,mu00)
        type (SGGFDTDINFO), intent(IN)   ::  sgg
        REAL (KIND=RKIND) , intent(inout) :: Eps00, Mu00
        REAL (KIND=RKIND),  DIMENSION(3,3) ::  sigma,epr,mur,sigmaM
        integer (kind=4) :: jmed
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      cluz=1.0_RKIND/sqrt(eps0*mu0)

      !calculate the coefficients
        do jmed=1,AniMed%NumMed
            sigma  = AniMed%info(jmed)%sigma  
            sigmam = AniMed%info(jmed)%sigmam 
            mur    = AniMed%info(jmed)%mur    
            epr    = AniMed%info(jmed)%epr 
            !copied from hirf_anis1.nb
            call CalculateCoeff(epr,mur,sigma,sigmam,sgg%dt,AniMed%info(jmed)%coeff)
        end do
        return
   end subroutine calc_anisotropicconstants




   !found by the mathematica notebook
   Subroutine CalculateCoeff(epr,mur,sigma,sigmam,dt,coeff)

      type (coeff_t), intent(out) :: coeff
      REAL (KIND=RKIND),  DIMENSION(3,3), intent(in)  ::  sigma,epr,mur,sigmaM
      REAL (KIND=RKIND_tiempo) :: dt

      coeff%eexx = ((-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2          &
      * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt      &
      - sigma(3,1)/2.0_RKIND) +((eps0 * epr(2,1))/dt - sigma(2,1)/2.0_RKIND) *((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2      &
      * eps0 * epr(3,2) + dt * sigma(3,2)) -(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(3,3) +      &
      dt * sigma(3,3))) +((eps0 * epr(1,1))/dt - sigma(1,1)/2.0_RKIND) *(-((2 * eps0 * epr(2,3) + dt * sigma(2,3)) *(2 *      &
      eps0 * epr(3,2) + dt * sigma(3,2))) +(2 * eps0 * epr(2,2) + dt * sigma(2,2)) *(2 * eps0 * epr(3,3) + dt      &
      * sigma(3,3))))/((-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +      &
      (2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt +      &
      sigma(3,1)/2.0_RKIND) -(-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt      &
      + sigma(3,2)/2.0_RKIND) +(-((2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt      &
      + sigma(3,3)/2.0_RKIND))

      coeff%eexy = (4 * dt * eps0 *(2 * eps0 *((epr(1,3) * epr(3,2) - epr(1,2) * epr(3,3)) * sigma(2,2) + epr      &
      (2,3) *(-(epr(3,2) * sigma(1,2)) + epr(1,2) * sigma(3,2)) + epr(2,2) *(epr(3,3) * sigma(1,2) - epr(1,3)      &
      * sigma(3,2))) + dt *(epr(3,2) *(sigma(1,3) * sigma(2,2) - sigma(1,2) * sigma(2,3)) + epr(2,2)      &
      *(-(sigma(1,3) * sigma(3,2)) + sigma(1,2) * sigma(3,3)) + epr(1,2) *(sigma(2,3) * sigma(3,2) - sigma(2,2)      &
      * sigma(3,3)))))/(8 * eps0**3 *(epr(1,3) *(epr(2,2) * epr(3,1) - epr(2,1) * epr(3,2)) + epr(1,2) *      &
      (-(epr(2,3) * epr(3,1)) + epr(2,1) * epr(3,3)) + epr(1,1) *(epr(2,3) * epr(3,2) - epr(2,2) * epr(3,3))) + 2      &
      * dt**2.0_RKIND * eps0 *(epr(3,1) * sigma(1,3) * sigma(2,2) + epr(3,3) *(sigma(1,2) * sigma(2,1) - sigma      &
      (1,1) * sigma(2,2)) - epr(3,1) * sigma(1,2) * sigma(2,3) + epr(3,2) *(-(sigma(1,3)      &
      * sigma(2,1)) + sigma(1,1) * sigma(2,3)) - epr(2,3) * sigma(1,2) * sigma(3,1) + epr(2,2) * sigma(1,3)      &
      * sigma(3,1) + epr(1,3) * sigma(2,2) * sigma(3,1) - epr(1,2) * sigma(2,3) * sigma(3,1) + epr(2,3) *      &
      sigma(1,1) * sigma(3,2) - epr(2,1) * sigma(1,3) * sigma(3,2) - epr(1,3) * sigma(2,1) * sigma(3,2) + epr(1,1)      &
      * sigma(2,3) * sigma(3,2) +(-(epr(2,2) * sigma(1,1)) + epr(2,1) * sigma(1,2) + epr(1,2) * sigma(2,1)      &
      - epr(1,1) * sigma(2,2)) * sigma(3,3)) + 4 * dt * eps0**2.0_RKIND *(epr(2,1)      &
      * epr(3,3) * sigma(1,2) - epr(2,1) * epr(3,2) * sigma(1,3) - epr(1,3) * epr(3,2) * sigma(2,1) + epr(1,2) * epr      &
      (3,3) * sigma(2,1) + epr(1,3) * epr(3,1) * sigma(2,2) - epr(1,1) * epr(3,3) * sigma(2,2) - epr(1,2) * epr      &
      (3,1) * sigma(2,3) + epr(1,1) * epr(3,2) * sigma(2,3) - epr(1,3) * epr(2,1) * sigma(3,2) + epr(2,3) *(epr      &
      (3,2) * sigma(1,1) - epr(3,1) * sigma(1,2) - epr(1,2) * sigma(3,1) + epr(1,1) * sigma(3,2)) + epr(1,2)      &
      * epr(2,1) * sigma(3,3) + epr(2,2) *(-(epr(3,3) * sigma(1,1)) + epr(3,1) * sigma(1,3) + epr(1,3) *      &
      sigma(3,1) - epr(1,1) * sigma(3,3))) + dt**3 *(sigma(1,3) *(sigma(2,2) * sigma(3,1) - sigma(2,1) * sigma      &
      (3,2)) + sigma(1,2) *(-(sigma(2,3) * sigma(3,1)) + sigma(2,1) * sigma(3,3)) + sigma(1,1) *(sigma(2,3) * sigma      &
      (3,2) - sigma(2,2) * sigma(3,3))))

      coeff%eexz = (4 * dt * eps0 *(2 * eps0 *((epr(1,3) * epr(3,2) - epr(1,2) * epr(3,3)) * sigma(2,3) + epr      &
      (2,3) *(-(epr(3,2) * sigma(1,3)) + epr(1,2) * sigma(3,3)) + epr(2,2) *(epr(3,3) * sigma(1,3) - epr(1,3)      &
      * sigma(3,3))) + dt *(epr(3,3) *(sigma(1,3) * sigma(2,2) - sigma(1,2) * sigma(2,3)) + epr(2,3)      &
      *(-(sigma(1,3) * sigma(3,2)) + sigma(1,2) * sigma(3,3)) + epr(1,3) *(sigma(2,3) * sigma(3,2) - sigma(2,2)      &
      * sigma(3,3)))))/(8 * eps0**3 *(epr(1,3) *(epr(2,2) * epr(3,1) - epr(2,1) * epr(3,2)) + epr(1,2) *      &
      (-(epr(2,3) * epr(3,1)) + epr(2,1) * epr(3,3)) + epr(1,1) *(epr(2,3) * epr(3,2) - epr(2,2) * epr(3,3))) + 2      &
      * dt**2.0_RKIND * eps0 *(epr(3,1) * sigma(1,3) * sigma(2,2) + epr(3,3) *(sigma(1,2) * sigma(2,1) - sigma      &
      (1,1) * sigma(2,2)) - epr(3,1) * sigma(1,2) * sigma(2,3) + epr(3,2) *(-(sigma(1,3)      &
      * sigma(2,1)) + sigma(1,1) * sigma(2,3)) - epr(2,3) * sigma(1,2) * sigma(3,1) + epr(2,2) * sigma(1,3)      &
      * sigma(3,1) + epr(1,3) * sigma(2,2) * sigma(3,1) - epr(1,2) * sigma(2,3) * sigma(3,1) + epr(2,3) *      &
      sigma(1,1) * sigma(3,2) - epr(2,1) * sigma(1,3) * sigma(3,2) - epr(1,3) * sigma(2,1) * sigma(3,2) + epr(1,1)      &
      * sigma(2,3) * sigma(3,2) +(-(epr(2,2) * sigma(1,1)) + epr(2,1) * sigma(1,2) + epr(1,2) * sigma(2,1)      &
      - epr(1,1) * sigma(2,2)) * sigma(3,3)) + 4 * dt * eps0**2.0_RKIND *(epr(2,1)      &
      * epr(3,3) * sigma(1,2) - epr(2,1) * epr(3,2) * sigma(1,3) - epr(1,3) * epr(3,2) * sigma(2,1) + epr(1,2) * epr      &
      (3,3) * sigma(2,1) + epr(1,3) * epr(3,1) * sigma(2,2) - epr(1,1) * epr(3,3) * sigma(2,2) - epr(1,2) * epr      &
      (3,1) * sigma(2,3) + epr(1,1) * epr(3,2) * sigma(2,3) - epr(1,3) * epr(2,1) * sigma(3,2) + epr(2,3) *(epr      &
      (3,2) * sigma(1,1) - epr(3,1) * sigma(1,2) - epr(1,2) * sigma(3,1) + epr(1,1) * sigma(3,2)) + epr(1,2)      &
      * epr(2,1) * sigma(3,3) + epr(2,2) *(-(epr(3,3) * sigma(1,1)) + epr(3,1) * sigma(1,3) + epr(1,3) *      &
      sigma(3,1) - epr(1,1) * sigma(3,3))) + dt**3 *(sigma(1,3) *(sigma(2,2) * sigma(3,1) - sigma(2,1) * sigma      &
      (3,2)) + sigma(1,2) *(-(sigma(2,3) * sigma(3,1)) + sigma(2,1) * sigma(3,3)) + sigma(1,1) *(sigma(2,3) * sigma      &
      (3,2) - sigma(2,2) * sigma(3,3))))

      coeff%eeyx = (4 * dt * eps0 *(2 * eps0 *((epr(1,3) * epr(3,1) - epr(1,1) * epr(3,3)) * sigma(2,1) + epr      &
      (2,3) *(-(epr(3,1) * sigma(1,1)) + epr(1,1) * sigma(3,1)) + epr(2,1) *(epr(3,3) * sigma(1,1) - epr(1,3)      &
      * sigma(3,1))) + dt *(epr(3,1) *(sigma(1,3) * sigma(2,1) - sigma(1,1) * sigma(2,3)) + epr(2,1)      &
      *(-(sigma(1,3) * sigma(3,1)) + sigma(1,1) * sigma(3,3)) + epr(1,1) *(sigma(2,3) * sigma(3,1) - sigma(2,1)      &
      * sigma(3,3)))))/(8 * eps0**3 *(epr(1,3) *(-(epr(2,2) * epr(3,1)) + epr(2,1) * epr(3,2)) + epr(1,2)      &
      *(epr(2,3) * epr(3,1) - epr(2,1) * epr(3,3)) + epr(1,1) *(-(epr(2,3) * epr(3,2)) + epr(2,2) * epr(3,3))) + 2      &
      * dt**2.0_RKIND * eps0 *(-(epr(3,1) * sigma(1,3) * sigma(2,2)) + epr(3,3) *(-(sigma(1,2) * sigma(2,1)) +      &
      sigma(1,1) * sigma(2,2)) + epr(3,1) * sigma(1,2) * sigma(2,3) + epr(3,2) *(sigma(1,3)      &
      * sigma(2,1) - sigma(1,1) * sigma(2,3)) + epr(2,3) * sigma(1,2) * sigma(3,1) - epr(2,2) * sigma(1,3)      &
      * sigma(3,1) - epr(1,3) * sigma(2,2) * sigma(3,1) + epr(1,2) * sigma(2,3) * sigma(3,1) - epr(2,3) *      &
      sigma(1,1) * sigma(3,2) + epr(2,1) * sigma(1,3) * sigma(3,2) + epr(1,3) * sigma(2,1) * sigma(3,2) - epr(1,1)      &
      * sigma(2,3) * sigma(3,2) +(epr(2,2) * sigma(1,1) - epr(2,1) * sigma(1,2) - epr(1,2) * sigma(2,1) +      &
      epr(1,1) * sigma(2,2)) * sigma(3,3)) + 4 * dt * eps0**2.0_RKIND *(-(epr(2,1) * epr(3,3) * sigma(1,2)) + epr(2,1)      &
      * epr(3,2) * sigma(1,3) + epr(1,3) * epr(3,2) * sigma(2,1) - epr(1,2)      &
      * epr(3,3) * sigma(2,1) - epr(1,3) * epr(3,1) * sigma(2,2) + epr(1,1) * epr(3,3) * sigma(2,2) + epr      &
      (1,2) * epr(3,1) * sigma(2,3) - epr(1,1) * epr(3,2) * sigma(2,3) + epr(1,3) * epr(2,1) * sigma(3,2) + epr      &
      (2,3) *(-(epr(3,2) * sigma(1,1)) + epr(3,1) * sigma(1,2) + epr(1,2) * sigma(3,1) - epr(1,1) * sigma(3,2)) -      &
      epr(1,2) * epr(2,1) * sigma(3,3) + epr(2,2) *(epr(3,3) * sigma(1,1) - epr(3,1) * sigma(1,3) - epr(1,3)      &
      * sigma(3,1) + epr(1,1) * sigma(3,3))) + dt**3 *(sigma(1,3) *(-(sigma(2,2) * sigma(3,1)) + sigma      &
      (2,1) * sigma(3,2)) + sigma(1,2) *(sigma(2,3) * sigma(3,1) - sigma(2,1) * sigma(3,3)) + sigma(1,1) *(-(sigma      &
      (2,3) * sigma(3,2)) + sigma(2,2) * sigma(3,3))))

      coeff%eeyy = (((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1)) -(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt      &
      - sigma(3,2)/2.0_RKIND) +((eps0 * epr(2,2))/dt - sigma(2,2)/2.0_RKIND) *(-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2      &
      * eps0 * epr(3,1) + dt * sigma(3,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(3,3) +      &
      dt * sigma(3,3))) +((eps0 * epr(1,2))/dt - sigma(1,2)/2.0_RKIND) *((2 * eps0 * epr(2,3) + dt * sigma(2,3)) *(2      &
      * eps0 * epr(3,1) + dt * sigma(3,1)) -(2 * eps0 * epr(2,1) + dt * sigma(2,1)) *(2 * eps0 * epr(3,3) +      &
      dt * sigma(3,3))))/((-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2      &
      * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt      &
      + sigma(3,1)/2.0_RKIND) -(-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0      &
      * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,3) + dt * sigma      &
      (2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0      &
      * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,2) + dt * sigma      &
      (2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%eeyz = (4 * dt * eps0 *(2 * eps0 *((-(epr(1,3) * epr(3,1)) + epr(1,1) * epr(3,3)) * sigma(2,3) +      &
      epr(2,3) *(epr(3,1) * sigma(1,3) - epr(1,1) * sigma(3,3)) + epr(2,1) *(-(epr(3,3) * sigma(1,3)) + epr(1,3)      &
      * sigma(3,3))) + dt *(epr(3,3) *(-(sigma(1,3) * sigma(2,1)) + sigma(1,1) * sigma(2,3)) + epr(2,3) *      &
      (sigma(1,3) * sigma(3,1) - sigma(1,1) * sigma(3,3)) + epr(1,3) *(-(sigma(2,3) * sigma(3,1)) + sigma(2,1)      &
      * sigma(3,3)))))/(8 * eps0**3 *(epr(1,3) *(epr(2,2) * epr(3,1) - epr(2,1) * epr(3,2)) + epr(1,2) *      &
      (-(epr(2,3) * epr(3,1)) + epr(2,1) * epr(3,3)) + epr(1,1) *(epr(2,3)      &
      * epr(3,2) - epr(2,2) * epr(3,3))) + 2 * dt**2.0_RKIND * eps0 *(epr(3,1) * sigma(1,3) * sigma(2,2) + epr(3,3) *(sigma      &
      (1,2) * sigma(2,1) - sigma(1,1) * sigma(2,2)) - epr(3,1) * sigma(1,2) * sigma(2,3) + epr(3,2) *(-(sigma(1,3)      &
      * sigma(2,1)) + sigma(1,1) * sigma(2,3)) - epr(2,3) * sigma(1,2) * sigma(3,1) + epr(2,2) * sigma(1,3)      &
      * sigma(3,1) + epr(1,3) * sigma(2,2) * sigma(3,1) - epr(1,2) * sigma(2,3) * sigma(3,1) + epr(2,3) *      &
      sigma(1,1) * sigma(3,2) - epr(2,1) * sigma(1,3) * sigma(3,2) - epr(1,3) * sigma(2,1) * sigma(3,2) + epr(1,1)      &
      * sigma(2,3) * sigma(3,2) +(-(epr(2,2) * sigma(1,1)) + epr(2,1) * sigma(1,2) + epr(1,2) * sigma(2,1)      &
      - epr(1,1) * sigma(2,2)) * sigma(3,3)) + 4 * dt * eps0**2.0_RKIND *(epr(2,1)      &
      * epr(3,3) * sigma(1,2) - epr(2,1) * epr(3,2) * sigma(1,3) - epr(1,3) * epr(3,2) * sigma(2,1) + epr(1,2) * epr      &
      (3,3) * sigma(2,1) + epr(1,3) * epr(3,1) * sigma(2,2) - epr(1,1) * epr(3,3) * sigma(2,2) - epr(1,2) * epr      &
      (3,1) * sigma(2,3) + epr(1,1) * epr(3,2) * sigma(2,3) - epr(1,3) * epr(2,1) * sigma(3,2) + epr(2,3) *(epr      &
      (3,2) * sigma(1,1) - epr(3,1) * sigma(1,2) - epr(1,2) * sigma(3,1) + epr(1,1) * sigma(3,2)) + epr(1,2)      &
      * epr(2,1) * sigma(3,3) + epr(2,2) *(-(epr(3,3) * sigma(1,1)) + epr(3,1) * sigma(1,3) + epr(1,3) *      &
      sigma(3,1) - epr(1,1) * sigma(3,3))) + dt**3 *(sigma(1,3) *(sigma(2,2) * sigma(3,1) - sigma(2,1) * sigma      &
      (3,2)) + sigma(1,2) *(-(sigma(2,3) * sigma(3,1)) + sigma(2,1) * sigma(3,3)) + sigma(1,1) *(sigma(2,3) * sigma      &
      (3,2) - sigma(2,2) * sigma(3,3))))

      coeff%eezx = (4 * dt * eps0 *(2 * eps0 *((epr(1,2) * epr(3,1) - epr(1,1) * epr(3,2)) * sigma(2,1) + epr      &
      (2,2) *(-(epr(3,1) * sigma(1,1)) + epr(1,1) * sigma(3,1)) + epr(2,1) *(epr(3,2) * sigma(1,1) - epr(1,2)      &
      * sigma(3,1))) + dt *(epr(3,1) *(sigma(1,2) * sigma(2,1) - sigma(1,1) * sigma(2,2)) + epr(2,1)      &
      *(-(sigma(1,2) * sigma(3,1)) + sigma(1,1) * sigma(3,2)) + epr(1,1) *(sigma(2,2) * sigma(3,1) - sigma(2,1)      &
      * sigma(3,2)))))/(8 * eps0**3 *(epr(1,3) *(epr(2,2) * epr(3,1) - epr(2,1) * epr(3,2)) + epr(1,2) *      &
      (-(epr(2,3) * epr(3,1)) + epr(2,1) * epr(3,3)) + epr(1,1) *(epr(2,3) * epr(3,2) - epr(2,2) * epr(3,3))) + 2      &
      * dt**2.0_RKIND * eps0 *(epr(3,1) * sigma(1,3) * sigma(2,2) + epr(3,3) *(sigma(1,2) * sigma(2,1) - sigma      &
      (1,1) * sigma(2,2)) - epr(3,1) * sigma(1,2) * sigma(2,3) + epr(3,2) *(-(sigma(1,3)      &
      * sigma(2,1)) + sigma(1,1) * sigma(2,3)) - epr(2,3) * sigma(1,2) * sigma(3,1) + epr(2,2) * sigma(1,3)      &
      * sigma(3,1) + epr(1,3) * sigma(2,2) * sigma(3,1) - epr(1,2) * sigma(2,3) * sigma(3,1) + epr(2,3) *      &
      sigma(1,1) * sigma(3,2) - epr(2,1) * sigma(1,3) * sigma(3,2) - epr(1,3) * sigma(2,1) * sigma(3,2) + epr(1,1)      &
      * sigma(2,3) * sigma(3,2) +(-(epr(2,2) * sigma(1,1)) + epr(2,1) * sigma(1,2) + epr(1,2) * sigma(2,1)      &
      - epr(1,1) * sigma(2,2)) * sigma(3,3)) + 4 * dt * eps0**2.0_RKIND *(epr(2,1)      &
      * epr(3,3) * sigma(1,2) - epr(2,1) * epr(3,2) * sigma(1,3) - epr(1,3) * epr(3,2) * sigma(2,1) + epr(1,2) * epr      &
      (3,3) * sigma(2,1) + epr(1,3) * epr(3,1) * sigma(2,2) - epr(1,1) * epr(3,3) * sigma(2,2) - epr(1,2) * epr      &
      (3,1) * sigma(2,3) + epr(1,1) * epr(3,2) * sigma(2,3) - epr(1,3) * epr(2,1) * sigma(3,2) + epr(2,3) *(epr      &
      (3,2) * sigma(1,1) - epr(3,1) * sigma(1,2) - epr(1,2) * sigma(3,1) + epr(1,1) * sigma(3,2)) + epr(1,2)      &
      * epr(2,1) * sigma(3,3) + epr(2,2) *(-(epr(3,3) * sigma(1,1)) + epr(3,1) * sigma(1,3) + epr(1,3) *      &
      sigma(3,1) - epr(1,1) * sigma(3,3))) + dt**3 *(sigma(1,3) *(sigma(2,2) * sigma(3,1) - sigma(2,1) * sigma      &
      (3,2)) + sigma(1,2) *(-(sigma(2,3) * sigma(3,1)) + sigma(2,1) * sigma(3,3)) + sigma(1,1) *(sigma(2,3) * sigma      &
      (3,2) - sigma(2,2) * sigma(3,3))))

      coeff%eezy = (4 * dt * eps0 *(2 * eps0 *((-(epr(1,2) * epr(3,1)) + epr(1,1) * epr(3,2)) * sigma(2,2) +      &
      epr(2,2) *(epr(3,1) * sigma(1,2) - epr(1,1) * sigma(3,2)) + epr(2,1) *(-(epr(3,2) * sigma(1,2)) + epr(1,2)      &
      * sigma(3,2))) + dt *(epr(3,2) *(-(sigma(1,2) * sigma(2,1)) + sigma(1,1) * sigma(2,2)) + epr(2,2) *      &
      (sigma(1,2) * sigma(3,1) - sigma(1,1) * sigma(3,2)) + epr(1,2) *(-(sigma(2,2) * sigma(3,1)) + sigma(2,1)      &
      * sigma(3,2)))))/(8 * eps0**3 *(epr(1,3) *(-(epr(2,2) * epr(3,1)) + epr(2,1) * epr(3,2)) + epr(1,2)      &
      *(epr(2,3) * epr(3,1) - epr(2,1) * epr(3,3)) + epr(1,1) *(-(epr(2,3) * epr(3,2)) + epr(2,2) * epr(3,3))) + 2      &
      * dt**2.0_RKIND * eps0 *(-(epr(3,1) * sigma(1,3) * sigma(2,2)) + epr(3,3) *(-(sigma(1,2) * sigma(2,1)) +      &
      sigma(1,1) * sigma(2,2)) + epr(3,1) * sigma(1,2) * sigma(2,3) + epr(3,2) *(sigma(1,3)      &
      * sigma(2,1) - sigma(1,1) * sigma(2,3)) + epr(2,3) * sigma(1,2) * sigma(3,1) - epr(2,2) * sigma(1,3)      &
      * sigma(3,1) - epr(1,3) * sigma(2,2) * sigma(3,1) + epr(1,2) * sigma(2,3) * sigma(3,1) - epr(2,3) *      &
      sigma(1,1) * sigma(3,2) + epr(2,1) * sigma(1,3) * sigma(3,2) + epr(1,3) * sigma(2,1) * sigma(3,2) - epr(1,1)      &
      * sigma(2,3) * sigma(3,2) +(epr(2,2) * sigma(1,1) - epr(2,1) * sigma(1,2) - epr(1,2) * sigma(2,1) +      &
      epr(1,1) * sigma(2,2)) * sigma(3,3)) + 4 * dt * eps0**2.0_RKIND *(-(epr(2,1) * epr(3,3) * sigma(1,2)) + epr(2,1)      &
      * epr(3,2) * sigma(1,3) + epr(1,3) * epr(3,2) * sigma(2,1) - epr(1,2)      &
      * epr(3,3) * sigma(2,1) - epr(1,3) * epr(3,1) * sigma(2,2) + epr(1,1) * epr(3,3) * sigma(2,2) + epr      &
      (1,2) * epr(3,1) * sigma(2,3) - epr(1,1) * epr(3,2) * sigma(2,3) + epr(1,3) * epr(2,1) * sigma(3,2) + epr      &
      (2,3) *(-(epr(3,2) * sigma(1,1)) + epr(3,1) * sigma(1,2) + epr(1,2) * sigma(3,1) - epr(1,1) * sigma(3,2)) -      &
      epr(1,2) * epr(2,1) * sigma(3,3) + epr(2,2) *(epr(3,3) * sigma(1,1) - epr(3,1)      &
      * sigma(1,3) - epr(1,3) * sigma(3,1) + epr(1,1) * sigma(3,3))) + dt**3 *(sigma(1,3) *(-(sigma(2,2) *      &
      sigma(3,1)) + sigma(2,1) * sigma(3,2)) + sigma(1,2) *(sigma(2,3) * sigma(3,1) - sigma(2,1) * sigma(3,3)) +      &
      sigma(1,1) *(-(sigma(2,3) * sigma(3,2)) + sigma(2,2) * sigma(3,3))))

      coeff%eezz = (((eps0 * epr(2,3))/dt - sigma(2,3)/2.0_RKIND) *((2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2      &
      * eps0 * epr(3,1) + dt * sigma(3,1)) -(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(3,2) + dt * sigma(3,2)))      &
      +((eps0 * epr(1,3))/dt - sigma(1,3)/2.0_RKIND) *(-((2 * eps0 * epr(2,2) + dt * sigma(2,2)) *(2 * eps0 * epr(3,1) +      &
      dt * sigma(3,1))) +(2 * eps0 * epr(2,1) + dt * sigma(2,1)) *(2 * eps0 * epr(3,2) + dt * sigma(3,2))) +(-((2      &
      * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) +      &
      dt * sigma(1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt - sigma(3,3)/2.0_RKIND))/((-((2      &
      * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) +      &
      dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2      &
      * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) +      &
      dt * sigma(1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2      &
      * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) +      &
      dt * sigma(1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehxx = (-((2 * eps0 * epr(2,3) + dt * sigma(2,3)) *(2 * eps0 * epr(3,2) + dt * sigma(3,2))) +(2      &
      * eps0 * epr(2,2) + dt * sigma(2,2)) *(2 * eps0 * epr(3,3) + dt * sigma(3,3)))/((-((2 * eps0 * epr(1,3)      &
      + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 *      &
      eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0      &
      * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0      &
      * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehxy = ((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(3,2) + dt * sigma(3,2)) -(2      &
      * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(3,3) + dt * sigma(3,3)))/((-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0 * epr(1,2) + dt * sigma      &
      (1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehxz = (-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2      &
      * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3)))/((-((2 * eps0 * epr(1,3)      &
      + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 *      &
      eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0      &
      * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0      &
      * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehyx = ((2 * eps0 * epr(2,3) + dt * sigma(2,3)) *(2 * eps0 * epr(3,1) + dt * sigma(3,1)) -(2      &
      * eps0 * epr(2,1) + dt * sigma(2,1)) *(2 * eps0 * epr(3,3) + dt * sigma(3,3)))/((-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0 * epr(1,2) + dt * sigma      &
      (1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehyy = (-((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(3,1) + dt * sigma(3,1))) +(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(3,3) + dt * sigma(3,3)))/((-((2 * eps0 * epr(1,3)      &
      + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 *      &
      eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0      &
      * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0      &
      * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehyz = ((2 * eps0 * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1)) -(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3)))/((-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0 * epr(1,2) + dt * sigma      &
      (1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehzx = (-((2 * eps0 * epr(2,2) + dt * sigma(2,2)) *(2 * eps0 * epr(3,1) + dt * sigma(3,1))) +(2      &
      * eps0 * epr(2,1) + dt * sigma(2,1)) *(2 * eps0 * epr(3,2) + dt * sigma(3,2)))/((-((2 * eps0 * epr(1,3)      &
      + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 *      &
      eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0      &
      * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0      &
      * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehzy = ((2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(3,1) + dt * sigma(3,1)) -(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(3,2) + dt * sigma(3,2)))/((-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0 * epr(1,3) + dt * sigma      &
      (1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0 * epr(1,2) + dt * sigma      &
      (1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr      &
      (2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%ehzz = (-((2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2      &
      * eps0 * epr(1,1) + dt * sigma(1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2)))/((-((2 * eps0 * epr(1,3)      &
      + dt * sigma(1,3)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) +(2 * eps0 * epr(1,2) + dt * sigma(1,2)) *(2 *      &
      eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,1))/dt + sigma(3,1)/2.0_RKIND) -(-((2 * eps0      &
      * epr(1,3) + dt * sigma(1,3)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,3) + dt * sigma(2,3))) *((eps0 * epr(3,2))/dt + sigma(3,2)/2.0_RKIND) +(-((2 * eps0      &
      * epr(1,2) + dt * sigma(1,2)) *(2 * eps0 * epr(2,1) + dt * sigma(2,1))) +(2 * eps0 * epr(1,1) + dt * sigma      &
      (1,1)) *(2 * eps0 * epr(2,2) + dt * sigma(2,2))) *((eps0 * epr(3,3))/dt + sigma(3,3)/2.0_RKIND))

      coeff%hhxx = ((-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2      &
      * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt -      &
      sigmam(3,1)/2.0_RKIND) +((mu0 * mur(2,1))/dt - sigmam(2,1)/2.0_RKIND) *((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2      &
      * mu0 * mur(3,2) + dt * sigmam(3,2)) -(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(3,3) +      &
      dt * sigmam(3,3))) +((mu0 * mur(1,1))/dt - sigmam(1,1)/2.0_RKIND) *(-((2 * mu0 * mur(2,3) + dt * sigmam(2,3)) *(2 *      &
      mu0 * mur(3,2) + dt * sigmam(3,2))) +(2 * mu0 * mur(2,2) + dt * sigmam(2,2)) *(2 * mu0 * mur(3,3) + dt      &
      * sigmam(3,3))))/((-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2)))      &
      +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt +      &
      sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3)))      &
      *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt      &
      * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2)))      &
      *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hhxy = (4 * dt * mu0 *(2 * mu0 *((mur(1,3) * mur(3,2) -      &
      mur(1,2) * mur(3,3)) * sigmam(2,2) + mur(2,3) *(-(mur(3,2) * sigmam(1,2)) + mur(1,2) * sigmam(3,2)) + mur      &
      (2,2) *(mur(3,3) * sigmam(1,2) - mur(1,3) * sigmam(3,2))) + dt *(mur(3,2) *(sigmam(1,3)      &
      * sigmam(2,2) - sigmam(1,2) * sigmam(2,3)) + mur(2,2) *(-(sigmam(1,3) * sigmam(3,2)) + sigmam(1,2) *      &
      sigmam(3,3)) + mur(1,2) *(sigmam(2,3) * sigmam(3,2) - sigmam(2,2) * sigmam(3,3)))))/(8 * mu0**3 *(mur(1,3)      &
      *(mur(2,2) * mur(3,1) - mur(2,1) * mur(3,2)) + mur(1,2) *(-(mur(2,3) * mur(3,1)) + mur(2,1) * mur(3,3)) + mur      &
      (1,1) *(mur(2,3) * mur(3,2) - mur(2,2) * mur(3,3))) + 2 * dt**2.0_RKIND * mu0 *(mur(3,1)      &
      * sigmam(1,3) * sigmam(2,2) + mur(3,3) *(sigmam(1,2) * sigmam(2,1) - sigmam(1,1) * sigmam(2,2)) - mur      &
      (3,1) * sigmam(1,2) * sigmam(2,3) + mur(3,2) *(-(sigmam(1,3) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,3)) -      &
      mur(2,3) * sigmam(1,2) * sigmam(3,1) + mur(2,2) * sigmam(1,3) * sigmam(3,1) + mur(1,3) * sigmam(2,2)      &
      * sigmam(3,1) - mur(1,2) * sigmam(2,3) * sigmam(3,1) + mur(2,3) * sigmam(1,1) * sigmam(3,2) - mur      &
      (2,1) * sigmam(1,3) * sigmam(3,2) - mur(1,3) * sigmam(2,1) * sigmam(3,2) + mur(1,1) * sigmam(2,3)      &
      * sigmam(3,2) +(-(mur(2,2) * sigmam(1,1)) + mur(2,1) * sigmam(1,2) + mur(1,2) * sigmam(2,1) - mur      &
      (1,1) * sigmam(2,2)) * sigmam(3,3)) + 4 * dt * mu0**2.0_RKIND *(mur(2,1) * mur(3,3) * sigmam(1,2) - mur(2,1) * mur      &
      (3,2) * sigmam(1,3) - mur(1,3) * mur(3,2) * sigmam(2,1) + mur(1,2) * mur(3,3) * sigmam(2,1) + mur(1,3)      &
      * mur(3,1) * sigmam(2,2) - mur(1,1) * mur(3,3) * sigmam(2,2) - mur(1,2) * mur(3,1) * sigmam(2,3) +      &
      mur(1,1) * mur(3,2) * sigmam(2,3) - mur(1,3) * mur(2,1) * sigmam(3,2) + mur(2,3) *(mur(3,2)      &
      * sigmam(1,1) - mur(3,1) * sigmam(1,2) - mur(1,2) * sigmam(3,1) + mur(1,1) * sigmam(3,2)) + mur(1,2)      &
      * mur(2,1) * sigmam(3,3) + mur(2,2) *(-(mur(3,3) * sigmam(1,1)) + mur(3,1) * sigmam(1,3) + mur(1,3)      &
      * sigmam(3,1) - mur(1,1) * sigmam(3,3))) + dt**3 *(sigmam(1,3) *(sigmam(2,2) * sigmam(3,1) - sigmam      &
      (2,1) * sigmam(3,2)) + sigmam(1,2) *(-(sigmam(2,3) * sigmam(3,1)) + sigmam(2,1) * sigmam(3,3)) + sigmam(1,1)      &
      *(sigmam(2,3) * sigmam(3,2) - sigmam(2,2) * sigmam(3,3))))

      coeff%hhxz = (4 * dt * mu0 *(2 * mu0 *((mur      &
      (1,3) * mur(3,2) - mur(1,2) * mur(3,3)) * sigmam(2,3) + mur(2,3) *(-(mur(3,2) * sigmam(1,3)) + mur(1,2)      &
      * sigmam(3,3)) + mur(2,2) *(mur(3,3) * sigmam(1,3) - mur(1,3)      &
      * sigmam(3,3))) + dt *(mur(3,3) *(sigmam(1,3) * sigmam(2,2) - sigmam(1,2) * sigmam(2,3)) + mur      &
      (2,3) *(-(sigmam(1,3) * sigmam(3,2)) + sigmam(1,2) * sigmam(3,3)) + mur(1,3) *(sigmam(2,3) * sigmam(3,2) -      &
      sigmam(2,2) * sigmam(3,3)))))/(8 * mu0**3 *(mur(1,3) *(mur(2,2) * mur(3,1) - mur(2,1) * mur(3,2)) + mur(1,2)      &
      *(-(mur(2,3) * mur(3,1)) + mur(2,1) * mur(3,3)) + mur(1,1) *(mur(2,3) * mur(3,2) - mur(2,2) * mur(3,3))) + 2      &
      * dt**2.0_RKIND * mu0 *(mur(3,1) * sigmam(1,3) * sigmam(2,2) + mur(3,3)      &
      *(sigmam(1,2) * sigmam(2,1) - sigmam(1,1) * sigmam(2,2)) - mur(3,1) * sigmam(1,2) * sigmam(2,3) + mur(3,2)      &
      *(-(sigmam(1,3) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,3)) - mur(2,3) * sigmam(1,2) * sigmam(3,1) + mur      &
      (2,2) * sigmam(1,3) * sigmam(3,1) + mur(1,3) * sigmam(2,2) * sigmam(3,1) - mur(1,2) * sigmam(2,3)      &
      * sigmam(3,1) + mur(2,3) * sigmam(1,1) * sigmam(3,2) - mur(2,1) * sigmam(1,3) * sigmam(3,2) - mur      &
      (1,3) * sigmam(2,1) * sigmam(3,2) + mur(1,1) * sigmam(2,3) * sigmam(3,2) +(-(mur(2,2) * sigmam(1,1)) + mur      &
      (2,1) * sigmam(1,2) + mur(1,2) * sigmam(2,1) - mur(1,1) * sigmam(2,2)) * sigmam(3,3)) + 4 * dt * mu0**2.0_RKIND *(mur      &
      (2,1) * mur(3,3) * sigmam(1,2) - mur(2,1) * mur(3,2) * sigmam(1,3) - mur(1,3) * mur(3,2)      &
      * sigmam(2,1) + mur(1,2) * mur(3,3) * sigmam(2,1) + mur(1,3) * mur(3,1) * sigmam(2,2) - mur(1,1) *      &
      mur(3,3) * sigmam(2,2) - mur(1,2) * mur(3,1) * sigmam(2,3) + mur(1,1) * mur(3,2) * sigmam(2,3) - mur(1,3)      &
      * mur(2,1) * sigmam(3,2) + mur(2,3) *(mur(3,2) * sigmam(1,1) - mur(3,1) * sigmam(1,2) - mur(1,2) *      &
      sigmam(3,1) + mur(1,1) * sigmam(3,2)) + mur(1,2) * mur(2,1) * sigmam(3,3) + mur(2,2) *(-(mur(3,3) * sigmam      &
      (1,1)) + mur(3,1) * sigmam(1,3) + mur(1,3) * sigmam(3,1) - mur(1,1) * sigmam(3,3))) + dt**3 *(sigmam(1,3)      &
      *(sigmam(2,2) * sigmam(3,1) - sigmam(2,1) * sigmam(3,2)) + sigmam(1,2) *(-(sigmam(2,3)      &
      * sigmam(3,1)) + sigmam(2,1) * sigmam(3,3)) + sigmam(1,1) *(sigmam(2,3)      &
      * sigmam(3,2) - sigmam(2,2) * sigmam(3,3))))

      coeff%hhyx = (4 * dt * mu0 *(2 * mu0 *((mur(1,3) * mur(3,1) - mur      &
      (1,1) * mur(3,3)) * sigmam(2,1) + mur(2,3) *(-(mur(3,1) * sigmam(1,1)) + mur(1,1) * sigmam(3,1)) + mur(2,1)      &
      *(mur(3,3) * sigmam(1,1) - mur(1,3) * sigmam(3,1))) + dt *(mur(3,1) *(sigmam(1,3) * sigmam(2,1) - sigmam      &
      (1,1) * sigmam(2,3)) + mur(2,1) *(-(sigmam(1,3) * sigmam(3,1)) + sigmam(1,1) * sigmam(3,3)) + mur(1,1)      &
      *(sigmam(2,3) * sigmam(3,1) - sigmam(2,1) * sigmam(3,3)))))/(8 * mu0**3 *(mur(1,3) *(-(mur(2,2) * mur(3,1)) +      &
      mur(2,1) * mur(3,2)) + mur(1,2) *(mur(2,3) * mur(3,1) - mur(2,1) * mur(3,3)) + mur(1,1) *(-(mur(2,3)      &
      * mur(3,2)) + mur(2,2) * mur(3,3))) + 2 * dt**2.0_RKIND * mu0 *(-(mur(3,1)      &
      * sigmam(1,3) * sigmam(2,2)) + mur(3,3) *(-(sigmam(1,2) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,2)) +      &
      mur(3,1) * sigmam(1,2) * sigmam(2,3) + mur(3,2) *(sigmam(1,3) * sigmam(2,1) - sigmam(1,1) * sigmam(2,3)) +      &
      mur(2,3) * sigmam(1,2) * sigmam(3,1) - mur(2,2) * sigmam(1,3) * sigmam(3,1) - mur(1,3)      &
      * sigmam(2,2) * sigmam(3,1) + mur(1,2) * sigmam(2,3) * sigmam(3,1) - mur(2,3) * sigmam(1,1) * sigmam(3,2) + mur      &
      (2,1) * sigmam(1,3) * sigmam(3,2) + mur(1,3) * sigmam(2,1) * sigmam(3,2) - mur(1,1) * sigmam(2,3) * sigmam      &
      (3,2) +(mur(2,2) * sigmam(1,1) - mur(2,1) * sigmam(1,2) - mur(1,2) * sigmam(2,1) + mur(1,1) * sigmam(2,2))      &
      * sigmam(3,3)) + 4 * dt * mu0**2.0_RKIND *(-(mur(2,1) * mur(3,3) * sigmam(1,2)) + mur(2,1) * mur(3,2) * sigmam      &
      (1,3) + mur(1,3) * mur(3,2) * sigmam(2,1) - mur(1,2) * mur(3,3) * sigmam(2,1) - mur(1,3) * mur(3,1) * sigmam      &
      (2,2) + mur(1,1) * mur(3,3) * sigmam(2,2) + mur(1,2) * mur(3,1) * sigmam(2,3) - mur(1,1) * mur(3,2) * sigmam      &
      (2,3) + mur(1,3) * mur(2,1) * sigmam(3,2) + mur(2,3) *(-(mur(3,2) * sigmam(1,1)) + mur(3,1) * sigmam(1,2) +      &
      mur(1,2) * sigmam(3,1) - mur(1,1) * sigmam(3,2)) - mur(1,2) * mur(2,1) * sigmam(3,3) + mur(2,2) *(mur(3,3)      &
      * sigmam(1,1) - mur(3,1) * sigmam(1,3) - mur(1,3) * sigmam(3,1) + mur(1,1) * sigmam(3,3))) + dt**3 *      &
      (sigmam(1,3) *(-(sigmam(2,2) * sigmam(3,1)) + sigmam(2,1) * sigmam(3,2)) + sigmam(1,2) *(sigmam(2,3)      &
      * sigmam(3,1) - sigmam(2,1) * sigmam(3,3)) + sigmam(1,1) *(-(sigmam(2,3) * sigmam(3,2)) + sigmam(2,2)      &
      * sigmam(3,3))))

      coeff%hhyy = (((2 * mu0      &
      * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1)) -(2 * mu0 * mur(1,1) + dt *      &
      sigmam(1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt - sigmam(3,2)/2.0_RKIND) +((mu0      &
      * mur(2,2))/dt - sigmam(2,2)/2.0_RKIND) *(-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(3,1) + dt      &
      * sigmam(3,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(3,3) + dt      &
      * sigmam(3,3))) +((mu0 * mur(1,2))/dt - sigmam(1,2)/2.0_RKIND) *((2 * mu0 * mur(2,3) + dt * sigmam(2,3)) *(2 * mu0 * mur  &
      (3,1) + dt * sigmam(3,1)) -(2 * mu0 * mur(2,1) + dt * sigmam(2,1)) *(2 * mu0 * mur(3,3) + dt * sigmam      &
      (3,3))))/((-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt +      &
      sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt +      &
      sigmam(3,2)/2.0_RKIND) +(-((2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt +      &
      sigmam(3,3)/2.0_RKIND))

      coeff%hhyz = (4 * dt * mu0 *(2 * mu0 *((-(mur(1,3) * mur(3,1)) + mur(1,1) * mur(3,3)) * sigmam(2,3) + mur      &
      (2,3) *(mur(3,1) * sigmam(1,3) - mur(1,1) * sigmam(3,3)) + mur(2,1) *(-(mur(3,3) * sigmam(1,3)) + mur(1,3)      &
      * sigmam(3,3))) + dt *(mur(3,3) *(-(sigmam(1,3) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,3)) + mur      &
      (2,3) *(sigmam(1,3) * sigmam(3,1) - sigmam(1,1) * sigmam(3,3)) + mur(1,3) *(-(sigmam(2,3) * sigmam(3,1)) +      &
      sigmam(2,1) * sigmam(3,3)))))/(8 * mu0**3 *(mur(1,3) *(mur(2,2) * mur(3,1) - mur(2,1)      &
      * mur(3,2)) + mur(1,2) *(-(mur(2,3) * mur(3,1)) + mur(2,1) * mur(3,3)) + mur(1,1) *(mur(2,3) * mur      &
      (3,2) - mur(2,2) * mur(3,3))) + 2 * dt**2.0_RKIND * mu0 *(mur(3,1) * sigmam(1,3) * sigmam(2,2) + mur(3,3) *(sigmam      &
      (1,2) * sigmam(2,1) - sigmam(1,1) * sigmam(2,2)) - mur(3,1) * sigmam(1,2) * sigmam(2,3) + mur(3,2) *      &
      (-(sigmam(1,3) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,3)) - mur(2,3)      &
      * sigmam(1,2) * sigmam(3,1) + mur(2,2) * sigmam(1,3) * sigmam(3,1) + mur(1,3) * sigmam(2,2) * sigmam(3,1) - mur      &
      (1,2) * sigmam(2,3) * sigmam(3,1) + mur(2,3) * sigmam(1,1) * sigmam(3,2) - mur(2,1) * sigmam(1,3) * sigmam      &
      (3,2) - mur(1,3) * sigmam(2,1) * sigmam(3,2) + mur(1,1) * sigmam(2,3) * sigmam(3,2) +(-(mur(2,2) * sigmam      &
      (1,1)) + mur(2,1) * sigmam(1,2) + mur(1,2) * sigmam(2,1) - mur(1,1) * sigmam(2,2)) * sigmam(3,3)) + 4      &
      * dt * mu0**2.0_RKIND *(mur(2,1) * mur(3,3) * sigmam(1,2) - mur(2,1) * mur(3,2) * sigmam(1,3) - mur(1,3) *      &
      mur(3,2) * sigmam(2,1) + mur(1,2) * mur(3,3) * sigmam(2,1) + mur(1,3) * mur(3,1) * sigmam(2,2) - mur(1,1)      &
      * mur(3,3) * sigmam(2,2) - mur(1,2) * mur(3,1) * sigmam(2,3) + mur(1,1) * mur(3,2) * sigmam(2,3) -      &
      mur(1,3) * mur(2,1) * sigmam(3,2) + mur(2,3) *(mur(3,2) * sigmam(1,1) - mur(3,1) * sigmam(1,2) - mur(1,2)      &
      * sigmam(3,1) + mur(1,1) * sigmam(3,2)) + mur(1,2) * mur(2,1) * sigmam(3,3) + mur(2,2) *(-(mur(3,3)      &
      * sigmam(1,1)) + mur(3,1) * sigmam(1,3) + mur(1,3) * sigmam(3,1) - mur(1,1) * sigmam(3,3))) + dt**3 *      &
      (sigmam(1,3) *(sigmam(2,2) * sigmam(3,1) - sigmam(2,1) * sigmam(3,2)) + sigmam(1,2) *(-(sigmam(2,3)      &
      * sigmam(3,1)) + sigmam(2,1) * sigmam(3,3)) + sigmam(1,1) *(sigmam(2,3) * sigmam(3,2) - sigmam(2,2) *      &
      sigmam(3,3))))

      coeff%hhzx = (4 * dt * mu0 *(2 * mu0 *((mur(1,2) * mur(3,1) - mur(1,1) * mur(3,2))      &
      * sigmam(2,1) + mur(2,2) *(-(mur(3,1) * sigmam(1,1)) + mur(1,1) * sigmam(3,1)) + mur(2,1) *(mur      &
      (3,2) * sigmam(1,1) - mur(1,2) * sigmam(3,1))) + dt *(mur(3,1) *(sigmam(1,2) * sigmam(2,1) - sigmam(1,1) *      &
      sigmam(2,2)) + mur(2,1) *(-(sigmam(1,2) * sigmam(3,1)) + sigmam(1,1) * sigmam(3,2)) + mur(1,1) *(sigmam(2,2)      &
      * sigmam(3,1) - sigmam(2,1) * sigmam(3,2)))))/(8      &
      * mu0**3 *(mur(1,3) *(mur(2,2) * mur(3,1) - mur(2,1) * mur(3,2)) + mur(1,2) *(-(mur(2,3) * mur      &
      (3,1)) + mur(2,1) * mur(3,3)) + mur(1,1) *(mur(2,3) * mur(3,2) - mur(2,2) * mur(3,3))) + 2 * dt**2.0_RKIND * mu0      &
      *(mur(3,1) * sigmam(1,3) * sigmam(2,2) + mur(3,3) *(sigmam(1,2) * sigmam(2,1) - sigmam(1,1)      &
      * sigmam(2,2)) - mur(3,1) * sigmam(1,2) * sigmam(2,3) + mur(3,2) *(-(sigmam(1,3) * sigmam(2,1)) +      &
      sigmam(1,1) * sigmam(2,3)) - mur(2,3) * sigmam(1,2) * sigmam(3,1) + mur(2,2) * sigmam(1,3) * sigmam(3,1) +      &
      mur(1,3) * sigmam(2,2) * sigmam(3,1) - mur(1,2) * sigmam(2,3) * sigmam(3,1) + mur(2,3)      &
      * sigmam(1,1) * sigmam(3,2) - mur(2,1) * sigmam(1,3) * sigmam(3,2) - mur(1,3) * sigmam(2,1) * sigmam(3,2) + mur      &
      (1,1) * sigmam(2,3) * sigmam(3,2) +(-(mur(2,2) * sigmam(1,1)) + mur(2,1) * sigmam(1,2) + mur(1,2) * sigmam      &
      (2,1) - mur(1,1) * sigmam(2,2)) * sigmam(3,3)) + 4 * dt * mu0**2.0_RKIND *(mur(2,1) * mur(3,3) * sigmam(1,2) - mur      &
      (2,1) * mur(3,2) * sigmam(1,3) - mur(1,3) * mur(3,2) * sigmam(2,1) + mur(1,2) * mur(3,3) * sigmam(2,1) + mur      &
      (1,3) * mur(3,1) * sigmam(2,2) - mur(1,1) * mur(3,3) * sigmam(2,2) - mur(1,2) * mur(3,1) * sigmam(2,3) + mur      &
      (1,1) * mur(3,2) * sigmam(2,3) - mur(1,3) * mur(2,1) * sigmam(3,2) + mur(2,3) *(mur(3,2)      &
      * sigmam(1,1) - mur(3,1) * sigmam(1,2) - mur(1,2) * sigmam(3,1) + mur(1,1) * sigmam(3,2)) + mur(1,2)      &
      * mur(2,1) * sigmam(3,3) + mur(2,2) *(-(mur(3,3) * sigmam(1,1)) + mur(3,1) * sigmam(1,3) + mur(1,3)      &
      * sigmam(3,1) - mur(1,1) * sigmam(3,3))) + dt**3 *(sigmam(1,3) *(sigmam(2,2) * sigmam(3,1) - sigmam      &
      (2,1) * sigmam(3,2)) + sigmam(1,2) *(-(sigmam(2,3) * sigmam(3,1)) + sigmam(2,1) * sigmam(3,3)) + sigmam(1,1)      &
      *(sigmam(2,3) * sigmam(3,2) - sigmam(2,2) * sigmam(3,3))))

      coeff%hhzy = (4 * dt * mu0 *(2 * mu0 *((-(mur      &
      (1,2) * mur(3,1)) + mur(1,1) * mur(3,2)) * sigmam(2,2) + mur(2,2) *(mur(3,1)      &
      * sigmam(1,2) - mur(1,1) * sigmam(3,2)) + mur(2,1) *(-(mur(3,2) * sigmam(1,2)) + mur(1,2)      &
      * sigmam(3,2))) + dt *(mur(3,2) *(-(sigmam(1,2) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,2)) + mur      &
      (2,2) *(sigmam(1,2) * sigmam(3,1) - sigmam(1,1) * sigmam(3,2)) + mur(1,2) *(-(sigmam(2,2) * sigmam(3,1)) +      &
      sigmam(2,1) * sigmam(3,2)))))/(8 * mu0**3 *(mur(1,3) *(-(mur(2,2)      &
      * mur(3,1)) + mur(2,1) * mur(3,2)) + mur(1,2) *(mur(2,3) * mur(3,1) - mur(2,1) * mur(3,3)) + mur      &
      (1,1) *(-(mur(2,3) * mur(3,2)) + mur(2,2) * mur(3,3))) + 2 * dt**2.0_RKIND * mu0 *(-(mur(3,1) * sigmam(1,3) * sigmam      &
      (2,2)) + mur(3,3) *(-(sigmam(1,2) * sigmam(2,1)) + sigmam(1,1) * sigmam(2,2)) + mur(3,1) * sigmam(1,2)      &
      * sigmam(2,3) + mur(3,2) *(sigmam(1,3) * sigmam(2,1) - sigmam(1,1) * sigmam(2,3)) + mur(2,3) * sigmam      &
      (1,2) * sigmam(3,1) - mur(2,2) * sigmam(1,3) * sigmam(3,1) - mur(1,3) * sigmam(2,2) * sigmam(3,1) + mur(1,2)      &
      * sigmam(2,3) * sigmam(3,1) - mur(2,3) * sigmam(1,1) * sigmam(3,2) + mur(2,1)      &
      * sigmam(1,3) * sigmam(3,2) + mur(1,3) * sigmam(2,1) * sigmam(3,2) - mur(1,1) * sigmam(2,3)      &
      * sigmam(3,2) +(mur(2,2) * sigmam(1,1) - mur(2,1) * sigmam(1,2) - mur(1,2) * sigmam(2,1) + mur(1,1) * sigmam      &
      (2,2)) * sigmam(3,3)) + 4 * dt * mu0**2.0_RKIND *(-(mur(2,1) * mur(3,3) * sigmam(1,2)) + mur(2,1) * mur(3,2) * sigmam      &
      (1,3) + mur(1,3) * mur(3,2) * sigmam(2,1) - mur(1,2) * mur(3,3) * sigmam(2,1) - mur(1,3) * mur(3,1) * sigmam      &
      (2,2) + mur(1,1) * mur(3,3) * sigmam(2,2) + mur(1,2) * mur(3,1) * sigmam(2,3) - mur(1,1) * mur(3,2) * sigmam      &
      (2,3) + mur(1,3) * mur(2,1) * sigmam(3,2) + mur(2,3) *(-(mur(3,2) * sigmam(1,1)) + mur(3,1) * sigmam(1,2) +      &
      mur(1,2) * sigmam(3,1) - mur(1,1) * sigmam(3,2)) - mur(1,2) * mur(2,1) * sigmam(3,3) + mur(2,2) *(mur(3,3)      &
      * sigmam(1,1) - mur(3,1) * sigmam(1,3) - mur(1,3) * sigmam(3,1) + mur(1,1) * sigmam(3,3))) + dt**3 *      &
      (sigmam(1,3) *(-(sigmam(2,2) * sigmam(3,1)) + sigmam(2,1) * sigmam(3,2)) + sigmam(1,2) *(sigmam(2,3)      &
      * sigmam(3,1) - sigmam(2,1) * sigmam(3,3)) + sigmam(1,1) *(-(sigmam(2,3) * sigmam(3,2)) + sigmam(2,2)      &
      * sigmam(3,3))))

      coeff%hhzz = (((mu0 * mur(2,3))/dt - sigmam(2,3)/2.0_RKIND) *((2 * mu0 * mur(1,2) + dt *      &
      sigmam(1,2)) *(2 * mu0 * mur(3,1) + dt * sigmam(3,1)) -(2 * mu0 * mur(1,1) + dt      &
      * sigmam(1,1)) *(2 * mu0 * mur(3,2) + dt * sigmam(3,2))) +((mu0 * mur(1,3))/dt - sigmam(1,3)/2.0_RKIND)      &
      *(-((2 * mu0 * mur(2,2) + dt * sigmam(2,2)) *(2 * mu0 * mur(3,1) + dt * sigmam(3,1))) +(2 * mu0 * mur(2,1) +      &
      dt * sigmam(2,1)) *(2 * mu0 * mur(3,2) + dt * sigmam(3,2))) +(-((2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 *      &
      mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,2) + dt      &
      * sigmam(2,2))) *((mu0 * mur(3,3))/dt - sigmam(3,3)/2.0_RKIND))/((-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *      &
      (2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,3) + dt      &
      * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam(1,3))      &
      *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,3) + dt      &
      * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt *      &
      sigmam(1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hexx = (-((2 * mu0 * mur(2,3) + dt * sigmam(2,3)) *(2 * mu0 * mur(3,2) + dt * sigmam(3,2))) +(2      &
      * mu0 * mur(2,2) + dt * sigmam(2,2)) *(2 * mu0 * mur(3,3) + dt * sigmam(3,3)))/((-((2 * mu0 * mur(1,3)      &
      + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 *      &
      mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0      &
      * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hexy = ((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(3,2) + dt * sigmam(3,2)) -(2      &
      * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(3,3) + dt * sigmam(3,3)))/((-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0 * mur(1,2) + dt * sigmam      &
      (1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hexz = (-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2      &
      * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3)))/((-((2 * mu0 * mur(1,3)      &
      + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 *      &
      mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0      &
      * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%heyx = ((2 * mu0 * mur(2,3) + dt * sigmam(2,3)) *(2 * mu0 * mur(3,1) + dt * sigmam(3,1)) -(2      &
      * mu0 * mur(2,1) + dt * sigmam(2,1)) *(2 * mu0 * mur(3,3) + dt * sigmam(3,3)))/((-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0 * mur(1,2) + dt * sigmam      &
      (1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%heyy = (-((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(3,1) + dt * sigmam(3,1))) +(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(3,3) + dt * sigmam(3,3)))/((-((2 * mu0 * mur(1,3)      &
      + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 *      &
      mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0      &
      * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%heyz = ((2 * mu0 * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1)) -(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3)))/((-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0 * mur(1,2) + dt * sigmam      &
      (1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hezx = (-((2 * mu0 * mur(2,2) + dt * sigmam(2,2)) *(2 * mu0 * mur(3,1) + dt * sigmam(3,1))) +(2      &
      * mu0 * mur(2,1) + dt * sigmam(2,1)) *(2 * mu0 * mur(3,2) + dt * sigmam(3,2)))/((-((2 * mu0 * mur(1,3)      &
      + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 *      &
      mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0      &
      * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hezy = ((2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(3,1) + dt * sigmam(3,1)) -(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(3,2) + dt * sigmam(3,2)))/((-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0 * mur(1,3) + dt * sigmam      &
      (1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0 * mur(1,2) + dt * sigmam      &
      (1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur      &
      (2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      coeff%hezz = (-((2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2      &
      * mu0 * mur(1,1) + dt * sigmam(1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2)))/((-((2 * mu0 * mur(1,3)      &
      + dt * sigmam(1,3)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) +(2 * mu0 * mur(1,2) + dt * sigmam(1,2)) *(2 *      &
      mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,1))/dt + sigmam(3,1)/2.0_RKIND) -(-((2 * mu0      &
      * mur(1,3) + dt * sigmam(1,3)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,3) + dt * sigmam(2,3))) *((mu0 * mur(3,2))/dt + sigmam(3,2)/2.0_RKIND) +(-((2 * mu0      &
      * mur(1,2) + dt * sigmam(1,2)) *(2 * mu0 * mur(2,1) + dt * sigmam(2,1))) +(2 * mu0 * mur(1,1) + dt * sigmam      &
      (1,1)) *(2 * mu0 * mur(2,2) + dt * sigmam(2,2))) *((mu0 * mur(3,3))/dt + sigmam(3,3)/2.0_RKIND))

      return

   end subroutine

end module Anisotropic
