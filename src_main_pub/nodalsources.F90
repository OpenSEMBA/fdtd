
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Nodal  Feeding modules
!  Creation date Date :  June0, 21, 2011
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module nodalsources

   use fdetypes
   USE REPORT

   IMPLICIT NONE
   private

   type XYZlimit_t_singlescaled
      integer (kind=4)  :: XI,XE,YI,YE,ZI,ZE
      REAL (KIND=RKIND)   :: amplitude
   end type


   type  ::  NodalLocal_t
      REAL (KIND=RKIND), pointer, dimension (:)  ::  evol
      REAL (KIND=RKIND)   :: deltaevol
      integer (kind=4) :: numus
      type (XYZlimit_t_singlescaled)   :: punto
      logical :: IsInitialValue
   end type NodalLocal_t


   type :: nodsou
      integer (kind=4) :: NumHard = 0 , NumSoft = 0
      type (NodalLocal_t), pointer, dimension(:) :: nodHard,nodSoft
   end type

   !!!!!variables locales

   type (nodsou), save, target :: Nodal_Ex,Nodal_Ey,Nodal_Ez
   type (nodsou), save, target :: Nodal_Hx,Nodal_Hy,Nodal_Hz

   public :: InitHopf,initNodalSources,AdvanceNodalE,AdvanceNodalH,DestroyNodal,nodsou,getnodal




contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes Nodal Source data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitnodalSources(sgg,layoutnumber,NumNodalSources,sggNodalSource,sggSweep,ThereAreNodalE,ThereAreNodalH)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      !!!
      integer, intent (in) :: NumNodalSources
      type (NodalSource_t), dimension(1:NumNodalSources),intent(in)           ::  sggNodalSource

      integer (kind=4):: layoutnumber,j,i
      logical, intent(out)  ::  ThereArenodalE,ThereArenodalH
      integer :: numNodalSoft_Ex,numNodalSoft_Ey,numNodalSoft_Ez, &
      numNodalSoft_Hx,numNodalSoft_Hy,numNodalSoft_Hz, &
      numNodalHard_Ex,numNodalHard_Ey,numNodalHard_Ez, &
      numNodalHard_Hx,numNodalHard_Hy,numNodalHard_Hz

      real (kind=rkind) :: amplit
      type (XYZlimit_t), dimension (1:6)    ::  sggSweep

      !!!

      numNodalSoft_Ex = 0
      numNodalSoft_Ey = 0
      numNodalSoft_Ez = 0
      numNodalSoft_Hx = 0
      numNodalSoft_Hy = 0
      numNodalSoft_Hz = 0
      numNodalHard_Ex = 0
      numNodalHard_Ey = 0
      numNodalHard_Ez = 0
      numNodalHard_Hx = 0
      numNodalHard_Hy = 0
      numNodalHard_Hz = 0

      ThereArenodalE=.false.
      ThereArenodalH=.false.
      
      do j=1,NumNodalSources
         if (sggNodalSource(j)%IsElec) then
            do i=1,sggNodalSource(j)%numpuntos
               if (sggNodalSource(j)%punto(i)%xc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Ex = numNodalHard_Ex  +1
                  else
                     numNodalSoft_Ex = numNodalSoft_Ex  +1
                  endif
               endif
               if (sggNodalSource(j)%punto(i)%yc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Ey = numNodalHard_Ey  +1
                  else
                     numNodalSoft_Ey = numNodalSoft_Ey  +1
                  endif
               endif
               if (sggNodalSource(j)%punto(i)%zc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Ez = numNodalHard_Ez  +1
                  else
                     numNodalSoft_Ez = numNodalSoft_Ez  +1
                  endif
               endif
            end do
         else
            do i=1,sggNodalSource(j)%numpuntos
               if (sggNodalSource(j)%punto(i)%xc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Hx = numNodalHard_Hx  +1
                  else
                     numNodalSoft_Hx = numNodalSoft_Hx  +1
                  endif
               endif
               if (sggNodalSource(j)%punto(i)%yc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Hy = numNodalHard_Hy  +1
                  else
                     numNodalSoft_Hy = numNodalSoft_Hy  +1
                  endif
               endif
               if (sggNodalSource(j)%punto(i)%zc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Hz = numNodalHard_Hz  +1
                  else
                     numNodalSoft_Hz = numNodalSoft_Hz  +1
                  endif
               endif
            end do
         endif
      end do


      if (NumNodalSoft_Ex+NumNodalSoft_Ey+NumNodalSoft_Ez /= 0) then
         ThereArenodalE =.true.
         ALLOCATE (Nodal_Ex%nodSoft(1:numNodalSoft_Ex), &
         Nodal_Ey%nodSoft(1:numNodalSoft_Ey), &
         Nodal_Ez%nodSoft(1:numNodalSoft_Ez))
      endif
      if (NumNodalHard_Ex+NumNodalHard_Ey+NumNodalHard_Ez /= 0) then
         ThereArenodalE =.true.
         ALLOCATE (Nodal_Ex%nodHard(1:numNodalHard_Ex), &
         Nodal_Ey%nodHard(1:numNodalHard_Ey), &
         Nodal_Ez%nodHard(1:numNodalHard_Ez))
      endif
      if (NumNodalSoft_Hx+NumNodalSoft_Hy+NumNodalSoft_Hz /= 0) then
         ThereArenodalH =.true.
         ALLOCATE (Nodal_Hx%nodSoft(1:numNodalSoft_Hx), &
         Nodal_Hy%nodSoft(1:numNodalSoft_Hy), &
         Nodal_Hz%nodSoft(1:numNodalSoft_Hz))
      endif
      if (NumNodalHard_Hx+NumNodalHard_Hy+NumNodalHard_Hz /= 0) then
         ThereArenodalH =.true.
         ALLOCATE (Nodal_Hx%nodHard(1:numNodalHard_Hx), &
         Nodal_Hy%nodHard(1:numNodalHard_Hy), &
         Nodal_Hz%nodHard(1:numNodalHard_Hz))
      endif


      Nodal_Ex%numHard = 0
      Nodal_Ey%numHard = 0
      Nodal_Ez%numHard = 0
      Nodal_Hx%numHard = 0
      Nodal_Hy%numHard = 0
      Nodal_Hz%numHard = 0
      !
      Nodal_Ex%numSoft = 0
      Nodal_Ey%numSoft = 0
      Nodal_Ez%numSoft = 0
      Nodal_Hx%numSoft = 0
      Nodal_Hy%numSoft = 0
      Nodal_Hz%numSoft = 0


      do j=1,NumNodalSources
         if (sggNodalSource(j)%IsElec) then
            do i=1,sggNodalSource(j)%numpuntos
               amplit = sggNodalSource(J)%punto(i)%xc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Ex,sggNodalSource(J),sggSweep(iEx),i,amplit)
               endif
               amplit = sggNodalSource(j)%punto(i)%yc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Ey,sggNodalSource(J),sggSweep(iEy),i,amplit)
               endif
               amplit = sggNodalSource(j)%punto(i)%zc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Ez,sggNodalSource(J),sggSweep(iEz),i,amplit)
               endif
            end do
         else !es magnetico
            do i=1,sggNodalSource(j)%numpuntos
               amplit = sggNodalSource(J)%punto(i)%xc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Hx,sggNodalSource(J),sggSweep(iHx),i,amplit)
               endif
               amplit = sggNodalSource(j)%punto(i)%yc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Hy,sggNodalSource(J),sggSweep(iHy),i,amplit)
               endif
               amplit = sggNodalSource(j)%punto(i)%zc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Hz,sggNodalSource(J),sggSweep(iHz),i,amplit)
               endif
            end do
         endif
      end do

      !print *,'tras nodal sources ',Nodal_Ex%numHard

      return


   contains

      subroutine createnodal(layoutnumber,dummy,sggdummy,sggSweep,index,amplit)

         type (nodsou), intent (inout) :: dummy
         type (NodalSource_t), intent (in), target :: sggdummy
         real (kind=rkind), intent(in) :: amplit
         integer (kind=4), intent(in) :: index
         integer (kind=4) :: layoutnumber,i,j,k
         type (XYZlimit_t)    ::  sggSweep

         character(len=BUFSIZE) :: buff
         
         
         if (sggdummy%IsHard) then
            dummy%numHard=dummy%numHard+1
            !
            dummy%nodHard(dummy%numHard)%IsInitialValue=sggdummy%IsInitialValue
            dummy%nodHard(dummy%numHard)%punto%XI = max(sggdummy%punto(index)%XI,sggSweep%XI)
            dummy%nodHard(dummy%numHard)%punto%XE = min(sggdummy%punto(index)%XE,sggSweep%XE)
            dummy%nodHard(dummy%numHard)%punto%YI = max(sggdummy%punto(index)%YI,sggSweep%YI)
            dummy%nodHard(dummy%numHard)%punto%YE = min(sggdummy%punto(index)%YE,sggSweep%YE)
            dummy%nodHard(dummy%numHard)%punto%ZI = max(sggdummy%punto(index)%ZI,sggSweep%ZI)
            dummy%nodHard(dummy%numHard)%punto%ZE = min(sggdummy%punto(index)%ZE,sggSweep%ZE)
            !
            dummy%nodHard(dummy%numHard)%punto%amplitude = amplit
            !Read the time evolution
            dummy%nodHard(dummy%numHard)%deltaevol=sggdummy%fichero%deltaSamples
            if (dummy%nodHard(dummy%numHard)%deltaevol > sgg%dt) then
               write (buff,'(a,e12.2e3)')  'WARNING: '//trim(adjustl(sggdummy%Fichero%Name))// &
               ' undersampled by a factor ',dummy%nodHard(dummy%numHard)%deltaevol/sgg%dt
               call WarnErrReport(buff)
            endif
            dummy%nodHard(dummy%numHard)%numus =  sggdummy%Fichero%NumSamples
            dummy%nodHard(dummy%numHard)%evol  => sggdummy%fichero%Samples
         else
            dummy%numSoft=dummy%numSoft+1
            !
            dummy%nodSoft(dummy%numSoft)%IsInitialValue=sggdummy%IsInitialValue
            dummy%nodSoft(dummy%numSoft)%punto%XI = max(sggdummy%punto(index)%XI,sggSweep%XI)
            dummy%nodSoft(dummy%numSoft)%punto%XE = min(sggdummy%punto(index)%XE,sggSweep%XE)
            dummy%nodSoft(dummy%numSoft)%punto%YI = max(sggdummy%punto(index)%YI,sggSweep%YI)
            dummy%nodSoft(dummy%numSoft)%punto%YE = min(sggdummy%punto(index)%YE,sggSweep%YE)
            dummy%nodSoft(dummy%numSoft)%punto%ZI = max(sggdummy%punto(index)%ZI,sggSweep%ZI)
            dummy%nodSoft(dummy%numSoft)%punto%ZE = min(sggdummy%punto(index)%ZE,sggSweep%ZE)
            !
            dummy%nodSoft(dummy%numSoft)%punto%amplitude = amplit
            !Read the time evolution
            dummy%nodSoft(dummy%numSoft)%deltaevol=sggdummy%fichero%deltaSamples
            if (dummy%nodSoft(dummy%numSoft)%deltaevol > sgg%dt) then
               write (buff,'(a,e12.2e3)')  'WARNING: '//trim(adjustl(sggdummy%Fichero%Name))// &
               ' undersampled by a factor ',dummy%nodSoft(dummy%numSoft)%deltaevol/sgg%dt
               call WarnErrReport(buff)
            endif
            dummy%nodSoft(dummy%numSoft)%numus =  sggdummy%Fichero%NumSamples
            dummy%nodSoft(dummy%numSoft)%evol  => sggdummy%fichero%Samples
         endif

         return

      end subroutine createnodal

   end subroutine InitnodalSources




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the E-field at n
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Evolution function to interpolate from the input file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   REAL (KIND=RKIND) function evolucion(t,dummy)
      REAL (KIND=RKIND) t,deltaevol
      integer (kind=4)  ::  numus
      integer (kind=8)  ::  nprev
      REAL (KIND=RKIND), pointer, dimension ( : )  ::  evol
      type (NodalLocal_t), intent (in) :: dummy

      if (dummy%IsInitialValue) then
        !!!evolucion=1.0_RKIND
        if (int(t/dummy%deltaevol)/=0) then
            print *,'error en initial values. '
            stop
        endif
        evolucion=dummy%evol(0) !should be 1 always
        return
      endif
      
      deltaevol = dummy%deltaevol
      evol => dummy%evol
      numus = dummy%numus

      evolucion=0.0_RKIND

      nprev=int(t/deltaevol)
      !first order interpolation
      if ((nprev+1 > numus).OR.(NPREV+1 <= 0)) then !SI NPREV<0 ES PORQUE SE HA DESBORADO EL ENTERO !BUG MIGEL 130614
         evolucion=0.0_RKIND !se asume que el fichero de entrada contiene una excitacion que se anula despues
      else
         evolucion=(evol(nprev+1)-evol(nprev))/deltaevol*((t)-nprev*deltaevol)+evol(nprev) !interpolacion lineal
      endif
      !second order !no advantages over first order
      !  if (nprev+2 > numus) then
      !      evolucion=0.0_RKIND !se asume que el fichero de entrada contiene una excitacion que se anula despues
      !  else
      !      evolucion=evol(nprev+2) * ( ((t)-nprev    *deltaevol) * ((t)-(nprev+1)*deltaevol) ) /(2.0_RKIND * deltaevol**2.0_RKIND ) - &
      !                evol(nprev+1) * ( ((t)-nprev    *deltaevol) * ((t)-(nprev+2)*deltaevol) ) /(   deltaevol**2.0_RKIND ) + &
      !                evol(nprev  ) * ( ((t)-(nprev+2)*deltaevol) * ((t)-(nprev+1)*deltaevol) ) /(2.0_RKIND * deltaevol**2.0_RKIND )
      !  endif


      return

   end function evolucion

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyNodal(sgg)
      type (SGGFDTDINFO), intent(INOUT)         ::  sgg


      if (Nodal_Ex%NumSoft+Nodal_Ey%NumSoft+Nodal_Ez%NumSoft /= 0) then
         if (associated(Nodal_Ex%nodSoft)) deALLOCATE (Nodal_Ex%nodSoft)
         if (associated(Nodal_Ey%nodSoft)) deALLOCATE (Nodal_Ey%nodSoft)
         if (associated(Nodal_Ez%nodSoft)) deALLOCATE (Nodal_Ez%nodSoft)
      endif
      if (Nodal_Ex%NumHard+Nodal_Ey%NumHard+Nodal_Ez%NumHard /= 0) then
         if (associated(Nodal_Ex%nodHard)) deALLOCATE (Nodal_Ex%nodHard)
         if (associated(Nodal_Ey%nodHard)) deALLOCATE (Nodal_Ey%nodHard)
         if (associated(Nodal_Ez%nodHard)) deALLOCATE (Nodal_Ez%nodHard)
      endif
      if (Nodal_Hx%NumSoft+Nodal_Hy%NumSoft+Nodal_Hz%NumSoft /= 0) then
         if (associated(Nodal_Hx%nodSoft)) deALLOCATE (Nodal_Hx%nodSoft)
         if (associated(Nodal_Hy%nodSoft)) deALLOCATE (Nodal_Hy%nodSoft)
         if (associated(Nodal_Hz%nodSoft)) deALLOCATE (Nodal_Hz%nodSoft)
      endif
      if (Nodal_Hx%NumHard+Nodal_Hy%NumHard+Nodal_Hz%NumHard /= 0) then
         if (associated(Nodal_Hx%nodHard)) deALLOCATE (Nodal_Hx%nodHard)
         if (associated(Nodal_Hy%nodHard)) deALLOCATE (Nodal_Hy%nodHard)
         if (associated(Nodal_Hz%nodHard)) deALLOCATE (Nodal_Hz%nodHard)
      endif


      if (associated(sgg%NodalSource)) deallocate (sgg%NodalSource)
   end subroutine DestroyNodal



   !**************************************************************************************************
   subroutine AdvancenodalE(sgg,sggMiEx, sggMiEy, sggMiEz,NumMedia,timeinstant, b, g2,Idxh,Idyh,Idzh,Ex,Ey,Ez,simu_devia)
      !---------------------------> inputs <----------------------------------------------------------
      type (SGGFDTDINFO), intent(IN)     , target  ::  sgg
      logical, intent(in) :: simu_devia
      integer, intent( IN)  ::  NumMedia, timeinstant
      !!!
      type( bounds_t), intent( IN)  ::  b
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent( IN)  ::  sggMiEx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent( IN)  ::  sggMiEy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent( IN)  ::  sggMiEz
      !--->
      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  g2
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN)  ::  Idxh
      real (kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN)  ::  Idyh
      real (kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN)  ::  Idzh
      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( INOUT)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( INOUT)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( INOUT)  ::  Ez

      !---------------------------> variables locales <-----------------------------------------------
      real (kind = RKIND)  ::  timei,amp
      integer  ::  i, j, k, i_m, j_m, k_m,ii,medio
      !---------------------------> empieza AdvancenodalE <---------------------------------------

      !!!
      !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219 
      !!!timei = (timeinstant +3) * sgg%dt !ORIGINAL sync
      timei = sgg%tiempo(timeinstant) 

      !
      barridonodalhardEx: do ii=1,Nodal_Ex%numHard
         if (Nodal_Ex%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardEx
         endif
         !     
         amp = Nodal_Ex%nodHard(ii)%punto%amplitude
         do k=Nodal_Ex%nodHard(ii)%punto%zi,Nodal_Ex%nodHard(ii)%punto%ze
            k_m = k - b%Ex%ZI
            do j=Nodal_Ex%nodHard(ii)%punto%yi,Nodal_Ex%nodHard(ii)%punto%ye
               j_m = j - b%Ex%YI
               do i=Nodal_Ex%nodHard(ii)%punto%xi,Nodal_Ex%nodHard(ii)%punto%xe
                  i_m = i - b%Ex%XI
                  medio = sggMiEx(i_m,j_m,k_m)
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medio)%Is%PEC) Ex(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Ex%nodHard(ii))
                  else
                        if (.not.sgg%Med(medio)%Is%PEC) Ex(i_m,j_m,k_m) = 0.0 !!!!!!
                  endif
               end do
            End do
         end do
      end do barridonodalhardEx
      !
      barridonodalsoftEx: do ii=1,Nodal_Ex%numSoft
         if (Nodal_Ex%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftEx
         endif
         !     
         amp = Nodal_Ex%nodSoft(ii)%punto%amplitude
         do k=Nodal_Ex%nodSoft(ii)%punto%zi,Nodal_Ex%nodSoft(ii)%punto%ze
            k_m = k - b%Ex%ZI
            do j=Nodal_Ex%nodSoft(ii)%punto%yi,Nodal_Ex%nodSoft(ii)%punto%ye
               j_m = j - b%Ex%YI
               do i=Nodal_Ex%nodSoft(ii)%punto%xi,Nodal_Ex%nodSoft(ii)%punto%xe
                  i_m = i - b%Ex%XI
                  medio = sggMiEx(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medio)%Is%PEC) Ex(i_m,j_m,k_m) = Ex(i_m,j_m,k_m)- G2(medio) * Idyh(j_m) * Idzh(k_m) * amp * evolucion(timei,Nodal_Ex%nodSoft(ii)) 
                  else
                       if (.not.sgg%Med(medio)%Is%PEC)  Ex(i_m,j_m,k_m) = Ex(i_m,j_m,k_m) !!!!!!
                  endif
               end do
            End do
         end do
      end do barridonodalsoftEx
      !
      !
      barridonodalhardEy: do ii=1,Nodal_Ey%numHard
         if (Nodal_Ey%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardEy
         endif
         !
         amp = Nodal_Ey%nodHard(ii)%punto%amplitude
         do k=Nodal_Ey%nodHard(ii)%punto%zi,Nodal_Ey%nodHard(ii)%punto%ze
            k_m = k - b%Ey%ZI
            do j=Nodal_Ey%nodHard(ii)%punto%yi,Nodal_Ey%nodHard(ii)%punto%ye
               j_m = j - b%Ey%YI
               do i=Nodal_Ey%nodHard(ii)%punto%xi,Nodal_Ey%nodHard(ii)%punto%xe
                  i_m = i - b%Ey%XI
                  medio = sggMiEy(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medio)%Is%PEC) Ey(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Ey%nodHard(ii))   
                  else
                        if (.not.sgg%Med(medio)%Is%PEC) Ey(i_m,j_m,k_m) = 0.0 !!!!!!
                  endif
               end do
            End do
         end do
      end do barridonodalhardEy
      !
      barridonodalsoftEy: do ii=1,Nodal_Ey%numSoft
         if (Nodal_Ey%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftEy
         endif
         !     
         amp = Nodal_Ey%nodSoft(ii)%punto%amplitude
         do k=Nodal_Ey%nodSoft(ii)%punto%zi,Nodal_Ey%nodSoft(ii)%punto%ze
            k_m = k - b%Ey%ZI
            do j=Nodal_Ey%nodSoft(ii)%punto%yi,Nodal_Ey%nodSoft(ii)%punto%ye
               j_m = j - b%Ey%YI
               do i=Nodal_Ey%nodSoft(ii)%punto%xi,Nodal_Ey%nodSoft(ii)%punto%xe
                  i_m = i - b%Ey%XI
                  medio = sggMiEy(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medio)%Is%PEC) Ey(i_m,j_m,k_m) = Ey(i_m,j_m,k_m)- G2(medio) * Idxh(i_m) * Idzh(k_m) * amp * evolucion(timei,Nodal_Ey%nodSoft(ii))   
                  else
                       if (.not.sgg%Med(medio)%Is%PEC)  Ey(i_m,j_m,k_m) = Ey(i_m,j_m,k_m) !!!!!!
                  endif
               end do
            End do
         end do
      end do barridonodalsoftEy

      barridonodalhardEz: do ii=1,Nodal_Ez%numHard
         if (Nodal_Ez%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardEz
         endif
         !
         amp = Nodal_Ez%nodHard(ii)%punto%amplitude
         do k=Nodal_Ez%nodHard(ii)%punto%zi,Nodal_Ez%nodHard(ii)%punto%ze
            k_m = k - b%Ez%ZI
            do j=Nodal_Ez%nodHard(ii)%punto%yi,Nodal_Ez%nodHard(ii)%punto%ye
               j_m = j - b%Ez%YI
               do i=Nodal_Ez%nodHard(ii)%punto%xi,Nodal_Ez%nodHard(ii)%punto%xe
                  i_m = i - b%Ez%XI
                  medio = sggMiEz(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medio)%Is%PEC) Ez(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Ez%nodHard(ii))  
                  else
                        if (.not.sgg%Med(medio)%Is%PEC) Ez(i_m,j_m,k_m) = 0.0 !!!!!!
                  endif
               end do
            End do
         end do
      end do barridonodalhardEz
      !
      barridonodalsoftEz: do ii=1,Nodal_Ez%numSoft
         if (Nodal_Ez%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftEz
         endif
         !     
         amp = Nodal_Ez%nodSoft(ii)%punto%amplitude
         do k=Nodal_Ez%nodSoft(ii)%punto%zi,Nodal_Ez%nodSoft(ii)%punto%ze
            k_m = k - b%Ez%ZI
            do j=Nodal_Ez%nodSoft(ii)%punto%yi,Nodal_Ez%nodSoft(ii)%punto%ye
               j_m = j - b%Ez%YI
               do i=Nodal_Ez%nodSoft(ii)%punto%xi,Nodal_Ez%nodSoft(ii)%punto%xe
                  i_m = i - b%Ez%XI
                  medio = sggMiEz(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medio)%Is%PEC) Ez(i_m,j_m,k_m) = Ez(i_m,j_m,k_m)- G2(medio) * Idyh(j_m) * Idxh(i_m) * amp * evolucion(timei,Nodal_Ez%nodSoft(ii))
                  else
                        if (.not.sgg%Med(medio)%Is%PEC) Ez(i_m,j_m,k_m) = Ez(i_m,j_m,k_m) !!!!!!
                  endif
               end do
            End do
         end do
      end do barridonodalsoftEz




      return

   endsubroutine AdvancenodalE
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the H-field at n+0.5_RKIND
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**************************************************************************************************
   subroutine AdvancenodalH(sgg,sggMiHx, sggMiHy, sggMiHz,NumMedia,timeinstant, b,gm2,Idxe,Idye,Idze,Hx,Hy,Hz,simu_devia)
      !---------------------------> inputs <----------------------------------------------------------
      type (SGGFDTDINFO), intent(IN)     , target  ::  sgg
      logical , intent(in) :: simu_devia !ojo untested con simu_devia este tipo de fuentes
      integer, intent( IN)  ::  NumMedia, timeinstant
      !!!
      type( bounds_t), intent( IN)  ::  b
      !--->
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent( IN)  ::  sggMiHx
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent( IN)  ::  sggMiHy
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent( IN)  ::  sggMiHz
      !--->
      real (kind = RKIND), dimension( 0 :  NumMedia), intent( IN)  ::  gm2
      !--->
      real (kind = RKIND), dimension( 0 :  b%dxh%NX-1), intent( IN)  ::  Idxe
      real (kind = RKIND), dimension( 0 :  b%dyh%NY-1), intent( IN)  ::  Idye
      real (kind = RKIND), dimension( 0 :  b%dzh%NZ-1), intent( IN)  ::  Idze

      !---------------------------> inputs/outputs <--------------------------------------------------
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( INOUT)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( INOUT)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( INOUT)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      real (kind = RKIND)  ::  timei,amp
      integer (kind=4)  ::  i, j, k, i_m, j_m, k_m,ii,medio
      real (kind = RKIND)  ::  GM2_1
      !!!
      if (simu_devia) then
          print *,'Devia H nodal/field sources untested. Aborting'
          stop
      endif
      GM2_1=GM2(1)
      !---------------------------> empieza AdvancenodalH <---------------------------------------
      
      timei = sgg%tiempo(timeinstant) + 0.5_RKIND  * sgg%dt
      !!!! deprecado en pscale y el+3 de la sincronia con ORIGINAL se jode para siempre 110219 
      !!! timei = ( timeinstant + 0.5_RKIND  +3.0_RKIND) * sgg%dt  !ORIGINAL sync


      barridonodalhardHx: do ii=1,Nodal_Hx%numHard
         if (Nodal_Hx%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardHx
         endif
         !     
         amp = Nodal_Hx%nodHard(ii)%punto%amplitude
         do k=Nodal_Hx%nodHard(ii)%punto%zi,Nodal_Hx%nodHard(ii)%punto%ze
            k_m = k - b%Hx%ZI
            do j=Nodal_Hx%nodHard(ii)%punto%yi,Nodal_Hx%nodHard(ii)%punto%ye
               j_m = j - b%Hx%YI
               do i=Nodal_Hx%nodHard(ii)%punto%xi,Nodal_Hx%nodHard(ii)%punto%xe
                  i_m = i - b%Hx%XI
                  medio = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medio)%Is%PMC) Hx(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Hx%nodHard(ii))
               end do
            End do
         end do
      end do barridonodalhardHx
      !
      barridonodalsoftHx: do ii=1,Nodal_Hx%numSoft
         if (Nodal_Hx%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftHx
         endif
         !     
         amp = Nodal_Hx%nodSoft(ii)%punto%amplitude
         do k=Nodal_Hx%nodSoft(ii)%punto%zi,Nodal_Hx%nodSoft(ii)%punto%ze
            k_m = k - b%Hx%ZI
            do j=Nodal_Hx%nodSoft(ii)%punto%yi,Nodal_Hx%nodSoft(ii)%punto%ye
               j_m = j - b%Hx%YI
               do i=Nodal_Hx%nodSoft(ii)%punto%xi,Nodal_Hx%nodSoft(ii)%punto%xe
                  i_m = i - b%Hx%XI
                  medio = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medio)%Is%PMC) Hx(i_m,j_m,k_m) = Hx(i_m,j_m,k_m)- Gm2(medio) * Idye(j_m) * Idze(k_m) * amp * evolucion(timei,Nodal_Hx%nodSoft(ii))
               end do
            End do
         end do
      end do barridonodalsoftHx
      !
      !
      barridonodalhardHy: do ii=1,Nodal_Hy%numHard
         if (Nodal_Hy%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardHy
         endif
         !
         amp = Nodal_Hy%nodHard(ii)%punto%amplitude
         do k=Nodal_Hy%nodHard(ii)%punto%zi,Nodal_Hy%nodHard(ii)%punto%ze
            k_m = k - b%Hy%ZI
            do j=Nodal_Hy%nodHard(ii)%punto%yi,Nodal_Hy%nodHard(ii)%punto%ye
               j_m = j - b%Hy%YI
               do i=Nodal_Hy%nodHard(ii)%punto%xi,Nodal_Hy%nodHard(ii)%punto%xe
                  i_m = i - b%Hy%XI
                  medio = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medio)%Is%PMC) Hy(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Hy%nodHard(ii))
               end do
            End do
         end do
      end do barridonodalhardHy
      !
      barridonodalsoftHy: do ii=1,Nodal_Hy%numSoft
         if (Nodal_Hy%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftHy
         endif
         !     
         amp = Nodal_Hy%nodSoft(ii)%punto%amplitude
         do k=Nodal_Hy%nodSoft(ii)%punto%zi,Nodal_Hy%nodSoft(ii)%punto%ze
            k_m = k - b%Hy%ZI
            do j=Nodal_Hy%nodSoft(ii)%punto%yi,Nodal_Hy%nodSoft(ii)%punto%ye
               j_m = j - b%Hy%YI
               do i=Nodal_Hy%nodSoft(ii)%punto%xi,Nodal_Hy%nodSoft(ii)%punto%xe
                  i_m = i - b%Hy%XI
                  medio = sggMiHy(i_m,j_m,k_m)
                  if (.not.sgg%Med(medio)%Is%PMC) Hy(i_m,j_m,k_m) = Hy(i_m,j_m,k_m)- Gm2(medio) * Idxe(i_m) * Idze(k_m) * amp * evolucion(timei,Nodal_Hy%nodSoft(ii))
               end do
            End do
         end do
      end do barridonodalsoftHy

      barridonodalhardHz: do ii=1,Nodal_Hz%numHard
         if (Nodal_Hz%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardHz
         endif
         !
         amp = Nodal_Hz%nodHard(ii)%punto%amplitude
         do k=Nodal_Hz%nodHard(ii)%punto%zi,Nodal_Hz%nodHard(ii)%punto%ze
            k_m = k - b%Hz%ZI
            do j=Nodal_Hz%nodHard(ii)%punto%yi,Nodal_Hz%nodHard(ii)%punto%ye
               j_m = j - b%Hz%YI
               do i=Nodal_Hz%nodHard(ii)%punto%xi,Nodal_Hz%nodHard(ii)%punto%xe
                  i_m = i - b%Hz%XI
                  medio = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medio)%Is%PMC) Hz(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Hz%nodHard(ii))
               end do
            End do
         end do
      end do barridonodalhardHz
      !
      barridonodalsoftHz: do ii=1,Nodal_Hz%numSoft
         if (Nodal_Hz%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftHz
         endif
         !     
         amp = Nodal_Hz%nodSoft(ii)%punto%amplitude
         do k=Nodal_Hz%nodSoft(ii)%punto%zi,Nodal_Hz%nodSoft(ii)%punto%ze
            k_m = k - b%Hz%ZI
            do j=Nodal_Hz%nodSoft(ii)%punto%yi,Nodal_Hz%nodSoft(ii)%punto%ye
               j_m = j - b%Hz%YI
               do i=Nodal_Hz%nodSoft(ii)%punto%xi,Nodal_Hz%nodSoft(ii)%punto%xe
                  i_m = i - b%Hz%XI
                  medio = sggMiHz(i_m,j_m,k_m)
                  if (.not.sgg%Med(medio)%Is%PMC) Hz(i_m,j_m,k_m) = Hz(i_m,j_m,k_m)- Gm2(medio) * Idye(j_m) * Idxe(i_m) * amp * evolucion(timei,Nodal_Hz%nodSoft(ii))
               end do
            End do
         end do
      end do barridonodalsoftHz




      return
   endsubroutine AdvancenodalH


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Function to publish the private output data (used in postprocess)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getnodal(rNodal_Ex,rNodal_Ey,rNodal_Ez,rNodal_Hx,rNodal_Hy,rNodal_Hz)

      type (nodsou), pointer :: rNodal_Ex ,rNodal_Ey ,rNodal_Ez
      type (nodsou), pointer :: rNodal_Hx ,rNodal_Hy ,rNodal_Hz

      rNodal_Ex  => Nodal_Ex
      rNodal_Ey  => Nodal_Ey
      rNodal_Ez  => Nodal_Ez
      rNodal_Hx  => Nodal_Hx
      rNodal_Hy  => Nodal_Hy
      rNodal_Hz  => Nodal_Hz


      return
   end subroutine

   
   subroutine InitHopf(sgg,NumNodalSources,sggNodalSource,sggSweep,ficherohopf)
      type (SGGFDTDINFO), intent(IN)         ::  sgg
      integer, intent (in) :: NumNodalSources
      type (NodalSource_t), dimension(1:NumNodalSources),intent(in)           ::  sggNodalSource
      type (XYZlimit_t), dimension (1:6)    ::  sggSweep
      character (LEN=BUFSIZE) :: ficherohopf
      return
   end subroutine InitHopf

END MODULE nodalsources
 