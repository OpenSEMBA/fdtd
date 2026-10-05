module nodalsources_m

   use FDETYPES_m
   use Report_m

   implicit none
   private

   type xyzlimit_singlescaled_t
      integer(kind=4) :: XI,XE,YI,YE,ZI,ZE
      real(kind=RKIND) :: amplitude
   end type


   type :: NodalLocal_t
      real(kind=RKIND), pointer, dimension(:) :: evol
      real(kind=RKIND) :: deltaevol
      integer(kind=4) :: numus
      type(xyzlimit_singlescaled_t) :: gridPoint
      logical :: IsInitialValue
   end type NodalLocal_t


   type, public :: nodsou_t
      integer(kind=4) :: NumHard = 0 , NumSoft = 0
      type(NodalLocal_t), pointer, dimension(:) :: nodHard,nodSoft
   end type

   !!!!!local variables

   type(nodsou_t), save, target :: Nodal_Ex,Nodal_Ey,Nodal_Ez
   type(nodsou_t), save, target :: Nodal_Hx,Nodal_Hy,Nodal_Hz

   public :: initNodalSources,AdvanceNodalE,AdvanceNodalH,DestroyNodal,getnodal




contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes Nodal Source data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitnodalSources(sgg,layoutnumber,NumNodalSources,sggNodalSource,sggSweep,ThereAreNodalE,ThereAreNodalH)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      !!!
      integer, intent(in) :: NumNodalSources
      type(NodalSource_t), dimension(1:NumNodalSources),intent(in) :: sggNodalSource

      integer(kind=4):: layoutnumber,j,i
      logical, intent(out) :: ThereArenodalE,ThereArenodalH
      integer :: numNodalSoft_Ex,numNodalSoft_Ey,numNodalSoft_Ez, &
      numNodalSoft_Hx,numNodalSoft_Hy,numNodalSoft_Hz, &
      numNodalHard_Ex,numNodalHard_Ey,numNodalHard_Ez, &
      numNodalHard_Hx,numNodalHard_Hy,numNodalHard_Hz

      real(kind=rkind) :: amplit
      type(XYZlimit_t), dimension(1:6) :: sggSweep

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
            do i=1,sggNodalSource(j)%numPoints
               if (sggNodalSource(j)%gridPoint(i)%xc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Ex = numNodalHard_Ex  +1
                  else
                     numNodalSoft_Ex = numNodalSoft_Ex  +1
                  end if
               end if
               if (sggNodalSource(j)%gridPoint(i)%yc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Ey = numNodalHard_Ey  +1
                  else
                     numNodalSoft_Ey = numNodalSoft_Ey  +1
                  end if
               end if
               if (sggNodalSource(j)%gridPoint(i)%zc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Ez = numNodalHard_Ez  +1
                  else
                     numNodalSoft_Ez = numNodalSoft_Ez  +1
                  end if
               end if
            end do
         else
            do i=1,sggNodalSource(j)%numPoints
               if (sggNodalSource(j)%gridPoint(i)%xc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Hx = numNodalHard_Hx  +1
                  else
                     numNodalSoft_Hx = numNodalSoft_Hx  +1
                  end if
               end if
               if (sggNodalSource(j)%gridPoint(i)%yc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Hy = numNodalHard_Hy  +1
                  else
                     numNodalSoft_Hy = numNodalSoft_Hy  +1
                  end if
               end if
               if (sggNodalSource(j)%gridPoint(i)%zc /= 0.0_RKIND) then
                  if (sggNodalSource(j)%IsHard) then
                     numNodalHard_Hz = numNodalHard_Hz  +1
                  else
                     numNodalSoft_Hz = numNodalSoft_Hz  +1
                  end if
               end if
            end do
         end if
      end do


      if (NumNodalSoft_Ex+NumNodalSoft_Ey+NumNodalSoft_Ez /= 0) then
         ThereArenodalE =.true.
        allocate(Nodal_Ex%nodSoft(1:numNodalSoft_Ex), &
         Nodal_Ey%nodSoft(1:numNodalSoft_Ey), &
         Nodal_Ez%nodSoft(1:numNodalSoft_Ez))
      end if
      if (NumNodalHard_Ex+NumNodalHard_Ey+NumNodalHard_Ez /= 0) then
         ThereArenodalE =.true.
        allocate(Nodal_Ex%nodHard(1:numNodalHard_Ex), &
         Nodal_Ey%nodHard(1:numNodalHard_Ey), &
         Nodal_Ez%nodHard(1:numNodalHard_Ez))
      end if
      if (NumNodalSoft_Hx+NumNodalSoft_Hy+NumNodalSoft_Hz /= 0) then
         ThereArenodalH =.true.
        allocate(Nodal_Hx%nodSoft(1:numNodalSoft_Hx), &
         Nodal_Hy%nodSoft(1:numNodalSoft_Hy), &
         Nodal_Hz%nodSoft(1:numNodalSoft_Hz))
      end if
      if (NumNodalHard_Hx+NumNodalHard_Hy+NumNodalHard_Hz /= 0) then
         ThereArenodalH =.true.
        allocate(Nodal_Hx%nodHard(1:numNodalHard_Hx), &
         Nodal_Hy%nodHard(1:numNodalHard_Hy), &
         Nodal_Hz%nodHard(1:numNodalHard_Hz))
      end if


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
            do i=1,sggNodalSource(j)%numPoints
               amplit = sggNodalSource(J)%gridPoint(i)%xc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Ex,sggNodalSource(J),sggSweep(IEX),i,amplit)
               end if
               amplit = sggNodalSource(j)%gridPoint(i)%yc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Ey,sggNodalSource(J),sggSweep(IEY),i,amplit)
               end if
               amplit = sggNodalSource(j)%gridPoint(i)%zc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Ez,sggNodalSource(J),sggSweep(IEZ),i,amplit)
               end if
            end do
         else !it is magnetic
            do i=1,sggNodalSource(j)%numPoints
               amplit = sggNodalSource(J)%gridPoint(i)%xc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Hx,sggNodalSource(J),sggSweep(IHX),i,amplit)
               end if
               amplit = sggNodalSource(j)%gridPoint(i)%yc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Hy,sggNodalSource(J),sggSweep(IHY),i,amplit)
               end if
               amplit = sggNodalSource(j)%gridPoint(i)%zc
               if (amplit /= 0.0_RKIND) then
                  call CreateNodal(layoutnumber,Nodal_Hz,sggNodalSource(J),sggSweep(IHZ),i,amplit)
               end if
            end do
         end if
      end do

      !print *,'tras nodal sources ',Nodal_Ex%numHard

      return


   contains

      subroutine createnodal(layoutnumber,dummy,sggdummy,sggSweep,elementIndex,amplit)

         type(nodsou_t), intent (inout) :: dummy
         type(NodalSource_t), intent(in), target :: sggdummy
         real(kind=rkind), intent(in) :: amplit
         integer(kind=4), intent(in) :: elementIndex
         integer(kind=4) :: layoutnumber,i,j,k
         type(XYZlimit_t) :: sggSweep

         character(len=BUFSIZE) :: buff
         
         
         if (sggdummy%IsHard) then
            dummy%numHard=dummy%numHard+1
            !
            dummy%nodHard(dummy%numHard)%IsInitialValue=sggdummy%IsInitialValue
            dummy%nodHard(dummy%numHard)%gridPoint%XI = max(sggdummy%gridPoint(elementIndex)%XI,sggSweep%XI)
            dummy%nodHard(dummy%numHard)%gridPoint%XE = min(sggdummy%gridPoint(elementIndex)%XE,sggSweep%XE)
            dummy%nodHard(dummy%numHard)%gridPoint%YI = max(sggdummy%gridPoint(elementIndex)%YI,sggSweep%YI)
            dummy%nodHard(dummy%numHard)%gridPoint%YE = min(sggdummy%gridPoint(elementIndex)%YE,sggSweep%YE)
            dummy%nodHard(dummy%numHard)%gridPoint%ZI = max(sggdummy%gridPoint(elementIndex)%ZI,sggSweep%ZI)
            dummy%nodHard(dummy%numHard)%gridPoint%ZE = min(sggdummy%gridPoint(elementIndex)%ZE,sggSweep%ZE)
            !
            dummy%nodHard(dummy%numHard)%gridPoint%amplitude = amplit
            !Read the time evolution
            dummy%nodHard(dummy%numHard)%deltaevol=sggdummy%sourceFile%deltaSamples
            if (dummy%nodHard(dummy%numHard)%deltaevol > sgg%dt) then
               write (buff,'(a,e12.2e3)')  'WARNING: '//trim(adjustl(sggdummy%sourceFile%Name))// &
               ' undersampled by a factor ',dummy%nodHard(dummy%numHard)%deltaevol/sgg%dt
               call WarnErrReport(buff)
            end if
            dummy%nodHard(dummy%numHard)%numus =  sggdummy%sourceFile%NumSamples
            dummy%nodHard(dummy%numHard)%evol  => sggdummy%sourceFile%Samples
         else
            dummy%numSoft=dummy%numSoft+1
            !
            dummy%nodSoft(dummy%numSoft)%IsInitialValue=sggdummy%IsInitialValue
            dummy%nodSoft(dummy%numSoft)%gridPoint%XI = max(sggdummy%gridPoint(elementIndex)%XI,sggSweep%XI)
            dummy%nodSoft(dummy%numSoft)%gridPoint%XE = min(sggdummy%gridPoint(elementIndex)%XE,sggSweep%XE)
            dummy%nodSoft(dummy%numSoft)%gridPoint%YI = max(sggdummy%gridPoint(elementIndex)%YI,sggSweep%YI)
            dummy%nodSoft(dummy%numSoft)%gridPoint%YE = min(sggdummy%gridPoint(elementIndex)%YE,sggSweep%YE)
            dummy%nodSoft(dummy%numSoft)%gridPoint%ZI = max(sggdummy%gridPoint(elementIndex)%ZI,sggSweep%ZI)
            dummy%nodSoft(dummy%numSoft)%gridPoint%ZE = min(sggdummy%gridPoint(elementIndex)%ZE,sggSweep%ZE)
            !
            dummy%nodSoft(dummy%numSoft)%gridPoint%amplitude = amplit
            !Read the time evolution
            dummy%nodSoft(dummy%numSoft)%deltaevol=sggdummy%sourceFile%deltaSamples
            if (dummy%nodSoft(dummy%numSoft)%deltaevol > sgg%dt) then
               write (buff,'(a,e12.2e3)')  'WARNING: '//trim(adjustl(sggdummy%sourceFile%Name))// &
               ' undersampled by a factor ',dummy%nodSoft(dummy%numSoft)%deltaevol/sgg%dt
               call WarnErrReport(buff)
            end if
            dummy%nodSoft(dummy%numSoft)%numus =  sggdummy%sourceFile%NumSamples
            dummy%nodSoft(dummy%numSoft)%evol  => sggdummy%sourceFile%Samples
         end if

         return

      end subroutine createnodal

   end subroutine InitnodalSources




   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the E-field at n
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Evolution function to interpolate from the input file
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   real(kind=RKIND) function evolucion(t,dummy)
      real(kind=RKIND) t,deltaevol
      integer(kind=4) :: numus
      integer(kind=8) :: nprev
      real(kind=RKIND), pointer, dimension(:) :: evol
      type(NodalLocal_t), intent(in) :: dummy

      if (dummy%IsInitialValue) then
        !!!evolucion=1.0_RKIND
        if (int(t/dummy%deltaevol)/=0) then
            print *,'error en initial values. '
            stop
        end if
        evolucion=dummy%evol(0) !should be 1 always
        return
      end if
      
      deltaevol = dummy%deltaevol
      evol => dummy%evol
      numus = dummy%numus

      evolucion=0.0_RKIND

      nprev=int(t/deltaevol)
      !first order interpolation
      if ((nprev+1 > numus).OR.(NPREV+1 <= 0)) then !IF NPREV<0 IT IS BECAUSE THE INTEGER HAS OVERFLOWED !BUG MIGEL 130614
         evolucion=0.0_RKIND !it is assumed that the input file contains an excitation that vanishes afterwards
      else
         evolucion=(evol(nprev+1)-evol(nprev))/deltaevol*((t)-nprev*deltaevol)+evol(nprev) !linear interpolation
      end if
      !second order !no advantages over first order
      !  if (nprev+2 > numus) then
      !      evolucion=0.0_RKIND !it is assumed that the input file contains an excitation that vanishes afterwards
      !  else
      !      evolucion=evol(nprev+2) * ( ((t)-nprev    *deltaevol) * ((t)-(nprev+1)*deltaevol) ) /(2.0_RKIND * deltaevol**2.0_RKIND ) - &
      !                evol(nprev+1) * ( ((t)-nprev    *deltaevol) * ((t)-(nprev+2)*deltaevol) ) /(   deltaevol**2.0_RKIND ) + &
      !                evol(nprev  ) * ( ((t)-(nprev+2)*deltaevol) * ((t)-(nprev+1)*deltaevol) ) /(2.0_RKIND * deltaevol**2.0_RKIND )
      !  end if


      return

   end function evolucion

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!  Free-up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine DestroyNodal(sgg)
      type(SGGFDTDINFO_t), intent(inout) :: sgg


      if (Nodal_Ex%NumSoft+Nodal_Ey%NumSoft+Nodal_Ez%NumSoft /= 0) then
         if (associated(Nodal_Ex%nodSoft)) deallocate(Nodal_Ex%nodSoft)
         if (associated(Nodal_Ey%nodSoft)) deallocate(Nodal_Ey%nodSoft)
         if (associated(Nodal_Ez%nodSoft)) deallocate(Nodal_Ez%nodSoft)
      end if
      if (Nodal_Ex%NumHard+Nodal_Ey%NumHard+Nodal_Ez%NumHard /= 0) then
         if (associated(Nodal_Ex%nodHard)) deallocate(Nodal_Ex%nodHard)
         if (associated(Nodal_Ey%nodHard)) deallocate(Nodal_Ey%nodHard)
         if (associated(Nodal_Ez%nodHard)) deallocate(Nodal_Ez%nodHard)
      end if
      if (Nodal_Hx%NumSoft+Nodal_Hy%NumSoft+Nodal_Hz%NumSoft /= 0) then
         if (associated(Nodal_Hx%nodSoft)) deallocate(Nodal_Hx%nodSoft)
         if (associated(Nodal_Hy%nodSoft)) deallocate(Nodal_Hy%nodSoft)
         if (associated(Nodal_Hz%nodSoft)) deallocate(Nodal_Hz%nodSoft)
      end if
      if (Nodal_Hx%NumHard+Nodal_Hy%NumHard+Nodal_Hz%NumHard /= 0) then
         if (associated(Nodal_Hx%nodHard)) deallocate(Nodal_Hx%nodHard)
         if (associated(Nodal_Hy%nodHard)) deallocate(Nodal_Hy%nodHard)
         if (associated(Nodal_Hz%nodHard)) deallocate(Nodal_Hz%nodHard)
      end if


      if (associated(sgg%NodalSource)) deallocate(sgg%NodalSource)
   end subroutine DestroyNodal



   !**************************************************************************************************
   subroutine AdvancenodalE(sgg,sggMiEx, sggMiEy, sggMiEz,NumMedia,timeinstant, b, g2,Idxh,Idyh,Idzh,Ex,Ey,Ez,simu_devia)
      !---------------------------> inputs <----------------------------------------------------------
      type(SGGFDTDINFO_t), intent(in)     , target  :: sgg
      logical, intent(in) :: simu_devia
      integer, intent(in) :: NumMedia, timeinstant
      !!!
      type(bounds_t), intent(in) :: b
      !--->
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiEx%NX-1, 0 :  b%sggMiEx%NY-1, 0 :  b%sggMiEx%NZ-1), intent(in) :: sggMiEx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiEy%NX-1, 0 :  b%sggMiEy%NY-1, 0 :  b%sggMiEy%NZ-1), intent(in) :: sggMiEy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiEz%NX-1, 0 :  b%sggMiEz%NY-1, 0 :  b%sggMiEz%NZ-1), intent(in) :: sggMiEz
      !--->
      real(kind = RKIND), dimension(0 :  NumMedia), intent(in) :: g2
      !--->
      real(kind = RKIND), dimension(0 :  b%dxh%NX-1), intent(in) :: Idxh
      real(kind = RKIND), dimension(0 :  b%dyh%NY-1), intent(in) :: Idyh
      real(kind = RKIND), dimension(0 :  b%dzh%NZ-1), intent(in) :: Idzh
      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent(inout) :: Ex
      real(kind = RKIND), dimension(0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent(inout) :: Ey
      real(kind = RKIND), dimension(0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent(inout) :: Ez

      !---------------------------> local variables <-----------------------------------------------
      real(kind = RKIND) :: timei,amp
      integer  :: i, j, k, i_m, j_m, k_m,ii,medium
      !---------------------------> starts AdvancenodalE <---------------------------------------

      !!!
      !!!! deprecated in pscale and the +3 of the synchronization with ORIGINAL is broken forever 110219 
      !!!timei = (timeinstant +3) * sgg%dt !ORIGINAL sync
      timei = sgg%time(timeinstant) 

      !
      barridonodalhardEx: do ii=1,Nodal_Ex%numHard
         if (Nodal_Ex%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardEx
         end if
         !     
         amp = Nodal_Ex%nodHard(ii)%gridPoint%amplitude
         do k=Nodal_Ex%nodHard(ii)%gridPoint%zi,Nodal_Ex%nodHard(ii)%gridPoint%ze
            k_m = k - b%Ex%ZI
            do j=Nodal_Ex%nodHard(ii)%gridPoint%yi,Nodal_Ex%nodHard(ii)%gridPoint%ye
               j_m = j - b%Ex%YI
               do i=Nodal_Ex%nodHard(ii)%gridPoint%xi,Nodal_Ex%nodHard(ii)%gridPoint%xe
                  i_m = i - b%Ex%XI
                  medium = sggMiEx(i_m,j_m,k_m)
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medium)%Is%PEC) Ex(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Ex%nodHard(ii))
                  else
                        if (.not.sgg%Med(medium)%Is%PEC) Ex(i_m,j_m,k_m) = 0.0 !!!!!!
                  end if
               end do
            end do
         end do
      end do barridonodalhardEx
      !
      barridonodalsoftEx: do ii=1,Nodal_Ex%numSoft
         if (Nodal_Ex%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftEx
         end if
         !     
         amp = Nodal_Ex%nodSoft(ii)%gridPoint%amplitude
         do k=Nodal_Ex%nodSoft(ii)%gridPoint%zi,Nodal_Ex%nodSoft(ii)%gridPoint%ze
            k_m = k - b%Ex%ZI
            do j=Nodal_Ex%nodSoft(ii)%gridPoint%yi,Nodal_Ex%nodSoft(ii)%gridPoint%ye
               j_m = j - b%Ex%YI
               do i=Nodal_Ex%nodSoft(ii)%gridPoint%xi,Nodal_Ex%nodSoft(ii)%gridPoint%xe
                  i_m = i - b%Ex%XI
                  medium = sggMiEx(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medium)%Is%PEC) Ex(i_m,j_m,k_m) = Ex(i_m,j_m,k_m)- G2(medium) * Idyh(j_m) * Idzh(k_m) * amp * evolucion(timei,Nodal_Ex%nodSoft(ii)) 
                  else
                       if (.not.sgg%Med(medium)%Is%PEC)  Ex(i_m,j_m,k_m) = Ex(i_m,j_m,k_m) !!!!!!
                  end if
               end do
            end do
         end do
      end do barridonodalsoftEx
      !
      !
      barridonodalhardEy: do ii=1,Nodal_Ey%numHard
         if (Nodal_Ey%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardEy
         end if
         !
         amp = Nodal_Ey%nodHard(ii)%gridPoint%amplitude
         do k=Nodal_Ey%nodHard(ii)%gridPoint%zi,Nodal_Ey%nodHard(ii)%gridPoint%ze
            k_m = k - b%Ey%ZI
            do j=Nodal_Ey%nodHard(ii)%gridPoint%yi,Nodal_Ey%nodHard(ii)%gridPoint%ye
               j_m = j - b%Ey%YI
               do i=Nodal_Ey%nodHard(ii)%gridPoint%xi,Nodal_Ey%nodHard(ii)%gridPoint%xe
                  i_m = i - b%Ey%XI
                  medium = sggMiEy(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medium)%Is%PEC) Ey(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Ey%nodHard(ii))   
                  else
                        if (.not.sgg%Med(medium)%Is%PEC) Ey(i_m,j_m,k_m) = 0.0 !!!!!!
                  end if
               end do
            end do
         end do
      end do barridonodalhardEy
      !
      barridonodalsoftEy: do ii=1,Nodal_Ey%numSoft
         if (Nodal_Ey%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftEy
         end if
         !     
         amp = Nodal_Ey%nodSoft(ii)%gridPoint%amplitude
         do k=Nodal_Ey%nodSoft(ii)%gridPoint%zi,Nodal_Ey%nodSoft(ii)%gridPoint%ze
            k_m = k - b%Ey%ZI
            do j=Nodal_Ey%nodSoft(ii)%gridPoint%yi,Nodal_Ey%nodSoft(ii)%gridPoint%ye
               j_m = j - b%Ey%YI
               do i=Nodal_Ey%nodSoft(ii)%gridPoint%xi,Nodal_Ey%nodSoft(ii)%gridPoint%xe
                  i_m = i - b%Ey%XI
                  medium = sggMiEy(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medium)%Is%PEC) Ey(i_m,j_m,k_m) = Ey(i_m,j_m,k_m)- G2(medium) * Idxh(i_m) * Idzh(k_m) * amp * evolucion(timei,Nodal_Ey%nodSoft(ii))   
                  else
                       if (.not.sgg%Med(medium)%Is%PEC)  Ey(i_m,j_m,k_m) = Ey(i_m,j_m,k_m) !!!!!!
                  end if
               end do
            end do
         end do
      end do barridonodalsoftEy

      barridonodalhardEz: do ii=1,Nodal_Ez%numHard
         if (Nodal_Ez%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardEz
         end if
         !
         amp = Nodal_Ez%nodHard(ii)%gridPoint%amplitude
         do k=Nodal_Ez%nodHard(ii)%gridPoint%zi,Nodal_Ez%nodHard(ii)%gridPoint%ze
            k_m = k - b%Ez%ZI
            do j=Nodal_Ez%nodHard(ii)%gridPoint%yi,Nodal_Ez%nodHard(ii)%gridPoint%ye
               j_m = j - b%Ez%YI
               do i=Nodal_Ez%nodHard(ii)%gridPoint%xi,Nodal_Ez%nodHard(ii)%gridPoint%xe
                  i_m = i - b%Ez%XI
                  medium = sggMiEz(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medium)%Is%PEC) Ez(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Ez%nodHard(ii))  
                  else
                        if (.not.sgg%Med(medium)%Is%PEC) Ez(i_m,j_m,k_m) = 0.0 !!!!!!
                  end if
               end do
            end do
         end do
      end do barridonodalhardEz
      !
      barridonodalsoftEz: do ii=1,Nodal_Ez%numSoft
         if (Nodal_Ez%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftEz
         end if
         !     
         amp = Nodal_Ez%nodSoft(ii)%gridPoint%amplitude
         do k=Nodal_Ez%nodSoft(ii)%gridPoint%zi,Nodal_Ez%nodSoft(ii)%gridPoint%ze
            k_m = k - b%Ez%ZI
            do j=Nodal_Ez%nodSoft(ii)%gridPoint%yi,Nodal_Ez%nodSoft(ii)%gridPoint%ye
               j_m = j - b%Ez%YI
               do i=Nodal_Ez%nodSoft(ii)%gridPoint%xi,Nodal_Ez%nodSoft(ii)%gridPoint%xe
                  i_m = i - b%Ez%XI
                  medium = sggMiEz(i_m,j_m,k_m)
                  
                  if (.not.simu_devia)   then !bug 280323 mdrc
                        if (.not.sgg%Med(medium)%Is%PEC) Ez(i_m,j_m,k_m) = Ez(i_m,j_m,k_m)- G2(medium) * Idyh(j_m) * Idxh(i_m) * amp * evolucion(timei,Nodal_Ez%nodSoft(ii))
                  else
                        if (.not.sgg%Med(medium)%Is%PEC) Ez(i_m,j_m,k_m) = Ez(i_m,j_m,k_m) !!!!!!
                  end if
               end do
            end do
         end do
      end do barridonodalsoftEz




      return

   end subroutine AdvancenodalE
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Feed the currents to illuminate the H-field at n+0.5_RKIND
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !**************************************************************************************************
   subroutine AdvancenodalH(sgg,sggMiHx, sggMiHy, sggMiHz,NumMedia,timeinstant, b,gm2,Idxe,Idye,Idze,Hx,Hy,Hz,simu_devia)
      !---------------------------> inputs <----------------------------------------------------------
      type(SGGFDTDINFO_t), intent(in)     , target  :: sgg
      logical , intent(in) :: simu_devia !note untested with simu_devia this type of sources
      integer, intent(in) :: NumMedia, timeinstant
      !!!
      type(bounds_t), intent(in) :: b
      !--->
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiHx%NX-1, 0 :  b%sggMiHx%NY-1, 0 :  b%sggMiHx%NZ-1), intent(in) :: sggMiHx
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiHy%NX-1, 0 :  b%sggMiHy%NY-1, 0 :  b%sggMiHy%NZ-1), intent(in) :: sggMiHy
      integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension(0 :  b%sggMiHz%NX-1, 0 :  b%sggMiHz%NY-1, 0 :  b%sggMiHz%NZ-1), intent(in) :: sggMiHz
      !--->
      real(kind = RKIND), dimension(0 :  NumMedia), intent(in) :: gm2
      !--->
      real(kind = RKIND), dimension(0 :  b%dxh%NX-1), intent(in) :: Idxe
      real(kind = RKIND), dimension(0 :  b%dyh%NY-1), intent(in) :: Idye
      real(kind = RKIND), dimension(0 :  b%dzh%NZ-1), intent(in) :: Idze

      !---------------------------> inputs/outputs <--------------------------------------------------
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(inout) :: Hx
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(inout) :: Hy
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(inout) :: Hz
      !---------------------------> local variables <-----------------------------------------------
      real(kind = RKIND) :: timei,amp
      integer(kind=4) :: i, j, k, i_m, j_m, k_m,ii,medium
      real(kind = RKIND) :: GM2_1
      !!!
      if (simu_devia) then
          print *,'Devia H nodal/field sources untested. Aborting'
          stop
      end if
      GM2_1=GM2(1)
      !---------------------------> starts AdvancenodalH <---------------------------------------
      
      timei = sgg%time(timeinstant) + 0.5_RKIND  * sgg%dt
      !!!! deprecated in pscale and the +3 of the synchronization with ORIGINAL is broken forever 110219 
      !!! timei = ( timeinstant + 0.5_RKIND  +3.0_RKIND) * sgg%dt  !ORIGINAL sync


      barridonodalhardHx: do ii=1,Nodal_Hx%numHard
         if (Nodal_Hx%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardHx
         end if
         !     
         amp = Nodal_Hx%nodHard(ii)%gridPoint%amplitude
         do k=Nodal_Hx%nodHard(ii)%gridPoint%zi,Nodal_Hx%nodHard(ii)%gridPoint%ze
            k_m = k - b%Hx%ZI
            do j=Nodal_Hx%nodHard(ii)%gridPoint%yi,Nodal_Hx%nodHard(ii)%gridPoint%ye
               j_m = j - b%Hx%YI
               do i=Nodal_Hx%nodHard(ii)%gridPoint%xi,Nodal_Hx%nodHard(ii)%gridPoint%xe
                  i_m = i - b%Hx%XI
                  medium = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medium)%Is%PMC) Hx(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Hx%nodHard(ii))
               end do
            end do
         end do
      end do barridonodalhardHx
      !
      barridonodalsoftHx: do ii=1,Nodal_Hx%numSoft
         if (Nodal_Hx%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftHx
         end if
         !     
         amp = Nodal_Hx%nodSoft(ii)%gridPoint%amplitude
         do k=Nodal_Hx%nodSoft(ii)%gridPoint%zi,Nodal_Hx%nodSoft(ii)%gridPoint%ze
            k_m = k - b%Hx%ZI
            do j=Nodal_Hx%nodSoft(ii)%gridPoint%yi,Nodal_Hx%nodSoft(ii)%gridPoint%ye
               j_m = j - b%Hx%YI
               do i=Nodal_Hx%nodSoft(ii)%gridPoint%xi,Nodal_Hx%nodSoft(ii)%gridPoint%xe
                  i_m = i - b%Hx%XI
                  medium = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medium)%Is%PMC) Hx(i_m,j_m,k_m) = Hx(i_m,j_m,k_m)- Gm2(medium) * Idye(j_m) * Idze(k_m) * amp * evolucion(timei,Nodal_Hx%nodSoft(ii))
               end do
            end do
         end do
      end do barridonodalsoftHx
      !
      !
      barridonodalhardHy: do ii=1,Nodal_Hy%numHard
         if (Nodal_Hy%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardHy
         end if
         !
         amp = Nodal_Hy%nodHard(ii)%gridPoint%amplitude
         do k=Nodal_Hy%nodHard(ii)%gridPoint%zi,Nodal_Hy%nodHard(ii)%gridPoint%ze
            k_m = k - b%Hy%ZI
            do j=Nodal_Hy%nodHard(ii)%gridPoint%yi,Nodal_Hy%nodHard(ii)%gridPoint%ye
               j_m = j - b%Hy%YI
               do i=Nodal_Hy%nodHard(ii)%gridPoint%xi,Nodal_Hy%nodHard(ii)%gridPoint%xe
                  i_m = i - b%Hy%XI
                  medium = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medium)%Is%PMC) Hy(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Hy%nodHard(ii))
               end do
            end do
         end do
      end do barridonodalhardHy
      !
      barridonodalsoftHy: do ii=1,Nodal_Hy%numSoft
         if (Nodal_Hy%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftHy
         end if
         !     
         amp = Nodal_Hy%nodSoft(ii)%gridPoint%amplitude
         do k=Nodal_Hy%nodSoft(ii)%gridPoint%zi,Nodal_Hy%nodSoft(ii)%gridPoint%ze
            k_m = k - b%Hy%ZI
            do j=Nodal_Hy%nodSoft(ii)%gridPoint%yi,Nodal_Hy%nodSoft(ii)%gridPoint%ye
               j_m = j - b%Hy%YI
               do i=Nodal_Hy%nodSoft(ii)%gridPoint%xi,Nodal_Hy%nodSoft(ii)%gridPoint%xe
                  i_m = i - b%Hy%XI
                  medium = sggMiHy(i_m,j_m,k_m)
                  if (.not.sgg%Med(medium)%Is%PMC) Hy(i_m,j_m,k_m) = Hy(i_m,j_m,k_m)- Gm2(medium) * Idxe(i_m) * Idze(k_m) * amp * evolucion(timei,Nodal_Hy%nodSoft(ii))
               end do
            end do
         end do
      end do barridonodalsoftHy

      barridonodalhardHz: do ii=1,Nodal_Hz%numHard
         if (Nodal_Hz%nodHard(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalhardHz
         end if
         !
         amp = Nodal_Hz%nodHard(ii)%gridPoint%amplitude
         do k=Nodal_Hz%nodHard(ii)%gridPoint%zi,Nodal_Hz%nodHard(ii)%gridPoint%ze
            k_m = k - b%Hz%ZI
            do j=Nodal_Hz%nodHard(ii)%gridPoint%yi,Nodal_Hz%nodHard(ii)%gridPoint%ye
               j_m = j - b%Hz%YI
               do i=Nodal_Hz%nodHard(ii)%gridPoint%xi,Nodal_Hz%nodHard(ii)%gridPoint%xe
                  i_m = i - b%Hz%XI
                  medium = sggMiHx(i_m,j_m,k_m)
                  if (.not.sgg%Med(medium)%Is%PMC) Hz(i_m,j_m,k_m) = amp * evolucion(timei,Nodal_Hz%nodHard(ii))
               end do
            end do
         end do
      end do barridonodalhardHz
      !
      barridonodalsoftHz: do ii=1,Nodal_Hz%numSoft
         if (Nodal_Hz%nodSoft(ii)%IsInitialValue .and. (timeinstant /= 0)) then
              cycle barridonodalsoftHz
         end if
         !     
         amp = Nodal_Hz%nodSoft(ii)%gridPoint%amplitude
         do k=Nodal_Hz%nodSoft(ii)%gridPoint%zi,Nodal_Hz%nodSoft(ii)%gridPoint%ze
            k_m = k - b%Hz%ZI
            do j=Nodal_Hz%nodSoft(ii)%gridPoint%yi,Nodal_Hz%nodSoft(ii)%gridPoint%ye
               j_m = j - b%Hz%YI
               do i=Nodal_Hz%nodSoft(ii)%gridPoint%xi,Nodal_Hz%nodSoft(ii)%gridPoint%xe
                  i_m = i - b%Hz%XI
                  medium = sggMiHz(i_m,j_m,k_m)
                  if (.not.sgg%Med(medium)%Is%PMC) Hz(i_m,j_m,k_m) = Hz(i_m,j_m,k_m)- Gm2(medium) * Idye(j_m) * Idxe(i_m) * amp * evolucion(timei,Nodal_Hz%nodSoft(ii))
               end do
            end do
         end do
      end do barridonodalsoftHz




      return
   end subroutine AdvancenodalH


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Function to publish the private output data (used in postprocess)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine getnodal(rNodal_Ex,rNodal_Ey,rNodal_Ez,rNodal_Hx,rNodal_Hy,rNodal_Hz)

      type(nodsou_t), pointer :: rNodal_Ex ,rNodal_Ey ,rNodal_Ez
      type(nodsou_t), pointer :: rNodal_Hx ,rNodal_Hy ,rNodal_Hz

      rNodal_Ex  => Nodal_Ex
      rNodal_Ey  => Nodal_Ey
      rNodal_Ez  => Nodal_Ez
      rNodal_Hx  => Nodal_Hx
      rNodal_Hy  => Nodal_Hy
      rNodal_Hz  => Nodal_Hz


      return
   end subroutine

end module nodalsources_m
 