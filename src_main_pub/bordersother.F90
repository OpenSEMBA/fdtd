
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Borders :  PML, PEC, PMC, Periodic handling.
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BORDERS_other
   use fdetypes
   implicit none
   private
   !
   public  :: InitOtherBorders, MinusCloneMagneticPMC,CloneMagneticPeriodic

contains
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes PEC and PML data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitOtherBorders(sgg,thereAre)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(logic_control_t), intent(inout) :: thereAre

      thereAre%PeriodicBorders=.false.
      if (sgg%Border%IsBackPeriodic.or.sgg%Border%IsFrontPeriodic.or.sgg%Border%IsLeftPeriodic.or.sgg%Border%IsRightPeriodic.or. &
      sgg%Border%IsUpPeriodic.or.sgg%Border%IsDownPeriodic) thereAre%PeriodicBorders=.true.

      thereAre%PMCBorders=.false.
      if (sgg%Border%IsBackPMC.or.sgg%Border%IsFrontPMC.or.sgg%Border%IsLeftPMC.or.sgg%Border%IsRightPMC.or. &
      sgg%Border%IsUpPMC.or.sgg%Border%IsDownPMC) thereAre%PMCBorders=.true.

      thereAre%PECBorders=.false.
      if (sgg%Border%IsBackPEC.or.sgg%Border%IsFrontPEC.or.sgg%Border%IsLeftPEC.or.sgg%Border%IsRightPEC.or. &
      sgg%Border%IsUpPEC.or.sgg%Border%IsDownPEC) thereAre%PECBorders=.true.
      return
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Mirrorizes the Magnetic fields one cell outside to be used by PMC conditions
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine MinusCloneMagneticPMC(sggalloc,sggBorder,Hx,Hy,Hz,c,layoutnumber,size)

      type(XYZlimit_t), dimension(1:6), intent(in)                      :: sggAlloc
      real(kind=RKIND)   , intent(inout) :: &
      Hx(sggalloc(iHx)%XI : sggalloc(iHx)%XE,sggalloc(iHx)%YI : sggalloc(iHx)%YE,sggalloc(iHx)%ZI : sggalloc(iHx)%ZE),&
      Hy(sggalloc(iHy)%XI : sggalloc(iHy)%XE,sggalloc(iHy)%YI : sggalloc(iHy)%YE,sggalloc(iHy)%ZI : sggalloc(iHy)%ZE),&
      Hz(sggalloc(iHz)%XI : sggalloc(iHz)%XE,sggalloc(iHz)%YI : sggalloc(iHz)%YE,sggalloc(iHz)%ZI : sggalloc(iHz)%ZE)

      type(XYZlimit_t), dimension(1:6) :: c
      integer , intent(in) :: layoutnumber,size
      type(Border_t), intent(in)                                         :: sggBorder

      !Hx Down
      if (sggBorder%IsDownPMC) then
         if (layoutnumber == 0)      Hx( : , : ,C(iHx)%ZI-1)=-Hx( : , : ,C(iHx)%ZI)
      end if
      !Hx Up
      if (sggBorder%IsUpPMC) then
         if (layoutnumber == size-1) Hx( : , : ,C(iHx)%ZE+1)=-Hx( : , : ,C(iHx)%ZE)
      end if
      !Hx Left
      if (sggBorder%IsLeftPMC) then
         Hx( : ,C(iHx)%YI-1, : )=-Hx( : ,C(iHx)%YI, : )
      end if
      !Hx Right
      if (sggBorder%IsRightPMC) then
         Hx( : ,C(iHx)%YE+1, : )=-Hx( : ,C(iHx)%YE, : )
      end if
      !Hy Back
      if (sggBorder%IsBackPMC) then
         Hy(C(iHy)%XI-1, : , : )=-Hy(C(iHy)%XI, : , : )
      end if
      !Hy Front
      if (sggBorder%IsFrontPMC) then
         Hy(C(iHy)%XE+1, : , : )=-Hy(C(iHy)%XE, : , : )
      end if
      !Hy Down
      if (sggBorder%IsDownPMC) then
         if (layoutnumber == 0)      Hy( : , : ,C(iHy)%ZI-1)=-Hy( : , : ,C(iHy)%ZI)
      end if
      !Hy Up
      if (sggBorder%IsUpPMC) then
         if (layoutnumber == size-1) Hy( : , : ,C(iHy)%ZE+1)=-Hy( : , : ,C(iHy)%ZE)
      end if
      !
      !Hz Down
      if (sggBorder%IsBackPMC) then
         Hz(C(iHz)%XI-1, : , : )=-Hz(C(iHz)%XI, : , : )
      end if
      !Hz Front
      if (sggBorder%IsFrontPMC) then
         Hz(C(iHz)%XE+1, : , : )=-Hz(C(iHz)%XE, : , : )
      end if
      !Hz Left
      if (sggBorder%IsLeftPMC) then
         Hz( : ,C(iHz)%YI-1, : )=-Hz( : ,C(iHz)%YI, : )
      end if
      !Hz Right
      if (sggBorder%IsRightPMC) then
         Hz( : ,C(iHz)%YE+1, : )=-Hz( : ,C(iHz)%YE, : )
      end if
      return
   end subroutine MinusCloneMagneticPMC



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Mirrorizes the Magnetic fields for Periodic
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine CloneMagneticPeriodic(sggalloc,sggBorder,Hx,Hy,Hz,c,layoutnumber,size)

      type(XYZlimit_t), dimension(1:6), intent(in)                      :: sggAlloc
      real(kind=RKIND)   , intent(inout) :: &
      Hx(sggalloc(iHx)%XI : sggalloc(iHx)%XE,sggalloc(iHx)%YI : sggalloc(iHx)%YE,sggalloc(iHx)%ZI : sggalloc(iHx)%ZE),&
      Hy(sggalloc(iHy)%XI : sggalloc(iHy)%XE,sggalloc(iHy)%YI : sggalloc(iHy)%YE,sggalloc(iHy)%ZI : sggalloc(iHy)%ZE),&
      Hz(sggalloc(iHz)%XI : sggalloc(iHz)%XE,sggalloc(iHz)%YI : sggalloc(iHz)%YE,sggalloc(iHz)%ZI : sggalloc(iHz)%ZE)

      type(XYZlimit_t), dimension(1:6) :: c
      integer(kind=4), intent(in) :: layoutnumber,size
      type(Border_t), intent(in)                                         :: sggBorder

      !Hx Down
      if (sggBorder%IsDownPeriodic) then
         if (layoutnumber == 0)      Hx( : , : ,C(iHx)%ZI-1) = Hx( : , : ,C(iHx)%ZE)
      end if
      !Hx Up
      if (sggBorder%IsUpPeriodic) then
         if (layoutnumber == size-1) Hx( : , : ,C(iHx)%ZE+1) = Hx( : , : ,C(iHx)%ZI)
      end if
      !Hx Left
      if (sggBorder%IsLeftPeriodic) then
         Hx( : ,C(iHx)%YI-1, : ) = Hx( : ,C(iHx)%YE, : )
      end if
      !Hx Right
      if (sggBorder%IsRightPeriodic) then
         Hx( : ,C(iHx)%YE+1, : ) = Hx( : ,C(iHx)%YI, : )
      end if
      !Hy Back
      if (sggBorder%IsBackPeriodic) then
         Hy(C(iHy)%XI-1, : , : ) = Hy(C(iHy)%XE, : , : )
      end if
      !Hy Front
      if (sggBorder%IsFrontPeriodic) then
         Hy(C(iHy)%XE+1, : , : ) = Hy(C(iHy)%XI, : , : )
      end if
      !Hy Down
      if (sggBorder%IsDownPeriodic) then
         if (layoutnumber == 0)      Hy( : , : ,C(iHy)%ZI-1) = Hy( : , : ,C(iHy)%ZE)
      end if
      !Hy Up
      if (sggBorder%IsUpPeriodic) then
         if (layoutnumber == size-1) Hy( : , : ,C(iHy)%ZE+1) = Hy( : , : ,C(iHy)%ZI)
      end if
      !
      !Hz Back
      if (sggBorder%IsBackPeriodic) then
         Hz(C(iHz)%XI-1, : , : ) = Hz(C(iHz)%XE, : , : )
      end if
      !Hz Front
      if (sggBorder%IsFrontPeriodic) then
         Hz(C(iHz)%XE+1, : , : ) = Hz(C(iHz)%XI, : , : )
      end if
      !Hz Left
      if (sggBorder%IsLeftPeriodic) then
         Hz( : ,C(iHz)%YI-1, : ) = Hz( : ,C(iHz)%YE, : )
      end if
      !Hz Right
      if (sggBorder%IsRightPeriodic) then
         Hz( : ,C(iHz)%YE+1, : ) = Hz( : ,C(iHz)%YI, : )
      end if
      return
   end subroutine CloneMagneticPeriodic

end module
