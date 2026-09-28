
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Borders :  PML, PEC, PMC, Periodic handling.
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module BORDERS_other_m
   use FDETYPES_m
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
   subroutine MinusCloneMagneticPMC(sggalloc,sggBorder,Hx,Hy,Hz,c,layoutnumber,num_procs)

      type(XYZlimit_t), dimension(1:6), intent(in)                      :: sggAlloc
      real(kind=RKIND)   , intent(inout) :: &
      Hx(sggalloc(IHX)%XI : sggalloc(IHX)%XE,sggalloc(IHX)%YI : sggalloc(IHX)%YE,sggalloc(IHX)%ZI : sggalloc(IHX)%ZE),&
      Hy(sggalloc(IHY)%XI : sggalloc(IHY)%XE,sggalloc(IHY)%YI : sggalloc(IHY)%YE,sggalloc(IHY)%ZI : sggalloc(IHY)%ZE),&
      Hz(sggalloc(IHZ)%XI : sggalloc(IHZ)%XE,sggalloc(IHZ)%YI : sggalloc(IHZ)%YE,sggalloc(IHZ)%ZI : sggalloc(IHZ)%ZE)

      type(XYZlimit_t), dimension(1:6) :: c
      integer , intent(in) :: layoutnumber,num_procs
      type(Border_t), intent(in)                                         :: sggBorder

      !Hx Down
      if (sggBorder%IsDownPMC) then
         if (layoutnumber == 0)      Hx( : , : ,C(IHX)%ZI-1)=-Hx( : , : ,C(IHX)%ZI)
      end if
      !Hx Up
      if (sggBorder%IsUpPMC) then
         if (layoutnumber == num_procs-1) Hx( : , : ,C(IHX)%ZE+1)=-Hx( : , : ,C(IHX)%ZE)
      end if
      !Hx Left
      if (sggBorder%IsLeftPMC) then
         Hx( : ,C(IHX)%YI-1, : )=-Hx( : ,C(IHX)%YI, : )
      end if
      !Hx Right
      if (sggBorder%IsRightPMC) then
         Hx( : ,C(IHX)%YE+1, : )=-Hx( : ,C(IHX)%YE, : )
      end if
      !Hy Back
      if (sggBorder%IsBackPMC) then
         Hy(C(IHY)%XI-1, : , : )=-Hy(C(IHY)%XI, : , : )
      end if
      !Hy Front
      if (sggBorder%IsFrontPMC) then
         Hy(C(IHY)%XE+1, : , : )=-Hy(C(IHY)%XE, : , : )
      end if
      !Hy Down
      if (sggBorder%IsDownPMC) then
         if (layoutnumber == 0)      Hy( : , : ,C(IHY)%ZI-1)=-Hy( : , : ,C(IHY)%ZI)
      end if
      !Hy Up
      if (sggBorder%IsUpPMC) then
         if (layoutnumber == num_procs-1) Hy( : , : ,C(IHY)%ZE+1)=-Hy( : , : ,C(IHY)%ZE)
      end if
      !
      !Hz Down
      if (sggBorder%IsBackPMC) then
         Hz(C(IHZ)%XI-1, : , : )=-Hz(C(IHZ)%XI, : , : )
      end if
      !Hz Front
      if (sggBorder%IsFrontPMC) then
         Hz(C(IHZ)%XE+1, : , : )=-Hz(C(IHZ)%XE, : , : )
      end if
      !Hz Left
      if (sggBorder%IsLeftPMC) then
         Hz( : ,C(IHZ)%YI-1, : )=-Hz( : ,C(IHZ)%YI, : )
      end if
      !Hz Right
      if (sggBorder%IsRightPMC) then
         Hz( : ,C(IHZ)%YE+1, : )=-Hz( : ,C(IHZ)%YE, : )
      end if
      return
   end subroutine MinusCloneMagneticPMC



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Mirrorizes the Magnetic fields for Periodic
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   subroutine CloneMagneticPeriodic(sggalloc,sggBorder,Hx,Hy,Hz,c,layoutnumber,num_procs)

      type(XYZlimit_t), dimension(1:6), intent(in)                      :: sggAlloc
      real(kind=RKIND)   , intent(inout) :: &
      Hx(sggalloc(IHX)%XI : sggalloc(IHX)%XE,sggalloc(IHX)%YI : sggalloc(IHX)%YE,sggalloc(IHX)%ZI : sggalloc(IHX)%ZE),&
      Hy(sggalloc(IHY)%XI : sggalloc(IHY)%XE,sggalloc(IHY)%YI : sggalloc(IHY)%YE,sggalloc(IHY)%ZI : sggalloc(IHY)%ZE),&
      Hz(sggalloc(IHZ)%XI : sggalloc(IHZ)%XE,sggalloc(IHZ)%YI : sggalloc(IHZ)%YE,sggalloc(IHZ)%ZI : sggalloc(IHZ)%ZE)

      type(XYZlimit_t), dimension(1:6) :: c
      integer(kind=4), intent(in) :: layoutnumber,num_procs
      type(Border_t), intent(in)                                         :: sggBorder

      !Hx Down
      if (sggBorder%IsDownPeriodic) then
         if (layoutnumber == 0)      Hx( : , : ,C(IHX)%ZI-1) = Hx( : , : ,C(IHX)%ZE)
      end if
      !Hx Up
      if (sggBorder%IsUpPeriodic) then
         if (layoutnumber == num_procs-1) Hx( : , : ,C(IHX)%ZE+1) = Hx( : , : ,C(IHX)%ZI)
      end if
      !Hx Left
      if (sggBorder%IsLeftPeriodic) then
         Hx( : ,C(IHX)%YI-1, : ) = Hx( : ,C(IHX)%YE, : )
      end if
      !Hx Right
      if (sggBorder%IsRightPeriodic) then
         Hx( : ,C(IHX)%YE+1, : ) = Hx( : ,C(IHX)%YI, : )
      end if
      !Hy Back
      if (sggBorder%IsBackPeriodic) then
         Hy(C(IHY)%XI-1, : , : ) = Hy(C(IHY)%XE, : , : )
      end if
      !Hy Front
      if (sggBorder%IsFrontPeriodic) then
         Hy(C(IHY)%XE+1, : , : ) = Hy(C(IHY)%XI, : , : )
      end if
      !Hy Down
      if (sggBorder%IsDownPeriodic) then
         if (layoutnumber == 0)      Hy( : , : ,C(IHY)%ZI-1) = Hy( : , : ,C(IHY)%ZE)
      end if
      !Hy Up
      if (sggBorder%IsUpPeriodic) then
         if (layoutnumber == num_procs-1) Hy( : , : ,C(IHY)%ZE+1) = Hy( : , : ,C(IHY)%ZI)
      end if
      !
      !Hz Back
      if (sggBorder%IsBackPeriodic) then
         Hz(C(IHZ)%XI-1, : , : ) = Hz(C(IHZ)%XE, : , : )
      end if
      !Hz Front
      if (sggBorder%IsFrontPeriodic) then
         Hz(C(IHZ)%XE+1, : , : ) = Hz(C(IHZ)%XI, : , : )
      end if
      !Hz Left
      if (sggBorder%IsLeftPeriodic) then
         Hz( : ,C(IHZ)%YI-1, : ) = Hz( : ,C(IHZ)%YE, : )
      end if
      !Hz Right
      if (sggBorder%IsRightPeriodic) then
         Hz( : ,C(IHZ)%YE+1, : ) = Hz( : ,C(IHZ)%YI, : )
      end if
      return
   end subroutine CloneMagneticPeriodic

end module
