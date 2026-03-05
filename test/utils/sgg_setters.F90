module mod_sggMethods
   use FDETYPES
   implicit none
   private
   

   public :: sgg_init

   public :: sgg_set_tiempo
   public :: sgg_set_dt
   public :: sgg_set_extraswitches

   public :: sgg_set_NumMedia
   public :: sgg_set_AllocMed
   public :: sgg_set_IniPMLMedia
   public :: sgg_set_EndPMLMedia

   public :: sgg_set_NumPlaneWaves
   public :: sgg_set_TimeSteps
   public :: sgg_set_InitialTimeStep

   public :: sgg_set_NumNodalSources
   public :: sgg_set_NumberRequest

   public :: sgg_set_LineX
   public :: sgg_set_LineY
   public :: sgg_set_LineZ
   public :: sgg_set_DX
   public :: sgg_set_DY
   public :: sgg_set_DZ

   public :: sgg_set_AllocDxI
   public :: sgg_set_AllocDyI
   public :: sgg_set_AllocDzI
   public :: sgg_set_AllocDxE
   public :: sgg_set_AllocDyE
   public :: sgg_set_AllocDzE

   public :: sgg_set_PlaneWave
   public :: sgg_set_Border
   public :: sgg_set_PML
   public :: sgg_set_Eshared
   public :: sgg_set_Hshared

   public :: sgg_set_Alloc
   public :: sgg_set_Sweep
   public :: sgg_set_SINPMLSweep

   public :: sgg_set_Med
   public :: sgg_set_NodalSource
   public :: sgg_set_Observation

   public :: sgg_set_thereAreMagneticMedia
   public :: sgg_set_thereArePMLMagneticMedia

   public :: sgg_set_nEntradaRoot
   public :: sgg_set_Punto

   public :: sgg_add_observation
contains
   subroutine sgg_init(obj, &
      tiempo, dt, extraswitches, &
      NumMedia, AllocMed, &
      IniPMLMedia, EndPMLMedia, &
      NumPlaneWaves, TimeSteps, InitialTimeStep, &
      NumNodalSources, NumberRequest, &
      thereAreMagneticMedia, thereArePMLMagneticMedia, &
      nEntradaRoot)

      implicit none

      type(SGGFDTDINFO), intent(inout) :: obj

      ! ===== Optional arguments =====
      real(kind=RKIND_tiempo), pointer, optional :: tiempo(:)
      real(kind=RKIND_tiempo), optional          :: dt
      character(len=*), optional          :: extraswitches

      integer(kind=SINGLE), optional :: NumMedia, AllocMed
      integer(kind=SINGLE), optional :: IniPMLMedia, EndPMLMedia
      integer(kind=SINGLE), optional :: NumPlaneWaves, TimeSteps, InitialTimeStep
      integer(kind=SINGLE), optional :: NumNodalSources, NumberRequest

      logical, optional :: thereAreMagneticMedia
      logical, optional :: thereArePMLMagneticMedia

      character(len=*), optional :: nEntradaRoot

      ! ===== Defaults =====

      nullify (obj%tiempo)
      obj%dt = 0.0_RKIND_tiempo
      obj%extraswitches = ""

      obj%NumMedia = 0_SINGLE
      obj%AllocMed = 0_SINGLE
      obj%IniPMLMedia = 0_SINGLE
      obj%EndPMLMedia = 0_SINGLE
      obj%NumPlaneWaves = 0_SINGLE
      obj%TimeSteps = 0_SINGLE
      obj%InitialTimeStep = 0_SINGLE
      obj%NumNodalSources = 0_SINGLE
      obj%NumberRequest = 0_SINGLE

      nullify (obj%LineX, obj%LineY, obj%LineZ)
      nullify (obj%DX, obj%DY, obj%DZ)

      obj%AllocDxI = 0_SINGLE
      obj%AllocDyI = 0_SINGLE
      obj%AllocDzI = 0_SINGLE
      obj%AllocDxE = 0_SINGLE
      obj%AllocDyE = 0_SINGLE
      obj%AllocDzE = 0_SINGLE

      nullify (obj%PlaneWave)
      nullify (obj%Med)
      nullify (obj%NodalSource)
      nullify (obj%Observation)

      obj%thereAreMagneticMedia = .false.
      obj%thereArePMLMagneticMedia = .false.

      obj%nEntradaRoot = ""

      ! NOTE:
      ! Derived-type components (Border, PML, Shared_t, XYZlimit_t, Punto)
      ! are automatically default-initialized if they define their own defaults.

      ! ===== Overrides from arguments =====

      if (present(tiempo)) obj%tiempo => tiempo
      if (present(dt)) obj%dt = dt
      if (present(extraswitches)) obj%extraswitches = extraswitches

      if (present(NumMedia)) obj%NumMedia = NumMedia
      if (present(AllocMed)) obj%AllocMed = AllocMed
      if (present(IniPMLMedia)) obj%IniPMLMedia = IniPMLMedia
      if (present(EndPMLMedia)) obj%EndPMLMedia = EndPMLMedia
      if (present(NumPlaneWaves)) obj%NumPlaneWaves = NumPlaneWaves
      if (present(TimeSteps)) obj%TimeSteps = TimeSteps
      if (present(InitialTimeStep)) obj%InitialTimeStep = InitialTimeStep
      if (present(NumNodalSources)) obj%NumNodalSources = NumNodalSources
      if (present(NumberRequest)) obj%NumberRequest = NumberRequest

      if (present(thereAreMagneticMedia)) &
         obj%thereAreMagneticMedia = thereAreMagneticMedia

      if (present(thereArePMLMagneticMedia)) &
         obj%thereArePMLMagneticMedia = thereArePMLMagneticMedia

      if (present(nEntradaRoot)) obj%nEntradaRoot = nEntradaRoot

   end subroutine sgg_init

   subroutine sgg_set_tiempo(sgg, tiempo)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND_tiempo), pointer :: tiempo(:)
      sgg%tiempo => tiempo
   end subroutine

   subroutine sgg_set_dt(sgg, dt)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND_tiempo), intent(in) :: dt
      sgg%dt = dt
   end subroutine

   subroutine sgg_set_extraswitches(sgg, extraswitches)
      type(SGGFDTDINFO), intent(inout) :: sgg
      character(len=*), intent(in) :: extraswitches
      sgg%extraswitches = extraswitches
   end subroutine

   subroutine sgg_set_NumMedia(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%NumMedia = newValue
   end subroutine

   subroutine sgg_set_AllocMed(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocMed = newValue
   end subroutine

   subroutine sgg_set_IniPMLMedia(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%IniPMLMedia = newValue
   end subroutine

   subroutine sgg_set_EndPMLMedia(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%EndPMLMedia = newValue
   end subroutine

   subroutine sgg_set_NumPlaneWaves(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%NumPlaneWaves = newValue
   end subroutine

   subroutine sgg_set_TimeSteps(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%TimeSteps = newValue
   end subroutine

   subroutine sgg_set_InitialTimeStep(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%InitialTimeStep = newValue
   end subroutine

   subroutine sgg_set_NumNodalSources(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%NumNodalSources = newValue
   end subroutine

   subroutine sgg_set_NumberRequest(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%NumberRequest = newValue
   end subroutine

   subroutine sgg_set_LineX(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND), pointer :: newValue(:)
      sgg%LineX => newValue
   end subroutine

   subroutine sgg_set_LineY(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND), pointer :: newValue(:)
      sgg%LineY => newValue
   end subroutine

   subroutine sgg_set_LineZ(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND), pointer :: newValue(:)
      sgg%LineZ => newValue
   end subroutine

   subroutine sgg_set_DX(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND), pointer :: newValue(:)
      sgg%DX => newValue
   end subroutine

   subroutine sgg_set_DY(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND), pointer :: newValue(:)
      sgg%DY => newValue
   end subroutine

   subroutine sgg_set_DZ(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      real(kind=RKIND), pointer :: newValue(:)
      sgg%DZ => newValue
   end subroutine

   subroutine sgg_set_AllocDxI(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocDxI = newValue
   end subroutine

   subroutine sgg_set_AllocDyI(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocDyI = newValue
   end subroutine

   subroutine sgg_set_AllocDzI(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocDzI = newValue
   end subroutine

   subroutine sgg_set_AllocDxE(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocDxE = newValue
   end subroutine

   subroutine sgg_set_AllocDyE(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocDyE = newValue
   end subroutine

   subroutine sgg_set_AllocDzE(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      integer(kind=SINGLE), intent(in) :: newValue
      sgg%AllocDzE = newValue
   end subroutine

   subroutine sgg_set_PlaneWave(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(planeonde_t), pointer :: newValue(:)
      sgg%PlaneWave => newValue
   end subroutine

   subroutine sgg_set_Med(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(MediaData_t), pointer :: newValue(:)
      sgg%Med => newValue
   end subroutine

   subroutine sgg_set_NodalSource(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(NodalSource_t), pointer :: newValue(:)
      sgg%NodalSource => newValue
   end subroutine

   subroutine sgg_set_Observation(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(obses_t), pointer :: newValue(:)
      sgg%Observation => newValue
   end subroutine

   subroutine sgg_set_Border(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(Border_t), intent(in) :: newValue
      sgg%Border = newValue
   end subroutine

   subroutine sgg_set_PML(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(PML_t), intent(in) :: newValue
      sgg%PML = newValue
   end subroutine

   subroutine sgg_set_Eshared(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(Shared_t), intent(in) :: newValue
      sgg%Eshared = newValue
   end subroutine

   subroutine sgg_set_Hshared(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(Shared_t), intent(in) :: newValue
      sgg%Hshared = newValue
   end subroutine

   subroutine sgg_set_Alloc(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(XYZlimit_t), intent(in) :: newValue(1:6)
      sgg%Alloc = newValue
   end subroutine

   subroutine sgg_set_Sweep(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(XYZlimit_t), intent(in) :: newValue(1:6)
      sgg%Sweep = newValue
   end subroutine

   subroutine sgg_set_SINPMLSweep(sgg, newValue)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(XYZlimit_t), intent(in) :: newValue(1:6)
      sgg%SINPMLSweep = newValue
   end subroutine

   subroutine sgg_set_thereAreMagneticMedia(sgg, value)
      type(SGGFDTDINFO), intent(inout) :: sgg
      logical, intent(in) :: value
      sgg%thereAreMagneticMedia = value
   end subroutine

   subroutine sgg_set_thereArePMLMagneticMedia(sgg, value)
      type(SGGFDTDINFO), intent(inout) :: sgg
      logical, intent(in) :: value
      sgg%thereArePMLMagneticMedia = value
   end subroutine

   subroutine sgg_set_nEntradaRoot(sgg, value)
      type(SGGFDTDINFO), intent(inout) :: sgg
      character(len=*), intent(in) :: value
      sgg%nEntradaRoot = value
   end subroutine

   subroutine sgg_set_Punto(sgg, value)
      type(SGGFDTDINFO), intent(inout) :: sgg
      type(coorsxyzP), intent(in) :: value
      sgg%Punto = value
   end subroutine

   subroutine sgg_add_observation(sgg, new_observation)
      implicit none

      type(SGGFDTDINFO), intent(inout) :: sgg
      type(Obses_t), intent(in), target :: new_observation

      type(Obses_t), dimension(:), pointer :: temp_obs
      integer :: old_size, new_size

      old_size = sgg%NumberRequest
      new_size = old_size + 1

      allocate (temp_obs(1:new_size))

      if (old_size > 0) then
         temp_obs(1:old_size) = sgg%Observation(1:old_size)
         deallocate (sgg%Observation)
      end if

      temp_obs(new_size) = new_observation

      sgg%Observation => temp_obs

      sgg%NumberRequest = new_size

   end subroutine sgg_add_observation
   
end module mod_sggMethods
