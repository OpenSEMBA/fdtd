
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Observation module to store the observed data
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Observa
  use FDETYPES

#ifdef CompileWithMPI
  use MPIcomm
#endif

  use wiresHolland_constants
  use HollandWires

#ifdef CompileWithMTLN
  use Wire_bundles_mtln_mod
  use mtln_solver_mod, mtln_solver_t => mtln_t
#endif
#ifdef CompileWithBerengerWires
  use WiresBerenger
#endif
#ifdef CompileWithSlantedWires
  use WiresSlanted
  use WiresSlanted_Types
  use WiresSlanted_Constants
#endif
  use report
  use farfield_m
  use nodalsources
!
  IMPLICIT NONE
  private

  type Serialized_t
    REAL(KIND=RKIND), pointer, dimension(:, :)   ::  valor, valor_x, valor_y, valor_z ! (step, valor)
    REAL(KIND=RKIND), pointer, dimension(:, :)   ::  valorE, valor_Ex, valor_Ey, valor_Ez ! (step, valor)
    REAL(KIND=RKIND), pointer, dimension(:, :)   ::  valorH, valor_Hx, valor_Hy, valor_Hz ! (step, valor)
    INTEGER(kind=4), POINTER, DIMENSION(:) :: eI, eJ, eK, currentType, sggmtag
    complex(kind=CKIND), dimension(:, :), allocatable  :: valorComplex_x, valorComplex_y, valorComplex_z
    complex(kind=CKIND), dimension(:, :), allocatable  :: valorComplex_Ex, valorComplex_Ey, valorComplex_Ez
    complex(kind=CKIND), dimension(:, :), allocatable  :: valorComplex_Hx, valorComplex_Hy, valorComplex_Hz

  contains
    procedure :: allocate_for_time_domain
    procedure :: allocate_for_frequency_domain
    procedure :: allocate_current_value
    procedure :: deallocate_for_time_domain
    procedure :: deallocate_for_frequency_domain
    procedure :: deallocate_current_value
  end type Serialized_t
  type item_t

    type(CurrentSegments), pointer  ::  segmento !segmento de hilo que se observa si lo hubiere

#ifdef CompileWithBerengerWires
    type(TSegment), pointer  ::  segmento_Berenger !segmento de hilo que se observa si lo hubiere
#endif
#ifdef CompileWithSlantedWires
    class(Segment), pointer  ::  segmento_Slanted !segmento de hilo que se observa si lo hubiere
#endif
    character(LEN=BUFSIZE)  ::  path
    integer(kind=4) :: unit, unitmaster !to store the unit of the file y en caso de singlefileginario el unitmaster que escribe
    integer(kind=4) :: columnas !number of columns in the output file
    REAL(KIND=RKIND), pointer, dimension(:)   ::  valor, valor2, valor3, valor4, valor5 !stored values at each time step !not read but calculate !210521 also store -edl+vdrop
    REAL(KIND=RKIND)    ::  valorsigno !just to store the sign of the current for the wires
    REAL(KIND=RKIND), pointer, dimension(:, :, :, :)   ::  valor3D !stored values at each time step !not read but calculate
    type(Serialized_t)  ::  Serialized !para almecenar valores serializados en volumenes en vez de bulk
    !freqdomain probles
    complex(kind=CKIND), dimension(:, :, :, :, :), allocatable  :: valor3DComplex !freqdomain probes
#ifdef CompileWithMPI
    integer(kind=4)      ::  MPISubcomm, MPIRoot, MPIGroupIndex
    integer(kind=4)      :: ZIorig, ZEorig
#endif
    integer(kind=4)      :: Xtrancos, Ytrancos, Ztrancos
    integer(kind=4)      :: XItrancos, YItrancos, ZItrancos
    integer(kind=4)      :: XEtrancos, YEtrancos, ZEtrancos
  end type

  type output_t
    type(item_t), dimension(:), pointer  ::  item !path con el output y sus valores
    integer(kind=4)      ::  Trancos
    logical  ::  SaveAll
    integer(kind=4) :: TimesWritten !to control the volumic probes
    !freqdomain probles
    integer(KIND=4) :: NumFreqs
    real(kind=Rkind), dimension(:), allocatable  :: Freq
    real(kind=Rkind) :: InitialFreq, FinalFreq, FreqStep
    complex(kind=CKIND), dimension(:), allocatable  :: auxExp_E, auxExp_H, dftEntrada   !para sondas freqdomain
  end type output_t

  type(Thinwires_t), pointer  ::  Hwireslocal
#ifdef CompileWithBerengerWires
  type(TWires), pointer  ::  Hwireslocal_Berenger
#endif
#ifdef CompileWithSlantedWires
  type(WiresData), pointer  ::  Hwireslocal_Slanted
#endif

#ifdef CompileWithMPI
  REAL(KIND=RKIND), pointer, dimension(:)   ::  valores, newvalores  !auxiliary for Bloque currents sync
#endif
!!!variables globales del modulo
  REAL(KIND=RKIND), save           ::  eps0, mu0
!!!
   !!!!!!!!!variables local

  REAL(KIND=RKIND), pointer, dimension(:), save  ::  InvEps, InvMu
  type(output_t), pointer, dimension(:), save  ::  output

  public InitObservation, FlushObservationFiles, UpdateObservation, DestroyObservation, CloseObservationFiles, unpacksinglefiles, &
    GetOutput, preprocess_observation
  public output_t, item_t, Serialized_t, dtft
#ifdef CompileWithMTLN
  public FlushMTLNObservationFiles
#endif
contains

  SUBROUTINE allocate_for_time_domain(this, numberOfSerialized)
    class(Serialized_t), intent(inout) :: this
    integer(kind=4) :: numberOfSerialized

    ALLOCATE (this%Valor(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_x(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_y(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_z(1, 1:numberOfSerialized))

    ALLOCATE (this%ValorE(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_Ex(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_Ey(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_Ez(1, 1:numberOfSerialized))

    ALLOCATE (this%ValorH(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_Hx(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_Hy(1, 1:numberOfSerialized))
    ALLOCATE (this%Valor_Hz(1, 1:numberOfSerialized))

    this%Valor = 0.
    this%Valor_x = 0.
    this%Valor_y = 0.
    this%Valor_z = 0.

    this%ValorE = 0.
    this%Valor_Ex = 0.
    this%Valor_Ey = 0.
    this%Valor_Ez = 0.

    this%ValorH = 0.
    this%Valor_Hx = 0.
    this%Valor_Hy = 0.
    this%Valor_Hz = 0.

  END SUBROUTINE

  SUBROUTINE deallocate_for_time_domain(this)
    class(Serialized_t), intent(inout) :: this

    DEALLOCATE (this%Valor)
    DEALLOCATE (this%Valor_x)
    DEALLOCATE (this%Valor_y)
    DEALLOCATE (this%Valor_z)

    DEALLOCATE (this%ValorE)
    DEALLOCATE (this%Valor_Ex)
    DEALLOCATE (this%Valor_Ey)
    DEALLOCATE (this%Valor_Ez)

    DEALLOCATE (this%ValorH)
    DEALLOCATE (this%Valor_Hx)
    DEALLOCATE (this%Valor_Hy)
    DEALLOCATE (this%Valor_Hz)

  END SUBROUTINE

  SUBROUTINE allocate_for_frequency_domain(this, numberOfSerialized)
    class(Serialized_t), intent(inout) :: this
    integer(kind=4) :: numberOfSerialized

    call this%allocate_for_time_domain(numberOfSerialized)

    ALLOCATE (this%ValorComplex_x(1, 1:numberOfSerialized))
    ALLOCATE (this%ValorComplex_y(1, 1:numberOfSerialized))
    ALLOCATE (this%ValorComplex_z(1, 1:numberOfSerialized))

    ALLOCATE (this%ValorComplex_Ex(1, 1:numberOfSerialized))
    ALLOCATE (this%ValorComplex_Ey(1, 1:numberOfSerialized))
    ALLOCATE (this%ValorComplex_Ez(1, 1:numberOfSerialized))

    ALLOCATE (this%ValorComplex_Hx(1, 1:numberOfSerialized))
    ALLOCATE (this%ValorComplex_Hy(1, 1:numberOfSerialized))
    ALLOCATE (this%ValorComplex_Hz(1, 1:numberOfSerialized))

    this%ValorComplex_x = 0.
    this%ValorComplex_y = 0.
    this%ValorComplex_z = 0.

    this%ValorComplex_Ex = 0.
    this%ValorComplex_Ey = 0.
    this%ValorComplex_Ez = 0.

    this%ValorComplex_Hx = 0.
    this%ValorComplex_Hy = 0.
    this%ValorComplex_Hz = 0.

  END SUBROUTINE

  SUBROUTINE deallocate_for_frequency_domain(this)
    class(Serialized_t), intent(inout) :: this
    call this%deallocate_for_time_domain()

    DEALLOCATE (this%ValorComplex_x)
    DEALLOCATE (this%ValorComplex_y)
    DEALLOCATE (this%ValorComplex_z)

    DEALLOCATE (this%ValorComplex_Ex)
    DEALLOCATE (this%ValorComplex_Ey)
    DEALLOCATE (this%ValorComplex_Ez)

    DEALLOCATE (this%ValorComplex_Hx)
    DEALLOCATE (this%ValorComplex_Hy)
    DEALLOCATE (this%ValorComplex_Hz)

  END SUBROUTINE

  SUBROUTINE allocate_current_value(this, numberOfSerialized)
    class(Serialized_t), intent(inout) :: this
    integer(kind=4) :: numberOfSerialized

    ALLOCATE (this%eI(1:numberOfSerialized))
    ALLOCATE (this%eJ(1:numberOfSerialized))
    ALLOCATE (this%eK(1:numberOfSerialized))

    ALLOCATE (this%currentType(1:numberOfSerialized))
    ALLOCATE (this%sggMtag(1:numberOfSerialized))

    this%eI = 0
    this%eJ = 0
    this%eK = 0

    this%currentType = 0
    this%sggMtag = 0
  END SUBROUTINE

  SUBROUTINE deallocate_current_value(this)
    class(Serialized_t), intent(inout) :: this
    DEALLOCATE (this%eI)
    DEALLOCATE (this%eJ)
    DEALLOCATE (this%eK)

    DEALLOCATE (this%currentType)
    DEALLOCATE (this%sggMtag)
  END SUBROUTINE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initializes observation stuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine preprocess_observation(observation, privateOutput, time, finaltimestep, dt, saveall)
    type(Obses_t), intent(inout) :: observation
    type(output_t), intent(inout) :: privateOutput
    real(kind=RKIND_tiempo), pointer, dimension(:), intent(in) :: time
    real(kind=RKIND_tiempo), intent(in) :: dt
    integer(kind=4) :: finaltimestep
    logical :: saveall

    observation%done = .false.
    observation%begun = .false.
    observation%flushed = .false.

    observation%TimeStep = max(observation%TimeStep, dt)

    if (10.0_RKIND*(observation%FinalTime - observation%InitialTime)/min(dt, observation%TimeStep) >= huge(1_4)) then
      observation%FinalTime = observation%InitialTime + min(dt, observation%TimeStep)*huge(1_4)/10.0_RKIND
    end if

    if (observation%InitialTime < observation%TimeStep) then
      observation%InitialTime = 0.0_RKIND !para que saque tambien el instante inicial
    end if

    if (observation%TimeStep > (observation%FinalTime - observation%InitialTime)) then

      if (observation%P(1)%what == mapvtk) then
        observation%FinalTime = 0.0_RKIND
        observation%InitialTime = 0.0_RKIND
      else
        observation%FinalTime = observation%InitialTime + observation%TimeStep
      end if
    end if

    observation%InitialTime = int(observation%InitialTime)
    observation%FinalTime = int(observation%FinalTime)
    observation%FreqStep = min(observation%FreqStep, 2.0_RKIND/dt)
    if ((observation%FreqStep > observation%FinalFreq - observation%InitialFreq) .or. (observation%FreqStep == 0)) then
      observation%FreqStep = observation%FinalFreq - observation%InitialFreq
      observation%FinalFreq = observation%InitialFreq + observation%FreqStep
    end if
    if (.not. observation%Volumic) then
      observation%Saveall = observation%Saveall .or. saveall
      privateOutput%SaveAll = observation%Saveall
    else
      privateOutput%SaveAll = .false.
      observation%Saveall = .false.
    end if
#ifdef miguelConformalStandAlone
    privateOutput%SaveAll = .false.
#endif
   if (observation%nP /= 0) then
      if (observation%P(1)%what == mapvtk) then
        privateOutput%SaveAll = .false.
        observation%Saveall = .false.
      end if
    end if

    if (observation%Saveall) then
      privateOutput%Trancos = 1
      observation%InitialTime = 0.0_RKIND
      observation%FinalTime = time(finaltimestep + 2)
    else
      privateOutput%Trancos = max(1, int(observation%TimeStep/dt))
      observation%InitialTime = max(0.0_RKIND, observation%InitialTime)
      observation%FinalTime = min(time(finaltimestep + 2), observation%FinalTime) !CLIPEA
      if (observation%FinalTime < observation%InitialTime) then
        observation%FinalTime = observation%InitialTime
      end if
    end if
    
!!!!
  end subroutine preprocess_observation

  subroutine eliminate_unnecesary_observation_points(observation_probe, output_item, sweep, SINPMLSweep, ZI, ZE, layoutnumber, size)
   type(item_t),  intent(inout) :: output_item
   type(observable_t),  intent(inout) :: observation_probe
   type(XYZlimit_t), dimension(1:6), intent(in) :: sweep, SINPMLSweep
   integer(kind=4), intent(in) :: ZI, ZE, layoutnumber, size
   integer(kind=4) :: field
   
   output_item%Xtrancos = observation_probe%Xtrancos
   output_item%Ytrancos = observation_probe%Ytrancos
   output_item%Ztrancos = observation_probe%Ztrancos

   output_item%XItrancos = int(observation_probe%XI/output_item%Xtrancos)
   output_item%YItrancos = int(observation_probe%YI/output_item%Ytrancos)
   output_item%ZItrancos = int(observation_probe%ZI/output_item%Ztrancos)

   output_item%XEtrancos = int(observation_probe%XE/output_item%Xtrancos)
   output_item%YEtrancos = int(observation_probe%YE/output_item%Ytrancos)
   output_item%ZEtrancos = int(observation_probe%ZE/output_item%Ztrancos)

   if (mod(observation_probe%XI,output_item%Xtrancos) /= 0) output_item%XItrancos=output_item%XItrancos+1
   if (mod(observation_probe%YI,output_item%Ytrancos) /= 0) output_item%YItrancos=output_item%YItrancos+1
   if (mod(observation_probe%ZI,output_item%Ztrancos) /= 0) output_item%ZItrancos=output_item%ZItrancos+1

#ifdef CompileWithMPI
   output_item%MPISubComm = -1 !just to void it
#endif

   field = observation_probe%What
   select case (field)
      case (iBloqueJx, iBloqueJy, iBloqueMx, iBloqueMy)
         call eliminate_observation_from_block(observation_probe, output_item, sweep, field)
      case (iEx, iVx, iEy, iVy, iHz, iBloqueMz, iJx, iJy, iQx, iQy)
         !in case of MPI the flushing is only cared by one of the sharing layouts
         !este es el unico caso en el que un punto es susceptible de ser escrito por dos layouts. Por eso se lo echo
         ! solo a uno de ellos: al de abajo (a menos que que sea el layout de mas arriba, en cuyo caso tiene que tratarlo el) !bug del itc2 con el pathx hasta el borde
         if (((observation_probe%ZI >= sweep(fieldo(field, 'Z'))%ZE) .and. (layoutnumber /= size - 1)) .or. &
            (observation_probe%ZI < sweep(fieldo(field, 'Z'))%ZI)) then
            observation_probe%What = Nothing !do not observe anything
         end if
      case (iEz, iVz, iJz, iQz, iBloqueJz, iHx, iHy)
         if ((observation_probe%ZI > sweep(fieldo(field, 'Z'))%ZE) .or. &
              (observation_probe%ZI < sweep(fieldo(field, 'Z'))%ZI)) then
            observation_probe%What = nothing !do not observe anything
         end if
      case (iExC, iEyC, iHzC, iMhC, iEzC, iHxC, iHyC, iMeC)
         call eliminate_observation_from_block(observation_probe, output_item, sweep, field)
      case (iCur, iCurX, iCurY, iCurZ, mapvtk)
         call eliminate_observation_from_current(observation_probe, output_item, sweep, field)
      case (FarField)
         call eliminate_observation_from_farfield(observation_probe, output_item, SINPMLSweep, field, ZI, ZE)
   end select
  end subroutine eliminate_unnecesary_observation_points

  subroutine eliminate_observation_from_block(observation_probe, output_item, sweep, field)
   type(item_t), intent(inout) :: output_item
   type(observable_t), intent(inout) :: observation_probe
   type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
   integer, intent(in) :: field

   if ((observation_probe%ZI > sweep(fieldo(field, 'Z'))%ZE) .or. &
   (observation_probe%ZE < sweep(fieldo(field, 'Z'))%ZI)) then
         observation_probe%What = nothing

#ifdef CompileWithMPI
      output_item%MPISubComm = -1 !just to void it
   else
      output_item%MPISubComm = 1
   end if
      output_item%MPIRoot = 0
      If ((observation_probe%ZI >= sweep(fieldo(field, 'Z'))%ZI) .and. &
         (observation_probe%ZI <= sweep(fieldo(field, 'Z'))%ZE)) then
         output_item%MPIRoot = layoutnumber
      end if
   !all of them must call the init routine even if they do not sync
   call MPIinitSubcomm(layoutnumber, size, &
         output_item%MPISubComm, output_item%MPIRoot, output_item%MPIGroupIndex)
#else
   end if
#endif

  end subroutine eliminate_observation_from_block

  subroutine eliminate_observation_from_electric_current(observation_probe, output_item, sweep, field)
   type(item_t), intent(inout) :: output_item
   type(observable_t), intent(inout) :: observation_probe
   type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
   integer, intent(in) :: field

   if ((observation_probe%ZI > sweep(fieldo(field, 'Z'))%ZE) .or. &
         (observation_probe%ZE < sweep(fieldo(field, 'Z'))%ZI)) then
      observation_probe%What = nothing
#ifdef CompileWithMPI
      output_item%MPISubComm = -1
   else
      output_item%MPISubComm = 1
   end if
   output_item%MPIRoot = 0
   If ((observation_probe%ZI >= sweep(fieldo(field, 'Z'))%ZI) .and. &
         (observation_probe%ZI <= sweep(fieldo(field, 'Z'))%ZE)) then
      output_item%MPIRoot = layoutnumber
   end if
   call MPIinitSubcomm(layoutnumber, size, &
      output_item%MPISubComm, output_item%MPIRoot, output_item%MPIGroupIndex)
#else
   end if
#endif

   end subroutine eliminate_observation_from_electric_current

   subroutine eliminate_observation_from_current(observation_probe, output_item, sweep, field)
   type(item_t), intent(inout) :: output_item
   type(observable_t), intent(inout) :: observation_probe
   type(XYZlimit_t), dimension(1:6), intent(in) :: sweep
   integer, intent(in) :: field

   if ((observation_probe%ZI >= sweep(iHz)%ZE) .or. &
         (observation_probe%ZE < sweep(iHZ)%ZI)) then
      observation_probe%What = nothing
#ifdef CompileWithMPI
      output_item%MPISubComm = -1
   else
      output_item%MPISubComm = 1
   end if
      !clipeo los finales porque luego tengo que interpolar y el MPI me puede molestar 06/07/15
   if ((field == icur) .or. (field == icurX) .or. (field == icurY) .or. (field == mapvtk)) then
      observation_probe%ZE = Min(observation_probe%ZE, sweep(iHx)%ZE)
   end if
   
   output_item%MPIRoot = 0
   If ((observation_probe%ZI >= sweep(fieldo(field, 'Z'))%ZI) .and. &
         (observation_probe%ZI <= sweep(fieldo(field, 'Z'))%ZE)) then
      output_item%MPIRoot = layoutnumber
   end if
   call MPIinitSubcomm(layoutnumber, size, &
      output(ii)%item(i)%MPISubComm, output(ii)%item(i)%MPIRoot, output(ii)%item(i)%MPIGroupIndex)
#else
   end if
#endif

   end subroutine eliminate_observation_from_current

   subroutine eliminate_observation_from_farfield(observation_probe, output_item, SINPMLSweep, field, ZI, ZE)
   type(item_t), intent(inout) :: output_item
   type(observable_t), intent(inout) :: observation_probe
   type(XYZlimit_t), dimension(1:6), intent(in) :: SINPMLSweep
   integer, intent(in) :: field
   INTEGER(kind=4), intent(in) :: ZI, ZE

   if ((ZI > SINPMLSweep(IHz)%ZE) .or. (ZE < SINPMLSweep(iHz)%ZI)) then   !MPI NO DUPLICAR CALCULOS
      observation_probe%What = nothing
#ifdef CompileWithMPI
      output_item%MPISubComm = -1 !just to void it
   else
      output_item%MPISubComm = 1
   end if
      output_item%MPIRoot = 0
   If ((observation_probe%ZI >= SINPMLSweep(iHz)%ZI) .and. &
         (observation_probe%ZI < SINPMLSweep(iHz)%ZE)) then
      output_item%MPIRoot = layoutnumber
   end if

   call MPIinitSubcomm(layoutnumber, size, &
         output_item%MPISubComm, output_item%MPIRoot, output_item%MPIGroupIndex)
#else
   end if
#endif
      
   end subroutine eliminate_observation_from_farfield

  subroutine InitObservation(sgg, media, tag_numbers, &
                  ThereAreObservation, ThereAreWires, ThereAreFarFields, resume, initialtimestep, finaltimestep, lastexecutedtime, &
                             nEntradaRoot, layoutnumber, size, saveall, singlefilewrite, wiresflavor, &
                             SINPML_fullsize, facesNF2FF, NF2FFDecim, eps00, mu00, simu_devia, mpidir, niapapostprocess, b)
    !solo lo precisa de entrada farfield
    type(media_matrices_t), intent(in) :: media
    type(bounds_t)  ::  b
    type(SGGFDTDINFO), intent(IN)         ::  sgg
    type(taglist_t) :: tag_numbers
    logical :: simu_devia, niapapostprocess
    REAL(KIND=RKIND)           ::  eps00, mu00
    !---------------------------> inputs <----------------------------------------------------------
    integer(kind=4), intent(in) :: layoutnumber, size, mpidir
    type(nf2ff_t) :: facesNF2FF
    type(limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
    character(len=*), INTENT(in) :: wiresflavor
    logical  ::  saveall, singlefilewrite, NF2FFDecim, INIT, GEOM, ASIGNA, electric, magnetic
    character(LEN=BUFSIZE)  ::  p1, p2
    real(kind=RKIND_tiempo) :: lastexecutedtime

    character(len=*), intent(in)  ::  nEntradaRoot

    integer (kind=4)  ::  i,field,ii,i1,j1,k1,n,i2,j2,k2,initialtimestep, finaltimestep,NO,NO2,iwi,iwj,compo,ntime,ntimeforvolumic,iff1,i0t
    integer(kind=4) :: Efield, HField
    logical, intent(inout)   ::  ThereAreObservation, ThereAreFarFields
    logical, intent(in)      ::  ThereAreWires, resume
    character(LEN=BUFSIZE)  ::  chari, charj, chark, chari2, charj2, chark2, charNO
    character(LEN=BUFSIZE)  ::  ext, extpoint, adum, prefix_field
    logical  ::  incident, errnofile, first
    REAL(KIND=RKIND)    ::  rdum, field1, field2
    REAL(KIND=RKIND_tiempo)    ::  at, dtevol, tiempo1, tiempo2
    integer(kind=4)  ::  unit, ndum, unitmaster, conta, III, JJJ, KKK, pozi, i1t, j1t, k1t
    character(LEN=BUFSIZE)  ::  whoami, whoamishort
    logical :: ok, existe, wrotemaster, found
    integer(kind=8)  :: memo, ntini, ntfin
    character(LEN=BUFSIZE) :: buff, path, buff2
#ifdef CompileWithMPI
    integer(kind=MPI_OFFSET_KIND) disp
    integer(kind=4)  ::  ierr
#endif
    logical :: Esborde
    integer(kind=4)  ::  imed, imed1, imed2, imed3, imed4, medium
    integer(kind=4)  ::  thefile !for file management
!for dft
    REAL(KIND=RKIND), allocatable, dimension(:) :: signal, fqPos
    REAL(KIND=RKIND_tiempo), allocatable, dimension(:) :: samplingtime
    complex(kind=CKIND), allocatable, dimension(:) :: fqValues
    integer(kind=4) :: timesteps, klk, fqlength
    integer :: my_iostat
!
    eps0 = eps00; mu0 = mu00; !chapuz para convertir la variables de paso en globales
!
      !!
    output => null()
#ifdef CompileWithMPI
    valores => null()
    newvalores => null()  !auxiliary for Bloque currents sync
#endif

    unitmaster = -1000 !!!no se bien. Lo pongo absurdo
    unit = 1000 !initial
    if (unit >= 2.0_RKIND**31.0_RKIND - 1.0_RKIND) then
      call stoponerror(layoutnumber, size, 'Excesive number of probes')
    end if
    !
    write (whoamishort, '(i5)') layoutnumber + 1
    write (whoami, '(a,i5,a,i5,a)') '(', layoutnumber + 1, '/', size, ') '

    call crea_gnuplot

    allocate (InvEps(0:sgg%NumMedia), InvMu(0:sgg%NumMedia))

    incident = .false.
    InvEps(0:sgg%NumMedia) = 1.0_RKIND/(Eps0*sgg%Med(0:sgg%NumMedia)%Epr)
    InvMu(0:sgg%NumMedia) = 1.0_RKIND/(Mu0*sgg%Med(0:sgg%NumMedia)%Mur)

    allocate (output(1:sgg%NumberRequest))
    output(1:sgg%NumberRequest)%Trancos = -1
    output(1:sgg%NumberRequest)%SaveAll = .false.
    output(1:sgg%NumberRequest)%TimesWritten = -1

    do ii = 1, sgg%NumberRequest
      call preprocess_observation(sgg%Observation(ii), output(ii), sgg%tiempo, finaltimestep, sgg%dt, saveall)
    end do

    do ii = 1, sgg%NumberRequest
      allocate (output(ii)%item(1:sgg%Observation(ii)%nP))
#ifdef CompileWithMPI
      do i = 1, sgg%Observation(ii)%nP
        output(ii)%item(i)%ZIorig = sgg%Observation(ii)%P(i)%ZI
        output(ii)%item(i)%ZEorig = sgg%Observation(ii)%P(i)%ZE
      end do
#endif
      output(ii)%TimesWritten = 0 !for volumic probes
    end do

   do ii = 1, sgg%NumberRequest
      do i = 1, sgg%Observation(ii)%nP
         call eliminate_unnecesary_observation_points(sgg%Observation(ii)%P(i), output(ii)%item(i), &
            sgg%Sweep, sgg%SINPMLSweep, sgg%Observation(ii)%P(1)%ZI, sgg%Observation(ii)%P(1)%ZE, layoutnumber, size)
      end do
   end do

      !
      ThereAreObservation = .false.
      ThereAreFarFields = .false.
      do ii = 1, sgg%NumberRequest
        do i = 1, sgg%Observation(ii)%nP
          field = sgg%observation(ii)%P(i)%what
          if (field /= nothing) ThereAreObservation = .true.
        end do
      end do
#ifdef CompileWithMTLN
      block
        type(mtln_solver_t), pointer :: mtln_solver
        integer :: i, j
        mtln_solver => GetSolverPtr()
        do i = 1, ubound(mtln_solver%bundles, 1)
          if (ubound(mtln_solver%bundles(i)%probes, 1) /= 0) then
            do j = 1, ubound(mtln_solver%bundles(i)%probes, 1)
              if (mtln_solver%bundles(i)%probes(j)%in_layer) ThereAreObservation = .true.
            end do
          end if
        end do
      end block
#endif
      !
      memo = 0
      !
      IF (ThereAreObservation) then
#ifdef CompileWithMPI
        allocate (valores(0:BuffObse), newvalores(0:BuffObse))
        valores = 0.0_RKIND
        newvalores = 0.0_RKIND
#endif
        if (sgg%NumPlaneWaves >= 1) incident = .true. !150419 creo que no se ha dado nunca >1 porque los RC tocan el numero de modos pero solo hay una planewave
        !Write also the incident fields in case there are plane waves (useful in SE calculations)

        open (119, file=trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt')
        write (119, '(a)') '!END'
        close (119, status='delete')
        my_iostat = 0
9138    if (my_iostat /= 0) write (*, fmt='(a)', advance='no') '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9138,'.',layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
         open (19,file=trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt',err=9138,iostat=my_iostat,status='new',action='write')

        if ((trim(adjustl(wiresflavor)) == 'holland') .or. &
            (trim(adjustl(wiresflavor)) == 'transition')) then
          if (Therearewires) Hwireslocal => GetHwires()
        end if
#ifdef CompileWithBerengerWires
        if (trim(adjustl(wiresflavor)) == 'berenger') then
          if (Therearewires) Hwireslocal_Berenger => GetHwires_Berenger()
        end if
#endif
#ifdef CompileWithSlantedWires
        if ((trim(adjustl(wiresflavor)) == 'slanted') .or. (trim(adjustl(wiresflavor)) == 'semistructured')) then
          if (Therearewires) Hwireslocal_Slanted => GetHwires_Slanted()
        end if
#endif

         !!!!!!!!!Comun a todas las sondas freqdomain
        do ii = 1, sgg%NumberRequest
          if (SGG%Observation(ii)%FreqDomain) then
            !
            !Count the frequencies that are going to be registered
            !
            output(ii)%InitialFreq = sgg%observation(ii)%InitialFreq
            output(ii)%FinalFreq = sgg%observation(ii)%FinalFreq
            output(ii)%FreqStep = sgg%observation(ii)%FreqStep
            !
            if (SGG%Observation(ii)%FreqStep /= 0) then
              output(ii)%NumFreqs = int(abs(SGG%Observation(ii)%FinalFreq - SGG%Observation(ii)%InitialFreq)/SGG%Observation(ii)%FreqStep) + 1
            else
              output(ii)%NumFreqs = 1 !default
            end if
            if ((output(ii)%NumFreqs < 0)) then
              Buff = 'Freq. range for Freq. probes invalid'
              call stoponerror(layoutnumber, size, Buff)
            end if
            if ((output(ii)%NumFreqs > 100000)) then
              Buff = 'Too many Freqs requested (>100000)'
              call stoponerror(layoutnumber, size, Buff)
            end if
            
            ALLOCATE (output(ii)%Freq(1:output(ii)%NumFreqs), &
                      output(ii)%auxExp_E(1:output(ii)%NumFreqs), &
                      output(ii)%auxExp_H(1:output(ii)%NumFreqs))
!
            pozi = index(sgg%observation(ii)%outputrequest, '_log_')
            if (pozi == 0) then

              do iff1 = 1, output(ii)%NumFreqs
                output(ii)%Freq(iff1) = output(ii)%InitialFreq + (iff1 - 1)*output(ii)%FreqStep

              END DO
            else !logaritmico
              output(ii)%InitialFreq = log10(output(ii)%InitialFreq)
              output(ii)%FinalFreq = log10(output(ii)%FinalFreq)
              output(ii)%FreqStep = abs((output(ii)%InitialFreq - output(ii)%FinalFreq)/(output(ii)%NumFreqs))

              do iff1 = 1, output(ii)%NumFreqs
                output(ii)%Freq(iff1) = 10.0_RKIND**(output(ii)%InitialFreq + (iff1 - 1)*output(ii)%FreqStep)

              END DO
            end if
            !!!!!!!!!!!!!!!!
            errnofile = .FALSE.
            !!! Revisa si existe transfer fucntion y almacena la info en dos arrays
            IF (SGG%Observation(ii)%Transfer) then
              ALLOCATE (output(ii)%dftEntrada(1:output(ii)%NumFreqs))
              output(ii)%dftEntrada = 0.0_RKIND
!
              INQUIRE (file=trim(adjustl(sgg%observation(ii)%FileNormalize)), EXIST=errnofile)
              IF (.NOT. errnofile) THEN
                Buff = trim(adjustl(sgg%observation(ii)%FileNormalize))//' NORMALIZATION FILE DOES NOT EXIST'
                CALL STOPONERROR(layoutnumber, size, Buff)
              END IF
!precuenta
              timesteps = 0
              OPEN (15, file=trim(adjustl(sgg%observation(ii)%FileNormalize)))
              do
                READ (15, *, end=998) tiempo1, field1
                timesteps = timesteps + 1
                continue
              end do
998           CLOSE (15)
              allocate (samplingTime(1:timesteps))
              allocate (signal(1:timesteps))
              signal = 0.0_RKIND
              samplingTime = 0.0_RKIND_tiempo
!
              !read the normalization file and find its DFT
              OPEN (15, file=trim(adjustl(sgg%observation(ii)%FileNormalize)))
              DO klk = 1, timesteps
                READ (15, *) samplingTime(klk), signal(klk)
              END DO
              CLOSE (15)
!
              !niapa quitar 200120 ojooo
              if (niapapostprocess) then
                print *, 'Correcting in observation ', timesteps, trim(adjustl(sgg%observation(ii)%FileNormalize))
                DO klk = 1, timesteps
                  samplingTime(klk) = real(klk*sgg%dt, RKIND_tiempo)
                end do
              end if
              !fin niapa

              fqLength = output(ii)%NumFreqs
              allocate (fqValues(1:fqLength), fqPos(1:fqLength))
              fqValues(1:fqLength) = 0.0_RKIND
              fqPos(1:fqLength) = output(ii)%Freq(1:fqLength)
              call dtft(fqValues, fqPos, fqLength, samplingTime, signal, timesteps)
              output(ii)%dftEntrada = fqValues
              deallocate (samplingTime, signal, fqValues, fqPos)
            END IF !DEL TRANSFER
          end if !del freqdomain
        end do
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do ii = 1, sgg%NumberRequest
          wrotemaster = .false.
          loop_ob: do i = 1, sgg%Observation(ii)%nP
            I1 = sgg%observation(ii)%P(i)%XI
            J1 = sgg%observation(ii)%P(i)%YI
            K1 = sgg%observation(ii)%P(i)%ZI
            I2 = sgg%observation(ii)%P(i)%XE
            J2 = sgg%observation(ii)%P(i)%YE
            K2 = sgg%observation(ii)%P(i)%ZE
            NO = sgg%observation(ii)%P(i)%NODE
            write (chari, '(i7)') i1
            write (charj, '(i7)') j1
            write (chark, '(i7)') k1
            field = sgg%observation(ii)%P(i)%what
            select case (field)
              !
            case (iEx, iEy, iEz, iVx, iVy, iVz, iJx, iJy, iJz, iQx, iQy, iQz, iHx, iHy, iHz, lineIntegral)
              !
              if (((field == iEx) .or. (field == iEy) .or. (field == iEz) .or. &
                   (field == iHx) .or. (field == iHy) .or. (field == iHz)) .and. &
                  (sgg%NumPlaneWaves >= 1)) then
                output(ii)%item(i)%columnas = 2
              elseif ((field == iJx) .or. (field == iJy) .or. (field == iJz)) then
                output(ii)%item(i)%columnas = 5 ! corriente -e*dl vplus vminus vplus-vminus
              elseif ((field == iVx) .or. (field == iVy) .or. (field == iVz)) then
                output(ii)%item(i)%columnas = 1
              else
                output(ii)%item(i)%columnas = 1
              end if
              !mpdir 190319       !desrotacion para que los nombres sean correctos
              if (mpidir == 3) then
                extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
                select case (field)
                case (iEx)
                  prefix_field = prefix(iEx)
                case (iEy)
                  prefix_field = prefix(iEy)
                case (iEz)
                  prefix_field = prefix(iEz)
                case (iJx)
                  prefix_field = prefix(iJx)
                case (iJy)
                  prefix_field = prefix(iJy)
                case (iJz)
                  prefix_field = prefix(iJz)
                case (iQx)
                  prefix_field = prefix(iQx)
                case (iQy)
                  prefix_field = prefix(iQy)
                case (iQz)
                  prefix_field = prefix(iQz)
                case (iVx)
                  prefix_field = prefix(iVx)
                case (iVy)
                  prefix_field = prefix(iVy)
                case (iVz)
                  prefix_field = prefix(iVz)
                case (iHx)
                  prefix_field = prefix(iHx)
                case (iHy)
                  prefix_field = prefix(iHy)
                case (iHz)
                  prefix_field = prefix(iHz)
                case default
                  prefix_field = prefix(field)
                end select
              elseif (mpidir == 2) then
                extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
                select case (field)
                case (iEx)
                  prefix_field = prefix(iEz)
                case (iEy)
                  prefix_field = prefix(iEx)
                case (iEz)
                  prefix_field = prefix(iEy)
                case (iJx)
                  prefix_field = prefix(iJz)
                case (iJy)
                  prefix_field = prefix(iJx)
                case (iJz)
                  prefix_field = prefix(iJy)
                case (iQx)
                  prefix_field = prefix(iQz)
                case (iQy)
                  prefix_field = prefix(iQx)
                case (iQz)
                  prefix_field = prefix(iQy)
                case (iVx)
                  prefix_field = prefix(iVz)
                case (iVy)
                  prefix_field = prefix(iVx)
                case (iVz)
                  prefix_field = prefix(iVy)
                case (iHx)
                  prefix_field = prefix(iHz)
                case (iHy)
                  prefix_field = prefix(iHx)
                case (iHz)
                  prefix_field = prefix(iHy)
                case default
                  prefix_field = prefix(field)
                end select
              elseif (mpidir == 1) then
                extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
                select case (field)
                case (iEx)
                  prefix_field = prefix(iEy)
                case (iEy)
                  prefix_field = prefix(iEz)
                case (iEz)
                  prefix_field = prefix(iEx)
                case (iJx)
                  prefix_field = prefix(iJy)
                case (iJy)
                  prefix_field = prefix(iJz)
                case (iJz)
                  prefix_field = prefix(iJx)
                case (iQx)
                  prefix_field = prefix(iQy)
                case (iQy)
                  prefix_field = prefix(iQz)
                case (iQz)
                  prefix_field = prefix(iQx)
                case (iVx)
                  prefix_field = prefix(iVy)
                case (iVy)
                  prefix_field = prefix(iVz)
                case (iVz)
                  prefix_field = prefix(iVx)
                case (iHx)
                  prefix_field = prefix(iHy)
                case (iHy)
                  prefix_field = prefix(iHz)
                case (iHz)
                  prefix_field = prefix(iHx)
                case default
                  prefix_field = prefix(field)
                end select
              else
                call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
              end if
              !
              if ((field == iJx) .or. (field == iJy) .or. (field == iJz)) then
                write (charNO, '(i7)') NO
                !append the number of the segment
                extpoint = trim(adjustl(extpoint))//'_s'//trim(adjustl(charNO))
              end if
              if ((field == iQx) .or. (field == iQy) .or. (field == iQz)) then
                write (charNO, '(i7)') NO
                !append the number of the segment
                extpoint = trim(adjustl(extpoint))//'_s'//trim(adjustl(charNO))
              end if
              ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
              !do not use layername since no two observations from different layers will overlap
              output(ii)%item(i)%path = &
                trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'
              !
              unit = unit + 1
              if (unit >= 2.0_RKIND**31.0_RKIND - 1.0_RKIND) then
                call stoponerror(layoutnumber, size, 'Excesive number of probes')
              end if
              output(ii)%item(i)%unit = unit
              !

                  !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
              call checkduplicatenames
                  !!!!!!

              my_iostat = 0
9235          if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9235,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
              write (19, '(a)', err=9235, iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))
              !
              memo = memo + rkind*BuffObse
              if (memo > MaxMemoryProbes) then
                call stoponerror(layoutnumber, size, 'Recompile: excesive memory for probes.'// &
                &                                   'Increase MaxMemoryProbes')
              end if
              allocate (output(ii)%item(i)%valor(0:BuffObse))
              output(ii)%item(i)%valor(0:BuffObse) = 0.0_RKIND

              if (field == iQx .or. field == iQy .or. field == iQz) then
                found = .false.
                do n = 1, HWireslocal%NumCurrentSegments
                  if ((HWireslocal%CurrentSegment(n)%origindex == no) .and. &
                      (HWireslocal%CurrentSegment(n)%i == i1) .and. &
                      (HWireslocal%CurrentSegment(n)%j == j1) .and. &
                      (HWireslocal%CurrentSegment(n)%k == k1) .and. &
                      (HWireslocal%CurrentSegment(n)%tipofield*10000 == field)) then
                    output(ii)%item(i)%segmento => HWireslocal%CurrentSegment(n)
                    if (output(ii)%item(i)%segmento%orientadoalreves) output(ii)%item(i)%valorsigno = -1
                    found = .true.
                  end if
                end do
                if ((.not. found) .and. ((field == iQx) .or. (field == iQy) .or. (field == iQz))) then
                  sgg%Observation(ii)%P(i)%What = nothing
                  write (buff, '(a,4i7,a)') 'ERROR: CHARGE probe ', no, i1, j1, k1, ' DOES NOT EXIST'
                  CALL WarnErrReport(buff, .true.)
                end if

              end if

              if ((trim(adjustl(wiresflavor)) == 'holland') .or. &
                  (trim(adjustl(wiresflavor)) == 'transition')) then
                found = .false.
                if ((Therearewires) .and. ((field == iJx) .or. (field == iJy) .or. (field == iJz))) then

                  memo = memo + 3*4*BuffObse
                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'Recompile: excesive memory for probes.'// &
                    &                                   'Increase MaxMemoryProbes')
                  end if
                  allocate ( &
                    output(ii)%item(i)%valor2(0:BuffObse), &
                    output(ii)%item(i)%valor3(0:BuffObse), &
                    output(ii)%item(i)%valor4(0:BuffObse), &
                    output(ii)%item(i)%valor5(0:BuffObse))
                  output(ii)%item(i)%valor2(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor3(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor4(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor5(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valorsigno = +1
                  !en caso de hilos se necesitan
                  !parsea los hilos
                  found = .false.
                  output(ii)%item(i)%segmento => HWireslocal%NullSegment !por si no encuentra el segmento por estar repetido !bug serio que impide discernir para segmentos
                  do n = 1, HWireslocal%NumCurrentSegments
                    if ((HWireslocal%CurrentSegment(n)%origindex == no) .and. & !el nodo tambien debe coincidir !bug serio que impide discernir para segmentos paralelos 12/09/13 cazado gracias a Model_unidos.nfde
                        (HWireslocal%CurrentSegment(n)%i == i1) .and. &
                        (HWireslocal%CurrentSegment(n)%j == j1) .and. &
                        (HWireslocal%CurrentSegment(n)%k == k1) .and. &
                        (HWireslocal%CurrentSegment(n)%tipofield*10 == field)) then
                      !I have chosen IJx=10 IEx, etc. !do not change
                      output(ii)%item(i)%segmento => HWireslocal%CurrentSegment(n)
                      if (output(ii)%item(i)%segmento%orientadoalreves) output(ii)%item(i)%valorsigno = -1
                      found = .true.
                    end if
                  end do
                  if (.not. found) then
                    !busca por si fuera un multirabito
                    buscarabono: do iwi = 1, Hwireslocal%NumDifferentWires
                      do iwj = 1, sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%numsegmentos
                                 if ((no==sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex).and.sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                          no2 = sgg%Med(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE
                          do n = 1, HWireslocal%NumCurrentSegments
                            if (HWireslocal%CurrentSegment(n)%origindex == no2) then !el nodo tambien debe coincidir aunque no necesariamente coordenadas ni campo porque se ha cortado el rabo original
                              !I have chosen IJx=10 IEx, etc. !do not change
                              output(ii)%item(i)%segmento => HWireslocal%CurrentSegment(n)
                              if (output(ii)%item(i)%segmento%orientadoalreves) output(ii)%item(i)%valorsigno = -1
                              found = .true.
                            end if
                          end do
                          found = .true.
                          exit buscarabono
                        end if
                      end do
                    end do buscarabono
                  end if
                end if
                if ((.not. found) .and. ((((field == iJx) .or. (field == iJy) .or. (field == iJz))))) then
                  sgg%Observation(ii)%P(i)%What = nothing
                  !ojoo 010423 para debugeo lbb1
                  write (buff, '(a,4i7,a)') 'ERROR: WIRE probe ', no, i1, j1, k1, ' DOES NOT EXIST'
                  CALL WarnErrReport(buff, .true.)
                end if
              end if

#ifdef CompileWithBerengerWires
              if (trim(adjustl(wiresflavor)) == 'berenger') then
                found = .false.
                if ((Therearewires) .and. ((field == iJx) .or. (field == iJy) .or. (field == iJz))) then

                  memo = memo + 3*4*BuffObse
                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'Recompile: excesive memory for probes.'// &
                    &                                   'Increase MaxMemoryProbes')
                  end if
                  allocate ( &
                    output(ii)%item(i)%valor2(0:BuffObse), &
                    output(ii)%item(i)%valor3(0:BuffObse), &
                    output(ii)%item(i)%valor4(0:BuffObse), &
                    output(ii)%item(i)%valor5(0:BuffObse))
                  output(ii)%item(i)%valor2(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor3(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor4(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor5(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valorsigno = +1
                  !en caso de hilos se necesitan
                  !parsea los hilos
                  found = .false.
                  do n = 1, Hwireslocal_Berenger%NumSegments
                    if (Hwireslocal_Berenger%Segments(n)%IndexSegment == no) then !solo miro el nodo. dama ya corrige esto la observation de Berenger
                      !bug dectectado por OLD 311019 con probe_issue.nfde
!                           if ((Hwireslocal_Berenger%Segments(n)%IndexSegment==no).and. & !el nodo tambien debe coincidir !bug serio que impide discernir para segmentos paralelos 12/09/13 cazado gracias a Model_unidos.nfde
!                           (Hwireslocal_Berenger%Segments(n)%ii==i1).and. &
!                           (Hwireslocal_Berenger%Segments(n)%ji==j1).and. &
!                           (Hwireslocal_Berenger%Segments(n)%ki==k1).and. &
!                           (Hwireslocal_Berenger%Segments(n)%orient*10==field))  then
                      !I have chosen IJx=10 IEx, etc. !do not change
                      output(ii)%item(i)%segmento_Berenger => Hwireslocal_Berenger%Segments(n)
                      if (output(ii)%item(i)%segmento_Berenger%orientadoalreves) output(ii)%item(i)%valorsigno = -1
                      found = .true.
                    end if
                  end do
                end if
                if ((.not. found) .and. ((((field == iJx) .or. (field == iJy) .or. (field == iJz))))) then
                  sgg%Observation(ii)%P(i)%What = nothing
                  write (buff, '(a,4i7,a)') 'ERROR: WIRE probe ', no, i1, j1, k1, ' DOES NOT EXIST'
                  CALL WarnErrReport(buff, .TRUE.)
                end if
              end if
#endif
#ifdef CompileWithSlantedWires
              if ((trim(adjustl(wiresflavor)) == 'slanted') .or. (trim(adjustl(wiresflavor)) == 'semistructured')) then
                found = .false.
                if ((Therearewires) .and. ((field == iJx) .or. (field == iJy) .or. (field == iJz))) then

                  memo = memo + 3*4*BuffObse
                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'Recompile: excesive memory for probes.'// &
                    &                                   'Increase MaxMemoryProbes')
                  end if
                  allocate ( &
                    output(ii)%item(i)%valor2(0:BuffObse), &
                    output(ii)%item(i)%valor3(0:BuffObse), &
                    output(ii)%item(i)%valor4(0:BuffObse), &
                    output(ii)%item(i)%valor5(0:BuffObse))
                  output(ii)%item(i)%valor2(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor3(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor4(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valor5(0:BuffObse) = 0.0_RKIND
                  output(ii)%item(i)%valorsigno = +1
                  !en caso de hilos se necesitan
                  !parsea los hilos
                  found = .false.
                  do n = 1, Hwireslocal_Slanted%NumSegments
                    if (Hwireslocal_Slanted%Segments(n)%ptr%Index == no) then
                      !I have chosen IJx=10 IEx, etc. !do not change
                      output(ii)%item(i)%segmento_Slanted => Hwireslocal_Slanted%Segments(n)%ptr
                      found = .true.
                    end if
                  end do
                end if
                !010423  creo que si no lo encuentra es porque el indice es el exterior bug lbb1 epg 0323
                if (.not. found) then
                  do n = 1, Hwireslocal_Slanted%NumSegments
                    if (Hwireslocal_Slanted%Segments(n)%ptr%elotroindice == no) then
                      output(ii)%item(i)%segmento_Slanted => Hwireslocal_Slanted%Segments(n)%ptr
                      found = .true.
                    end if
                  end do
                end if
                if ((.not. found) .and. ((((field == iJx) .or. (field == iJy) .or. (field == iJz))))) then
                  sgg%Observation(ii)%P(i)%What = nothing

                  !ojoo 010423 para debugeo lbb1
                  write (buff, '(a,4i7,a)') 'ERROR: WIRE probe ', no, i1, j1, k1, ' DOES NOT EXIST'
                  CALL WarnErrReport(buff, .true.)
                end if
              end if
#endif

                  !!!!!!!!!!!!!!!
              !erase pre-existing data unless this is a resuming simulation
              if (.not. resume) then
                !
                if (singlefilewrite) then
                  if (.not. wrotemaster) then
                    wrotemaster = .true.
                    unitmaster = output(ii)%item(i)%unit
                    output(ii)%item(i)%unitmaster = unitmaster
           open (unitmaster, recl=1000, file=trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin')
                    write (unitmaster, *) '!END'
                    close (unitmaster, status='delete')

                    my_iostat = 0
9238                if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9238,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin'
                           open (unitmaster,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin',form='unformatted',err=9238,iostat=my_iostat,status='new',action='write')
                  else
                    output(ii)%item(i)%unitmaster = unitmaster
                  end if
                else
                  !
                  open (output(ii)%item(i)%unit, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)))
                  write (output(ii)%item(i)%unit, *) '!END'
                  close (output(ii)%item(i)%unit, status='delete')
                  my_iostat = 0
8766              if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8766,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
                        open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),err=8766,iostat=my_iostat,status='new',action='write')
                  write (output(ii)%item(i)%unit, '(a)') trim(adjustl(' t'//'              '// &
                                         trim(adjustl(output(ii)%item(i)%path))//'       '//trim(adjustl(suffix(field, incident)))))
                end if
                !
                !                        write(output(ii)%item(i)%unit,'(5a)') ' t','              ', &
                     !!                            trim(adjustl(prefix(field)))//trim(adjustl(extpoint)),'   ',trim(adjustl(suffix(field,incident)))
                !                            trim(adjustl(output(ii)%item(i)%path)),'       ',trim(adjustl(suffix(field,incident)))
                !
              else !wipe out duplicate data after non synchronous data and field resuming
                !
                if (singlefilewrite) then
                  if (.not. wrotemaster) then
                    wrotemaster = .true.
                    unitmaster = output(ii)%item(i)%unit
                    output(ii)%item(i)%unitmaster = unitmaster
open (unitmaster, recl=1000, file=trim(adjustl(output(ii)%item(i)%path))//'_master'//'_'//trim(adjustl(whoamishort))//'_master.bin')
                    write (unitmaster, *) '!END'
                    close (unitmaster, status='delete')
                    my_iostat = 0
9239                if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9239,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin'
                           open (unitmaster,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin',form='unformatted',err=9239,iostat=my_iostat,status='new',action='write')
                  else
                    output(ii)%item(i)%unitmaster = unitmaster
                  end if
                else
                  inquire (file=trim(adjustl(output(ii)%item(i)%path)), exist=existe)
                  if (.not. existe) then
    call stoponerror(layoutnumber, size, 'Data files for resuming non existent (Ex, etc.) '//trim(adjustl(output(ii)%item(i)%path)))
                  end if
                  !
                  open (output(ii)%item(i)%unit, recl=1000, access='sequential', file=trim(adjustl(output(ii)%item(i)%path)))
                  read (output(ii)%item(i)%unit, '(a)') adum !first line contains characters
                  cutting: do
                    read (output(ii)%item(i)%unit, *, end=678) at
                    if (rdum > lastexecutedtime) then
                     print '(i4,a,a,2e19.9e3)', quienmpi, 'Cutting 1 ', trim(adjustl(output(ii)%item(i)%path)), at, lastexecutedtime
                      backspace (output(ii)%item(i)%unit)
                      endfile (output(ii)%item(i)%unit)
                      exit cutting
                    end if
                  end do cutting
678               continue
                  close (output(ii)%item(i)%unit)
                  open (output(ii)%item(i)%unit, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)), position='append')
                end if

              end if
            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
              output(ii)%item(i)%columnas = 1
              write (chari2, '(i7)') i2
              write (charj2, '(i7)') j2
              write (chark2, '(i7)') k2
              !mpidir 190319     !desrotacion para que los nombres sean correctos
              if (mpidir == 3) then
                extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
                           trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                select case (field)
                case (iBloqueJx)
                  prefix_field = prefix(iBloqueJx)
                case (iBloqueJy)
                  prefix_field = prefix(iBloqueJy)
                case (iBloqueJz)
                  prefix_field = prefix(iBloqueJz)
                case (iBloqueMx)
                  prefix_field = prefix(iBloqueMx)
                case (iBloqueMy)
                  prefix_field = prefix(iBloqueMy)
                case (iBloqueMz)
                  prefix_field = prefix(iBloqueMz)
                end select
              elseif (mpidir == 2) then
                extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
                           trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                select case (field)
                case (iBloqueJx)
                  prefix_field = prefix(iBloqueJz)
                case (iBloqueJy)
                  prefix_field = prefix(iBloqueJx)
                case (iBloqueJz)
                  prefix_field = prefix(iBloqueJy)
                case (iBloqueMx)
                  prefix_field = prefix(iBloqueMz)
                case (iBloqueMy)
                  prefix_field = prefix(iBloqueMx)
                case (iBloqueMz)
                  prefix_field = prefix(iBloqueMy)
                end select
              elseif (mpidir == 1) then
                extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
                           trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                select case (field)
                case (iBloqueJx)
                  prefix_field = prefix(iBloqueJy)
                case (iBloqueJy)
                  prefix_field = prefix(iBloqueJz)
                case (iBloqueJz)
                  prefix_field = prefix(iBloqueJx)
                case (iBloqueMx)
                  prefix_field = prefix(iBloqueMy)
                case (iBloqueMy)
                  prefix_field = prefix(iBloqueMz)
                case (iBloqueMz)
                  prefix_field = prefix(iBloqueMx)
                end select
              else
                call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
              end if
              !
              ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
              !
              output(ii)%item(i)%path = &
                trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'

              !
              unit = unit + 1
              if (unit >= 2.0_RKIND**31.0_RKIND - 1.0_RKIND) then
                call stoponerror(layoutnumber, size, 'Excesive number of probes')
              end if
              output(ii)%item(i)%unit = unit
              !
                  !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
              call checkduplicatenames
                  !!!!!!

              memo = memo + rkind*BuffObse
              if (memo > MaxMemoryProbes) then
                call stoponerror(layoutnumber, size, 'ERROR: Recompile: excesive memory for the probes.'// &
                &                                   'Recompile increasing MaxMemoryProbes')
              end if
              allocate (output(ii)%item(i)%valor(0:BuffObse))
              output(ii)%item(i)%valor(0:BuffObse) = 0.0_RKIND
              !readjust correctly the calculation region (for currents crossing the MPI region)
              select case (field)
              case (iBloqueJx, iBloqueJy)
                sgg%observation(ii)%P(i)%ZI = max(sgg%Sweep(fieldo(field, 'Z'))%ZI + 1, sgg%observation(ii)%P(i)%ZI)
                sgg%observation(ii)%P(i)%ZE = min(sgg%Sweep(fieldo(field, 'Z'))%ZE, sgg%observation(ii)%P(i)%ZE)
              case (iBloqueMx, iBloqueMy, iBloqueJz, iBloqueMz)
                sgg%observation(ii)%P(i)%ZI = max(sgg%Sweep(fieldo(field, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
                sgg%observation(ii)%P(i)%ZE = min(sgg%Sweep(fieldo(field, 'Z'))%ZE, sgg%observation(ii)%P(i)%ZE)
              end select
              !
#ifdef CompileWithMPI
              if ((layoutnumber == output(ii)%item(i)%MPIRoot) .or. &
                  (field == iBloqueJz) .or. (field == iBloqueMz)) then !only the master
#endif
                my_iostat = 0
9837            if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9837,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                write (19, '(a)', err=9837, iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))
                !erase pre-existing data unless this is a resuming simulation
                if (.not. resume) then
                  open (output(ii)%item(i)%unit, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)))
                  write (output(ii)%item(i)%unit, *) '!END'
                  close (output(ii)%item(i)%unit, status='DELETE')
                  my_iostat = 0
9838              if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9838,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
                     open (output(ii)%item(i)%unit,recl=1000,file= trim(adjustl(output(ii)%item(i)%path)),err=9838,iostat=my_iostat,status='new',action='write')
                  write (output(ii)%item(i)%unit, '(a)') trim(adjustl(' t'//'              '// &
                                                                      trim(adjustl(output(ii)%item(i)%path))//'       '// &
                                                                      trim(adjustl(suffix(field, incident)))))
                  !ojo a esto no le he corregido lo del mpidir de a 190319 porque estaba comentado
                  !
                  !                            write(output(ii)%item(i)%unit,'(5a)') ' t','              ', &
                  !                               !trim(adjustl(prefix(field)))//trim(adjustl(extpoint)),'   ',&
                  !                                trim(adjustl(output(ii)%item(i)%path)),'       ',&
                  !                                trim(adjustl(suffix(field,incident)))
                  !
                  !wipe out duplicate data after non synchronous data and field resuming
                else
                  inquire (file=trim(adjustl(output(ii)%item(i)%path)), exist=existe)
                  if (.not. existe) then
call stoponerror(layoutnumber, size, 'Data files for resuming non existent (Bloque, etc.) '//trim(adjustl(output(ii)%item(i)%path)))
                  end if
                  open (output(ii)%item(i)%unit, recl=1000, access='sequential', &
                        file=trim(adjustl(output(ii)%item(i)%path)))
                  read (output(ii)%item(i)%unit, '(a)') adum !first line contains characters
                  cutting2: do
                    read (output(ii)%item(i)%unit, *, end=679) at
                    if (at > lastexecutedtime) then
                     print '(i4,a,a,2e19.9e3)', quienmpi, 'Cutting 2 ', trim(adjustl(output(ii)%item(i)%path)), at, lastexecutedtime
                      backspace (output(ii)%item(i)%unit)
                      endfile (output(ii)%item(i)%unit)
                      exit cutting2
                    end if
                  end do cutting2
679               continue
                  close (output(ii)%item(i)%unit)
                  open (output(ii)%item(i)%unit, recl=1000, file= &
                        trim(adjustl(output(ii)%item(i)%path)), position='append')
                end if
#ifdef CompileWithMPI
              end if
#endif
              !voluminc Bloque current probes along edges (wires, surfaces
            case (iCur, iCurX, iCurY, iCurZ, mapvtk)
              if (sgg%Observation(ii)%Volumic) then !they are necssaryly
                if (sgg%Observation(ii)%nP /= 1) then
                  call stoponerror(layoutnumber, size, 'ERROR! More than a volumic probe per group')
                end if
                !readjut correctly the calculation region
                !de momento sere conservador 20/2/14 por lo que truene el MPI luego quitare el -1 si acaso !!!! !a priori puedo necesitar el HZ(alloc+1) para calcular las Bloque currents pero de momento me estoy quieto
                sgg%observation(ii)%P(i)%ZI = max(sgg%Sweep(iEx)%ZI, sgg%observation(ii)%P(i)%ZI) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
                sgg%observation(ii)%P(i)%ZE = min(sgg%Sweep(iEx)%ZE, sgg%observation(ii)%P(i)%ZE) !ojo estaba sweep(iEz) para ser conservador...puede dar problemas!! 03/07/15
                !solo acepto que P(1:1) !!!
                write (chari, '(i7)') sgg%observation(ii)%P(i)%XI
                write (charj, '(i7)') sgg%observation(ii)%P(i)%YI
                write (chark, '(i7)') sgg%observation(ii)%P(i)%ZI
                write (chari2, '(i7)') sgg%observation(ii)%P(i)%XE
                write (charj2, '(i7)') sgg%observation(ii)%P(i)%YE
                write (chark2, '(i7)') sgg%observation(ii)%P(i)%ZE
                !mpidir 190319     !desrotacion para que los nombres sean correctos
                if (mpidir == 3) then
                  extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
                             trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                  select case (field)
                  case (iCurX)
                    prefix_field = prefix(iCurX)
                  case (iCurY)
                    prefix_field = prefix(iCurY)
                  case (iCurZ)
                    prefix_field = prefix(iCurZ)
                  case default
                    prefix_field = prefix(field)
                  end select
                elseif (mpidir == 2) then
                  extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
                             trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                  select case (field)
                  case (iCurX)
                    prefix_field = prefix(iCurZ)
                  case (iCurY)
                    prefix_field = prefix(iCurX)
                  case (iCurZ)
                    prefix_field = prefix(iCurY)
                  case default
                    prefix_field = prefix(field)
                  end select
                elseif (mpidir == 1) then
                  extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
                             trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                  select case (field)
                  case (iCurX)
                    prefix_field = prefix(iCurY)
                  case (iCurY)
                    prefix_field = prefix(iCurZ)
                  case (iCurZ)
                    prefix_field = prefix(iCurX)
                  case default
                    prefix_field = prefix(field)
                  end select
                else
                  call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
                end if
                !
                ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                !cada mpi su nombre
                output(ii)%item(i)%path = trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.bin'

                !
                unit = unit + 1
                if (unit >= 2.0_RKIND**31.0_RKIND - 1.0_RKIND) then
                  call stoponerror(layoutnumber, size, 'Excesive number of probes')
                end if
                output(ii)%item(i)%unit = unit
                output(ii)%item(i)%columnas = 0 !numero de items a escribir

                     !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                call checkduplicatenames
                     !!!!!!
                !precontaje

                conta = 0
                do kkk = sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                  do jjj = sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                    do iii = sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                      if (field /= mapvtk) then
                        do Efield = iEx, iEz

                          if (isWithinBounds(Efield, iii, jjj, kkk)) then
                            if (isThinWire(Efield, iii, jjj, kkk)) then
                              conta = conta + 1
                            end if

                            if (.not. isMediaVacuum(Efield, iii, jjj, kkk) .and. &
                                .not. isSplitOrAdvanced(Efield, iii, jjj, kkk)) then
                              conta = conta + 1
                            end if
                          end if
                        end do

                      else !si es mapvtk
                        !si es mapvtk y si no es vacio
                        do Efield = iEx, iEz
                          call assignMedia(imed, imed1, imed2, imed3, imed4, Efield, iii, jjj, kkk)
                          call contabordes(sgg, imed, imed1, imed2, imed3, imed4, EsBorde, SINPML_fullsize, Efield, iii, jjj, kkk)
                          if (EsBorde) then
                            conta = conta + 1
                          end if
                        end do

                      end if
                    end do
                  end do
                end do

                     !!!
                if (field == mapvtk) then
                  INIT = .TRUE.; geom = .false.; asigna = .false.; magnetic = .false.; electric = .true.

  call nodalvtk(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,media%sggMtag,tag_numbers, &
                                init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)

        call wirebundlesvtk(sgg, init, geom, asigna, conta, i, ii, output, Ntimeforvolumic, wiresflavor, media%sggMtag, tag_numbers)
                end if
                     !!!
                do kkk = sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                  do jjj = sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                    do iii = sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                      if (field /= mapvtk) then
                        ! count PEC surfaces
                        do Hfield = iHx, iHz
                          if ((isPECorSurface(Hfield, iii, jjj, kkk) .or. field == blockCurrent(Hfield)) &
                              .and. isWithinBounds(Hfield, iii, jjj, kkk)) then
                            conta = conta + 1
                          end if
                        end do
                      else
                        ! si es mapvtk y si no es vacio
                        do Hfield = iHx, iHz
                          ! count media surfaces
                          if (.not. isMediaVacuum(Hfield, iii, jjj, kkk) .and. &
                              .not. isPML(Hfield, iii, jjj, kkk) .and. &
                              isWithinBounds(Hfield, iii, jjj, kkk)) then
                            conta = conta + 1
                          end if
                          ! count negative vacuum tag numbers
                          if (tag_numbers%getFaceTag(Hfield, iii, jjj, kkk) < 0 .and. &
                              (btest(iabs(tag_numbers%getFaceTag(Hfield, iii, jjj, kkk)), Hfield - 1)) .and. &
                              .not. isPML(Hfield, iii, jjj, kkk) .and. &
                              isWithinBounds(Hfield, iii, jjj, kkk)) then
                            conta = conta + 1
                          end if
                        end do
                      end if
                    end do
                  end do
                end do
                     !!!
                if (field == mapvtk) then
                  INIT = .false.; geom = .false.; asigna = .false.; magnetic = .true.; electric = .false.
  call nodalvtk(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,media%sggMtag,tag_numbers, &
                                init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)

                end if
                     !!!
                output(ii)%item(i)%columnas = conta
                !print *,' ---- ALLOC DE output(ii)%item(i)%columnas=conta ', II,I,CONTA
                !allocateo

                !
                IF (SGG%Observation(ii)%TimeDomain) THEN
                  !ojo por si algun dia esto molestara a Cray
                  !replico los ifs de transferencia y escritura
                  ntini = 0
                  ntfin = 0
                  first = .true.
                  do Ntime = initialtimestep, finaltimestep
                    at = sgg%tiempo(Ntime)
                    Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                    if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                      Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                        if (first) then
                          ntini = ntimeforvolumic
                          first = .false.
                        end if
                        ntfin = ntimeforvolumic
                      end if
                    end if
                  end do
                        !!!                            ntinI=0
                        !!!                            ntfin=min(int((finaltimestep*sgg%dt-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos), &
                        !!!                                      int((sgg%Observation(ii)%FinalTIME-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos))+1

                  memo = memo + RKIND*(ntfin - ntini)*output(ii)%item(i)%columnas + 16*output(ii)%item(i)%columnas ! 4 integers de 4 bytes

                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'ERROR: Recompile: excesive memory for the 3D probes.'// &
                    &                                   'Recompile increasing MaxMemoryProbes')
                  end if

                  allocate (output(ii)%item(i)%Serialized%valor(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  !almaceno tambien los vectores
                  allocate (output(ii)%item(i)%Serialized%valor_x(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  allocate (output(ii)%item(i)%Serialized%valor_y(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  allocate (output(ii)%item(i)%Serialized%valor_z(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  output(ii)%item(i)%Serialized%Valor = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_x = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_y = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_z = 0.0_RKIND
                  !electric

                  allocate (output(ii)%item(i)%Serialized%valorE(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  !almaceno tambien los vectores
                  allocate (output(ii)%item(i)%Serialized%valor_Ex(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  allocate (output(ii)%item(i)%Serialized%valor_Ey(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  allocate (output(ii)%item(i)%Serialized%valor_Ez(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  output(ii)%item(i)%Serialized%ValorE = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_Ex = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_Ey = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_Ez = 0.0_RKIND
                  !magnetic

                  allocate (output(ii)%item(i)%Serialized%valorH(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  !almaceno tambien los vectores
                  allocate (output(ii)%item(i)%Serialized%valor_Hx(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  allocate (output(ii)%item(i)%Serialized%valor_Hy(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  allocate (output(ii)%item(i)%Serialized%valor_Hz(ntinI:ntfin, 1:output(ii)%item(i)%columnas))
                  output(ii)%item(i)%Serialized%ValorH = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_Hx = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_Hy = 0.0_RKIND
                  output(ii)%item(i)%Serialized%Valor_Hz = 0.0_RKIND
                ELSEIF (SGG%Observation(ii)%FreqDomain) THEN
                  memo = memo + RKIND*output(ii)%NumFreqs*output(ii)%item(i)%columnas + 16*output(ii)%item(i)%columnas ! 4 integers de 4 bytes
                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'ERROR: Recompile: excesive memory for the probes.'// &
                    &                                   'Recompile increasing MaxMemoryProbes')
                  end if
                  !almaceno tambien los vectores
                  allocate (output(ii)%item(i)%Serialized%valorComplex_x(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                  allocate (output(ii)%item(i)%Serialized%valorComplex_y(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                  allocate (output(ii)%item(i)%Serialized%valorComplex_z(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes

                  output(ii)%item(i)%Serialized%ValorComplex_x = (0.0_RKIND, 0.0_RKIND)
                  output(ii)%item(i)%Serialized%ValorComplex_y = (0.0_RKIND, 0.0_RKIND)
                  output(ii)%item(i)%Serialized%ValorComplex_z = (0.0_RKIND, 0.0_RKIND)
                  !ELECTRIC

                  !almaceno tambien los vectores
                  allocate (output(ii)%item(i)%Serialized%valorComplex_Ex(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                  allocate (output(ii)%item(i)%Serialized%valorComplex_Ey(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                  allocate (output(ii)%item(i)%Serialized%valorComplex_Ez(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes

                  output(ii)%item(i)%Serialized%ValorComplex_Ex = (0.0_RKIND, 0.0_RKIND)
                  output(ii)%item(i)%Serialized%ValorComplex_Ey = (0.0_RKIND, 0.0_RKIND)
                  output(ii)%item(i)%Serialized%ValorComplex_Ez = (0.0_RKIND, 0.0_RKIND)
                  !MAGNETIC

                  !almaceno tambien los vectores
                  allocate (output(ii)%item(i)%Serialized%valorComplex_Hx(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                  allocate (output(ii)%item(i)%Serialized%valorComplex_Hy(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes
                  allocate (output(ii)%item(i)%Serialized%valorComplex_Hz(1:output(ii)%NumFreqs, 1:output(ii)%item(i)%columnas)) !dos posibles componentes

                  output(ii)%item(i)%Serialized%ValorComplex_Hx = (0.0_RKIND, 0.0_RKIND)
                  output(ii)%item(i)%Serialized%ValorComplex_Hy = (0.0_RKIND, 0.0_RKIND)
                  output(ii)%item(i)%Serialized%ValorComplex_Hz = (0.0_RKIND, 0.0_RKIND)
                end if

                allocate (output(ii)%item(i)%Serialized%eI(1:output(ii)%item(i)%columnas)) !he aprovechado la variable columnas!!jeje...comentario hecho el 120120 por bug vtk OLD
                allocate (output(ii)%item(i)%Serialized%eJ(1:output(ii)%item(i)%columnas))
                allocate (output(ii)%item(i)%Serialized%eK(1:output(ii)%item(i)%columnas))
                allocate (output(ii)%item(i)%Serialized%currentType(1:output(ii)%item(i)%columnas))
                allocate (output(ii)%item(i)%Serialized%sggMtag(1:output(ii)%item(i)%columnas))

                !relleno info geometrica edge
                conta = 0
                do kkk = sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                  do jjj = sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                    do iii = sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                      if (field /= mapvtk) then
                        do Efield = iEx, iEz
                          if (isWithinBounds(Efield, iii, jjj, kkk)) then
                            if (isThinWire(Efield, iii, jjj, kkk)) then
                              conta = conta + 1
                              call writeSerialized(ii, i, conta, iii, jjj, kkk, &
                                                   currentType(Efield), &
                                                   iabs(tag_numbers%getEdgeTag(Efield, iii, jjj, kkk)))
                            end if

                            if (.not. isMediaVacuum(Efield, iii, jjj, kkk) .and. &
                                .not. isSplitOrAdvanced(Efield, iii, jjj, kkk)) then
                              conta = conta + 1
                              call writeSerialized(ii, i, conta, iii, jjj, kkk, &
                                                   currentType(Efield), &
                                                   iabs(tag_numbers%getEdgeTag(Efield, iii, jjj, kkk)))
                            end if

                          end if
                        end do

                      else !si es mapvtk
                        !si es mapvtk y si no es vacio
                        do Efield = iEx, iEz
                          call assignMedia(imed, imed1, imed2, imed3, imed4, Efield, iii, jjj, kkk)
                          call contabordes(sgg, imed, imed1, imed2, imed3, imed4, EsBorde, SINPML_fullsize, Efield, iii, jjj, kkk)
                          if (EsBorde) then
                            conta = conta + 1
                            call writeSerialized(ii, i, conta, iii, jjj, kkk, &
                                                 currentType(Efield), &
                                                 tag_numbers%getEdgeTag(Efield, iii, jjj, kkk))
                          end if
                        end do

                      end if
                      !
                    end do
                  end do
                end do

                     !!!
                if (field == mapvtk) then
                  INIT = .false.; geom = .true.; asigna = .false.; magnetic = .false.; electric = .true.
  call nodalvtk(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,media%sggMtag,tag_numbers, &
                                init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)
        call wirebundlesvtk(sgg, init, geom, asigna, conta, i, ii, output, Ntimeforvolumic, wiresflavor, media%sggMtag, tag_numbers)
                end if
                     !!!
                do kkk = sgg%Observation(ii)%P(i)%ZI, sgg%Observation(ii)%P(i)%ZE
                  do jjj = sgg%observation(ii)%P(i)%YI, sgg%observation(ii)%P(i)%YE
                    do iii = sgg%observation(ii)%P(i)%XI, sgg%observation(ii)%P(i)%XE
                      if (field /= mapvtk) then
                        do Hfield = iHx, iHz
                          if ((isPECorSurface(Hfield, iii, jjj, kkk) .or. &
                               field == blockCurrent(Hfield)) &
                              .and. isWithinBounds(Hfield, iii, jjj, kkk)) then
                            conta = conta + 1
                            call writeSerialized(ii, i, conta, iii, jjj, kkk, &
                                                 currentType(Hfield), &
                                                 iabs(tag_numbers%getFaceTag(Hfield, iii, jjj, kkk)))
                          end if
                        end do
                      else !mapvtk y si no es vacio, asimilo la salida a corrientes iBloqueJ? para que vtk.f90 los escriba en quads
                        do Hfield = iHx, iHz
                          if (.not. isMediaVacuum(Hfield, iii, jjj, kkk) .and. &
                              .not. isPML(Hfield, iii, jjj, kkk) .and. &
                              isWithinBounds(Hfield, iii, jjj, kkk)) then
                            conta = conta + 1
                            call writeSerialized(ii, i, conta, iii, jjj, kkk, &
                                                 currentType(Hfield), &
                                                 tag_numbers%getFaceTag(Hfield, iii, jjj, kkk))
                          end if

                          if (tag_numbers%getFaceTag(Hfield, iii, jjj, kkk) < 0 .and. &
                              (btest(iabs(tag_numbers%getFaceTag(Hfield, iii, jjj, kkk)), Hfield - 1)) .and. &
                              .not. isPML(Hfield, iii, jjj, kkk) .and. &
                              isWithinBounds(Hfield, iii, jjj, kkk)) then
                            conta = conta + 1
                            call writeSerialized(ii, i, conta, iii, jjj, kkk, &
                                                 currentType(Hfield), &
                                                 tag_numbers%getFaceTag(Hfield, iii, jjj, kkk))
                          end if
                        end do
                        !   endif
                      end if
                      !
                    end do
                  end do
                end do

                     !!!
                if (field == mapvtk) then
                  INIT = .false.; geom = .true.; asigna = .false.; magnetic = .true.; electric = .false.
  call nodalvtk(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,media%sggMtag,tag_numbers, &
                                init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)

                end if
                     !!!
                my_iostat = 0
9137            if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9137,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                write (19, '(a)', err=9137, iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))

                !erase pre-existing data unless this is a resuming simulation

                if (.not. resume) then
                  if (SGG%Observation(ii)%TimeDomain) then
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted')
                    write (output(ii)%item(i)%unit) '!END'
                    close (output(ii)%item(i)%unit, status='DELETE')
                    my_iostat = 0
9240                if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a,2x,i5)',9240,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path)),output(ii)%item(i)%unit
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted',err=9240,iostat=my_iostat,status='new',action='write') !

                    write (output(ii)%item(i)%unit) output(ii)%item(i)%columnas
                    do conta = 1, output(ii)%item(i)%columnas
                      write (output(ii)%item(i)%unit) &
                        output(ii)%item(i)%Serialized%eI(conta), &
                        output(ii)%item(i)%Serialized%eJ(conta), &
                        output(ii)%item(i)%Serialized%eK(conta), &
                        output(ii)%item(i)%Serialized%currentType(conta), &
                        output(ii)%item(i)%Serialized%sggMtag(conta) !added to resuming file 121020
                    end do
                  elseif (SGG%Observation(ii)%FreqDomain) then
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted')
                    write (output(ii)%item(i)%unit) '!END'
                    close (output(ii)%item(i)%unit, status='DELETE')
                  end if !no need to keep it open
                  !wipe out duplicate data after non synchronous data and field resuming !later
                else !SE RESUMEA
                  inquire (file=trim(adjustl(output(ii)%item(i)%path)), exist=existe)
                  if (.not. existe) then
      call stoponerror(layoutnumber, size, 'Data files for resuming non existent (Volume) '//trim(adjustl(output(ii)%item(i)%path)))
                  end if
                  open (output(ii)%item(i)%unit, access='sequential', file=trim(adjustl(output(ii)%item(i)%path)), &
                        form='unformatted')

                  if ((SGG%Observation(ii)%TimeDomain) .and. (sgg%observation(ii)%P(1)%what /= mapvtk)) then
                    !
                    read (output(ii)%item(i)%unit) ndum
              if (output(ii)%item(i)%columnas /= ndum) call stoponerror(layoutnumber, size, 'BUGGYError reading resuming files () ')
                    do conta = 1, output(ii)%item(i)%columnas
                      read (output(ii)%item(i)%unit) ndum, ndum, ndum, ndum, ndum
                    end do
                    cutting3: do
                      read (output(ii)%item(i)%unit, end=699) at
          if (output(ii)%item(i)%columnas /= 0) read (output(ii)%item(i)%unit, end=699) (rdum, conta=1, output(ii)%item(i)%columnas)
                      output(ii)%TimesWritten = output(ii)%TimesWritten + 1
                      if (at > lastexecutedtime) then
                     print '(i4,a,a,2e19.9e3)', quienmpi, 'Cutting 3 ', trim(adjustl(output(ii)%item(i)%path)), at, lastexecutedtime
                        backspace (output(ii)%item(i)%unit)
                        endfile (output(ii)%item(i)%unit)
                        exit cutting3
                      end if
                    end do cutting3
699                 continue
                           !!            backspace(output(ii)%item(i)%unit) !machaco el timeswritten proque solo puede haber uno al final
                    close (output(ii)%item(i)%unit)
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), position='append', &
                          form='unformatted')
                  elseif ((SGG%Observation(ii)%FreqDomain) .and. (sgg%observation(ii)%P(1)%what /= mapvtk)) then
                    !
                    read (output(ii)%item(i)%unit) ndum
                    do conta = 1, output(ii)%item(i)%columnas
                      read (output(ii)%item(i)%unit) ndum, ndum, ndum, ndum, ndum
                    end do
                    read (output(ii)%item(i)%unit) at
                    if (int(at/sgg%dt) /= (initialtimestep - 1)) then
 write (buff,*) nint(at/sgg%dt), initialtimestep-1,' Data files for resuming 3D freq domain probes might be corrupt. Continuing....'
                      call print11(layoutnumber, buff)
                    end if
                    do N = 1, output(ii)%NumFreqs
                      read (output(ii)%item(i)%unit, end=6919) rdum
                      if (output(ii)%item(i)%columnas /= 0) then
                        do conta = 1, output(ii)%item(i)%columnas
                          read (output(ii)%item(i)%unit, end=6917) &
                            output(ii)%item(i)%Serialized%ValorComplex_x(N, conta), &
                            output(ii)%item(i)%Serialized%ValorComplex_y(N, conta), &
                            output(ii)%item(i)%Serialized%ValorComplex_z(N, conta)
6917                      read (output(ii)%item(i)%unit, end=6918) &
                            output(ii)%item(i)%Serialized%ValorComplex_Ex(N, conta), &
                            output(ii)%item(i)%Serialized%ValorComplex_Ey(N, conta), &
                            output(ii)%item(i)%Serialized%ValorComplex_Ez(N, conta)
6918                      read (output(ii)%item(i)%unit, end=6919) &
                            output(ii)%item(i)%Serialized%ValorComplex_Hx(N, conta), &
                            output(ii)%item(i)%Serialized%ValorComplex_Hy(N, conta), &
                            output(ii)%item(i)%Serialized%ValorComplex_Hz(N, conta)
                        end do
                      end if
                      if (SGG%Observation(ii)%transfer) then
  output(ii)%item(i)%Serialized%ValorComplex_x = output(ii)%item(i)%Serialized%ValorComplex_x*output(ii)%dftEntrada(N) !desnormaliza
  output(ii)%item(i)%Serialized%ValorComplex_y = output(ii)%item(i)%Serialized%ValorComplex_y*output(ii)%dftEntrada(N) !desnormaliza
  output(ii)%item(i)%Serialized%ValorComplex_z = output(ii)%item(i)%Serialized%ValorComplex_z*output(ii)%dftEntrada(N) !desnormaliza

output(ii)%item(i)%Serialized%ValorComplex_Ex = output(ii)%item(i)%Serialized%ValorComplex_Ex*output(ii)%dftEntrada(N) !desnormaliza
output(ii)%item(i)%Serialized%ValorComplex_Ey = output(ii)%item(i)%Serialized%ValorComplex_Ey*output(ii)%dftEntrada(N) !desnormaliza
output(ii)%item(i)%Serialized%ValorComplex_Ez = output(ii)%item(i)%Serialized%ValorComplex_Ez*output(ii)%dftEntrada(N) !desnormaliza

output(ii)%item(i)%Serialized%ValorComplex_Hx = output(ii)%item(i)%Serialized%ValorComplex_Hx*output(ii)%dftEntrada(N) !desnormaliza
output(ii)%item(i)%Serialized%ValorComplex_Hy = output(ii)%item(i)%Serialized%ValorComplex_Hy*output(ii)%dftEntrada(N) !desnormaliza
output(ii)%item(i)%Serialized%ValorComplex_Hz = output(ii)%item(i)%Serialized%ValorComplex_Hz*output(ii)%dftEntrada(N) !desnormaliza
                      end if
                    end do
6919                continue
                    close (output(ii)%item(i)%unit, status='delete')
                  end if
                end if
              end if
                  !!!!!!!!!!!!!!!!!!fin vtk
              !Volumic probes
            case (iMEC, iMHC, iExC, iEyC, iEzC, iHxC, iHyC, iHzC)
              if (sgg%Observation(ii)%Volumic) then !they are necssaryly
                if (sgg%Observation(ii)%nP /= 1) then
                  call stoponerror(layoutnumber, size, 'ERROR! More than a volumic probe per group')
                end if
                !readjust correctly the calculation region
                select case (field)
                case (iExC, iEyC, iHzC, iMhC)
                  sgg%observation(ii)%P(i)%ZI = max(sgg%Sweep(fieldo(field, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
                  sgg%observation(ii)%P(i)%ZE = min(sgg%Sweep(fieldo(field, 'Z'))%ZE - 1, sgg%observation(ii)%P(i)%ZE)
                case (iEzC, iHxC, iHyC, iMeC)
                  sgg%observation(ii)%P(i)%ZI = max(sgg%Sweep(fieldo(field, 'Z'))%ZI, sgg%observation(ii)%P(i)%ZI)
                  sgg%observation(ii)%P(i)%ZE = min(sgg%Sweep(fieldo(field, 'Z'))%ZE, sgg%observation(ii)%P(i)%ZE)
                end select
                !solo acepto que P(1:1) !!!
                write (chari, '(i7)') sgg%observation(ii)%P(i)%XI
                write (charj, '(i7)') sgg%observation(ii)%P(i)%YI
                write (chark, '(i7)') sgg%observation(ii)%P(i)%ZI
                write (chari2, '(i7)') sgg%observation(ii)%P(i)%XE
                write (charj2, '(i7)') sgg%observation(ii)%P(i)%YE
                write (chark2, '(i7)') sgg%observation(ii)%P(i)%ZE

                !mpidir 190319      !desrotacion para que los nombres sean correctos
                if (mpidir == 3) then
                  extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
                             trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                  select case (field)
                  case (iExC)
                    prefix_field = prefix(iExC)
                  case (iEyC)
                    prefix_field = prefix(iEyC)
                  case (iEzC)
                    prefix_field = prefix(iEzC)
                  case (iHxC)
                    prefix_field = prefix(iHxC)
                  case (iHyC)
                    prefix_field = prefix(iHyC)
                  case (iHzC)
                    prefix_field = prefix(iHzC)
                  case default
                    prefix_field = prefix(field)
                  end select
                elseif (mpidir == 2) then
                  extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
                             trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                  select case (field)
                  case (iExC)
                    prefix_field = prefix(iEzC)
                  case (iEyC)
                    prefix_field = prefix(iExC)
                  case (iEzC)
                    prefix_field = prefix(iEyC)
                  case (iHxC)
                    prefix_field = prefix(iHzC)
                  case (iHyC)
                    prefix_field = prefix(iHxC)
                  case (iHzC)
                    prefix_field = prefix(iHyC)
                  case default
                    prefix_field = prefix(field)
                  end select
                elseif (mpidir == 1) then
                  extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
                             trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                  select case (field)
                  case (iExC)
                    prefix_field = prefix(iEyC)
                  case (iEyC)
                    prefix_field = prefix(iEzC)
                  case (iEzC)
                    prefix_field = prefix(iExC)
                  case (iHxC)
                    prefix_field = prefix(iHyC)
                  case (iHyC)
                    prefix_field = prefix(iHzC)
                  case (iHzC)
                    prefix_field = prefix(iHxC)
                  case default
                    prefix_field = prefix(field)
                  end select
                else
                  call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
                end if
                !
                !
                ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
                !cada mpi su nombre
                output(ii)%item(i)%path = trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.bin'

                !
                unit = unit + 1
                if (unit >= 2.0_RKIND**31.0_RKIND - 1.0_RKIND) then
                  call stoponerror(layoutnumber, size, 'Excesive number of probes')
                end if
                output(ii)%item(i)%unit = unit

                     !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
                call checkduplicatenames
                     !!!!!!
                !
                output(ii)%item(i)%columnas = (sgg%Observation(ii)%P(i)%XE - sgg%Observation(ii)%P(i)%XI + 1)* &
                                              (sgg%observation(ii)%P(i)%YE - sgg%observation(ii)%P(i)%YI + 1)* &
                                              (sgg%observation(ii)%P(i)%ZE - sgg%observation(ii)%P(i)%ZI + 1)
                !
                !ojo por si algun dia esto molestara a Cray
                IF (SGG%Observation(ii)%TimeDomain) THEN

                  !replico los ifs de transferencia y escritura
                  ntini = 0
                  ntfin = 0
                  first = .true.
                  do Ntime = initialtimestep, finaltimestep
                    at = sgg%tiempo(Ntime)
                    Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                    if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                      Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                      if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2))) then
                        if (first) then
                          ntini = ntimeforvolumic
                          first = .false.
                        end if
                        ntfin = ntimeforvolumic
                      end if
                    end if
                  end do

                        !!!                             ntinI=0
                        !!!                             ntfin=min(int((finaltimestep*sgg%dt-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos), &
                        !!!                                       int((sgg%Observation(ii)%FinalTIME-sgg%OBSERVATION(ii)%InitialTime)/sgg%dt/output(ii)%Trancos))+1

                  memo = memo + RKIND* &
                         (ntfin - ntini + 1)* &
                         (output(ii)%item(i)%XEtrancos - output(ii)%item(i)%XItrancos + 1)* &
                         (output(ii)%item(i)%YEtrancos - output(ii)%item(i)%YItrancos + 1)* &
                         (output(ii)%item(i)%ZEtrancos - output(ii)%item(i)%ZItrancos + 1)

                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'ERROR: Recompile: excesive memory for the 3D probes.'// &
                    &                                   'Recompile increasing MaxMemoryProbes')
                  end if

                  allocate (output(ii)%item(i)%valor3D(ntinI:ntfin, &
                                                       output(ii)%item(i)%XItrancos:output(ii)%item(i)%XEtrancos, &
                                                       output(ii)%item(i)%YItrancos:output(ii)%item(i)%YEtrancos, &
                                                       output(ii)%item(i)%ZItrancos:output(ii)%item(i)%ZEtrancos))
                  output(ii)%item(i)%valor3D = 0.0_RKIND

                ELSEIF (SGG%Observation(ii)%FreqDomain) THEN
                  memo = memo + RKIND*output(ii)%NumFreqs*output(ii)%item(i)%columnas + 16*output(ii)%item(i)%columnas ! 4 integers de 4 bytes
                  if (memo > MaxMemoryProbes) then
                    call stoponerror(layoutnumber, size, 'ERROR: Recompile: excesive memory for the probes.'// &
                    &                                   'Recompile increasing MaxMemoryProbes')
                  end if

                  allocate (output(ii)%item(i)%valor3DComplex(1:output(ii)%NumFreqs, 1:3, & !tres posibles componentes
                                                              output(ii)%item(i)%XItrancos:output(ii)%item(i)%XEtrancos, &
                                                              output(ii)%item(i)%YItrancos:output(ii)%item(i)%YEtrancos, &
                                                              output(ii)%item(i)%ZItrancos:output(ii)%item(i)%ZEtrancos))
                  output(ii)%item(i)%valor3DComplex = (0.0_RKIND, 0.0_RKIND)
                end if
                !
                my_iostat = 0
9234            if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9234,layoutnumber,trim(adjustl(nEntradaRoot))//'_Outputrequests_'//trim(adjustl(whoamishort))//'.txt'
                write (19, '(a)', err=9234, iostat=my_iostat) trim(adjustl(output(ii)%item(i)%path))
                !erase pre-existing data unless this is a resuming simulation

                if (.not. resume) then
                  if (SGG%Observation(ii)%TimeDomain) then
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted')
                    write (output(ii)%item(i)%unit) '!END'
                    close (output(ii)%item(i)%unit, status='DELETE')
                    my_iostat = 0
9271                if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9271,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
                           open ( output(ii)%item(i)%unit,file= trim(adjustl(output(ii)%item(i)%path)),form='unformatted',err=9271,iostat=my_iostat,status='new',action='write')
                    write (output(ii)%item(i)%unit) &
                      output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos, &
                      output(ii)%item(i)%YItrancos, output(ii)%item(i)%YEtrancos, &
                      output(ii)%item(i)%ZItrancos, output(ii)%item(i)%ZEtrancos
                          !!!&      sgg%observation(ii)%P(i)%xI,sgg%observation(ii)%P(i)%xE, &
                          !!!&      sgg%observation(ii)%P(i)%YI,sgg%observation(ii)%P(i)%YE, &
                          !!!&      sgg%observation(ii)%P(i)%zI,sgg%observation(ii)%P(i)%ZE
                    !wipe out duplicate data after non synchronous data and field resuming !later

                  elseif (SGG%Observation(ii)%FreqDomain) then
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted')
                    write (output(ii)%item(i)%unit) '!END'
                    close (output(ii)%item(i)%unit, status='DELETE')
                  end if !no need to keep it open
                else
                  if (SGG%Observation(ii)%TimeDomain) then
                    inquire (file=trim(adjustl(output(ii)%item(i)%path)), exist=existe)
                    if (.not. existe) then
call stoponerror(layoutnumber,size,'Data files for resuming non existent (volume xdmf...) '//trim(adjustl(output(ii)%item(i)%path)))
                    end if
                    open (output(ii)%item(i)%unit, access='sequential', file=trim(adjustl(output(ii)%item(i)%path)), &
                          form='unformatted')
                    read (output(ii)%item(i)%unit) ndum, ndum, ndum, ndum, ndum, ndum
                    cutting4: do
                      read (output(ii)%item(i)%unit, end=6999) at
                      do k1 = output(ii)%item(i)%ZItrancos, output(ii)%item(i)%ZEtrancos !sgg%Observation(ii)%P(i)%ZI,sgg%Observation(ii)%P(i)%ZE
                        do j1 = output(ii)%item(i)%YItrancos, output(ii)%item(i)%YEtrancos !sgg%Observation(ii)%P(i)%YI,sgg%Observation(ii)%P(i)%YE
                          read (output(ii)%item(i)%unit, end=6999) (rdum, &
                          &            i1=output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos) ! sgg%Observation(ii)%P(i)%XI,sgg%Observation(ii)%P(i)%XE)
                        end do
                      end do
                      output(ii)%TimesWritten = output(ii)%TimesWritten + 1
                      if (at > lastexecutedtime) then
                     print '(i4,a,a,2e19.9e3)', quienmpi, 'Cutting 4 ', trim(adjustl(output(ii)%item(i)%path)), at, lastexecutedtime
                        backspace (output(ii)%item(i)%unit)
                        endfile (output(ii)%item(i)%unit)
                        exit cutting4
                      end if
                    end do cutting4
6999                continue
                           !!!    backspace(output(ii)%item(i)%unit) !machaco el timeswritten proque solo puede haber uno al final
                    close (output(ii)%item(i)%unit)
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), position='append', &
                          form='unformatted')
                  elseif (SGG%Observation(ii)%FreqDomain) then
                    open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted')
                    read (output(ii)%item(i)%unit) ndum, ndum, ndum, ndum, ndum, ndum
                    read (output(ii)%item(i)%unit) at
                    if (nint(at/sgg%dt) /= (initialtimestep - 1)) then  !estaba mal? ponia initialstep) sin el -1 !261119. lo he cambiado
 write (buff,*) nint(at/sgg%dt), initialtimestep-1,' Data files for resuming 3D freq domain probes might be corrupt. Continuing....'
                      call print11(layoutnumber, buff)
                    end if

                    DO N = 1, output(ii)%NumFreqs
                      read (output(ii)%item(i)%unit) rdum
                      do compo = 1, 3
                      do k1t = output(ii)%item(i)%ZItrancos, output(ii)%item(i)%ZEtrancos
                        do j1t = output(ii)%item(i)%YItrancos, output(ii)%item(i)%YEtrancos
                          read (output(ii)%item(i)%unit) (output(ii)%item(i)%valor3DComplex(N, compo, i1t, j1t, k1t), &
                       &                                i1t=output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos)
                        end do
                      end do
                      end do
if (sgg%Observation(ii)%Transfer) output(ii)%item(i)%valor3DComplex = output(ii)%item(i)%valor3DComplex*output(ii)%dftEntrada(n) !desnormaliza
                    end do
                    close (output(ii)%item(i)%unit, status='delete')
                  end if
                end if
              end if

            case (farfield)
              ThereAreFarFields = .true.
              !
              write (chari, '(i7)') sgg%observation(ii)%P(1)%XI
              write (charj, '(i7)') sgg%observation(ii)%P(1)%YI
              write (chark, '(i7)') sgg%observation(ii)%P(1)%ZI
              write (chari2, '(i7)') sgg%observation(ii)%P(1)%XE
              write (charj2, '(i7)') sgg%observation(ii)%P(1)%YE
              write (chark2, '(i7)') sgg%observation(ii)%P(1)%ZE
              !mpidir 190319      !desrotacion para que los nombres sean correctos
              if (mpidir == 3) then
                extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
                           trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                prefix_field = prefix(field)
              elseif (mpidir == 2) then
                extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
                           trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                prefix_field = prefix(field)
              elseif (mpidir == 1) then
                extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
                           trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                prefix_field = prefix(field)
              else
                call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
              end if
              !
              !
              ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
              output(ii)%item(i)%path = trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))// &
                                        trim(adjustl(extpoint))//'.dat'

              output(ii)%item(i)%columnas = 1
              !
              unit = unit + 1
              if (unit >= 2.0_RKIND**31.0_RKIND - 1.0_RKIND) then
                call stoponerror(layoutnumber, size, 'Excesive number of probes')
              end if
              output(ii)%item(i)%unit = unit
                  !!!busca nombres de ficheros por duplicado y resuelve la duplicidad
              call checkduplicatenames
                  !!!!!!
              !inicializacion especifica del farfield
      call InitFarField(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,layoutnumber,size, &
                                b, resume, &
                                output(ii)%item(i)%unit, &
                                output(ii)%item(i)%path, &
                                sgg%observation(ii)%P(1)%XI, &
                                sgg%observation(ii)%P(1)%XE, &
                                sgg%observation(ii)%P(1)%YI, &
                                sgg%observation(ii)%P(1)%YE, &
                                sgg%observation(ii)%P(1)%ZI, &
                                sgg%observation(ii)%P(1)%ZE, &
                                sgg%observation(ii)%InitialFreq, &
                                sgg%observation(ii)%FinalFreq, &
                                sgg%observation(ii)%FreqStep, &
                                sgg%observation(ii)%phiStart, &
                                sgg%observation(ii)%phiStop, &
                                sgg%observation(ii)%phiStep, &
                                sgg%observation(ii)%thetaStart, &
                                sgg%observation(ii)%thetaStop, &
                                sgg%observation(ii)%thetaStep, &
                                sgg%observation(ii)%FileNormalize, SINPML_fullsize, facesNF2FF, NF2FFDecim &
#ifdef CompileWithMPI
                                , output(ii)%item(i)%MPISubComm, output(ii)%item(i)%MPIRoot &
#endif
                                , eps0, mu0)
              !no es necesario hacer wipe out pq en DF se van machacando
            end select
          end do loop_ob
            !!!!        endif !del time domain !NO ES PRECISO 25/02/14
        end do !del ii=1,numberrequest

        write (19, '(a)') '!END '
        close (19)

      end if

      return

      contains

      subroutine writeSerialized(i_out, i_item, conta, i, j, k, current, tag)
        integer(kind=4) :: i_out, i_item, conta, i, j, k, current
        integer(kind=IKINDMTAG) :: tag
        output(i_out)%item(i_item)%Serialized%eI(conta) = i
        output(i_out)%item(i_item)%Serialized%eJ(conta) = j
        output(i_out)%item(i_item)%Serialized%eK(conta) = k
        output(i_out)%item(i_item)%Serialized%currentType(conta) = current
        output(i_out)%item(i_item)%Serialized%sggMtag(conta) = tag
      end subroutine

      integer function blockCurrent(field)
        integer(kind=4) :: field
        select case (field)
        case (iHx)
          blockCurrent = iCurX
        case (iHy)
          blockCurrent = iCurY
        case (iHz)
          blockCurrent = iCurZ
        case default
          call StopOnError(layoutnumber, size, 'field is not H field')
        end select
      end function

      integer function currentType(field)
        integer(kind=4) :: field
        select case (field)
        case (iEx)
          currentType = iJx
        case (iEy)
          currentType = iJy
        case (iEz)
          currentType = iJz
        case (iHx)
          currentType = iBloqueJx
        case (iHy)
          currentType = iBloqueJy
        case (iHz)
          currentType = iBloqueJz
        case default
          call StopOnError(layoutnumber, size, 'field is not a E or H field')
        end select
      end function

      function getMedia(field, i, j, k) result(res)
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: res
        integer(kind=4) :: field, i, j, k
        select case (field)
        case (iEx)
          res = media%sggMiEx(i, j, k)
        case (iEy)
          res = media%sggMiEy(i, j, k)
        case (iEz)
          res = media%sggMiEz(i, j, k)
        case (iHx)
          res = media%sggMiHx(i, j, k)
        case (iHy)
          res = media%sggMiHy(i, j, k)
        case (iHz)
          res = media%sggMiHz(i, j, k)
        case default
          call StopOnError(layoutnumber, size, 'Unrecognized field')
        end select
      end function

      logical function isMediaVacuum(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media, vacuum = 1
        media = getMedia(field, i, j, k)
        isMediaVacuum = (media == vacuum)
      end function

      logical function isSplitOrAdvanced(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isSplitOrAdvanced = sgg%med(media)%is%split_and_useless .or. &
                            sgg%med(media)%is%already_YEEadvanced_byconformal

      end function

      logical function isThinWireWithinBounds(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        isThinWireWithinBounds = isThinWire(field, i, j, k) .and. &
                                 isWithinBounds(field, i, j, k)
      end function

      logical function isPECorSurface(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isPECorSurface = sgg%med(media)%is%PEC .or. &
                         sgg%med(media)%is%Surface
      end function

      logical function isThinWire(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isThinWire = sgg%Med(media)%is%ThinWire
      end function

      logical function isPML(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isPML = sgg%med(media)%is%PML
      end function

      logical function isPEC(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isPEC = sgg%med(media)%is%PEC
      end function

      logical function isWithinBounds(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        isWithinBounds = (i <= SINPML_fullsize(field)%XE) .and. &
                         (j <= SINPML_fullsize(field)%YE) .and. &
                         (k <= SINPML_fullsize(field)%ZE)
      end function

      subroutine assignMedia(m, m1, m2, m3, m4, dir, i, j, k)
        integer(kind=4), intent(inout) :: m, m1, m2, m3, m4
        integer(kind=4) :: dir, i, j, k
        m = getMedia(dir, i, j, k)
        m1 = getMedia(4 + modulo(dir, 3), i, j, k)
        m2 = getMedia(4 + modulo(dir, 3), i - merge(1, 0, dir == iEy), j - merge(1, 0, dir == iEz), k - merge(1, 0, dir == iEx))
        m3 = getMedia(4 + modulo(dir + 1, 3), i, j, k)
        m4 = getMedia(4 + modulo(dir + 1, 3), i - merge(1, 0, dir == iEz), j - merge(1, 0, dir == iEx), k - merge(1, 0, dir == iEy))
      end subroutine

      subroutine checkduplicatenames
        integer(Kind=4) :: n_ii, n_i, off
        p1 = output(ii)%item(i)%path
        do n_ii = 1, ii
          off = sgg%Observation(n_ii)%nP
          if (n_ii == ii) off = i - 1
          do n_i = 1, off
            if (sgg%Observation(n_ii)%P(n_i)%What /= nothing) then
              p2 = output(n_ii)%item(n_i)%path
              if (trim(adjustl(p1)) == trim(adjustl(p2))) then
                write (charNO, '(i7)') output(ii)%item(i)%unit
                output(ii)%item(i)%path = trim(adjustl(output(ii)%item(i)%path))//'_duplicate_'//trim(adjustl(charNO))//'.dat'
              end if
            end if
          end do
        end do
        return
      end subroutine

      !*************************

      subroutine crea_gnuplot

        buff2 = trim(adjustl(nEntradaRoot))//'_gnuplot.pl'
        thefile = openfile_mpi(layoutnumber, buff2)

         !!!!!
        conta = 0
        do ii = 1, sgg%NumberRequest
          do i = 1, sgg%Observation(ii)%nP
            I1 = sgg%observation(ii)%P(i)%XI
            J1 = sgg%observation(ii)%P(i)%YI
            K1 = sgg%observation(ii)%P(i)%ZI
            I2 = sgg%observation(ii)%P(i)%XE
            J2 = sgg%observation(ii)%P(i)%YE
            K2 = sgg%observation(ii)%P(i)%ZE
            NO = sgg%observation(ii)%P(i)%NODE
            write (chari, '(i7)') i1
            write (charj, '(i7)') j1
            write (chark, '(i7)') k1
            field = sgg%observation(ii)%P(i)%What
            select case (field)
            case (iEx, iEy, iEz, iVx, iVy, iVz, iJx, iJy, iJz, iHx, iHy, iHz)
              conta = conta + 1
              !mpidir 190319 !desrotacion para que los nombres sean correctos
              if (mpidir == 3) then
                extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
                select case (field)
                case (iEx)
                  prefix_field = prefix(iEx)
                case (iEy)
                  prefix_field = prefix(iEy)
                case (iEz)
                  prefix_field = prefix(iEz)
                case (iJx)
                  prefix_field = prefix(iJx)
                case (iJy)
                  prefix_field = prefix(iJy)
                case (iJz)
                  prefix_field = prefix(iJz)
                case (iVx)
                  prefix_field = prefix(iVx)
                case (iVy)
                  prefix_field = prefix(iVy)
                case (iVz)
                  prefix_field = prefix(iVz)
                case (iHx)
                  prefix_field = prefix(iHx)
                case (iHy)
                  prefix_field = prefix(iHy)
                case (iHz)
                  prefix_field = prefix(iHz)
                case default
                  prefix_field = prefix(field)
                end select
              elseif (mpidir == 2) then
                extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
                select case (field)
                case (iEx)
                  prefix_field = prefix(iEz)
                case (iEy)
                  prefix_field = prefix(iEx)
                case (iEz)
                  prefix_field = prefix(iEy)
                case (iJx)
                  prefix_field = prefix(iJz)
                case (iJy)
                  prefix_field = prefix(iJx)
                case (iJz)
                  prefix_field = prefix(iJy)
                case (iVx)
                  prefix_field = prefix(iVz)
                case (iVy)
                  prefix_field = prefix(iVx)
                case (iVz)
                  prefix_field = prefix(iVy)
                case (iHx)
                  prefix_field = prefix(iHz)
                case (iHy)
                  prefix_field = prefix(iHx)
                case (iHz)
                  prefix_field = prefix(iHy)
                case default
                  prefix_field = prefix(field)
                end select
              elseif (mpidir == 1) then
                extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
                select case (field)
                case (iEx)
                  prefix_field = prefix(iEy)
                case (iEy)
                  prefix_field = prefix(iEz)
                case (iEz)
                  prefix_field = prefix(iEx)
                case (iJx)
                  prefix_field = prefix(iJy)
                case (iJy)
                  prefix_field = prefix(iJz)
                case (iJz)
                  prefix_field = prefix(iJx)
                case (iVx)
                  prefix_field = prefix(iVy)
                case (iVy)
                  prefix_field = prefix(iVz)
                case (iVz)
                  prefix_field = prefix(iVx)
                case (iHx)
                  prefix_field = prefix(iHy)
                case (iHy)
                  prefix_field = prefix(iHz)
                case (iHz)
                  prefix_field = prefix(iHx)
                case default
                  prefix_field = prefix(field)
                end select
              else
                call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
              end if
              !
              if ((field == iJx) .or. (field == iJy) .or. (field == iJz)) then
                write (charNO, '(i7)') NO
                !append the number of the segment
                extpoint = trim(adjustl(extpoint))//'_s'//trim(adjustl(charNO))
              end if
              ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
              !do not use layername since no two observations from different layers will overlap
              path = CHAR(39)//trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'//CHAR(39)
              write (buff2, *) 'set term x11 persist', conta + 10000*layoutnumber
              call writefile_mpi(layoutnumber, thefile, buff2)
              write (buff2, *) 'plot ', trim(adjustl(path)), ' using 1:2 every 1::2 with lines '
              call writefile_mpi(layoutnumber, thefile, buff2)
            case (iBloqueJx, iBloqueJy, iBloqueJz, iBloqueMx, iBloqueMy, iBloqueMz)
              conta = conta + 1
              write (chari2, '(i7)') i2
              write (charj2, '(i7)') j2
              write (chark2, '(i7)') k2
              !mpidir 190319    !desrotacion para que los nombres sean correctos
              if (mpidir == 3) then
                extpoint = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))//'__'// &
                           trim(adjustl(chari2))//'_'//trim(adjustl(charj2))//'_'//trim(adjustl(chark2))
                select case (field)
                case (iBloqueJX)
                  prefix_field = prefix(iBloqueJX)
                case (iBloqueJY)
                  prefix_field = prefix(iBloqueJY)
                case (iBloqueJZ)
                  prefix_field = prefix(iBloqueJZ)
                case (iBloqueMX)
                  prefix_field = prefix(iBloqueMX)
                case (iBloqueMY)
                  prefix_field = prefix(iBloqueMY)
                case (iBloqueMZ)
                  prefix_field = prefix(iBloqueMZ)
                case default
                  prefix_field = prefix(field)
                end select
              elseif (mpidir == 2) then
                extpoint = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))//'__'// &
                           trim(adjustl(charj2))//'_'//trim(adjustl(chark2))//'_'//trim(adjustl(chari2))
                select case (field)
                case (iBloqueJX)
                  prefix_field = prefix(iBloqueJZ)
                case (iBloqueJY)
                  prefix_field = prefix(iBloqueJX)
                case (iBloqueJZ)
                  prefix_field = prefix(iBloqueJY)
                case (iBloqueMX)
                  prefix_field = prefix(iBloqueMZ)
                case (iBloqueMY)
                  prefix_field = prefix(iBloqueMX)
                case (iBloqueMZ)
                  prefix_field = prefix(iBloqueMY)
                case default
                  prefix_field = prefix(field)
                end select
              elseif (mpidir == 1) then
                extpoint = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))//'__'// &
                           trim(adjustl(chark2))//'_'//trim(adjustl(chari2))//'_'//trim(adjustl(charj2))
                select case (field)
                case (iBloqueJX)
                  prefix_field = prefix(iBloqueJY)
                case (iBloqueJY)
                  prefix_field = prefix(iBloqueJZ)
                case (iBloqueJZ)
                  prefix_field = prefix(iBloqueJX)
                case (iBloqueMX)
                  prefix_field = prefix(iBloqueMY)
                case (iBloqueMY)
                  prefix_field = prefix(iBloqueMZ)
                case (iBloqueMZ)
                  prefix_field = prefix(iBloqueMX)
                case default
                  prefix_field = prefix(field)
                end select
              else
                call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
              end if
              !
              !
              ext = trim(adjustl(nEntradaRoot))//'_'//trim(adjustl(sgg%observation(ii)%outputrequest))
              !
              path = CHAR(39)//trim(adjustl(ext))//'_'//trim(adjustl(prefix_field))//trim(adjustl(extpoint))//'.dat'//CHAR(39)
              write (buff2, *) 'set term x11 persist', conta + 10000*layoutnumber
              call writefile_mpi(layoutnumber, thefile, buff2)
              write (buff2, *) 'plot ', trim(adjustl(path)), ' using 1:2 every 1::2 with lines '
              call writefile_mpi(layoutnumber, thefile, buff2)
            end select
          end do
        end do !del ii=1,numberrequest

        buff2 = trim(adjustl(nEntradaRoot))//'_gnuplot.pl'
        call closefile_mpi(layoutnumber, size, buff2, thefile)

        return
      end subroutine crea_gnuplot

    end subroutine InitObservation

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Closes observation stuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine CloseObservationFiles(sgg, layoutnumber, size, singlefilewrite, initialtimestep, lastexecutedtime, resume)
      type(SGGFDTDINFO), intent(IN)         ::  sgg
      integer(kind=4)  ::  i, ii, layoutnumber, field, initialtimestep, unidad, size, idum
      logical :: singlefilewrite, resume, incident, existe, wrotemaster
      REAL(KIND=RKIND)    ::  rdum1, rdum2, rdum3, rdum4, rdum5, rdum6, rdum
      REAL(KIND=RKIND_tiempo)    ::  lastexecutedtime
      character(LEN=BUFSIZE) :: chdum
      character(LEN=BUFSIZE)  ::  whoamishort
      integer :: my_iostat
      real(kind=RKIND_tiempo) :: at
      !!!

      write (whoamishort, '(i5)') layoutnumber + 1
      !
      if (sgg%NumPlaneWaves >= 1) incident = .true.
      do ii = 1, sgg%NumberRequest
        if (SGG%Observation(ii)%TimeDomain) then
          if (.not. SGG%Observation(ii)%Volumic) then
            wrotemaster = .false.
            do i = 1, sgg%Observation(ii)%nP
              if ((sgg%observation(ii)%P(i)%What /= nothing) .and. &
                  (SGG%Observation(ii)%P(1)%what /= farfield)) then !el farfield se cierra y abre de su propio modo
                field = sgg%observation(ii)%P(i)%what
                if (singlefilewrite .and. ((field == iEx) .or. (field == iEy) .or. (field == iEz) .or. &
                                           (field == iVx) .or. (field == iVy) .or. (field == iVz) .or. &
                                           (field == iJx) .or. (field == iJy) .or. (field == iJz) .or. &
                                           (field == iHx) .or. (field == iHy) .or. (field == iHz))) then
                  if (.not. wrotemaster) then
                    wrotemaster = .true.
                    close (output(ii)%item(i)%unitmaster)
                           open  (output(ii)%item(i)%unitmaster-1,file= trim(adjustl(output(ii)%item(i)%path))//'_'//trim(adjustl(whoamishort))//'_master.bin',form='unformatted')
                  else
                    rewind (output(ii)%item(i)%unitmaster - 1)
                  end if
                  if (.not. resume) then
                    unidad = output(ii)%item(i)%unit
                    open (unidad, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)))
                    write (unidad, *) '!END'
                    close (unidad, status='delete')
                    my_iostat = 0
9242                if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9242,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
     open (unidad, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)), err=9242, iostat=my_iostat, status='new', action='write')
                    write (unidad, '(a)') trim(adjustl(' t'//'              '// &
                                         trim(adjustl(output(ii)%item(i)%path))//'       '//trim(adjustl(suffix(field, incident)))))
                  else
                    inquire (file=trim(adjustl(output(ii)%item(i)%path)), exist=existe)
                    if (.not. existe) then
                              call stoponerror(layoutnumber,size,'Data files for resuming non existent (generic closing) '//trim(adjustl(output(ii)%item(i)%path)))
                    end if
                    !
                    open (output(ii)%item(i)%unit, recl=1000, access='sequential', file=trim(adjustl(output(ii)%item(i)%path)))
                    read (output(ii)%item(i)%unit, '(a)') chdum !first line contains characters
                    cutting: do
                      read (output(ii)%item(i)%unit, *, end=678) at
                      if (at > lastexecutedtime) then
                     print '(i4,a,a,2e19.9e3)', quienmpi, 'Cutting 5 ', trim(adjustl(output(ii)%item(i)%path)), at, lastexecutedtime
                        backspace (output(ii)%item(i)%unit)
                        endfile (output(ii)%item(i)%unit)
                        exit cutting
                      end if
                    end do cutting
678                 continue
                    close (output(ii)%item(i)%unit)
                    !
                    open (output(ii)%item(i)%unit, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)), position='append')
                  end if
                  !
                  do
                    select case (field)
                    case (iHx, iHy, iHz, iEx, iEy, iEz)
                      if (incident) then
                        read (output(ii)%item(i)%unitmaster - 1, end=777) idum, rdum1, rdum2, rdum3
                        if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit, fmt) rdum1, rdum2, rdum3
                      else
                        read (output(ii)%item(i)%unitmaster - 1, end=777) idum, rdum1, rdum2
                        if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit, fmt) rdum1, rdum2
                      end if
                    case (iJx, iJy, iJz)
                      read (output(ii)%item(i)%unitmaster - 1, end=777) idum, rdum1, rdum2, rdum3, rdum4, rdum5, rdum6
                  if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit, fmt) rdum1, rdum2, rdum3, rdum4, rdum5, rdum6
                    case default
                      read (output(ii)%item(i)%unitmaster - 1, end=777) idum, rdum1, rdum2 !caso de los votalges ivx, etc
                      if (idum == output(ii)%item(i)%unit) write (output(ii)%item(i)%unit, fmt) rdum1, rdum2
                    end select
                  end do
777               close (output(ii)%item(i)%unit)
                else
                  close (output(ii)%item(i)%unit)
                end if
              end if
            end do
          else !sondas volumicas
            do i = 1, sgg%Observation(ii)%nP
              if ((sgg%observation(ii)%P(i)%What /= nothing) .and. &
                  (SGG%Observation(ii)%P(1)%what /= farfield)) then
                !write(output(ii)%item(i)%unit) output(ii)%TimesWritten !el farfield se cierra y abre de su propio modo
                endfile (output(ii)%item(i)%unit)
                close (output(ii)%item(i)%unit)
              end if
            end do
               !!!!
          end if
        elseif (SGG%Observation(ii)%FreqDomain) then
          continue !nothing to do since the freq domain is updated by opening, writing and closing each time !esot hace que no se pueda hacer cutting (si hay fallos y restarteos)
        end if
        wrotemaster = .false.
        do i = 1, sgg%Observation(ii)%nP
          if (SGG%Observation(ii)%TimeDomain) then
            if ((sgg%observation(ii)%P(i)%What /= nothing) .and. &
                (SGG%Observation(ii)%P(1)%what /= farfield)) then !el farfield se cierra y abre de su propio modo
              field = sgg%observation(ii)%P(i)%what
              if (singlefilewrite .and. ((field == iEx) .or. (field == iEy) .or. (field == iEz) .or. &
                                         (field == iVx) .or. (field == iVy) .or. (field == iVz) .or. &
                                         (field == iJx) .or. (field == iJy) .or. (field == iJz) .or. &
                                         (field == iHx) .or. (field == iHy) .or. (field == iHz))) then
                if (.not. wrotemaster) then
                  wrotemaster = .true.
                  close (output(ii)%item(i)%unitmaster - 1, status='delete') !el binario no se precisa para nada
                end if
              end if
            end if
          end if !no hay singlewrite para sondas freqdomain
        end do
      end do

      return
    end subroutine CloseObservationFiles

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Upacks .bin files observation stuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine UnpackSingleFiles(sgg, layoutnumber, size, singlefilewrite, initialtimestep, resume)
      type(SGGFDTDINFO), intent(IN)         ::  sgg
      integer(kind=4)  ::  i, ii, layoutnumber, field, initialtimestep, unidad, size, idum
      logical :: singlefilewrite, resume, incident, existe, wrotemaster
      REAL(KIND=RKIND)    ::  rdum1, rdum2, rdum3, rdum4, rdum5, rdum6, rdum
      character(LEN=BUFSIZE) :: chdum
      character(LEN=BUFSIZE)  ::  whoamishort
      integer :: my_iostat
      !!!
      write (whoamishort, '(i5)') layoutnumber + 1
      !
      if (sgg%NumPlaneWaves >= 1) incident = .true.
      do ii = 1, sgg%NumberRequest
        if (SGG%Observation(ii)%TimeDomain) then
          if (.not. SGG%Observation(ii)%Volumic) then
            wrotemaster = .false.
            do i = 1, sgg%Observation(ii)%nP
              if ((sgg%observation(ii)%P(i)%What /= nothing) .and. &
                  (SGG%Observation(ii)%P(1)%what /= farfield)) then !el farfield se cierra y abre de su propio modo
                field = sgg%observation(ii)%P(i)%what
                if (singlefilewrite .and. ((field == iEx) .or. (field == iEy) .or. (field == iEz) .or. &
                                           (field == iVx) .or. (field == iVy) .or. (field == iVz) .or. &
                                           (field == iJx) .or. (field == iJy) .or. (field == iJz) .or. &
                                           (field == iHx) .or. (field == iHy) .or. (field == iHz))) then
                  rewind (output(ii)%item(i)%unitmaster)
                  unidad = 35
                  open (unidad, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)))
                  write (unidad, *) '!END'
                  close (unidad, status='delete')
                  my_iostat = 0
9243              if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9243,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
     open (unidad, recl=1000, file=trim(adjustl(output(ii)%item(i)%path)), err=9243, iostat=my_iostat, status='new', action='write')
                  write (unidad, '(a)') trim(adjustl(' t'//'              '// &
                                         trim(adjustl(output(ii)%item(i)%path))//'       '//trim(adjustl(suffix(field, incident)))))
                  !
                  do
                    select case (field)
                    case (iHx, iHy, iHz, iEx, iEy, iEz)
                      if (incident) then
                        read (output(ii)%item(i)%unitmaster, end=7778) idum, rdum1, rdum2, rdum3
                        if (idum == output(ii)%item(i)%unit) write (unidad, fmt) rdum1, rdum2, rdum3
                      else
                        read (output(ii)%item(i)%unitmaster, end=7778) idum, rdum1, rdum2
                        if (idum == output(ii)%item(i)%unit) write (unidad, fmt) rdum1, rdum2
                      end if
                    case (iJx, iJy, iJz)
                      read (output(ii)%item(i)%unitmaster, end=7778) idum, rdum1, rdum2, rdum3, rdum4, rdum5, rdum6
                      if (idum == output(ii)%item(i)%unit) write (unidad, fmt) rdum1, rdum2, rdum3, rdum4, rdum5, rdum6
                    case default
                      read (output(ii)%item(i)%unitmaster, end=7778) idum, rdum1, rdum2 !caso de los votalges ivx, etc
                      if (idum == output(ii)%item(i)%unit) write (unidad, fmt) rdum1, rdum2
                    end select
                  end do
7778              close (unidad)
                end if
              end if
            end do
          end if
        end if
      end do

      return
    end subroutine UnpackSingleFiles

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Updates the observed values. A nodal average is used for each field
   !!! The Wire modules uses its own updating procedure
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine UpdateObservation(sgg, media, tag_numbers, &
                      nTime, nInit, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh, wiresflavor, SINPML_fullsize, wirecrank, &
                                 noconformalmapvtk, b)
      !solo lo precisa de entrada farfield
      type(bounds_t)  ::  b
      logical :: noconformalmapvtk
      type(SGGFDTDINFO), intent(IN)         ::  sgg
      type(media_matrices_t), intent(in) :: media
      type(taglist_t) :: tag_numbers
      !---------------------------> inputs <----------------------------------------------------------
      type(limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      integer, intent(IN)  ::  nTime, nInit
      REAL(KIND=RKIND), intent(in), target     :: &
        Ex(sgg%alloc(iEx)%XI:sgg%alloc(iEx)%XE, sgg%alloc(iEx)%YI:sgg%alloc(iEx)%YE, sgg%alloc(iEx)%ZI:sgg%alloc(iEx)%ZE), &
        Ey(sgg%alloc(iEy)%XI:sgg%alloc(iEy)%XE, sgg%alloc(iEy)%YI:sgg%alloc(iEy)%YE, sgg%alloc(iEy)%ZI:sgg%alloc(iEy)%ZE), &
        Ez(sgg%alloc(iEz)%XI:sgg%alloc(iEz)%XE, sgg%alloc(iEz)%YI:sgg%alloc(iEz)%YE, sgg%alloc(iEz)%ZI:sgg%alloc(iEz)%ZE), &
        Hx(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE, sgg%alloc(iHx)%YI:sgg%alloc(iHx)%YE, sgg%alloc(iHx)%ZI:sgg%alloc(iHx)%ZE), &
        Hy(sgg%alloc(iHy)%XI:sgg%alloc(iHy)%XE, sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE, sgg%alloc(iHy)%ZI:sgg%alloc(iHy)%ZE), &
        Hz(sgg%alloc(iHz)%XI:sgg%alloc(iHz)%XE, sgg%alloc(iHz)%YI:sgg%alloc(iHz)%YE, sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)
      !--->
      REAL(KIND=RKIND), dimension(:), intent(in)   :: dxh(sgg%ALLOC(iEx)%XI:sgg%ALLOC(iEx)%XE), &
                                                      dyh(sgg%ALLOC(iEy)%YI:sgg%ALLOC(iEy)%YE), &
                                                      dzh(sgg%ALLOC(iEz)%ZI:sgg%ALLOC(iEz)%ZE), &
                                                      dxe(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE), &
                                                      dye(sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE), &
                                                      dze(sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)

      !---------------------------> variables locales <-----------------------------------------------
   integer( kind = 4)  ::  i, ii, i1, i2, j1, j2, k1, k2, i1_m, i2_m, j1_m, j2_m, k1_m, k2_m, field,jjx,jjy,jjz,if1,i1t,j1t,k1t,iff1
      integer(kind=4)  :: Efield, HField
      integer(kind=4)  ::  iii, kkk, jjj, jjj_m, iii_m, kkk_m, NtimeforVolumic, imed, imed1, imed2, imed3, imed4, medium
      logical :: esborde, wirecrank
      REAL(KIND=RKIND_tiempo)    ::  at
      real(kind=RKIND) :: jx, jy, jz, jdir, jdir1, jdir2
      complex(kind=ckind) :: z_cplx
      integer(kind=4) :: conta !para realmente dar tangenciales de campos en los medios superficiales
      character(len=*), INTENT(in) :: wiresflavor

      type(CurrentSegments), pointer  ::  segmDumm !segmento de hilo que se observa si lo hubiere
      !
#ifdef CompileWithBerengerWires
      type(TSegment), pointer  ::  segmDumm_Berenger !segmento de hilo que se observa si lo hubiere
#endif
      !
#ifdef CompileWithSlantedWires
      class(Segment), pointer  ::  segmDumm_Slanted !segmento de hilo que se observa si lo hubiere
#endif

      logical ::  INIT, GEOM, ASIGNA, electric, magnetic

      at = -1; jx = -1; jy = -1; jz = -1; jdir = -1; jdir1 = -1; jdir2 = -1  !para que gfortran no me diga que no las inicializo

      !---------------------------> empieza UpdateObservation <---------------------------------------
      do ii = 1, sgg%NumberRequest
        loop_obser: do i = 1, sgg%Observation(ii)%nP
          field = SGG%Observation(ii)%P(i)%what
          if (field /= nothing) then
            I1 = SGG%Observation(ii)%P(i)%XI
            J1 = SGG%Observation(ii)%P(i)%YI
            K1 = SGG%Observation(ii)%P(i)%ZI
            I2 = SGG%Observation(ii)%P(i)%XE !ojo estos no se usan salvo en Bloques y Volumics
            J2 = SGG%Observation(ii)%P(i)%YE
            K2 = SGG%Observation(ii)%P(i)%ZE
            !--->
            if (SGG%Observation(ii)%TimeDomain) then
              selectcase (field)
              case (iEx)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                output(ii)%item(i)%valor(nTime - nInit) = Ex(i1_m, j1_m, k1_m)
              case (iEy)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                output(ii)%item(i)%valor(nTime - nInit) = Ey(i1_m, j1_m, k1_m)
              case (iEz)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                output(ii)%item(i)%valor(nTime - nInit) = Ez(i1_m, j1_m, k1_m)
              case (iHx)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                output(ii)%item(i)%valor(nTime - nInit) = Hx(i1_m, j1_m, k1_m)
              case (iHy)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                output(ii)%item(i)%valor(nTime - nInit) = Hy(i1_m, j1_m, k1_m)
              case (iHz)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                output(ii)%item(i)%valor(nTime - nInit) = Hz(i1_m, j1_m, k1_m)
              case (iBloqueJx)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                k1_m = K1
                k2_m = K2
                do JJJ = j1, j2
                  JJJ_m = JJJ
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (Hy(i1_m, JJJ_m, k1_m - 1) - Hy(i1_m, JJJ_m, k2_m))*dyh(JJJ_m)
                end do
                !--->
                i1_m = I1
                j1_m = J1
                j2_m = J2
                do KKK = k1, k2
                  KKK_m = KKK
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (-Hz(i1_m, j1_m - 1, KKK_m) + Hz(i1_m, j2_m, KKK_m))*dzh(KKK_m)
                end do
              case (iBloqueJy)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                i2_m = I2
                do KKK = k1, k2
                  KKK_m = KKK
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (-Hz(i2_m, j1_m, KKK_m) + Hz(i1_m - 1, j1_m, KKK_m))*dzh(KKK_m)
                end do
                !--->
                j1_m = J1
                k1_m = K1
                k2_m = K2
                do III = i1, i2
                  III_m = III
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (Hx(III_m, j1_m, k2_m) - Hx(III_m, j1_m, k1_m - 1))*dxh(III_m)
                end do
              case (iBloqueJz)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                j1_m = J1
                k1_m = K1
                j2_m = J2
                do III = i1, i2
                  III_m = III
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (Hx(III_m, j1_m - 1, k1_m) - Hx(III_m, j2_m, k1_m))*dxh(III_m)
                end do
                !--->
                i1_m = I1
                k1_m = K1
                i2_m = I2
                do JJJ = j1, j2
                  JJJ_m = JJJ
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (-Hy(i1_m - 1, JJJ_m, k1_m) + Hy(i2_m, JJJ_m, k1_m))*dyh(JJJ_m)
                end do
              case (iBloqueMx)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                k1_m = K1
                k2_m = K2
                do JJJ = j1, j2
                  JJJ_m = JJJ
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (-Ey(i1_m, JJJ_m, k1_m) + Ey(i1_m, JJJ_m, k2_m + 1))*dye(JJJ_m)
                end do
                !--->
                i1_m = I1
                j1_m = J1
                j2_m = J2
                do KKK = k1, k2
                  KKK_m = KKK
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (Ez(i1_m, j1_m, KKK_m) - Ez(i1_m, j2_m + 1, KKK_m))*dze(KKK_m)
                end do
              case (iBloqueMy)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                i1_m = I1
                j1_m = J1
                i2_m = I2
                do KKK = k1, k2
                  KKK_m = KKK
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (Ez(i2_m + 1, j1_m, KKK_m) - Ez(i1_m, j1_m, KKK_m))*dze(KKK_m)
                end do
                !--->
                j1_m = J1
                k1_m = K1
                k2_m = K2
                do III = i1, i2
                  III_m = III
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (-Ex(III_m, j1_m, k2_m + 1) + Ex(III_m, j1_m, k1_m))*dxe(III_m)
                end do
              case (iBloqueMz)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                j1_m = J1
                k1_m = K1
                j2_m = J2
                do III = i1, i2
                  III_m = III
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (-Ex(III_m, j1_m, k1_m) + Ex(III_m, j2_m + 1, k1_m))*dxe(III_m)
                end do
                !--->
                i1_m = I1
                k1_m = K1
                i2_m = I2
                do JJJ = j1, j2
                  JJJ_m = JJJ
                  !--->
                  output(ii)%item(i)%valor(nTime - nInit) = &
                    output(ii)%item(i)%valor(nTime - nInit) + &
                    (Ey(i1_m, JJJ_m, k1_m) - Ey(i2_m + 1, JJJ_m, k1_m))*dye(JJJ_m)
                end do

              case (lineIntegral)
                block
                  integer(kind=4) :: lidx, lx, ly, lz, lor
                  type(direction_t), dimension(:), allocatable :: line
                  line = sgg%observation(ii)%P(i)%line
                  output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value

                  do lidx = 1, ubound(line, 1) - lbound(line, 1) + 1
                    lor = line(lidx)%orientation
                    lx = line(lidx)%x
                    ly = line(lidx)%y
                    lz = line(lidx)%z
                    select case (abs(lor))
                    case (iEx)
                      output(ii)%item(i)%valor(nTime - nInit) = output(ii)%item(i)%valor(nTime - nInit) + &
                                                                Ex(lx, ly, lz)*sign(1, lor)*dxe(lx)
                    case (iEy)
                      output(ii)%item(i)%valor(nTime - nInit) = output(ii)%item(i)%valor(nTime - nInit) + &
                                                                Ey(lx, ly, lz)*sign(1, lor)*dye(ly)
                    case (iEz)
                      output(ii)%item(i)%valor(nTime - nInit) = output(ii)%item(i)%valor(nTime - nInit) + &
                                                                Ez(lx, ly, lz)*sign(1, lor)*dze(lz)
                    end select
                  end do
                end block

              case (iQx, iQy, iQz)
                output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                SegmDumm => output(ii)%item(i)%Segmento
                output(ii)%item(i)%valor(nTime - nInit) = SegmDumm%ChargeMinus%ChargePresent

              case (iJx, iJy, iJz)
                if ((trim(adjustl(wiresflavor)) == 'holland') .or. &
                    (trim(adjustl(wiresflavor)) == 'transition')) then
                  output(ii)%item(i)%valor(nTime - nInit) = 0.0_RKIND !wipe value
                  SegmDumm => output(ii)%item(i)%Segmento
                  if (wirecrank) then !no hay que promediar nada porque estan co-locados en tiempo
                    output(ii)%item(i)%valor(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                                              SegmDumm%Currentpast
                    output(ii)%item(i)%valor2(nTime - nInit) = -SegmDumm%Efield_wire2main*SegmDumm%delta
                    output(ii)%item(i)%valor3(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                          (((SegmDumm%ChargePlus%ChargePresent)))*SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
                    output(ii)%item(i)%valor4(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                         (((SegmDumm%ChargeMinus%ChargePresent)))*SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
      output(ii)%item(i)%valor5(nTime - nInit) = output(ii)%item(i)%valor3(nTime - nInit) - output(ii)%item(i)%valor4(nTime - nInit)

                  else
!!saco el potencial calculado con E*delta !051115 !!y aniado el vdrop antinugo porque la Z se hacien bien con este 030719
                    output(ii)%item(i)%valor(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                                              SegmDumm%currentpast
                    output(ii)%item(i)%valor2(nTime - nInit) = -SegmDumm%Efield_wire2main*SegmDumm%delta
                    output(ii)%item(i)%valor3(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                               (((SegmDumm%ChargePlus%ChargePresent + SegmDumm%ChargePlus%ChargePast))/2.0_RKIND)* &
                                                               SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
                    output(ii)%item(i)%valor4(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                             (((SegmDumm%ChargeMinus%ChargePresent + SegmDumm%ChargeMinus%ChargePast))/2.0_RKIND)* &
                                                               SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
      output(ii)%item(i)%valor5(nTime - nInit) = output(ii)%item(i)%valor3(nTime - nInit) - output(ii)%item(i)%valor4(nTime - nInit)

                  end if
                        !!!!!!!!!!!!!!!!!!
                end if

#ifdef CompileWithBerengerWires
                if (trim(adjustl(wiresflavor)) == 'berenger') then
                  SegmDumm_Berenger => output(ii)%item(i)%Segmento_Berenger
                  output(ii)%item(i)%valor(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                                            SegmDumm_Berenger%Currentpast
                  output(ii)%item(i)%valor2(nTime - nInit) = -SegmDumm_Berenger%field*SegmDumm_Berenger%dl
                  output(ii)%item(i)%valor3(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                                  (((SegmDumm_Berenger%ChargePlus + SegmDumm_Berenger%ChargePlusPast))/2.0_RKIND)* &
                                                  SegmDumm_Berenger%L*(InvMu(SegmDumm_Berenger%imed)*InvEps(SegmDumm_Berenger%imed))
                  output(ii)%item(i)%valor4(nTime - nInit) = output(ii)%item(i)%valorsigno* &
                                                (((SegmDumm_Berenger%ChargeMinus + SegmDumm_Berenger%ChargeMinusPast))/2.0_RKIND)* &
                                                  SegmDumm_Berenger%L*(InvMu(SegmDumm_Berenger%imed)*InvEps(SegmDumm_Berenger%imed))
      output(ii)%item(i)%valor5(nTime - nInit) = output(ii)%item(i)%valor3(nTime - nInit) - output(ii)%item(i)%valor4(nTime - nInit)
                end if
#endif
#ifdef CompileWithSlantedWires
                if ((trim(adjustl(wiresflavor)) == 'slanted') .or. (trim(adjustl(wiresflavor)) == 'semistructured')) then !del wiresflavor
                  SegmDumm_Slanted => output(ii)%item(i)%Segmento_Slanted
                  output(ii)%item(i)%valor(nTime - nInit) = & !ojo: slanted ya los orienta bien y no hay que multiplicar por valorsigno
                    SegmDumm_Slanted%Currentpast
                  output(ii)%item(i)%valor2(nTime - nInit) = -SegmDumm_Slanted%field*SegmDumm_Slanted%dl
                  output(ii)%item(i)%valor3(nTime - nInit) = &
                    (((SegmDumm_Slanted%Voltage(iPlus)%ptr%Voltage + SegmDumm_Slanted%Voltage(iPlus)%ptr%VoltagePast))/2.0_RKIND)
                  output(ii)%item(i)%valor4(nTime - nInit) = &
                    (((SegmDumm_Slanted%Voltage(iMinus)%ptr%Voltage + SegmDumm_Slanted%Voltage(iMinus)%ptr%VoltagePast))/2.0_RKIND)
      output(ii)%item(i)%valor5(nTime - nInit) = output(ii)%item(i)%valor3(nTime - nInit) - output(ii)%item(i)%valor4(nTime - nInit)
                end if
#endif
                !Volumic probes
              case (iExC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) = Ex(III_m, JJJ_m, KKK_m)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
              case (iEyC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) = Ey(III_m, JJJ_m, KKK_m)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
              case (iEzC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) = Ez(III_m, JJJ_m, KKK_m)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
!por aqui voy con los i1t, j1t, k1t
              case (iHxC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) = Hx(III_m, JJJ_m, KKK_m)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
              case (iHyC)
                AT = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) = Hy(III_m, JJJ_m, KKK_m)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
              case (iHzC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) = Hz(III_m, JJJ_m, KKK_m)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
                !
              case (iMEC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) =  &
                          &    sqrt(Ex(III_m, JJJ_m, KKK_m)**2.0_RKIND + Ey(III_m, JJJ_m, KKK_m)**2.0_RKIND + &
                          Ez(III_m, JJJ_m, KKK_m)**2.0_RKIND)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
              case (iMHC)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    do KKK = k1, k2
                    if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                      k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                      KKK_m = KKK
                      do JJJ = j1, j2
                      if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                        j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                        JJJ_m = JJJ
                        do III = i1, i2
                        if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                          i1t = int(iii/output(ii)%item(i)%Xtrancos)
                          III_m = III
                          output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t) =  &
                          &    sqrt(Hx(III_m, JJJ_m, KKK_m)**2.0_RKIND + Hy(III_m, JJJ_m, KKK_m)**2.0_RKIND + &
                          Hz(III_m, JJJ_m, KKK_m)**2.0_RKIND)
                        end if
                        end do
                      end if
                      end do
                    end if
                    end do
                  end if
                end if
                !sondas corriente 20/2/14
              case (iCur, iCurX, iCurY, iCurZ, mapvtk)
                at = sgg%tiempo(ntime)
                if (at > sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done = .true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                Ntimeforvolumic = Ntime !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    conta = 0 !en el mismo orden que se allocatearon
                    do KKK = k1, k2
                      do JJJ = j1, j2
                        do III = i1, i2
                          !saca  current a lo largo del edge con las sondas icur
                          if (field /= mapvtk) then
                            ! refactoring done. Needs tests
                            do Efield = iEx, iEz
                              if (isThinWire(Efield, iii, jjj, kkk) .and. isWithinBounds(Efield, iii, jjj, kkk)) then
                                conta = conta + 1
                                jdir = computeJ(EField, iii, jjj, kkk)
                               output(ii)%item(i)%Serialized%valor_x(Ntimeforvolumic, conta) = merge(jdir, 0.0_RKIND, Efield == iEx)
                               output(ii)%item(i)%Serialized%valor_y(Ntimeforvolumic, conta) = merge(jdir, 0.0_RKIND, Efield == iEy)
                               output(ii)%item(i)%Serialized%valor_z(Ntimeforvolumic, conta) = merge(jdir, 0.0_RKIND, Efield == iEz)
                                             output( ii)%item( i)%Serialized%valor_Ex(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEx,Efield) 
                                             output( ii)%item( i)%Serialized%valor_Ey(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEy,Efield)
                                             output( ii)%item( i)%Serialized%valor_Ez(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEz,Efield)
                                             output( ii)%item( i)%Serialized%valor_Hx(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHx,Efield)
                                             output( ii)%item( i)%Serialized%valor_Hy(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHy,Efield)
                                             output( ii)%item( i)%Serialized%valor_Hz(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHz,Efield)
                              end if

                              if (.not. isMediaVacuum(Efield, iii, jjj, kkk) .and. &
                                  .not. isSplitOrAdvanced(Efield, iii, jjj, kkk) .and. &
                                  isWithinBounds(Efield, iii, jjj, kkk)) then
                                conta = conta + 1
                                jdir = computeJ(EField, iii, jjj, kkk)
                               output(ii)%item(i)%Serialized%valor_x(Ntimeforvolumic, conta) = merge(jdir, 0.0_RKIND, Efield == iEx)
                               output(ii)%item(i)%Serialized%valor_y(Ntimeforvolumic, conta) = merge(jdir, 0.0_RKIND, Efield == iEy)
                               output(ii)%item(i)%Serialized%valor_z(Ntimeforvolumic, conta) = merge(jdir, 0.0_RKIND, Efield == iEz)
                                             output( ii)%item( i)%Serialized%valor_Ex(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEx,Efield) 
                                             output( ii)%item( i)%Serialized%valor_Ey(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEy,Efield)
                                             output( ii)%item( i)%Serialized%valor_Ez(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEz,Efield)
                                             output( ii)%item( i)%Serialized%valor_Hx(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHx,Efield)
                                             output( ii)%item( i)%Serialized%valor_Hy(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHy,Efield)
                                             output( ii)%item( i)%Serialized%valor_Hz(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHz,Efield)
                              end if
                            end do

                          else !si es mapvtk

                            do Efield = iEx, iEz
                              call assignMedia(imed, imed1, imed2, imed3, imed4, Efield, iii, jjj, kkk)
                            call contabordes(sgg, imed, imed1, imed2, imed3, imed4, EsBorde, SINPML_fullsize, Efield, iii, jjj, kkk)
                              if (esBorde) then
                                conta = conta + 1
                                output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = &
                                  assignEdgeMediaType(Efield, iii, jjj, kkk)
                              end if
                            end do

                          end if
                        end do
                      end do
                    end do

                           !!!
                    if (field == mapvtk) then
                      INIT = .false.; geom = .false.; asigna = .true.; magnetic = .false.; electric = .true.
  call nodalvtk(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,media%sggMtag,tag_numbers, &
                                    init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)

        call wirebundlesvtk(sgg, init, geom, asigna, conta, i, ii, output, Ntimeforvolumic, wiresflavor, media%sggMtag, tag_numbers)
                    end if
                           !!!
                    do KKK = k1, k2
                      do JJJ = j1, j2
                        do III = i1, i2
                          !saca current en surfaces 0124
                          if (field /= mapvtk) then
                            do HField = iHx, iHz
                              if ((isPECorSurface(Hfield, iii, jjj, kkk) .or. field == blockCurrent(Hfield)) .and. &
                                  isWithinBounds(Hfield, iii, jjj, kkk)) then
                                conta = conta + 1
                                jdir1 = computeJ1(HField, iii, jjj, kkk)
                                jdir2 = computeJ2(HField, iii, jjj, kkk)

 output(ii)%item(i)%Serialized%valor_x(Ntimeforvolumic, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHz), Hfield == iHx)
 output(ii)%item(i)%Serialized%valor_y(Ntimeforvolumic, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHx), Hfield == iHy)
 output(ii)%item(i)%Serialized%valor_z(Ntimeforvolumic, conta) = merge(0.0_RKIND, merge(jdir1, jdir2, HField == iHy), Hfield == iHz)

                                                   output( ii)%item( i)%Serialized%valor_Ex(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEx,HField)
                                                   output( ii)%item( i)%Serialized%valor_Ey(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEy,HField)
                                                   output( ii)%item( i)%Serialized%valor_Ez(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iEz,HField)
                                                   output( ii)%item( i)%Serialized%valor_Hx(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHx,HField)
                                                   output( ii)%item( i)%Serialized%valor_Hy(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHy,HField)
                                                   output( ii)%item( i)%Serialized%valor_Hz(Ntimeforvolumic,conta) = interpolate_field_atwhere(sgg,Ex,Ey,Ez,Hx,Hy,Hz,iii, jjj, kkk, iHz,HField)
                              end if
                            end do
                          else
                            do Hfield = iHx, iHz
                              if (surfaceIsMedia(Hfield, iii, jjj, kkk)) then
                                conta = conta + 1
                                output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = &
                                  assignSurfaceMediaType(Hfield, iii, jjj, kkk)
                              end if
                              ! faces or edges?
                              if (tag_numbers%getFaceTag(Hfield, iii, jjj, kkk) < 0 .and. &
                                  (btest(iabs(tag_numbers%getFaceTag(Hfield, iii, jjj, kkk)), Hfield - 1)) .and. &
                                  .not. isPML(Hfield, iii, jjj, kkk) .and. isWithinBounds(Hfield, iii, jjj, kkk)) then
                                conta = conta + 1
                                call updateJ(Hfield, jdir)
                                output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = jdir
                              end if

                            end do
                          end if !del if mapvtk
                          !
                        end do
                      end do
                    end do

                           !!!
                    if (field == mapvtk) then
                      INIT = .false.; geom = .false.; asigna = .true.; magnetic = .true.; electric = .false.
  call nodalvtk(sgg,media%sggMiEx,media%sggMiEy,media%sggMiEz,media%sggMiHx,media%sggMiHy,media%sggMiHz,media%sggMtag,tag_numbers, &
                                    init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)
                    end if
                           !!!
                           !!!!!!!!!!!!esto dara problemas en los angulos y aristas donde porque ahi sacara la Bloque current en Hx!!!! 19/2/14
                  end if
                end if
                     !!!!!!!!fin sondas corriente
              case (FarField)
                call UpdateFarField(ntime, b, Ex, Ey, Ez, Hx, Hy, Hz)
              end select
                  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FREQMAIN
                  !!!!!!!!!!!!!!!!!!!!
                  !!!!!!!!!!!!!!!!!!!!
            elseif (SGG%Observation(ii)%FreqDomain) then
              at = sgg%tiempo(ntime)
!!!!! permit scaling

              do iff1 = 1, output(ii)%NumFreqs
                output(ii)%auxExp_E(iff1) = sgg%dt*(1.0E0_RKIND, 0.0E0_RKIND)*Exp(mcpi2*output(ii)%Freq(iff1)*at)   !el dt deberia ser algun tipo de promedio pero no me complico permit scaling 211118
                output(ii)%auxExp_H(iff1) = output(ii)%auxExp_E(iff1)*Exp(mcpi2*output(ii)%Freq(iff1)*sgg%dt*0.5_RKIND)
              end do
!!!
              select case (field)
!!!las freqdomain NUNCA ESTAN DONE
              case (iMEC, iExC, iEyC, iEzC) !como los tengo que guardar todas las componentes cartesianas solo variara el output final en el .xdmf
                !                    if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                do KKK = k1, k2
                if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                  k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                  KKK_m = KKK
                  do JJJ = j1, j2
                  if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                    j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                    JJJ_m = JJJ
                    do III = i1, i2
                    if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                      i1t = int(iii/output(ii)%item(i)%Xtrancos)
                      III_m = III
                      do if1 = 1, output(ii)%NumFreqs
                                 !!!
                        output(ii)%item(i)%valor3DComplex(if1, 1, i1t, j1t, k1t) = output(ii)%item(i)%valor3DComplex(if1, 1, i1t, j1t, k1t) + &
                                                                                   output(ii)%auxExp_E(if1)*Ex(III_m, JJJ_m, KKK_m)
                        output(ii)%item(i)%valor3DComplex(if1, 2, i1t, j1t, k1t) = output(ii)%item(i)%valor3DComplex(if1, 2, i1t, j1t, k1t) + &
                                                                                   output(ii)%auxExp_E(if1)*Ey(III_m, JJJ_m, KKK_m)
                        output(ii)%item(i)%valor3DComplex(if1, 3, i1t, j1t, k1t) = output(ii)%item(i)%valor3DComplex(if1, 3, i1t, j1t, k1t) + &
                                                                                   output(ii)%auxExp_E(if1)*Ez(III_m, JJJ_m, KKK_m)
                      end do
                    end if
                    end do
                  end if
                  end do
                end if
                end do
              case (iMHC, iHxC, iHyC, iHzC)
!                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                do KKK = k1, k2
                if (mod(KKK, output(ii)%item(i)%Ztrancos) == 0) then
                  k1t = int(kkk/output(ii)%item(i)%Ztrancos)
                  KKK_m = KKK
                  do JJJ = j1, j2
                  if (mod(jjj, output(ii)%item(i)%Ytrancos) == 0) then
                    j1t = int(jjj/output(ii)%item(i)%Ytrancos)
                    JJJ_m = JJJ
                    do III = i1, i2
                    if (mod(iii, output(ii)%item(i)%Xtrancos) == 0) then
                      i1t = int(iii/output(ii)%item(i)%Xtrancos)
                      III_m = III
                      do if1 = 1, output(ii)%NumFreqs
                                 !!!
             output(ii)%item(i)%valor3DComplex(if1, 1, i1t, j1t, k1t) = output(ii)%item(i)%valor3DComplex(if1, 1, i1t, j1t, k1t) + &
                                                                                   output(ii)%auxExp_H(if1)*Hx(III_m, JJJ_m, KKK_m)
             output(ii)%item(i)%valor3DComplex(if1, 2, i1t, j1t, k1t) = output(ii)%item(i)%valor3DComplex(if1, 2, i1t, j1t, k1t) + &
                                                                                   output(ii)%auxExp_H(if1)*Hy(III_m, JJJ_m, KKK_m)
             output(ii)%item(i)%valor3DComplex(if1, 3, i1t, j1t, k1t) = output(ii)%item(i)%valor3DComplex(if1, 3, i1t, j1t, k1t) + &
                                                                                   output(ii)%auxExp_H(if1)*Hz(III_m, JJJ_m, KKK_m)
                      end do
                    end if
                    end do
                  end if
                  end do
                end if
                end do
                !sondas corriente 20/2/14
              case (iCur, iCurX, iCurY, iCurZ)
!                     if (at > sgg%OBSERVATION(ii)%FinalTime+sgg%dt/2.0_RKIND) sgg%OBSERVATION(ii)%Done=.true.
                if (at >= sgg%OBSERVATION(ii)%InitialTime) sgg%OBSERVATION(ii)%Begun = .true.
                conta = 0 !en el mismo orden que se allocatearon
                do KKK = k1, k2
                  do JJJ = j1, j2
                    do III = i1, i2
                      !saca bul current a lo largo del edgje con las sondas icur
                      do Hfield = iHx, iHz
                        if (isThinWire(Hfield, iii, jjj, kkk) .and. isWithinBounds(Hfield, iii, jjj, kkk)) then
                          conta = conta + 1
                          jdir = computeJ(Hfield, iii, jjj, kkk)
                          do if1 = 1, output(ii)%NumFreqs
                            output(ii)%item(i)%Serialized%valorComplex_x(conta, if1) = &
              merge(output(ii)%item(i)%Serialized%valorComplex_x(conta, if1) + output(ii)%auxExp_H(if1)*jdir, z_cplx, Hfield == iHx)
                            output(ii)%item(i)%Serialized%valorComplex_y(conta, if1) = &
              merge(output(ii)%item(i)%Serialized%valorComplex_y(conta, if1) + output(ii)%auxExp_H(if1)*jdir, z_cplx, Hfield == iHy)
                            output(ii)%item(i)%Serialized%valorComplex_z(conta, if1) = &
              merge(output(ii)%item(i)%Serialized%valorComplex_z(conta, if1) + output(ii)%auxExp_H(if1)*jdir, z_cplx, Hfield == iHz)

                            output(ii)%item(i)%Serialized%valorComplex_Ex(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ex(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEx, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Ey(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ey(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEy, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Ez(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ez(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEz, Hfield)

                            output(ii)%item(i)%Serialized%valorComplex_Hx(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hx(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHx, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Hy(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hy(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHy, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Hz(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hz(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHz, Hfield)
                          end do
                        end if

                        if (.not. isMediaVacuum(Hfield, iii, jjj, kkk) .and. &
                            .not. isSplitOrAdvanced(Hfield, iii, jjj, kkk) .and. &
                            isWithinBounds(Hfield, iii, jjj, kkk)) then
                          conta = conta + 1
                          jdir = computeJ(Hfield, iii, jjj, kkk)
                          do if1 = 1, output(ii)%NumFreqs
                            output(ii)%item(i)%Serialized%valorComplex_x(conta, if1) = &
              merge(output(ii)%item(i)%Serialized%valorComplex_x(conta, if1) + output(ii)%auxExp_H(if1)*jdir, z_cplx, Hfield == iHx)
                            output(ii)%item(i)%Serialized%valorComplex_y(conta, if1) = &
              merge(output(ii)%item(i)%Serialized%valorComplex_y(conta, if1) + output(ii)%auxExp_H(if1)*jdir, z_cplx, Hfield == iHy)
                            output(ii)%item(i)%Serialized%valorComplex_z(conta, if1) = &
              merge(output(ii)%item(i)%Serialized%valorComplex_z(conta, if1) + output(ii)%auxExp_H(if1)*jdir, z_cplx, Hfield == iHz)

                            output(ii)%item(i)%Serialized%valorComplex_Ex(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ex(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEx, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Ey(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ey(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEy, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Ez(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ez(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEz, Hfield)

                            output(ii)%item(i)%Serialized%valorComplex_Hx(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hx(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHx, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Hy(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hy(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHy, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Hz(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hz(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHz, Hfield)
                          end do
                        end if
                      end do
                    end do
                  end do
                end do
                do KKK = k1, k2
                  do JJJ = j1, j2
                    do III = i1, i2
                      !saca current en surfaces 0124
                      do HField = iHx, iHz
                        if ((isPECorSurface(Hfield, iii, jjj, kkk) .or. &
                             field == blockCurrent(Hfield)) .and. &
                            isWithinBounds(Hfield, iii, jjj, kkk)) then

                          conta = conta + 1
                          jdir1 = computeJ1(HField, iii, jjj, kkk)
                          jdir2 = computeJ2(HField, iii, jjj, kkk)
                          do if1 = 1, output(ii)%NumFreqs

                            output(ii)%item(i)%Serialized%valorComplex_x(conta, if1) = &
                                             merge(z_cplx, output( ii)%item( i)%Serialized%valorComplex_x(conta, if1)+output( ii)%auxExp_H(if1)*merge(jdir1, jdir2, HField == iHz), Hfield == iHx)
                            output(ii)%item(i)%Serialized%valorComplex_y(conta, if1) = &
                                             merge(z_cplx, output( ii)%item( i)%Serialized%valorComplex_y(conta, if1)+output( ii)%auxExp_H(if1)*merge(jdir1, jdir2, HField == iHx), Hfield == iHy)
                            output(ii)%item(i)%Serialized%valorComplex_z(conta, if1) = &
                                             merge(z_cplx, output( ii)%item( i)%Serialized%valorComplex_z(conta, if1)+output( ii)%auxExp_H(if1)*merge(jdir1, jdir2, HField == iHy), Hfield == iHz)

                            output(ii)%item(i)%Serialized%valorComplex_Ex(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ex(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEx, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Ey(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ey(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEy, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Ez(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Ez(conta, if1) + &
                         output(ii)%auxExp_E(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iEz, Hfield)

                            output(ii)%item(i)%Serialized%valorComplex_Hx(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hx(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHx, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Hy(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hy(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHy, Hfield)
                            output(ii)%item(i)%Serialized%valorComplex_Hz(conta, if1) = &
                              output(ii)%item(i)%Serialized%valorComplex_Hz(conta, if1) + &
                         output(ii)%auxExp_H(if1)*interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, iii, jjj, kkk, iHz, Hfield)
                          end do
                        end if
                      end do
                    end do
                  end do
                end do
                     !!!!!!!!!!!!esto dara problemas en los angulos y aristas donde porque ahi sacara la Bloque current en Hx!!!! 19/2/14

                     !!!!!!!!fin sondas corriente
              end select

            end if !time domain
          end if !del nothing
        end do loop_obser
      end do
      !---------------------------> acaba UpdateObservation <-----------------------------------------
      return

    contains

      logical function isPECorSurface(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isPECorSurface = sgg%med(media)%is%PEC .or. &
                         sgg%med(media)%is%Surface
      end function

      integer function blockCurrent(field)
        integer(kind=4) :: field
        select case (field)
        case (iHx)
          blockCurrent = iCurX
        case (iHy)
          blockCurrent = iCurY
        case (iHz)
          blockCurrent = iCurZ
        end select
      end function

      function computeJ1(f, i, j, k) result(res)
        integer(kind=4) :: f, i, j, k, c
        real(kind=rkind) :: res

        c = mod(f - 2, 3) + 4

         res = - ((getDelta(c, i, j, k)* getField(c,i,j,k)                               + getDelta(c, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx))*getField(c, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx))) - &
                  (getDelta(c, i, j, k)* getField(c,i-u(f,iHx),j-u(f,iHy),k-u(f,iHz))    + getDelta(c, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx))*getField(c, i-u(f,iHx)+u(f,iHy), j-u(f,iHy)+u(f,iHz), k-u(f,iHz)+u(f,iHx))) + &
                   getDelta(f, i, j, k)*(getField(f,i-u(f,iHy),j-u(f,iHz),k-u(f,iHx))    -                                                 getField(f, i+u(f,iHy), j+u(f,iHz), k+u(f,iHx))))

      end function
      function computeJ2(f, i, j, k) result(res)
        integer(kind=4) :: f, i, j, k, c
        real(kind=rkind) :: res

        c = mod(f - 3, 3) + 4

         res =   (getDelta(c, i, j, k)* getField(c,i,j,k)                               + getDelta(c, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy))*getField(c, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy))) - &
                 (getDelta(c, i, j, k)* getField(c,i-u(f,iHx),j-u(f,iHy),k-u(f,iHz))    + getDelta(c, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy))*getField(c, i-u(f,iHx)+u(f,iHz),j-u(f,iHy)+u(f,iHx),k-u(f,iHz)+u(f,iHy))) + &
                  getDelta(f, i, j, k)*(getField(f,i-u(f,iHz),j-u(f,iHx),k-u(f,iHy))    -                                                 getField(f, i+u(f,iHz), j+u(f,iHx), k+u(f,iHy)))
      end function

      integer function u(field1, field2)
        integer(kind=4) :: field1, field2
        if (field1 == field2) then
          u = 1
        else
          u = 0
        end if
      end function

      logical function isSplitOrAdvanced(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isSplitOrAdvanced = sgg%med(media)%is%split_and_useless .or. &
                            sgg%med(media)%is%already_YEEadvanced_byconformal

      end function

      function getDelta(field, i, j, k) result(res)
        integer(kind=4) :: field, i, j, k
        real(kind=rkind) :: res
        select case (field)
        case (iEx, iHx)
          res = dxh(i)
        case (iEy, iHy)
          res = dyh(j)
        case (iEz, iHz)
          res = dzh(k)
        end select
      end function

      function computeJ(field, i, j, k) result(res)
        integer(kind=4) :: field, i, j, k, i1, i2, j1, j2, k1, k2
        real(kind=rkind) :: res
        i1 = i - merge(1, 0, 1 + mod(field + 1, 3) == iEx)
        j1 = j - merge(1, 0, 1 + mod(field + 1, 3) == iEy)
        k1 = k - merge(1, 0, 1 + mod(field + 1, 3) == iEz)
        i2 = i - merge(1, 0, 1 + mod(field, 3) == iEx)
        j2 = j - merge(1, 0, 1 + mod(field, 3) == iEy)
        k2 = k - merge(1, 0, 1 + mod(field, 3) == iEz)
 res =  getDelta(1+mod(field,3)  , i, j, k) * (-getField(1+mod(field,3) + 3,  i,j,k) + getField(1+mod(field,3) + 3,  i1,j1,k1))  + &
             getDelta(1+mod(field+1,3), i, j, k) * ( getField(1+mod(field+1,3) + 3,i,j,k) - getField(1+mod(field+1,3) + 3,i2,j2,k2))
      end function

      function getField(field, i, j, k) result(res)
        real(kind=rkind) :: res
        integer(kind=4) :: field, i, j, k
        select case (field)
        case (iEx)
          res = Ex(i, j, k)
        case (iEy)
          res = Ey(i, j, k)
        case (iEz)
          res = Ez(i, j, k)
        case (iHx)
          res = Hx(i, j, k)
        case (iHy)
          res = Hy(i, j, k)
        case (iHz)
          res = Hz(i, j, k)
        end select
      end function

      function getMedia(field, i, j, k) result(res)
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: res
        integer(kind=4) :: field, i, j, k
        select case (field)
        case (iEx)
          res = media%sggMiEx(i, j, k)
        case (iEy)
          res = media%sggMiEy(i, j, k)
        case (iEz)
          res = media%sggMiEz(i, j, k)
        case (iHx)
          res = media%sggMiHx(i, j, k)
        case (iHy)
          res = media%sggMiHy(i, j, k)
        case (iHz)
          res = media%sggMiHz(i, j, k)
        end select
      end function

      logical function isMediaVacuum(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media, vacuum = 1
        media = getMedia(field, i, j, k)
        isMediaVacuum = (media == vacuum)
      end function

      logical function isWithinBounds(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        isWithinBounds = (i <= SINPML_fullsize(field)%XE) .and. &
                         (j <= SINPML_fullsize(field)%YE) .and. &
                         (k <= SINPML_fullsize(field)%ZE)
      end function

      logical function isThinWire(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isThinWire = sgg%Med(media)%is%ThinWire
      end function

      logical function isPML(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        isPML = sgg%med(media)%is%PML
      end function

      logical function isSGBCorMultiport(media)
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        isSGBCorMultiport = (sgg%Med(media)%is%SGBC) .or. &
                            (sgg%Med(media)%is%multiport) .or. &
                            (sgg%Med(media)%is%anismultiport)
      end function

      logical function isDispersive(media)
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        isDispersive = (sgg%Med(media)%is%edispersive) .or. &
                       (sgg%Med(media)%is%EDispersiveANIS) .or. &
                       (sgg%Med(media)%is%mDispersive) .or. &
                       (sgg%Med(media)%is%mDispersiveANIS)
      end function

      logical function surfaceIsMedia(field, i, j, k)
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: surfaceMedia, vacuum = 1
        integer(kind=4) :: field, i, j, k
        surfaceMedia = getMedia(field, i, j, k)
        surfaceIsMedia = ((surfaceMedia /= vacuum) .and. &
                          .not. sgg%med(surfaceMedia)%is%PML .and. &
                          isWithinBounds(field, i, j, k))
      end function

      function assignSurfaceMediaType(field, i, j, k) result(res)
        integer(kind=4) :: field, i, j, k
        real(kind=RKIND) :: res
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        if ((media == 0) .or. (sgg%Med(media)%is%Pec)) then
          res = 0
        elseif (sgg%Med(media)%is%thinwire) then
          CALL StopOnError(0, 1, 'ERROR: A magnetic field cannot be a thin-wire')
        elseif (isSGBCorMultiport(media)) then
          res = 300 + media
        elseif (isDispersive(media)) then
          res = 100 + media
        elseif ((sgg%Med(media)%is%Dielectric) .or. &
                (sgg%Med(media)%is%Anisotropic)) then
          res = 200 + media
        elseif (sgg%Med(media)%is%thinslot) then
          res = 400 + media
        elseif ((sgg%Med(media)%is%already_YEEadvanced_byconformal) .and. (.not. noconformalmapvtk)) then
          res = 5
        elseif ((sgg%Med(media)%is%split_and_useless) .and. (.not. noconformalmapvtk)) then
          res = 6
        else
          res = -1
        end if
      end function

      function assignEdgeMediaType(field, i, j, k) result(res)
        integer(kind=4) :: field, i, j, k
        real(kind=RKIND) :: res
        integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: media
        media = getMedia(field, i, j, k)
        if ((sgg%Med(media)%is%already_YEEadvanced_byconformal) .and. (.not. noconformalmapvtk)) then
          res = 5.5
        elseif ((sgg%Med(media)%is%split_and_useless) .and. (.not. noconformalmapvtk)) then
          res = 6.5
        elseif (sgg%Med(media)%is%thinwire) then
          if (collidesWithNonThinWire(field, i, j, k)) then
            res = 8
          else
            res = 7
          end if
        elseif ((media == 0) .or. (sgg%Med(media)%is%Pec)) then
          res = 0.5_RKIND
        elseif (isSGBCorMultiport(media)) then
          res = 3.5
        elseif (isDispersive(media)) then
          res = 1.5
        elseif ((sgg%Med(media)%is%Dielectric) .or. &
                (sgg%Med(media)%is%Anisotropic)) then
          res = 2.5
        elseif (sgg%Med(media)%is%thinslot) then
          res = 4.5
        else
          res = -0.5_RKIND
        end if

      end function

      subroutine assignMedia(m, m1, m2, m3, m4, dir, i, j, k)
        integer(kind=4), intent(inout) :: m, m1, m2, m3, m4
        integer(kind=4) :: dir, i, j, k
        m = getMedia(dir, i, j, k)
        m1 = getMedia(4 + modulo(dir, 3), i, j, k)
        m2 = getMedia(4 + modulo(dir, 3), i - merge(1, 0, dir == iEy), j - merge(1, 0, dir == iEz), k - merge(1, 0, dir == iEx))
        m3 = getMedia(4 + modulo(dir + 1, 3), i, j, k)
        m4 = getMedia(4 + modulo(dir + 1, 3), i - merge(1, 0, dir == iEz), j - merge(1, 0, dir == iEx), k - merge(1, 0, dir == iEy))
      end subroutine

      logical function collidesWithNonThinWire(field, i, j, k)
        integer(kind=4) :: field, i, j, k
        integer(kind=4) :: idx(6)

        collidesWithNonThinWire = .false.
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field,3),i,j,k)/=1 .and. .not.   sgg%med(getMedia(1+mod(field,3),i,j,k))%is%thinWire)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field+1,3),i,j,k)/=1 .and. .not. sgg%med(getMedia(1+mod(field+1,3),i,j,k))%is%thinWire)
        idx = assignIndices1(field, i, j, k)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field,3),idx(1),idx(2),idx(3))/=1 .and. .not.   sgg%med(getMedia(1+mod(field,3),idx(1),idx(2),idx(3)))%is%thinWire)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field+1,3),idx(4),idx(5),idx(6))/=1 .and. .not. sgg%med(getMedia(1+mod(field+1,3),idx(4),idx(5),idx(6)))%is%thinWire)
        idx = assignIndices2(field, i, j, k)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field,3),idx(1),idx(2),idx(3))/=1 .and. .not.   sgg%med(getMedia(1+mod(field,3),idx(1),idx(2),idx(3)))%is%thinWire)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field+1,3),idx(4),idx(5),idx(6))/=1 .and. .not. sgg%med(getMedia(1+mod(field+1,3),idx(4),idx(5),idx(6)))%is%thinWire)
        idx = assignIndices3(field, i, j, k)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field,3),idx(1),idx(2),idx(3))/=1 .and. .not.   sgg%med(getMedia(1+mod(field,3),idx(1),idx(2),idx(3)))%is%thinWire)
      collidesWithNonThinWire = collidesWithNonThinWire .or. (getMedia(1+mod(field+1,3),idx(4),idx(5),idx(6))/=1 .and. .not. sgg%med(getMedia(1+mod(field+1,3),idx(4),idx(5),idx(6)))%is%thinWire)

      end function

      function assignIndices1(field, i, j, k) result(res)
        integer(kind=4), intent(in)  :: field, i, j, k
        integer(kind=4) :: res(6)
        res(1) = i - merge(1, 0, 1 + mod(field, 3) == iEx)
        res(2) = j - merge(1, 0, 1 + mod(field, 3) == iEy)
        res(3) = k - merge(1, 0, 1 + mod(field, 3) == iEz)
        res(4) = i - merge(1, 0, 1 + mod(field + 1, 3) == iEx)
        res(5) = j - merge(1, 0, 1 + mod(field + 1, 3) == iEy)
        res(6) = k - merge(1, 0, 1 + mod(field + 1, 3) == iEz)
      end function
      function assignIndices2(field, i, j, k) result(res)
        integer(kind=4), intent(in)  :: field, i, j, k
        integer(kind=4) :: res(6)
        res(1) = i + merge(1, 0, field == iEx)
        res(2) = j + merge(1, 0, field == iEy)
        res(3) = k + merge(1, 0, field == iEz)
        res(4) = i + merge(1, 0, field == iEx)
        res(5) = j + merge(1, 0, field == iEy)
        res(6) = k + merge(1, 0, field == iEz)
      end function
      function assignIndices3(field, i, j, k) result(res)
        integer(kind=4), intent(in)  :: field, i, j, k
        integer(kind=4) :: res(6)
        res(1) = i + merge(1, 0, field == iEx) - merge(1, 0, 1 + mod(field, 3) == iEx)
        res(2) = j + merge(1, 0, field == iEy) - merge(1, 0, 1 + mod(field, 3) == iEy)
        res(3) = k + merge(1, 0, field == iEz) - merge(1, 0, 1 + mod(field, 3) == iEz)
        res(4) = i + merge(1, 0, field == iEx) - merge(1, 0, 1 + mod(field + 1, 3) == iEx)
        res(5) = j + merge(1, 0, field == iEy) - merge(1, 0, 1 + mod(field + 1, 3) == iEy)
        res(6) = k + merge(1, 0, field == iEz) - merge(1, 0, 1 + mod(field + 1, 3) == iEz)
      end function

      subroutine updateJ(field, jdir)
        integer(kind=4) :: field
        real(kind=rkind) :: jdir
        select case (field)
        case (iEx)
          Jx = -100
          jdir = jx
        case (iEy)
          Jy = -100
          jdir = jy
        case (iEz)
          Jz = -100
          jdir = jz
        end select
      end subroutine

    end subroutine UpdateObservation

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Flushes the observed magnitudes to disk
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine FlushObservationFiles(sgg,nInit,FinalInstant,layoutnumber,size, dxe,dye,dze,dxh,dyh,dzh,b,singlefilewrite,facesNF2FF,flushff)
      USE ILUMINA !is needed to also calculate the incident field in the observed points
      !solo lo precisa de entrada farfield
      type(bounds_t)  ::  b
      !
      type(nf2ff_t) :: facesNF2FF
      !!!
      !
      type(SGGFDTDINFO), intent(IN)         ::  sgg
      REAL(KIND=RKIND), dimension(:), intent(in)   :: dxh(sgg%ALLOC(iEx)%XI:sgg%ALLOC(iEx)%XE), &
                                                      dyh(sgg%ALLOC(iEy)%YI:sgg%ALLOC(iEy)%YE), &
                                                      dzh(sgg%ALLOC(iEz)%ZI:sgg%ALLOC(iEz)%ZE), &
                                                      dxe(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE), &
                                                      dye(sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE), &
                                                      dze(sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)
      integer(kind=4), intent(in) :: layoutnumber, size
      integer(kind=4)  ::  nInit, FinalInstant, unidad, compo, conta
      integer(kind=4)  ::  i, field, N, ii, i1, j1, k1, Ntimeforvolumic, dummy_jjj, i1t, j1t, k1t, i0t
      logical  ::  incident, singlefilewrite, flushff, ISyaopen
      REAL(KIND=RKIND_tiempo)  ::  at

      logical :: ok
      logical :: called_fromobservation, dummy_logical
      integer :: my_iostat

      character(LEN=BUFSIZE)  ::  whoami
      !!!
      write (whoami, '(a,i5,a,i5,a)') '(', layoutnumber + 1, '/', size, ') '
      called_fromobservation = .true.
      !!!ojo dummy_logical lo dejo que incid( lo fije a still_planewave_time sin tocarlo aqui para no afectar al principal 210419
      dummy_jjj = 1  !no es preciso fijarlo, porque a incid( se le pasa  called_fromobservation
      !Write also the incident fields in case there are plane waves (useful in SE calculations)
      incident = .false.
      if (sgg%NumPlaneWaves >= 1) incident = .true.
      do ii = 1, sgg%NumberRequest
        loop_obb: do i = 1, sgg%Observation(ii)%nP
          field = sgg%observation(ii)%P(i)%what
          if (SGG%Observation(ii)%TimeDomain) then
            select case (field)
            case (iEx, iEy, iEz, iJx, iJy, iJz, iQx, iQy, iQz, iHx, iHy, iHz, lineIntegral)
              I1 = sgg%OBSERVATION(ii)%P(i)%XI
              J1 = sgg%OBSERVATION(ii)%P(i)%YI
              K1 = sgg%OBSERVATION(ii)%P(i)%ZI
              !
              !                nInit=max(nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt),nInit)
              DO N = nInit, FinalInstant !!! bug octubre'14 mur1.nfde. Quitado --->>>,output(ii)%Trancos  !save only for the requested data
                if (mod(n, output(ii)%Trancos) == 0) then  !save only for the requested data
                  select case (field)
                  case (iHx, iHy, iHz)
                    at = sgg%tiempo(N) + sgg%dt/2.0_RKIND_tiempo
                  case (iEx, iEy, iEz, iJx, iJy, iJz, iQx, iQy, iQz, lineIntegral)
                    at = sgg%tiempo(N)
                  end select
                 if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND)) .or. &
                      output(ii)%saveall) then

                    select case (field)
                    case (iEx, iEy, iEz)
                      !
                      if (singlefilewrite) then
                        unidad = output(ii)%item(i)%unitmaster
                        if (incident) then
                          WRITE (unidad) output(ii)%item(i)%unit, at, output(ii)%item(i)%valor(n - nInit), &
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                          !by hand para clavarlo
                        else
                          WRITE (unidad) output(ii)%item(i)%unit, at, output(ii)%item(i)%valor(n - nInit)
                        end if
                      else
                        unidad = output(ii)%item(i)%unit
                        if (incident) then
#ifdef CompileWithReal16
                          WRITE (unidad, *) at, output(ii)%item(i)%valor(n - nInit), &
Incid(sgg, dummy_jjj, field, at*1.0_RKIND + 0.0_RKIND*sgg%dt, i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
#ifdef CompileWithmMiguelStandaloneObservation
                          WRITE (unidad, *) at, output(ii)%item(i)%valor(n - nInit), &
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
                          WRITE (unidad, fmt) at, output(ii)%item(i)%valor(n - nInit), &
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                          !by hand para clavarlo
#endif
#endif
                        else
                          WRITE (unidad, fmt) at, output(ii)%item(i)%valor(n - nInit)
                        end if
                      end if
                      !
                    case (iHx, iHy, iHz)
                      !
                      if (singlefilewrite) then
                        unidad = output(ii)%item(i)%unitmaster
                        if (incident) then
WRITE (unidad) output(ii)%item(i)%unit, at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                          !by hand para clavarlo
                        else
                          WRITE (unidad) output(ii)%item(i)%unit, at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
                        end if
                      else
                        unidad = output(ii)%item(i)%unit
                        if (incident) then
#ifdef CompileWithReal16
                          WRITE (unidad, *) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
#ifdef CompileWithmMiguelStandaloneObservation
                          WRITE (unidad, *) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
#else
                          WRITE (unidad, fmt) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit), &  !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
Incid(sgg, dummy_jjj, field, real(at + 0.0_RKIND*sgg%dt, RKIND), i1, j1, k1, dummy_logical, called_fromobservation) !quitado el 3 de ORIGINAL sync para pscale bien sincronizado
                          !by hand para clavarlo
#endif
#endif
                        else
                          WRITE (unidad, fmt) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
                        end if
                      end if
                      !
                    case (iJx, iJy, iJz)
                      if (singlefilewrite) then
                        unidad = output(ii)%item(i)%unitmaster
                        WRITE (unidad) output(ii)%item(i)%unit, at, &
                          output(ii)%item(i)%valor(n - nInit), &
                          output(ii)%item(i)%valor2(n - nInit), & !saco el valor2 -e*dl
                          output(ii)%item(i)%valor3(n - nInit), & ! VpluS
                          output(ii)%item(i)%valor4(n - nInit), & ! Vminus
                          output(ii)%item(i)%valor5(n - nInit) ! vplus-vminus
                      else
                        unidad = output(ii)%item(i)%unit
                        WRITE (unidad, fmt) at, output(ii)%item(i)%valor(n - nInit), &
                          output(ii)%item(i)%valor2(n - nInit), & !saco el valor2 -e*dl
                          output(ii)%item(i)%valor3(n - nInit), & ! VPLUS
                          output(ii)%item(i)%valor4(n - nInit), & ! Vminus
                          output(ii)%item(i)%valor5(n - nInit) ! vplus-vminus
                      end if
                    case (iQx, iQy, iQz)
                      if (singlefilewrite) then
                        unidad = output(ii)%item(i)%unitmaster
                        WRITE (unidad) output(ii)%item(i)%unit, at, &
                          output(ii)%item(i)%valor(n - nInit) ! node charge
                      else
                        unidad = output(ii)%item(i)%unit
                        WRITE (unidad, fmt) at, output(ii)%item(i)%valor(n - nInit) ! node charge
                      end if

                    case (lineIntegral)
                      if (singlefilewrite) then
                        unidad = output(ii)%item(i)%unitmaster
                        WRITE (unidad) output(ii)%item(i)%unit, at, &
                          output(ii)%item(i)%valor(n - nInit) ! e*dl sum along line
                      else
                        unidad = output(ii)%item(i)%unit
                        WRITE (unidad, fmt) at, output(ii)%item(i)%valor(n - nInit) ! e*dl sum along line
                      end if
                    end select
                  end if
                end if
              END DO
              if (singlefilewrite .and. ((field == iEx) .or. (field == iEy) .or. (field == iEz) .or. &
                                         (field == iVx) .or. (field == iVy) .or. (field == iVz) .or. &
                                         (field == iJx) .or. (field == iJy) .or. (field == iJz) .or. &
                                         (field == iQx) .or. (field == iQy) .or. (field == iQz) .or. &
                                         (field == iHx) .or. (field == iHy) .or. (field == iHz))) then
                unidad = output(ii)%item(i)%unitmaster
              else
                unidad = output(ii)%item(i)%unit
              end if
              call flush (unidad)
              !
            case (iBloqueMx, iBloqueMz, iBloqueMy, iBloqueJx, iBloqueJz, iBloqueJy)
#ifdef CompileWithMPI
              if ((field == iBloqueMx) .or. (field == iBloqueMy) .or. (field == iBloqueJx) .or. (field == iBloqueJy)) then
                if (output(ii)%item(i)%MPISubComm /= -1) then !solo si alguien tiene que hacerlo
                  valores(0:BuffObse) = output(ii)%item(i)%valor(0:BuffObse)
                  call MPIupdateBloques(layoutnumber, valores, newvalores, &
                                        output(ii)%item(i)%MPISubComm)
                  output(ii)%item(i)%valor(0:BuffObse) = newvalores(0:BuffObse)
                end if
              end if

              if ((layoutnumber == output(ii)%item(i)%MPIRoot) .or. &
                  (field == iBloqueJz) .or. (field == iBloqueMz)) then !only the master
#endif
                !                    nInit=max(nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt),nInit)
                DO N = nInit, FinalInstant
                  if (mod(n, output(ii)%Trancos) == 0) then  !save only for the requested data
                    select case (field)
                    case (iBloqueMx, iBloqueMz, iBloqueMy)
                      at = sgg%tiempo(N) + sgg%dt/2.0_RKIND_tiempo
                    case (iBloqueJx, iBloqueJz, iBloqueJy)
                      at = sgg%tiempo(N)
                    end select
                 if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND)) .or. &
                        output(ii)%saveall) then
                      select case (field)
                      case (iBloqueMx, iBloqueMz, iBloqueMy)
#ifdef CompileWithReal16
                        WRITE (output(ii)%item(i)%unit, *) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
                      case (iBloqueJx, iBloqueJz, iBloqueJy)
                        WRITE (output(ii)%item(i)%unit, *) at, output(ii)%item(i)%valor(n - nInit)
#else
#ifdef CompileWithmMiguelStandaloneObservation
                        WRITE (output(ii)%item(i)%unit, *) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
                      case (iBloqueJx, iBloqueJz, iBloqueJy)
                        WRITE (output(ii)%item(i)%unit, *) at, output(ii)%item(i)%valor(n - nInit)
#else
                        WRITE (output(ii)%item(i)%unit, fmt) at - sgg%dt/2.0_RKIND_tiempo, output(ii)%item(i)%valor(n - nInit) !SOLO A EFECTOS DE SALIDA EN FICHERO CHAPUZ SGG MAIL OLD 070916
                      case (iBloqueJx, iBloqueJz, iBloqueJy)
                        WRITE (output(ii)%item(i)%unit, fmt) at, output(ii)%item(i)%valor(n - nInit)
#endif
#endif
                      end select
                    end if
                  end if
                END DO
                call flush (output(ii)%item(i)%unit)
                !--->

                !--->
#ifdef CompileWithMPI
              end if
#endif

            case (FarField) !no emplear tiempo calculando rcs por el camino solo al final
              at = sgg%tiempo(FinalInstant)
              if (flushFF) call FlushFarfield(layoutnumber, size, b, dxe, dye, dze, dxh, dyh, dzh, facesNF2FF, at)
            case (iMHC, iHxC, iHyC, iHzC, iMEC, iExC, iEyC, iEzC, icur, iCurX, iCurY, iCurZ, mapvtk)
              DO N = nInit, FinalInstant
                at = sgg%tiempo(N)
                Ntimeforvolumic = N !!!-nint(0.4999999+sgg%OBSERVATION(ii)%InitialTime/sgg%dt)
                if (mod(Ntimeforvolumic, output(ii)%Trancos) == 0) then
                  Ntimeforvolumic = Ntimeforvolumic/output(ii)%Trancos
                  if (((at >= sgg%OBSERVATION(ii)%InitialTime) .and. (at <= sgg%OBSERVATION(ii)%FinalTime + sgg%dt/2.0_RKIND))) then
                    !assumo que todos son electricos o magneticos en una probe Volumic para calcular el tiempo !logico
                    output(ii)%TimesWritten = output(ii)%TimesWritten + 1
                    select case (field)
                    case (iMHC, iHxC, iHyC, iHzC, iMEC, iExC, iEyC, iEzC)
                      write (output(ii)%item(i)%unit) at
                      do k1t = output(ii)%item(i)%ZItrancos, output(ii)%item(i)%ZEtrancos
                        do j1t = output(ii)%item(i)%YItrancos, output(ii)%item(i)%YEtrancos
                          write (output(ii)%item(i)%unit) (output(ii)%item(i)%valor3D(Ntimeforvolumic, i1t, j1t, k1t), &
                          &                                i1t=output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos)
                        end do
                      end do
                    case (icur, iCurX, iCurY, iCurZ, mapvtk)
                      write (output(ii)%item(i)%unit) at
                      if (output(ii)%item(i)%columnas /= 0) then
                        do i1 = 1, output(ii)%item(i)%columnas
                          write (output(ii)%item(i)%unit) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_x(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_y(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_z(Ntimeforvolumic, i1)
                          write (output(ii)%item(i)%unit) output(ii)%item(i)%Serialized%valorE(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_Ex(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_Ey(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_Ez(Ntimeforvolumic, i1)
                          write (output(ii)%item(i)%unit) output(ii)%item(i)%Serialized%valorH(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_Hx(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_Hy(Ntimeforvolumic, i1), &
                            output(ii)%item(i)%Serialized%valor_Hz(Ntimeforvolumic, i1)
                        end do
                      end if
                    end select
                  end if
                end if
              end do
                  !!!!!!!!!!!!!!!!!!!
              !
              call flush (output(ii)%item(i)%unit)
            end select

          elseif (SGG%Observation(ii)%FreqDomain) then !only volumic probes are handled in freq domain in this way
            select case (field)
            case (iMHC, iHxC, iHyC, iHzC, iMEC, iExC, iEyC, iEzC)
              at = sgg%tiempo(FinalInstant)
                  !!!!assumo que todos son electricos o magneticos en una probe Volumic para calcular el tiempo !logico
              output(ii)%TimesWritten = output(ii)%NumFreqs  !util para leer el numero exacto de freq points
              INQUIRE (file=trim(adjustl(output(ii)%item(i)%path)), OPENED=ISyaopen)
              if (isyaopen) close (output(ii)%item(i)%unit, status='delete')
              !
              my_iostat = 0
9244          if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9244,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
         open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted', err=9244, iostat=my_iostat)
              write (output(ii)%item(i)%unit) &
                output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos, &
                output(ii)%item(i)%YItrancos, output(ii)%item(i)%YEtrancos, &
                output(ii)%item(i)%ZItrancos, output(ii)%item(i)%ZEtrancos
                !!! &      sgg%observation(ii)%P(i)%xI,sgg%observation(ii)%P(i)%xE, &
                !!! &      sgg%observation(ii)%P(i)%YI,sgg%observation(ii)%P(i)%YE, &
                !!! &      sgg%observation(ii)%P(i)%zI,sgg%observation(ii)%P(i)%ZE
              write (output(ii)%item(i)%unit) at !deteccion errores dft si se resumea a partir de instantes posteriores al ultimo escrito
              DO N = 1, output(ii)%NumFreqs
                write (output(ii)%item(i)%unit) output(ii)%Freq(n)
                do compo = 1, 3
                  do k1t = output(ii)%item(i)%ZItrancos, output(ii)%item(i)%ZEtrancos
                    do j1t = output(ii)%item(i)%YItrancos, output(ii)%item(i)%YEtrancos
                      IF (SGG%Observation(ii)%Transfer) then
             write (output(ii)%item(i)%unit) (output(ii)%item(i)%valor3DComplex(N, compo, i1t, j1t, k1t)/output(ii)%dftEntrada(n), &
                       &                                             i1t=output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos)
                      else !solo la transformada sin normalizar
                        write (output(ii)%item(i)%unit) (output(ii)%item(i)%valor3DComplex(N, compo, i1t, j1t, k1t), &
          &                                             i1t=output(ii)%item(i)%XItrancos, output(ii)%item(i)%XEtrancos)
                      end if
                    end do
                  end do
                end do
              end do
              close (output(ii)%item(i)%unit)
            case (icur, iCurX, iCurY, iCurZ)  !!!quitadp de aqui el mapvtk porque nunca puede estar en frecuencia!!!! 050216
              at = sgg%tiempo(FinalInstant)
              output(ii)%TimesWritten = output(ii)%NumFreqs  !util para leer el numero exacto de freq points
              INQUIRE (file=trim(adjustl(output(ii)%item(i)%path)), OPENED=ISyaopen)
              if (isyaopen) close (output(ii)%item(i)%unit, status='delete')
              !
              my_iostat = 0
9245          if (my_iostat /= 0) write (*, fmt='(a)', advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',9245,'.',layoutnumber,trim(adjustl(output(ii)%item(i)%path))
         open (output(ii)%item(i)%unit, file=trim(adjustl(output(ii)%item(i)%path)), form='unformatted', err=9245, iostat=my_iostat)
              write (output(ii)%item(i)%unit) output(ii)%item(i)%columnas
              do conta = 1, output(ii)%item(i)%columnas
                write (output(ii)%item(i)%unit) &
                  output(ii)%item(i)%Serialized%eI(conta), &
                  output(ii)%item(i)%Serialized%eJ(conta), &
                  output(ii)%item(i)%Serialized%eK(conta), &
                  output(ii)%item(i)%Serialized%currentType(conta), &
                  output(ii)%item(i)%Serialized%sggMtag(conta) !added to resuming file 121020
              end do
              write (output(ii)%item(i)%unit) at !deteccion errores dft si se resumea a partir de instantes posteriores al ultimo escrito
              DO N = 1, output(ii)%NumFreqs
                write (output(ii)%item(i)%unit) output(ii)%Freq(n)
                IF (SGG%Observation(ii)%Transfer) then
                  if (output(ii)%item(i)%columnas /= 0) then
                    do i1 = 1, output(ii)%item(i)%columnas
                      write (output(ii)%item(i)%unit) &
                        output(ii)%item(i)%Serialized%valorComplex_x(N, i1)/output(ii)%dftEntrada(n), &
                        output(ii)%item(i)%Serialized%valorComplex_y(N, i1)/output(ii)%dftEntrada(n), &
                        output(ii)%item(i)%Serialized%valorComplex_z(N, i1)/output(ii)%dftEntrada(n)
                    end do
                  end if
                ELSE
                  if (output(ii)%item(i)%columnas /= 0) then

                    do i1 = 1, output(ii)%item(i)%columnas
                      write (output(ii)%item(i)%unit) &
                        output(ii)%item(i)%Serialized%valorComplex_x(N, i1), &
                        output(ii)%item(i)%Serialized%valorComplex_y(N, i1), &
                        output(ii)%item(i)%Serialized%valorComplex_z(N, i1)
                    end do
                  end if

                END IF
              end do
              close (output(ii)%item(i)%unit)
            end select
            !
          end if !del time domain
        end do loop_obb
      end do

      nInit = FinalInstant + 1
      !voids valor
      do ii = 1, sgg%NumberRequest
        do i = 1, sgg%Observation(ii)%nP
          if (SGG%Observation(ii)%TimeDomain) then
            field = sgg%observation(ii)%P(i)%what
            select case (field)
              !estas sondas no se anulan tras escribir
              !case (iCur,iCurX,iCurY,iCurZ,mapvtk)
              !    output(ii)%item(i)%Serialized%valor = 0.0_RKIND
              !case (iMHC,iHxC,iHyC,iHzC,iMEC,iExC,iEyC,iEzC)
              !    output(ii)%item(i)%valor3D = 0.0_RKIND
            case (iQx, iQy, iQz)
              output(ii)%item(i)%valor(0:BuffObse) = 0.0_RKIND
            case (iJx, iJy, iJz)
              output(ii)%item(i)%valor(0:BuffObse) = 0.0_RKIND
              output(ii)%item(i)%valor2(0:BuffObse) = 0.0_RKIND
              output(ii)%item(i)%valor3(0:BuffObse) = 0.0_RKIND
              output(ii)%item(i)%valor4(0:BuffObse) = 0.0_RKIND
              output(ii)%item(i)%valor5(0:BuffObse) = 0.0_RKIND
            case (iBloqueMx, iBloqueMz, iBloqueMy, iBloqueJx, iBloqueJz, iBloqueJy, iEx, iEy, iEz, iHx, iHy, iHz)
              output(ii)%item(i)%valor(0:BuffObse) = 0.0_RKIND
            end select
          end if
        end do
      end do

      return
    end subroutine FlushObservationFiles

#ifdef CompileWithMTLN
    subroutine FlushMTLNObservationFiles(nEntradaRoot, mtlnProblem)
      character(len=*), intent(in)  ::  nEntradaRoot
      logical, intent(in) :: mtlnProblem
      type(mtln_solver_t), pointer :: mtln_solver
      integer :: i, j, k, n
      integer :: unit
      character(len=bufsize)  ::  temp
      character(len=bufsize)  ::  path
      character(len=:), allocatable :: buffer
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif

      mtln_solver => GetSolverPtr()
      unit = 2000
      do i = 1, size(mtln_solver%bundles)
        do j = 1, size(mtln_solver%bundles(i)%probes)
          path = trim(trim(nEntradaRoot)//"_"//trim(mtln_solver%bundles(i)%probes(j)%name)//".dat")
          open (unit=unit, file=trim(path))
          write (*, *) 'name: ', trim(mtln_solver%bundles(i)%probes(j)%name)
          buffer = "time"

          do k = 1, size(mtln_solver%bundles(i)%probes(j)%val, 2)
            write (temp, *) k
            buffer = buffer//" "//"conductor_"//trim(adjustl(temp))
          end do
          write (unit, *) trim(buffer)
          do k = 1, size(mtln_solver%bundles(i)%probes(j)%t)
            buffer = ""
            write (temp, *) mtln_solver%bundles(i)%probes(j)%t(k)
            buffer = buffer//trim(temp)
            do n = 1, size(mtln_solver%bundles(i)%probes(j)%val, 2)
              write (temp, *) mtln_solver%bundles(i)%probes(j)%val(k, n)
              buffer = buffer//" "//trim(temp)
            end do
            write (unit, *) trim(buffer)
          end do
          close (unit)
          unit = unit + 1
        end do
      end do

    end subroutine
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Free up memory
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine DestroyObservation(sgg)
      type(SGGFDTDINFO), intent(INOUT)         ::  sgg
      integer(kind=4)  ::  ii, i, field

#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif

#ifdef CompileWithMPI
      if (associated(valores)) deallocate (valores, newvalores)
#endif
      do ii = 1, sgg%NumberRequest
        if (SGG%Observation(ii)%Transfer) deallocate (output(ii)%dftEntrada)
        if (SGG%Observation(ii)%FreqDomain) deallocate (output(ii)%auxExp_E, output(ii)%auxExp_H, output(ii)%Freq)
        DO i = 1, sgg%Observation(ii)%nP
          field = sgg%observation(ii)%P(i)%what
          select case (field)
          case (iQx, iQy, iQz)
            deallocate (output(ii)%item(i)%valor)
          case (iJx, iJy, iJz)
            deallocate (output(ii)%item(i)%valor)
            deallocate (output(ii)%item(i)%valor2, output(ii)%item(i)%valor3, output(ii)%item(i)%valor4, output(ii)%item(i)%valor5)  !en caso de hilos se necesitan
          case (iBloqueJx, iBloqueJy, iBloqueMx, iBloqueMy)
            deallocate (output(ii)%item(i)%valor)
          case (lineIntegral)
            deallocate (output(ii)%item(i)%valor)
#ifdef CompileWithMPI
            if (output(ii)%item(i)%MPISubComm /= -1) then
              call MPI_Group_free(output(ii)%item(i)%MPIgroupindex, ierr)
            end if
#endif
          case (iMHC, iHxC, iHyC, iHzC, iMEC, iExC, iEyC, iEzC)
#ifdef CompileWithMPI
            if (output(ii)%item(i)%MPISubComm /= -1) then
              call MPI_Group_free(output(ii)%item(i)%MPIgroupindex, ierr)
            end if
#endif
            if (SGG%Observation(ii)%TimeDomain) deallocate (output(ii)%item(i)%valor3D)
            if (SGG%Observation(ii)%FreqDomain) deallocate (output(ii)%item(i)%valor3DComplex)
          case (iCur, iCurX, iCurY, iCurZ, mapvtk) !!!
#ifdef CompileWithMPI
            if (output(ii)%item(i)%MPISubComm /= -1) then
              call MPI_Group_free(output(ii)%item(i)%MPIgroupindex, ierr)
            end if
#endif
            if (SGG%Observation(ii)%TimeDomain) deallocate (output(ii)%item(i)%Serialized%valor, &
                                                            output(ii)%item(i)%Serialized%valor_x, &
                                                            output(ii)%item(i)%Serialized%valor_y, &
                                                            output(ii)%item(i)%Serialized%valor_z)
            if (SGG%Observation(ii)%FreqDomain) then
              deallocate (output(ii)%item(i)%Serialized%valorComplex_x)
              deallocate (output(ii)%item(i)%Serialized%valorComplex_y)
              deallocate (output(ii)%item(i)%Serialized%valorComplex_z)
            end if
            deallocate (output(ii)%item(i)%Serialized%eI)
            deallocate (output(ii)%item(i)%Serialized%eJ)
            deallocate (output(ii)%item(i)%Serialized%eK)
            deallocate (output(ii)%item(i)%Serialized%currentType)
            deallocate (output(ii)%item(i)%Serialized%sggMtag)

          case (iBloqueMz, iBloqueJz, iEx, iEy, iEz, iHx, iHy, iHz)
            deallocate (output(ii)%item(i)%valor)
          case (farfield)
            call DestroyFarField
#ifdef CompileWithMPI
            if (output(ii)%item(i)%MPISubComm /= -1) then
              call MPI_Group_free(output(ii)%item(i)%MPIgroupindex, ierr)
            end if
#endif
          end select
        end do
        if (associated(sgg%Observation(ii)%P)) deallocate (sgg%Observation(ii)%P)
        if (associated(output(ii)%item)) deallocate (output(ii)%item)

      end do

      if (associated(sgg%Observation)) deallocate (sgg%Observation)
      if (associated(output)) deallocate (output)

    end subroutine

   !!!!!!!!!!!!!!!!!!!!!

    function prefix(campo) result(ext)
      integer(kind=4)  ::  campo
      character(len=BUFSIZE)  ::  ext

      select case (campo)
      case (iEx)
        ext = 'Ex_'
      case (iEy)
        ext = 'Ey_'
      case (iEz)
        ext = 'Ez_'
      case (iVx)
        ext = 'Vx_'
      case (iVy)
        ext = 'Vy_'
      case (iVz)
        ext = 'Vz_'
      case (iHx)
        ext = 'Hx_'
      case (iHy)
        ext = 'Hy_'
      case (iHz)
        ext = 'Hz_'
      case (iBloqueJx)
        ext = 'Jx_'
      case (iBloqueJy)
        ext = 'Jy_'
      case (iBloqueJz)
        ext = 'Jz_'
      case (iBloqueMx)
        ext = 'Mx_'
      case (iBloqueMy)
        ext = 'My_'
      case (iBloqueMz)
        ext = 'Mz_'
      case (iJx)
        ext = 'Wx_'
      case (iJy)
        ext = 'Wy_'
      case (iJz)
        ext = 'Wz_'
      case (iQx)
        ext = 'Qx_'
      case (iQy)
        ext = 'Qy_'
      case (iQz)
        ext = 'Qz_'
      case (iExC)
        ext = 'ExC_'
      case (iEyC)
        ext = 'EyC_'
      case (iEzC)
        ext = 'EzC_'
      case (iHxC)
        ext = 'HxC_'
      case (iHyC)
        ext = 'HyC_'
      case (iHzC)
        ext = 'HzC_'
      case (iMEC)
        ext = 'ME_'
      case (iMHC)
        ext = 'MH_'
      case (iCur)
        ext = 'BC_'
      case (mapvtk)
        ext = 'MAP_'
      case (iCurX)
        ext = 'BCX_'
      case (iCurY)
        ext = 'BCY_'
      case (iCurZ)
        ext = 'BCZ_'
      case (farfield)
        ext = 'FF_'
      case (lineIntegral)
        ext = 'LI_'
      end select

      return

    end function prefix

    function suffix(campo, incid) result(ext)
      integer(kind=4)  ::  campo
      character(LEN=BUFSIZE)  ::  ext
      logical  ::  incid

      ext = ' '

      select case (campo)
      case (iEx, iEy, iEz, iHx, iHy, iHz)
        if (incid) ext = 'incid'
      case (iJx, iJy, iJz)
        ext = ' -E*dl Vplus Vminus Vplus-Vminus'
      end select

      return

    end function suffix

    function fieldo(field, dir) result(fieldo2)
      integer  ::  fieldo2, field
      character(len=1) :: dir
      fieldo2 = -1
      select case (field)
      case (iEx, iEy, iEz, iHx, iHy, iHz)
        fieldo2 = field
      case (iJx, iVx, iBloqueJx, iExC, iQx)
        fieldo2 = iEx
      case (iJy, iVy, iBloqueJy, iEyC, iQy)
        fieldo2 = iEy
      case (iJz, iVz, iBloqueJz, iEzC, iQz)
        fieldo2 = iEz
      case (iBloqueMx, iHxC)
        fieldo2 = iHx
      case (iBloqueMy, iHyC)
        fieldo2 = iHy
      case (iBloqueMz, iHzC)
        fieldo2 = iHz
      case (iMEC)
        select case (dir)
        CASE ('X', 'x')
          fieldo2 = iEx
        CASE ('Y', 'y')
          fieldo2 = iEY
        CASE ('Z', 'z')
          fieldo2 = iEz
        END SELECT
      case (iMHC)
        select case (dir)
        CASE ('X', 'x')
          fieldo2 = ihx
        CASE ('Y', 'y')
          fieldo2 = iHY
        CASE ('Z', 'z')
          fieldo2 = iHz
        END SELECT
      case (iCur, iCurX, icurY, icurZ, mapvtk)  !los pongo en efield para evitar problemas con el MPI
        select case (dir)
        CASE ('X', 'x')
          fieldo2 = iEx
        CASE ('Y', 'y')
          fieldo2 = iEY
        CASE ('Z', 'z')
          fieldo2 = iEz
        END SELECT
      end select
    end function

   !!!cuenta los bordes adyacentes
    subroutine contabordes(sgg, imed, imed1, imed2, imed3, imed4, EsBorde, SINPML_fullsize, campo, iii, jjj, kkk)
      type(SGGFDTDINFO), intent(IN)       ::  sgg
      type(limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      integer(Kind=4) imed, imed1, imed2, imed3, imed4, contaborde, campo, iii, jjj, kkk
      logical :: esborde
      !!!!
      esborde = .false.
      contaborde = 0
      !esta primera opcion solo considera bordes los externos
      IF (imed /= 1) THEN
        !    if     (sgg%med(imed )%is%SGBC) then
        !        if (sgg%med(imed1)%is%SGBC) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed2)%is%SGBC) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed3)%is%SGBC) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed4)%is%SGBC) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !   elseif  (sgg%med(imed )%is%Multiport) then
        !        if (sgg%med(imed1)%is%Multiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed2)%is%Multiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed3)%is%Multiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed4)%is%Multiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !    elseif (sgg%med(imed )%is%AnisMultiport) then
        !        if (sgg%med(imed1)%is%AnisMultiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed1)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed2)%is%AnisMultiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed2)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed3)%is%AnisMultiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed3)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !        if (sgg%med(imed4)%is%AnisMultiport) THEN
        !            if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))==trim(adjustl(sgg%Med(imed4)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
        !        endif
        !    else
        !        if (imed==imed1) contaborde=contaborde+1
        !        if (imed==imed2) contaborde=contaborde+1
        !        if (imed==imed3) contaborde=contaborde+1
        !        if (imed==imed4) contaborde=contaborde+1
        !    endif
        !    if (contaborde <=1) esborde=.true.
         !!!!alternativa
        if (sgg%med(imed)%is%SGBC) then
          if (sgg%med(imed1)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed1 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed2)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed2 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed3)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed3 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed4)%is%SGBC) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed4 /= 1) THEN
            contaborde = contaborde + 1
          end if
        elseif (sgg%med(imed)%is%Multiport) then
          if (sgg%med(imed1)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed1)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed1 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed2)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed2)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed2 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed3)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed3)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed3 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed4)%is%Multiport) THEN
               if (trim(adjustl(sgg%Med(imed )%Multiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed4)%Multiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed4 /= 1) THEN
            contaborde = contaborde + 1
          end if
        elseif (sgg%med(imed)%is%AnisMultiport) then
          if (sgg%med(imed1)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed1)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed1 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed2)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed2)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed2 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed3)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed3)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed3 /= 1) THEN
            contaborde = contaborde + 1
          end if
          if (sgg%med(imed4)%is%AnisMultiport) THEN
               if (trim(adjustl(sgg%Med(imed )%AnisMultiport(1)%MultiportFileZ11))/=trim(adjustl(sgg%Med(imed4)%AnisMultiport(1)%MultiportFileZ11)) ) contaborde=contaborde+1
          elseif (imed4 /= 1) THEN
            contaborde = contaborde + 1
          end if
        else
          if ((imed /= imed1) .and. (imed1 /= 1)) contaborde = contaborde + 1
          if ((imed /= imed2) .and. (imed2 /= 1)) contaborde = contaborde + 1
          if ((imed /= imed3) .and. (imed3 /= 1)) contaborde = contaborde + 1
          if ((imed /= imed4) .and. (imed4 /= 1)) contaborde = contaborde + 1
        end if
        if ((imed1 == 1) .and. (imed2 == 1) .and. (imed3 == 1) .and. (imed4 /= 1)) esborde = .true. !un borde con vacion
        if ((imed2 == 1) .and. (imed3 == 1) .and. (imed4 == 1) .and. (imed1 /= 1)) esborde = .true. !un borde con vacion
        if ((imed3 == 1) .and. (imed4 == 1) .and. (imed1 == 1) .and. (imed2 /= 1)) esborde = .true. !un borde con vacion
        if ((imed4 == 1) .and. (imed1 == 1) .and. (imed2 == 1) .and. (imed3 /= 1)) esborde = .true. !un borde con vacion
        if (contaborde > 0) esborde = .true.
        if ((imed1 == 1) .and. (imed2 == 1) .and. (imed3 == 1) .and. (imed4 == 1)) esborde = .true. !un segmento aislado
         !!!!!! Fin de la alternativa para que considere bordes tambien ejes internos
         !!!!!!!!!!!!!!
         !!!!!!!!!!!!!
        if ((iii > SINPML_fullsize(campo)%XE) .or. (jjj > SINPML_fullsize(campo)%YE) .or. (kkk > SINPML_fullsize(campo)%ZE)) then
          esborde = .false.
        end if
        if ((iii < SINPML_fullsize(campo)%XI) .or. (jjj < SINPML_fullsize(campo)%Yi) .or. (kkk < SINPML_fullsize(campo)%Zi)) then
          esborde = .false.
        end if
      else
        esborde = .false.
      END IF !DEL IMED1
      return
    end subroutine contabordes

    subroutine nodalvtk(sgg, sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz, sggMtag, tag_numbers, &
                        init, geom, asigna, electric, magnetic, conta, i, ii, output, Ntimeforvolumic)
      type(SGGFDTDINFO), intent(IN)         ::  sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)

      type(output_t), pointer, dimension(:)  ::  output
      integer(kind=4), intent(IN) :: i, ii, Ntimeforvolumic

      logical geom, asigNa, init, electric, magnetic
      integer(kind=4) conta, sweep, ni, nj, nk, i_m, j_m, k_m, IMED
      type(taglist_t) :: tag_numbers
      integer(KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in)   :: &
        sggMiEx(sgg%alloc(iEx)%XI:sgg%alloc(iEx)%XE, sgg%alloc(iEx)%YI:sgg%alloc(iEx)%YE, sgg%alloc(iEx)%ZI:sgg%alloc(iEx)%ZE), &
        sggMiEy(sgg%alloc(iEy)%XI:sgg%alloc(iEy)%XE, sgg%alloc(iEy)%YI:sgg%alloc(iEy)%YE, sgg%alloc(iEy)%ZI:sgg%alloc(iEy)%ZE), &
        sggMiEz(sgg%alloc(iEz)%XI:sgg%alloc(iEz)%XE, sgg%alloc(iEz)%YI:sgg%alloc(iEz)%YE, sgg%alloc(iEz)%ZI:sgg%alloc(iEz)%ZE), &
        sggMiHx(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE, sgg%alloc(iHx)%YI:sgg%alloc(iHx)%YE, sgg%alloc(iHx)%ZI:sgg%alloc(iHx)%ZE), &
        sggMiHy(sgg%alloc(iHy)%XI:sgg%alloc(iHy)%XE, sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE, sgg%alloc(iHy)%ZI:sgg%alloc(iHy)%ZE), &
        sggMiHz(sgg%alloc(iHz)%XI:sgg%alloc(iHz)%XE, sgg%alloc(iHz)%YI:sgg%alloc(iHz)%YE, sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)

      !to fetch info of nodal sources for the vtkmap
      type(nodsou), pointer, save :: rNodal_Ex, rNodal_Ey, rNodal_Ez
      type(nodsou), pointer, save :: rNodal_Hx, rNodal_Hy, rNodal_Hz
      !!!!!!!!!!!!

      if (init) call getnodal(rNodal_Ex, rNodal_Ey, rNodal_Ez, rNodal_Hx, rNodal_Hy, rNodal_Hz)
      IF (ELECTRIC) THEN
        do sweep = 1, rNodal_Ex%numHard
          do nk = rNodal_Ex%nodHard(sweep)%punto%zi, rNodal_Ex%nodHard(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Ex%nodHard(sweep)%punto%yi, rNodal_Ex%nodHard(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Ex%nodHard(sweep)%punto%xi, rNodal_Ex%nodHard(sweep)%punto%xe
                i_m = ni
                imed = sggMiEx(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PEC) then
                  conta = conta + 1
                  if (geom) then
                    ! print *,'-------antes de GEOM ',II,I, conta
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iJx
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%x(ni, nj, nk))
                    ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                  end if
                  if (asigna) then
                    ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                    output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 8.5
                    ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                  end if
                end if
              end do
            End do
          end do
        end do
        !
        do sweep = 1, rNodal_Ex%numSoft
          do nk = rNodal_Ex%nodSoft(sweep)%punto%zi, rNodal_Ex%nodSoft(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Ex%nodSoft(sweep)%punto%yi, rNodal_Ex%nodSoft(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Ex%nodSoft(sweep)%punto%xi, rNodal_Ex%nodSoft(sweep)%punto%xe
                i_m = ni
                imed = sggMiEx(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PEC) then
                  conta = conta + 1
                  if (geom) then
                    ! print *,'-------antes de GEOM ',II,I, conta
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iJx
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%x(ni, nj, nk))
                    ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                  end if
                  if (asigna) then
                    ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                    output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 8.5
                    ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                  end if
                end if
              end do
            End do
          end do
        end do
        !
        !
        do sweep = 1, rNodal_Ey%numHard
          do nk = rNodal_Ey%nodHard(sweep)%punto%zi, rNodal_Ey%nodHard(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Ey%nodHard(sweep)%punto%yi, rNodal_Ey%nodHard(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Ey%nodHard(sweep)%punto%xi, rNodal_Ey%nodHard(sweep)%punto%xe
                i_m = ni
                imed = sggMiEy(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PEC) then
                  conta = conta + 1
                  if (geom) then
                    ! print *,'-------antes de GEOM ',II,I, conta
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iJy
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%y(ni, nj, nk))
                    ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                  end if
                  if (asigna) then
                    ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                    output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 8.5
                    ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                  end if
                end if
              end do
            End do
          end do
        end do
        !
        do sweep = 1, rNodal_Ey%numSoft
          do nk = rNodal_Ey%nodSoft(sweep)%punto%zi, rNodal_Ey%nodSoft(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Ey%nodSoft(sweep)%punto%yi, rNodal_Ey%nodSoft(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Ey%nodSoft(sweep)%punto%xi, rNodal_Ey%nodSoft(sweep)%punto%xe
                i_m = ni
                imed = sggMiEy(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PEC) then
                  conta = conta + 1
                  if (geom) then
                    ! print *,'-------antes de GEOM ',II,I, conta
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iJy
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%y(ni, nj, nk))
                    ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                  end if
                  if (asigna) then
                    ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                    output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 8.5
                    ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                  end if
                end if
              end do
            End do
          end do
        end do

        do sweep = 1, rNodal_Ez%numHard
          do nk = rNodal_Ez%nodHard(sweep)%punto%zi, rNodal_Ez%nodHard(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Ez%nodHard(sweep)%punto%yi, rNodal_Ez%nodHard(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Ez%nodHard(sweep)%punto%xi, rNodal_Ez%nodHard(sweep)%punto%xe
                i_m = ni
                imed = sggMiEz(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PEC) then
                  conta = conta + 1
                  if (geom) then
                    ! print *,'-------antes de GEOM ',II,I, conta
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iJz
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%z(ni, nj, nk))
                    ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                  end if
                  if (asigna) then
                    ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                    output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 8.5
                    ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                  end if
                end if
              end do
            End do
          end do
        end do
        !
        do sweep = 1, rNodal_Ez%numSoft
          do nk = rNodal_Ez%nodSoft(sweep)%punto%zi, rNodal_Ez%nodSoft(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Ez%nodSoft(sweep)%punto%yi, rNodal_Ez%nodSoft(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Ez%nodSoft(sweep)%punto%xi, rNodal_Ez%nodSoft(sweep)%punto%xe
                i_m = ni
                imed = sggMiEz(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PEC) then
                  conta = conta + 1
                  if (geom) then
                    ! print *,'-------antes de GEOM ',II,I, conta
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iJz
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%z(ni, nj, nk))
                    ! print *,'-------tras  GEOM',output(ii)%item(i)%Serialized%eI(conta),output(ii)%item(i)%Serialized%currentType(conta)
                  end if
                  if (asigna) then
                    ! print *,'-------antes de asigna ', Ntimeforvolumic,conta
                    output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 8.5
                    ! print *,'-------tras  de asigna ',output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta)
                  end if
                end if
              end do
            End do
          end do
        end do
      END IF !DEL ELECTRIC

      IF (MAGNETIC) THEN

        do sweep = 1, rNodal_Hx%numHard
          do nk = rNodal_Hx%nodHard(sweep)%punto%zi, rNodal_Hx%nodHard(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Hx%nodHard(sweep)%punto%yi, rNodal_Hx%nodHard(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Hx%nodHard(sweep)%punto%xi, rNodal_Hx%nodHard(sweep)%punto%xe
                i_m = ni
                imed = sggMiHx(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PMC) then
                  conta = conta + 1
                  if (geom) then
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iBloqueJx
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%face%x(ni, nj, nk))
                  end if
                  if (asigna) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 9.0
                end if
              end do
            End do
          end do
        end do
        !
        do sweep = 1, rNodal_Hx%numSoft
          do nk = rNodal_Hx%nodSoft(sweep)%punto%zi, rNodal_Hx%nodSoft(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Hx%nodSoft(sweep)%punto%yi, rNodal_Hx%nodSoft(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Hx%nodSoft(sweep)%punto%xi, rNodal_Hx%nodSoft(sweep)%punto%xe
                i_m = ni
                imed = sggMiHx(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PMC) then
                  conta = conta + 1
                  if (geom) then
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iBloqueJx
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%face%x(ni, nj, nk))
                  end if
                  if (asigna) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 9.0
                end if
              end do
            End do
          end do
        end do
        !
        !
        do sweep = 1, rNodal_Hy%numHard
          do nk = rNodal_Hy%nodHard(sweep)%punto%zi, rNodal_Hy%nodHard(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Hy%nodHard(sweep)%punto%yi, rNodal_Hy%nodHard(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Hy%nodHard(sweep)%punto%xi, rNodal_Hy%nodHard(sweep)%punto%xe
                i_m = ni
                imed = sggMiHx(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PMC) then
                  conta = conta + 1
                  if (geom) then
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iBloqueJy
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%face%y(ni, nj, nk))
                  end if
                  if (asigna) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 9.0
                end if
              end do
            End do
          end do
        end do
        !
        do sweep = 1, rNodal_Hy%numSoft
          do nk = rNodal_Hy%nodSoft(sweep)%punto%zi, rNodal_Hy%nodSoft(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Hy%nodSoft(sweep)%punto%yi, rNodal_Hy%nodSoft(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Hy%nodSoft(sweep)%punto%xi, rNodal_Hy%nodSoft(sweep)%punto%xe
                i_m = ni
                imed = sggMiHy(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PMC) then
                  conta = conta + 1
                  if (geom) then
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iBloqueJy
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%face%y(ni, nj, nk))
                  end if
                  if (asigna) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 9.0
                end if
              end do
            End do
          end do
        end do

        do sweep = 1, rNodal_Hz%numHard
          do nk = rNodal_Hz%nodHard(sweep)%punto%zi, rNodal_Hz%nodHard(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Hz%nodHard(sweep)%punto%yi, rNodal_Hz%nodHard(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Hz%nodHard(sweep)%punto%xi, rNodal_Hz%nodHard(sweep)%punto%xe
                i_m = ni
                imed = sggMiHx(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PMC) then
                  conta = conta + 1
                  if (geom) then
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iBloqueJz
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%face%z(ni, nj, nk))
                  end if
                  if (asigna) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 9.0
                end if
              end do
            End do
          end do
        end do
        !
        do sweep = 1, rNodal_Hz%numSoft
          do nk = rNodal_Hz%nodSoft(sweep)%punto%zi, rNodal_Hz%nodSoft(sweep)%punto%ze
            k_m = nk
            do nj = rNodal_Hz%nodSoft(sweep)%punto%yi, rNodal_Hz%nodSoft(sweep)%punto%ye
              j_m = nj
              do ni = rNodal_Hz%nodSoft(sweep)%punto%xi, rNodal_Hz%nodSoft(sweep)%punto%xe
                i_m = ni
                imed = sggMiHz(i_m, j_m, k_m)
                if (.not. sgg%Med(imed)%Is%PMC) then
                  conta = conta + 1
                  if (geom) then
                    output(ii)%item(i)%Serialized%eI(conta) = ni
                    output(ii)%item(i)%Serialized%eJ(conta) = nj
                    output(ii)%item(i)%Serialized%eK(conta) = nk
                    output(ii)%item(i)%Serialized%currentType(conta) = iBloqueJz
                    output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%face%z(ni, nj, nk))
                  end if
                  if (asigna) output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 9.0
                end if
              end do
            End do
          end do
        end do
      END IF

      !!!!!!

      !print *,'----tras nodalvtk init,geom,asigna,electric,magnetic,conta,i,ii ',init,geom,asigna,electric,magnetic,conta,i,ii
      return
    end subroutine

    subroutine wirebundlesvtk(sgg, init, geom, asigna, conta, i, ii, output, Ntimeforvolumic, wiresflavor, sggMtag, tag_numbers)

      type(SGGFDTDINFO), intent(IN)   :: sgg
      INTEGER (KIND=IKINDMTAG), intent(in) :: sggMtag  (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, sgg%Alloc(iHy)%YI:sgg%Alloc(iHy)%YE, sgg%Alloc(iHz)%ZI:sgg%Alloc(iHz)%ZE)
      type(taglist_t) :: tag_numbers
      type(output_t), pointer, dimension(:)  ::  output
      integer(kind=4), intent(IN) :: i, ii, Ntimeforvolumic
      logical, intent(IN) :: geom, asigNa, init
      integer(kind=4) conta, ni, nj, nk, n
      character(len=*), INTENT(in) :: wiresflavor
      integer(kind=4), SAVE :: MINIMED

      type(Thinwires_t), pointer, save  ::  Hwireslocal
      !
#ifdef CompileWithBerengerWires
      type(TWires), pointer, save  ::  Hwireslocal_Berenger
#endif
#ifdef CompileWithSlantedWires
      !
      type(WiresData), pointer, save  ::  Hwireslocal_Slanted
#endif

      !print *,'----antes wires init,geom,asigna,conta,i,ii',init,geom,asigna,conta,i,ii
      if (init) then
        if ((trim(adjustl(wiresflavor)) == 'holland') .or. &
            (trim(adjustl(wiresflavor)) == 'transition')) then
          Hwireslocal => GetHwires()
        end if
#ifdef CompileWithBerengerWires
        if (trim(adjustl(wiresflavor)) == 'berenger') then
          Hwireslocal_Berenger => GetHwires_Berenger()
        end if
#endif
#ifdef CompileWithSlantedWires
        if ((trim(adjustl(wiresflavor)) == 'slanted') .or. (trim(adjustl(wiresflavor)) == 'semistructured')) then
          Hwireslocal_Slanted => GetHwires_Slanted()
        end if
#endif
      end if
#ifdef CompileWithBerengerWires
      if (trim(adjustl(wiresflavor)) == 'berenger') then

        !parsea los hilos
        if (geom) then
          MINIMED = 2**12
          do n = 1, Hwireslocal_Berenger%NumSegments
            conta = conta + 1
            minimed = MIN(MINIMED, Hwireslocal_Berenger%Segments(n)%imeD)
            ni = Hwireslocal_Berenger%Segments(n)%ii
            nj = Hwireslocal_Berenger%Segments(n)%ji
            nk = Hwireslocal_Berenger%Segments(n)%ki

            output(ii)%item(i)%Serialized%eI(conta) = ni
            output(ii)%item(i)%Serialized%eJ(conta) = nj
            output(ii)%item(i)%Serialized%eK(conta) = nk

            select case (Hwireslocal_Berenger%Segments(n)%orient)
            case (iEx)
              output(ii)%item(i)%Serialized%currentType(conta) = iJx
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%x(ni, nj, nk))
            case (iEy)
              output(ii)%item(i)%Serialized%currentType(conta) = iJy
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%y(ni, nj, nk))
            case (iEz)
              output(ii)%item(i)%Serialized%currentType(conta) = iJz
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%z(ni, nj, nk))
            end select
          end do

        elseif (asigna) then
          do n = 1, Hwireslocal_Berenger%NumSegments
            conta = conta + 1
            if ((Hwireslocal_Berenger%Segments(n)%Is_LeftEnd) .or. &
                (Hwireslocal_Berenger%Segments(n)%Is_RightEnd)) then
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 10
            else
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 20 + Hwireslocal_Berenger%Segments(n)%imed - MINIMED
            end if
          end do
        else
          do n = 1, Hwireslocal_Berenger%NumSegments
            conta = conta + 1
          end do
        end if
      end if
#endif
      if ((trim(adjustl(wiresflavor)) == 'holland') .or. &
          (trim(adjustl(wiresflavor)) == 'transition')) then
        if (geom) then
          MINIMED = 2**30
          do n = 1, Hwireslocal%NumCurrentSegments
            conta = conta + 1
            minimed = MIN(MINIMED, Hwireslocal%CurrentSegment(n)%indexmed)
            ni = Hwireslocal%CurrentSegment(n)%i
            nj = Hwireslocal%CurrentSegment(n)%j
            nk = Hwireslocal%CurrentSegment(n)%k

            output(ii)%item(i)%Serialized%eI(conta) = ni
            output(ii)%item(i)%Serialized%eJ(conta) = nj
            output(ii)%item(i)%Serialized%eK(conta) = nk

            select case (Hwireslocal%CurrentSegment(n)%tipofield)
            case (iEx)
              output(ii)%item(i)%Serialized%currentType(conta) = iJx
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%x(ni, nj, nk))
            case (iEy)
              output(ii)%item(i)%Serialized%currentType(conta) = iJy
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%y(ni, nj, nk))
            case (iEz)
              output(ii)%item(i)%Serialized%currentType(conta) = iJz
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%z(ni, nj, nk))
            end select
          end do

        elseif (asigna) then
          do n = 1, Hwireslocal%NumCurrentSegments
            conta = conta + 1
            if ((Hwireslocal%CurrentSegment(n)%Is_LeftEnd) .or. &
                (Hwireslocal%CurrentSegment(n)%Is_RightEnd)) then
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 10
            elseif (Hwireslocal%CurrentSegment(n)%IsEnd_norLeft_norRight) then
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 11
            else
              !  output( ii)%item( i)%Serialized%valor(Ntimeforvolumic,conta) = 20 + Hwireslocal%CurrentSegment(n)%indexmed-MINIMED
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 20 + Hwireslocal%CurrentSegment(n)%NumParallel
            end if
          end do
        else
          do n = 1, Hwireslocal%NumCurrentSegments
            conta = conta + 1
          end do
        end if
      end if
#ifdef CompileWithSlantedWires
      if ((trim(adjustl(wiresflavor)) == 'slanted') .or. (trim(adjustl(wiresflavor)) == 'semistructured')) then
        !parsea los hilos
        if (geom) then
          MINIMED = 2**12
          do n = 1, Hwireslocal_Slanted%NumSegmentsStr
            conta = conta + 1
            minimed = MIN(MINIMED, Hwireslocal_Slanted%SegmentsStr(n)%imeD)
            ni = Hwireslocal_Slanted%SegmentsStr(n)%cell(iX)
            nj = Hwireslocal_Slanted%SegmentsStr(n)%cell(iY)
            nk = Hwireslocal_Slanted%SegmentsStr(n)%cell(iZ)

            output(ii)%item(i)%Serialized%eI(conta) = ni
            output(ii)%item(i)%Serialized%eJ(conta) = nj
            output(ii)%item(i)%Serialized%eK(conta) = nk

            select case (Hwireslocal_Slanted%SegmentsStr(n)%orient)
            case (iEx)
              output(ii)%item(i)%Serialized%currentType(conta) = iJx
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%x(ni, nj, nk))
            case (iEy)
              output(ii)%item(i)%Serialized%currentType(conta) = iJy
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%y(ni, nj, nk))
            case (iEz)
              output(ii)%item(i)%Serialized%currentType(conta) = iJz
              output(ii)%item(i)%Serialized%sggMtag(conta) = iabs(tag_numbers%edge%z(ni, nj, nk))
            end select
          end do

        elseif (asigna) then
          do n = 1, Hwireslocal_Slanted%NumSegmentsStr
            conta = conta + 1
            if ((Hwireslocal_Slanted%SegmentsStr(n)%IsExt(iBeg)) .or. &
                (Hwireslocal_Slanted%SegmentsStr(n)%IsExt(iEnd))) then
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 10
            else
              output(ii)%item(i)%Serialized%valor(Ntimeforvolumic, conta) = 20 + Hwireslocal_Slanted%SegmentsStr(n)%imed - MINIMED
            end if
          end do
        else
          do n = 1, Hwireslocal_Slanted%NumSegmentsStr
            conta = conta + 1
          end do
        end if

      end if
#endif

      !print *,'----tras wires init,geom,asigna,conta,i,ii',init,geom,asigna,conta,i,ii

      return
    end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !Function to publish the private output data (used in postprocess)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function GetOutput() result(r)
      type(output_t), pointer, dimension(:)  ::  r

      r => output
      return
    end function

    ! ========================================================================
    ! PURPOSE:
    ! Performs Discrete Time Fourier Transform in a signal given in sig
    ! sampled at times st.
    ! The frequency values are stored in fqVal for frequencies given in fq.
    ! This subroutine is efficient when fqSize << sigSize.
    ! ========================================================================
    subroutine dtft(fqVal, fq, fqSize, st, sig, sigSize)
      implicit none
      ! ----------- Variables --------------------------------------------------
      ! In & out vars.
      ! Size of the input frequencies and signal vectors.
      integer, intent(in) :: fqSize, sigSize
      ! Fourier transform value.
      complex(kind=CKIND), intent(out), dimension(fqSize) :: fqVal
      ! Vector of frequencies to compute the values.
      real(kind=RKIND), intent(in), dimension(fqSize) ::  fq
      ! Input signal.
      real(kind=RKIND), intent(in), dimension(sigSize) :: sig
      ! Input signal sampling time.
      real(kind=RKIND_tiempo), intent(in), dimension(sigSize) :: st
      ! Aux. vars.
      integer :: i, j; 
      ! ----------- Body -------------------------------------------------------
      ! Checks that frequencies are below the sf.
      fqval = 0.0_RKIND

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j)
#endif
      do i = 1, fqSize
        do j = 1, sigSize - 1 !algun delta promedio habra que tomar permit scaling 211118
          fqVal(i) = fqVal(i) + abs(st(j + 1) - st(j))/2.0_RKIND*sig(j)*exp(mcPI2*fq(i)*st(j)); !nosisosampleada 200120 !ojo valor absoluto delta
!!            fqVal(i) = fqVal(i) + abs(st(2)-st(1))/2.0_RKIND * sig(j) * exp(mcPI2 * fq(i) * (j-1)*(st(2)-st(1))); !isosampleada 200120 !ojo valor absoluto delta
!!            if ((i==1).and.(j==sigsize-1)) write (634,*) j,abs(st(2)-st(1)),abs(st(j+1)-st(j))
        end do
      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

#ifdef CompileWithOpenMP
!$OMP PARALLEL DO DEFAULT(SHARED) private (i,j)
#endif
      do i = 1, fqSize
        j = sigSize  !algun delta promedio habra que tomar permit scaling 211118
        !no debia mucho tener impacto por ser el ultimo timpicamente pequenio, pero....
        fqVal(i) = fqVal(i) + abs(st(j - 1) - st(j))/2.0_RKIND*sig(j)*exp(mcPI2*fq(i)*st(j)); !nosisosampleada 200120 !ojo valor absoluto delta  !ojo valor absoluto delta pq el ultimo cambiaba elsigno 14111'
!!         fqVal(i) = fqVal(i) + abs(st(2)-st(1))/2.0_RKIND * sig(j) * exp(mcPI2 * fq(i) * (j-1)*(st(2)-st(1))); !isosampleada 200120 !ojo valor absoluto delta

      end do
#ifdef CompileWithOpenMP
!$OMP END PARALLEL DO
#endif

      fqVal = fqVal*2.0_RKIND !BUG HIRAI ENERGIA DOBLE PARSEVAL  mail 24/07/19

    end subroutine

    real(kind=RKIND) function interpolate_field_atwhere(sgg, Ex, Ey, Ez, Hx, Hy, Hz, i, j, k, field, atwhere) result(interp)

      type(SGGFDTDINFO), intent(IN) :: sgg
      REAL(KIND=RKIND), intent(in), target :: &
        Ex(sgg%alloc(iEx)%XI:sgg%alloc(iEx)%XE, sgg%alloc(iEx)%YI:sgg%alloc(iEx)%YE, sgg%alloc(iEx)%ZI:sgg%alloc(iEx)%ZE), &
        Ey(sgg%alloc(iEy)%XI:sgg%alloc(iEy)%XE, sgg%alloc(iEy)%YI:sgg%alloc(iEy)%YE, sgg%alloc(iEy)%ZI:sgg%alloc(iEy)%ZE), &
        Ez(sgg%alloc(iEz)%XI:sgg%alloc(iEz)%XE, sgg%alloc(iEz)%YI:sgg%alloc(iEz)%YE, sgg%alloc(iEz)%ZI:sgg%alloc(iEz)%ZE), &
        Hx(sgg%alloc(iHx)%XI:sgg%alloc(iHx)%XE, sgg%alloc(iHx)%YI:sgg%alloc(iHx)%YE, sgg%alloc(iHx)%ZI:sgg%alloc(iHx)%ZE), &
        Hy(sgg%alloc(iHy)%XI:sgg%alloc(iHy)%XE, sgg%alloc(iHy)%YI:sgg%alloc(iHy)%YE, sgg%alloc(iHy)%ZI:sgg%alloc(iHy)%ZE), &
        Hz(sgg%alloc(iHz)%XI:sgg%alloc(iHz)%XE, sgg%alloc(iHz)%YI:sgg%alloc(iHz)%YE, sgg%alloc(iHz)%ZI:sgg%alloc(iHz)%ZE)

      integer, intent(in) :: i, j, k
      integer, intent(in) :: field, atwhere
   !! real (kind=RKIND) :: interp

      ! Index variables for each field
      integer :: im1_Ex, ip1_Ex, jm1_Ex, jp1_Ex, km1_Ex, kp1_Ex
      integer :: im1_Ey, ip1_Ey, jm1_Ey, jp1_Ey, km1_Ey, kp1_Ey
      integer :: im1_Ez, ip1_Ez, jm1_Ez, jp1_Ez, km1_Ez, kp1_Ez
      integer :: im1_Hx, ip1_Hx, jm1_Hx, jp1_Hx, km1_Hx, kp1_Hx
      integer :: im1_Hy, ip1_Hy, jm1_Hy, jp1_Hy, km1_Hy, kp1_Hy
      integer :: im1_Hz, ip1_Hz, jm1_Hz, jp1_Hz, km1_Hz, kp1_Hz

      ! Compute indices for Ex
      im1_Ex = Max(i - 1, sgg%alloc(iEx)%XI)
      ip1_Ex = Min(i + 1, sgg%alloc(iEx)%XE)
      jm1_Ex = Max(j - 1, sgg%alloc(iEx)%YI)
      jp1_Ex = Min(j + 1, sgg%alloc(iEx)%YE)
      km1_Ex = Max(k - 1, sgg%alloc(iEx)%ZI)
      kp1_Ex = Min(k + 1, sgg%alloc(iEx)%ZE)

      ! Compute indices for Ey
      im1_Ey = Max(i - 1, sgg%alloc(iEy)%XI)
      ip1_Ey = Min(i + 1, sgg%alloc(iEy)%XE)
      jm1_Ey = Max(j - 1, sgg%alloc(iEy)%YI)
      jp1_Ey = Min(j + 1, sgg%alloc(iEy)%YE)
      km1_Ey = Max(k - 1, sgg%alloc(iEy)%ZI)
      kp1_Ey = Min(k + 1, sgg%alloc(iEy)%ZE)

      ! Compute indices for Ez
      im1_Ez = Max(i - 1, sgg%alloc(iEz)%XI)
      ip1_Ez = Min(i + 1, sgg%alloc(iEz)%XE)
      jm1_Ez = Max(j - 1, sgg%alloc(iEz)%YI)
      jp1_Ez = Min(j + 1, sgg%alloc(iEz)%YE)
      km1_Ez = Max(k - 1, sgg%alloc(iEz)%ZI)
      kp1_Ez = Min(k + 1, sgg%alloc(iEz)%ZE)

      ! Compute indices for Hx
      im1_Hx = Max(i - 1, sgg%alloc(iHx)%XI)
      ip1_Hx = Min(i + 1, sgg%alloc(iHx)%XE)
      jm1_Hx = Max(j - 1, sgg%alloc(iHx)%YI)
      jp1_Hx = Min(j + 1, sgg%alloc(iHx)%YE)
      km1_Hx = Max(k - 1, sgg%alloc(iHx)%ZI)
      kp1_Hx = Min(k + 1, sgg%alloc(iHx)%ZE)

      ! Compute indices for Hy
      im1_Hy = Max(i - 1, sgg%alloc(iHy)%XI)
      ip1_Hy = Min(i + 1, sgg%alloc(iHy)%XE)
      jm1_Hy = Max(j - 1, sgg%alloc(iHy)%YI)
      jp1_Hy = Min(j + 1, sgg%alloc(iHy)%YE)
      km1_Hy = Max(k - 1, sgg%alloc(iHy)%ZI)
      kp1_Hy = Min(k + 1, sgg%alloc(iHy)%ZE)

      ! Compute indices for Hz
      im1_Hz = Max(i - 1, sgg%alloc(iHz)%XI)
      ip1_Hz = Min(i + 1, sgg%alloc(iHz)%XE)
      jm1_Hz = Max(j - 1, sgg%alloc(iHz)%YI)
      jp1_Hz = Min(j + 1, sgg%alloc(iHz)%YE)
      km1_Hz = Max(k - 1, sgg%alloc(iHz)%ZI)
      kp1_Hz = Min(k + 1, sgg%alloc(iHz)%ZE)

      ! Initialize output
      interp = 0.0

      ! Electric field interpolation at various positions
      if (atwhere == iEx) then
        if (field == iEx) then
          interp = Ex(i, j, k)
        elseif (field == iEy) then
          interp = (Ey(i, j, k) + Ey(i, jm1_Ey, k) + Ey(ip1_Ey, j, k) + Ey(ip1_Ey, jm1_Ey, k))/4.0
        elseif (field == iEz) then
          interp = (Ez(i, j, k) + Ez(i, j, km1_Ez) + Ez(ip1_Ez, j, k) + Ez(ip1_Ez, j, km1_Ez))/4.0
        elseif (field == iHx) then
          interp = (Hx(i, j, k) + Hx(i, j, km1_Hx) + Hx(i, jm1_Hx, k) + Hx(i, jm1_Hx, km1_Hx) + &
                    Hx(ip1_Hx, j, k) + Hx(ip1_Hx, j, km1_Hx) + Hx(ip1_Hx, jm1_Hx, k) + Hx(ip1_Hx, jm1_Hx, km1_Hx))/8.0
        elseif (field == iHy) then
          interp = (Hy(i, j, k) + Hy(i, j, km1_Hy))/2.0
        elseif (field == iHz) then
          interp = (Hz(i, j, k) + Hz(i, jm1_Hz, k))/2.0
        end if
      elseif (atwhere == iEy) then
        if (field == iEx) then
          interp = (Ex(i, j, k) + Ex(im1_Ex, j, k) + Ex(i, jp1_Ex, k) + Ex(im1_Ex, jp1_Ex, k))/4.0
        elseif (field == iEy) then
          interp = Ey(i, j, k)
        elseif (field == iEz) then
          interp = (Ez(i, j, k) + Ez(i, j, km1_Ez) + Ez(i, jp1_Ez, k) + Ez(i, jp1_Ez, km1_Ez))/4.0
        elseif (field == iHx) then
          interp = (Hx(i, j, k) + Hx(i, j, km1_Hx))/2.0
        elseif (field == iHy) then
          interp = (Hy(i, j, k) + Hy(im1_Hy, j, k) + Hy(i, j, km1_Hy) + Hy(im1_Hy, j, km1_Hy) + &
                    Hy(i, jp1_Hy, k) + Hy(im1_Hy, jp1_Hy, k) + Hy(i, jp1_Hy, km1_Hy) + Hy(im1_Hy, jp1_Hy, km1_Hy))/8.0
        elseif (field == iHz) then
          interp = (Hz(i, j, k) + Hz(im1_Hz, j, k))/2.0
        end if
      elseif (atwhere == iEz) then
        if (field == iEx) then
          interp = (Ex(i, j, k) + Ex(im1_Ex, j, k) + Ex(i, j, kp1_Ex) + Ex(im1_Ex, j, kp1_Ex))/4.0
        elseif (field == iEy) then
          interp = (Ey(i, j, k) + Ey(i, jm1_Ey, k) + Ey(i, j, kp1_Ey) + Ey(i, jm1_Ey, kp1_Ey))/4.0
        elseif (field == iEz) then
          interp = Ez(i, j, k)
        elseif (field == iHx) then
          interp = (Hx(i, j, k) + Hx(i, jm1_Hx, k))/2.0
        elseif (field == iHy) then
          interp = (Hy(i, j, k) + Hy(im1_Hy, j, k))/2.0
        elseif (field == iHz) then
          interp = (Hz(i, j, k) + Hz(i, jm1_Hz, k) + Hz(im1_Hz, j, k) + Hz(im1_Hz, jm1_Hz, k) + &
                    Hz(i, j, kp1_Hz) + Hz(i, jm1_Hz, kp1_Hz) + Hz(im1_Hz, j, kp1_Hz) + Hz(im1_Hz, jm1_Hz, kp1_Hz))/8.0
        end if
      end if

      ! Magnetic field interpolation at various positions
      if (atwhere == iHx) then
        if (field == iEx) then
          interp = (Ex(i, j, k) + Ex(im1_Ex, j, k) + Ex(i, jp1_Ex, k) + Ex(im1_Ex, jp1_Ex, k) + &
                    Ex(i, j, kp1_Ex) + Ex(im1_Ex, j, kp1_Ex) + Ex(i, jp1_Ex, kp1_Ex) + Ex(im1_Ex, jp1_Ex, kp1_Ex))/8.0
        elseif (field == iEy) then
          interp = (Ey(i, j, k) + Ey(i, j, kp1_Ey))/2.0
        elseif (field == iEz) then
          interp = (Ez(i, j, k) + Ez(i, jp1_Ez, k))/2.0
        elseif (field == iHx) then
          interp = Hx(i, j, k)
        elseif (field == iHy) then
          interp = (Hy(i, j, k) + Hy(im1_Hy, j, k) + Hy(i, jp1_Hy, k) + Hy(im1_Hy, jp1_Hy, k))/4.0
        elseif (field == iHz) then
          interp = (Hz(i, j, k) + Hz(im1_Hz, j, k) + Hz(i, j, kp1_Hz) + Hz(im1_Hz, j, kp1_Hz))/4.0
        end if
      elseif (atwhere == iHy) then
        if (field == iEx) then
          interp = (Ex(i, j, k) + Ex(i, j, kp1_Ex))/2.0
        elseif (field == iEy) then
          interp = (Ey(i, j, k) + Ey(ip1_Ey, j, k) + Ey(i, jm1_Ey, k) + Ey(ip1_Ey, jm1_Ey, k) + &
                    Ey(i, j, kp1_Ey) + Ey(ip1_Ey, j, kp1_Ey) + Ey(i, jm1_Ey, kp1_Ey) + Ey(ip1_Ey, jm1_Ey, kp1_Ey))/8.0
        elseif (field == iEz) then
          interp = (Ez(ip1_Ez, j, k) + Ez(i, j, k))/2.0
        elseif (field == iHy) then
          interp = Hy(i, j, k)
        elseif (field == iHx) then
          interp = (Hx(i, j, k) + Hx(i, jm1_Hx, k) + Hx(ip1_Hx, j, k) + Hx(ip1_Hx, jm1_Hx, k))/4.0
        elseif (field == iHz) then
          interp = (Hz(i, j, k) + Hz(i, jm1_Hz, k) + Hz(i, j, kp1_Hz) + Hz(i, jm1_Hz, kp1_Hz))/4.0
        end if
      elseif (atwhere == iHz) then
        if (field == iEx) then
          interp = (Ex(i, j, k) + Ex(i, jp1_Ex, k))/2.0
        elseif (field == iEy) then
          interp = (Ey(i, j, k) + Ey(ip1_Ey, j, k))/2.0
        elseif (field == iEz) then
          interp = (Ez(i, j, k) + Ez(i, jp1_Ez, k) + Ez(i, j, km1_Ez) + Ez(i, jp1_Ez, km1_Ez) + &
                    Ez(ip1_Ez, j, k) + Ez(ip1_Ez, jp1_Ez, k) + Ez(ip1_Ez, j, km1_Ez) + Ez(ip1_Ez, jp1_Ez, km1_Ez))/8.0
        elseif (field == iHz) then
          interp = Hz(i, j, k)
        elseif (field == iHx) then
          interp = (Hx(i, j, k) + Hx(ip1_Hx, j, k) + Hx(i, j, km1_Hx) + Hx(ip1_Hx, j, km1_Hx))/4.0
        elseif (field == iHy) then
          interp = (Hy(i, j, k) + Hy(i, jp1_Hy, k) + Hy(i, j, km1_Hy) + Hy(i, jp1_Hy, km1_Hy))/4.0
        end if
      end if

    end function interpolate_field_atwhere

  end module Observa

