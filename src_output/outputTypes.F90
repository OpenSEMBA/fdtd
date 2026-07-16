module outputTypes_m
   use FDETYPES_m
   use HollandWires_m
   use wiresHolland_constants_m
#ifdef CompileWithBerengerWires
   use WiresBerenger
#endif
#ifdef CompileWithSlantedWires
   use WiresSlanted
   use WiresSlanted_Types
   use WiresSlanted_Constants
#endif
   implicit none

!=====================================================
! Parameters & constants
!=====================================================
   integer, parameter :: UNDEFINED_DOMAIN  = -1
   integer, parameter :: TIME_DOMAIN       =  0
   integer, parameter :: FREQUENCY_DOMAIN  =  1
   integer, parameter :: BOTH_DOMAIN       =  2


   character(len=4), parameter :: binaryExtension = '.bin'
   character(len=4), parameter :: pvdExtension = '.pvd'
   character(len=4), parameter :: datFileExtension = '.dat'
   character(len=4), parameter :: vtkFileExtension = '.vtk'
   character(len=4), parameter :: vtuFileExtension = '.vtu'
   character(len=2), parameter :: timeExtension    = 'tm'
    character(len=2), parameter :: frequencyExtension = 'fq'
    character(len=1), parameter :: wordseparation = '_'

    integer, parameter :: OUTPUT_ARTIFACT_UNDEFINED = 0
    integer, parameter :: OUTPUT_ARTIFACT_TEXT = 1
    integer, parameter :: OUTPUT_ARTIFACT_BINARY = 2
    integer, parameter :: OUTPUT_ARTIFACT_METADATA = 3
    integer, parameter :: OUTPUT_ARTIFACT_VISUALISATION_METADATA = 4
    integer, parameter :: OUTPUT_ARTIFACT_VISUALISATION_DATA = 5
    integer, parameter :: OUTPUT_ARTIFACT_GEOMETRY = 6

    integer, parameter :: OUTPUT_LIFECYCLE_DECLARED = 0
    integer, parameter :: OUTPUT_LIFECYCLE_ACTIVE = 1
    integer, parameter :: OUTPUT_LIFECYCLE_FINALISING = 2
    integer, parameter :: OUTPUT_LIFECYCLE_COMPLETE = 3
    integer, parameter :: OUTPUT_LIFECYCLE_FAILED = 4

    integer, parameter :: BINARY_ENDIAN_UNSPECIFIED = 0
    integer, parameter :: BINARY_ENDIAN_LITTLE = 1
    integer, parameter :: BINARY_ENDIAN_BIG = 2

    integer, parameter :: BINARY_NUMERIC_UNSPECIFIED = 0
    integer, parameter :: BINARY_NUMERIC_REAL32 = 1
    integer, parameter :: BINARY_NUMERIC_REAL64 = 2
    integer, parameter :: BINARY_NUMERIC_INT32 = 3
    integer, parameter :: BINARY_NUMERIC_INT64 = 4

    integer, parameter :: BINARY_COMPLEX_UNSPECIFIED = 0
    integer, parameter :: BINARY_COMPLEX_REAL_IMAG = 1
    integer, parameter :: BINARY_COMPLEX_MAGNITUDE_PHASE = 2

!=====================================================
! Basic helper / geometry types
!=====================================================
    type :: cell_coordinate_t
       integer(kind=SINGLE) :: x, y, z
    end type cell_coordinate_t

    type :: output_artifact_t
       integer :: kind = OUTPUT_ARTIFACT_UNDEFINED
       character(len=BUFSIZE) :: relative_path = ''
       logical :: required = .true.
       integer :: byte_order = BINARY_ENDIAN_UNSPECIFIED
       integer :: numeric_representation = BINARY_NUMERIC_UNSPECIFIED
       integer :: complex_representation = BINARY_COMPLEX_UNSPECIFIED
       integer(kind=8) :: record_bytes = 0
    end type output_artifact_t

     type :: output_lifecycle_t
        integer :: state = OUTPUT_LIFECYCLE_DECLARED
        character(len=BUFSIZE) :: diagnostic = ''
     end type output_lifecycle_t

     type :: probe_output_ownership_t
        integer, allocatable :: participant_ranks(:)
        integer :: scalar_writer_rank = -1
     end type probe_output_ownership_t

     type :: probe_metadata_t
       integer :: schema_version = 1
       character(len=BUFSIZE) :: probe_id = ''
       character(len=BUFSIZE) :: quantity = ''
       type(cell_coordinate_t) :: lower_bound
       type(cell_coordinate_t) :: upper_bound
        integer :: domain_type = UNDEFINED_DOMAIN
        type(output_lifecycle_t) :: lifecycle
        type(probe_output_ownership_t) :: ownership
        type(output_artifact_t), allocatable :: artifacts(:)
     end type probe_metadata_t

   type :: domain_t
      real(kind=RKIND_tiempo) :: tstart = 0.0_RKIND_tiempo
      real(kind=RKIND_tiempo) :: tstop  = 0.0_RKIND_tiempo
      real(kind=RKIND_tiempo) :: tstep  = 0.0_RKIND_tiempo
      real(kind=RKIND)        :: fstart = 0.0_RKIND
      real(kind=RKIND)        :: fstop  = 0.0_RKIND
      real(kind=RKIND)        :: fstep
      integer(kind=SINGLE)    :: fnum = 0
      integer(kind=SINGLE)    :: domainType = UNDEFINED_DOMAIN
      logical                 :: logarithmicSpacing = .false.
   end type domain_t

   type :: spheric_domain_t
      real(kind=RKIND) :: phiStart   = 0.0_RKIND
      real(kind=RKIND) :: phiStop    = 0.0_RKIND
      real(kind=RKIND) :: phiStep    = 0.0_RKIND
      real(kind=RKIND) :: thetaStart = 0.0_RKIND
      real(kind=RKIND) :: thetaStop  = 0.0_RKIND
      real(kind=RKIND) :: thetastep  = 0.0_RKIND
   end type spheric_domain_t

!=====================================================
! Field & current data containers
!=====================================================
   type :: field_data_t
      real(kind=RKIND), pointer, dimension(:, :, :), contiguous :: x => NULL()
      real(kind=RKIND), pointer, dimension(:, :, :), contiguous :: y => NULL()
      real(kind=RKIND), pointer, dimension(:, :, :), contiguous :: z => NULL()
      real(kind=RKIND), pointer, dimension(:), contiguous :: deltaX => NULL()
      real(kind=RKIND), pointer, dimension(:), contiguous :: deltaY => NULL()
      real(kind=RKIND), pointer, dimension(:), contiguous :: deltaZ => NULL()
   end type field_data_t

   type :: fields_reference_t
      type(field_data_t) :: E
      type(field_data_t) :: H
   end type fields_reference_t

   type :: current_values_t
      real(kind=RKIND) :: current = 0.0_RKIND
      real(kind=RKIND) :: deltaVoltage = 0.0_RKIND
      real(kind=RKIND) :: plusVoltage  = 0.0_RKIND
      real(kind=RKIND) :: minusVoltage = 0.0_RKIND
      real(kind=RKIND) :: voltageDiference = 0.0_RKIND
   end type current_values_t

!=====================================================
! Abstract probe hierarchy
!=====================================================
    type :: abstract_probe_t
       integer(kind=SINGLE)      :: columnas
       type(domain_t)            :: domain
       type(cell_coordinate_t)   :: mainCoords
       integer(kind=SINGLE)      :: component
       character(len=BUFSIZE)    :: path
   end type abstract_probe_t
   type, extends(abstract_probe_t) :: abstract_time_probe_t
      character(len=BUFSIZE) :: filePathTime
      integer(kind=SINGLE) :: nTime = 0_SINGLE
      integer(kind=SINGLE) :: nTimesFlushed = 0_SINGLE !times alredy writen in disk
      real(kind=RKIND_tiempo), allocatable :: timeStep(:)
   end type abstract_time_probe_t

   type, extends(abstract_probe_t) :: abstract_frequency_probe_t
      character(len=BUFSIZE) :: filePathFreq
      integer(kind=SINGLE) :: nFreq = 0_SINGLE
      real(kind=RKIND), allocatable    :: frequencySlice(:)
      complex(kind=CKIND), allocatable :: auxExp_E(:), auxExp_H(:)
   end type abstract_frequency_probe_t

   type, extends(abstract_probe_t) :: abstract_time_frequency_probe_t
      character(len=BUFSIZE) :: filePathTime, filePathFreq
      integer(kind=SINGLE) :: nTime = 0_SINGLE, nFreq = 0_SINGLE
      real(kind=RKIND_tiempo), allocatable :: timeStep(:)
      real(kind=RKIND), allocatable        :: frequencySlice(:)
      complex(kind=CKIND), allocatable     :: auxExp_E(:), auxExp_H(:)
   end type abstract_time_frequency_probe_t

!=====================================================
! Concrete probe types
!=====================================================
   type, extends(abstract_probe_t) :: mapvtk_output_t
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE), allocatable :: coords(:, :)
      integer(kind=SINGLE), allocatable :: currentType(:)
      integer(kind=SINGLE), allocatable :: materialTag(:)
       integer :: nPoints = -1
       type(output_artifact_t) :: artifacts(2)
   end type mapvtk_output_t
   type, extends(abstract_time_frequency_probe_t) :: point_probe_output_t
      real(kind=RKIND), allocatable :: valueForTime(:)
      complex(kind=CKIND), allocatable :: valueForFreq(:)
      type(output_artifact_t) :: artifacts(2)
   end type point_probe_output_t

    type, extends(abstract_time_probe_t) :: wire_charge_probe_output_t
      integer(kind=SINGLE) :: sign = +1
      real(kind=RKIND), allocatable :: chargeValue(:)
       type(CurrentSegments_t), pointer :: segment
       type(output_artifact_t) :: artifacts(2)
   end type wire_charge_probe_output_t

    type, extends(abstract_time_probe_t) :: wire_current_probe_output_t
      integer(kind=SINGLE) :: sign = +1
      type(current_values_t) :: currentValues(BuffObse)
       type(CurrentSegments_t), pointer :: segment
       type(output_artifact_t) :: artifacts(2)
#ifdef CompileWithBerengerWires
      type(TSegment), pointer :: segmentBerenger
#endif
#ifdef CompileWithSlantedWires
      class(Segment), pointer :: segmentSlanted
#endif
   end type wire_current_probe_output_t

    type, extends(abstract_time_probe_t) :: bulk_current_probe_output_t
      !Binary format: timeStamp, Val. Total register size: 16
      type(cell_coordinate_t) :: auxCoords
       real(kind=RKIND), allocatable :: valueForTime(:)
       type(output_artifact_t) :: artifacts(2)
   end type bulk_current_probe_output_t

    type, extends(abstract_frequency_probe_t) :: far_field_probe_output_t
      type(spheric_domain_t)  :: sphericRange
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE)    :: nPoints = -1
      integer(kind=SINGLE), allocatable :: coords(:, :)
       complex(kind=CKIND), allocatable :: valueForFreq(:, :)
       type(output_artifact_t), allocatable :: artifacts(:)
   end type far_field_probe_output_t

   type, extends(abstract_time_probe_t) :: movie_probe_output_t
      !Binary format: timeStamp, x, y, z, xVal, yVal, zVal. Total register size: 44
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE)    :: nPoints = -1
      integer(kind=SINGLE), allocatable :: coords(:, :)     !(3, coordIdx)
      real(kind=RKIND), allocatable :: xValueForTime(:, :)  !(time, coordIdx)
      real(kind=RKIND), allocatable :: yValueForTime(:, :)  !(time, coordIdx)
      real(kind=RKIND), allocatable :: zValueForTime(:, :)  !(time, coordIdx)
      character(len=BUFSIZE) :: filesPath
   end type movie_probe_output_t

   type, extends(abstract_frequency_probe_t) :: frequency_slice_probe_output_t
      !Binary format: frequencySlice, x, y, z, xVal, yVal, zVal. Total register size: 44
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE)    :: nPoints = -1
      integer(kind=SINGLE), allocatable :: coords(:, :)        !(3, coordIdx)
      complex(kind=CKIND), allocatable :: xValueForFreq(:, :)  !(time, coordIdx)
      complex(kind=CKIND), allocatable :: yValueForFreq(:, :)  !(time, coordIdx)
      complex(kind=CKIND), allocatable :: zValueForFreq(:, :)  !(time, coordIdx)
      character(len=BUFSIZE) :: filesPath
   end type frequency_slice_probe_output_t

!=====================================================
! High-level aggregation types
!=====================================================
   type :: solver_output_t
      integer(kind=SINGLE) :: outputID = -1
      type(point_probe_output_t), allocatable :: pointProbe
      type(wire_current_probe_output_t), allocatable :: wireCurrentProbe
      type(wire_charge_probe_output_t), allocatable  :: wireChargeProbe
      type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe
      type(movie_probe_output_t), allocatable         :: movieProbe
      type(frequency_slice_probe_output_t), allocatable :: frequencySliceProbe
      type(far_field_probe_output_t), allocatable     :: farFieldOutput
      type(mapvtk_output_t), allocatable              :: mapvtkOutput
#ifdef CompileWithMPI
      integer(kind=4) :: MPISubcomm, MPIRoot, MPIGroupIndex
      integer(kind=4) :: ZIorig, ZEorig
#endif
   end type solver_output_t

   type :: problem_info_t
      type(media_matrices_t), pointer :: geometryToMaterialData
      type(limit_t), pointer :: problemDimension(:)
      type(bounds_t), pointer :: simulationBounds
      type(MediaData_t), pointer :: materialList(:)
      type(taglist_t), pointer :: materialTag
      real(RKIND), pointer, dimension(:) :: xSteps
      real(RKIND), pointer, dimension(:) :: ySteps
      real(RKIND), pointer, dimension(:) :: zSteps
   end type problem_info_t

contains

   subroutine declare_output_artifacts(metadata, paths, kinds)
      type(probe_metadata_t), intent(inout) :: metadata
      character(len=*), intent(in) :: paths(:)
      integer, intent(in) :: kinds(:)
      integer :: i

      if (allocated(metadata%artifacts)) deallocate(metadata%artifacts)
      allocate(metadata%artifacts(size(paths)))
      do i = 1, size(paths)
         metadata%artifacts(i)%kind = kinds(i)
         metadata%artifacts(i)%relative_path = paths(i)
      end do
      metadata%lifecycle%state = OUTPUT_LIFECYCLE_DECLARED
   end subroutine declare_output_artifacts

   subroutine declare_probe_artifacts(artifacts, paths, kinds)
      type(output_artifact_t), intent(out) :: artifacts(:)
      character(len=*), intent(in) :: paths(:)
      integer, intent(in) :: kinds(:)
      integer :: i

      artifacts = output_artifact_t()
      do i = 1, size(paths)
         artifacts(i)%kind = kinds(i)
         artifacts(i)%relative_path = paths(i)
      end do
   end subroutine declare_probe_artifacts

   pure logical function output_lifecycle_is_terminal(lifecycle)
      type(output_lifecycle_t), intent(in) :: lifecycle

      output_lifecycle_is_terminal = lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE .or. &
                                     lifecycle%state == OUTPUT_LIFECYCLE_FAILED
   end function output_lifecycle_is_terminal

   pure logical function probe_metadata_is_complete(metadata)
      type(probe_metadata_t), intent(in) :: metadata
      integer :: i

      probe_metadata_is_complete = metadata%lifecycle%state == OUTPUT_LIFECYCLE_COMPLETE
      if (.not. probe_metadata_is_complete) return

      if (.not. allocated(metadata%artifacts)) return
      do i = 1, size(metadata%artifacts)
         if (metadata%artifacts(i)%required .and. metadata%artifacts(i)%kind == OUTPUT_ARTIFACT_UNDEFINED) then
            probe_metadata_is_complete = .false.
            return
         end if
      end do
   end function probe_metadata_is_complete

end module outputTypes_m
