module outputTypes
   use FDETYPES
   use HollandWires
   use wiresHolland_constants
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

   character(len=4), parameter :: datFileExtension = '.dat'
   character(len=4), parameter :: timeExtension    = 'tm'
   character(len=4), parameter :: frequencyExtension = 'fq'

!=====================================================
! Basic helper / geometry types
!=====================================================
   type :: cell_coordinate_t
      integer(kind=SINGLE) :: x, y, z
   end type cell_coordinate_t

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
      integer(kind=SINGLE) :: fileUnitTime
      integer(kind=SINGLE) :: nTime
      real(kind=RKIND_tiempo), allocatable :: timeStep(:)
   end type abstract_time_probe_t

   type, extends(abstract_probe_t) :: abstract_frequency_probe_t
      integer(kind=SINGLE) :: fileUnitFreq
      integer(kind=SINGLE) :: nFreq
      real(kind=RKIND), allocatable    :: frequencySlice(:)
      complex(kind=CKIND), allocatable :: auxExp_E(:), auxExp_H(:)
   end type abstract_frequency_probe_t

   type, extends(abstract_probe_t) :: abstract_time_frequency_probe_t
      integer(kind=SINGLE) :: fileUnitTime, fileUnitFreq
      integer(kind=SINGLE) :: nTime, nFreq
      real(kind=RKIND_tiempo), allocatable :: timeStep(:)
      real(kind=RKIND), allocatable        :: frequencySlice(:)
      complex(kind=CKIND), allocatable     :: auxExp_E(:), auxExp_H(:)
   end type abstract_time_frequency_probe_t

!=====================================================
! Concrete probe types
!=====================================================
   type, extends(abstract_time_frequency_probe_t) :: point_probe_output_t
      real(kind=RKIND), allocatable :: valueForTime(:)
      complex(kind=CKIND), allocatable :: valueForFreq(:)
   end type point_probe_output_t

   type, extends(abstract_time_probe_t) :: wire_charge_probe_output_t
      integer(kind=SINGLE) :: sign = +1
      real(kind=RKIND), allocatable :: chargeValue(:)
      type(CurrentSegments), pointer :: segment
   end type wire_charge_probe_output_t

   type, extends(abstract_time_probe_t) :: wire_current_probe_output_t
      integer(kind=SINGLE) :: sign = +1
      type(current_values_t) :: currentValues(BuffObse)
      type(CurrentSegments), pointer :: segment
#ifdef CompileWithBerengerWires
      type(TSegment), pointer :: segmentBerenger
#endif
#ifdef CompileWithSlantedWires
      class(Segment), pointer :: segmentSlanted
#endif
   end type wire_current_probe_output_t

   type, extends(abstract_time_probe_t) :: bulk_current_probe_output_t
      type(cell_coordinate_t) :: auxCoords
      real(kind=RKIND), allocatable :: valueForTime(:)
   end type bulk_current_probe_output_t

   type, extends(abstract_frequency_probe_t) :: far_field_probe_output_t
      type(spheric_domain_t)  :: sphericRange
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE)    :: nPoints
      integer(kind=SINGLE), allocatable :: coords(:, :)
      complex(kind=CKIND), allocatable :: valueForFreq(:, :)
   end type far_field_probe_output_t

   type, extends(abstract_time_probe_t) :: movie_probe_output_t
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE)    :: nPoints
      integer(kind=SINGLE), allocatable :: coords(:, :)
      real(kind=RKIND), allocatable :: xValueForTime(:, :)
      real(kind=RKIND), allocatable :: yValueForTime(:, :)
      real(kind=RKIND), allocatable :: zValueForTime(:, :)
   end type movie_probe_output_t

   type, extends(abstract_frequency_probe_t) :: frequency_slice_probe_output_t
      type(cell_coordinate_t) :: auxCoords
      integer(kind=SINGLE)    :: nPoints
      integer(kind=SINGLE), allocatable :: coords(:, :)
      complex(kind=CKIND), allocatable :: xValueForFreq(:, :)
      complex(kind=CKIND), allocatable :: yValueForFreq(:, :)
      complex(kind=CKIND), allocatable :: zValueForFreq(:, :)
   end type frequency_slice_probe_output_t

!=====================================================
! High-level aggregation types
!=====================================================
   type :: solver_output_t
      integer(kind=SINGLE) :: outputID
      type(point_probe_output_t), allocatable :: pointProbe
      type(wire_current_probe_output_t), allocatable :: wireCurrentProbe
      type(wire_charge_probe_output_t), allocatable  :: wireChargeProbe
      type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe
      type(movie_probe_output_t), allocatable         :: movieProbe
      type(frequency_slice_probe_output_t), allocatable :: frequencySliceProbe
      type(far_field_probe_output_t), allocatable     :: farFieldOutput
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
   end type problem_info_t

contains

end module outputTypes
