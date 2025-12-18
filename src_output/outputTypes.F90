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

   integer, parameter :: UNDEFINED_DOMAIN = -1
   integer, parameter :: TIME_DOMAIN = 0
   integer, parameter :: FREQUENCY_DOMAIN = 1
   integer, parameter :: BOTH_DOMAIN = 2

   character(len=4), parameter :: datFileExtension = '.dat'
   character(len=4), parameter :: timeExtension = 'tm'
   character(len=4), parameter :: frequencyExtension = 'fq'

   type :: domain_t
      real(kind=RKIND_tiempo) :: tstart = 0.0_RKIND_tiempo, tstop = 0.0_RKIND_tiempo, tstep = 0.0_RKIND_tiempo
      real(kind=RKIND)        :: fstart = 0.0_RKIND, fstop = 0.0_RKIND, fstep
      integer(kind=SINGLE)    :: fnum = 0
      integer(kind=SINGLE)    :: domainType = UNDEFINED_DOMAIN
      logical                 :: logarithmicSpacing = .false.
   end type domain_t

   type spheric_domain_t
      real(kind=RKIND) :: phiStart = 0.0_RKIND, phiStop = 0.0_RKIND, phiStep = 0.0_RKIND
      real(kind=RKIND) :: thetaStart = 0.0_RKIND, thetaStop = 0.0_RKIND, thetastep = 0.0_RKIND
   end type

   type cell_coordinate_t
      integer(kind=SINGLE) :: x,y,z
   end type cell_coordinate_t

   type field_data_t
      real(kind=RKIND), pointer, dimension(:, :, :), contiguous :: x => NULL()
      real(kind=RKIND), pointer, dimension(:, :, :), contiguous :: y => NULL()
      real(kind=RKIND), pointer, dimension(:, :, :), contiguous :: z => NULL()
      real(kind=RKIND), pointer, dimension(:), contiguous :: deltaX => NULL()
      real(kind=RKIND), pointer, dimension(:), contiguous :: deltaY => NULL()
      real(kind=RKIND), pointer, dimension(:), contiguous :: deltaZ => NULL()
   end type field_data_t

   type fields_reference_t
      type(field_data_t) :: E
      type(field_data_t) :: H
   end type fields_reference_t

   type point_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE !reference and field
      type(domain_t) :: domain
      type(cell_coordinate_t) :: coordinates
      integer(kind=SINGLE) :: fileUnitTime, fileUnitFreq
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE, nFreq = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      real(kind=RKIND), dimension(BuffObse) :: valueForTime = 0.0_RKIND

      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      complex(kind=CKIND), dimension(:), allocatable :: valueForFreq
      complex(kind=CKIND), dimension(:), allocatable :: auxExp_E
      complex(kind=CKIND), dimension(:), allocatable :: auxExp_H
   end type point_probe_output_t

   type wire_charge_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE
      integer(kind=SINGLE) :: fileUnitTime
      type(domain_t) :: domain
      type(cell_coordinate_t) :: coordinates
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: chargeComponent
      integer(kind=SINGLE) :: sign = +1

      type(CurrentSegments), pointer :: segment

      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      real(kind=RKIND), dimension(BuffObse) :: chargeValue
   end type wire_charge_probe_output_t

   type current_values_t
      real(kind=RKIND) :: current = 0.0_RKIND, deltaVoltage = 0.0_RKIND
      real(kind=RKIND) :: plusVoltage = 0.0_RKIND, minusVoltage = 0.0_RKIND, voltageDiference = 0.0_RKIND
   end type

   type wire_current_probe_output_t
      integer(kind=SINGLE) :: columnas = 6_SINGLE !reference, corriente, -e*dl, vplus, vminus, vplus-vminus
      integer(kind=SINGLE) :: fileUnitTime
      type(domain_t) :: domain
      type(cell_coordinate_t) :: coordinates
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: currentComponent
      integer(kind=SINGLE) :: sign = +1

      type(CurrentSegments), pointer :: segment
#ifdef CompileWithBerengerWires
      type(TSegment), pointer  :: segmentBerenger
#endif
#ifdef CompileWithSlantedWires
      class(Segment), pointer  :: segmentSlanted
#endif

      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      type(current_values_t), dimension(BuffObse) :: currentValues
   end type wire_current_probe_output_t

   type bulk_current_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE !reference and field
      integer(kind=SINGLE) :: fileUnitTime
      type(domain_t) :: domain
      type(cell_coordinate_t) :: lowerBound
      type(cell_coordinate_t) :: upperBound
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      real(kind=RKIND), dimension(BuffObse) :: valueForTime = 0.0_RKIND

   end type bulk_current_probe_output_t

   type volumic_current_probe_t
      integer(kind=SINGLE) :: columnas = 4_SINGLE !reference and current components
      type(domain_t) :: domain
      type(cell_coordinate_t) :: lowerBound
      type(cell_coordinate_t) :: upperBound
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent

      !Intent storage order:
      !(:) == (timeinstance) => timeValue
      !(:,:) == (timeInstance, componentId) => escalar

      !Time Domain (requires first allocation)
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(:), allocatable :: timeStep
      real(kind=RKIND), dimension(:, :), allocatable :: xValueForTime
      real(kind=RKIND), dimension(:, :), allocatable :: yValueForTime
      real(kind=RKIND), dimension(:, :), allocatable :: zValueForTime

      !Intent storage order:
      !(:) == (frquencyinstance) => timeValue
      !(:,:) == (frquencyinstance, componentId) => escalar

      !Frequency Domain (requires first allocation)
      integer(kind=SINGLE) :: nFreq = 0_SINGLE
      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      complex(kind=CKIND), dimension(:, :), allocatable :: xValueForFreq
      complex(kind=CKIND), dimension(:, :), allocatable :: yValueForFreq
      complex(kind=CKIND), dimension(:, :), allocatable :: zValueForFreq
      complex(kind=CKIND), dimension(:), allocatable :: auxExp_E
      complex(kind=CKIND), dimension(:), allocatable :: auxExp_H

   end type volumic_current_probe_t

   type volumic_field_probe_output_t
   !!!!!Pending
   end type volumic_field_probe_output_t
   type line_integral_probe_output_t
    !!!!!Pending
   end type line_integral_probe_output_t
   type far_field_probe_output_t
      integer(kind=SINGLE) :: fileUnitFreq
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: columnas = 6_SINGLE !reference and current components
      type(domain_t) :: domain
      type(spheric_domain_t) :: sphericRange
      type(cell_coordinate_t) :: lowerBound
      type(cell_coordinate_t) :: upperBound
      character(len=BUFSIZE) :: path

      integer(kind=SINGLE) :: nMeasuredElements = 0_SINGLE
      integer(kind=SINGLE), dimension(:,:), allocatable :: coords
      integer(kind=SINGLE) :: nFreq = 0_SINGLE
      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      complex(kind=CKIND), dimension(:, :), allocatable :: valueForFreq
   end type far_field_probe_output_t
   type movie_probe_output_t
      integer(kind=SINGLE) :: PDVUnit
      integer(kind=SINGLE) :: columnas = 4_SINGLE !reference and current components
      type(domain_t) :: domain
      type(cell_coordinate_t) :: lowerBound
      type(cell_coordinate_t) :: upperBound
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent

      integer(kind=SINGLE) :: nMeasuredElements = 0_SINGLE
      integer(kind=SINGLE), dimension(:,:), allocatable :: coords

      !Intent storage order:
      !(:) == (timeinstance) => timeValue
      !(:,:) == (timeInstance, componentId) => escalar

      !Time Domain (requires first allocation)
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(:), allocatable :: timeStep
      real(kind=RKIND), dimension(:, :), allocatable :: xValueForTime
      real(kind=RKIND), dimension(:, :), allocatable :: yValueForTime
      real(kind=RKIND), dimension(:, :), allocatable :: zValueForTime
   end type movie_probe_output_t
   type frequency_slice_probe_output_t
      integer(kind=SINGLE) :: PDVUnit
      integer(kind=SINGLE) :: columnas = 4_SINGLE !reference and current components
      type(domain_t) :: domain
      type(cell_coordinate_t) :: lowerBound
      type(cell_coordinate_t) :: upperBound
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent

      integer(kind=SINGLE) :: nMeasuredElements = 0_SINGLE
      integer(kind=SINGLE), dimension(:,:), allocatable :: coords

      !Intent storage order:
      !(:) == (frquencyinstance) => timeValue
      !(:,:) == (frquencyinstance, componentId) => escalar

      !Frequency Domain (requires first allocation)
      integer(kind=SINGLE) :: nFreq = 0_SINGLE
      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      complex(kind=CKIND), dimension(:, :), allocatable :: xValueForFreq
      complex(kind=CKIND), dimension(:, :), allocatable :: yValueForFreq
      complex(kind=CKIND), dimension(:, :), allocatable :: zValueForFreq
      complex(kind=CKIND), dimension(:), allocatable :: auxExp_E
      complex(kind=CKIND), dimension(:), allocatable :: auxExp_H
   end type frequency_slice_probe_output_t

contains

end module outputTypes
