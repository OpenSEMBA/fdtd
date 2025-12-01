module mod_wireChargeProbeOutput
  use FDETYPES
  use mod_domain
  use wiresHolland_constants
  use HollandWires
  implicit none
  type wire_charge_probe_output_t
      integer(kind=SINGLE) :: columnas = 6_SINGLE !reference, corriente, -e*dl, vplus, vminus, vplus-vminus
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: chargeComponent
      integer(kind=SINGLE) :: sign = +1

      type(CurrentSegments), pointer :: segment


      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      type(rkind), dimension(BuffObse) :: chargeValue
   end type wire_current_probe_output_t
contains

  subroutine init_wire_charge_probe_output(this, iCoord, jCoord, kCoord, node, field, domain, media, outputTypeExtension, mpidir,)
    type(wire_charge_probe_output_t), intent(out) :: this 
    integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord, node
    integer(kind=SINGLE), intent(in) :: field, mpidir
    character(len=BUFSIZE), intent(in) :: outputTypeExtension
    type(domain_t), intent(in) :: domain

    type(Thinwires_t), pointer  ::  Hwireslocal
    if (trim(adjustl(wiresflavor))=='holland' .or. trim(adjustl(wiresflavor))=='transition') Hwireslocal => GetHwires()

    call find_segment()

    this%xCoord = iCoord
    this%yCoord = jCoord
    this%zCoord = kCoord

    this%chargeComponent = field

    this%domain = domain
    this%path = get_output_path()

    contains
    subroutine find_segment()
      logical :: found = .false.
      do n = 1, HWireslocal%NumCurrentSegments
        currentSegment => HWireslocal%CurrentSegment(n)
        if ((currentSegment%origindex == node) .and. &
            (currentSegment%i == iCoord) .and. (currentSegment%j == jCoord) .and. (currentSegment%k == kCoord) .and. &
            (currentSegment%tipofield*10000 == field)) then
              found = .true.
              this%segment => currentSegment
              if (currentSegment%orientadoalreves) this%sign = -1
        end if
      end do
      if (.not. found) then
            write (buff, '(a,4i7,a)') 'ERROR: CHARGE probe ', node, iCoord, jCoord, kCoord, ' DOES NOT EXIST'
            CALL WarnErrReport(buff, .true.)
         end if
      end subroutine find_segment
  end subroutine init_wire_charge_probe_output

  subroutine update_wire_charge_probe_output(this, step)
    type(wire_charge_probe_output_t), intent(inout) :: this 
    real(kind=RKIND_tiempo), intent(in) :: step
    type(CurrentSegments), pointer  ::  segmDumm

    this%serializedTimeSize = this%serializedTimeSize + 1
    this%timeStep(this%serializedTimeSize) = step
    SegmDumm => this%segment
    this%chargeValue(this%serializedTimeSize) = SegmDumm%ChargeMinus%ChargePresent
  end subroutine update_wire_charge_probe_output
  
end module wire_charge_probe_output_t