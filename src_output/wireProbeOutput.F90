module mod_wireProbeOutput
   use FDETYPES
   USE mod_UTILS
   use Report
   use outputTypes
   use mod_outputUtils
   use wiresHolland_constants
   use HollandWires

   implicit none
   private

   !===========================
   ! Public interface
   !===========================
   public :: init_wire_current_probe_output
   public :: init_wire_charge_probe_output
   public :: create_wire_current_probe_output
   public :: create_wire_charge_probe_output
   public :: update_wire_current_probe_output
   public :: update_wire_charge_probe_output
   public :: flush_wire_current_probe_output
   public :: flush_wire_charge_probe_output
   !===========================

   !===========================
   ! Private interface
   !===========================
   private :: find_current_segment
   private :: find_charge_segment
   private :: build_output_path
   private :: probe_bounds_extension
   private :: clear_current_time_data
   private :: clear_charge_time_data
   private :: update_current_holland

#ifdef CompileWithBerengerWires
   private :: update_current_berenger
#endif

#ifdef CompileWithSlantedWires
   private :: update_current_slanted
#endif
   !===========================

   contains
   !======================================================================
   ! INITIALIZATION
   !======================================================================
   subroutine init_wire_current_probe_output(this, coordinates, node, field, domain, media, &
                                             outputTypeExtension, mpidir, wiresflavor)
      type(wire_current_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in)            :: coordinates
      type(domain_t), intent(in)                     :: domain
      type(MediaData_t), intent(in)                  :: media(:)
      integer(kind=SINGLE), intent(in)               :: node, field, mpidir
      character(len=*), intent(in)                   :: outputTypeExtension, wiresflavor

      this%mainCoords = coordinates
      this%component = field
      this%domain    = domain
      this%sign      = 1

      call find_current_segment(this, node, field, media, wiresflavor)
      this%path = build_output_path(outputTypeExtension, field, node, mpidir, coordinates)

      call alloc_and_init(this%timeStep, BuffObse, 0.0_RKIND_tiempo)

   end subroutine init_wire_current_probe_output


   subroutine init_wire_charge_probe_output(this, coordinates, node, field, domain, &
                                            outputTypeExtension, mpidir, wiresflavor)
      type(wire_charge_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in)           :: coordinates
      type(domain_t), intent(in)                    :: domain
      integer(kind=SINGLE), intent(in)              :: node, field, mpidir
      character(len=*), intent(in)                  :: outputTypeExtension, wiresflavor

      this%mainCoords = coordinates
      this%component = field
      this%domain    = domain
      this%sign      = 1

      call find_charge_segment(this, node, field, wiresflavor)
      this%path = build_output_path(outputTypeExtension, field, node, mpidir, coordinates)

      call alloc_and_init(this%timeStep, BuffObse, 0.0_RKIND_tiempo)
      call alloc_and_init(this%chargeValue, BuffObse, 0.0_RKIND)

   end subroutine init_wire_charge_probe_output

   !======================================================================
   ! FILE CREATION
   !======================================================================
   subroutine create_wire_current_probe_output(this)
      type(wire_current_probe_output_t), intent(inout) :: this
      integer(kind=SINGLE) :: err
      call create_or_clear_file(trim(this%path)//'_'//timeExtension//'_'//datFileExtension, &
                                this%fileUnitTime, err)
   end subroutine

   subroutine create_wire_charge_probe_output(this)
      type(wire_charge_probe_output_t), intent(inout) :: this
      integer(kind=SINGLE) :: err
      call create_or_clear_file(trim(this%path)//'_'//timeExtension//'_'//datFileExtension, &
                                this%fileUnitTime, err)
   end subroutine

   !======================================================================
   ! UPDATE
   !======================================================================
   subroutine update_wire_current_probe_output(this, step, control, InvEps, InvMu)
      type(wire_current_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(sim_control_t), intent(in)     :: control
      real(kind=RKIND), intent(in)        :: InvEps(:), InvMu(:)

      this%nTime = this%nTime + 1
      this%timeStep(this%nTime) = step

      select case (trim(control%wiresflavor))
      case ('holland','transition')
         call update_current_holland(this, control, InvEps, InvMu)
#ifdef CompileWithBerengerWires
      case ('berenger')
         call update_current_berenger(this, InvEps, InvMu)
#endif
#ifdef CompileWithSlantedWires
      case ('slanted','semistructured')
         call update_current_slanted(this)
#endif
      end select
   end subroutine


   subroutine update_wire_charge_probe_output(this, step)
      type(wire_charge_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step

      this%nTime = this%nTime + 1
      this%timeStep(this%nTime) = step
      this%chargeValue(this%nTime) = this%segment%ChargeMinus%ChargePresent
   end subroutine

   !======================================================================
   ! FLUSH
   !======================================================================
   subroutine flush_wire_current_probe_output(this)
      type(wire_current_probe_output_t), intent(inout) :: this
      integer :: i

      open(this%fileUnitTime, file=trim(this%path)//'_'//timeExtension//'_'//datFileExtension, &
           status='old', position='append')

      do i = 1, this%nTime
         write(this%fileUnitTime, fmt) this%timeStep(i), &
            this%currentValues(i)%current, &
            this%currentValues(i)%deltaVoltage, &
            this%currentValues(i)%plusVoltage, &
            this%currentValues(i)%minusVoltage, &
            this%currentValues(i)%voltageDiference
      end do
      close(this%fileUnitTime)

      call clear_current_time_data(this)
   end subroutine


   subroutine flush_wire_charge_probe_output(this)
      type(wire_charge_probe_output_t), intent(inout) :: this
      integer :: i

      open(this%fileUnitTime, file=trim(this%path)//'_'//timeExtension//'_'//datFileExtension, &
           status='old', position='append')

      do i = 1, this%nTime
         write(this%fileUnitTime, fmt) this%timeStep(i), this%chargeValue(i)
      end do
      close(this%fileUnitTime)

      call clear_charge_time_data(this)
   end subroutine

   subroutine find_current_segment(this, node, field, media, wiresflavor)
      type(wire_current_probe_output_t), intent(inout) :: this
      type(MediaData_t), intent(in) :: media(:)
      integer(kind=SINGLE), intent(in) :: node, field
      character(len=*), intent(in) :: wiresflavor

      type(Thinwires_t), pointer :: Hwireslocal
      type(CurrentSegments), pointer :: seg
#ifdef CompileWithBerengerWires
      type(TWires), pointer :: Hwireslocal_B
#endif
#ifdef CompileWithSlantedWires
      type(WiresData), pointer :: Hwireslocal_S
#endif

      integer :: n, iwi, iwj, node2
      logical :: found
      character(len=BUFSIZE) :: buff

      found = .false.
      this%sign = 1

      select case (trim(adjustl(wiresflavor)))
      case ('holland','transition')
         Hwireslocal => GetHwires()
         this%segment => Hwireslocal%NullSegment

         do n = 1, Hwireslocal%NumCurrentSegments
            seg => Hwireslocal%CurrentSegment(n)
            if (seg%origindex == node .and. &
                seg%i == iCoord .and. seg%j == jCoord .and. seg%k == kCoord .and. &
                seg%tipofield*10 == field) then
               found = .true.
               this%segment => seg
               if (seg%orientadoalreves) this%sign = -1
               exit
            end if
         end do

#ifdef CompileWithBerengerWires
      case ('berenger')
         Hwireslocal_B => GetHwires_Berenger()
         do n = 1, Hwireslocal_B%NumSegments
            if (Hwireslocal_B%Segments(n)%IndexSegment == node) then
               found = .true.
               this%segmentBerenger => Hwireslocal_B%Segments(n)
               if (Hwireslocal_B%Segments(n)%orientadoalreves) this%sign = -1
               exit
            end if
         end do
#endif

#ifdef CompileWithSlantedWires
      case ('slanted','semistructured')
         Hwireslocal_S => GetHwires_Slanted()
         do n = 1, Hwireslocal_S%NumSegments
            if (Hwireslocal_S%Segments(n)%ptr%Index == node) then
               found = .true.
               this%segmentSlanted => Hwireslocal_S%Segments(n)%ptr
               exit
            end if
         end do
#endif
      end select

      ! --- multirabo fallback (Holland only)
      if (.not. found .and. trim(adjustl(wiresflavor)) /= 'berenger') then
         buscarabono: do iwi = 1, Hwireslocal%NumDifferentWires
            do iwj = 1, media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%numsegmentos
               if (node == media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex .and. &
                   media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then

                  node2 = media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE
                  do n = 1, Hwireslocal%NumCurrentSegments
                     seg => Hwireslocal%CurrentSegment(n)
                     if (seg%origindex == node2) then
                        found = .true.
                        this%segment => seg
                        if (seg%orientadoalreves) this%sign = -1
                        exit buscarabono
                     end if
                  end do
               end if
            end do
         end do buscarabono
      end if

      if (.not. found) then
         write(buff,'(a,4i7,a)') 'ERROR: WIRE probe ',node,iCoord,jCoord,kCoord,' DOES NOT EXIST'
         call WarnErrReport(buff,.true.)
      end if
   end subroutine find_current_segment

   subroutine find_charge_segment(this, node, field, wiresflavor)
      type(wire_charge_probe_output_t), intent(inout) :: this
      integer(kind=SINGLE), intent(in) :: node, field
      character(len=*), intent(in) :: wiresflavor

      type(Thinwires_t), pointer :: Hwireslocal
      type(CurrentSegments), pointer :: seg
      integer :: n
      logical :: found
      character(len=BUFSIZE) :: buff

      found = .false.
      this%sign = 1

      if (trim(adjustl(wiresflavor)) /= 'holland' .and. &
          trim(adjustl(wiresflavor)) /= 'transition') then
         call WarnErrReport('Charge probes only supported for holland wires', .true.)
         return
      end if

      Hwireslocal => GetHwires()

      do n = 1, Hwireslocal%NumCurrentSegments
         seg => Hwireslocal%CurrentSegment(n)
         if (seg%origindex == node .and. &
             seg%i == iCoord .and. seg%j == jCoord .and. seg%k == kCoord .and. &
             seg%tipofield*10000 == field) then
            found = .true.
            this%segment => seg
            if (seg%orientadoalreves) this%sign = -1
            exit
         end if
      end do

      if (.not. found) then
         write(buff,'(a,4i7,a)') 'ERROR: CHARGE probe ',node,iCoord,jCoord,kCoord,' DOES NOT EXIST'
         call WarnErrReport(buff,.true.)
      end if
   end subroutine find_charge_segment

   function probe_bounds_extension(mpidir, coords) result(ext)
      integer(kind=SINGLE), intent(in) :: mpidir
      type(cell_coordinate_t), intent(in) :: coords
      character(len=BUFSIZE) :: ext
      character(len=BUFSIZE) :: ci, cj, ck

      write(ci,'(i7)') coords%x
      write(cj,'(i7)') coords%y
      write(ck,'(i7)') coords%z

#if CompileWithMPI
      select case (mpidir)
      case (3)
         ext = trim(adjustl(ci))//'_'//trim(adjustl(cj))//'_'//trim(adjustl(ck))
      case (2)
         ext = trim(adjustl(cj))//'_'//trim(adjustl(ck))//'_'//trim(adjustl(ci))
      case (1)
         ext = trim(adjustl(ck))//'_'//trim(adjustl(ci))//'_'//trim(adjustl(cj))
      case default
         call stoponerror(0,0,'Buggy error in mpidir.')
      end select
#else
      ext = trim(adjustl(ci))//'_'//trim(adjustl(cj))//'_'//trim(adjustl(ck))
#endif
   end function probe_bounds_extension

   function build_output_path(outExt, field, node, mpidir, coords) result(path)
      character(len=*), intent(in) :: outExt
      integer(kind=SINGLE), intent(in) :: field, node, mpidir
      type(cell_coordinate_t), intent(in) :: coords
      character(len=BUFSIZE) :: path
      character(len=BUFSIZE) :: nodeStr, fieldExt, boundsExt

      write(nodeStr,'(i7)') node
      fieldExt  = get_prefix_extension(field, mpidir)
      boundsExt = probe_bounds_extension(mpidir, coords)

      path = trim(outExt)//'_'//trim(fieldExt)//'_'// &
             trim(boundsExt)//'_s'//trim(adjustl(nodeStr))
   end function build_output_path

   subroutine clear_current_time_data(this)
      type(wire_current_probe_output_t), intent(inout) :: this

      this%timeStep = 0.0_RKIND_tiempo
      this%currentValues%current          = 0.0_RKIND
      this%currentValues%deltaVoltage     = 0.0_RKIND
      this%currentValues%plusVoltage      = 0.0_RKIND
      this%currentValues%minusVoltage     = 0.0_RKIND
      this%currentValues%voltageDiference = 0.0_RKIND
      this%nTime = 0
   end subroutine clear_current_time_data

   subroutine clear_charge_time_data(this)
      type(wire_charge_probe_output_t), intent(inout) :: this

      this%timeStep    = 0.0_RKIND_tiempo
      this%chargeValue = 0.0_RKIND
      this%nTime = 0
   end subroutine clear_charge_time_data

   subroutine update_current_holland(this, control, InvEps, InvMu)
      type(wire_current_probe_output_t), intent(inout) :: this
      type(sim_control_t), intent(in) :: control
      real(kind=RKIND), intent(in) :: InvEps(:), InvMu(:)

      type(CurrentSegments), pointer :: seg

      seg => this%segment

      this%currentValues(this%nTime)%current = &
         this%sign * seg%currentpast

      this%currentValues(this%nTime)%deltaVoltage = &
         - seg%Efield_wire2main * seg%delta

      if (control%wirecrank) then
         this%currentValues(this%nTime)%plusVoltage = this%sign * &
            (seg%ChargePlus%ChargePresent) * seg%Lind * &
            (InvMu(seg%indexmed) * InvEps(seg%indexmed))

         this%currentValues(this%nTime)%minusVoltage = this%sign * &
            (seg%ChargeMinus%ChargePresent) * seg%Lind * &
            (InvMu(seg%indexmed) * InvEps(seg%indexmed))
      else
         this%currentValues(this%nTime)%plusVoltage = this%sign * &
            ((seg%ChargePlus%ChargePresent + seg%ChargePlus%ChargePast) / 2.0_RKIND) * &
            seg%Lind * (InvMu(seg%indexmed) * InvEps(seg%indexmed))

         this%currentValues(this%nTime)%minusVoltage = this%sign * &
            ((seg%ChargeMinus%ChargePresent + seg%ChargeMinus%ChargePast) / 2.0_RKIND) * &
            seg%Lind * (InvMu(seg%indexmed) * InvEps(seg%indexmed))
      end if

      this%currentValues(this%nTime)%voltageDiference = &
         this%currentValues(this%nTime)%plusVoltage - &
         this%currentValues(this%nTime)%minusVoltage
   end subroutine update_current_holland

#ifdef CompileWithBerengerWires
   subroutine update_current_berenger(this, InvEps, InvMu)
      type(wire_current_probe_output_t), intent(inout) :: this
      real(kind=RKIND), intent(in) :: InvEps(:), InvMu(:)

      type(TSegment), pointer :: seg

      seg => this%segmentBerenger

      this%currentValues(this%nTime)%current = &
         this%sign * seg%currentpast

      this%currentValues(this%nTime)%deltaVoltage = &
         - seg%field * seg%dl

      this%currentValues(this%nTime)%plusVoltage = this%sign * &
         ((seg%ChargePlus + seg%ChargePlusPast) / 2.0_RKIND) * &
         seg%L * (InvMu(seg%imed) * InvEps(seg%imed))

      this%currentValues(this%nTime)%minusVoltage = this%sign * &
         ((seg%ChargeMinus + seg%ChargeMinusPast) / 2.0_RKIND) * &
         seg%L * (InvMu(seg%imed) * InvEps(seg%imed))

      this%currentValues(this%nTime)%voltageDiference = &
         this%currentValues(this%nTime)%plusVoltage - &
         this%currentValues(this%nTime)%minusVoltage
   end subroutine update_current_berenger
#endif

#ifdef CompileWithSlantedWires
   subroutine update_current_slanted(this)
      type(wire_current_probe_output_t), intent(inout) :: this

      class(Segment), pointer :: seg

      seg => this%segmentSlanted

      this%currentValues(this%nTime)%current = &
         seg%Currentpast

      this%currentValues(this%nTime)%deltaVoltage = &
         - seg%field * seg%dl

      this%currentValues(this%nTime)%plusVoltage = &
         (seg%Voltage(iPlus)%ptr%Voltage + &
          seg%Voltage(iPlus)%ptr%VoltagePast) / 2.0_RKIND

      this%currentValues(this%nTime)%minusVoltage = &
         (seg%Voltage(iMinus)%ptr%Voltage + &
          seg%Voltage(iMinus)%ptr%VoltagePast) / 2.0_RKIND

      this%currentValues(this%nTime)%voltageDiference = &
         this%currentValues(this%nTime)%plusVoltage - &
         this%currentValues(this%nTime)%minusVoltage
   end subroutine update_current_slanted
#endif

end module mod_wireProbeOutput
