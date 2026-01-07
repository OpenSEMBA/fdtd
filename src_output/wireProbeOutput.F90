module mod_wireProbeOutput
   use FDETYPES
   use Report
   use outputTypes
   use mod_outputUtils
   use wiresHolland_constants
   use HollandWires

   implicit none

   private

   !===========================
   !  Public interface summary
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

contains
  subroutine init_wire_current_probe_output(this, coordinates, node, field, domain, media, outputTypeExtension, mpidir, wiresflavor)
      type(wire_current_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: node
      integer(kind=SINGLE), intent(in) :: field, mpidir
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      character(len=*), intent(in) :: wiresflavor
      type(domain_t), intent(in) :: domain
      type(MediaData_t), dimension(:), intent(in) :: media

      type(cell_coordinate_t) :: coordinates

      type(Thinwires_t), pointer  ::  Hwireslocal
#ifdef CompileWithBerengerWires
      type(TWires), pointer  ::  Hwireslocal_Berenger
#endif
#ifdef CompileWithSlantedWires
      type(WiresData), pointer  ::  Hwireslocal_Slanted
#endif

      select case (trim(adjustl(wiresflavor)))
      case ('holland', 'transition'); Hwireslocal => GetHwires()
#ifdef CompileWithBerengerWires
      case ('berenger'); Hwireslocal_Berenger => GetHwires_Berenger()
#endif
#ifdef CompileWithSlantedWires
      case ('slanted', 'semistructured'); Hwireslocal_Slanted => GetHwires_Slanted()
#endif
      end select

      call find_segment()

      this%mainCoords = coordinates

      this%component = field

      this%domain = domain
      this%path = get_output_path()

   contains
      subroutine find_segment()
         integer(kind=SINGLE) :: n, iwi, iwj, node2
         type(CurrentSegments), pointer :: currentSegment
         logical :: found = .false.
         character(len=BUFSIZE) :: buff

         select case (trim(adjustl(wiresflavor)))
         case ('holland', 'transition')
            this%segment => HWireslocal%NullSegment
            do n = 1, HWireslocal%NumCurrentSegments
               currentSegment => HWireslocal%CurrentSegment(n)
               if ((currentSegment%origindex == node) .and. &
                   (currentSegment%i == iCoord) .and. (currentSegment%j == jCoord) .and. (currentSegment%k == kCoord) .and. &
                   (currentSegment%tipofield*10 == field)) then
                  found = .true.
                  this%segment => currentSegment
                  if (currentSegment%orientadoalreves) this%sign = -1
               end if
            end do
#ifdef CompileWithBerengerWires
         case ('berenger')
            do n = 1, Hwireslocal_Berenger%NumSegments
               currentSegment => Hwireslocal_Berenger%Segments(n)
               if (currentSegment%IndexSegment == node) then
                  found = .true.
                  this%segmentBerenger => currentSegment
                  if (currentSegment%orientadoalreves) this%sign = -1
               end if
            end do
#endif
#ifdef CompileWithSlantedWires
         case ('slanted', 'semistructured')
            do n = 1, Hwireslocal_Slanted%NumSegments
               currentSegment => Hwireslocal_Slanted%Segments(n)
               if (currentSegment%ptr%Index == node) then
                  found = .true.
                  this%segmentSlanted => currentSegment%ptr
               end if
            end do
#endif
         end select

         if (.not. found) then
            select case (trim(adjustl(wiresflavor)))
            case ('holland', 'transition')
               buscarabono: do iwi = 1, Hwireslocal%NumDifferentWires
                  do iwj = 1, media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%numsegmentos
                     if ((node == media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%origindex) .and. &
                         media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multirabo) then
                        node2 = media(Hwireslocal%WireTipoMedio(iwi))%wire(1)%segm(iwj)%multiraboDE
                        do n = 1, HWireslocal%NumCurrentSegments
                           currentSegment => HWireslocal%CurrentSegment(n)
                           if (currentSegment%origindex == node2) then
                              found = .true.
                              this%segment => currentSegment
                              if (currentSegment%orientadoalreves) this%sign = -1
                           end if
                        end do
                        exit buscarabono
                     end if
                  end do
               end do buscarabono
#ifdef CompileWithSlantedWires
            case ('slanted', 'semistructured')
               do n = 1, Hwireslocal_Slanted%NumSegments
                  currentSegment => Hwireslocal_Slanted%Segments(n)
                  if (currentSegment%ptr%elotroindice == node) then
                     found = .true.
                     this%segmentSlanted => currentSegment%ptr
                  end if
               end do
#endif
            end select
         end if

         if (.not. found) then
            write (buff, '(a,4i7,a)') 'ERROR: WIRE probe ', node, iCoord, jCoord, kCoord, ' DOES NOT EXIST'
            CALL WarnErrReport(buff, .true.)
         end if
      end subroutine find_segment

      function get_output_path() result(outputPath)
         character(len=BUFSIZE) :: outputPath
         character(len=BUFSIZE)  ::  charNO
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension, prefixNodeExtension

         write (charNO, '(i7)') node
         prefixNodeExtension = 's'//trim(adjustl(charNO))
         probeBoundsExtension = get_probe_bounds_extension()
         prefixFieldExtension = get_prefix_extension(field, mpidir)

         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_' &
            //trim(adjustl(probeBoundsExtension))//'_'//trim(adjustl(prefixNodeExtension))
         return
      end function get_output_path

      function get_probe_bounds_extension() result(ext)
         character(len=BUFSIZE) :: ext
         character(len=BUFSIZE)  ::  chari, charj, chark

         write (chari, '(i7)') iCoord
         write (charj, '(i7)') jCoord
         write (chark, '(i7)') kCoord

#if CompileWithMPI
         if (mpidir == 3) then
            ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
         elseif (mpidir == 2) then
            ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
         elseif (mpidir == 1) then
            ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
         else
            call stoponerror('Buggy error in mpidir. ')
         end if
#else
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
#endif

         return
      end function get_probe_bounds_extension

   end subroutine init_wire_current_probe_output

   subroutine init_wire_charge_probe_output(this, coordinates, node, field, domain, outputTypeExtension, mpidir, wiresflavor)
      type(wire_charge_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) ::  node
      integer(kind=SINGLE), intent(in) :: field, mpidir
      character(len=*), intent(in) :: outputTypeExtension, wiresflavor
      type(domain_t), intent(in) :: domain

      type(Thinwires_t), pointer  ::  Hwireslocal
      type(CurrentSegments), pointer :: currentSegment
      type(cell_coordinate_t) :: coordinates
      character(len=BUFSIZE) :: buff
      integer(kind=SINGLE) :: n
      if (trim(adjustl(wiresflavor)) == 'holland' .or. trim(adjustl(wiresflavor)) == 'transition') Hwireslocal => GetHwires()

      call find_segment()

      this%mainCoords = coordinates

      this%component = field

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

      function get_output_path() result(outputPath)
         character(len=BUFSIZE) :: outputPath
         character(len=BUFSIZE)  ::  charNO
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension, prefixNodeExtension

         write (charNO, '(i7)') node
         prefixNodeExtension = 's'//trim(adjustl(charNO))
         probeBoundsExtension = get_probe_bounds_extension()
         prefixFieldExtension = get_prefix_extension(field, mpidir)

         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_' &
            //trim(adjustl(probeBoundsExtension))//'_'//trim(adjustl(prefixNodeExtension))
         return
      end function get_output_path

      function get_probe_bounds_extension() result(ext)
         character(len=BUFSIZE) :: ext
         character(len=BUFSIZE)  ::  chari, charj, chark

         write (chari, '(i7)') iCoord
         write (charj, '(i7)') jCoord
         write (chark, '(i7)') kCoord

#if CompileWithMPI
         if (mpidir == 3) then
            ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
         elseif (mpidir == 2) then
            ext = trim(adjustl(charj))//'_'//trim(adjustl(chark))//'_'//trim(adjustl(chari))
         elseif (mpidir == 1) then
            ext = trim(adjustl(chark))//'_'//trim(adjustl(chari))//'_'//trim(adjustl(charj))
         else
            call stoponerror('Buggy error in mpidir. ')
         end if
#else
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
#endif

         return
      end function get_probe_bounds_extension

   end subroutine init_wire_charge_probe_output

   subroutine create_wire_current_probe_output(this)
      type(wire_current_probe_output_t), intent(inout) :: this
      character(len=BUFSIZE) :: file_time
      integer(kind=SINGLE) :: err
      err = 0

      file_time = trim(adjustl(this%path))//'_'// &
                  trim(adjustl(timeExtension))//'_'// &
                  trim(adjustl(datFileExtension))
      call create_or_clear_file(file_time, this%fileUnitTime, err)
   end subroutine create_wire_current_probe_output

   subroutine create_wire_charge_probe_output(this)
      character(len=BUFSIZE) :: file_time
      type(wire_charge_probe_output_t), intent(inout) :: this
      integer(kind=SINGLE) :: err
      err = 0

      file_time = trim(adjustl(this%path))//'_'// &
                  trim(adjustl(timeExtension))//'_'// &
                  trim(adjustl(datFileExtension))
      call create_or_clear_file(file_time, this%fileUnitTime, err)
   end subroutine create_wire_charge_probe_output

   subroutine update_wire_current_probe_output(this, step, wiresflavor, wirecrank, InvEps, InvMu)
      type(wire_current_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      character(len=*), intent(in) :: wiresflavor
      logical, intent(in) :: wirecrank
      real(KIND=RKIND), pointer, dimension(:), intent(in) :: InvEps, InvMu

      type(CurrentSegments), pointer  ::  segmDumm
#ifdef CompileWithBerengerWires
      type(TSegment), pointer  ::  segmDumm_Berenger
#endif
#ifdef CompileWithSlantedWires
      class(Segment), pointer  ::  segmDumm_Slanted
#endif

      select case (trim(adjustl(wiresflavor)))
      case ('holland', 'transition')
         this%nTime = this%nTime + 1
         this%timeStep(this%nTime) = step
         SegmDumm => this%segment

         this%currentValues(this%nTime)%current = this%sign*SegmDumm%currentpast
         this%currentValues(this%nTime)%deltaVoltage = -SegmDumm%Efield_wire2main*SegmDumm%delta

         if (wirecrank) then
            this%currentValues(this%nTime)%plusVoltage = this%sign* &
                          (((SegmDumm%ChargePlus%ChargePresent)))*SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
            this%currentValues(this%nTime)%minusVoltage = this%sign* &
                         (((SegmDumm%ChargeMinus%ChargePresent)))*SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
         else
            this%currentValues(this%nTime)%plusVoltage = this%sign* &
                                               (((SegmDumm%ChargePlus%ChargePresent + SegmDumm%ChargePlus%ChargePast))/2.0_RKIND)* &
                                                                  SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
            this%currentValues(this%nTime)%minusVoltage = this%sign* &
                                             (((SegmDumm%ChargeMinus%ChargePresent + SegmDumm%ChargeMinus%ChargePast))/2.0_RKIND)* &
                                                                  SegmDumm%Lind*(InvMu(SegmDumm%indexmed)*InvEps(SegmDumm%indexmed))
         end if

         this%currentValues(this%nTime)%voltageDiference = &
            this%currentValues(this%nTime)%plusVoltage - this%currentValues(this%nTime)%minusVoltage

#ifdef CompileWithBerengerWires
      case ('berenger')
         this%nTime = this%nTime + 1
         this%timeStep(this%nTime) = step
         SegmDumm_Berenger => this%segmentBerenger

         this%currentValues(this%nTime)%current = this%sign*SegmDumm_Berenger%currentpast
         this%currentValues(this%nTime)%deltaVoltage = -SegmDumm_Berenger%field*SegmDumm_Berenger%dl

         this%currentValues(this%nTime)%plusVoltage = this%sign* &
                                                  (((SegmDumm_Berenger%ChargePlus + SegmDumm_Berenger%ChargePlusPast))/2.0_RKIND)* &
                                                  SegmDumm_Berenger%L*(InvMu(SegmDumm_Berenger%imed)*InvEps(SegmDumm_Berenger%imed))
         this%currentValues(this%nTime)%minusVoltage = this%sign* &
                                                (((SegmDumm_Berenger%ChargeMinus + SegmDumm_Berenger%ChargeMinusPast))/2.0_RKIND)* &
                                                  SegmDumm_Berenger%L*(InvMu(SegmDumm_Berenger%imed)*InvEps(SegmDumm_Berenger%imed))
         this%currentValues(this%nTime)%voltageDiference = &
            this%currentValues(this%nTime)%plusVoltage - this%currentValues(this%nTime)%minusVoltage

#endif
#ifdef CompileWithSlantedWires
      case ('slanted', 'semistructured')
         this%nTime = this%nTime + 1
         this%timeStep(this%nTime) = step
         SegmDumm_Slanted => this%segmentSlanted

         this%currentValues(this%nTime)%current = SegmDumm_Slanted%Currentpast !ojo: slanted ya los orienta bien y no hay que multiplicar por valorsigno
         this%currentValues(this%nTime)%deltaVoltage = -SegmDumm_Slanted%field*SegmDumm_Slanted%dl
         this%currentValues(this%nTime)%plusVoltage = &
            (((SegmDumm_Slanted%Voltage(iPlus)%ptr%Voltage + SegmDumm_Slanted%Voltage(iPlus)%ptr%VoltagePast))/2.0_RKIND)
         this%currentValues(this%nTime)%minusVoltage = &
            (((SegmDumm_Slanted%Voltage(iMinus)%ptr%Voltage + SegmDumm_Slanted%Voltage(iMinus)%ptr%VoltagePast))/2.0_RKIND)
         this%currentValues(this%nTime)%voltageDiference = &
            this%currentValues(this%nTime)%plusVoltage - this%currentValues(this%nTime)%minusVoltage
#endif
      end select

   end subroutine

   subroutine update_wire_charge_probe_output(this, step)
      type(wire_charge_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(CurrentSegments), pointer  ::  segmDumm

      this%nTime = this%nTime + 1
      this%timeStep(this%nTime) = step
      SegmDumm => this%segment
      this%chargeValue(this%nTime) = SegmDumm%ChargeMinus%ChargePresent
   end subroutine update_wire_charge_probe_output

   subroutine flush_wire_current_probe_output(this)
      type(wire_current_probe_output_t), intent(inout) :: this
      character(len=BUFSIZE) :: filename
      integer :: i

      filename = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
      open (unit=this%fileUnitTime, file=filename, status="old", action="write", position="append")

      do i = 1, this%nTime
            write (this%fileUnitTime, fmt) this%timeStep(i), &
            this%currentValues%current, &
            this%currentValues%deltaVoltage, &
            this%currentValues%plusVoltage, &
            this%currentValues%minusVoltage, &
            this%currentValues%voltageDiference
      end do
      close (this%fileUnitTime)

      call clear_time_data()
      contains 
      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_tiempo

         this%currentValues%current = 0.0_RKIND
         this%currentValues%deltaVoltage = 0.0_RKIND
         this%currentValues%plusVoltage = 0.0_RKIND
         this%currentValues%minusVoltage = 0.0_RKIND
         this%currentValues%voltageDiference = 0.0_RKIND

         this%nTime = 0
      end subroutine clear_time_data
   end subroutine flush_wire_current_probe_output

   subroutine flush_wire_charge_probe_output(this)
      type(wire_charge_probe_output_t), intent(inout) :: this
      character(len=BUFSIZE) :: filename
      integer :: i

      filename = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
      open (unit=this%fileUnitTime, file=filename, status="old", action="write", position="append")

      do i = 1, this%nTime
            write (this%fileUnitTime, fmt) this%timeStep(i), &
            this%chargeValue
      end do
      close (this%fileUnitTime)
            call clear_time_data()
      contains 
      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_tiempo

         this%chargeValue = 0.0_RKIND

         this%nTime = 0
      end subroutine clear_time_data
   end subroutine flush_wire_charge_probe_output
end module mod_wireProbeOutput
