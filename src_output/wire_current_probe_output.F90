module mod_wireCurrentProbeOutput
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   use wiresHolland_constants
   use HollandWires
   implicit none

   type current_values_t
      real(kind=RKIND) :: current = 0.0_RKIND, deltaVoltage = 0.0_RKIND
      real(kind=RKIND) :: plusVoltage = 0.0_RKIND, minusVoltage = 0.0_RKIND, voltageDiference = 0.0_RKIND
   end type
   type wire_current_probe_output_t
      integer(kind=SINGLE) :: columnas = 6_SINGLE !reference, corriente, -e*dl, vplus, vminus, vplus-vminus
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: currentComponent
      integer(kind=SINGLE) :: sign = +1
      type(CurrentSegments), pointer :: segment

      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      type(current_values_t), dimension(BuffObse) :: currentValues
   end type wire_current_probe_output_t

contains
   subroutine init_wire_current_probe_output(this, iCoord, jCoord, kCoord, node, field, domain, media, outputTypeExtension, mpidir, wiresflavor)
      type(wire_current_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord, node
      integer(kind=SINGLE), intent(in) :: field, mpidir
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      character(len=*), intent(in) :: wiresflavor
      type(domain_t), intent(in) :: domain
      type(MediaData_t), pointer, dimension(:), intent(in) :: media

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

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%currentComponent = field

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


end module mod_wireCurrentProbeOutput
