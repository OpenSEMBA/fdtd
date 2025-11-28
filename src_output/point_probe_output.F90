module mod_pointProbeOutput
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   implicit none

   type point_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE !reference and field
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE, nFreq = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(BuffObse) :: timeStep = 0.0_RKIND
      real(kind=RKIND), dimension(BuffObse) :: valueForTime = 0.0_RKIND

      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      real(kind=CKIND), dimension(:), allocatable :: valueForFreq
   end type point_probe_output_t

contains
   subroutine init_point_probe_output(this, iCoord, jCoord, kCoord, field, domain, outputTypeExtension, mpidir)
      type(point_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      integer(kind=SINGLE) :: i

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%fieldComponent = field

      this%domain = domain
      this%path = get_output_path()

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         allocate (this%frequencySlice(this%domain%fnum))
         allocate (this%valueForFreq(this%domain%fnum))
         do i = 1, this%nFreq
            call init_frequency_slice(this%frequencySlice, this%domain)
         end do
         this%valueForFreq = (0.0_RKIND, 0.0_RKIND)
      end if

   contains
      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_probe_bounds_extension()
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//trim(adjustl(probeBoundsExtension))
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
   end subroutine init_point_probe_output

   subroutine update_point_probe_output(this, step, field)
      type(point_probe_output_t), intent(inout) :: this
      real(kind=RKIND), pointer, dimension(:, :, :) :: field
      real(kind=RKIND_tiempo) :: step
      integer(kind=SINGLE) :: iter

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%serializedTimeSize = this%serializedTimeSize + 1
         this%timeStep(this%serializedTimeSize) = step
         this%valueForTime(this%serializedTimeSize) = field(this%xCoord, this%yCoord, this%zCoord)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         do iter = 1, this%nFreq
            this%valueForFreq(iter) = &
   this%valueForFreq(iter) + field(this%xCoord, this%yCoord, this%zCoord) !*get_auxExp(this%frequencySlice(iter), this%fieldComponent)
         end do
      end if
   end subroutine update_point_probe_output

   subroutine flush_point_probe_output(this)
      type(point_probe_output_t), intent(inout) :: this

      integer(kind=SINGLE) :: timeUnitFile, frequencyUnitFile, status
      character(len=BUFSIZE) :: timeFileName, frequencyFileName
      integer(kind=SINGLE) :: i

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         timeFileName = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         timeUnitFile = FILE_UNIT + 1

         status = open_file(timeUnitFile, timeFileName)
         if (status /= 0) call stoponerror('Failed to open timeDomainFile. ')

         do i = 1, this%serializedTimeSize
            write (timeUnitFile, '(F12.4, 2X, F12.4)') this%timeStep(i), this%valueForTime(i)
         end do

         status = close_file(timeUnitFile)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         frequencyFileName = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         frequencyUnitFile = FILE_UNIT + 2

         OPEN (UNIT=frequencyUnitFile, FILE=frequencyFileName, STATUS='REPLACE', ACTION='WRITE', iostat=status)
         if (status /= 0) call stoponerror('Failed to open frequencyDomainFile. ')

         do i = 1, this%nFreq
            write (frequencyUnitFile, '(F12.4, 2X, F12.4)') this%frequencySlice(i), this%valueForFreq(i)
         end do

         status = close_file(frequencyUnitFile)
      end if
   end subroutine flush_point_probe_output

   subroutine delete_point_probe_output()
      !TODO
   end subroutine delete_point_probe_output
end module
