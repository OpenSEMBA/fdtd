module output
   use FDETYPES
   use mod_domain
   use mod_outputUtils
   implicit none
   character(len=4) :: datFileExtension = '.dat', timeExtension = 'tm', frequencyExtension = 'fq'
   integer(kind=SINGLE) :: MAX_SERIALIZED_COUNT = 500, FILE_UNIT = 400

   type solver_output_t
      type(point_probe_output_t), allocatable :: pointProbe
      type(wire_probe_output_t), allocatable :: wireProbe
      type(bulk_current_probe_output_t), allocatable :: bulkCurrentProbe
      type(far_field_t), allocatable :: farField
      type(time_movie_output_t), allocatable :: timeMovie
      type(frequency_slice_output_t), allocatable :: frequencySlice
   end type solver_output_t

   type point_probe_output_t
      integer(kind=SINGLE) :: columnas = 2_SINGLE
      type(domain_t) :: domain
      integer(kind=SINGLE) :: xCoord, yCoord, zCoord
      character(len=BUFSIZE) :: path
      integer(kind=SINGLE) :: fieldComponent
      integer(kind=SINGLE) :: serializedTimeSize = 0_SINGLE, nFreq = 0_SINGLE
      real(kind=RKIND_tiempo), dimension(MAX_SERIALIZED_COUNT), allocatable :: timeStep
      real(kind=RKIND), dimension(MAX_SERIALIZED_COUNT), allocatable :: valueForTime
      real(kind=RKIND), dimension(:), allocatable :: frequencySlice
      real(kind=CKIND), dimension(:), allocatable :: valueForFreq
   end type point_probe_output_t

   interface init_solver_output
      module procedure &
         init_point_probe_output, &
         init_wire_probe_output, &
         init_bulk_current_probe_output, &
         init_far_field, &
         initime_movie_output, &
         init_frequency_slice_output
   end interface

   interface update_solver_output
      module procedure &
         update_point_probe_output, &
         update_wire_probe_output, &
         update_bulk_current_probe_output, &
         update_far_field, &
         updateime_movie_output, &
         update_frequency_slice_output
   end interface

   interface flush_solver_output
      module procedure &
         flush_point_probe_output, &
         flush_wire_probe_output, &
         flush_bulk_current_probe_output, &
         flush_far_field, &
         flushime_movie_output, &
         flush_frequency_slice_output
   end interface

   interface delete_solver_output
      module procedure &
         delete_point_probe_output, &
         delete_wire_probe_output, &
         delete_bulk_current_probe_output, &
         delete_far_field, &
         deleteime_movie_output, &
         delete_frequency_slice_output
   end interface
contains

   subroutine init_point_probe_output(this, iCoord, jCoord, kCoord, field, domain, outputTypeExtension, mpidir)
      type(point_probe_output_t), intent(out) :: this
      integer(kind=SINGLE), intent(in) :: iCoord, jCoord, kCoord
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension

      this%xCoord = iCoord
      this%yCoord = jCoord
      this%zCoord = kCoord

      this%domain = domain
      this%path = get_output_path(outputTypeExtension, iCoord, jCoord, kCoord, field, mpidir)

      if (any(this%domain%domainType=(/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         this%nFreq = this%domain%fnum
         allocate (this%frequencySlice(this%domain%fnum))
         allocate (this%valueForFreq(this%domain%fnum))
      end if

   contains
      function get_output_path() result(outputPath)
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
            call stoponerror(layoutnumber, size, 'Buggy error in mpidir. ')
         end if
#else
         ext = trim(adjustl(chari))//'_'//trim(adjustl(charj))//'_'//trim(adjustl(chark))
#endif

         return
      end function get_probe_bounds_extension
   end subroutine init_point_probe_output

   subroutine update_point_probe_output(this, step)
      type(point_probe_output_t), intent(inout) :: this
      real(kind=RKIND), pointer, dimension(:, :, :) :: field
      real(kind=RKIND_tiempo) :: step

      field => get_field_component(this%fieldComponent)

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         this%serializedTimeSize = this%serializedTimeSize + 1
         this%timeStep(this%serializedTimeSize) = step
         this%valueForTime(this%serializedTimeSize) = field(i, j, k)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         do iter = 1, this%nFreq
            this%valueForFreq(iter) = &
               this%valueForFreq(iter) + field(i, j, k)*get_auxExp(this%frequencySlice(iter), this%fieldComponent)
         end do
      end if
   end subroutine update_point_probe_output

   subroutine flush_point_probe_output(this)
      type(point_probe_output_t), intent(inout) :: this

      integer(kind=SINGLE) :: timeUnitFile, frequencyUnitFile, status
      character(len=BUFSIZE) :: timeFileName, frequencyFileName

      if (any(this%domain%domainType == (/TIME_DOMAIN, BOTH_DOMAIN/))) then
         timeFileName = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         timeUnitFile = FILE_UNIT + 1

         status = open_file(timeUnitFile, timeFileName)
         if (status /= 0) call stoponerror()

         do i = 1, this%serializedTimeSize
            write (timeUnitFile, '(F12.4, 2X, F12.4)') this%timeStep(i), this%valueForTime(i)
         end do

         status = close_file(timeUnitFile)
      end if

      if (any(this%domain%domainType == (/FREQUENCY_DOMAIN, BOTH_DOMAIN/))) then
         frequencyFileName = trim(adjustl(this%path))//'_'//trim(adjustl(timeExtension))//'_'//trim(adjustl(datFileExtension))
         frequencyUnitFile = FILE_UNIT + 2

         OPEN (UNIT=frequencyUnitFile, FILE=frequencyFileName, STATUS='REPLACE', ACTION='WRITE', iostat=status)
         if (status /= 0) call stoponerror()

         do i = 1, this%nFreq
            write (frequencyUnitFile, '(F12.4, 2X, F12.4)') this%frequencySlice(i), this%valueForFreq(i)
         end do

         status = close_file(frequencyUnitFile)
      end if
   end subroutine flush_point_probe_output

end module output
