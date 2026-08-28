module bulkProbeOutput_m
   use FDETYPES_m
   use utils_m
   use allocationUtils_m, only: alloc_and_init
   use outputTypes_m
   use outputUtils_m
#ifdef CompileWithMPI
   use mpi
#endif
   implicit none
   private

   public :: init_bulk_probe_output, update_bulk_probe_output, flush_bulk_probe_output

contains

   subroutine init_bulk_probe_output(this, lowerBound, upperBound, field, domain, outputTypeExtension, mpidir)
      type(bulk_current_probe_output_t), intent(out) :: this
      type(cell_coordinate_t), intent(in) :: lowerBound, upperBound
      integer(kind=SINGLE), intent(in) :: mpidir, field
      character(len=BUFSIZE), intent(in) :: outputTypeExtension
      type(domain_t), intent(in) :: domain

      character(len=BUFSIZE) :: artifact_paths(1)
      character(len=13) :: data_header
      integer :: artifact_kinds(1)
#ifdef CompileWithMPI
      integer :: mpi_rank, ierr
#endif

      this%mainCoords = lowerBound
      this%auxCoords = upperBound
      this%component = field

      this%domain = domain
      this%path = get_output_path()

#ifdef CompileWithMPI
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
      this%isWriter = mpi_rank == 0
#endif

      call alloc_and_init(this%timeStep, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND_tiempo)
      call alloc_and_init(this%valueForTime, OUTPUT_TIME_BUFFER_SIZE, 0.0_RKIND)
      artifact_paths(1) = trim(this%path)//'_'//timeExtension//datFileExtension
      artifact_kinds = OUTPUT_ARTIFACT_TEXT
      call declare_probe_artifacts(this%artifacts, artifact_paths, artifact_kinds)
      this%filePathTime = this%artifacts(1)%relative_path
      select case (field)
      case (iBloqueMx, iBloqueMy, iBloqueMz)
         data_header = 't circulation'
      case default
         data_header = 't current'
      end select
      if (this%isWriter) then
         call create_data_file(this%filePathTime, this%path, timeExtension, datFileExtension, data_header)
      end if

   contains

      function get_output_path() result(outputPath)
         character(len=BUFSIZE)  :: probeBoundsExtension, prefixFieldExtension
         character(len=BUFSIZE) :: outputPath
         probeBoundsExtension = get_coordinates_extension(this%mainCoords, this%auxCoords, mpidir)
         prefixFieldExtension = get_prefix_extension(field, mpidir)
         outputPath = &
            trim(adjustl(outputTypeExtension))//'_'//trim(adjustl(prefixFieldExtension))//'_'//trim(adjustl(probeBoundsExtension))
         return
      end function get_output_path

   end subroutine init_bulk_probe_output

   subroutine update_bulk_probe_output(this, step, field)
      type(bulk_current_probe_output_t), intent(inout) :: this
      real(kind=RKIND_tiempo), intent(in) :: step
      type(field_data_t), intent(in) :: field

      integer(kind=SINGLE) :: i1_m, i2_m, j1_m, j2_m, k1_m, k2_m
      integer(kind=SINGLE) :: i1, i2, j1, j2, k1, k2
      integer(kind=SINGLE) :: iii, jjj, kkk
#ifdef CompileWithMPI
      integer :: ierr
      real(kind=RKIND) :: localValue
#endif

      real(kind=RKIND), pointer, dimension(:, :, :) :: xF, yF, zF
      real(kind=RKIND), pointer, dimension(:) :: dx, dy, dz

      i1_m = this%mainCoords%x
      j1_m = this%mainCoords%y
      k1_m = this%mainCoords%z
      i2_m = this%auxCoords%x
      j2_m = this%auxCoords%y
      k2_m = this%auxCoords%z

      i1 = i1_m
      j1 = j1_m
      k1 = k1_m
      i2 = i2_m
      j2 = j2_m
      k2 = k2_m

      xF => field%x
      yF => field%y
      zF => field%z
      dx => field%deltaX
      dy => field%deltaY
      dz => field%deltaZ

      this%nTime = this%nTime + 1
      this%timeStep(this%nTime) = step
      this%valueForTime(this%nTime) = 0.0_RKIND !Clear uninitialized value
      selectcase (this%component)
      case (iBloqueJx)
         do JJJ = j1, j2
            if (k1_m - 1 >= lbound(yF, 3) .and. k1_m - 1 <= ubound(yF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + yF(i1_m, JJJ, k1_m - 1)*dy(JJJ)
            end if
            if (k2_m >= lbound(yF, 3) .and. k2_m <= ubound(yF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) - yF(i1_m, JJJ, k2_m)*dy(JJJ)
            end if
         end do
         do KKK = max(k1, lbound(zF, 3)), min(k2, ubound(zF, 3))
            this%valueForTime(this%nTime) = &
               this%valueForTime(this%nTime) + &
               (-zF(i1_m, j1_m - 1, KKK) + zF(i1_m, j2_m, KKK))*dz(KKK)
         end do

      case (iBloqueJy)
         do KKK = max(k1, lbound(zF, 3)), min(k2, ubound(zF, 3))
            this%valueForTime(this%nTime) = &
               this%valueForTime(this%nTime) + &
               (-zF(i2_m, j1_m, KKK) + zF(i1_m - 1, j1_m, KKK))*dz(KKK)
         end do
         do III = i1, i2
            if (k2_m >= lbound(xF, 3) .and. k2_m <= ubound(xF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + xF(III, j1_m, k2_m)*dx(III)
            end if
            if (k1_m - 1 >= lbound(xF, 3) .and. k1_m - 1 <= ubound(xF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) - xF(III, j1_m, k1_m - 1)*dx(III)
            end if
         end do

      case (iBloqueJz)
         if (k1_m >= lbound(xF, 3) .and. k1_m <= ubound(xF, 3)) then
            do III = i1, i2
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + &
                                               (xF(III, j1_m - 1, k1_m) - xF(III, j2_m, k1_m))*dx(III)
            end do
         end if
         if (k1_m >= lbound(yF, 3) .and. k1_m <= ubound(yF, 3)) then
            do JJJ = j1, j2
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + &
                                               (-yF(i1_m - 1, JJJ, k1_m) + yF(i2_m, JJJ, k1_m))*dy(JJJ)
            end do
         end if

      case (iBloqueMx)
         do JJJ = j1, j2
            if (k1_m >= lbound(yF, 3) .and. k1_m <= ubound(yF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) - yF(i1_m, JJJ, k1_m)*dy(JJJ)
            end if
            if (k2_m + 1 >= lbound(yF, 3) .and. k2_m + 1 <= ubound(yF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + yF(i1_m, JJJ, k2_m + 1)*dy(JJJ)
            end if
         end do
         do KKK = max(k1, lbound(zF, 3)), min(k2, ubound(zF, 3))
            this%valueForTime(this%nTime) = &
               this%valueForTime(this%nTime) + &
               (zF(i1_m, j1_m, KKK) - zF(i1_m, j2_m + 1, KKK))*dz(KKK)
         end do

      case (iBloqueMy)
         do KKK = max(k1, lbound(zF, 3)), min(k2, ubound(zF, 3))
            this%valueForTime(this%nTime) = &
               this%valueForTime(this%nTime) + &
               (zF(i2_m + 1, j1_m, KKK) - zF(i1_m, j1_m, KKK))*dz(KKK)
         end do
         do III = i1, i2
            if (k2_m + 1 >= lbound(xF, 3) .and. k2_m + 1 <= ubound(xF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) - xF(III, j1_m, k2_m + 1)*dx(III)
            end if
            if (k1_m >= lbound(xF, 3) .and. k1_m <= ubound(xF, 3)) then
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + xF(III, j1_m, k1_m)*dx(III)
            end if
         end do

      case (iBloqueMz)
         if (k1_m >= lbound(xF, 3) .and. k1_m <= ubound(xF, 3)) then
            do III = i1, i2
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + &
                                               (-xF(III, j1_m, k1_m) + xF(III, j2_m + 1, k1_m))*dx(III)
            end do
         end if
         if (k1_m >= lbound(yF, 3) .and. k1_m <= ubound(yF, 3)) then
            do JJJ = j1, j2
               this%valueForTime(this%nTime) = this%valueForTime(this%nTime) + &
                                               (yF(i1_m, JJJ, k1_m) - yF(i2_m + 1, JJJ, k1_m))*dy(JJJ)
            end do
         end if

      end select

#ifdef CompileWithMPI
      localValue = this%valueForTime(this%nTime)
      call MPI_Allreduce(localValue, this%valueForTime(this%nTime), 1, REALSIZE, MPI_SUM, MPI_COMM_WORLD, ierr)
#endif

   end subroutine update_bulk_probe_output

   subroutine flush_bulk_probe_output(this)
      type(bulk_current_probe_output_t), intent(inout) :: this
      integer :: i
      integer :: unit
      if (this%nTime <= 0) then
         print *, "No data to write."
         return
      end if

      if (this%isWriter) then
         open (newunit=unit, file=this%filePathTime, status="old", action="write", position="append")

         do i = 1, this%nTime
            write (unit, fmt) this%timeStep(i), this%valueForTime(i)
         end do

         close (unit)
      end if
      call clear_time_data()
   contains
      subroutine clear_time_data()
         this%timeStep = 0.0_RKIND_tiempo
         this%valueForTime = 0.0_RKIND

         this%nTime = 0
      end subroutine clear_time_data
   end subroutine flush_bulk_probe_output

end module bulkProbeOutput_m
