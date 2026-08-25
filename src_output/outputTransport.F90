module outputTransport_m
   use, intrinsic :: iso_fortran_env, only: real64
#ifdef CompileWithMPI
   use mpi
#endif
   implicit none
   private

   integer, parameter, public :: OUTPUT_TRANSPORT_SUCCESS = 0
   integer, parameter, public :: OUTPUT_TRANSPORT_INVALID_CONTEXT = 1
   integer, parameter, public :: OUTPUT_TRANSPORT_RUNTIME_FAILURE = 2

   type, public :: output_transport_t
      integer :: rank = 0
      integer :: rank_count = 1
      integer :: root_rank = 0
      integer :: communicator = 0
      integer :: real64_datatype = 0
   end type output_transport_t

   public :: init_output_transport
   public :: transfer_flush_batch

contains

   subroutine init_output_transport(transport, root_rank, status, communicator)
      type(output_transport_t), intent(out) :: transport
      integer, intent(in), optional :: root_rank
      integer, intent(out) :: status
      integer, intent(in), optional :: communicator
#ifdef CompileWithMPI
      integer :: ierr
#endif

      transport = output_transport_t()
      if (present(root_rank)) transport%root_rank = root_rank

#ifdef CompileWithMPI
      transport%communicator = MPI_COMM_WORLD
      if (present(communicator)) transport%communicator = communicator
      call MPI_Type_match_size(MPI_TYPECLASS_REAL, storage_size(0.0_real64) / 8, &
                               transport%real64_datatype, ierr)
      if (ierr /= MPI_SUCCESS) then
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
      call MPI_Comm_rank(transport%communicator, transport%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
      call MPI_Comm_size(transport%communicator, transport%rank_count, ierr)
      if (ierr /= MPI_SUCCESS) then
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
#endif

      if (.not. valid_transport(transport)) then
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if
      status = OUTPUT_TRANSPORT_SUCCESS
   end subroutine init_output_transport

   subroutine transfer_flush_batch(transport, local_batch, gathered_batch, counts, displacements, status)
      type(output_transport_t), intent(in) :: transport
      real(real64), intent(in) :: local_batch(:)
      real(real64), allocatable, intent(out) :: gathered_batch(:)
      integer, allocatable, intent(out) :: counts(:), displacements(:)
      integer, intent(out) :: status
      integer :: i, total_count
#ifdef CompileWithMPI
      integer :: ierr
#endif

      if (allocated(gathered_batch)) deallocate(gathered_batch)
      if (allocated(counts)) deallocate(counts)
      if (allocated(displacements)) deallocate(displacements)
      if (.not. valid_transport(transport)) then
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if

      allocate(counts(transport%rank_count), displacements(transport%rank_count))
      counts = 0
      displacements = 0

#ifdef CompileWithMPI
      call MPI_Gather(size(local_batch), 1, MPI_INTEGER, counts, 1, MPI_INTEGER, transport%root_rank, &
                      transport%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         call clear_flush_batch(gathered_batch, counts, displacements)
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
#else
      counts(1) = size(local_batch)
#endif

      if (transport%rank == transport%root_rank) then
         do i = 2, transport%rank_count
            displacements(i) = displacements(i - 1) + counts(i - 1)
         end do
         total_count = sum(counts)
         allocate(gathered_batch(total_count))
      else
         allocate(gathered_batch(0))
      end if

#ifdef CompileWithMPI
      call MPI_Gatherv(local_batch, size(local_batch), transport%real64_datatype, gathered_batch, counts, displacements, &
                       transport%real64_datatype, transport%root_rank, transport%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
         call clear_flush_batch(gathered_batch, counts, displacements)
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
#else
      gathered_batch = local_batch
#endif

      status = OUTPUT_TRANSPORT_SUCCESS
   end subroutine transfer_flush_batch

   pure logical function valid_transport(transport)
      type(output_transport_t), intent(in) :: transport

      valid_transport = transport%rank_count >= 1 .and. transport%rank >= 0 .and. &
                        transport%rank < transport%rank_count .and. transport%root_rank >= 0 .and. &
                        transport%root_rank < transport%rank_count
   end function valid_transport

   subroutine clear_flush_batch(gathered_batch, counts, displacements)
      real(real64), allocatable, intent(inout) :: gathered_batch(:)
      integer, allocatable, intent(inout) :: counts(:), displacements(:)

      if (allocated(gathered_batch)) deallocate(gathered_batch)
      if (allocated(counts)) deallocate(counts)
      if (allocated(displacements)) deallocate(displacements)
   end subroutine clear_flush_batch

end module outputTransport_m
