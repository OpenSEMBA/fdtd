module outputTransport_m
   use, intrinsic :: iso_fortran_env, only: int64, real64
#ifdef CompileWithMPI
   use mpi
#endif
   implicit none
   private

   integer, parameter, public :: OUTPUT_TRANSPORT_SUCCESS = 0
   integer, parameter, public :: OUTPUT_TRANSPORT_INVALID_CONTEXT = 1
   integer, parameter, public :: OUTPUT_TRANSPORT_RUNTIME_FAILURE = 2
   integer, parameter, public :: OUTPUT_TRANSPORT_SIZE_OVERFLOW = 3
   integer, parameter :: RELAY_TAG = 7301

   type, public :: output_transport_t
      integer :: rank = 0
      integer :: rank_count = 1
      integer :: root_rank = 0
      integer :: communicator = 0
      integer :: real64_datatype = 0
   end type output_transport_t

   public :: init_output_transport
   public :: transfer_flush_batch
   public :: relay_flush_batch

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
      call MPI_Type_match_size(MPI_TYPECLASS_REAL, storage_size(0.0_real64)/8, &
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

      if (allocated(gathered_batch)) deallocate (gathered_batch)
      if (allocated(counts)) deallocate (counts)
      if (allocated(displacements)) deallocate (displacements)
      if (.not. valid_transport(transport)) then
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if

      allocate (counts(transport%rank_count), displacements(transport%rank_count))
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

#ifdef CompileWithMPI
      if (transport%rank == transport%root_rank) then
          ! MPI_Gatherv counts and displacements are default integers.
          if (sum(int(counts, int64)) > int(huge(0), int64)) then
             status = OUTPUT_TRANSPORT_SIZE_OVERFLOW
          else
             status = OUTPUT_TRANSPORT_SUCCESS
          end if
      end if
      call MPI_Bcast(status, 1, MPI_INTEGER, transport%root_rank, transport%communicator, ierr)
      if (ierr /= MPI_SUCCESS) then
          call clear_flush_batch(gathered_batch, counts, displacements)
          status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
          return
      end if
      if (status /= OUTPUT_TRANSPORT_SUCCESS) then
          call clear_flush_batch(gathered_batch, counts, displacements)
          return
      end if
#endif

      if (transport%rank == transport%root_rank) then
          do i = 2, transport%rank_count
             displacements(i) = displacements(i - 1) + counts(i - 1)
          end do
         total_count = sum(counts)
         allocate (gathered_batch(total_count))
      else
         allocate (gathered_batch(0))
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

   ! Moves the batch of source_rank to the root, so the root holds one rank's batch at a time.
   ! Every rank must call it with the same source_rank.
   subroutine relay_flush_batch(transport, source_rank, local_batch, batch, status)
      type(output_transport_t), intent(in) :: transport
      integer, intent(in) :: source_rank
      real(real64), intent(in) :: local_batch(:)
      real(real64), allocatable, intent(out) :: batch(:)
      integer, intent(out) :: status
#ifdef CompileWithMPI
      integer :: batch_size, ierr
#endif

      if (.not. valid_transport(transport) .or. source_rank < 0 .or. source_rank >= transport%rank_count) then
         allocate (batch(0))
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if
      status = OUTPUT_TRANSPORT_SUCCESS
      if (source_rank == transport%root_rank) then
         if (transport%rank == transport%root_rank) then
            batch = local_batch
         else
            allocate (batch(0))
         end if
         return
      end if

#ifdef CompileWithMPI
      ierr = MPI_SUCCESS
      if (transport%rank == source_rank) then
         allocate (batch(0))
         call MPI_Send(size(local_batch), 1, MPI_INTEGER, transport%root_rank, RELAY_TAG, &
                       transport%communicator, ierr)
         if (ierr == MPI_SUCCESS) call MPI_Send(local_batch, size(local_batch), transport%real64_datatype, &
                                                transport%root_rank, RELAY_TAG, transport%communicator, ierr)
      else if (transport%rank == transport%root_rank) then
         call MPI_Recv(batch_size, 1, MPI_INTEGER, source_rank, RELAY_TAG, transport%communicator, &
                       MPI_STATUS_IGNORE, ierr)
         if (ierr /= MPI_SUCCESS) batch_size = 0
         allocate (batch(batch_size))
         if (ierr == MPI_SUCCESS) call MPI_Recv(batch, batch_size, transport%real64_datatype, source_rank, &
                                                RELAY_TAG, transport%communicator, MPI_STATUS_IGNORE, ierr)
      else
         allocate (batch(0))
      end if
      if (ierr /= MPI_SUCCESS) status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
#else
      allocate (batch(0))
#endif
   end subroutine relay_flush_batch

   pure logical function valid_transport(transport)
      type(output_transport_t), intent(in) :: transport

      valid_transport = transport%rank_count >= 1 .and. transport%rank >= 0 .and. &
                        transport%rank < transport%rank_count .and. transport%root_rank >= 0 .and. &
                        transport%root_rank < transport%rank_count
   end function valid_transport

   subroutine clear_flush_batch(gathered_batch, counts, displacements)
      real(real64), allocatable, intent(inout) :: gathered_batch(:)
      integer, allocatable, intent(inout) :: counts(:), displacements(:)

      if (allocated(gathered_batch)) deallocate (gathered_batch)
      if (allocated(counts)) deallocate (counts)
      if (allocated(displacements)) deallocate (displacements)
   end subroutine clear_flush_batch

end module outputTransport_m
