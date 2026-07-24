module outputTransport_m
   use FDETYPES_m, only: RKIND
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
      logical :: distributed = .false.
   end type output_transport_t

   public :: init_output_transport
   public :: gather_point_eligibility
   public :: reduce_scalar_batch
   public :: transfer_flush_batch

contains

   subroutine init_output_transport(transport, root_rank, status)
      type(output_transport_t), intent(out) :: transport
      integer, intent(in), optional :: root_rank
      integer, intent(out) :: status
#ifdef CompileWithMPI
      integer :: ierr
#endif

      transport = output_transport_t()
      if (present(root_rank)) transport%root_rank = root_rank

#ifdef CompileWithMPI
      call MPI_Comm_rank(MPI_COMM_WORLD, transport%rank, ierr)
      if (ierr /= MPI_SUCCESS) then
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
      call MPI_Comm_size(MPI_COMM_WORLD, transport%rank_count, ierr)
      if (ierr /= MPI_SUCCESS) then
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
      transport%distributed = transport%rank_count > 1
#endif

      if (.not. valid_transport(transport)) then
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if
      status = OUTPUT_TRANSPORT_SUCCESS
   end subroutine init_output_transport

   subroutine gather_point_eligibility(transport, local_eligible, rank_is_eligible, status)
      type(output_transport_t), intent(in) :: transport
      logical, intent(in) :: local_eligible
      logical, allocatable, intent(out) :: rank_is_eligible(:)
      integer, intent(out) :: status
#ifdef CompileWithMPI
      integer :: ierr
#endif

      if (allocated(rank_is_eligible)) deallocate(rank_is_eligible)
      if (.not. valid_transport(transport)) then
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if
      allocate(rank_is_eligible(transport%rank_count))

#ifdef CompileWithMPI
      call MPI_Allgather(local_eligible, 1, MPI_LOGICAL, rank_is_eligible, 1, MPI_LOGICAL, MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         deallocate(rank_is_eligible)
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
#else
      rank_is_eligible(1) = local_eligible
#endif

      status = OUTPUT_TRANSPORT_SUCCESS
   end subroutine gather_point_eligibility

   subroutine reduce_scalar_batch(transport, local_values, canonical_values, status)
      type(output_transport_t), intent(in) :: transport
      real(kind=RKIND), intent(in) :: local_values(:)
      real(kind=RKIND), allocatable, intent(out) :: canonical_values(:)
      integer, intent(out) :: status
#ifdef CompileWithMPI
      integer :: ierr
#endif

      if (allocated(canonical_values)) deallocate(canonical_values)
      if (.not. valid_transport(transport)) then
         status = OUTPUT_TRANSPORT_INVALID_CONTEXT
         return
      end if
      allocate(canonical_values(size(local_values)))
      canonical_values = 0.0_RKIND

#ifdef CompileWithMPI
      call MPI_Reduce(local_values, canonical_values, size(local_values), mpi_rkind(), MPI_SUM, transport%root_rank, &
                      MPI_COMM_WORLD, ierr)
      if (ierr /= MPI_SUCCESS) then
         deallocate(canonical_values)
         status = OUTPUT_TRANSPORT_RUNTIME_FAILURE
         return
      end if
#else
      canonical_values = local_values
#endif

      status = OUTPUT_TRANSPORT_SUCCESS
   end subroutine reduce_scalar_batch

   subroutine transfer_flush_batch(transport, local_batch, gathered_batch, counts, displacements, status)
      type(output_transport_t), intent(in) :: transport
      real(kind=RKIND), intent(in) :: local_batch(:)
      real(kind=RKIND), allocatable, intent(out) :: gathered_batch(:)
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
                      MPI_COMM_WORLD, ierr)
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
      call MPI_Gatherv(local_batch, size(local_batch), mpi_rkind(), gathered_batch, counts, displacements, mpi_rkind(), &
                       transport%root_rank, MPI_COMM_WORLD, ierr)
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
      real(kind=RKIND), allocatable, intent(inout) :: gathered_batch(:)
      integer, allocatable, intent(inout) :: counts(:), displacements(:)

      if (allocated(gathered_batch)) deallocate(gathered_batch)
      if (allocated(counts)) deallocate(counts)
      if (allocated(displacements)) deallocate(displacements)
   end subroutine clear_flush_batch

#ifdef CompileWithMPI
   integer function mpi_rkind()
#ifdef CompileWithReal8
      mpi_rkind = MPI_DOUBLE_PRECISION
#else
      mpi_rkind = MPI_REAL
#endif
   end function mpi_rkind
#endif

end module outputTransport_m
