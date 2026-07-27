program test_output_transport
   use mpi
   use FDETYPES_m, only: RKIND
    use outputTransport_m, only: output_transport_t, init_output_transport, &
                                 gather_point_eligibility, reduce_scalar_batch, transfer_flush_batch, &
                                 OUTPUT_TRANSPORT_SUCCESS
    use lineProbeOutput_m, only: reduce_line_probe_sample
   implicit none

   integer, parameter :: root_rank = 0
   integer :: ierr, rank, rank_count, status, failures, i
   integer, allocatable :: counts(:), displacements(:)
   logical :: local_eligible
   logical, allocatable :: eligibility(:)
    real(kind=RKIND) :: local_scalars(2), canonical_line_value
   real(kind=RKIND), allocatable :: reduced_scalars(:), local_batch(:), gathered_batch(:)
   type(output_transport_t) :: transport

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, rank_count, ierr)

   failures = 0
   call init_output_transport(transport, root_rank, status)
   if (status /= OUTPUT_TRANSPORT_SUCCESS) failures = failures + 1

   local_eligible = rank == 0 .or. rank == rank_count - 1
   call gather_point_eligibility(transport, local_eligible, eligibility, status)
   if (status /= OUTPUT_TRANSPORT_SUCCESS .or. size(eligibility) /= rank_count) failures = failures + 1
   do i = 1, rank_count
      if (eligibility(i) .neqv. (i == 1 .or. i == rank_count)) failures = failures + 1
   end do

   local_scalars = [real(rank + 1, RKIND), real(2 * (rank + 1), RKIND)]
   call reduce_scalar_batch(transport, local_scalars, reduced_scalars, status)
   if (status /= OUTPUT_TRANSPORT_SUCCESS) failures = failures + 1
    if (rank == root_rank) then
      if (any(abs(reduced_scalars - [real(rank_count * (rank_count + 1) / 2, RKIND), &
                                     real(rank_count * (rank_count + 1), RKIND)]) > 0.0_RKIND)) then
         failures = failures + 1
    end if

    call reduce_line_probe_sample(transport, real(rank + 1, RKIND), canonical_line_value, status)
    if (status /= OUTPUT_TRANSPORT_SUCCESS) failures = failures + 1
    if (rank == root_rank .and. canonical_line_value /= real(rank_count * (rank_count + 1) / 2, RKIND)) then
       failures = failures + 1
    end if
   end if

   allocate(local_batch(rank + 1))
   do i = 1, size(local_batch)
      local_batch(i) = real(10 * rank + i, RKIND)
   end do
   call transfer_flush_batch(transport, local_batch, gathered_batch, counts, displacements, status)
   if (status /= OUTPUT_TRANSPORT_SUCCESS) failures = failures + 1
   if (rank == root_rank) then
      if (size(counts) /= rank_count .or. size(displacements) /= rank_count) failures = failures + 1
      do i = 1, rank_count
         if (counts(i) /= i) failures = failures + 1
         if (displacements(i) /= (i - 1) * i / 2) failures = failures + 1
      end do
      do i = 1, size(gathered_batch)
         if (gathered_batch(i) /= real(10 * rank_for_index(i) + local_index(i), RKIND)) failures = failures + 1
      end do
   end if

   call MPI_Allreduce(MPI_IN_PLACE, failures, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
   call MPI_Finalize(ierr)
   if (failures /= 0) error stop failures

contains

   pure integer function rank_for_index(index)
      integer, intent(in) :: index

      rank_for_index = 0
      do while ((rank_for_index + 1) * (rank_for_index + 2) / 2 < index)
         rank_for_index = rank_for_index + 1
      end do
   end function rank_for_index

   pure integer function local_index(index)
      integer, intent(in) :: index
      integer :: owner

      owner = rank_for_index(index)
      local_index = index - owner * (owner + 1) / 2
   end function local_index

end program test_output_transport
