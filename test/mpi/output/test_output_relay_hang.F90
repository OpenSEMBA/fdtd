program test_output_relay_hang
   ! Manual reproducer: rank 1 stalls without sending, so root blocks in relay_flush_batch.
   ! Never returns by design. Not registered in CTest.
   use, intrinsic :: iso_fortran_env, only: real64
   use mpi
   use outputTransport_m, only: output_transport_t, init_output_transport, relay_flush_batch
   implicit none

   integer :: ierr, rank, status
   real(real64), allocatable :: relayed_batch(:)
   real(real64) :: local_batch(2)
   type(output_transport_t) :: transport

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

   call init_output_transport(transport, 0, status)
   local_batch = real(rank, real64)

   if (rank == 1) then
      write (*, '(a)') 'rank 1: simulating a stalled process (never sends)'
      flush (6)
      call execute_command_line('sleep 600')
   else if (rank == 0) then
      write (*, '(a)') 'rank 0: waiting on relay from rank 1 (should block)'
      flush (6)
      call relay_flush_batch(transport, 1, local_batch, relayed_batch, status)
      write (*, '(a,i0)') 'rank 0: relay returned, status=', status
   end if

   call MPI_Finalize(ierr)
end program test_output_relay_hang
