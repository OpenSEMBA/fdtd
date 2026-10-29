program test_sync_timeout
   ! Mode "ok": a normal exchange and a timed barrier must complete.
   ! Mode "stall": rank 1 never answers, so rank 0 must abort after SEMBA_FDTD_SYNC_TIMEOUT_S.
   use, intrinsic :: iso_fortran_env, only: real64
   use mpi
   use MPIcomm_m, only: barrier_timed, wait_requests_timed
   use Report_m, only: OnPrint
   implicit none

   integer :: ierr, rank
   integer, dimension(1) :: req
   real(real64) :: buf
   character(len=16) :: mode

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call OnPrint
   call get_command_argument(1, mode)
   buf = 0.0_real64

   select case (trim(mode))
   case ('ok')
      if (rank == 0) call MPI_Irecv(buf, 1, MPI_DOUBLE_PRECISION, 1, 7, MPI_COMM_WORLD, req(1), ierr)
      if (rank == 1) then
         buf = 1.0_real64
         call MPI_Isend(buf, 1, MPI_DOUBLE_PRECISION, 0, 7, MPI_COMM_WORLD, req(1), ierr)
      end if
      call wait_requests_timed(req)
      call barrier_timed(MPI_COMM_WORLD)
      if (rank == 0 .and. buf /= 1.0_real64) error stop 1
   case ('stall')
      if (rank == 0) call MPI_Irecv(buf, 1, MPI_DOUBLE_PRECISION, 1, 7, MPI_COMM_WORLD, req(1), ierr)
      if (rank == 1) call execute_command_line('sleep 120')
      if (rank == 0) call wait_requests_timed(req)
   case default
      error stop 2
   end select

   call MPI_Finalize(ierr)
end program test_sync_timeout
