program test_memory_budget_mpi
   use mpi
   use FDETYPES_m, only: RKIND
   use memoryBudget_m, only: evaluate_memory_budget, memory_verdict_t
   use systemMemory_m, only: system_memory_t
   implicit none

   type(system_memory_t) :: status
   type(memory_verdict_t) :: verdict
   integer :: rank, size, ierr
   integer(kind=8) :: pending

   call MPI_Init(ierr)
   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
   call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

   if (size < 2) then
      if (rank == 0) print *, 'test_memory_budget_mpi requires at least 2 ranks'
      call MPI_Finalize(ierr)
      stop 1
   end if

   status = system_memory_t()
   status%availableBytes = 1000000_8
   status%availableKnown = .true.
   status%residentBytes = 1000_8
   status%residentKnown = .true.
   status%virtualBytes = 1000_8
   status%virtualKnown = .true.
   status%processLimitBytes = 0_8
   status%limitKnown = .true.

   pending = 100_8
   if (rank == 1) pending = 2000000_8
   call evaluate_memory_budget(pending, status, MPI_COMM_WORLD, verdict)
   if (verdict%fits) then
      if (rank == 0) print *, 'test_memory_budget_mpi FAILED: expected insufficient memory'
      call MPI_Finalize(ierr)
      stop 1
   end if
   if (verdict%offenderRank /= 1) then
      if (rank == 0) print *, 'test_memory_budget_mpi FAILED: wrong offender rank', verdict%offenderRank
      call MPI_Finalize(ierr)
      stop 1
   end if

   pending = 100_8
   call evaluate_memory_budget(pending, status, MPI_COMM_WORLD, verdict)
   if (.not. verdict%fits) then
      if (rank == 0) print *, 'test_memory_budget_mpi FAILED: expected a fitting run'
      call MPI_Finalize(ierr)
      stop 1
   end if

   call MPI_Finalize(ierr)
end program test_memory_budget_mpi
