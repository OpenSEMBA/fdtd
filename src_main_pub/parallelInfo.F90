module parallelInfo_m
   use Report_m, only: print11
   use systemMemory_m, only: system_host_name
#ifdef CompileWithMPI
   use mpi
#endif
#if defined(_OPENMP)
   use omp_lib, only: omp_get_max_threads
#endif
   implicit none
   private

   public :: collect_parallel_info
   public :: report_parallel_info
   public :: parallel_info_text

   type, public :: parallel_info_t
      integer(kind=4) :: mpiProcesses = 1_4
      integer(kind=4) :: mpiNodes = 1_4
      integer(kind=4) :: ranksPerNodeMin = 1_4
      integer(kind=4) :: ranksPerNodeMax = 1_4
      integer(kind=4) :: threadsMin = 1_4
      integer(kind=4) :: threadsMax = 1_4
      logical :: mpiEnabled = .false.
      logical :: openmpEnabled = .false.
      character(len=256) :: rootHost = ''
      character(len=64) :: ompNumThreads = ''
   end type parallel_info_t

contains

   ! Collective on MPI_COMM_WORLD when MPI is compiled in.
   subroutine collect_parallel_info(info)
      type(parallel_info_t), intent(out) :: info
      integer(kind=4) :: ierr
      integer(kind=4) :: mpiSize
      integer(kind=4) :: mpiRank
      integer(kind=4) :: nodeComm
      integer(kind=4) :: nodeRank
      integer(kind=4) :: nodeSize
      integer(kind=4) :: nodeFlag
      integer(kind=4) :: threads
      integer(kind=4) :: nameLength
      character(len=256) :: processorName
      logical :: nodeCommCreated

      info = parallel_info_t()
      info%mpiEnabled = .false.
      info%openmpEnabled = .false.

#ifdef CompileWithMPI
      info%mpiEnabled = .true.
      call MPI_Comm_size(MPI_COMM_WORLD, mpiSize, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpiRank, ierr)
      info%mpiProcesses = mpiSize

      nodeCommCreated = .true.
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0_4, MPI_INFO_NULL, nodeComm, ierr)
      if (ierr /= MPI_SUCCESS) then
         nodeCommCreated = .false.
         nodeComm = MPI_COMM_WORLD
      end if
      call MPI_Comm_rank(nodeComm, nodeRank, ierr)
      call MPI_Comm_size(nodeComm, nodeSize, ierr)

      nodeFlag = merge(1_4, 0_4, nodeRank == 0_4)
      call MPI_Allreduce(nodeFlag, info%mpiNodes, 1_4, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(nodeSize, info%ranksPerNodeMin, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(nodeSize, info%ranksPerNodeMax, 1_4, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

      if (mpiRank == 0_4) then
         processorName = ''
         nameLength = 0_4
         call MPI_Get_processor_name(processorName, nameLength, ierr)
         if (ierr == MPI_SUCCESS .and. nameLength > 0_4) then
            info%rootHost = processorName(1:nameLength)
         else
            info%rootHost = system_host_name()
         end if
      end if

      if (nodeCommCreated) call MPI_Comm_free(nodeComm, ierr)
#else
      info%mpiProcesses = 1_4
      info%mpiNodes = 1_4
      info%ranksPerNodeMin = 1_4
      info%ranksPerNodeMax = 1_4
      info%rootHost = system_host_name()
#endif

      threads = 1_4
#if defined(_OPENMP)
      info%openmpEnabled = .true.
      threads = int(omp_get_max_threads(), 4)
#endif
#ifdef CompileWithMPI
      call MPI_Allreduce(threads, info%threadsMin, 1_4, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_Allreduce(threads, info%threadsMax, 1_4, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
#else
      info%threadsMin = threads
      info%threadsMax = threads
#endif

      info%ompNumThreads = ''
      call get_environment_variable('OMP_NUM_THREADS', info%ompNumThreads)
      info%ompNumThreads = trim(adjustl(info%ompNumThreads))
   end subroutine collect_parallel_info

   subroutine report_parallel_info(info, layoutnumber)
      type(parallel_info_t), intent(in) :: info
      integer(kind=4), intent(in) :: layoutnumber

      if (layoutnumber /= 0_4) return
      call print11(layoutnumber, '&'//trim(parallel_info_text(info)))
   end subroutine report_parallel_info

   function parallel_info_text(info) result(text)
      type(parallel_info_t), intent(in) :: info
      character(len=1024) :: text
      character(len=256) :: line

      text = 'Parallel configuration:'

      if (info%mpiEnabled) then
         write(line, '(a,i0)') 'MPI processes: ', info%mpiProcesses
         call append_line(text, trim(line))
         if (info%ranksPerNodeMin == info%ranksPerNodeMax) then
            write(line, '(a,i0,a,i0,a)') 'MPI nodes: ', info%mpiNodes, ' (', &
                                        info%ranksPerNodeMin, ' ranks per node)'
         else
            write(line, '(a,i0,a,i0,a,i0,a)') 'MPI nodes: ', info%mpiNodes, ' (', &
                                             info%ranksPerNodeMin, '-', info%ranksPerNodeMax, ' ranks per node)'
         end if
         call append_line(text, trim(line))
         if (len_trim(info%rootHost) > 0) then
            call append_line(text, 'MPI rank 0 host: '//trim(info%rootHost))
         end if
      else
         call append_line(text, 'MPI processes: 1 (serial)')
         call append_line(text, 'MPI nodes: 1')
      end if

      if (.not. info%openmpEnabled) then
         call append_line(text, 'OpenMP threads per process: 1 (OpenMP not enabled)')
      else if (info%threadsMin == info%threadsMax) then
         write(line, '(a,i0)') 'OpenMP threads per process: ', info%threadsMax
         if (len_trim(info%ompNumThreads) > 0) line = trim(line)//' (OMP_NUM_THREADS='// &
                                                          trim(info%ompNumThreads)//')'
         call append_line(text, trim(line))
      else
         write(line, '(a,i0,a,i0)') 'OpenMP threads per process: ', info%threadsMin, &
                                    '-', info%threadsMax
         call append_line(text, trim(line))
      end if

      text = trim(text)
   end function parallel_info_text

   subroutine append_line(text, line)
      character(len=*), intent(inout) :: text
      character(len=*), intent(in) :: line
      text = trim(text)//achar(10)//trim(line)
   end subroutine append_line

end module parallelInfo_m
