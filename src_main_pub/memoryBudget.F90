module memoryBudget_m
   use FDETYPES_m
   use Report_m, only: log_msg, LOG_INFO, LOG_WARN, LOG_FATAL
   use systemMemory_m, only: system_memory_t, read_system_memory, system_host_name
   implicit none
   private

   public :: estimate_solver_memory
   public :: evaluate_memory_budget
   public :: report_memory_verdict
   public :: check_memory_budget
   public :: format_memory_amount
   public :: format_memory_verdict_message

   integer, parameter :: HOST_NAME_LENGTH = 256
   integer(kind=8), parameter :: UNKNOWN_MEMORY = huge(1_8)/4_8

   type, public :: memory_estimate_t
      integer(kind=8) :: fieldBytes = 0_8
      integer(kind=8) :: distanceBytes = 0_8
      integer(kind=8) :: mediaCopyBytes = 0_8
      integer(kind=8) :: tagCopyBytes = 0_8
      integer(kind=8) :: totalBytes = 0_8
   end type memory_estimate_t

   type, public :: memory_verdict_t
      logical :: checked = .false.
      logical :: fits = .true.
      integer(kind=8) :: requiredBytes = 0_8
      integer(kind=8) :: availableBytes = 0_8
      integer(kind=8) :: shortfallBytes = 0_8
      integer(kind=4) :: offenderRank = 0_4
      integer(kind=4) :: reportingRank = 0_4
      integer(kind=4) :: commSize = 1_4
      character(len=HOST_NAME_LENGTH) :: offenderNode = ''
   end type memory_verdict_t

contains

   ! Estimates the memory still to be allocated by the solver at the point
   ! where it is called.  Copies are only included when the caller is still
   ! before solver construction.
   subroutine estimate_solver_memory(sgg, media, tags, includeCopies, estimate)
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(media_matrices_t), intent(in) :: media
      type(taglist_t), intent(in) :: tags
      logical, intent(in) :: includeCopies
      type(memory_estimate_t), intent(out) :: estimate

      integer(kind=8) :: bytesPerReal
      integer(kind=8) :: bytesPerMedia
      integer(kind=8) :: bytesPerTag

      estimate = memory_estimate_t()
      bytesPerReal = int(storage_size(1.0_RKIND)/8, 8)
      bytesPerMedia = int(INTEGERSIZEOFMEDIAMATRICES, 8)
      bytesPerTag = int(IKINDMTAG, 8)

      estimate%fieldBytes = bytesPerReal * (cells_of(sgg%Alloc(IEX)) + cells_of(sgg%Alloc(IEY)) &
                                           + cells_of(sgg%Alloc(IEZ)) + cells_of(sgg%Alloc(IHX)) &
                                           + cells_of(sgg%Alloc(IHY)) + cells_of(sgg%Alloc(IHZ)))

      ! dx/Idx pairs for the electric and magnetic grids along each axis.
      estimate%distanceBytes = 2_8 * bytesPerReal * &
         (int(sgg%Alloc(IHX)%XE - sgg%Alloc(IHX)%XI + 1, 8) + &
          int(sgg%Alloc(IHY)%YE - sgg%Alloc(IHY)%YI + 1, 8) + &
          int(sgg%Alloc(IHZ)%ZE - sgg%Alloc(IHZ)%ZI + 1, 8) + &
          int(sgg%Alloc(IEX)%XE - sgg%Alloc(IEX)%XI + 1, 8) + &
          int(sgg%Alloc(IEY)%YE - sgg%Alloc(IEY)%YI + 1, 8) + &
          int(sgg%Alloc(IEZ)%ZE - sgg%Alloc(IEZ)%ZI + 1, 8))

      if (includeCopies) then
         estimate%mediaCopyBytes = media_copy_bytes(media, bytesPerMedia, bytesPerTag)
         estimate%tagCopyBytes = tag_copy_bytes(tags, bytesPerTag)
      end if

      estimate%totalBytes = estimate%fieldBytes + estimate%distanceBytes + &
                            estimate%mediaCopyBytes + estimate%tagCopyBytes
   end subroutine estimate_solver_memory

   function media_copy_bytes(media, bytesPerMedia, bytesPerTag) result(bytes)
      type(media_matrices_t), intent(in) :: media
      integer(kind=8), intent(in) :: bytesPerMedia
      integer(kind=8), intent(in) :: bytesPerTag
      integer(kind=8) :: bytes

      bytes = 0_8
      if (allocated(media%sggMiNo)) then
         bytes = bytes + bytesPerMedia * (size(media%sggMiNo, kind=8) + size(media%sggMiEx, kind=8) &
                                         + size(media%sggMiEy, kind=8) + size(media%sggMiEz, kind=8) &
                                         + size(media%sggMiHx, kind=8) + size(media%sggMiHy, kind=8) &
                                         + size(media%sggMiHz, kind=8))
      end if
      if (allocated(media%sggMtag)) bytes = bytes + bytesPerTag * size(media%sggMtag, kind=8)
   end function media_copy_bytes

   function tag_copy_bytes(tags, bytesPerTag) result(bytes)
      type(taglist_t), intent(in) :: tags
      integer(kind=8), intent(in) :: bytesPerTag
      integer(kind=8) :: bytes

      bytes = 0_8
      bytes = bytes + allocated_int3_bytes(tags%edge%x, bytesPerTag)
      bytes = bytes + allocated_int3_bytes(tags%edge%y, bytesPerTag)
      bytes = bytes + allocated_int3_bytes(tags%edge%z, bytesPerTag)
      bytes = bytes + allocated_int3_bytes(tags%face%x, bytesPerTag)
      bytes = bytes + allocated_int3_bytes(tags%face%y, bytesPerTag)
      bytes = bytes + allocated_int3_bytes(tags%face%z, bytesPerTag)
   end function tag_copy_bytes

   function allocated_int3_bytes(array, bytesPerElement) result(bytes)
      integer(kind=IKINDMTAG), allocatable, intent(in) :: array(:, :, :)
      integer(kind=8), intent(in) :: bytesPerElement
      integer(kind=8) :: bytes

      if (allocated(array)) then
         bytes = bytesPerElement * size(array, kind=8)
      else
         bytes = 0_8
      end if
   end function allocated_int3_bytes

   pure integer(kind=8) function cells_of(lim) result(elements)
      type(XYZlimit_t), intent(in) :: lim
      elements = int(lim%XE - lim%XI + 1, 8) * int(lim%YE - lim%YI + 1, 8) * int(lim%ZE - lim%ZI + 1, 8)
   end function cells_of

   ! Collective check of the memory still needed by this rank.
   subroutine evaluate_memory_budget(pendingBytes, status, comm, verdict)
      integer(kind=8), intent(in) :: pendingBytes
      type(system_memory_t), intent(in) :: status
      integer(kind=4), intent(in) :: comm
      type(memory_verdict_t), intent(out) :: verdict

#ifdef CompileWithMPI
      integer(kind=4) :: ierr
      integer(kind=4) :: commRank
      integer(kind=4) :: commSize
      integer(kind=4) :: nodeComm
      integer(kind=4) :: nodeRank
      integer(kind=4) :: nodeSize
      integer(kind=4) :: candidate
      integer(kind=8) :: localAvailable
      integer(kind=8) :: nodeAvailable
      integer(kind=8) :: localHeadroom
      integer(kind=8) :: localShortfall
      integer(kind=8) :: maxShortfall
      integer(kind=8) :: shareBytes(2)
      character(len=HOST_NAME_LENGTH) :: hostName
      logical :: splitSucceeded
      logical :: localKnown, checkKnown

      verdict = memory_verdict_t()
      call MPI_Comm_rank(comm, commRank, ierr)
      call MPI_Comm_size(comm, commSize, ierr)
      verdict%reportingRank = commRank
      verdict%commSize = commSize

      splitSucceeded = .true.
      call MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0_4, MPI_INFO_NULL, nodeComm, ierr)
      if (ierr /= MPI_SUCCESS) then
         splitSucceeded = .false.
         nodeComm = comm
      end if
      call MPI_Comm_rank(nodeComm, nodeRank, ierr)
      call MPI_Comm_size(nodeComm, nodeSize, ierr)

      if (status%availableKnown) then
         localAvailable = status%availableBytes
      else
         localAvailable = UNKNOWN_MEMORY
      end if
      call MPI_Allreduce(localAvailable, nodeAvailable, 1_4, MPI_INTEGER8, MPI_MIN, nodeComm, ierr)

      localKnown = status%availableKnown
      call MPI_Allreduce(localKnown, checkKnown, 1_4, MPI_LOGICAL, MPI_LAND, comm, ierr)
      if (.not. checkKnown) then
         verdict%checked = .false.
         verdict%fits = .true.
         if (splitSucceeded) call MPI_Comm_free(nodeComm, ierr)
         return
      end if

      localHeadroom = 0_8
      if (nodeAvailable /= UNKNOWN_MEMORY) then
         localHeadroom = nodeAvailable / max(1_4, nodeSize)
      end if
      if (status%limitKnown .and. status%processLimitBytes > 0_8) then
         localHeadroom = min(localHeadroom, max(0_8, status%processLimitBytes - process_used_bytes(status)))
      end if

      verdict%checked = nodeAvailable /= UNKNOWN_MEMORY
      verdict%requiredBytes = pendingBytes
      verdict%availableBytes = localHeadroom
      verdict%shortfallBytes = max(0_8, pendingBytes - localHeadroom)
      verdict%offenderRank = commRank

      call MPI_Allreduce(verdict%shortfallBytes, maxShortfall, 1_4, MPI_INTEGER8, MPI_MAX, comm, ierr)
      if (maxShortfall <= 0_8) then
         verdict%fits = .true.
         verdict%shortfallBytes = 0_8
      else
         verdict%fits = .false.
         candidate = merge(commRank, -1_4, verdict%shortfallBytes == maxShortfall)
         call MPI_Allreduce(candidate, verdict%offenderRank, 1_4, MPI_INTEGER, MPI_MAX, comm, ierr)

         hostName = system_host_name()
         call MPI_Bcast(hostName, HOST_NAME_LENGTH, MPI_CHARACTER, verdict%offenderRank, comm, ierr)
         verdict%offenderNode = hostName

         shareBytes = [verdict%requiredBytes, verdict%availableBytes]
         call MPI_Bcast(shareBytes, 2_4, MPI_INTEGER8, verdict%offenderRank, comm, ierr)
         verdict%requiredBytes = shareBytes(1)
         verdict%availableBytes = shareBytes(2)
         verdict%shortfallBytes = max(0_8, verdict%requiredBytes - verdict%availableBytes)
      end if

      if (splitSucceeded) then
         call MPI_Comm_free(nodeComm, ierr)
      end if
#else
      integer(kind=8) :: headroom

      verdict = memory_verdict_t()
      verdict%commSize = 1_4
      if (status%availableKnown) then
         headroom = status%availableBytes
         if (status%limitKnown .and. status%processLimitBytes > 0_8) then
            headroom = min(headroom, max(0_8, status%processLimitBytes - process_used_bytes(status)))
         end if
         verdict%checked = .true.
         verdict%requiredBytes = pendingBytes
         verdict%availableBytes = headroom
         verdict%offenderRank = 0_4
         verdict%offenderNode = system_host_name()
         verdict%fits = pendingBytes <= headroom
         verdict%shortfallBytes = max(0_8, pendingBytes - headroom)
      end if
#endif
   end subroutine evaluate_memory_budget

   integer(kind=8) function process_used_bytes(status) result(bytes)
      type(system_memory_t), intent(in) :: status

      if (status%virtualKnown) then
         bytes = status%virtualBytes
      else if (status%residentKnown) then
         bytes = status%residentBytes
      else
         bytes = 0_8
      end if
   end function process_used_bytes

   subroutine check_memory_budget(pendingBytes, context)
      integer(kind=8), intent(in) :: pendingBytes
      character(len=*), intent(in) :: context

      type(system_memory_t) :: status
      type(memory_verdict_t) :: verdict
      integer(kind=4) :: comm

      call read_system_memory(status)
#ifdef CompileWithMPI
      comm = SUBCOMM_MPI
      if (comm == 0_4) comm = MPI_COMM_WORLD
#else
      comm = 0_4
#endif
      call evaluate_memory_budget(pendingBytes, status, comm, verdict)
      call report_memory_verdict(verdict, context)
   end subroutine check_memory_budget

   subroutine report_memory_verdict(verdict, context)
      type(memory_verdict_t), intent(in) :: verdict
      character(len=*), intent(in) :: context

      character(len=BUFSIZE) :: message

      if (.not. verdict%checked) then
#ifdef CompileWithMPI
         if (verdict%reportingRank == 0_4) then
#endif
            call log_msg(LOG_WARN, 'MEMORY', &
                         'memory check skipped ('//trim(context)//'): available memory could not be determined')
#ifdef CompileWithMPI
         end if
#endif
         return
      end if

      if (verdict%fits) then
#ifdef CompileWithMPI
         if (verdict%reportingRank == 0_4) then
#endif
            message = 'memory check passed ('//trim(context)//'): required '// &
                      format_memory_amount(verdict%requiredBytes)//', available '// &
                      format_memory_amount(verdict%availableBytes)
            call log_msg(LOG_INFO, 'MEMORY', trim(message))
#ifdef CompileWithMPI
         end if
#endif
         return
      end if

#ifdef CompileWithMPI
      if (verdict%reportingRank /= 0_4) return
#endif
      message = format_memory_verdict_message(verdict, context)
      call log_msg(LOG_FATAL, 'MEMORY', trim(message))
   end subroutine report_memory_verdict

   function format_memory_verdict_message(verdict, context) result(message)
      type(memory_verdict_t), intent(in) :: verdict
      character(len=*), intent(in) :: context
      character(len=BUFSIZE) :: message
      character(len=16) :: rankText
      character(len=64) :: suffix

      if (verdict%commSize <= 1_4) then
         message = 'Not enough memory to run the simulation.'//achar(10)// &
                   '  Required:  '//trim(format_memory_amount(verdict%requiredBytes))//achar(10)// &
                   '  Available:  '//trim(format_memory_amount(verdict%availableBytes))
      else
         write(rankText, '(i0)') verdict%offenderRank
         suffix = trim(context)
         if (len_trim(suffix) > 0) suffix = ' ['//trim(suffix)//']'
         message = 'Not enough memory on rank '//trim(rankText)//' ('// &
                   trim(verdict%offenderNode)//').'//trim(suffix)//achar(10)// &
                   '  Required by process:  '//trim(format_memory_amount(verdict%requiredBytes))//achar(10)// &
                   '  Available to process:  '//trim(format_memory_amount(verdict%availableBytes))
      end if
   end function format_memory_verdict_message

   function format_memory_amount(bytes) result(text)
      integer(kind=8), intent(in) :: bytes
      character(len=:), allocatable :: text
      character(len=32) :: buffer
      real(kind=8) :: value

      value = real(bytes, 8)
      if (abs(value) >= 1073741824.0_8) then
         write(buffer, '(f12.1,1x,a)') value/1073741824.0_8, 'GB'
      else if (abs(value) >= 1048576.0_8) then
         write(buffer, '(f12.1,1x,a)') value/1048576.0_8, 'MB'
      else
         write(buffer, '(f12.1,1x,a)') value/1024.0_8, 'KB'
      end if
      text = trim(adjustl(buffer))
   end function format_memory_amount

end module memoryBudget_m
