module systemMemory_m
   use, intrinsic :: iso_c_binding, only: c_char, c_int, c_int64_t, c_null_char, c_size_t
   implicit none
   private

   public :: read_system_memory
   public :: system_host_name

   integer, parameter :: HOST_NAME_LENGTH = 256

   ! Snapshot of the memory available to this process.  A negative C result is
   ! translated into a zero value plus a "known" flag set to false.
   type, public :: system_memory_t
      integer(kind=8) :: totalBytes = 0_8
      integer(kind=8) :: availableBytes = 0_8
      integer(kind=8) :: residentBytes = 0_8
      integer(kind=8) :: virtualBytes = 0_8
      integer(kind=8) :: processLimitBytes = 0_8
      logical :: totalKnown = .false.
      logical :: availableKnown = .false.
      logical :: residentKnown = .false.
      logical :: virtualKnown = .false.
      logical :: limitKnown = .false.
   end type system_memory_t

   interface
      function fdtd_memory_total_bytes() bind(C, name='fdtd_memory_total_bytes') result(bytes)
         import :: c_int64_t
         integer(c_int64_t) :: bytes
      end function fdtd_memory_total_bytes

      function fdtd_memory_available_bytes() bind(C, name='fdtd_memory_available_bytes') result(bytes)
         import :: c_int64_t
         integer(c_int64_t) :: bytes
      end function fdtd_memory_available_bytes

      function fdtd_process_resident_bytes() bind(C, name='fdtd_process_resident_bytes') result(bytes)
         import :: c_int64_t
         integer(c_int64_t) :: bytes
      end function fdtd_process_resident_bytes

      function fdtd_process_virtual_bytes() bind(C, name='fdtd_process_virtual_bytes') result(bytes)
         import :: c_int64_t
         integer(c_int64_t) :: bytes
      end function fdtd_process_virtual_bytes

      function fdtd_process_memory_limit_bytes() bind(C, name='fdtd_process_memory_limit_bytes') result(bytes)
         import :: c_int64_t
         integer(c_int64_t) :: bytes
      end function fdtd_process_memory_limit_bytes

      function fdtd_host_name(buffer, size) bind(C, name='fdtd_host_name') result(status)
         import :: c_char, c_int, c_size_t
         character(kind=c_char), dimension(*) :: buffer
         integer(c_size_t), value :: size
         integer(c_int) :: status
      end function fdtd_host_name
   end interface

contains

   ! Reads a consistent snapshot of the OS memory state.
   subroutine read_system_memory(status)
      type(system_memory_t), intent(out) :: status
      integer(c_int64_t) :: bytes

      status = system_memory_t()

      bytes = fdtd_memory_total_bytes()
      if (bytes >= 0_c_int64_t) then
         status%totalBytes = int(bytes, 8)
         status%totalKnown = .true.
      end if

      bytes = fdtd_memory_available_bytes()
      if (bytes >= 0_c_int64_t) then
         status%availableBytes = int(bytes, 8)
         status%availableKnown = .true.
      end if

      bytes = fdtd_process_resident_bytes()
      if (bytes >= 0_c_int64_t) then
         status%residentBytes = int(bytes, 8)
         status%residentKnown = .true.
      end if

      bytes = fdtd_process_virtual_bytes()
      if (bytes >= 0_c_int64_t) then
         status%virtualBytes = int(bytes, 8)
         status%virtualKnown = .true.
      end if

      bytes = fdtd_process_memory_limit_bytes()
      if (bytes >= 0_c_int64_t) then
         status%processLimitBytes = int(bytes, 8)
         status%limitKnown = .true.
      end if
   end subroutine read_system_memory

   ! Returns the machine name used in the fatal message.
   function system_host_name() result(name)
      character(len=HOST_NAME_LENGTH) :: name
      character(kind=c_char), dimension(HOST_NAME_LENGTH) :: buffer
      integer :: i
      integer :: status

      buffer = c_null_char
      status = fdtd_host_name(buffer, int(HOST_NAME_LENGTH, c_size_t))
      if (status /= 0) then
         name = 'unknown'
         return
      end if

      name = ' '
      do i = 1, HOST_NAME_LENGTH
         if (buffer(i) == c_null_char) exit
         name(i:i) = buffer(i)
      end do
      name = trim(adjustl(name))
      if (len_trim(name) == 0) name = 'unknown'
   end function system_host_name

end module systemMemory_m
