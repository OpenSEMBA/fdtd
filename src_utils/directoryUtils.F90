module directoryUtils_m
   use FDETYPES_m
   use iso_c_binding, only: c_char, c_int, c_null_char
   implicit none
   private

   public :: add_extension
   public :: create_folder
   public :: remove_extension
   public :: folder_exists
   public :: join_path
   public :: get_last_component
   public :: remove_folder
   public :: file_exists
   public :: delete_file
   public :: list_files
   public :: create_file_with_path
   public :: get_path_separator

   abstract interface
      integer(c_int) function native_path_operation(path) bind(C)
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: path(*)
      end function native_path_operation
   end interface

   interface
      integer(c_int) function fdtd_create_directories(path) bind(C, name='fdtd_create_directories')
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: path(*)
      end function fdtd_create_directories

      integer(c_int) function fdtd_remove_tree(path) bind(C, name='fdtd_remove_tree')
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: path(*)
      end function fdtd_remove_tree

      integer(c_int) function fdtd_delete_file(path) bind(C, name='fdtd_delete_file')
         import :: c_char, c_int
         character(kind=c_char), intent(in) :: path(*)
      end function fdtd_delete_file
   end interface

contains

   !------------------------------------------------------------
   ! Add an extension to a filename
   !------------------------------------------------------------
   function add_extension(filename, ext) result(fullname)
       character(len=*), intent(in) :: filename, ext
       character(len=:), allocatable :: fullname

       fullname = trim(filename) // trim(ext)
   end function add_extension

   !------------------------------------------------------------
   ! Check if a folder exists
   !------------------------------------------------------------
   function folder_exists(path) result(exists)
      character(len=*), intent(in) :: path
      logical :: exists
      character(len=BUFSIZE) :: p

      p = trim(path)
      if (index(p, '\') > 0) then
         p = trim(p)//"\"
      else
         p = trim(p)//"/"
      end if
#ifdef GNUCompiler
      inquire (file=p, exist=exists)
#else
      inquire (directory=p, exist=exists)
#endif
   end function folder_exists

   !------------------------------------------------------------
   ! Remove the final extension from a filename
   !------------------------------------------------------------
   function remove_extension(filename) result(base)
       character(len=*), intent(in) :: filename
       character(len=BUFSIZE) :: base
       integer :: last_dot, n, i

       base = trim(filename)
       n = len_trim(base)
       last_dot = 0

       do i = n, 1, -1
           if (base(i:i) == '.') then
               last_dot = i
               exit
           end if
       end do

       if (last_dot > 0) then
           base = base(:last_dot-1)
       end if
   end function remove_extension

   !------------------------------------------------------------
   ! Join two path components into one (simplified)
   !------------------------------------------------------------
   function join_path(base, child) result(fullpath)
      character(len=*), intent(in) :: base, child
      character(len=:), allocatable :: fullpath
      character(len=1) :: sep
      integer :: n


      sep = get_path_separator()


      fullpath = trim(base)
      n = len_trim(fullpath)


      if (n > 0) then
      if (fullpath(n:n) /= sep) fullpath = fullpath // sep
      end if


      fullpath = fullpath // trim(child)
      end function join_path

   !------------------------------------------------------------
   ! Get the last component of a path (file or folder)
   !------------------------------------------------------------
   function get_last_component(path) result(component)
      character(len=*), intent(in) :: path
      character(len=BUFSIZE) :: component
      integer :: last_slash, n

      n = len_trim(path)
      component = path(:n)

      if (n > 0) then
         if (component(n:n) == get_path_separator()) component = component(:n - 1)
      end if

      last_slash = scan(component, get_path_separator(),.TRUE.)

      if (last_slash > 0) then
         component = component(last_slash + 1:)
      end if
   end function get_last_component

   !------------------------------------------------------------
   ! Create a folder (portable)
   !------------------------------------------------------------
   subroutine create_folder(path, ios)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios

      if (folder_exists(path)) then
         ios = 0
         return
      end if

       call call_native_path_operation(path, fdtd_create_directories, ios)
   end subroutine create_folder

   !------------------------------------------------------------
   ! Remove a folder
   !------------------------------------------------------------
   subroutine remove_folder(path, ios)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios

      if (.not. folder_exists(path)) then
         ios = 0
         return
      end if

       call call_native_path_operation(path, fdtd_remove_tree, ios)

   end subroutine remove_folder

   !------------------------------------------------------------
   ! Check if a file exists
   !------------------------------------------------------------
   function file_exists(path) result(exists)
      character(len=*), intent(in) :: path
      logical :: exists

      inquire (file=trim(path), exist=exists)
   end function file_exists

   !------------------------------------------------------------
   ! Delete a file
   !------------------------------------------------------------
   subroutine delete_file(path, ios)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ios

      if (.not. file_exists(path)) then
         ios = 0
         return
      end if

       call call_native_path_operation(path, fdtd_delete_file, ios)

   end subroutine delete_file

   !------------------------------------------------------------
   ! List files in a folder (simple)
   !------------------------------------------------------------
   subroutine list_files(path, files, nfiles, ios)
      character(len=*), intent(in) :: path
      character(len=BUFSIZE), dimension(:), intent(out) :: files
      integer, intent(out) :: nfiles
      integer, intent(out) :: ios

      character(len=BUFSIZE) :: cmd
      character(len=BUFSIZE) :: line
      integer :: i
      integer :: unit

      nfiles = 0
      ios = 0

      if (.not. folder_exists(path)) then
         ios = 1
         return
      end if

#ifdef __WIN32__
      cmd = 'dir /B "'//trim(path)//'"'
#else
      cmd = 'ls -1 "'//trim(path)//'"'
#endif

      open (newunit=unit, file=cmd, action='read', status='old', iostat=ios)

      if (ios /= 0) return

      do
         read (unit, '(A)', iostat=ios) line
         if (ios /= 0) exit
         i = i + 1
         if (i > size(files)) then
            ios = 2
            exit
         end if
         files(i) = adjustl(trim(line))
      end do

      nfiles = i
      close (unit)

   end subroutine list_files

   !------------------------------------------------------------
   ! Create a file, creating its folder if needed
   !------------------------------------------------------------
    subroutine create_file_with_path(fullpath, ios)
      character(len=*), intent(in) :: fullpath
      integer, intent(out) :: ios
      integer :: unit

      character(len=BUFSIZE) :: folder
      integer :: pos

      ios = 0
      ! Find last slash or backslash
       pos = scan(trim(fullpath), '/\\', back=.true.)

      if (pos > 0) then
         folder = adjustl(fullpath(:pos - 1))
         call create_folder(trim(folder), ios)
         if (ios /= 0) return
      end if
      open (newunit=unit, file=trim(adjustl(fullpath)), status='replace', iostat=ios)
      if (ios == 0) close (unit)

    end subroutine create_file_with_path

   subroutine call_native_path_operation(path, operation, ios)
      character(len=*), intent(in) :: path
      procedure(native_path_operation) :: operation
      integer, intent(out) :: ios
      character(kind=c_char), allocatable :: c_path(:)
      integer :: i, path_length

      path_length = len_trim(path)
      allocate(c_path(path_length + 1))
      do i = 1, path_length
         c_path(i) = path(i:i)
      end do
      c_path(path_length + 1) = c_null_char
      ios = operation(c_path)
   end subroutine call_native_path_operation

    function get_path_separator() result(sep)
      character(len=1) :: sep

#ifdef __WIN32__
      sep = '\'
#else
      sep = '/'
#endif

   end function get_path_separator

end module directoryUtils_m
