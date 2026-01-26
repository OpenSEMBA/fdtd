module mod_directoryUtils
   implicit none
   private

   public :: create_folder
   public :: folder_exists
   public :: remove_folder
   public :: file_exists
   public :: delete_file
   public :: list_files
   public :: create_file_with_path
   public :: get_path_separator

contains

   !------------------------------------------------------------
   ! Check if a folder exists
   !------------------------------------------------------------
   function folder_exists(path) result(exists)
      character(len=*), intent(in) :: path
      logical :: exists
      character(len=256) :: p

      p = trim(path)
      if (index(p, '\') > 0) then
         p = trim(p)//"\"
      else
         p = trim(p)//"/"
      end if

      inquire (file=p, exist=exists)
   end function folder_exists

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

#ifdef _WIN32
      call execute_command_line("mkdir """//trim(path)//"""", exitstat=ios)
#else
      call execute_command_line("mkdir -p "//trim(path), exitstat=ios)
#endif
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

#ifdef _WIN32
      call execute_command_line("rmdir /S /Q """//trim(path)//"""", exitstat=ios)
#else
      call execute_command_line("rm -rf "//trim(path), exitstat=ios)
#endif

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

#ifdef _WIN32
      call execute_command_line("del /Q """//trim(path)//"""", exitstat=ios)
#else
      call execute_command_line("rm -f "//trim(path), exitstat=ios)
#endif

   end subroutine delete_file

   !------------------------------------------------------------
   ! List files in a folder (simple)
   !------------------------------------------------------------
   subroutine list_files(path, files, nfiles, ios)
      character(len=*), intent(in) :: path
      character(len=256), dimension(:), intent(out) :: files
      integer, intent(out) :: nfiles
      integer, intent(out) :: ios

      character(len=512) :: cmd
      character(len=512) :: line
      integer :: i
      integer :: unit

      nfiles = 0
      ios = 0

      if (.not. folder_exists(path)) then
         ios = 1
         return
      end if

#ifdef _WIN32
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

      character(len=512) :: folder
      integer :: pos

      ios = 0

      ! Find last slash or backslash
      pos = max(index(fullpath, '/'), index(fullpath, '\'))

      if (pos > 0) then
         folder = adjustl(fullpath(:pos - 1))
         call create_folder(trim(folder), ios)
         if (ios /= 0) return
      end if

      open (newunit=unit, file=trim(fullpath), status='replace', action='write', iostat=ios)
      if (ios == 0) close (unit)

   end subroutine create_file_with_path

   function get_path_separator() result(sep)
      character(len=1) :: sep

#ifdef _WIN32
      sep = '\'
#else
      sep = '/'
#endif

   end function get_path_separator

end module mod_directoryUtils
