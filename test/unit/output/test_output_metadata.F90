integer function test_atomic_file_replacement() bind(c) result(err)
   ! Verifies replacement leaves the target complete and removes the temporary file.
   use directoryUtils_m, only: atomic_replace_file, create_file_with_path, file_exists, &
                               delete_file
   use assertionTools_m, only: assert_integer_equal, assert_true, assert_string_equal
   implicit none

   character(len=*), parameter :: TARGETVALUE = 'testing atomic replacement/result.json'
   character(len=*), parameter :: TEMPORARY = 'testing atomic replacement/result.json.tmp'
   character(len=32) :: line
   integer :: ios, unit

   err = 0
   call create_file_with_path(targetValue, ios)
   call create_file_with_path(TEMPORARY, ios)
   open(newunit=unit, file=TEMPORARY, status='old', action='write', position='rewind', iostat=ios)
   write(unit, '(a)') 'complete'
   close(unit)

   call atomic_replace_file(TEMPORARY, targetValue, ios)
   err = err + assert_integer_equal(ios, 0, 'Atomic replacement failed')
   err = err + assert_true(file_exists(targetValue), 'Replacement target does not exist')
   err = err + assert_true(.not. file_exists(TEMPORARY), 'Temporary file remains after replacement')

   open(newunit=unit, file=targetValue, status='old', action='read', iostat=ios)
   read(unit, '(a)', iostat=ios) line
   close(unit)
   err = err + assert_string_equal(line, 'complete', 'Replacement target is incomplete')
   call delete_file(targetValue, ios)
end function test_atomic_file_replacement

integer function test_keep_previous_run_files() bind(c) result(err)
   ! Verifies a resumed run moves the files of the previous run aside instead of overwriting them.
   use directoryUtils_m, only: create_file_with_path, file_exists, folder_exists, remove_folder
   use outputUtils_m, only: keep_previous_run_files
   use assertionTools_m, only: assert_true, assert_string_equal
   implicit none

   character(len=*), parameter :: FOLDER = 'testing keep previous run'
   character(len=32) :: line
   integer :: ios, unit, run

   err = 0
   call remove_folder(FOLDER, ios)
   do run = 1, 2
      call create_file_with_path(FOLDER//'/probe.h5', ios)
      call create_file_with_path(FOLDER//'/probe.xdmf', ios)
      open(newunit=unit, file=FOLDER//'/probe.h5', status='old', action='write', position='rewind', iostat=ios)
      write(unit, '(a,i0)') 'run ', run
      close(unit)
      call keep_previous_run_files(FOLDER, 'probe')
      err = err + assert_true(.not. file_exists(FOLDER//'/probe.h5'), 'Previous run file was not moved')
   end do

   do run = 1, 2
      open(newunit=unit, file=FOLDER//'/before_resume_'//achar(48 + run)//'/probe.h5', status='old', action='read', iostat=ios)
      line = ''
      read(unit, '(a)', iostat=ios) line
      close(unit)
      err = err + assert_string_equal(line, 'run '//achar(48 + run), 'Previous run content was lost')
      err = err + assert_true(file_exists(FOLDER//'/before_resume_'//achar(48 + run)//'/probe.xdmf'), 'XDMF was not moved')
   end do
   call remove_folder(FOLDER, ios)
end function test_keep_previous_run_files
