integer function test_atomic_file_replacement() bind(c) result(err)
   ! Verifies replacement leaves the target complete and removes the temporary file.
   use directoryUtils_m, only: atomic_replace_file, create_file_with_path, file_exists, &
                               delete_file
   use assertionTools_m, only: assert_integer_equal, assert_true, assert_string_equal
   implicit none

   character(len=*), parameter :: targetValue = 'testing atomic replacement/result.json'
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
