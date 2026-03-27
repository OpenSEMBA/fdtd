integer function test_unrecognized_key_at_root() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools
   use Report_m

   implicit none

   character(len=*), parameter :: filename = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'unrecognized_key_root.fdtd.json'
   type(parser_t) :: parser
   type(Parseador_t) :: pr
   err = 0

   call resetFatalError()
   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (.not. isFatalError()) then
      call testFails(err, 'Expected fatal error for unrecognized root key, but none was raised.')
   end if
   call resetFatalError()

end function

integer function test_unrecognized_key_in_section() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools
   use Report_m

   implicit none

   character(len=*), parameter :: filename = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'unrecognized_key_section.fdtd.json'
   type(parser_t) :: parser
   type(Parseador_t) :: pr
   err = 0

   call resetFatalError()
   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (.not. isFatalError()) then
      call testFails(err, 'Expected fatal error for unrecognized key in section, but none was raised.')
   end if
   call resetFatalError()

end function
