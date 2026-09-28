integer function test_read_nodal_source_zero_length() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*), parameter :: FILENAME = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'nodal_source_zero_length.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser
   err = 0

   call resetFatalError()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   if (.not. isFatalError()) &
      call testFails(err, 'Expected a fatal error for nodal source defined over a zero-length interval')

end function

integer function test_read_nodal_source_non_line_interval() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*), parameter :: FILENAME = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'nodal_source_non_line_interval.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser
   err = 0

   call resetFatalError()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   if (.not. isFatalError()) &
      call testFails(err, 'Expected a fatal error for nodal source defined over non-line intervals')

end function

integer function test_read_nodal_source_one_cell_interval() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*), parameter :: FILENAME = PATH_TO_TEST_DATA// &
      'cases/nodalSource/nodal-source-with-movie.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser
   err = 0

   call resetFatalError()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   if (isFatalError()) &
      call testFails(err, 'Did not expect a fatal error for a one-cell nodal source')
   call expect_eq_int(err, 1, problem%nodSrc%n_nodSrc, 'Expected one parsed nodal source')

end function
