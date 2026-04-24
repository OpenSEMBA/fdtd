integer function test_read_holland1981_short() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools
   use mtln_types_m, only: TERMINATION_OPEN

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'holland1981_short.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser

   err = 0

   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   ! MTLN parsing must keep parity with legacy wire handling for this case.
   call expect_eq_int(err, 2, size(problem%mtln%networks), 'Expected two MTLN terminal networks')
   call expect_eq_int(err, TERMINATION_OPEN, &
      problem%mtln%networks(1)%connections(1)%nodes(1)%termination%termination_type, &
      'Wire short terminal should be normalized to open (initial side)')
   call expect_eq_int(err, TERMINATION_OPEN, &
      problem%mtln%networks(2)%connections(1)%nodes(1)%termination%termination_type, &
      'Wire short terminal should be normalized to open (end side)')

end function
