integer function test_read_mtln_no_mtln_compilation() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   ! This function is used in a death test (EXPECT_DEATH).
   ! It is expected to abort via error stop when reading a JSON with multiwires
   ! in a binary compiled without MTLN support.
   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'mtln.fdtd.json'
   type(parser_t) :: parser
   type(Parseador_t) :: problem
   err = 0

   parser = parser_t(filename)
   problem = parser%readProblemDescription()

end function
