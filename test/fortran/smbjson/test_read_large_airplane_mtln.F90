integer function test_read_large_airplane_mtln() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'large_airplane_mtln.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   
   err = 0

   parser = parser_t(filename)

   pr = parser%readProblemDescription()
   
contains
   
end function

