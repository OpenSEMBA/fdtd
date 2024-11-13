integer function test_read_large_airplane_mtln() bind (C) result(err)
   use smbjson
   use smbjson_testingTools

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'cases/large_airplane_mtln.fdtd.json'
   type(Parseador) :: pr
   type(parser_t) :: parser
   
   err = 0

   parser = parser_t(filename)

   pr = parser%readProblemDescription()
   
contains
   
end function

