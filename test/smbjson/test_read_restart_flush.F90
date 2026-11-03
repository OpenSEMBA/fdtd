integer function test_read_restart_flush_minutes() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use NFDETypes_m

   implicit none

   character(len=*), parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'planewave_restart_flush.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (pr%general%restartFlushMinutes /= 7) &
      call testFails(err, 'general.restartFlushMinutes should be read as 7')

end function

integer function test_read_restart_flush_default() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use NFDETypes_m

   implicit none

   character(len=*), parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'planewave.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (pr%general%restartFlushMinutes /= 0) &
      call testFails(err, 'general.restartFlushMinutes should default to 0 (no periodic checkpoint)')

end function
