integer function test_read_wire_probe_defaults() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'wireProbeDefaults.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser
   integer :: i, found
   err = 0
   found = 0

   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   ! A wire probe which does not specify a "field" entry must default to current.
#ifdef CompileWithMTLN
   do i = 1, size(problem%mtln%probes)
      if (problem%mtln%probes(i)%probe_name == "no_field_probe") then
         found = found + 1
         call expect_eq_int(err, PROBE_TYPE_CURRENT, problem%mtln%probes(i)%probe_type)
      end if
   end do
#else
   do i = 1, size(problem%Sonda%collection)
      if (problem%Sonda%collection(i)%outputrequest == "no_field_probe") then
         found = found + 1
         call expect_eq_int(err, NP_COR_WIRECURRENT, problem%Sonda%collection(i)%cordinates(1)%Or)
      end if
   end do
#endif
   call expect_eq_int(err, 1, found)

end function
