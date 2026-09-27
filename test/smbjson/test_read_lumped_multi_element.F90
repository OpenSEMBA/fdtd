integer function test_read_lumped_multi_element() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'lumped_multi_element.fdtd.json'
   type(Parseador_t) :: problem
   type(parser_t) :: parser
   integer :: i, nResistors
   logical :: foundA, foundB, foundC
   err = 0

   parser = parser_t(filename)
   problem = parser%readProblemDescription()

   ! A single lumped association listing three elementIds must yield three
   ! lumped components, one per element.
   if (problem%dielRegs%nLins /= 3) then
      print *, 'test_read_lumped_multi_element FAILED: expected 3 lins, got', problem%dielRegs%nLins
      err = err + 1
   end if

   nResistors = 0
   foundA = .false.
   foundB = .false.
   foundC = .false.
   do i = 1, problem%dielRegs%nLins
      if (problem%dielRegs%lins(i)%resistor .and. problem%dielRegs%lins(i)%R == 100.0_RKIND) then
         nResistors = nResistors + 1
         if (trim(adjustl(problem%dielRegs%lins(i)%c2P(1)%tag)) == '100ohm_resistor@lumped_line_a') foundA = .true.
         if (trim(adjustl(problem%dielRegs%lins(i)%c2P(1)%tag)) == '100ohm_resistor@lumped_line_b') foundB = .true.
         if (trim(adjustl(problem%dielRegs%lins(i)%c2P(1)%tag)) == '100ohm_resistor@lumped_line_c') foundC = .true.
      end if
   end do

   if (nResistors /= 3) then
      print *, 'test_read_lumped_multi_element FAILED: expected 3 resistors, got', nResistors
      err = err + 1
   end if

   if (.not. (foundA .and. foundB .and. foundC)) then
      print *, 'test_read_lumped_multi_element FAILED: missing per-element lumped tags'
      err = err + 1
   end if

   if (err == 0) print *, 'test_read_lumped_multi_element PASSED'
end function
