integer function test_spice_connectors() bind (C) result(err)
   use smbjson
   use smbjson_testingTools
   use mtln_solver_mod, mtln_solver_t => mtln_t

   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'cases/spice_connectors.fdtd.json'
   type(Parseador) :: problem
   type(parser_t) :: parser
   type(mtln_solver_t) :: solver
   integer :: i
   
   err  = 0    
   parser = parser_t(filename)
   problem = parser%readProblemDescription()
   solver = mtlnCtor(problem%mtln)
   call solver%run()
   write(*,*) 'here'
   open(unit = 1, file =  'wire_panel.txt')
   write(1,*) 'time probe_1_I probe_1_V'
   do i = 1, size(solver%bundles(1)%probes(1)%t)
       write(1,*) solver%bundles(1)%probes(1)%t(i)," ", &
               solver%bundles(1)%probes(1)%val(i,1) ," ", &
               solver%bundles(1)%probes(2)%val(i,1)
   ! lines depending on the number and type of probes
   end do

end function