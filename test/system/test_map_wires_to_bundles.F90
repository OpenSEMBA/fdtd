integer function test_map_wires_to_bundles() bind (C) result(err)
   use smbjson
   use system_testingTools_mod
   use mtln_solver_mod, mtln_solver_t => mtln_t
   
   implicit none

   character(len=*),parameter :: filename = PATH_TO_TEST_DATA//'system/mtln.fdtd.json'
   type(Parseador) :: problem, expected
   type(parser_t) :: parser
   type(mtln_solver_t) :: solver
   logical :: areSame
   err = 0

   ! expected = expectedProblemDescription()
   parser = parser_t(filename)
   problem = parser%readProblemDescription()
   solver = mtlnCtor(problem%mtln)
   ! solver = mtlnCtor(problem%mtln)
   ! call expect_eq(err, expected, problem)
end function

