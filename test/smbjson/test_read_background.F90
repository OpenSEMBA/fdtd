integer function test_read_background_defaults() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools
   use NFDETypes_m

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'planewave.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (pr%mats%mats(1)%eps /= EPSILON_VACUUM) &
      call testFails(err, 'Default background permittivity should be EPSILON_VACUUM')
   if (pr%mats%mats(1)%mu /= MU_VACUUM) &
      call testFails(err, 'Default background permeability should be MU_VACUUM')

end function

integer function test_read_background_set() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools
   use NFDETypes_m

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'background.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   real :: eps, mu
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   eps = 1.7708e-11
   mu  = 2.5133e-6

   if (.not. expect_near(real(pr%mats%mats(1)%eps), eps, eps * 1e-4)) &
      call testFails(err, 'Background permittivity not set correctly')
   if (.not. expect_near(real(pr%mats%mats(1)%mu),  mu,  mu  * 1e-4)) &
      call testFails(err, 'Background permeability not set correctly')

end function
