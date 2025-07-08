integer function test_init_solver() bind (C) result(err)
   use SEMBA_FDTD_mod
   use system_testingTools_mod   
   implicit none
   type(semba_fdtd_t) :: semba
   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'holland1981.fdtd.json'

   call semba%init("-i "//filename)
   call semba%launch()
   call semba%end()
end function