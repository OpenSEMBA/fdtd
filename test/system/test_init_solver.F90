integer function test_init_solver() bind (C) result(err)
   use SEMBA_FDTD_mod
   
   implicit none
   type(semba_fdtd_t) :: semba
   write(*,*) 'solver'
   call semba%init()
end function

