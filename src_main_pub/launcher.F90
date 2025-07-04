program SEMBA_FDTD_launcher
   use SEMBA_FDTD_mod
   implicit none

   type(semba_fdtd_t) :: semba
   call semba%init()
   write(*,*) 'solver'
end program SEMBA_FDTD_launcher


