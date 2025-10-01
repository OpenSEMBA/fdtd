integer function test_init_solver() bind (C) result(err)
   use SEMBA_FDTD_mod
   use system_testingTools_mod   
   implicit none
   type(semba_fdtd_t) :: semba
   type(solver_t) :: solver
   real(kind=RKIND) :: field_value(1.0)
   err = 0
   call chdir("./test/system/")

   call semba%init("-i init_solver.fdtd.json")
   solver = semba%create_solver()
   call solver%init()
   call solver%set_field_value(iEx, [2,4], [2,2], [2,2], field_value)
   call solver%step()
   if (solver%get_field_value(iHy, 2,2,2) == 0) err = err + 1
   if (solver%get_field_value(iHz, 2,2,2) == 0) err = err + 1

   call solver%destroy_and_deallocate()
   call chdir("../../")
end function


integer function test_rank_remapping() bind (C) result(err)
   use SEMBA_FDTD_mod
   use system_testingTools_mod   
   implicit none

   type(semba_fdtd_t) :: semba
   type(solver_t) :: solver

   err = 0
   call chdir("./test/system/")

   call semba%init("-i init_solver.fdtd.json")
   solver = semba%create_solver()
   call solver%init()
   call solver%set_field_value(iHy, [2,2], [2,2], [2,2], 1.0)
   call solver%advanceEx(solver%media%sggMiEx)
   if (solver%get_field_value(iEx, 2,2,2) /= -33.8822708) err = err + 1
   if (solver%get_field_value(iEx, 2,2,3) /= 33.8822708) err = err + 1

   call solver%destroy_and_deallocate()

   call chdir("../../")

end function