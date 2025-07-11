integer function test_init_solver() bind (C) result(err)
   use SEMBA_FDTD_mod
   use system_testingTools_mod   
   implicit none
   type(semba_fdtd_t) :: semba
   type(solver_t) :: solver

   err = 0
   call chdir("./test/system/")

   call semba%init("-i init_solver.fdtd.json")
   call solver%init_control(semba%l, semba%maxSourceValue, semba%time_desdelanzamiento)
   call solver%init(semba%sgg,semba%eps0, semba%mu0, semba%sggMiNo,& 
                    semba%sggMiEx,semba%sggMiEy,semba%sggMiEz,& 
                    semba%sggMiHx,semba%sggMiHy,semba%sggMiHz, & 
                    semba%sggMtag, semba%SINPML_fullsize, semba%fullsize, semba%tag_numbers)
   call solver%set_field_value(iEx, [2,4], [2,2], [2,2], 1.0)
   call solver%step(semba%sgg, semba%eps0, semba%mu0, semba%SINPML_FULLSIZE, semba%tag_numbers)
   if (solver%get_field_value(iHy, 2,2,2) == 0) err = err + 1
   if (solver%get_field_value(iHz, 2,2,2) == 0) err = err + 1

   call chdir("../../")
end function