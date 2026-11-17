integer function test_read_lorentz_material() bind(C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use NFDETypes_m

   implicit none

   character(len=*), parameter :: FILENAME = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'lorentz_slab.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   real(kind=RKIND), parameter :: TOL = 1.0e-5_RKIND
   complex(kind=CKIND) :: expected_pole1, expected_pole2, expected_residue1, expected_residue2

   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   ! The same material is associated with a volume, a surface and a line, so
   ! one dispersive entry is expected for each geometry type.
   call expect_eq_int(err, 1, pr%frqDepMats%nVols, 'frqDepMats should have one volume')
   call expect_eq_int(err, 1, pr%frqDepMats%nSurfs, 'frqDepMats should have one surface')
   call expect_eq_int(err, 1, pr%frqDepMats%nLins, 'frqDepMats should have one line')
   call expect_eq_int(err, 1, pr%frqDepMats%nVols_max, 'frqDepMats nVols_max')
   call expect_eq_int(err, 1, pr%frqDepMats%nSurfs_max, 'frqDepMats nSurfs_max')
   call expect_eq_int(err, 1, pr%frqDepMats%nLins_max, 'frqDepMats nLins_max')
   call expect_eq_int(err, 1, pr%frqDepMats%n_c_max, 'frqDepMats n_c_max')

   ! eps_r(w) = epsInf + omegaP^2 / (omega1^2 - w^2 + j*gamma*w) with
   ! epsInf = 2, omegaP = 3c, omega1 = 2c and gamma = 0.3c.
   ! The solver receives the equivalent first-order pole pair.
   expected_pole1 = cmplx(4.49688687e7_RKIND, 5.97896205325447083e8_RKIND, CKIND)
   expected_pole2 = cmplx(4.49688687e7_RKIND, -5.97896205325447083e8_RKIND, CKIND)
   expected_residue1 = cmplx(0.0_RKIND, 5.98931083300393435e-3_RKIND, CKIND)
   expected_residue2 = cmplx(0.0_RKIND, -5.98931083300393435e-3_RKIND, CKIND)

   call expect_lorentz(err, pr%frqDepMats%Vols(1), 0, 3, 0, 3, 20, 29, 0, &
      'plasma@volume', expected_pole1, expected_pole2, expected_residue1, expected_residue2, TOL)
   call expect_lorentz(err, pr%frqDepMats%Surfs(1), 0, 3, 0, 3, 10, 10, 3, &
      'plasma@surface', expected_pole1, expected_pole2, expected_residue1, expected_residue2, TOL)
   call expect_lorentz(err, pr%frqDepMats%Lins(1), 0, 0, 0, 0, 5, 14, 3, &
      'plasma@line', expected_pole1, expected_pole2, expected_residue1, expected_residue2, TOL)

contains
   subroutine expect_lorentz(err, entry, xi, xe, yi, ye, zi, ze, orientation, tag, &
                             pole1, pole2, residue1, residue2, tol)
      integer, intent(inout) :: err
      type(FreqDepenMaterial_t), intent(in) :: entry
      integer, intent(in) :: xi, xe, yi, ye, zi, ze, orientation
      character(len=*), intent(in) :: tag
      complex(kind=CKIND), intent(in) :: pole1, pole2, residue1, residue2
      real(kind=RKIND), intent(in) :: tol

      if (entry%n_c /= 1) then
         call testFails(err, 'Lorentz entry should have one coordinate')
         return
      end if
      if (entry%c(1)%Xi /= xi .or. entry%c(1)%Xe /= xe .or. &
          entry%c(1)%Yi /= yi .or. entry%c(1)%Ye /= ye .or. &
          entry%c(1)%Zi /= zi .or. entry%c(1)%Ze /= ze) then
         call testFails(err, 'Lorentz entry has unexpected coordinates')
      end if
      if (entry%c(1)%Or /= orientation) call testFails(err, 'Lorentz entry has unexpected orientation')
      if (trim(entry%c(1)%tag) /= tag) call testFails(err, 'Lorentz entry has unexpected tag')

      if (entry%K11 /= 2) call testFails(err, 'Lorentz entry should have two poles')
      if (abs(entry%eps11 - 2.0_RKIND*EPSILON_VACUUM) > tol*EPSILON_VACUUM) &
         call testFails(err, 'Lorentz epsInf should scale the vacuum permittivity')
      if (abs(entry%mu11 - MU_VACUUM) > tol*MU_VACUUM) &
         call testFails(err, 'Lorentz permeability should be the vacuum value')
      if (entry%sigma11 /= 0.0_RKIND) call testFails(err, 'Lorentz conductivity should be zero')

      if (abs(entry%b11(1) - pole1) > tol*abs(pole1)) call testFails(err, 'Lorentz first pole')
      if (abs(entry%b11(2) - pole2) > tol*abs(pole2)) call testFails(err, 'Lorentz second pole')
      if (abs(entry%a11(1) - residue1) > tol*abs(residue1)) call testFails(err, 'Lorentz first residue')
      if (abs(entry%a11(2) - residue2) > tol*abs(residue2)) call testFails(err, 'Lorentz second residue')
   end subroutine
end function
