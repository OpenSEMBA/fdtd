integer function test_read_nodal_source_resistance_per_meter() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA// &
      'cases/nodalSource/nodalSource.fdtd.json'

   err = 0
   call expect_cable_resistance(err, filename, 10000.0_RKIND)
contains
   subroutine expect_cable_resistance(err, filename, expectedResistance)
      integer, intent(inout) :: err
      character(len=*), intent(in) :: filename
      real(kind=RKIND), intent(in) :: expectedResistance
      real(kind=RKIND) :: tol
      type(Parseador_t) :: problem
      type(parser_t) :: parser

      tol = max(abs(expectedResistance) * 1.0e-6_RKIND, 1.0e-6_RKIND)

      parser = parser_t(filename)
      problem = parser%readProblemDescription()

#ifdef CompileWithMTLN
      call expect_eq_int(err, 1, size(problem%mtln%cables), 'Expected one parsed MTLN cable')
      select type(cable => problem%mtln%cables(1)%ptr)
      type is (unshielded_multiwire_t)
         call expect_eq_int(err, 10, size(cable%step_size), 'Expected 10 cable steps')
         if (abs(cable%resistance_per_meter(1,1) - expectedResistance) > tol) then
            call testFails(err, 'Parsed cable resistance per meter is not correct')
         end if
      class default
         call testFails(err, 'Expected an unshielded multiwire cable')
      end select
#else
      call expect_eq_int(err, 1, problem%tWires%n_tw, 'Expected one parsed thin wire')
      call expect_eq_int(err, 10, problem%tWires%tw(1)%n_twc, 'Expected 10 thin-wire cells')
      if (abs(problem%tWires%tw(1)%res - expectedResistance) > tol) then
         call testFails(err, 'Parsed thin-wire resistance per meter is not correct')
      end if
#endif
   end subroutine expect_cable_resistance
end function

integer function test_read_nodal_source_total_resistance() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA// &
      'cases/nodalSource/nodalSource_totalResistance.fdtd.json'

   err = 0
   call expect_cable_resistance(err, filename, 10000.0_RKIND)
contains
   subroutine expect_cable_resistance(err, filename, expectedResistance)
      integer, intent(inout) :: err
      character(len=*), intent(in) :: filename
      real(kind=RKIND), intent(in) :: expectedResistance
      real(kind=RKIND) :: tol
      type(Parseador_t) :: problem
      type(parser_t) :: parser

      tol = max(abs(expectedResistance) * 1.0e-6_RKIND, 1.0e-6_RKIND)

      parser = parser_t(filename)
      problem = parser%readProblemDescription()

#ifdef CompileWithMTLN
      call expect_eq_int(err, 1, size(problem%mtln%cables), 'Expected one parsed MTLN cable')
      select type(cable => problem%mtln%cables(1)%ptr)
      type is (unshielded_multiwire_t)
         call expect_eq_int(err, 10, size(cable%step_size), 'Expected 10 cable steps')
         if (abs(cable%resistance_per_meter(1,1) - expectedResistance) > tol) then
            call testFails(err, 'Parsed totalResistance override did not produce the expected resistance per meter')
         end if
      class default
         call testFails(err, 'Expected an unshielded multiwire cable')
      end select
#else
      call expect_eq_int(err, 1, problem%tWires%n_tw, 'Expected one parsed thin wire')
      call expect_eq_int(err, 10, problem%tWires%tw(1)%n_twc, 'Expected 10 thin-wire cells')
      if (abs(problem%tWires%tw(1)%res - expectedResistance) > tol) then
         call testFails(err, 'Parsed totalResistance override did not produce the expected resistance per meter')
      end if
#endif
   end subroutine expect_cable_resistance
end function