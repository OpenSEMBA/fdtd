integer function test_read_thinSlot_surface() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'thinSlot_surface.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   type(coords_t) :: c
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (pr%tSlots%n_tg /= 1) then
      call testFails(err, 'Expected one thin slot group')
      return
   end if

   if (abs(pr%tSlots%tg(1)%width - 1e-3_RKIND) > 1e-12_RKIND) then
      call testFails(err, 'Unexpected surface thin-slot width')
      return
   end if

   ! Surface-defined slots stay as raw intervals until preprocessing.
   if (pr%tSlots%tg(1)%n_tgc /= 0) then
      call testFails(err, 'Surface thin slot must not be expanded by the parser')
      return
   end if
   if (.not. allocated(pr%tSlots%tg(1)%surfaceCoords)) then
      call testFails(err, 'Expected raw surface thin-slot intervals')
      return
   end if
   if (size(pr%tSlots%tg(1)%surfaceCoords) /= 1) then
      call testFails(err, 'Expected one raw surface thin-slot interval')
      return
   end if

   c = pr%tSlots%tg(1)%surfaceCoords(1)
   if (c%Xi /= 5 .or. c%Xe /= 5 .or. &
       c%Yi /= 2 .or. c%Ye /= 4 .or. &
       c%Zi /= 5 .or. c%Ze /= 5 .or. &
       c%Or /= +IEZ) then
      call testFails(err, 'Unexpected raw surface thin-slot interval')
      return
   end if
   if (trim(adjustl(c%tag)) /= 'gap@surface-slot') then
      call testFails(err, 'Unexpected surface thin-slot tag')
      return
   end if

end function

integer function test_read_thinSlot_surface_negativeNormal() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*), parameter :: filename = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'thinSlot_surface_negative.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   type(coords_t) :: c
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (pr%tSlots%n_tg /= 1) then
      call testFails(err, 'Expected one thin slot group')
      return
   end if
   if (.not. allocated(pr%tSlots%tg(1)%surfaceCoords)) then
      call testFails(err, 'Expected raw surface thin-slot intervals')
      return
   end if
   if (size(pr%tSlots%tg(1)%surfaceCoords) /= 1) then
      call testFails(err, 'Expected one raw surface thin-slot interval')
      return
   end if

   c = pr%tSlots%tg(1)%surfaceCoords(1)
   if (c%Xi /= 5 .or. c%Xe /= 5 .or. &
       c%Yi /= 2 .or. c%Ye /= 4 .or. &
       c%Zi /= 5 .or. c%Ze /= 5 .or. &
       c%Or /= -IEZ) then
      call testFails(err, 'Unexpected negative-normal surface thin-slot interval')
      return
   end if

end function

integer function test_read_thinSlot_zeroLength() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m
   use Report_m, only: isFatalError, resetFatalError

   implicit none

   character(len=*), parameter :: filename = &
      PATH_TO_TEST_DATA//INPUT_EXAMPLES//'thinSlot_zeroLength.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   err = 0

   call resetFatalError()
   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (isFatalError()) then
      call testFails(err, 'A zero-length thin slot must not be a fatal error')
   end if
   if (pr%tSlots%n_tg /= 1) then
      call testFails(err, 'Expected one zero-length thin slot group')
   end if
   if (pr%tSlots%tg(1)%n_tgc /= 0) then
      call testFails(err, 'Expected zero components for a zero-length thin slot')
   end if
   if (allocated(pr%tSlots%tg(1)%surfaceCoords)) then
      call testFails(err, 'Zero-length thin slot must not keep surface intervals')
   end if
   call resetFatalError()

end function
