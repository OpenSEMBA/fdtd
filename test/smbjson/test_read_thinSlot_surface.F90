integer function test_read_thinSlot_surface() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools_m

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'thinSlot_surface.fdtd.json'
   type(Parseador_t) :: pr
   type(parser_t) :: parser
   integer :: i
   err = 0

   parser = parser_t(filename)
   pr = parser%readProblemDescription()

   if (pr%tSlots%n_tg /= 1) then
      call testFails(err, 'Expected one thin slot group')
      return
   end if

   ! The surface interval spans three cells along Y and one cell along X.
   if (pr%tSlots%tg(1)%n_tgc /= 3) then
      call testFails(err, 'Expected three surface thin-slot components')
      return
   end if

   if (abs(pr%tSlots%tg(1)%width - 1e-3_RKIND) > 1e-12_RKIND) then
      call testFails(err, 'Unexpected surface thin-slot width')
      return
   end if

   do i = 1, 3
      if (pr%tSlots%tg(1)%tgc(i)%i /= 5 .or. &
          pr%tSlots%tg(1)%tgc(i)%j /= 1 + i .or. &
          pr%tSlots%tg(1)%tgc(i)%k /= 5 .or. &
          pr%tSlots%tg(1)%tgc(i)%dir /= IEY .or. &
          pr%tSlots%tg(1)%tgc(i)%normal /= +IEZ .or. &
          .not. pr%tSlots%tg(1)%tgc(i)%isSurface) then
         call testFails(err, 'Unexpected surface thin-slot component')
         return
      end if
      if (trim(adjustl(pr%tSlots%tg(1)%tgc(i)%tag)) /= 'gap@surface-slot') then
         call testFails(err, 'Unexpected surface thin-slot tag')
         return
      end if
   end do

end function
