integer function test_read_thinSlot_ez() bind (C) result(err)
   use smbjson_m
   use smbjson_testingTools

   implicit none

   character(len=*), parameter :: filename = PATH_TO_TEST_DATA//INPUT_EXAMPLES//'thinSlot_ez.fdtd.json'
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

   ! Interval [[5,5,8],[5,5,14]], 6 Ez linels
   if (pr%tSlots%tg(1)%n_tgc /= 6) then
      call testFails(err, 'Expected six Ez thin-slot components')
      return
   end if

   if (abs(pr%tSlots%tg(1)%width - 1e-3_RKIND) > 1e-12_RKIND) then
      call testFails(err, 'Unexpected thin-slot width')
      return
   end if

   do i = 1, 6
      if (pr%tSlots%tg(1)%tgc(i)%i /= 5 .or. &
          pr%tSlots%tg(1)%tgc(i)%j /= 5 .or. &
          pr%tSlots%tg(1)%tgc(i)%k /= 7 + i .or. &
          pr%tSlots%tg(1)%tgc(i)%dir /= iEz) then
         call testFails(err, 'Unexpected Ez thin-slot component coordinates')
         return
      end if
      if (trim(adjustl(pr%tSlots%tg(1)%tgc(i)%tag)) /= 'gap@vertical-slot') then
         call testFails(err, 'Unexpected Ez thin-slot tag')
         return
      end if
   end do

end function
