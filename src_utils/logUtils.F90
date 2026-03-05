module mod_logUtils
   implicit none

contains
   subroutine printMessage(layoutNumber, message)
      implicit none
      integer, intent(in) :: layoutnumber
      character(len=*), intent(in) :: message
      logical :: printea
      printea = .true.

      ! Print into console
      if (printea) then
         write (*, '(a)') trim(adjustl(message))
      end if

      ! Print into unitFile 11
      if (layoutnumber == 0) then
         write (11, '(a)') trim(adjustl(message))
      end if
   end subroutine printMessage

   subroutine printMessageWithSeparator(layoutnumber, message)
      implicit none
      integer, intent(in) :: layoutnumber
      character(len=*), intent(in) :: message
      logical :: printea
      character(len=50) :: SEPARADOR
      SEPARADOR = repeat('_', 50)
      printea = .true.

      ! Print into console
      if (printea) then
         write (*, '(a)') SEPARADOR
         write (*, '(a)') trim(adjustl(message))
         write (*, '(a)') SEPARADOR
      end if

      ! Print into unitFile 11
      if (layoutnumber == 0) then
         write (11, '(a)') SEPARADOR
         write (11, '(a)') trim(adjustl(message))
         write (11, '(a)') SEPARADOR
      end if

   end subroutine printMessageWithSeparator

   subroutine printMessageWithEndingSeparator(layoutNumber, message)
      implicit none
      integer, intent(in) :: layoutnumber
      character(len=*), intent(in) :: message
      logical :: printea
      character(len=50) :: SEPARADOR
      SEPARADOR = repeat('_', 50)
      printea = .true.

      ! Print into console
      if (printea) then
         write (*, '(a)') SEPARADOR
         write (*, '(a)') trim(adjustl(message))
         write (*, '(a)') SEPARADOR
      end if

      ! Print into unitFile 11
      if (layoutnumber == 0) then
         write (11, '(a)') SEPARADOR
         write (11, '(a)') trim(adjustl(message))
         write (11, '(a)') SEPARADOR
      end if
   end subroutine printMessageWithEndingSeparator

   subroutine printSeparator(layoutnumber)
      implicit none
      integer, intent(in) :: layoutnumber
      logical :: printea
      character(len=50) :: SEPARADOR
      SEPARADOR = repeat('_', 50)
      printea = .true.

      ! Print into console
      if (printea) then
         write (*, '(a)') SEPARADOR
      end if

      ! Print into unitFile 11
      if (layoutnumber == 0) then
         write (11, '(a)') SEPARADOR
      end if
   end subroutine printSeparator
end module mod_logUtils
