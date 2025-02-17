
    
Module Getargs
   USE NFDETYPES , ONLY: BUFSIZE
   implicit none
   private 

   public getcommandargument,commandargumentcount

contains

   subroutine getcommandargument(chain2,posic,argum,length,status)
      character (LEN=BUFSIZE)  ::chain2, argum, argument
      integer (kind=4)  :: length, status, posic
      integer (kind=4) :: n

      CALL getarg(posic, argument)
      
      argum = argument
      status=0

      !100615 para evitar el crlf del .sh
      if ( (argum(1:1) ==char(10)) .or. (argum(1:1) ==char(13)) .or. (argum(1:1)==char( 0)) ) then
         argum=''
         return
      endif

      return
   end subroutine

   function commandargumentcount(chain2)
      character (LEN=BUFSIZE)  ::  chain2
      integer (kind=4)  ::  status,n,commandargumentcount
      n = command_argument_count()  
      status=0
      commandargumentcount=n
      return
   end function
end module
