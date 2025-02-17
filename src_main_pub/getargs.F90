
    
Module Getargs
   USE NFDETYPES , ONLY: BUFSIZE
   implicit none
      character(LEN=10) :: sembaExecutableLabel = "semba-fdtd"
   private 

   public getcommandargument,commandargumentcount

contains

   subroutine getcommandargument(chain2,posic,argum,length,status)
      character (LEN=BUFSIZE)  ::  chain2, argum
      character(LEN=1) :: a
      integer (kind=4)  :: length, status, posic, argumentStart, argumentEnd, endPathExecutableIndex
      integer (kind=4) :: i,j,n

      !length is unused
      !first remove all multiple blanks
      CALL removeDoubleWhiteSpaces(chain2)

      !remove binary path from string
      !First, find the semba-fdtd list index
      endPathExecutableIndex = index(chain2, sembaExecutableLabel) + len(sembaExecutableLabel)

      status=0
      argumentStart=0
      argumentEnd=0

      if (posic == 1) then
         argum=trim(adjustl(chain2(:endPathExecutableIndex)))
         return
      endif
      chain2=' '//trim(adjustl(chain2))//' '
      n=1




      findStart :  do i=endPathExecutableIndex  ,len(trim(adjustl(chain2)))+2
         if (chain2(i : i)==' ') n=n+1
         if (n==posic) then
            do j=i+1,len(trim(adjustl(chain2)))
               if (chain2(j : j)/=' ') argumentStart=j
               exit findStart
            end do
         endif
      end do findStart


      findEnd :  do i=argumentStart + 1 ,len(trim(adjustl(chain2)))+2
         a = chain2(i : i)
         if (chain2(i : i)==' ') then
            argumentEnd=i-1
            continue
            exit findEnd
         endif
      end do findEnd
      if (argumentStart*argumentEnd==0) status=1
      argum=trim(adjustl(chain2(argumentStart : argumentEnd)))

      !100615 para evitar el crlf del .sh
      if ( (argum(1:1) ==char(10)) .or. (argum(1:1) ==char(13)) .or. (argum(1:1)==char( 0)) ) then
         argum=''
         return
      endif

      return
   end subroutine

   function commandargumentcount(chain2)
      character (LEN=BUFSIZE)  ::  chain2
      integer (kind=4)  ::  status,n,i,j,commandargumentcount, endPathExecutableIndex

      !!if (chain2(1 : 5)=='-----') then
      !!    n=command_argument_count()
      !!else
      !length is unused
      !first remove all multiple blanks
      CALL removeDoubleWhiteSpaces(chain2)

      endPathExecutableIndex = index(chain2, sembaExecutableLabel) + len(sembaExecutableLabel)
      n=1
      do i=endPathExecutableIndex  ,len(trim(adjustl(chain2)))
         if (chain2(i : i)==' ') then
            n=n+1
         endif
      end do
      status=0
      !!!!!!endif
      commandargumentcount=n
      return
   end function

   subroutine removeDoubleWhiteSpaces(chain2)
      integer(kind=4) :: i, j
      character(LEN=BUFSIZE) :: chain2

      do i=1,len(trim(adjustl(chain2)))
         if (chain2(i : i)==' ') then
            rebus: do j=i+1,len(trim(adjustl(chain2)))
               if (chain2(j : j)/=' ') then
                  chain2(i+1 :)=chain2(j :)
                  exit rebus
               endif
            end do rebus
         endif
      end do
   
      end subroutine

end module
