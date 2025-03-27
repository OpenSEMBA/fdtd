
    
Module Getargs
   USE NFDETYPES , ONLY: BUFSIZE
   implicit none
   private 

   public getcommandargument,commandargumentcount

contains

   subroutine getcommandargument(chain2,posic,argum,length,status)
      character (LEN=BUFSIZE)  ::chain2, argum, argument, binaryPath
      integer (kind=4)  :: length, status, posic, binaryPathLenght, argumentStart, argumentEnd
      integer (kind=4) :: n, i, j

      CALL removeDoubleWhiteSpaces(chain2)

      CALL getarg(0, binaryPath)
      binaryPathLenght = len(trim(adjustl(binaryPath)))

      if (posic == 1) then 
         argum = binaryPath
         return
      end if
      status=0
      argumentStart=0
      argumentEnd=0
      n=1
      findStart : do i = binaryPathLenght, len(trim(adjustl(chain2)))+1
         if(chain2(i:i)== ' ') n=n+1
         if(n==posic) then
            do j=i+1,len(trim(adjustl(chain2)))
               if (chain2(j : j)/=' ') argumentStart=j
               exit findStart
            end do
         endif
      end do findStart

      findEnd :  do i=argumentStart + 1 ,len(trim(adjustl(chain2)))+2
         if (chain2(i : i)==' ') then
            argumentEnd=i-1
            exit findEnd
         endif
      end do findEnd

      if (argumentStart+argumentEnd==0) status=1
      argum=trim(adjustl(chain2(argumentStart : argumentEnd)))

      !100615 para evitar el crlf del .sh
      if ( (argum(1:1) ==char(10)) .or. (argum(1:1) ==char(13)) .or. (argum(1:1)==char( 0)) ) then
         argum=''
         return
      endif

      return
   end subroutine

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

   function commandargumentcount(chain2)
      character (LEN=BUFSIZE)  ::  chain2, binaryPath
      integer (kind=4)  ::  status,n,commandargumentcount, binaryPathLenght, i

      CALL removeDoubleWhiteSpaces(chain2)

      CALL getarg(0, binaryPath)
      binaryPathLenght = len(trim(adjustl(binaryPath)))

      n = 1

      do i=binaryPathLenght  ,len(trim(adjustl(chain2)))
         if (chain2(i : i)==' ') then
            n=n+1
         endif 
      end do

      status=0
      commandargumentcount=n
      return
   end function
end module
