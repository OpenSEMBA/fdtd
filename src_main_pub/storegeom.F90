!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module to handle the storing of the geometry in ascii files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module storeData_m
   use FDETYPES_m
   !
   implicit none
   private
   !
   integer(kind=4), parameter, private :: BLOCK_SIZE = 1024
   public store_geomData
   !
contains
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Stores the geometrical data given by the parser into disk
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine store_geomData (sgg,media, fileFDE)
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) :: INTJ
      type(media_matrices_t), intent(in) :: media
      type(SGGFDTDINFO_t), intent(in) :: sgg
      integer(kind=4) :: i, j, k, fieldIndex, q
      character(len=*), intent(in) :: fileFDE
      !Writes an ASCII map of the media matrix for each field component
      open(20, FILE=trim(adjustl(fileFDE))//'_MapEx.txt')
      open(21, FILE=trim(adjustl(fileFDE))//'_MapEy.txt')
      open(22, FILE=trim(adjustl(fileFDE))//'_MapEz.txt')
      open(23, FILE=trim(adjustl(fileFDE))//'_MapHx.txt')
      open(24, FILE=trim(adjustl(fileFDE))//'_MapHy.txt')
      open(25, FILE=trim(adjustl(fileFDE))//'_MapHz.txt')
      do fieldIndex = 1, 6
         i = 19 + fieldIndex
         q = 19 + fieldIndex
         write(q,*) '____ 1-Sustrato, -n PML_______'
         do j = 0, sgg%NumMedia
            INTJ=J
            write(q,*) '_____________________________'
            write(q,*) 'MEDIO :  ', chartranslate (Intj)
            write(q,*) 'Priority ', sgg%Med(j)%Priority
            write(q,*) 'Epr ', sgg%Med(j)%Epr
            write(q,*) 'Sigma ', sgg%Med(j)%Sigma
            write(q,*) 'Mur ', sgg%Med(j)%Mur
            write(q,*) 'Is PML ', sgg%Med(j)%Is%PML
            write(q,*) 'Is PEC ', sgg%Med(j)%Is%PEC
            write(q,*) 'SigmaM ', sgg%Med(j)%SigmaM
            write(q,*) 'Is ThinWIRE ', sgg%Med(j)%Is%ThinWire
            write(q,*) 'Is SlantedWIRE ', sgg%Med(j)%Is%SlantedWire
            write(q,*) 'Is EDispersive ', sgg%Med(j)%Is%EDispersive
            write(q,*) 'Is MDispersive ', sgg%Med(j)%Is%MDispersive
            write(q,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
            write(q,*) 'Is SGBC ', sgg%Med(j)%Is%SGBC
            write(q,*) 'Is Lossy ', sgg%Med(j)%Is%Lossy
            write(q,*) 'Is Multiport ', sgg%Med(j)%Is%multiport
            write(q,*) 'Is AnisMultiport ', sgg%Med(j)%Is%anismultiport
            write(q,*) 'Is MultiportPadding ', sgg%Med(j)%Is%multiportpadding
            write(q,*) 'Is Dielectric ', sgg%Med(j)%Is%DIELECTRIC
            write(q,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
            write(q,*) 'Is Anisotropic ', sgg%Med(j)%Is%Anisotropic
            write(q,*) 'Is Needed ', sgg%Med(j)%Is%Needed
            write(q,*) 'Is already_YEEadvanced_byconformal ', sgg%Med(j)%Is%already_YEEadvanced_byconformal
            write(q,*) 'Is split_and_useless ', sgg%Med(j)%Is%split_and_useless
            write(q,*) 'Is Volume ', sgg%Med(j)%Is%Volume
            write(q,*) 'Is Surface ', sgg%Med(j)%Is%Surface
            write(q,*) 'Is Line ', sgg%Med(j)%Is%Line
         end do
         !
         write(i,*) fieldIndex, ' con PML IINIC, IFIN ', sgg%sweep(fieldIndex)%XI, sgg%sweep(fieldIndex)%XE
         write(i,*) fieldIndex, ' con PML JINIC, JFIN ', sgg%sweep(fieldIndex)%YI, sgg%sweep(fieldIndex)%YE
         write(i,*) fieldIndex, ' con PML KINIC, KFIN ', sgg%sweep(fieldIndex)%ZI, sgg%sweep(fieldIndex)%ZE
         write(i,*) fieldIndex, ' sin PML IINIC, IFIN ', sgg%SINPMLsweep(fieldIndex)%XI, sgg%SINPMLsweep(fieldIndex)%XE
         write(i,*) fieldIndex, ' sin PML JINIC, JFIN ', sgg%SINPMLsweep(fieldIndex)%YI, sgg%SINPMLsweep(fieldIndex)%YE
         write(i,*) fieldIndex, ' sin PML KINIC, KFIN ', sgg%SINPMLsweep(fieldIndex)%ZI, sgg%SINPMLsweep(fieldIndex)%ZE
         !
         do k = sgg%sweep(fieldIndex)%ZI, sgg%sweep(fieldIndex)%ZE
            i = 19 + fieldIndex
            write(i, '(A)') '_______________________________________________________________________'
            write(i,*) '!!!!!!** k=', k
            write(19+fieldIndex, '(A,400a)') 'I=  |', ('0123456789', i=sgg%Alloc(fieldIndex)%XI, sgg%Alloc(fieldIndex)%XE+10, 10)
            write(19+fieldIndex, '(A)') 'J______________________________________________________________________'
            do j = sgg%sweep(fieldIndex)%YE, sgg%sweep(fieldIndex)%YI, - 1
               select case (fieldIndex)
                case (iEx)
                  write(19+fieldIndex, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEx(i, j, k)), i=sgg%sweep(fieldIndex)%XI, &
                  & sgg%sweep(fieldIndex)%XE)
                case (iEy)
                  write(19+fieldIndex, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEy(i, j, k)), i=sgg%sweep(fieldIndex)%XI, &
                  & sgg%sweep(fieldIndex)%XE)
                case (IEZ)
                  write(19+fieldIndex, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEz(i, j, k)), i=sgg%sweep(fieldIndex)%XI, &
                  & sgg%sweep(fieldIndex)%XE)
                case (IHX)
                  write(19+fieldIndex, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHx(i, j, k)), i=sgg%sweep(fieldIndex)%XI, &
                  & sgg%sweep(fieldIndex)%XE)
                case (IHY)
                  write(19+fieldIndex, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHy(i, j, k)), i=sgg%sweep(fieldIndex)%XI, &
                  & sgg%sweep(fieldIndex)%XE)
                case (IHZ)
                  write(19+fieldIndex, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHz(i, j, k)), i=sgg%sweep(fieldIndex)%XI, &
                  & sgg%sweep(fieldIndex)%XE)
               end select
            end do
         end do
      end do
      do i = 20, 25
         close (i)
      end do
      !
      return
      !
   contains
      !
      !Function to translate media indexes into characters for the mapping files
      !
      function chartranslate (entero) result (chara)
         integer(kind=INTEGERSIZEOFMEDIAMATRICES) entero
         character(len=1) chara
         if (entero == 1) then
            chara = '_'
         else if (entero == 0) then
            chara = '0'
         else if (entero ==-1) then
            chara = '#'
         else
            chara = char (48+Abs(entero))
         end if
         return
      end function chartranslate
      !
   end subroutine store_geomData
   !
end module storeData_m
!
