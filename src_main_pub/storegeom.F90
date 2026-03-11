!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module to handle the storing of the geometry in ascii files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module storeData
   use fdetypes
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
      integer(kind=4) :: i, j, k, campo, q
      character(len=*), intent(in) :: fileFDE
      !Writes an ASCII map of the media matrix for each field component
      open(20, FILE=trim(adjustl(fileFDE))//'_MapEx.txt')
      open(21, FILE=trim(adjustl(fileFDE))//'_MapEy.txt')
      open(22, FILE=trim(adjustl(fileFDE))//'_MapEz.txt')
      open(23, FILE=trim(adjustl(fileFDE))//'_MapHx.txt')
      open(24, FILE=trim(adjustl(fileFDE))//'_MapHy.txt')
      open(25, FILE=trim(adjustl(fileFDE))//'_MapHz.txt')
      do campo = 1, 6
         i = 19 + campo
         q = 19 + campo
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
            write(q,*) 'Is Dielectric ', sgg%Med(j)%Is%dielectric
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
         write(i,*) campo, ' con PML IINIC, IFIN ', sgg%sweep(campo)%XI, sgg%sweep(campo)%XE
         write(i,*) campo, ' con PML JINIC, JFIN ', sgg%sweep(campo)%YI, sgg%sweep(campo)%YE
         write(i,*) campo, ' con PML KINIC, KFIN ', sgg%sweep(campo)%ZI, sgg%sweep(campo)%ZE
         write(i,*) campo, ' sin PML IINIC, IFIN ', sgg%SINPMLsweep(campo)%XI, sgg%SINPMLsweep(campo)%XE
         write(i,*) campo, ' sin PML JINIC, JFIN ', sgg%SINPMLsweep(campo)%YI, sgg%SINPMLsweep(campo)%YE
         write(i,*) campo, ' sin PML KINIC, KFIN ', sgg%SINPMLsweep(campo)%ZI, sgg%SINPMLsweep(campo)%ZE
         !
         do k = sgg%sweep(campo)%ZI, sgg%sweep(campo)%ZE
            i = 19 + campo
            write(i, '(A)') '_______________________________________________________________________'
            write(i,*) '!!!!!!** k=', k
            write(19+campo, '(A,400a)') 'I=  |', ('0123456789', i=sgg%Alloc(campo)%XI, sgg%Alloc(campo)%XE+10, 10)
            write(19+campo, '(A)') 'J______________________________________________________________________'
            do j = sgg%sweep(campo)%YE, sgg%sweep(campo)%YI, - 1
               SELECT CASE (campo)
                CASE (iEx)
                  write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEx(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iEy)
                  write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEy(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iEz)
                  write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEz(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iHx)
                  write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHx(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iHy)
                  write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHy(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iHz)
                  write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHz(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
               end select
            end do
         end do
      end do
      do i = 20, 25
         CLOSE (i)
      end do
      !
      return
      !
   contains
      !
      !Function to translate media indexes into characters for the mapping files
      !
      function chartranslate (entero) RESULT (chara)
         integer(kind=INTEGERSIZEOFMEDIAMATRICES) entero
         character(len=1) chara
         if (entero == 1) then
            chara = '_'
         ELSE if (entero == 0) then
            chara = '0'
         ELSE if (entero ==-1) then
            chara = '#'
         ELSE
            chara = char (48+Abs(entero))
         end if
         return
      end function chartranslate
      !
   end subroutine store_geomData
   !
end module storeData
!
