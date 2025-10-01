
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module to handle the storing of the geometry in ascii files
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE storeData
   USE fdetypes
   !
   IMPLICIT NONE
   PRIVATE
   !
   INTEGER (KIND=4), PARAMETER, PRIVATE :: BLOCK_SIZE = 1024
   PUBLIC store_geomData
   !
CONTAINS
   !
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Stores the geometrical data given by the parser into disk
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   SUBROUTINE store_geomData (sgg,media, fileFDE)
      INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) :: INTJ
      type(media_matrices_t), intent(in) :: media
      type (SGGFDTDINFO), intent(IN) :: sgg
      INTEGER (KIND=4) :: i, j, k, campo, q
      CHARACTER (LEN=*), INTENT (IN) :: fileFDE
      !Writes an ASCII map of the media matrix for each field component
      OPEN (20, FILE=trim(adjustl(fileFDE))//'_MapEx.txt')
      OPEN (21, FILE=trim(adjustl(fileFDE))//'_MapEy.txt')
      OPEN (22, FILE=trim(adjustl(fileFDE))//'_MapEz.txt')
      OPEN (23, FILE=trim(adjustl(fileFDE))//'_MapHx.txt')
      OPEN (24, FILE=trim(adjustl(fileFDE))//'_MapHy.txt')
      OPEN (25, FILE=trim(adjustl(fileFDE))//'_MapHz.txt')
      DO campo = 1, 6
         i = 19 + campo
         q = 19 + campo
         WRITE (q,*) '____ 1-Sustrato, -n PML_______'
         DO j = 0, sgg%NumMedia
            INTJ=J
            WRITE (q,*) '_____________________________'
            WRITE (q,*) 'MEDIO :  ', chartranslate (Intj)
            WRITE (q,*) 'Priority ', sgg%Med(j)%Priority
            WRITE (q,*) 'Epr ', sgg%Med(j)%Epr
            WRITE (q,*) 'Sigma ', sgg%Med(j)%Sigma
            WRITE (q,*) 'Mur ', sgg%Med(j)%Mur
            WRITE (q,*) 'Is PML ', sgg%Med(j)%Is%PML
            WRITE (q,*) 'Is PEC ', sgg%Med(j)%Is%PEC
            WRITE (q,*) 'SigmaM ', sgg%Med(j)%SigmaM
            WRITE (q,*) 'Is ThinWIRE ', sgg%Med(j)%Is%ThinWire
            WRITE (q,*) 'Is SlantedWIRE ', sgg%Med(j)%Is%SlantedWire
            WRITE (q,*) 'Is EDispersive ', sgg%Med(j)%Is%EDispersive
            WRITE (q,*) 'Is MDispersive ', sgg%Med(j)%Is%MDispersive
            WRITE (q,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
            WRITE (q,*) 'Is SGBC ', sgg%Med(j)%Is%SGBC
            WRITE (q,*) 'Is Lossy ', sgg%Med(j)%Is%Lossy
            WRITE (q,*) 'Is Multiport ', sgg%Med(j)%Is%multiport
            WRITE (q,*) 'Is AnisMultiport ', sgg%Med(j)%Is%anismultiport
            WRITE (q,*) 'Is MultiportPadding ', sgg%Med(j)%Is%multiportpadding
            WRITE (q,*) 'Is Dielectric ', sgg%Med(j)%Is%dielectric
            WRITE (q,*) 'Is ThinSlot ', sgg%Med(j)%Is%ThinSlot
            WRITE (q,*) 'Is Anisotropic ', sgg%Med(j)%Is%Anisotropic
            WRITE (q,*) 'Is Needed ', sgg%Med(j)%Is%Needed
            WRITE (q,*) 'Is already_YEEadvanced_byconformal ', sgg%Med(j)%Is%already_YEEadvanced_byconformal
            WRITE (q,*) 'Is split_and_useless ', sgg%Med(j)%Is%split_and_useless
            WRITE (q,*) 'Is Volume ', sgg%Med(j)%Is%Volume
            WRITE (q,*) 'Is Surface ', sgg%Med(j)%Is%Surface
            WRITE (q,*) 'Is Line ', sgg%Med(j)%Is%Line
         END DO
         !
         WRITE (i,*) campo, ' con PML IINIC, IFIN ', sgg%sweep(campo)%XI, sgg%sweep(campo)%XE
         WRITE (i,*) campo, ' con PML JINIC, JFIN ', sgg%sweep(campo)%YI, sgg%sweep(campo)%YE
         WRITE (i,*) campo, ' con PML KINIC, KFIN ', sgg%sweep(campo)%ZI, sgg%sweep(campo)%ZE
         WRITE (i,*) campo, ' sin PML IINIC, IFIN ', sgg%SINPMLsweep(campo)%XI, sgg%SINPMLsweep(campo)%XE
         WRITE (i,*) campo, ' sin PML JINIC, JFIN ', sgg%SINPMLsweep(campo)%YI, sgg%SINPMLsweep(campo)%YE
         WRITE (i,*) campo, ' sin PML KINIC, KFIN ', sgg%SINPMLsweep(campo)%ZI, sgg%SINPMLsweep(campo)%ZE
         !
         DO k = sgg%sweep(campo)%ZI, sgg%sweep(campo)%ZE
            i = 19 + campo
            WRITE (i, '(A)') '_______________________________________________________________________'
            WRITE (i,*) '!!!!!!** k=', k
            WRITE (19+campo, '(A,400a)') 'I=  |', ('0123456789', i=sgg%Alloc(campo)%XI, sgg%Alloc(campo)%XE+10, 10)
            WRITE (19+campo, '(A)') 'J______________________________________________________________________'
            DO j = sgg%sweep(campo)%YE, sgg%sweep(campo)%YI, - 1
               SELECT CASE (campo)
                CASE (iEx)
                  WRITE (19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEx(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iEy)
                  WRITE (19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEy(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iEz)
                  WRITE (19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiEz(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iHx)
                  WRITE (19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHx(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iHy)
                  WRITE (19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHy(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
                CASE (iHz)
                  WRITE (19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(media%sggMiHz(i, j, k)), i=sgg%sweep(campo)%XI, &
                  & sgg%sweep(campo)%XE)
               END SELECT
            END DO
         END DO
      END DO
      DO i = 20, 25
         CLOSE (i)
      END DO
      !
      RETURN
      !
   CONTAINS
      !
      !Function to translate media indexes into characters for the mapping files
      !
      FUNCTION chartranslate (entero) RESULT (chara)
         INTEGER (KIND=INTEGERSIZEOFMEDIAMATRICES) entero
         CHARACTER (LEN=1) chara
         IF (entero == 1) THEN
            chara = '_'
         ELSE IF (entero == 0) THEN
            chara = '0'
         ELSE IF (entero ==-1) THEN
            chara = '#'
         ELSE
            chara = char (48+Abs(entero))
         END IF
         RETURN
      END FUNCTION chartranslate
      !
   END SUBROUTINE store_geomData
   !
END MODULE storeData
!
