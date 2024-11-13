!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MIT License
! 
! Copyright (c) 2023 University of Granada
! 
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
! 
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
! 
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DMMA
   USE FDETYPES
   IMPLICIT NONE
   PRIVATE
   !
   PUBLIC dmma_thin_Slot
   !
CONTAINS
   !
   !subroutine dmma_thin_Slot (incx,incy,incz,dir,orientacion,direccion,&
   !                           thickness,efm,ufm,epr,mur,epse,mue)
   !
   SUBROUTINE dmma_thin_Slot (incx, incy, incz, dir, orientacion, direccion, thickness, efm, ufm, epse, mue,eps0,mu0)
      !
      !
      !----------------------------------------------------------------------------------------------------
      ! This subroutine implement the DMMA thin Slot aproximation.
      ! It computes the equivalent epsilon and mu to be assigned to cells covered by the Slot.
      ! Imput data:
      !   incx        -> Cell size in x direction (meters)
      !   incy        -> Cell size in y direction (meters)
      !   incz        -> Cell size in z direction (meters)
      !   dir         -> Coordinates of the direction of incidence of the plane wave (meters)
      !   orientacion       -> orientacion ('XY','XZ','YZ') that contains the Slot line
      !   direccion   -> direccion ('X','Y','Z') of the Slot line
      !   thickness   -> Lower size of the Slot (meters)
      !   efm         -> Relative epsilon of the filling media
      !   ufm         -> Relative mu of the filling media
      !   epr        -> Relative espsilon of the media in rank 2 tensor format
      !   mur         -> Relative mu of the media in rank 2 tensor format
      ! Output data:
      !   epse        -> Equivalent espsilon of the Slot line in rank 2 tensor format
      !   mue         -> Equivalent mu of the Slot line in rank 2 tensor format
      !----------------------------------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL (KIND=RKIND)           ::  eps0,mu0
      REAL (KIND=RKIND), INTENT (IN) :: incx, incy, incz
      INTEGER (KIND=4), INTENT (IN) :: orientacion, direccion
      REAL (KIND=RKIND), INTENT (IN) :: thickness
      REAL (KIND=RKIND), INTENT (IN) :: efm, ufm
      REAL (KIND=RKIND), INTENT (IN), DIMENSION (3) :: dir
      !real (kind=RKIND),intent(in),dimension(3,3) :: epr
      !real (kind=RKIND),intent(in),dimension(3,3) :: mur
      REAL (KIND=RKIND), INTENT (OUT), DIMENSION (3, 3) :: epse, mue
      !
      ! Maximum frequency allowed by the cell size
      REAL (KIND=RKIND) :: maxfreq
      ! Absolut epsilon and mu of the filling media
      REAL (KIND=RKIND) :: eabs, uabs
      ! Speed of light in the filling media
      REAL (KIND=RKIND) :: cfm
      REAL (KIND=RKIND) :: omega
      ! Capacitance of the Slot line
      REAL (KIND=RKIND) :: cap,E0,U0
      !
      E0=EPS0
      U0=MU0

      eabs = e0 * efm
      uabs = u0 * ufm
      !
      !SGG
      epse = 0.0_RKIND
      mue = 0.0_RKIND
      epse (1, 1) = efm
      epse (2, 2) = efm
      epse (3, 3) = efm
      mue (1, 1) = ufm
      mue (2, 2) = ufm
      mue (3, 3) = ufm
      !
      cfm = 1.0_RKIND /  Sqrt (eabs*uabs) !si lo tomo relativo a la direccion de incidencia puede ser cfm=0.0_RKIND y se jode el logaritmo.
      !asi que lo tomo fijo !2012 bug articulo1_tgap_sgg_stair
      IF (orientacion == iEz) THEN
         !        cfm = Abs (dir(3)) / Sqrt (eabs*uabs)
         maxfreq = cfm / (incz*10.0)
         omega = 2.0_RKIND * pi * maxfreq
         !!!!OLD pre 2011
         !  cap=(4.232*eabs)/pi-(2.0_RKIND *eabs*thickness)/(pi*cfm)*(log(omega*thickness/cfm)-1.0_RKIND)
         !2011 mathem
         cap = eabs * (0.9918536053486919-0.3183098861837907*Log((omega*thickness)/cfm))
         IF (direccion == iEx) THEN
            epse (2, 2) = (incy/incz) * (cap/eabs)
            mue (3, 3) = 1.0_RKIND / epse (2, 2)
         END IF
         IF (direccion == iEy) THEN
            epse (1, 1) = (incx/incz) * (cap/eabs)
            mue (3, 3) = 1.0_RKIND / epse (1, 1)
         END IF
      END IF
      !
      IF (orientacion == iEy) THEN
         !        cfm = Abs (dir(2)) / Sqrt (eabs*uabs)
         maxfreq = cfm / (incy*10.0)
         omega = 2.0_RKIND * pi * maxfreq
         !  cap=(4.232*eabs)/pi-(2.0_RKIND *eabs*thickness)/(pi*cfm)*(log(omega*thickness/cfm)-1.0_RKIND)
         !
         !2011 mathem
         cap = eabs * (0.9918536053486919-0.3183098861837907*Log((omega*thickness)/cfm))
         IF (direccion == iEx) THEN
            epse (3, 3) = (incz/incy) * (cap/eabs)
            mue (2, 2) = 1.0_RKIND / epse (3, 3)
         END IF
         IF (direccion == iEz) THEN
            epse (1, 1) = (incx/incy) * (cap/eabs)
            mue (2, 2) = 1.0_RKIND / epse (1, 1)
         END IF
      END IF
      !
      IF (orientacion == iEx) THEN
         !        cfm = Abs (dir(1)) / Sqrt (eabs*uabs)
         maxfreq = cfm / (incx*10.0)
         omega = 2.0_RKIND * pi * maxfreq
         !  cap=(4.232*eabs)/pi-(2.0_RKIND *eabs*thickness)/(pi*cfm)*(log(omega*thickness/cfm)-1.0_RKIND)
         !
         !2011 mathem
         cap = eabs * (0.9918536053486919-0.3183098861837907*Log((omega*thickness)/cfm))
         IF (direccion == iEy) THEN
            epse (3, 3) = (incz/incx) * (cap/eabs)
            mue (1, 1) = 1.0_RKIND / epse (3, 3)
         END IF
         IF (direccion == iEz) THEN
            epse (2, 2) = (incy/incx) * (cap/eabs)
            mue (1, 1) = 1.0_RKIND / epse (2, 2)
         END IF
      END IF
      !
      RETURN
   END SUBROUTINE dmma_thin_Slot
   !
   !
   !
   !

END MODULE
