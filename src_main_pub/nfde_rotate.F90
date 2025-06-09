MODULE nfde_rotate_m
    !       
   USE NFDETypes

   !
   IMPLICIT NONE
   !
   PUBLIC

CONTAINS
    
   SUBROUTINE nfde_rotate (this,mpidir) 
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir
      print *,'-----> Warning novel 2024 rotate routine......'
      call  rotate_generateSpaceSteps                (this, mpidir)
      call  rotate_generateCurrent_Field_Sources     (this, mpidir)
      call  rotate_generatePlaneWaves                (this, mpidir)
      call  rotate_generateBoxSources                (this, mpidir)
      call  rotate_generateFronteras                 (this, mpidir)
      call  rotate_generatePECs                      (this, mpidir)
      call  rotate_generatePMCs                      (this, mpidir)
      call  rotate_generateNONMetals                 (this, mpidir)
      call  rotate_generateANISOTROPICs              (this, mpidir)
      call  rotate_generateThinWires                 (this, mpidir)
      call  rotate_generateSlantedWires              (this, mpidir)
      call  rotate_generateThinSlots                 (this, mpidir)
      call  rotate_generateLossyThinSurface          (this, mpidir)
      call  rotate_generateFDMs                      (this, mpidir)
      call rotate_generateSONDAs                     (this, mpidir)
      call rotate_generateMasSondas                  (this, mpidir)
      call rotate_generateBloqueProbes               (this, mpidir)
      call rotate_generateVolumicProbes              (this, mpidir)



      RETURN
   END SUBROUTINE nfde_rotate

   SUBROUTINE rotate_generateSpaceSteps (this, mpidir)
      TYPE (Parseador), INTENT (INOUT) :: this          
      INTEGER (KIND=4) ::  mpidir
      TYPE (Desplazamiento), POINTER :: old_despl => NULL ()      
      TYPE (MatrizMedios), POINTER :: old_matriz => NULL ()
      integer (kind=4) :: oxi,oyi,ozi
      REAL (KIND=RK) :: roxi,royi,rozi
      REAL (KIND=RK), dimension(:),POINTER :: poxi, poyi, pozi
      
      !!! MPI ROTATE           
      allocate(old_despl,source=this%despl)
      allocate(old_matriz,source=this%matriz)

      !X->Y->Z->X
      IF (MPIDIR==2 ) THEN
         OXI=old_matriz%totalX
         OYI=old_matriz%totalY
         OZI=old_matriz%totalZ
         !
         this%matriz%totalX=OZI
         this%matriz%totalY=OXI
         this%matriz%totalZ=OYI
         !
         OXI=old_despl%nX
         OYI=old_despl%nY
         OZI=old_despl%nZ
         !
         this%despl%nX=OZI
         this%despl%nY=OXI
         this%despl%nZ=OYI
         !
         OXI=old_despl%mX1
         OYI=old_despl%mY1
         OZI=old_despl%mZ1
         !
         this%despl%mX1=OZI
         this%despl%mY1=OXI
         this%despl%mZ1=OYI
         !
         OXI=old_despl%mX2
         OYI=old_despl%mY2
         OZI=old_despl%mZ2
         !
         this%despl%mX2=OZI
         this%despl%mY2=OXI
         this%despl%mZ2=OYI
         !
         rOXI=old_despl%originX
         rOYI=old_despl%originY
         rOZI=old_despl%originZ
         !
         this%despl%originX=rOZI
         this%despl%originY=rOXI
         this%despl%originZ=rOYI
         !         
         poxi => old_despl%desX 
         poyi => old_despl%desY 
         pozi => old_despl%desZ 
         !                  
         this%despl%desX => pozi
         this%despl%desY => poxi
         this%despl%desZ => poyi
      !X->Z->Y->X
      ELSEIF (MPIDIR==1 ) THEN
         OXI=old_matriz%totalX
         OYI=old_matriz%totalY
         OZI=old_matriz%totalZ
         !
         this%matriz%totalX=OYI
         this%matriz%totalY=OZI
         this%matriz%totalZ=OXI
         !
         OXI=old_despl%nX
         OYI=old_despl%nY
         OZI=old_despl%nZ
         !
         this%despl%nX=OYI
         this%despl%nY=OZI
         this%despl%nZ=OXI
         !
         OXI=old_despl%mX1
         OYI=old_despl%mY1
         OZI=old_despl%mZ1
         !
         this%despl%mX1=OYI
         this%despl%mY1=OZI
         this%despl%mZ1=OXI
         !
         OXI=old_despl%mX2
         OYI=old_despl%mY2
         OZI=old_despl%mZ2
         !
         this%despl%mX2=OYI
         this%despl%mY2=OZI
         this%despl%mZ2=OXI
         !
         rOXI=old_despl%originX
         rOYI=old_despl%originY
         rOZI=old_despl%originZ
         !
         this%despl%originX=rOYI
         this%despl%originY=rOZI
         this%despl%originZ=rOXI
         !
         poxi => old_despl%desX 
         poyi => old_despl%desY 
         pozi => old_despl%desZ 
         !
         this%despl%desX => poyi
         this%despl%desY => pozi
         this%despl%desZ => poxi

      ENDIF
      !!!!!!!!! fin rotacion
      deallocate (old_despl,old_matriz)
      RETURN
   END SUBROUTINE rotate_generateSpaceSteps
   
   SUBROUTINE rotate_generateCurrent_Field_Sources (this,mpidir)  
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama2,tama3,tama,i,ii
      
      tama = this%nodsrc%n_nodSrc   
      DO i = 1, tama
          tama2 = this%nodsrc%NodalSource(i)%n_c1P
          DO ii = 1, tama2
            CALL ROTATEMPI_SCALED(mpidir,this%nodsrc%NodalSource(i)%c1P(ii))
          end do           
          tama3 = this%nodsrc%NodalSource(i)%n_c2P
          DO ii = 1, tama3
            CALL ROTATEMPI_SCALED(mpidir,this%nodsrc%NodalSource(i)%c2P(ii))
          end do         
      end do         
      RETURN
   END SUBROUTINE rotate_generateCurrent_Field_Sources
  
   SUBROUTINE rotate_generatePlaneWaves (this, mpidir)
      TYPE (Parseador), INTENT (INOUT) :: this   
      TYPE (PlaneWaves),  pointer :: old_plnSrc => null ( )
      INTEGER (KIND=4) ::  mpidir
      integer (kind=4) :: oxi,oyi,ozi,oxe,oye,oze
      REAL (KIND=RK) :: theta,phi,alpha,beta     
      integer (kind=4) :: tama,i
      
      !!! MPI ROTATE         
      
      tama = (this%plnSrc%nc)     
      allocate(old_plnSrc,source=this%plnSrc)  
      do i=1,tama
      
          theta = old_plnSrc%collection(i)%theta 
          phi   = old_plnSrc%collection(i)%phi    
          alpha = old_plnSrc%collection(i)%alpha 
          beta  = old_plnSrc%collection(i)%beta 

          IF (MPIDIR==2 ) THEN
             OXI=  old_plnSrc%collection(i)%coor1 (1)
             OXE=  old_plnSrc%collection(i)%coor2 (1)
             OYI=  old_plnSrc%collection(i)%coor1 (2)
             OYE=  old_plnSrc%collection(i)%coor2 (2)
             OZI=  old_plnSrc%collection(i)%coor1 (3)
             OZE=  old_plnSrc%collection(i)%coor2 (3)
             !
             this%plnSrc%collection(i)%coor1 (1) =OZI
             this%plnSrc%collection(i)%coor2 (1) =OZE
             this%plnSrc%collection(i)%coor1 (2) =OXI
             this%plnSrc%collection(i)%coor2 (2) =OXE
             this%plnSrc%collection(i)%coor1 (3) =OYI
             this%plnSrc%collection(i)%coor2 (3) =OYE

             this%plnSrc%collection(i)%theta = atan2(Sqrt(Cos(theta)**2.0_RKIND+ Cos(phi)**2*Sin(theta)**2),Sin(phi)*Sin(theta))
             this%plnSrc%collection(i)%phi =   atan2(Cos(phi)*Sin(theta),Cos(theta))
             this%plnSrc%collection(i)%alpha = atan2(Sqrt(Cos(alpha)**2.0_RKIND+ Cos(beta)**2*Sin(alpha)**2),Sin(beta)*Sin(alpha))
             this%plnSrc%collection(i)%beta =  atan2(Cos(beta)*Sin(alpha),Cos(alpha))

          ELSEIF (MPIDIR==1 ) THEN
             OXI=  old_plnSrc%collection(i)%coor1 (1)
             OXE=  old_plnSrc%collection(i)%coor2 (1)
             OYI=  old_plnSrc%collection(i)%coor1 (2)
             OYE=  old_plnSrc%collection(i)%coor2 (2)
             OZI=  old_plnSrc%collection(i)%coor1 (3)
             OZE=  old_plnSrc%collection(i)%coor2 (3)
             !
             this%plnSrc%collection(i)%coor1 (1) =OYI
             this%plnSrc%collection(i)%coor2 (1) =OYE
             this%plnSrc%collection(i)%coor1 (2) =OZI
             this%plnSrc%collection(i)%coor2 (2) =OZE
             this%plnSrc%collection(i)%coor1 (3) =OXI
             this%plnSrc%collection(i)%coor2 (3) =OXE
                                  
             this%plnSrc%collection(i)%theta = atan2(Sqrt(Cos(theta)**2.0_RKIND+ Sin(phi)**2*Sin(theta)**2),Cos(phi)*Sin(theta))
             this%plnSrc%collection(i)%phi =   atan2(Cos(theta),Sin(phi)*Sin(theta))
             this%plnSrc%collection(i)%alpha = atan2(Sqrt(Cos(alpha)**2.0_RKIND+ Sin(beta)**2*Sin(alpha)**2),Cos(beta)*Sin(alpha))
             this%plnSrc%collection(i)%beta =  atan2(Cos(alpha),Sin(beta)*Sin(alpha))
          ENDIF
          !!!!!
          
        end do          
        deallocate(old_plnSrc)

      RETURN

   END SUBROUTINE rotate_generatePlaneWaves
   
   SUBROUTINE rotate_generateBoxSources (this, mpidir) 
      TYPE (Parseador), INTENT (INOUT) :: this   
      TYPE (Boxes),  pointer :: old_boxSrc => null ( )
      INTEGER (KIND=4) ::  mpidir
      integer (kind=4) :: oxi,oyi,ozi,oxe,oye,oze
      integer (kind=4) :: tama,i
      
         
      
      tama = (this%boxSrc%nvols)   
      allocate(old_boxSrc,source=this%boxSrc)                      
      do i=1,tama
      
          !MPI  ROTATE BOX
          IF (MPIDIR==2 ) THEN
             OXI=  old_boxSrc%vols(i)%coor1 (1)
             OXE=  old_boxSrc%vols(i)%coor2 (1)
             OYI=  old_boxSrc%vols(i)%coor1 (2)
             OYE=  old_boxSrc%vols(i)%coor2 (2)
             OZI=  old_boxSrc%vols(i)%coor1 (3)
             OZE=  old_boxSrc%vols(i)%coor2 (3)
             !
             this%boxSrc%vols(i)%coor1 (1) =OZI
             this%boxSrc%vols(i)%coor2 (1) =OZE
             this%boxSrc%vols(i)%coor1 (2) =OXI
             this%boxSrc%vols(i)%coor2 (2) =OXE
             this%boxSrc%vols(i)%coor1 (3) =OYI
             this%boxSrc%vols(i)%coor2 (3) =OYE
          ELSEIF (MPIDIR==1 ) THEN
             OXI=  old_boxSrc%vols(i)%coor1 (1)
             OXE=  old_boxSrc%vols(i)%coor2 (1)
             OYI=  old_boxSrc%vols(i)%coor1 (2)
             OYE=  old_boxSrc%vols(i)%coor2 (2)
             OZI=  old_boxSrc%vols(i)%coor1 (3)
             OZE=  old_boxSrc%vols(i)%coor2 (3)
             !
             this%boxSrc%vols(i)%coor1 (1) =OYI
             this%boxSrc%vols(i)%coor2 (1) =OYE
             this%boxSrc%vols(i)%coor1 (2) =OZI
             this%boxSrc%vols(i)%coor2 (2) =OZE
             this%boxSrc%vols(i)%coor1 (3) =OXI
             this%boxSrc%vols(i)%coor2 (3) =OXE
          ENDIF 
      end do      
     deallocate(old_boxSrc)
      !!!!!
      RETURN
   END SUBROUTINE rotate_generateBoxSources
   
   SUBROUTINE rotate_generateFronteras (this,mpidir)
      TYPE (Parseador), INTENT (INOUT) :: this   
      INTEGER (KIND=4) ::  mpidir       
      integer (kind=4) :: oxl,oxu,oyl,oyu,ozl,ozu  
      TYPE (FronteraPML) :: OPML_XL,OPML_XU,OPML_YL,OPML_YU,OPML_ZL,OPML_ZU
      
      !!! MPI ROTATE
      IF (MPIDIR==2 ) THEN
         OXL=this%front%tipofrontera(1)
         OXU=this%front%tipofrontera(2)
         OYL=this%front%tipofrontera(3)
         OYU=this%front%tipofrontera(4)
         OZL=this%front%tipofrontera(5)
         OZU=this%front%tipofrontera(6)
         !
         this%front%tipofrontera(1) = OZL
         this%front%tipofrontera(2) = OZU
         this%front%tipofrontera(3) = OXL
         this%front%tipofrontera(4) = OXU
         this%front%tipofrontera(5) = OYL
         this%front%tipofrontera(6) = OYU
         !
         OPML_XL%orden=this%front%propiedadesPML(1)%orden
         OPML_XU%orden=this%front%propiedadesPML(2)%orden
         OPML_YL%orden=this%front%propiedadesPML(3)%orden
         OPML_YU%orden=this%front%propiedadesPML(4)%orden
         OPML_ZL%orden=this%front%propiedadesPML(5)%orden
         OPML_ZU%orden=this%front%propiedadesPML(6)%orden
         !
         this%front%propiedadesPML(1)%orden = OPML_ZL%orden
         this%front%propiedadesPML(2)%orden = OPML_ZU%orden
         this%front%propiedadesPML(3)%orden = OPML_XL%orden
         this%front%propiedadesPML(4)%orden = OPML_XU%orden
         this%front%propiedadesPML(5)%orden = OPML_YL%orden
         this%front%propiedadesPML(6)%orden = OPML_YU%orden
         !
         !
         OPML_XL%refl=this%front%propiedadesPML(1)%refl
         OPML_XU%refl=this%front%propiedadesPML(2)%refl
         OPML_YL%refl=this%front%propiedadesPML(3)%refl
         OPML_YU%refl=this%front%propiedadesPML(4)%refl
         OPML_ZL%refl=this%front%propiedadesPML(5)%refl
         OPML_ZU%refl=this%front%propiedadesPML(6)%refl
         !
         this%front%propiedadesPML(1)%refl = OPML_ZL%refl
         this%front%propiedadesPML(2)%refl = OPML_ZU%refl
         this%front%propiedadesPML(3)%refl = OPML_XL%refl
         this%front%propiedadesPML(4)%refl = OPML_XU%refl
         this%front%propiedadesPML(5)%refl = OPML_YL%refl
         this%front%propiedadesPML(6)%refl = OPML_YU%refl
         !
         !
         OPML_XL%numCapas=this%front%propiedadesPML(1)%numCapas
         OPML_XU%numCapas=this%front%propiedadesPML(2)%numCapas
         OPML_YL%numCapas=this%front%propiedadesPML(3)%numCapas
         OPML_YU%numCapas=this%front%propiedadesPML(4)%numCapas
         OPML_ZL%numCapas=this%front%propiedadesPML(5)%numCapas
         OPML_ZU%numCapas=this%front%propiedadesPML(6)%numCapas
         !
         this%front%propiedadesPML(1)%numCapas = OPML_ZL%numCapas
         this%front%propiedadesPML(2)%numCapas = OPML_ZU%numCapas
         this%front%propiedadesPML(3)%numCapas = OPML_XL%numCapas
         this%front%propiedadesPML(4)%numCapas = OPML_XU%numCapas
         this%front%propiedadesPML(5)%numCapas = OPML_YL%numCapas
         this%front%propiedadesPML(6)%numCapas = OPML_YU%numCapas

      ELSEIF (MPIDIR==1 ) THEN
         OXL=this%front%tipofrontera(1)
         OXU=this%front%tipofrontera(2)
         OYL=this%front%tipofrontera(3)
         OYU=this%front%tipofrontera(4)
         OZL=this%front%tipofrontera(5)
         OZU=this%front%tipofrontera(6)
         !
         this%front%tipofrontera(1) = OYL
         this%front%tipofrontera(2) = OYU
         this%front%tipofrontera(3) = OZL
         this%front%tipofrontera(4) = OZU
         this%front%tipofrontera(5) = OXL
         this%front%tipofrontera(6) = OXU
         !
         OPML_XL%orden=this%front%propiedadesPML(1)%orden
         OPML_XU%orden=this%front%propiedadesPML(2)%orden
         OPML_YL%orden=this%front%propiedadesPML(3)%orden
         OPML_YU%orden=this%front%propiedadesPML(4)%orden
         OPML_ZL%orden=this%front%propiedadesPML(5)%orden
         OPML_ZU%orden=this%front%propiedadesPML(6)%orden
         !
         this%front%propiedadesPML(1)%orden = OPML_YL%orden
         this%front%propiedadesPML(2)%orden = OPML_YU%orden
         this%front%propiedadesPML(3)%orden = OPML_ZL%orden
         this%front%propiedadesPML(4)%orden = OPML_ZU%orden
         this%front%propiedadesPML(5)%orden = OPML_XL%orden
         this%front%propiedadesPML(6)%orden = OPML_XU%orden
         !
         OPML_XL%refl=this%front%propiedadesPML(1)%refl
         OPML_XU%refl=this%front%propiedadesPML(2)%refl
         OPML_YL%refl=this%front%propiedadesPML(3)%refl
         OPML_YU%refl=this%front%propiedadesPML(4)%refl
         OPML_ZL%refl=this%front%propiedadesPML(5)%refl
         OPML_ZU%refl=this%front%propiedadesPML(6)%refl
         !
         this%front%propiedadesPML(1)%refl = OPML_YL%refl
         this%front%propiedadesPML(2)%refl = OPML_YU%refl
         this%front%propiedadesPML(3)%refl = OPML_ZL%refl
         this%front%propiedadesPML(4)%refl = OPML_ZU%refl
         this%front%propiedadesPML(5)%refl = OPML_XL%refl
         this%front%propiedadesPML(6)%refl = OPML_XU%refl
         !
         OPML_XL%numCapas=this%front%propiedadesPML(1)%numCapas
         OPML_XU%numCapas=this%front%propiedadesPML(2)%numCapas
         OPML_YL%numCapas=this%front%propiedadesPML(3)%numCapas
         OPML_YU%numCapas=this%front%propiedadesPML(4)%numCapas
         OPML_ZL%numCapas=this%front%propiedadesPML(5)%numCapas
         OPML_ZU%numCapas=this%front%propiedadesPML(6)%numCapas
         !
         this%front%propiedadesPML(1)%numCapas = OPML_YL%numCapas
         this%front%propiedadesPML(2)%numCapas = OPML_YU%numCapas
         this%front%propiedadesPML(3)%numCapas = OPML_ZL%numCapas
         this%front%propiedadesPML(4)%numCapas = OPML_ZU%numCapas
         this%front%propiedadesPML(5)%numCapas = OPML_XL%numCapas
         this%front%propiedadesPML(6)%numCapas = OPML_XU%numCapas
      ENDIF
      !!!!!!!!! fin rotacion

      !END ROTATE FRONTERAS
      return
   END SUBROUTINE rotate_generateFronteras
   
   SUBROUTINE rotate_generatePECs (this,mpidir)   
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama,i
      
      tama = (this%pecregs%nvols)       
      DO i = 1, tama
        CALL ROTATEMPI(mpidir,this%pecRegs%Vols(i))    
      end do
      
      tama = (this%pecregs%nsurfs)    
      DO i = 1, tama
        CALL ROTATEMPI(mpidir,this%pecRegs%Surfs(i)) 
      end do
      
      tama = (this%pecregs%nlins)
      DO i = 1, tama
        CALL ROTATEMPI(mpidir,this%pecRegs%Lins(i))
      end do
      RETURN
   END SUBROUTINE rotate_generatePECs
   
   SUBROUTINE rotate_generatePMCs (this,mpidir)    
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama,i
      
      tama = (this%pmcregs%nvols)       
      DO i = 1, tama
        CALL ROTATEMPI(mpidir,this%pmcRegs%Vols(i))    
      end do
      
      tama = (this%pmcregs%nsurfs)    
      DO i = 1, tama
        CALL ROTATEMPI(mpidir,this%pmcRegs%Surfs(i)) 
      end do
      
      tama = (this%pmcregs%nlins)
      DO i = 1, tama
        CALL ROTATEMPI(mpidir,this%pmcRegs%Lins(i))
      end do
      RETURN
   END SUBROUTINE rotate_generatePMCs
   
   SUBROUTINE rotate_generateNONMetals (this,mpidir)       
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama2,tama3,tama,i,ii
                  
      !volumes
      tama = (this%DielRegs%nvols)
      DO i = 1, tama            
         tama2 = (this%DielRegs%vols(i)%n_c1P)
         DO ii = 1, tama2
            CALL ROTATEMPI(mpidir,this%DielRegs%vols(i)%C1P(ii)) 
         end do      
         this%DielRegs%vols(i)%DiodOrI = this%DielRegs%vols(i)%c1P(tama2)%Or    !UPDATE diodos POR SI HAY ROTACION ojo es un chapuz solo valido para diodos de 1p o 2p
         
         tama3 = (this%DielRegs%vols(i)%n_c2P)  
         DO ii = 1, tama3
            CALL ROTATEMPI(mpidir,this%DielRegs%vols(i)%C2P(ii)) 
         end do        
         this%DielRegs%vols(i)%DiodOrI = this%DielRegs%vols(i)%c2P(tama3)%Or   !UPDATE diodos POR SI HAY ROTACION. ojo es un chapuz solo valido para diodos de 1p o 2p
      end do
      !surfaces
      tama = (this%DielRegs%nsurfs)
      DO i = 1, tama            
         tama2 = (this%DielRegs%surfs(i)%n_c1P)
         DO ii = 1, tama2
            CALL ROTATEMPI(mpidir,this%DielRegs%surfs(i)%C1P(ii)) 
         end do      
         this%DielRegs%surfs(i)%DiodOrI = this%DielRegs%surfs(i)%c1P(tama2)%Or    !UPDATE diodos POR SI HAY ROTACION ojo es un chapuz solo valido para diodos de 1p o 2p
         
         tama3 = (this%DielRegs%surfs(i)%n_c2P)  
         DO ii = 1, tama3
            CALL ROTATEMPI(mpidir,this%DielRegs%surfs(i)%C2P(ii)) 
         end do        
         this%DielRegs%surfs(i)%DiodOrI = this%DielRegs%surfs(i)%c2P(tama3)%Or   !UPDATE diodos POR SI HAY ROTACION. ojo es un chapuz solo valido para diodos de 1p o 2p
      end do
      !lines
      tama = (this%DielRegs%nlins)
      DO i = 1, tama            
         tama2 = (this%DielRegs%lins(i)%n_c1P)
         DO ii = 1, tama2
            CALL ROTATEMPI(mpidir,this%DielRegs%lins(i)%C1P(ii)) 
         end do      
         this%DielRegs%lins(i)%DiodOrI = this%DielRegs%lins(i)%c1P(tama2)%Or    !UPDATE diodos POR SI HAY ROTACION ojo es un chapuz solo valido para diodos de 1p o 2p
         
         tama3 = (this%DielRegs%lins(i)%n_c2P)  
         DO ii = 1, tama3
            CALL ROTATEMPI(mpidir,this%DielRegs%lins(i)%C2P(ii)) 
         end do        
         this%DielRegs%lins(i)%DiodOrI = this%DielRegs%lins(i)%c2P(tama3)%Or   !UPDATE diodos POR SI HAY ROTACION. ojo es un chapuz solo valido para diodos de 1p o 2p
      end do
      RETURN
   END SUBROUTINE rotate_generateNONMetals
   
   SUBROUTINE rotate_generateANISOTROPICs (this,mpidir)         
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      if ((mpidir/=1).and.(this%ANIMATS%nvols+this%ANIMATS%nsurfs+this%ANIMATS%nlins/=0)) then
           print *,'Rotations in anisotropic unsupported'
           stop
           return
      endif
   !si algun dia lo hubiera es un cut and paste de rotate_generateNONMetals y a la que hay que aniadir la rotacion de la matriz de medios 
   RETURN
   END SUBROUTINE rotate_generateANISOTROPICs
   
   SUBROUTINE rotate_generateThinWires (this,mpidir)     
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      TYPE (ThinWires), POINTER :: old_tWires    
      integer (kind=4) :: tama,tama2,i,ii
      
      tama = this%twires%n_tw
      allocate(old_twires,source=this%twires)
      do i=1, tama       
         tama2 = this%twires%TW(i)%N_TWC
         DO ii = 1, tama2
      !!!ROTATE THINWIRE
             IF (MPIDIR==2 ) THEN
                   this%twires%tw(i)%tWc(ii)%i = old_twires%tw(i)%twc(ii)%K
                   this%twires%tw(i)%tWc(ii)%j = old_twires%tw(i)%twc(ii)%I
                   this%twires%tw(i)%tWc(ii)%K = old_twires%tw(i)%twc(ii)%J

                   SELECT CASE (old_twires%tw(i)%tWc(ii)%d)
                    CASE (iEx)
                      this%twires%tw(i)%tWc(ii)%d = iEy
                    CASE (iEY)
                      this%twires%tw(i)%tWc(ii)%d = iEz
                    CASE (iEZ)
                      this%twires%tw(i)%tWc(ii)%d = iEx
                   END SELECT
            ELSEIF (MPIDIR==1 ) THEN
                   this%twires%tw(i)%tWc(ii)%i = old_twires%tw(ii)%twc(i)%J
                   this%twires%tw(i)%tWc(ii)%j = old_twires%tw(ii)%twc(i)%K
                   this%twires%tw(i)%tWc(ii)%K = old_twires%tw(ii)%twc(i)%I
      
                   SELECT CASE (old_twires%tw(i)%tWc(ii)%d)
                    CASE (iEx)
                      this%twires%tw(i)%tWc(ii)%d = iEz
                    CASE (iEY)
                      this%twires%tw(i)%tWc(ii)%d = iEx
                    CASE (iEZ)
                      this%twires%tw(i)%tWc(ii)%d = iEy
                   END SELECT
             ENDIF
      !!!FIN  
         end do
      end do
      
      deallocate(old_twires)
     
      RETURN
   END SUBROUTINE rotate_generateThinWires
   
   SUBROUTINE rotate_generateSlantedWires (this,mpidir)      
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      TYPE (SlantedWires), POINTER :: old_swires    
      integer (kind=4) :: tama,tama2,i,ii
      
      tama = this%swires%n_sw
      allocate(old_swires,source=this%swires)
      do i=1, tama       
         tama2 = this%swires%sW(i)%N_SWC
         do ii=1,tama2
             !!!ROTATE THINWIRE
             IF (MPIDIR==2 ) THEN
                      this%swires%sw(i)%swc(ii)%x = old_swires%sw(i)%swc(ii)%z
                      this%swires%sw(i)%swc(ii)%y = old_swires%sw(i)%swc(ii)%x
                      this%swires%sw(i)%swc(ii)%z = old_swires%sw(i)%swc(ii)%y
             ELSEIF (MPIDIR==1 ) THEN
                      this%swires%sw(i)%swc(ii)%x = old_swires%sw(i)%swc(ii)%y
                      this%swires%sw(i)%swc(ii)%y = old_swires%sw(i)%swc(ii)%z
                      this%swires%sw(i)%swc(ii)%z = old_swires%sw(i)%swc(ii)%x
             ENDIF
         end do
         !!!FIN
      end do   
      allocate(old_swires)

      RETURN
   END SUBROUTINE rotate_generateSlantedWires
   
   SUBROUTINE rotate_generateThinSlots (this,mpidir)     
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      TYPE (ThinSlots), POINTER :: old_tSlots    
      integer (kind=4) :: tama,tama2,i,ii
      
      tama = this%tSlots%n_Tg
      allocate(old_tSlots,source=this%tSlots)
      do i=1, tama       
         tama2 = this%tSlots%Tg(i)%N_Tgc
         DO ii = 1, tama2
              !!!ROTATE THIN SLOT
              IF (MPIDIR==2 ) THEN
                       this%tSlots%Tg(i)%TgC(ii)%i = old_tsLots%Tg(i)%TgC(ii)%K
                       this%tSlots%Tg(i)%TgC(ii)%j = old_tsLots%Tg(i)%TgC(ii)%I
                       this%tSlots%Tg(i)%TgC(ii)%K = old_tsLots%Tg(i)%TgC(ii)%J
                       SELECT CASE (old_tsLots%Tg(i)%TgC(ii)%dir)
                        CASE (iEx)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEy
                        CASE (iEY)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEz
                        CASE (iEZ)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEx
                       END SELECT
              ELSEIF (MPIDIR==1 ) THEN
                       this%tSlots%Tg(i)%TgC(ii)%i = old_tsLots%Tg(i)%TgC(ii)%J
                       this%tSlots%Tg(i)%TgC(ii)%j = old_tsLots%Tg(i)%TgC(ii)%K
                       this%tSlots%Tg(i)%TgC(ii)%K = old_tsLots%Tg(i)%TgC(ii)%I

                       SELECT CASE (old_tsLots%Tg(i)%TgC(ii)%dir)
                        CASE (iEx)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEz
                        CASE (iEY)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEx
                        CASE (iEZ)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEy
                       END SELECT
              ENDIF
        end do
      end do      
      RETURN
      
   END SUBROUTINE rotate_generateThinSlots
!!   
   SUBROUTINE rotate_generateLossyThinSurface (this,mpidir)    
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama2,tama,i,ii
                                         
      tama = this%LossyThinSurfs%length
      DO i = 1, tama
         tama2 = this%LossyThinSurfs%cs(i)%nc
          DO ii = 1, tama2
             CALL ROTATEMPI(mpidir,this%LossyThinSurfs%cs(i)%C(ii)) 
          end do   
      end do    

      RETURN
   END SUBROUTINE rotate_generateLossyThinSurface

   SUBROUTINE rotate_generateFDMs (this,mpidir) 
      
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama,tama2,i,ii
      
      tama = (this%FRQDEPMATS%nvols)       
      DO i = 1, tama   
         tama2 = this%FRQDEPMATS%vols(i)%n_C
         DO ii = 1, tama2
            CALL ROTATEMPI(mpidir,this%FRQDEPMATS%Vols(i)%c(ii))  
         end do
      end do
      
      tama = (this%FRQDEPMATS%nsurfs)    
      DO i = 1, tama 
         tama2 = this%FRQDEPMATS%surfs(i)%n_C
         DO ii = 1, tama2
            CALL ROTATEMPI(mpidir,this%FRQDEPMATS%Surfs(i)%c(ii))   
         end do
      end do
      
      tama = (this%FRQDEPMATS%nlins)
      DO i = 1, tama 
         tama2 = this%FRQDEPMATS%Lins(i)%n_C
         DO ii = 1, tama2
            CALL ROTATEMPI(mpidir,this%FRQDEPMATS%Lins(i)%c(ii))   
         end do
      end do
      RETURN
      
   END SUBROUTINE rotate_generateFDMs
   
   SUBROUTINE rotate_generateSONDAs (this,mpidir) 
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir    
      integer (kind=4) :: tama,tama2,tama3,i,ii,iii  
      TYPE (FarField_Sonda), POINTER :: old_FarField => NULL ()
      TYPE (Electric_Sonda), POINTER :: old_Electric => NULL ()
      TYPE (Magnetic_Sonda), POINTER :: old_Magnetic => NULL ()
      REAL(KIND=RK) :: THETASTART,THETASTOP,PHISTART,PHISTOP
      
      tama = this%oldSONDA%n_probes        
      ! tres posibilidades FarField, Electric,Magnetic
      DO i = 1, tama      
         tama2 = (this%oldSONDA%probes(i)%n_FarField)  
         DO ii = 1, tama2                                 
            allocate (old_FarField,source=this%oldSONDA%probes(i)%FarField(ii))
            thetastart=old_FarField%probe%thetastart
            thetastop =old_FarField%probe%thetastop   
            phistart  =old_FarField%probe%phistart
            phistop   =old_FarField%probe%phistop
            !!!mpirotate angulos farfield .... las coordenadas se rotan luego
            IF (MPIDIR==2 ) THEN
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastart = atan2(Sqrt(Cos(thetastart)**2.0_RKIND+ Cos(phistart)**2*Sin(thetastart)**2),Sin(phistart)*Sin(thetastart))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastart   = atan2(Cos(phistart)*Sin(thetastart),Cos(thetastart))      
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastop = atan2(Sqrt(Cos(thetastop)**2.0_RKIND+ Cos(phistop)**2*Sin(thetastop)**2),Sin(phistop)*Sin(thetastop))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%phistop   = atan2(Cos(phistop)*Sin(thetastop),Cos(thetastop))
            ELSEIF (MPIDIR==1 ) THEN
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastart = atan2(Sqrt(Cos(thetastart)**2.0_RKIND+ Sin(phistart)**2*Sin(thetastart)**2),Cos(phistart)*Sin(thetastart))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastart   = atan2(Cos(thetastart),Sin(phistart)*Sin(thetastart))    
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastop = atan2(Sqrt(Cos(thetastop)**2.0_RKIND+ Sin(phistop)**2*Sin(thetastop)**2),Cos(phistop)*Sin(thetastop))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%phistop   = atan2(Cos(thetastop),Sin(phistop)*Sin(thetastop))
            ENDIF        
            tama3 = (this%oldSONDA%probes(i)%FarField(ii)%probe%n_cord)    
            DO iii = 1, tama3
              !!!ROTATE MPI
              IF (MPIDIR==2 ) THEN
                 this%oldSONDA%probes(i)%FarField(ii)%probe%i(iii) = old_FarField%probe%K(iii)
                 this%oldSONDA%probes(i)%FarField(ii)%probe%j(iii) = old_FarField%probe%I(iii)
                 this%oldSONDA%probes(i)%FarField(ii)%probe%K(iii) = old_FarField%probe%J(iii)
              ELSEIF (MPIDIR==1 ) THEN                    
                 this%oldSONDA%probes(i)%FarField(ii)%probe%i(iii) = old_FarField%probe%J(iii)
                 this%oldSONDA%probes(i)%FarField(ii)%probe%j(iii) = old_FarField%probe%K(iii)
                 this%oldSONDA%probes(i)%FarField(ii)%probe%K(iii) = old_FarField%probe%I(iii)
              ENDIF
            end do             
            deallocate (old_FarField)
         end do
      end do
      !     
      DO i = 1, tama      
         tama2 = (this%oldSONDA%probes(i)%n_Electric)  
         DO ii = 1, tama2                                 
            allocate (old_Electric,source=this%oldSONDA%probes(i)%Electric(ii))        
            tama3 = (this%oldSONDA%probes(i)%Electric(ii)%probe%n_cord)    
            DO iii = 1, tama3
              !!!ROTATE MPI
              IF (MPIDIR==2 ) THEN
                 this%oldSONDA%probes(i)%Electric(ii)%probe%i(iii) = old_Electric%probe%K(iii)
                 this%oldSONDA%probes(i)%Electric(ii)%probe%j(iii) = old_Electric%probe%I(iii)
                 this%oldSONDA%probes(i)%Electric(ii)%probe%K(iii) = old_Electric%probe%J(iii)
              ELSEIF (MPIDIR==1 ) THEN                    
                 this%oldSONDA%probes(i)%Electric(ii)%probe%i(iii) = old_Electric%probe%J(iii)
                 this%oldSONDA%probes(i)%Electric(ii)%probe%j(iii) = old_Electric%probe%K(iii)
                 this%oldSONDA%probes(i)%Electric(ii)%probe%K(iii) = old_Electric%probe%I(iii)
              ENDIF
            end do             
            deallocate (old_Electric)
         end do
       end do
      !
      DO i = 1, tama      
         tama2 = (this%oldSONDA%probes(i)%n_Magnetic)  
         DO ii = 1, tama2                                 
            allocate (old_Magnetic,source=this%oldSONDA%probes(i)%Magnetic(ii))        
            tama3 = (this%oldSONDA%probes(i)%Magnetic(ii)%probe%n_cord)    
            DO iii = 1, tama3
              !!!ROTATE MPI
              IF (MPIDIR==2 ) THEN
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%i(iii) = old_Magnetic%probe%K(iii)
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%j(iii) = old_Magnetic%probe%I(iii)
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%K(iii) = old_Magnetic%probe%J(iii)
              ELSEIF (MPIDIR==1 ) THEN                    
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%i(iii) = old_Magnetic%probe%J(iii)
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%j(iii) = old_Magnetic%probe%K(iii)
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%K(iii) = old_Magnetic%probe%I(iii)
              ENDIF
            end do             
            deallocate (old_Magnetic)
         end do
       end do
      !
       return
   END SUBROUTINE rotate_generateSONDAs
   
   SUBROUTINE rotate_generateMasSondas (this,mpidir)     
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir                          
      integer (kind=4) :: tama,tama2,i,ii  
      integer (kind=4) :: oxi,oyi,ozi,oxe,oye,oze,oor,TXI,TYI,TZI  
      TYPE (Coords), POINTER :: old_MasSonda => NULL ()
      
      tama = this%Sonda%length       
      ! tres posibilidades FarField, Electric,Magnetic
      DO i = 1, tama      
         tama2 = (this%Sonda%collection(i)%len_cor)    
         DO ii = 1, tama2                                               
              allocate (old_MasSonda,source=this%Sonda%collection(i)%cordinates(ii))                    
              OXI=old_MasSonda%XI
              OXE=old_MasSonda%XE
              OYI=old_MasSonda%YI
              OYE=old_MasSonda%YE
              OZI=old_MasSonda%ZI
              OZE=old_MasSonda%ZE
              OOR=old_MasSonda%OR          
              TXI=old_MasSonda%Xtrancos      
              TYI=old_MasSonda%Ytrancos
              TZI=old_MasSonda%Ztrancos
              IF ((OOR/=NP_COR_EX).AND.(OOR/=NP_COR_EY).AND.(OOR/=NP_COR_EZ).AND. &
              (OOR/=NP_COR_HX).AND.(OOR/=NP_COR_HY).AND.(OOR/=NP_COR_HZ)) RETURN
              !!LAS IW Y LAS VG NO SE ROTAN
              IF (MPIDIR==2 ) THEN
                 this%Sonda%collection(i)%cordinates(ii)%XI=OZI   
                 this%Sonda%collection(i)%cordinates(ii)%XE=OZE
                 this%Sonda%collection(i)%cordinates(ii)%Xtrancos=TZI
                 
                 this%Sonda%collection(i)%cordinates(ii)%YI=OXI 
                 this%Sonda%collection(i)%cordinates(ii)%YE=OXE
                 this%Sonda%collection(i)%cordinates(ii)%Ytrancos=TXI
                 
                 this%Sonda%collection(i)%cordinates(ii)%ZI=OYI 
                 this%Sonda%collection(i)%cordinates(ii)%ZE=OYE
                 this%Sonda%collection(i)%cordinates(ii)%Ztrancos=TYI
                 
                 IF (OOR== NP_COR_EX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EY
                 IF (OOR== NP_COR_EY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EZ
                 IF (OOR== NP_COR_EZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EX
                           
                 IF (OOR== NP_COR_hX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_HY
                 IF (OOR== NP_COR_hY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_HZ
                 IF (OOR== NP_COR_hZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_HX
              ELSEIF (MPIDIR==1 ) THEN
                 this%Sonda%collection(i)%cordinates(ii)%XI=OYI
                 this%Sonda%collection(i)%cordinates(ii)%XE=OYE 
                 this%Sonda%collection(i)%cordinates(ii)%Xtrancos= TYI
                 
                 this%Sonda%collection(i)%cordinates(ii)%YI=OZI
                 this%Sonda%collection(i)%cordinates(ii)%YE=OZE               
                 this%Sonda%collection(i)%cordinates(ii)%Ytrancos=TZI
                 
                 this%Sonda%collection(i)%cordinates(ii)%ZI=OXI
                 this%Sonda%collection(i)%cordinates(ii)%ZE=OXE 
                 this%Sonda%collection(i)%cordinates(ii)%Ztrancos=TXI
                 
                 IF (OOR== NP_COR_EX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EZ
                 IF (OOR== NP_COR_EY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EX
                 IF (OOR== NP_COR_EZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EY
                           
                 IF (OOR== NP_COR_HX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_hZ
                 IF (OOR== NP_COR_HY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_hX
                 IF (OOR== NP_COR_HZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_hY
              ENDIF                                              
              deallocate (old_MasSonda)
         end do
      end do
      RETURN
   END SUBROUTINE rotate_generateMasSondas
   
   SUBROUTINE rotate_generateBloqueProbes (this,mpidir)    
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir                          
      integer (kind=4) :: tama,i  
      integer (kind=4) :: oxi,oyi,ozi,oxe,oye,oze  
      TYPE (BloqueProbe), POINTER :: old_BloqueProbe => NULL ()
      
      tama = this%BloquePRB%N_BP   
      DO i = 1, tama      
          allocate(old_BloqueProbe,source=this%BloquePRB%BP(i))
          !MPI  ROTATE Bloque CURRENT
          IF (MPIDIR==2 ) THEN
             OXI=  old_BloqueProbe%i1
             OXE=  old_BloqueProbe%i2
             OYI=  old_BloqueProbe%j1
             OYE=  old_BloqueProbe%j2
             OZI=  old_BloqueProbe%k1
             OZE=  old_BloqueProbe%k2
             !
             this%BloquePRB%BP(i)%i1 =OZI
             this%BloquePRB%BP(i)%i2 =OZE
             this%BloquePRB%BP(i)%j1 =OXI
             this%BloquePRB%BP(i)%j2 =OXE
             this%BloquePRB%BP(i)%k1 =OYI
             this%BloquePRB%BP(i)%k2 =OYE
             SELECT CASE (this%BloquePRB%BP(i)%nml)
              CASE (iEx)
                this%BloquePRB%BP(i)%nml =  iEy
              CASE (iEy)
                this%BloquePRB%BP(i)%nml =  iEz
              CASE (iEz)
                this%BloquePRB%BP(i)%nml =  iEx
              CASE DEFAULT
             END SELECT
          ELSEIF (MPIDIR==1 ) THEN
             OXI=  old_BloqueProbe%i1
             OXE=  old_BloqueProbe%i2
             OYI=  old_BloqueProbe%j1
             OYE=  old_BloqueProbe%j2
             OZI=  old_BloqueProbe%k1
             OZE=  old_BloqueProbe%k2
             !
             this%BloquePRB%BP(i)%i1 =OYI
             this%BloquePRB%BP(i)%i2 =OYE
             this%BloquePRB%BP(i)%j1 =OZI
             this%BloquePRB%BP(i)%j2 =OZE
             this%BloquePRB%BP(i)%k1 =OXI
             this%BloquePRB%BP(i)%k2 =OXE
             SELECT CASE (this%BloquePRB%BP(i)%nml)
              CASE (iEx)
                this%BloquePRB%BP(i)%nml =  iEz
              CASE (iEy)
                this%BloquePRB%BP(i)%nml =  iEx
              CASE (iEz)
                this%BloquePRB%BP(i)%nml =  iEy
              CASE DEFAULT
             END SELECT
          ENDIF           
          deallocate(old_BloqueProbe)
      end do
      
          
      !!!!FIN ROTATE
      RETURN
   END SUBROUTINE rotate_generateBloqueProbes
!!   
   SUBROUTINE rotate_generateVolumicProbes(this,mpidir)    
      TYPE (Parseador), INTENT (INOUT) :: this
      INTEGER (KIND=4) ::  mpidir                          
      integer (kind=4) :: tama,tama2,i,ii  
      integer (kind=4) :: oxi,oyi,ozi,oxe,oye,oze,oor,TXI,TYI,TZI  
      TYPE (Coords), POINTER :: old_Coordinates => NULL ()
      
      tama = this%Sonda%length       
      ! tres posibilidades FarField, Electric,Magnetic
      DO i = 1, tama      
         tama2 = (this%Sonda%collection(i)%len_cor)    
         DO ii = 1, tama2                                               
              allocate (old_Coordinates,source=this%VolPrb%collection(i)%cordinates(ii))    
              OXI=old_Coordinates%XI
              OXE=old_Coordinates%XE        
              TXI=old_Coordinates%Xtrancos  
              OYI=old_Coordinates%YI
              OYE=old_Coordinates%YE     
              TYI=old_Coordinates%Ytrancos
              OZI=old_Coordinates%ZI
              OZE=old_Coordinates%ZE  
              TZI=old_Coordinates%Ztrancos    
              OOR=old_Coordinates%OR    

              OXI=old_Coordinates%XI
              OXE=old_Coordinates%XE
              OYI=old_Coordinates%YI
              OYE=old_Coordinates%YE
              OZI=old_Coordinates%ZI
              OZE=old_Coordinates%ZE
              OOR=old_Coordinates%OR
              !TRANCOS
              TXI=old_Coordinates%XTRANCOS
              TYI=old_Coordinates%YTRANCOS
              TZI=old_Coordinates%ZTRANCOS
              !
              IF ((OOR/=iExC).AND.(OOR/=iEyC).AND.(OOR/=iEzC).AND. &
              (OOR/=iHxC).AND.(OOR/=iHyC).AND.(OOR/=iHzC).AND. &
              (OOR/=iCurX).AND.(OOR/=iCurY).AND.(OOR/=iCurZ).AND. &
               (OOR/=iMEC).AND.(OOR/=iMHC).AND.(OOR/=iCur)) RETURN
              !!LAS IW Y LAS VG NO SE ROTAN.
        !!las imec, imhc e icur no le afecta el oor
              IF (MPIDIR==2 ) THEN
                 this%VolPrb%collection(i)%cordinates(ii)%XI=OZI
                 this%VolPrb%collection(i)%cordinates(ii)%XE=OZE
                 this%VolPrb%collection(i)%cordinates(ii)%YI=OXI
                 this%VolPrb%collection(i)%cordinates(ii)%YE=OXE
                 this%VolPrb%collection(i)%cordinates(ii)%ZI=OYI
                 this%VolPrb%collection(i)%cordinates(ii)%ZE=OYE
        !    
                 this%VolPrb%collection(i)%cordinates(ii)%XTRANCOS=TZI
                 this%VolPrb%collection(i)%cordinates(ii)%YTRANCOS=TXI
                 this%VolPrb%collection(i)%cordinates(ii)%ZTRANCOS=TYI
                 !
                 IF (OOR== iEXc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iEYc
                 IF (OOR== iEYc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iEZc
                 IF (OOR== iEZc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iEXc
                     
                 IF (OOR== ihXc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iHYc
                 IF (OOR== ihYc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iHZc
                 IF (OOR== ihZc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iHXc
                           
                 IF (OOR== iCurX)  this%VolPrb%collection(i)%cordinates(ii)%OR= iCurY
                 IF (OOR== iCurY)  this%VolPrb%collection(i)%cordinates(ii)%OR= iCurZ
                 IF (OOR== iCurZ)  this%VolPrb%collection(i)%cordinates(ii)%OR= iCurX
              ELSEIF (MPIDIR==1 ) THEN
                 this%VolPrb%collection(i)%cordinates(ii)%XI=OYI
                 this%VolPrb%collection(i)%cordinates(ii)%XE=OYE
                 this%VolPrb%collection(i)%cordinates(ii)%YI=OZI
                 this%VolPrb%collection(i)%cordinates(ii)%YE=OZE
                 this%VolPrb%collection(i)%cordinates(ii)%ZI=OXI
                 this%VolPrb%collection(i)%cordinates(ii)%ZE=OXE
        !    
                 this%VolPrb%collection(i)%cordinates(ii)%XTRANCOS=TYI
                 this%VolPrb%collection(i)%cordinates(ii)%YTRANCOS=TZI
                 this%VolPrb%collection(i)%cordinates(ii)%ZTRANCOS=TXI
                 !
                 IF (OOR== iEXc)  this%VolPrb%collection(i)%cordinates(ii)%OR= iEZc
                 IF (OOR== iEYc)  this%VolPrb%collection(i)%cordinates(ii)%OR= iEXc
                 IF (OOR== iEZc)  this%VolPrb%collection(i)%cordinates(ii)%OR= iEYc
                    
                 IF (OOR== iHXc)  this%VolPrb%collection(i)%cordinates(ii)%OR= ihZc
                 IF (OOR== iHYc)  this%VolPrb%collection(i)%cordinates(ii)%OR= ihXc
                 IF (OOR== iHZc)  this%VolPrb%collection(i)%cordinates(ii)%OR= ihYc
                         
                 IF (OOR== iCurX) this%VolPrb%collection(i)%cordinates(ii)%OR= iCurZ
                 IF (OOR== iCurY) this%VolPrb%collection(i)%cordinates(ii)%OR= iCurX
                 IF (OOR== iCurZ) this%VolPrb%collection(i)%cordinates(ii)%OR= iCurY
              ENDIF                                                  
              deallocate (old_Coordinates)
         end do
      end do
      
      RETURN
   END SUBROUTINE rotate_generateVolumicProbes

!!                                                    
!!!---------------------------------------------------->
!!!---------------------------------------------------->
!!!---------------------------------------------------->
!!
   SUBROUTINE ROTATEMPI(mpidir,COORDEN)         
      INTEGER (KIND=4) ::  mpidir    
      TYPE (coords), INTENT (INOUT) :: COORDEN
      INTEGER (KIND=4)   :: OXI, OXE, OYI, OYE,OZI, OZE, OOR,TXI,TYI,TZI
      OXI=COORDEN%XI
      OXE=COORDEN%XE        
      TXI=COORDEN%Xtrancos  
      OYI=COORDEN%YI
      OYE=COORDEN%YE     
      TYI=COORDEN%Ytrancos
      OZI=COORDEN%ZI
      OZE=COORDEN%ZE  
      TZI=COORDEN%Ztrancos    
      OOR=COORDEN%OR    
      IF (MPIDIR==2 ) THEN
         COORDEN%XI=OZI
         COORDEN%XE=OZE      
         COORDEN%Xtrancos =TZI    
         COORDEN%YI=OXI
         COORDEN%YE=OXE  
         COORDEN%Ytrancos =TXI 
         COORDEN%ZI=OYI
         COORDEN%ZE=OYE     
         COORDEN%Ztrancos =TYI 
         IF (OOR== iEx) COORDEN%OR= iEy
         IF (OOR==-iEx) COORDEN%OR=-iEy
         IF (OOR== iEy) COORDEN%OR= iEz
         IF (OOR==-iEy) COORDEN%OR=-iEz
         IF (OOR== iEz) COORDEN%OR= iEx
         IF (OOR==-iEz) COORDEN%OR=-iEx
      ELSEIF (MPIDIR==1 ) THEN
         COORDEN%XI=OYI
         COORDEN%XE=OYE         
         COORDEN%Xtrancos =TYI    
         COORDEN%YI=OZI
         COORDEN%YE=OZE   
         COORDEN%Ytrancos =TZI 
         COORDEN%ZI=OXI
         COORDEN%ZE=OXE 
         COORDEN%Ztrancos =TXI 
         IF (OOR== iEx) COORDEN%OR= iEz
         IF (OOR==-iEx) COORDEN%OR=-iEz
         IF (OOR== iEy) COORDEN%OR= iEx
         IF (OOR==-iEy) COORDEN%OR=-iEx
         IF (OOR== iEz) COORDEN%OR= iEy
         IF (OOR==-iEz) COORDEN%OR=-iEy
      ENDIF
      RETURN
   END SUBROUTINE ROTATEMPI

   SUBROUTINE ROTATEMPI_SCALED(mpidir,COORDEN)        
      INTEGER (KIND=4) ::  mpidir    
      TYPE (coords_SCALED), INTENT (INOUT) :: COORDEN
      INTEGER (KIND=4)   :: OXI, OXE, OYI, OYE,OZI, OZE,oor
      REAL (KIND=RK) :: OXC,OYC,OZC
      OXI=COORDEN%XI
      OXE=COORDEN%XE
      OYI=COORDEN%YI
      OYE=COORDEN%YE
      OZI=COORDEN%ZI
      OZE=COORDEN%ZE
      OXC=COORDEN%XC
      OYC=COORDEN%YC  
      OZC=COORDEN%ZC
      OOr=COORDEN%or
      IF (MPIDIR==2 ) THEN
         COORDEN%XI=OZI
         COORDEN%XE=OZE
         COORDEN%YI=OXI
         COORDEN%YE=OXE
         COORDEN%ZI=OYI
         COORDEN%ZE=OYE
         COORDEN%XC=OzC
         COORDEN%YC=OxC
         COORDEN%ZC=OyC  
         IF (OOR== iEx) COORDEN%OR= iEy
         IF (OOR==-iEx) COORDEN%OR=-iEy
         IF (OOR== iEy) COORDEN%OR= iEz
         IF (OOR==-iEy) COORDEN%OR=-iEz
         IF (OOR== iEz) COORDEN%OR= iEx
         IF (OOR==-iEz) COORDEN%OR=-iEx

      ELSEIF (MPIDIR==1 ) THEN
         COORDEN%XI=OYI
         COORDEN%XE=OYE
         COORDEN%YI=OZI
         COORDEN%YE=OZE
         COORDEN%ZI=OXI
         COORDEN%ZE=OXE
         COORDEN%XC=OyC
         COORDEN%YC=OzC
         COORDEN%ZC=OxC         
         IF (OOR== iEx) COORDEN%OR= iEz
         IF (OOR==-iEx) COORDEN%OR=-iEz
         IF (OOR== iEy) COORDEN%OR= iEx
         IF (OOR==-iEy) COORDEN%OR=-iEx
         IF (OOR== iEz) COORDEN%OR= iEy
         IF (OOR==-iEz) COORDEN%OR=-iEy

      ENDIF
      RETURN
   END SUBROUTINE ROTATEMPI_SCALED
     

!!   

END MODULE nfde_rotate_m