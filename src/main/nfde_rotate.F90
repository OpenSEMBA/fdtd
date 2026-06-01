module nfde_rotate_m
    !       
   use NFDETypes_m

   !
   implicit none
   !
   public

contains
    
   subroutine nfde_rotate (this,mpidir) 
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir
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
      call  rotate_generateSONDAs                     (this, mpidir)
      call  rotate_generateMasSondas                  (this, mpidir)
      call  rotate_generateBloqueProbes               (this, mpidir)
      call  rotate_generateVolumicProbes              (this, mpidir)
#ifdef CompileWithMTLN
      call  rotate_mtln                               (this, mpidir)
#endif

      return
   end subroutine nfde_rotate

#ifdef CompileWithMTLN
   subroutine rotate_mtln(this, mpidir)
      type(Parseador_t), intent(inout) :: this          
      integer(kind=4) :: mpidir
      type(mtln_t), pointer :: old_mtln => null()
      integer(kind=4) :: i, j
      integer(kind=4) :: x, y, z, or

      allocate(old_mtln, source = this%mtln)
      do i = 1, size(old_mtln%cables)
         do j = 1, size(old_mtln%cables(i)%ptr%segments)
            x = old_mtln%cables(i)%ptr%segments(j)%x
            y = old_mtln%cables(i)%ptr%segments(j)%y
            z = old_mtln%cables(i)%ptr%segments(j)%z
            or = old_mtln%cables(i)%ptr%segments(j)%orientation
            if (mpidir == 2) then 
               this%mtln%cables(i)%ptr%segments(j)%x = z
               this%mtln%cables(i)%ptr%segments(j)%y = x
               this%mtln%cables(i)%ptr%segments(j)%z = y
               select case(abs(or))
               case(1)
                  this%mtln%cables(i)%ptr%segments(j)%orientation = sign(2, or)
               case(2)
                  this%mtln%cables(i)%ptr%segments(j)%orientation = sign(3, or)
               case(3)
                  this%mtln%cables(i)%ptr%segments(j)%orientation = sign(1, or)
               end select
            else if (mpidir == 1) then 
               this%mtln%cables(i)%ptr%segments(j)%x = y
               this%mtln%cables(i)%ptr%segments(j)%y = z
               this%mtln%cables(i)%ptr%segments(j)%z = x
               select case(abs(old_mtln%cables(i)%ptr%segments(j)%orientation))
               case(1)
                  this%mtln%cables(i)%ptr%segments(j)%orientation = sign(3, or)
               case(2)
                  this%mtln%cables(i)%ptr%segments(j)%orientation = sign(1, or)
               case(3)
                  this%mtln%cables(i)%ptr%segments(j)%orientation = sign(2, or)
               end select
            end if
         end do
      end do
      deallocate(old_mtln)
   end subroutine
#endif


   subroutine rotate_generateSpaceSteps (this, mpidir)
      type(Parseador_t), intent(inout) :: this          
      integer(kind=4) :: mpidir
      type(Desplazamiento_t), pointer :: old_despl => NULL ()      
      type(MatrizMedios_t), pointer :: old_matriz => NULL ()
      integer(kind=4) :: oxi,oyi,ozi
      real(kind=RK) :: roxi,royi,rozi
      real(kind=RK), dimension(:),pointer :: poxi, poyi, pozi
      
      !!! MPI ROTATE           
      allocate(old_despl,source=this%despl)
      allocate(old_matriz,source=this%matriz)

      !X->Y->Z->X
      if (MPIDIR==2 ) then
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
      ELSEIF (MPIDIR==1 ) then
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

      end if
      !!!!!!!!! fin rotacion
      deallocate(old_despl,old_matriz)
      return
   end subroutine rotate_generateSpaceSteps
   
   subroutine rotate_generateCurrent_Field_Sources (this,mpidir)  
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama2,tama3,tama,i,ii
      
      tama = this%nodsrc%n_nodSrc   
      do i = 1, tama
          tama2 = this%nodsrc%NodalSource(i)%n_c1P
          do ii = 1, tama2
            call ROTATEMPI_SCALED(mpidir,this%nodsrc%NodalSource(i)%c1P(ii))
          end do           
          tama3 = this%nodsrc%NodalSource(i)%n_c2P
          do ii = 1, tama3
            call ROTATEMPI_SCALED(mpidir,this%nodsrc%NodalSource(i)%c2P(ii))
          end do         
      end do         
      return
   end subroutine rotate_generateCurrent_Field_Sources
  
   subroutine rotate_generatePlaneWaves (this, mpidir)
      type(Parseador_t), intent(inout) :: this   
      type(PlaneWaves_t),  pointer :: old_plnSrc => null ( )
      integer(kind=4) :: mpidir
      integer(kind=4) :: oxi,oyi,ozi,oxe,oye,oze
      real(kind=RK) :: theta,phi,alpha,beta     
      integer(kind=4) :: tama,i
      
      !!! MPI ROTATE         
      
      tama = (this%plnSrc%nc)     
      allocate(old_plnSrc,source=this%plnSrc)  
      do i=1,tama
      
          theta = old_plnSrc%collection(i)%theta 
          phi   = old_plnSrc%collection(i)%phi    
          alpha = old_plnSrc%collection(i)%alpha 
          beta  = old_plnSrc%collection(i)%beta 

          if (MPIDIR==2 ) then
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

          ELSEIF (MPIDIR==1 ) then
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
          end if
          !!!!!
          
        end do          
        deallocate(old_plnSrc)

      return

   end subroutine rotate_generatePlaneWaves
   
   subroutine rotate_generateBoxSources (this, mpidir) 
      type(Parseador_t), intent(inout) :: this   
      type(Boxes_t),  pointer :: old_boxSrc => null ( )
      integer(kind=4) :: mpidir
      integer(kind=4) :: oxi,oyi,ozi,oxe,oye,oze
      integer(kind=4) :: tama,i
      
         
      
      tama = (this%boxSrc%nvols)   
      allocate(old_boxSrc,source=this%boxSrc)                      
      do i=1,tama
      
          !MPI  ROTATE BOX
          if (MPIDIR==2 ) then
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
          ELSEIF (MPIDIR==1 ) then
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
          end if 
      end do      
     deallocate(old_boxSrc)
      !!!!!
      return
   end subroutine rotate_generateBoxSources
   
   subroutine rotate_generateFronteras (this,mpidir)
      type(Parseador_t), intent(inout) :: this   
      integer(kind=4) :: mpidir       
      integer(kind=4) :: oxl,oxu,oyl,oyu,ozl,ozu  
      type(FronteraPML_t) :: OPML_XL,OPML_XU,OPML_YL,OPML_YU,OPML_ZL,OPML_ZU
      
      !!! MPI ROTATE
      if (MPIDIR==2 ) then
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

      ELSEIF (MPIDIR==1 ) then
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
      end if
      !!!!!!!!! fin rotacion

      !END ROTATE FRONTERAS
      return
   end subroutine rotate_generateFronteras
   
   subroutine rotate_generatePECs (this,mpidir)   
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama,i
      
      tama = (this%pecregs%nvols)       
      do i = 1, tama
        call ROTATEMPI(mpidir,this%pecRegs%Vols(i))    
      end do
      
      tama = (this%pecregs%nsurfs)    
      do i = 1, tama
        call ROTATEMPI(mpidir,this%pecRegs%Surfs(i)) 
      end do
      
      tama = (this%pecregs%nlins)
      do i = 1, tama
        call ROTATEMPI(mpidir,this%pecRegs%Lins(i))
      end do
      return
   end subroutine rotate_generatePECs
   
   subroutine rotate_generatePMCs (this,mpidir)    
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama,i
      
      tama = (this%pmcregs%nvols)       
      do i = 1, tama
        call ROTATEMPI(mpidir,this%pmcRegs%Vols(i))    
      end do
      
      tama = (this%pmcregs%nsurfs)    
      do i = 1, tama
        call ROTATEMPI(mpidir,this%pmcRegs%Surfs(i)) 
      end do
      
      tama = (this%pmcregs%nlins)
      do i = 1, tama
        call ROTATEMPI(mpidir,this%pmcRegs%Lins(i))
      end do
      return
   end subroutine rotate_generatePMCs
   
   subroutine rotate_generateNONMetals (this,mpidir)       
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama2,tama3,tama,i,ii
                  
      !volumes
      tama = (this%DielRegs%nvols)
      do i = 1, tama            
         tama2 = (this%DielRegs%vols(i)%n_c1P)
         do ii = 1, tama2
            call ROTATEMPI(mpidir,this%DielRegs%vols(i)%C1P(ii)) 
         end do
         if (tama2 > 0) then
            this%DielRegs%vols(i)%DiodOrI = this%DielRegs%vols(i)%c1P(tama2)%Or    !UPDATE diodos POR SI HAY ROTACION ojo es un chapuz solo valido para diodos de 1p o 2p
         end if
         tama3 = (this%DielRegs%vols(i)%n_c2P)  
         do ii = 1, tama3
            call ROTATEMPI(mpidir,this%DielRegs%vols(i)%C2P(ii)) 
         end do
         if (tama3 > 0) then 
            this%DielRegs%vols(i)%DiodOrI = this%DielRegs%vols(i)%c2P(tama3)%Or   !UPDATE diodos POR SI HAY ROTACION. ojo es un chapuz solo valido para diodos de 1p o 2p
         end if
      end do
      !surfaces
      tama = (this%DielRegs%nsurfs)
      do i = 1, tama            
         tama2 = (this%DielRegs%surfs(i)%n_c1P)
         do ii = 1, tama2
            call ROTATEMPI(mpidir,this%DielRegs%surfs(i)%C1P(ii)) 
         end do
         if (tama2 > 0) then      
            this%DielRegs%surfs(i)%DiodOrI = this%DielRegs%surfs(i)%c1P(tama2)%Or    !UPDATE diodos POR SI HAY ROTACION ojo es un chapuz solo valido para diodos de 1p o 2p
         end if
         tama3 = (this%DielRegs%surfs(i)%n_c2P)  
         do ii = 1, tama3
            call ROTATEMPI(mpidir,this%DielRegs%surfs(i)%C2P(ii)) 
         end do
         if (tama3 > 0) then    
            this%DielRegs%surfs(i)%DiodOrI = this%DielRegs%surfs(i)%c2P(tama3)%Or   !UPDATE diodos POR SI HAY ROTACION. ojo es un chapuz solo valido para diodos de 1p o 2p
         end if
      end do
      !lines
      tama = (this%DielRegs%nlins)
      do i = 1, tama            
         tama2 = (this%DielRegs%lins(i)%n_c1P)
         do ii = 1, tama2
            call ROTATEMPI(mpidir,this%DielRegs%lins(i)%C1P(ii)) 
         end do
         if (tama2 > 0) then       
            this%DielRegs%lins(i)%DiodOrI = this%DielRegs%lins(i)%c1P(tama2)%Or    !UPDATE diodos POR SI HAY ROTACION ojo es un chapuz solo valido para diodos de 1p o 2p
         end if

         tama3 = (this%DielRegs%lins(i)%n_c2P)  
         do ii = 1, tama3
            call ROTATEMPI(mpidir,this%DielRegs%lins(i)%C2P(ii)) 
         end do
         if (tama3 > 0) then        
            this%DielRegs%lins(i)%DiodOrI = this%DielRegs%lins(i)%c2P(tama3)%Or   !UPDATE diodos POR SI HAY ROTACION. ojo es un chapuz solo valido para diodos de 1p o 2p
         end if
      end do
      return
   end subroutine rotate_generateNONMetals
   
   subroutine rotate_generateANISOTROPICs (this,mpidir)         
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      if ((mpidir/=1).and.(this%ANIMATS%nvols+this%ANIMATS%nsurfs+this%ANIMATS%nlins/=0)) then
           print *,'Rotations in anisotropic unsupported'
           stop
           return
      end if
   !si algun dia lo hubiera es un cut and paste de rotate_generateNONMetals y a la que hay que aniadir la rotacion de la matriz de medios 
   return
   end subroutine rotate_generateANISOTROPICs
   
   subroutine rotate_generateThinWires (this,mpidir)     
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir     
      integer(kind=4) :: tama,tama2,i,ii
      integer(kind=4) :: oldx, oldy, oldz
      
      tama = this%twires%n_tw
      do i=1, tama       
         tama2 = this%twires%TW(i)%N_TWC
         do ii = 1, tama2
      !!!ROTATE THINWIRE
             if (MPIDIR==2 ) then
                   oldx = this%twires%tw(i)%tWc(ii)%i
                   oldy = this%twires%tw(i)%tWc(ii)%j
                   oldz = this%twires%tw(i)%tWc(ii)%K
                   this%twires%tw(i)%tWc(ii)%i = oldz
                   this%twires%tw(i)%tWc(ii)%j = oldx
                   this%twires%tw(i)%tWc(ii)%K = oldy

                   SELECT CASE (this%twires%tw(i)%tWc(ii)%d)
                    CASE (iEx)
                      this%twires%tw(i)%tWc(ii)%d = iEy
                    CASE (iEY)
                      this%twires%tw(i)%tWc(ii)%d = iEz
                    CASE (iEZ)
                      this%twires%tw(i)%tWc(ii)%d = iEx
                   end select
            ELSEIF (MPIDIR==1 ) then
                      oldx = this%twires%tw(i)%tWc(ii)%i
                      oldy = this%twires%tw(i)%tWc(ii)%j
                      oldz = this%twires%tw(i)%tWc(ii)%K
                      this%twires%tw(i)%tWc(ii)%i = oldy
                      this%twires%tw(i)%tWc(ii)%j = oldz
                      this%twires%tw(i)%tWc(ii)%K = oldx
      
                   SELECT CASE (this%twires%tw(i)%tWc(ii)%d)
                    CASE (iEx)
                      this%twires%tw(i)%tWc(ii)%d = iEz
                    CASE (iEY)
                      this%twires%tw(i)%tWc(ii)%d = iEx
                    CASE (iEZ)
                      this%twires%tw(i)%tWc(ii)%d = iEy
                   end select
             end if
      !!!FIN  
         end do
      end do
     
      return
   end subroutine rotate_generateThinWires
   
   subroutine rotate_generateSlantedWires (this,mpidir)      
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      type(SlantedWiresInfo_t), pointer :: old_swires    
      real(kind=8) :: oldx, oldy, oldz
      integer(kind=4) :: tama,tama2,i,ii
      
      tama = this%swires%n_sw
      allocate(old_swires,source=this%swires)
      do i=1, tama       
         tama2 = this%swires%sW(i)%N_SWC
         do ii=1,tama2
             !!!ROTATE THINWIRE
             if (MPIDIR==2 ) then
                      oldx = this%swires%sw(i)%swc(ii)%x
                      oldy = this%swires%sw(i)%swc(ii)%y
                      oldz = this%swires%sw(i)%swc(ii)%z

                      this%swires%sw(i)%swc(ii)%x = oldz
                      this%swires%sw(i)%swc(ii)%y = oldx
                      this%swires%sw(i)%swc(ii)%z = oldy
             ELSEIF (MPIDIR==1 ) then
                      oldx = this%swires%sw(i)%swc(ii)%x
                      oldy = this%swires%sw(i)%swc(ii)%y
                      oldz = this%swires%sw(i)%swc(ii)%z

                      this%swires%sw(i)%swc(ii)%x = oldy
                      this%swires%sw(i)%swc(ii)%y = oldz
                      this%swires%sw(i)%swc(ii)%z = oldx
             end if
         end do
         !!!FIN
      end do   
      allocate(old_swires)

      return
   end subroutine rotate_generateSlantedWires
   
   subroutine rotate_generateThinSlots (this,mpidir)     
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      type(ThinSlots_t), pointer :: old_tSlots    
      integer(kind=4) :: tama,tama2,i,ii
      integer(kind=4) :: oldx, oldy, oldz
      
      tama = this%tSlots%n_Tg
      allocate(old_tSlots,source=this%tSlots)
      do i=1, tama       
         tama2 = this%tSlots%Tg(i)%N_Tgc
         do ii = 1, tama2
              !!!ROTATE THIN SLOT
              if (MPIDIR==2 ) then
                       oldx = this%tSlots%Tg(i)%TgC(ii)%i
                       oldy = this%tSlots%Tg(i)%TgC(ii)%j
                       oldz = this%tSlots%Tg(i)%TgC(ii)%K

                       this%tSlots%Tg(i)%TgC(ii)%i = oldz
                       this%tSlots%Tg(i)%TgC(ii)%j = oldx
                       this%tSlots%Tg(i)%TgC(ii)%K = oldy
                       SELECT CASE (old_tsLots%Tg(i)%TgC(ii)%dir)
                        CASE (iEx)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEy
                        CASE (iEY)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEz
                        CASE (iEZ)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEx
                       end select
              ELSEIF (MPIDIR==1 ) then
                       oldx = this%tSlots%Tg(i)%TgC(ii)%i
                       oldy = this%tSlots%Tg(i)%TgC(ii)%j
                       oldz = this%tSlots%Tg(i)%TgC(ii)%K
                       
                       this%tSlots%Tg(i)%TgC(ii)%i = oldy
                       this%tSlots%Tg(i)%TgC(ii)%j = oldz
                       this%tSlots%Tg(i)%TgC(ii)%K = oldx

                       SELECT CASE (old_tsLots%Tg(i)%TgC(ii)%dir)
                        CASE (iEx)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEz
                        CASE (iEY)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEx
                        CASE (iEZ)
                          this%tSlots%Tg(i)%TgC(ii)%dir = iEy
                       end select
              end if
        end do
      end do      
      return
      
   end subroutine rotate_generateThinSlots
!!   
   subroutine rotate_generateLossyThinSurface (this,mpidir)    
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama2,tama,i,ii
                                         
      tama = this%LossyThinSurfs%length
      do i = 1, tama
         tama2 = this%LossyThinSurfs%cs(i)%nc
          do ii = 1, tama2
             call ROTATEMPI(mpidir,this%LossyThinSurfs%cs(i)%C(ii)) 
          end do   
      end do    

      return
   end subroutine rotate_generateLossyThinSurface

   subroutine rotate_generateFDMs (this,mpidir) 
      
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama,tama2,i,ii
      
      tama = (this%FRQDEPMATS%nvols)       
      do i = 1, tama   
         tama2 = this%FRQDEPMATS%vols(i)%n_C
         do ii = 1, tama2
            call ROTATEMPI(mpidir,this%FRQDEPMATS%Vols(i)%c(ii))  
         end do
         call rotate_freq_depend_material_properties(mpidir,this%FRQDEPMATS%Vols(i))
      end do
      
      tama = (this%FRQDEPMATS%nsurfs)    
      do i = 1, tama 
         tama2 = this%FRQDEPMATS%surfs(i)%n_C
         do ii = 1, tama2
            call ROTATEMPI(mpidir,this%FRQDEPMATS%Surfs(i)%c(ii))   
         end do
         call rotate_freq_depend_material_properties(mpidir,this%FRQDEPMATS%Surfs(i))
      end do
      
      tama = (this%FRQDEPMATS%nlins)
      do i = 1, tama 
         tama2 = this%FRQDEPMATS%Lins(i)%n_C
         do ii = 1, tama2
            call ROTATEMPI(mpidir,this%FRQDEPMATS%Lins(i)%c(ii))   
         end do
         call rotate_freq_depend_material_properties(mpidir,this%FRQDEPMATS%Lins(i))
      end do
      return
      
   end subroutine rotate_generateFDMs
   
   subroutine rotate_generateSONDAs (this,mpidir) 
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir    
      integer(kind=4) :: tama,tama2,tama3,i,ii,iii  
      type(FarField_Sonda_t), pointer :: old_FarField => NULL ()
      type(Electric_Sonda_t), pointer :: old_Electric => NULL ()
      type(Magnetic_Sonda_t), pointer :: old_Magnetic => NULL ()
      real(kind=RK) :: THETASTART,THETASTOP,PHISTART,PHISTOP
      integer :: iox, ioy, ioz
      
      tama = this%oldSONDA%n_probes        
      ! tres posibilidades FarField, Electric,Magnetic
      do i = 1, tama      
         tama2 = (this%oldSONDA%probes(i)%n_FarField)  
         do ii = 1, tama2                                 
            allocate (old_FarField,source=this%oldSONDA%probes(i)%FarField(ii))
            thetastart=old_FarField%probe%thetastart
            thetastop =old_FarField%probe%thetastop   
            phistart  =old_FarField%probe%phistart
            phistop   =old_FarField%probe%phistop
            !!!mpirotate angulos farfield .... las coordenadas se rotan luego
            if (MPIDIR==2 ) then
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastart = atan2(Sqrt(Cos(thetastart)**2.0_RKIND+ Cos(phistart)**2*Sin(thetastart)**2),Sin(phistart)*Sin(thetastart))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%phistart = atan2(Cos(phistart)*Sin(thetastart),Cos(thetastart))      
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastop = atan2(Sqrt(Cos(thetastop)**2.0_RKIND+ Cos(phistop)**2*Sin(thetastop)**2),Sin(phistop)*Sin(thetastop))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%phistop   = atan2(Cos(phistop)*Sin(thetastop),Cos(thetastop))
            ELSEIF (MPIDIR==1 ) then
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastart = atan2(Sqrt(Cos(thetastart)**2.0_RKIND+ Sin(phistart)**2*Sin(thetastart)**2),Cos(phistart)*Sin(thetastart))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%phistart   = atan2(Cos(thetastart),Sin(phistart)*Sin(thetastart))    
                   this%oldSONDA%probes(i)%FarField(ii)%probe%thetastop = atan2(Sqrt(Cos(thetastop)**2.0_RKIND+ Sin(phistop)**2*Sin(thetastop)**2),Cos(phistop)*Sin(thetastop))
                   this%oldSONDA%probes(i)%FarField(ii)%probe%phistop   = atan2(Cos(thetastop),Sin(phistop)*Sin(thetastop))
            end if        
            tama3 = (this%oldSONDA%probes(i)%FarField(ii)%probe%n_cord)    
            do iii = 1, tama3
              !!!ROTATE MPI
              if (MPIDIR==2 ) then
                  iox = this%oldSONDA%probes(i)%FarField(ii)%probe%i(iii)
                  ioy = this%oldSONDA%probes(i)%FarField(ii)%probe%j(iii)
                  ioz = this%oldSONDA%probes(i)%FarField(ii)%probe%K(iii)
                 
                 this%oldSONDA%probes(i)%FarField(ii)%probe%i(iii) = ioz
                 this%oldSONDA%probes(i)%FarField(ii)%probe%j(iii) = iox
                 this%oldSONDA%probes(i)%FarField(ii)%probe%K(iii) = ioy
              ELSEIF (MPIDIR==1 ) then                    
                 iox = this%oldSONDA%probes(i)%FarField(ii)%probe%i(iii)
                  ioy = this%oldSONDA%probes(i)%FarField(ii)%probe%j(iii)
                  ioz = this%oldSONDA%probes(i)%FarField(ii)%probe%K(iii)
                 
                 this%oldSONDA%probes(i)%FarField(ii)%probe%i(iii) = ioy
                 this%oldSONDA%probes(i)%FarField(ii)%probe%j(iii) = ioz
                 this%oldSONDA%probes(i)%FarField(ii)%probe%K(iii) = iox
              end if
            end do             
            deallocate(old_FarField)
         end do
      end do
      !     
      do i = 1, tama      
         tama2 = (this%oldSONDA%probes(i)%n_Electric)  
         do ii = 1, tama2                                 
            allocate (old_Electric,source=this%oldSONDA%probes(i)%Electric(ii))        
            tama3 = (this%oldSONDA%probes(i)%Electric(ii)%probe%n_cord)    
            do iii = 1, tama3
              !!!ROTATE MPI
              if (MPIDIR==2 ) then
                 iox = this%oldSONDA%probes(i)%Electric(ii)%probe%i(iii)
                  ioy = this%oldSONDA%probes(i)%Electric(ii)%probe%j(iii)
                  ioz = this%oldSONDA%probes(i)%Electric(ii)%probe%K(iii)
                 
                 this%oldSONDA%probes(i)%Electric(ii)%probe%i(iii) = ioz
                 this%oldSONDA%probes(i)%Electric(ii)%probe%j(iii) = iox
                 this%oldSONDA%probes(i)%Electric(ii)%probe%K(iii) = ioy
              ELSEIF (MPIDIR==1 ) then                    
                 iox = this%oldSONDA%probes(i)%Electric(ii)%probe%i(iii)
                  ioy = this%oldSONDA%probes(i)%Electric(ii)%probe%j(iii)
                  ioz = this%oldSONDA%probes(i)%Electric(ii)%probe%K(iii)
                 
                 this%oldSONDA%probes(i)%Electric(ii)%probe%i(iii) = ioy
                 this%oldSONDA%probes(i)%Electric(ii)%probe%j(iii) = ioz
                 this%oldSONDA%probes(i)%Electric(ii)%probe%K(iii) = iox
              end if
            end do             
            deallocate(old_Electric)
         end do
       end do
      !
      do i = 1, tama      
         tama2 = (this%oldSONDA%probes(i)%n_Magnetic)  
         do ii = 1, tama2                                 
            allocate (old_Magnetic,source=this%oldSONDA%probes(i)%Magnetic(ii))        
            tama3 = (this%oldSONDA%probes(i)%Magnetic(ii)%probe%n_cord)    
            do iii = 1, tama3
              !!!ROTATE MPI
              if (MPIDIR==2 ) then
                 iox = this%oldSONDA%probes(i)%Magnetic(ii)%probe%i(iii)
                  ioy = this%oldSONDA%probes(i)%Magnetic(ii)%probe%j(iii)
                  ioz = this%oldSONDA%probes(i)%Magnetic(ii)%probe%K(iii)
                 
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%i(iii) = ioz
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%j(iii) = iox
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%K(iii) = ioy
              ELSEIF (MPIDIR==1 ) then                    
                 iox = this%oldSONDA%probes(i)%Magnetic(ii)%probe%i(iii)
                  ioy = this%oldSONDA%probes(i)%Magnetic(ii)%probe%j(iii)
                  ioz = this%oldSONDA%probes(i)%Magnetic(ii)%probe%K(iii)
                 
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%i(iii) = ioy
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%j(iii) = ioz
                 this%oldSONDA%probes(i)%Magnetic(ii)%probe%K(iii) = iox
              end if
            end do             
            deallocate(old_Magnetic)
         end do
       end do
      !
       return
   end subroutine rotate_generateSONDAs
   
   subroutine rotate_generateMasSondas (this,mpidir)     
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir                          
      integer(kind=4) :: tama,tama2,i,ii  
      integer(kind=4) :: oxi,oyi,ozi,oxe,oye,oze,oor,TXI,TYI,TZI  
      type(coords_t), pointer :: old_MasSonda => NULL ()
      
      tama = this%Sonda%length       
      ! tres posibilidades FarField, Electric,Magnetic
      do i = 1, tama      
         tama2 = (this%Sonda%collection(i)%len_cor)    
         do ii = 1, tama2                                               
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
              if ((OOR/=NP_COR_EX).AND.(OOR/=NP_COR_EY).AND.(OOR/=NP_COR_EZ).AND. &
              (OOR/=NP_COR_HX).AND.(OOR/=NP_COR_HY).AND.(OOR/=NP_COR_HZ)) return
              !!LAS IW Y LAS VG NO SE ROTAN
              if (MPIDIR==2 ) then
                 this%Sonda%collection(i)%cordinates(ii)%XI=OZI   
                 this%Sonda%collection(i)%cordinates(ii)%XE=OZE
                 this%Sonda%collection(i)%cordinates(ii)%Xtrancos=TZI
                 
                 this%Sonda%collection(i)%cordinates(ii)%YI=OXI 
                 this%Sonda%collection(i)%cordinates(ii)%YE=OXE
                 this%Sonda%collection(i)%cordinates(ii)%Ytrancos=TXI
                 
                 this%Sonda%collection(i)%cordinates(ii)%ZI=OYI 
                 this%Sonda%collection(i)%cordinates(ii)%ZE=OYE
                 this%Sonda%collection(i)%cordinates(ii)%Ztrancos=TYI
                 
                 if (OOR== NP_COR_EX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EY
                 if (OOR== NP_COR_EY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EZ
                 if (OOR== NP_COR_EZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EX
                           
                 if (OOR== NP_COR_hX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_HY
                 if (OOR== NP_COR_hY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_HZ
                 if (OOR== NP_COR_hZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_HX
              ELSEIF (MPIDIR==1 ) then
                 this%Sonda%collection(i)%cordinates(ii)%XI=OYI
                 this%Sonda%collection(i)%cordinates(ii)%XE=OYE 
                 this%Sonda%collection(i)%cordinates(ii)%Xtrancos= TYI
                 
                 this%Sonda%collection(i)%cordinates(ii)%YI=OZI
                 this%Sonda%collection(i)%cordinates(ii)%YE=OZE               
                 this%Sonda%collection(i)%cordinates(ii)%Ytrancos=TZI
                 
                 this%Sonda%collection(i)%cordinates(ii)%ZI=OXI
                 this%Sonda%collection(i)%cordinates(ii)%ZE=OXE 
                 this%Sonda%collection(i)%cordinates(ii)%Ztrancos=TXI
                 
                 if (OOR== NP_COR_EX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EZ
                 if (OOR== NP_COR_EY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EX
                 if (OOR== NP_COR_EZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_EY
                           
                 if (OOR== NP_COR_HX) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_hZ
                 if (OOR== NP_COR_HY) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_hX
                 if (OOR== NP_COR_HZ) this%Sonda%collection(i)%cordinates(ii)%OR= NP_COR_hY
              end if                                              
              deallocate(old_MasSonda)
         end do
      end do
      return
   end subroutine rotate_generateMasSondas
   
   subroutine rotate_generateBloqueProbes (this,mpidir)    
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir                          
      integer(kind=4) :: tama,i  
      integer(kind=4) :: oxi,oyi,ozi,oxe,oye,oze  
      type(BloqueProbe_t), pointer :: old_BloqueProbe => NULL ()
      
      tama = this%BloquePRB%N_BP   
      do i = 1, tama      
          allocate(old_BloqueProbe,source=this%BloquePRB%BP(i))
          !MPI  ROTATE Bloque CURRENT
          if (MPIDIR==2 ) then
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
             end select
          ELSEIF (MPIDIR==1 ) then
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
             end select
          end if           
          deallocate(old_BloqueProbe)
      end do
      
          
      !!!!FIN ROTATE
      return
   end subroutine rotate_generateBloqueProbes
!!   
   subroutine rotate_generateVolumicProbes(this,mpidir)    
      type(Parseador_t), intent(inout) :: this
      integer(kind=4) :: mpidir                          
      integer(kind=4) :: tama,tama2,i,ii  
      integer(kind=4) :: oxi,oyi,ozi,oxe,oye,oze,oor,TXI,TYI,TZI  
      type(coords_t), pointer :: old_Coordinates => NULL ()
      
      tama = this%VolPrb%length
      ! tres posibilidades FarField, Electric,Magnetic
      do i = 1, tama      
         tama2 = (this%VolPrb%collection(i)%len_cor)    
         do ii = 1, tama2                                               
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
              if ((OOR/=iExC).AND.(OOR/=iEyC).AND.(OOR/=iEzC).AND. &
              (OOR/=iHxC).AND.(OOR/=iHyC).AND.(OOR/=iHzC).AND. &
              (OOR/=iCurX).AND.(OOR/=iCurY).AND.(OOR/=iCurZ).AND. &
               (OOR/=iMEC).AND.(OOR/=iMHC).AND.(OOR/=iCur)) return
              !!LAS IW Y LAS VG NO SE ROTAN.
        !!las imec, imhc e icur no le afecta el oor
              if (MPIDIR==2 ) then
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
                 if (OOR== iEXc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iEYc
                 if (OOR== iEYc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iEZc
                 if (OOR== iEZc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iEXc
                     
                 if (OOR== ihXc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iHYc
                 if (OOR== ihYc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iHZc
                 if (OOR== ihZc)   this%VolPrb%collection(i)%cordinates(ii)%OR= iHXc
                           
                 if (OOR== iCurX)  this%VolPrb%collection(i)%cordinates(ii)%OR= iCurY
                 if (OOR== iCurY)  this%VolPrb%collection(i)%cordinates(ii)%OR= iCurZ
                 if (OOR== iCurZ)  this%VolPrb%collection(i)%cordinates(ii)%OR= iCurX
              ELSEIF (MPIDIR==1 ) then
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
                 if (OOR== iEXc)  this%VolPrb%collection(i)%cordinates(ii)%OR= iEZc
                 if (OOR== iEYc)  this%VolPrb%collection(i)%cordinates(ii)%OR= iEXc
                 if (OOR== iEZc)  this%VolPrb%collection(i)%cordinates(ii)%OR= iEYc
                    
                 if (OOR== iHXc)  this%VolPrb%collection(i)%cordinates(ii)%OR= ihZc
                 if (OOR== iHYc)  this%VolPrb%collection(i)%cordinates(ii)%OR= ihXc
                 if (OOR== iHZc)  this%VolPrb%collection(i)%cordinates(ii)%OR= ihYc
                         
                 if (OOR== iCurX) this%VolPrb%collection(i)%cordinates(ii)%OR= iCurZ
                 if (OOR== iCurY) this%VolPrb%collection(i)%cordinates(ii)%OR= iCurX
                 if (OOR== iCurZ) this%VolPrb%collection(i)%cordinates(ii)%OR= iCurY
              end if                                                  
              deallocate(old_Coordinates)
         end do
      end do
      
      return
   end subroutine rotate_generateVolumicProbes

!!                                                    
!!!---------------------------------------------------->
!!!---------------------------------------------------->
!!!---------------------------------------------------->
!!
   subroutine ROTATEMPI(mpidir,COORDEN)         
      integer(kind=4) :: mpidir    
      type(coords_t), intent(inout) :: COORDEN
      integer(kind=4) :: OXI, OXE, OYI, OYE,OZI, OZE, OOR,TXI,TYI,TZI
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
      if (MPIDIR==2 ) then
         COORDEN%XI=OZI
         COORDEN%XE=OZE      
         COORDEN%Xtrancos =TZI    
         COORDEN%YI=OXI
         COORDEN%YE=OXE  
         COORDEN%Ytrancos =TXI 
         COORDEN%ZI=OYI
         COORDEN%ZE=OYE     
         COORDEN%Ztrancos =TYI 
         if (OOR== iEx) COORDEN%OR= iEy
         if (OOR==-iEx) COORDEN%OR=-iEy
         if (OOR== iEy) COORDEN%OR= iEz
         if (OOR==-iEy) COORDEN%OR=-iEz
         if (OOR== iEz) COORDEN%OR= iEx
         if (OOR==-iEz) COORDEN%OR=-iEx
      ELSEIF (MPIDIR==1 ) then
         COORDEN%XI=OYI
         COORDEN%XE=OYE         
         COORDEN%Xtrancos =TYI    
         COORDEN%YI=OZI
         COORDEN%YE=OZE   
         COORDEN%Ytrancos =TZI 
         COORDEN%ZI=OXI
         COORDEN%ZE=OXE 
         COORDEN%Ztrancos =TXI 
         if (OOR== iEx) COORDEN%OR= iEz
         if (OOR==-iEx) COORDEN%OR=-iEz
         if (OOR== iEy) COORDEN%OR= iEx
         if (OOR==-iEy) COORDEN%OR=-iEx
         if (OOR== iEz) COORDEN%OR= iEy
         if (OOR==-iEz) COORDEN%OR=-iEy
      end if
      return
   end subroutine ROTATEMPI

   subroutine ROTATEMPI_SCALED(mpidir,COORDEN)        
      integer(kind=4) :: mpidir    
      type(coords_scaled_t), intent(inout) :: COORDEN
      integer(kind=4) :: OXI, OXE, OYI, OYE,OZI, OZE,oor
      real(kind=RK) :: OXC,OYC,OZC
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
      if (MPIDIR==2 ) then
         COORDEN%XI=OZI
         COORDEN%XE=OZE
         COORDEN%YI=OXI
         COORDEN%YE=OXE
         COORDEN%ZI=OYI
         COORDEN%ZE=OYE
         COORDEN%XC=OzC
         COORDEN%YC=OxC
         COORDEN%ZC=OyC  
         if (OOR== iEx) COORDEN%OR= iEy
         if (OOR==-iEx) COORDEN%OR=-iEy
         if (OOR== iEy) COORDEN%OR= iEz
         if (OOR==-iEy) COORDEN%OR=-iEz
         if (OOR== iEz) COORDEN%OR= iEx
         if (OOR==-iEz) COORDEN%OR=-iEx

      ELSEIF (MPIDIR==1 ) then
         COORDEN%XI=OYI
         COORDEN%XE=OYE
         COORDEN%YI=OZI
         COORDEN%YE=OZE
         COORDEN%ZI=OXI
         COORDEN%ZE=OXE
         COORDEN%XC=OyC
         COORDEN%YC=OzC
         COORDEN%ZC=OxC         
         if (OOR== iEx) COORDEN%OR= iEz
         if (OOR==-iEx) COORDEN%OR=-iEz
         if (OOR== iEy) COORDEN%OR= iEx
         if (OOR==-iEy) COORDEN%OR=-iEx
         if (OOR== iEz) COORDEN%OR= iEy
         if (OOR==-iEz) COORDEN%OR=-iEy

      end if
      return
   end subroutine ROTATEMPI_SCALED
     
   subroutine rotate_freq_depend_material_properties(mpidir, freqDepMat)
      integer(kind=4), intent(in) :: mpidir    
      type(FreqDepenMaterial_t), intent(inout) :: freqDepMat
      complex, dimension(:), pointer :: po11, po12, po13, po22, po23, po33
      real(kind=RK) :: ro11, ro12, ro13, ro22, ro23, ro33
      integer(kind=4) :: io11, io12, io13, io22, io23, io33
!!    
      if (mpidir==2) then
         !Rotate a matrix
         po11 => freqDepMat%a11
         po12 => freqDepMat%a12
         po13 => freqDepMat%a13
         po22 => freqDepMat%a22
         po23 => freqDepMat%a23
         po33 => freqDepMat%a33

         freqDepMat%a11 => po33 
         freqDepMat%a12 => po23 
         freqDepMat%a13 => po12 
         freqDepMat%a22 => po11 
         freqDepMat%a23 => po13 
         freqDepMat%a33 => po22
         
         !Rotate am matrix
         po11 => freqDepMat%am11
         po12 => freqDepMat%am12
         po13 => freqDepMat%am13
         po22 => freqDepMat%am22
         po23 => freqDepMat%am23
         po33 => freqDepMat%am33

         freqDepMat%am11 => po33 
         freqDepMat%am12 => po23 
         freqDepMat%am13 => po12 
         freqDepMat%am22 => po11 
         freqDepMat%am23 => po13 
         freqDepMat%am33 => po22

         !Rotate b matrix
         po11 => freqDepMat%b11
         po12 => freqDepMat%b12
         po13 => freqDepMat%b13
         po22 => freqDepMat%b22
         po23 => freqDepMat%b23
         po33 => freqDepMat%b33

         freqDepMat%b11 => po33 
         freqDepMat%b12 => po23 
         freqDepMat%b13 => po12 
         freqDepMat%b22 => po11 
         freqDepMat%b23 => po13 
         freqDepMat%b33 => po22
         
         !Rotate bm matrix
         po11 => freqDepMat%bm11
         po12 => freqDepMat%bm12
         po13 => freqDepMat%bm13
         po22 => freqDepMat%bm22
         po23 => freqDepMat%bm23
         po33 => freqDepMat%bm33

         freqDepMat%bm11 => po33 
         freqDepMat%bm12 => po23 
         freqDepMat%bm13 => po12 
         freqDepMat%bm22 => po11 
         freqDepMat%bm23 => po13 
         freqDepMat%bm33 => po22

         !Rotate eps matrix
         ro11 = freqDepMat%eps11
         ro12 = freqDepMat%eps12
         ro13 = freqDepMat%eps13
         ro22 = freqDepMat%eps22
         ro23 = freqDepMat%eps23
         ro33 = freqDepMat%eps33

         freqDepMat%eps11 = ro33 
         freqDepMat%eps12 = ro23 
         freqDepMat%eps13 = ro12 
         freqDepMat%eps22 = ro11 
         freqDepMat%eps23 = ro13 
         freqDepMat%eps33 = ro22

         !Rotate mu matrix
         ro11 = freqDepMat%mu11
         ro12 = freqDepMat%mu12
         ro13 = freqDepMat%mu13
         ro22 = freqDepMat%mu22
         ro23 = freqDepMat%mu23
         ro33 = freqDepMat%mu33

         freqDepMat%mu11 = ro33 
         freqDepMat%mu12 = ro23 
         freqDepMat%mu13 = ro12 
         freqDepMat%mu22 = ro11 
         freqDepMat%mu23 = ro13 
         freqDepMat%mu33 = ro22

         !Rotate sigma matrix
         ro11 = freqDepMat%sigma11
         ro12 = freqDepMat%sigma12
         ro13 = freqDepMat%sigma13
         ro22 = freqDepMat%sigma22
         ro23 = freqDepMat%sigma23
         ro33 = freqDepMat%sigma33

         freqDepMat%sigma11 = ro33 
         freqDepMat%sigma12 = ro23 
         freqDepMat%sigma13 = ro12 
         freqDepMat%sigma22 = ro11 
         freqDepMat%sigma23 = ro13 
         freqDepMat%sigma33 = ro22

         !Rotate sigmam matrix
         ro11 = freqDepMat%sigmam11
         ro12 = freqDepMat%sigmam12
         ro13 = freqDepMat%sigmam13
         ro22 = freqDepMat%sigmam22
         ro23 = freqDepMat%sigmam23
         ro33 = freqDepMat%sigmam33

         freqDepMat%sigmam11 = ro33 
         freqDepMat%sigmam12 = ro23 
         freqDepMat%sigmam13 = ro12 
         freqDepMat%sigmam22 = ro11 
         freqDepMat%sigmam23 = ro13 
         freqDepMat%sigmam33 = ro22

         !Rotate K matrix
         io11 = freqDepMat%k11
         io12 = freqDepMat%k12
         io13 = freqDepMat%k13
         io22 = freqDepMat%k22
         io23 = freqDepMat%k23
         io33 = freqDepMat%k33

         freqDepMat%k11 = io33 
         freqDepMat%k12 = io23 
         freqDepMat%k13 = io12 
         freqDepMat%k22 = io11 
         freqDepMat%k23 = io13 
         freqDepMat%k33 = io22
         
         !Rotate Km matrix
         io11 = freqDepMat%km11
         io12 = freqDepMat%km12
         io13 = freqDepMat%km13
         io22 = freqDepMat%km22
         io23 = freqDepMat%km23
         io33 = freqDepMat%km33

         freqDepMat%km11 = io33 
         freqDepMat%km12 = io23 
         freqDepMat%km13 = io12 
         freqDepMat%km22 = io11 
         freqDepMat%km23 = io13 
         freqDepMat%km33 = io22
      end if 
      if (mpidir==1) then
         !Rotate a matrix
         po11 => freqDepMat%a11
         po12 => freqDepMat%a12
         po13 => freqDepMat%a13
         po22 => freqDepMat%a22
         po23 => freqDepMat%a23
         po33 => freqDepMat%a33

         freqDepMat%a11 => po22
         freqDepMat%a12 => po13 
         freqDepMat%a13 => po23 
         freqDepMat%a22 => po33 
         freqDepMat%a23 => po12
         freqDepMat%a33 => po11
         
         !Rotate am matrix
         po11 => freqDepMat%am11
         po12 => freqDepMat%am12
         po13 => freqDepMat%am13
         po22 => freqDepMat%am22
         po23 => freqDepMat%am23
         po33 => freqDepMat%am33

         freqDepMat%am11 => po22 
         freqDepMat%am12 => po13 
         freqDepMat%am13 => po23 
         freqDepMat%am22 => po33 
         freqDepMat%am23 => po12 
         freqDepMat%am33 => po11

         !Rotate b matrix
         po11 => freqDepMat%b11
         po12 => freqDepMat%b12
         po13 => freqDepMat%b13
         po22 => freqDepMat%b22
         po23 => freqDepMat%b23
         po33 => freqDepMat%b33

         freqDepMat%b11 => po22 
         freqDepMat%b12 => po13 
         freqDepMat%b13 => po23 
         freqDepMat%b22 => po33 
         freqDepMat%b23 => po12 
         freqDepMat%b33 => po11
         
         !Rotate bm matrix
         po11 => freqDepMat%bm11
         po12 => freqDepMat%bm12
         po13 => freqDepMat%bm13
         po22 => freqDepMat%bm22
         po23 => freqDepMat%bm23
         po33 => freqDepMat%bm33

         freqDepMat%bm11 => po22 
         freqDepMat%bm12 => po13 
         freqDepMat%bm13 => po23 
         freqDepMat%bm22 => po33 
         freqDepMat%bm23 => po12 
         freqDepMat%bm33 => po11

         !Rotate eps matrix
         ro11 = freqDepMat%eps11
         ro12 = freqDepMat%eps12
         ro13 = freqDepMat%eps13
         ro22 = freqDepMat%eps22
         ro23 = freqDepMat%eps23
         ro33 = freqDepMat%eps33

         freqDepMat%eps11 = ro22 
         freqDepMat%eps12 = ro13 
         freqDepMat%eps13 = ro23 
         freqDepMat%eps22 = ro33 
         freqDepMat%eps23 = ro12 
         freqDepMat%eps33 = ro11

         !Rotate mu matrix
         ro11 = freqDepMat%mu11
         ro12 = freqDepMat%mu12
         ro13 = freqDepMat%mu13
         ro22 = freqDepMat%mu22
         ro23 = freqDepMat%mu23
         ro33 = freqDepMat%mu33

         freqDepMat%mu11 = ro22 
         freqDepMat%mu12 = ro13 
         freqDepMat%mu13 = ro23 
         freqDepMat%mu22 = ro33 
         freqDepMat%mu23 = ro12 
         freqDepMat%mu33 = ro11

         !Rotate sigma matrix
         ro11 = freqDepMat%sigma11
         ro12 = freqDepMat%sigma12
         ro13 = freqDepMat%sigma13
         ro22 = freqDepMat%sigma22
         ro23 = freqDepMat%sigma23
         ro33 = freqDepMat%sigma33

         freqDepMat%sigma11 = ro22 
         freqDepMat%sigma12 = ro13 
         freqDepMat%sigma13 = ro23 
         freqDepMat%sigma22 = ro33 
         freqDepMat%sigma23 = ro12 
         freqDepMat%sigma33 = ro11

         !Rotate sigmam matrix
         ro11 = freqDepMat%sigmam11
         ro12 = freqDepMat%sigmam12
         ro13 = freqDepMat%sigmam13
         ro22 = freqDepMat%sigmam22
         ro23 = freqDepMat%sigmam23
         ro33 = freqDepMat%sigmam33

         freqDepMat%sigmam11 = ro22 
         freqDepMat%sigmam12 = ro13 
         freqDepMat%sigmam13 = ro23 
         freqDepMat%sigmam22 = ro33 
         freqDepMat%sigmam23 = ro12 
         freqDepMat%sigmam33 = ro11

         !Rotate K matrix
         io11 = freqDepMat%k11
         io12 = freqDepMat%k12
         io13 = freqDepMat%k13
         io22 = freqDepMat%k22
         io23 = freqDepMat%k23
         io33 = freqDepMat%k33

         freqDepMat%k11 = io22 
         freqDepMat%k12 = io13 
         freqDepMat%k13 = io23 
         freqDepMat%k22 = io33 
         freqDepMat%k23 = io12 
         freqDepMat%k33 = io11
         
         !Rotate Km matrix
         io11 = freqDepMat%km11
         io12 = freqDepMat%km12
         io13 = freqDepMat%km13
         io22 = freqDepMat%km22
         io23 = freqDepMat%km23
         io33 = freqDepMat%km33

         freqDepMat%km11 = io22 
         freqDepMat%km12 = io13 
         freqDepMat%km13 = io23 
         freqDepMat%km22 = io33 
         freqDepMat%km23 = io12 
         freqDepMat%km33 = io11
      end if 
   end subroutine


end module nfde_rotate_m