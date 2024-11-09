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
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module to handle the resuming of a problem
!  Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module resuming

   Use Report

   use fdetypes



   !Thin metals

#ifdef CompileWithSGBC
#ifdef CompileWithStochastic
   use SGBC_stoch
#else
   use SGBC_NOstoch
#endif  
#endif
   use PMLbodies

   use Lumped
#ifdef CompileWithNIBC
   use Multiports
#endif

   !EDispersives
#ifdef CompileWithEDispersives
   use EDispersives
   use MDispersives
#endif

   use farfield_m

   !Wires Thin Module
#ifdef CompileWithWires
   use HollandWires
#endif
#ifdef CompileWithBerengerWires
   use WiresBerenger
#ifdef CompileWithMPI
   use WiresBerenger_MPI
#endif
#endif   
#ifdef CompileWithSlantedWires
   use WiresSlanted
#endif

   !Plane Wave Module
   use Ilumina
   !Observation Module
   use Observa
   !PMC and PML Module
   use Borders_CPML
   use Borders_MUR



#ifdef CompileWithMPI
   use MPIComm
#endif
#ifdef CompileWithStochastic
   use MPI_stochastic
#endif


   implicit none
   private

   
!!!variables globales del modulo
   REAL (KIND=RKIND), save           ::  zvac,cluz
   REAL (KIND=RKIND), save           ::  eps0,mu0
!!!   
   integer (kind=4), parameter, private  ::  BLOCK_SIZE = 1024
   public ReadFields,flush_and_save_resume


contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Read the main stepping program fields from isk for resuming simulation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ReadFields(sggalloc,lastexecutedtimestep,lastexecutedtime,ultimodt,eps00,mu00,Ex,Ey,Ez,Hx,Hy,Hz)

      type (XYZlimit_t), dimension(1:6)  ::  sggalloc
      REAL (KIND=RKIND)   , intent(inout)      :: &
      Ex(sggalloc(iEx)%XI : sggalloc(iEx)%XE,sggalloc(iEx)%YI : sggalloc(iEx)%YE,sggalloc(iEx)%ZI : sggalloc(iEx)%ZE),&
      Ey(sggalloc(iEy)%XI : sggalloc(iEy)%XE,sggalloc(iEy)%YI : sggalloc(iEy)%YE,sggalloc(iEy)%ZI : sggalloc(iEy)%ZE),&
      Ez(sggalloc(iEz)%XI : sggalloc(iEz)%XE,sggalloc(iEz)%YI : sggalloc(iEz)%YE,sggalloc(iEz)%ZI : sggalloc(iEz)%ZE),&
      Hx(sggalloc(iHx)%XI : sggalloc(iHx)%XE,sggalloc(iHx)%YI : sggalloc(iHx)%YE,sggalloc(iHx)%ZI : sggalloc(iHx)%ZE),&
      Hy(sggalloc(iHy)%XI : sggalloc(iHy)%XE,sggalloc(iHy)%YI : sggalloc(iHy)%YE,sggalloc(iHy)%ZI : sggalloc(iHy)%ZE),&
      Hz(sggalloc(iHz)%XI : sggalloc(iHz)%XE,sggalloc(iHz)%YI : sggalloc(iHz)%YE,sggalloc(iHz)%ZI : sggalloc(iHz)%ZE)
      REAL (KIND=RKIND_tiempo) :: lastexecutedtime,ultimodt
      REAL (KIND=RKIND) :: eps00,mu00
      integer (kind=4)  ::  lastexecutedtimestep,i,j,k,i_block,n_block,ini,fin

      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      cluz=1.0_RKIND/sqrt(mu0*eps0)

      READ (14) lastexecutedtimestep,lastexecutedtime,ultimodt,eps0,mu0
      Do k=sggalloc(iEx)%ZI,sggalloc(iEx)%ZE
         Do j=sggalloc(iEx)%YI,sggalloc(iEx)%YE
            n_block = int(((sggalloc(iEx)%XE) - (sggalloc(iEx)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iEx)%XI
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               READ (14) ( Ex(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            End do
            READ (14) ( Ex(i,j,k), i = ini, sggalloc(iEx)%XE)
         End do
      End do
      Do k=sggalloc(iEy)%ZI,sggalloc(iEy)%ZE
         Do j=sggalloc(iEy)%YI,sggalloc(iEy)%YE
            n_block = int(((sggalloc(iEy)%XE) - (sggalloc(iEy)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iEy)%XI
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               READ (14) ( Ey(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            End do
            READ (14) ( Ey(i,j,k), i = ini, sggalloc(iEy)%XE)
         End do
      End do
      Do k=sggalloc(iEz)%ZI,sggalloc(iEz)%ZE
         Do j=sggalloc(iEz)%YI,sggalloc(iEz)%YE
            n_block = int(((sggalloc(iEz)%XE) - (sggalloc(iEz)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iEz)%XI
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               READ (14) ( Ez(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            End do
            READ (14) ( Ez(i,j,k), i = ini, sggalloc(iEz)%XE)
         End do
      End do
      Do k=sggalloc(iHx)%ZI,sggalloc(iHx)%ZE
         Do j=sggalloc(iHx)%YI,sggalloc(iHx)%YE
            n_block = int(((sggalloc(iHx)%XE) - (sggalloc(iHx)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iHx)%XI
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               READ (14) ( Hx(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            End do
            READ (14) ( Hx(i,j,k), i = ini, sggalloc(iHx)%XE)
         End do
      End do
      Do k=sggalloc(iHy)%ZI,sggalloc(iHy)%ZE
         Do j=sggalloc(iHy)%YI,sggalloc(iHy)%YE
            n_block = int(((sggalloc(iHy)%XE) - (sggalloc(iHy)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iHy)%XI
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               READ (14) ( Hy(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            End do
            READ (14) ( Hy(i,j,k), i = ini, sggalloc(iHy)%XE)
         End do
      End do
      Do k=sggalloc(iHz)%ZI,sggalloc(iHz)%ZE
         Do j=sggalloc(iHz)%YI,sggalloc(iHz)%YE
            n_block = int(((sggalloc(iHz)%XE) - (sggalloc(iHz)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iHz)%XI
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               READ (14) ( Hz(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            End do
            READ (14) ( Hz(i,j,k), i = ini, sggalloc(iHz)%XE)
         End do
      End do

      return
   end subroutine




   !---------------------------------------------------->
   !**************************************************************************************************
   subroutine flush_and_save_resume(sgg, b, layoutnumber, size, nentradaroot, nresumeable2, thereare, fin,eps00,mu00, everflushed,  &
   Ex, Ey, Ez, Hx, Hy, Hz,wiresflavor,simu_devia,stochastic)
      logical :: simu_devia,stochastic
      type (SGGFDTDINFO), intent(IN)    :: sgg
      !---------------------------> inputs <----------------------------------------------------------
      character(len=*), INTENT(in) :: wiresflavor
      integer (kind=4) :: ierr
      type( bounds_t), intent( IN)  ::  b
      integer( kind = 4), intent( IN)  ::  layoutnumber, size
      !--->
      character( LEN=*), intent( IN)  ::  nresumeable2, nEntradaRoot
      type( logic_control), intent( IN)  ::  thereare
      integer( kind=4), intent( IN)  ::  fin
      logical :: existe
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
      !--->
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hz
      !---------------------------> output <----------------------------------------------------------
      logical, intent( OUT)  ::  everflushed
      !---------------------------> variables locales <-----------------------------------------------
      character (LEN=BUFSIZE)  ::  whoami
      character (LEN=BUFSIZE)     ::  dubuf
      real (kind = RKIND) :: eps00,mu00
      !---------------------------> empieza flush_and_save_resume <-----------------------------------
      integer :: my_iostat
      
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      cluz=1.0_RKIND/sqrt(mu0*eps0)
      
      write( whoami, '(a,i5,a,i5,a)') '(', layoutnumber+1, '/', size,') '
      everflushed = .TRUE.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!  Flush observation data to disk
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !SYNC TO DISK ENERGY AND REPORTING FILES
      IF (layoutnumber == 0) THEN
         call flush(11)
         call flush(10)
      endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  Open unit 14 and store the fields of each module for resuming pruposesdata
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
#ifdef CompileWithOldSaving
      inquire (file = trim(adjustl( nresumeable2)),exist=existe)
      if (existe) then
         my_iostat=0
8766     if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8766,'.',layoutnumber,trim(adjustl( nresumeable2))//'.old'
         open ( 14, file = trim(adjustl( nresumeable2))//'.old', form = 'formatted',err=8766,iostat=my_iostat)
         write( 14, '(a)',err=634) '!END'
         close ( 14, status = 'delete',err=634)
         call rename(trim(adjustl( nresumeable2)),trim(adjustl( nresumeable2))//'.old')
      endif
#endif
      !
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      my_iostat=0
8776  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8776,'.',layoutnumber,trim(adjustl( nresumeable2))//'.old'
      open ( 14, file = trim(adjustl( nresumeable2)), form = 'formatted',err=8776,iostat=my_iostat)
      write( 14, '(a)',err=634) '!END'
      close ( 14, status = 'delete',err=634)
      !
      my_iostat=0
8777  if(my_iostat /= 0) write(*,fmt='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8777,'.',layoutnumber,trim(adjustl( nresumeable2))//'.old'
      open ( 14, file = trim(adjustl( nresumeable2)), form = 'unformatted',err=8777,iostat=my_iostat,status='new',action='write')
      !--->
      call StoreFields(sgg,fin,eps0,mu0, b, Ex, Ey, Ez, Hx, Hy, Hz)
      !this module data !warning the calling order must be the same that the calling to the init routines
      if( Thereare%PMLBorders)       call StoreFieldsCPMLBorders
      If (Thereare%PMLbodies)        call StorefieldsPMLbodies
      if( Thereare%MURBorders)       call StoreFieldsMURBorders
#ifdef CompileWithMPI
      !do an update of the currents to later read the currents OK
      if (size>1)  then
#ifdef CompileWithWires
         if ((trim(adjustl(wiresflavor))=='holland') .or. &
             (trim(adjustl(wiresflavor))=='transition')) then
             if ((size>1).and.(thereare%wires))   then
                call newFlushWiresMPI(layoutnumber,size)
             endif
#ifdef CompileWithStochastic
             if (stochastic)  then
                call syncstoch_mpi_wires(simu_devia,layoutnumber,size)
             endif
#endif             
             endif
#endif
#ifdef CompileWithBerengerWires
         if (trim(adjustl(wiresflavor))=='berenger') then
            call FlushWiresMPI_Berenger(layoutnumber,size)
         endif
#endif
      endif
      

#endif
      if( Thereare%Wires)       then
#ifdef CompileWithWires
         if ((trim(adjustl(wiresflavor))=='holland') .or. &
             (trim(adjustl(wiresflavor))=='transition')) then
            call StoreFieldsWires
         endif
#endif
#ifdef CompileWithBerengerWires
         if (trim(adjustl(wiresflavor))=='berenger') then
            call StoreFieldsWires_Berenger
         endif
#endif
#ifdef CompileWithSlantedWires
         if((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
            call StoreFieldsWires_Slanted
         endif
#endif
      endif

      
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
      if (stochastic)  then
         call syncstoch_mpi_lumped(simu_devia,layoutnumber,size)
      endif
#endif    
#endif    
      if (ThereAre%Lumpeds) call StoreFieldsLumpeds(stochastic)
      
#ifdef CompileWithSGBC
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
      if (stochastic)  then
         call syncstoch_mpi_SGBCs(simu_devia,layoutnumber,size)
      endif
#endif    
#endif    
      if( Thereare%SGBCs)       then
          call StoreFieldsSGBCs(stochastic)
      endif
      
#endif
#ifdef CompileWithNIBC
      if( Thereare%Multiports)       call StoreFieldsMultiports
#endif

#ifdef CompileWithEDispersives
      if( Thereare%EDispersives)     call StoreFieldsEDispersives
      if( Thereare%MDispersives)     call StoreFieldsMDispersives
#endif
      if( Thereare%PlaneWaveBoxes)     call StorePlaneWaves(sgg)
      if( Thereare%FarFields)       call StoreFarFields(b)  !called at initobservation
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      close (14,err=634)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'RESUMING FLUSHSAVEANDRESUME: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
      return
   endsubroutine flush_and_save_resume
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Flush the main stepping program fields to disk after simulation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine StoreFields( sgg,finaltimestep,eps0,mu0, b, Ex, Ey, Ez, Hx, Hy, Hz)
      !---------------------------> inputs <----------------------------------------------------------
      type (SGGFDTDINFO), intent(IN)   ::  sgg
      type( bounds_t), intent( IN)  ::  b
      integer( kind = 4), intent(IN)  ::  finaltimestep
      !--->
      real (kind = RKIND), dimension( 0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent( IN)  ::  Ex
      real (kind = RKIND), dimension( 0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent( IN)  ::  Ey
      real (kind = RKIND), dimension( 0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent( IN)  ::  Ez
      !--->
      real (kind = RKIND), dimension( 0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent( IN)  ::  Hx
      real (kind = RKIND), dimension( 0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent( IN)  ::  Hy
      real (kind = RKIND), dimension( 0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent( IN)  ::  Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer( kind = 4)  ::  i, j, k, i_block, n_block, ini, fin
      real (kind = RKIND) :: eps0,mu0,cluz,zvac
      !---------------------------> empieza StoreFields <---------------------------------------------
      write(14,err=634) finaltimestep,sgg%tiempo(finaltimestep),sgg%dt,eps0,mu0
      !--->
      do k = 0, b%Ex%NZ-1
         do j = 0, b%Ex%NY-1
            n_block = int( b%Ex%NX / BLOCK_SIZE)
            ini = 0
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) ( Ex( i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            enddo
            write(14,err=634) ( Ex( i, j, k), i = ini, b%Ex%NX-1)
         enddo
      enddo
      !--->
      do k = 0, b%Ey%NZ-1
         do j= 0, b%Ey%NY-1
            n_block = int( b%Ey%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) ( Ey( i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            enddo
            write(14,err=634) ( Ey(i,j,k), i = ini, b%Ey%NX-1)
         enddo
      enddo
      !--->
      do k = 0, b%Ez%NZ-1
         do j = 0, b%Ez%NY-1
            n_block = int( b%Ez%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) ( Ez( i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            enddo
            write(14,err=634) ( Ez( i, j, k), i = ini, b%Ez%NX-1)
         enddo
      enddo
      !--->
      do k = 0, b%Hx%NZ-1
         do j = 0, b%Hx%NY-1
            n_block = int( b%Hx%NX / BLOCK_SIZE)
            ini = 0
            Do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) ( Hx( i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            enddo
            write(14,err=634) ( Hx( i, j, k), i = ini, b%Hx%NX-1)
         enddo
      enddo
      !--->
      do k = 0, b%Hy%NZ-1
         do j = 0, b%Hy%NY-1
            n_block = int( b%Hy%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) ( Hy( i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            enddo
            write(14,err=634) ( Hy( i, j, k), i = ini, b%Hy%NX-1)
         enddo
      enddo
      !--->
      do k = 0, b%Hz%NZ-1
         do j = 0, b%Hz%NY-1
            n_block = int( b%Hz%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) ( Hz( i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            enddo
            write(14,err=634) ( Hz( i, j, k), i = ini, b%Hz%NX-1)
         enddo
      enddo

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'RESUMING STOREFIELDS: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine


end module

