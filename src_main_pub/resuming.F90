
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Module to handle the resuming of a problem
!  Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module resuming_m

   use Report_m

   use FDETYPES_m
#ifdef CompileWithStochastic
   use SGBC_stoch
#else
   use SGBC_nostoch_m
#endif  
   use PMLbodies_m
   use Lumped_m
#ifdef CompileWithNIBC
   use Multiports
#endif
   use EDispersives_m
   use Mdispersives_m
   use farfield_m
   use HollandWires_m
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
   use ilumina_m
   !PMC and PML Module
   use BORDERS_CPML_m
   use BORDERS_MUR_m



#ifdef CompileWithMPI
   use MPIcomm_m
#endif
#ifdef CompileWithStochastic
   use MPI_stochastic
#endif


   implicit none
   private

   
!!!variables globales del modulo
   real(kind=RKIND), save           :: zvac,cluz
   real(kind=RKIND), save           :: eps0,mu0
!!!   
   integer(kind=4), parameter, private  :: BLOCK_SIZE = 1024
   public ReadFields,flush_and_save_resume


contains


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Read the main stepping program fields from isk for resuming simulation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine ReadFields(sggalloc,lastexecutedtimestep,lastexecutedtime,ultimodt,eps00,mu00,Ex,Ey,Ez,Hx,Hy,Hz)

      type(XYZlimit_t), dimension(1:6) :: sggalloc
      real(kind=RKIND)   , intent(inout) :: &
      Ex(sggalloc(iEx)%XI : sggalloc(iEx)%XE,sggalloc(iEx)%YI : sggalloc(iEx)%YE,sggalloc(iEx)%ZI : sggalloc(iEx)%ZE),&
      Ey(sggalloc(iEy)%XI : sggalloc(iEy)%XE,sggalloc(iEy)%YI : sggalloc(iEy)%YE,sggalloc(iEy)%ZI : sggalloc(iEy)%ZE),&
      Ez(sggalloc(IEZ)%XI : sggalloc(IEZ)%XE,sggalloc(IEZ)%YI : sggalloc(IEZ)%YE,sggalloc(IEZ)%ZI : sggalloc(IEZ)%ZE),&
      Hx(sggalloc(IHX)%XI : sggalloc(IHX)%XE,sggalloc(IHX)%YI : sggalloc(IHX)%YE,sggalloc(IHX)%ZI : sggalloc(IHX)%ZE),&
      Hy(sggalloc(IHY)%XI : sggalloc(IHY)%XE,sggalloc(IHY)%YI : sggalloc(IHY)%YE,sggalloc(IHY)%ZI : sggalloc(IHY)%ZE),&
      Hz(sggalloc(IHZ)%XI : sggalloc(IHZ)%XE,sggalloc(IHZ)%YI : sggalloc(IHZ)%YE,sggalloc(IHZ)%ZI : sggalloc(IHZ)%ZE)
      real(kind=RKIND_TIEMPO) :: lastexecutedtime,ultimodt
      real(kind=RKIND) :: eps00,mu00
      integer(kind=4) :: lastexecutedtimestep,i,j,k,i_block,n_block,ini,fin

      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      cluz=1.0_RKIND/sqrt(mu0*eps0)

      read (14) lastexecutedtimestep,lastexecutedtime,ultimodt,eps0,mu0
      do k=sggalloc(iEx)%ZI,sggalloc(iEx)%ZE
         do j=sggalloc(iEx)%YI,sggalloc(iEx)%YE
            n_block = int(((sggalloc(iEx)%XE) - (sggalloc(iEx)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iEx)%XI
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               read (14) (Ex(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            read (14) (Ex(i,j,k), i = ini, sggalloc(iEx)%XE)
         end do
      end do
      do k=sggalloc(iEy)%ZI,sggalloc(iEy)%ZE
         do j=sggalloc(iEy)%YI,sggalloc(iEy)%YE
            n_block = int(((sggalloc(iEy)%XE) - (sggalloc(iEy)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(iEy)%XI
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               read (14) (Ey(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            read (14) (Ey(i,j,k), i = ini, sggalloc(iEy)%XE)
         end do
      end do
      do k=sggalloc(IEZ)%ZI,sggalloc(IEZ)%ZE
         do j=sggalloc(IEZ)%YI,sggalloc(IEZ)%YE
            n_block = int(((sggalloc(IEZ)%XE) - (sggalloc(IEZ)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(IEZ)%XI
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               read (14) (Ez(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            read (14) (Ez(i,j,k), i = ini, sggalloc(IEZ)%XE)
         end do
      end do
      do k=sggalloc(IHX)%ZI,sggalloc(IHX)%ZE
         do j=sggalloc(IHX)%YI,sggalloc(IHX)%YE
            n_block = int(((sggalloc(IHX)%XE) - (sggalloc(IHX)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(IHX)%XI
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               read (14) (Hx(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            read (14) (Hx(i,j,k), i = ini, sggalloc(IHX)%XE)
         end do
      end do
      do k=sggalloc(IHY)%ZI,sggalloc(IHY)%ZE
         do j=sggalloc(IHY)%YI,sggalloc(IHY)%YE
            n_block = int(((sggalloc(IHY)%XE) - (sggalloc(IHY)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(IHY)%XI
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               read (14) (Hy(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            read (14) (Hy(i,j,k), i = ini, sggalloc(IHY)%XE)
         end do
      end do
      do k=sggalloc(IHZ)%ZI,sggalloc(IHZ)%ZE
         do j=sggalloc(IHZ)%YI,sggalloc(IHZ)%YE
            n_block = int(((sggalloc(IHZ)%XE) - (sggalloc(IHZ)%XI) + 1) / BLOCK_SIZE)
            ini = sggalloc(IHZ)%XI
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               read (14) (Hz(i,j,k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            read (14) (Hz(i,j,k), i = ini, sggalloc(IHZ)%XE)
         end do
      end do

      return
   end subroutine




   !---------------------------------------------------->
   !**************************************************************************************************
   subroutine flush_and_save_resume(sgg, b, layoutnumber, num_procs, nentradaroot, nresumeable2, thereare, fin,eps00,mu00, everflushed,  &
   Ex, Ey, Ez, Hx, Hy, Hz,wiresflavor,simu_devia,stochastic)
      logical :: simu_devia,stochastic
      type(SGGFDTDINFO_t), intent(in) :: sgg
      !---------------------------> inputs <----------------------------------------------------------
      character(len=*), intent(in) :: wiresflavor
      integer(kind=4) :: ierr
      type(bounds_t), intent(in) :: b
      integer(kind = 4), intent(in) :: layoutnumber, num_procs
      !--->
      character(LEN=*), intent(in) :: nresumeable2, nEntradaRoot
      type(logic_control_t), intent(in) :: thereare
      integer(kind=4), intent(in) :: fin
      logical :: existe
      !--->
      real(kind = RKIND), dimension(0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent(in) :: Ex
      real(kind = RKIND), dimension(0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent(in) :: Ey
      real(kind = RKIND), dimension(0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent(in) :: Ez
      !--->
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(in) :: Hx
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(in) :: Hy
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(in) :: Hz
      !---------------------------> output <----------------------------------------------------------
      logical, intent(out) :: everflushed
      !---------------------------> variables locales <-----------------------------------------------
      character(len=BUFSIZE) :: whoami
      character(len=BUFSIZE) :: dubuf
      real(kind = RKIND) :: eps00,mu00
      !---------------------------> empieza flush_and_save_resume <-----------------------------------
      integer :: my_iostat
      
      eps0=eps00; mu0=mu00; !chapuz para convertir la variables de paso en globales
      zvac=sqrt(mu0/eps0)
      cluz=1.0_RKIND/sqrt(mu0*eps0)
      
      write(whoami, '(a,i5,a,i5,a)') '(', layoutnumber+1, '/', num_procs,') '
      everflushed = .TRUE.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!  Flush observation data to disk
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !SYNC TO DISK ENERGY AND REPORTING FILES
      if (layoutnumber == 0) then
         call flush(11)
         call flush(10)
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  Open unit 14 and store the fields of each module for resuming pruposesdata
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !
#ifdef CompileWithOldSaving
      inquire (file = trim(adjustl(nresumeable2)),exist=existe)
      if (existe) then
         my_iostat=0
8766     if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8766,'.',layoutnumber,trim(adjustl( nresumeable2))//'.old'
         open (14, file = trim(adjustl(nresumeable2))//'.old', form = 'formatted',err=8766,iostat=my_iostat)
         write(14, '(a)',err=634) '!END'
         close (14, status = 'delete',err=634)
         call rename(trim(adjustl(nresumeable2)),trim(adjustl(nresumeable2))//'.old')
      end if
#endif
      !
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      my_iostat=0
8776  if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8776,'.',layoutnumber,trim(adjustl( nresumeable2))//'.old'
      open (14, file = trim(adjustl(nresumeable2)), form = 'formatted',err=8776,iostat=my_iostat)
      write(14, '(a)',err=634) '!END'
      close (14, status = 'delete',err=634)
      !
      my_iostat=0
8777  if(my_iostat /= 0) write(*,FMT='(a)',advance='no'), '.' !!if(my_iostat /= 0) print '(i5,a1,i4,2x,a)',8777,'.',layoutnumber,trim(adjustl( nresumeable2))//'.old'
      open (14, file = trim(adjustl(nresumeable2)), form = 'unformatted',err=8777,iostat=my_iostat,status='new',action='write')
      !--->
      call StoreFields(sgg,fin,eps0,mu0, b, Ex, Ey, Ez, Hx, Hy, Hz)
      !this module data !warning the calling order must be the same that the calling to the init routines
      if(Thereare%PMLBorders)       call StoreFieldsCPMLBorders
      if (Thereare%PMLbodies)        call StorefieldsPMLbodies
      if(Thereare%MURBorders)       call StoreFieldsMURBorders
#ifdef CompileWithMPI
      !do an update of the currents to later read the currents OK
      if (num_procs>1)  then
         if ((trim(adjustl(wiresflavor))=='holland') .or. &
             (trim(adjustl(wiresflavor))=='transition')) then
             if ((num_procs>1).and.(thereare%wires))   then
                call newFlushWiresMPI(layoutnumber,num_procs)
             end if
#ifdef CompileWithStochastic
             if (stochastic)  then
                call syncstoch_mpi_wires(simu_devia,layoutnumber,num_procs)
             end if
#endif             
             end if
#ifdef CompileWithBerengerWires
         if (trim(adjustl(wiresflavor))=='berenger') then
            call FlushWiresMPI_Berenger(layoutnumber,num_procs)
         end if
#endif
      end if
      

#endif
      if(Thereare%Wires)       then
         if ((trim(adjustl(wiresflavor))=='holland') .or. &
             (trim(adjustl(wiresflavor))=='transition')) then
            call StoreFieldsWires
         end if
#ifdef CompileWithBerengerWires
         if (trim(adjustl(wiresflavor))=='berenger') then
            call StoreFieldsWires_Berenger
         end if
#endif
#ifdef CompileWithSlantedWires
         if((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
            call StoreFieldsWires_Slanted
         end if
#endif
      end if

      
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
      if (stochastic)  then
         call syncstoch_mpi_lumped(simu_devia,layoutnumber,num_procs)
      end if
#endif    
#endif    
      if (ThereAre%Lumpeds) call StoreFieldsLumpeds(stochastic)
      
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
      if (stochastic)  then
         call syncstoch_mpi_SGBCs(simu_devia,layoutnumber,num_procs)
      end if
#endif    
#endif    
      if(Thereare%SGBCs)       then
          call StoreFieldsSGBCs(stochastic)
      end if      
#ifdef CompileWithNIBC
      if(Thereare%Multiports)       call StoreFieldsMultiports
#endif
      if(Thereare%EDispersives)     call StoreFieldsEDispersives
      if(Thereare%MDispersives)     call StoreFieldsMDispersives
      if(Thereare%PlaneWaveBoxes)     call StorePlaneWaves(sgg)
      if(Thereare%FarFields)       call StoreFarFields(b)  !called at initobservation
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
   end subroutine flush_and_save_resume
   !**************************************************************************************************
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Flush the main stepping program fields to disk after simulation
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine StoreFields(sgg,finaltimestep,eps0,mu0, b, Ex, Ey, Ez, Hx, Hy, Hz)
      !---------------------------> inputs <----------------------------------------------------------
      type(SGGFDTDINFO_t), intent(in) :: sgg
      type(bounds_t), intent(in) :: b
      integer(kind = 4), intent(in) :: finaltimestep
      !--->
      real(kind = RKIND), dimension(0 :  b%Ex%NX-1, 0 :  b%Ex%NY-1, 0 :  b%Ex%NZ-1), intent(in) :: Ex
      real(kind = RKIND), dimension(0 :  b%Ey%NX-1, 0 :  b%Ey%NY-1, 0 :  b%Ey%NZ-1), intent(in) :: Ey
      real(kind = RKIND), dimension(0 :  b%Ez%NX-1, 0 :  b%Ez%NY-1, 0 :  b%Ez%NZ-1), intent(in) :: Ez
      !--->
      real(kind = RKIND), dimension(0 :  b%Hx%NX-1, 0 :  b%Hx%NY-1, 0 :  b%Hx%NZ-1), intent(in) :: Hx
      real(kind = RKIND), dimension(0 :  b%Hy%NX-1, 0 :  b%Hy%NY-1, 0 :  b%Hy%NZ-1), intent(in) :: Hy
      real(kind = RKIND), dimension(0 :  b%Hz%NX-1, 0 :  b%Hz%NY-1, 0 :  b%Hz%NZ-1), intent(in) :: Hz
      !---------------------------> variables locales <-----------------------------------------------
      integer(kind = 4) :: i, j, k, i_block, n_block, ini, fin
      real(kind = RKIND) :: eps0,mu0,cluz,zvac
      !---------------------------> empieza StoreFields <---------------------------------------------
      write(14,err=634) finaltimestep,sgg%tiempo(finaltimestep),sgg%dt,eps0,mu0
      !--->
      do k = 0, b%Ex%NZ-1
         do j = 0, b%Ex%NY-1
            n_block = int(b%Ex%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) (Ex(i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            write(14,err=634) (Ex(i, j, k), i = ini, b%Ex%NX-1)
         end do
      end do
      !--->
      do k = 0, b%Ey%NZ-1
         do j= 0, b%Ey%NY-1
            n_block = int(b%Ey%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) (Ey(i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            write(14,err=634) (Ey(i,j,k), i = ini, b%Ey%NX-1)
         end do
      end do
      !--->
      do k = 0, b%Ez%NZ-1
         do j = 0, b%Ez%NY-1
            n_block = int(b%Ez%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) (Ez(i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            write(14,err=634) (Ez(i, j, k), i = ini, b%Ez%NX-1)
         end do
      end do
      !--->
      do k = 0, b%Hx%NZ-1
         do j = 0, b%Hx%NY-1
            n_block = int(b%Hx%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) (Hx(i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            write(14,err=634) (Hx(i, j, k), i = ini, b%Hx%NX-1)
         end do
      end do
      !--->
      do k = 0, b%Hy%NZ-1
         do j = 0, b%Hy%NY-1
            n_block = int(b%Hy%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) (Hy(i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            write(14,err=634) (Hy(i, j, k), i = ini, b%Hy%NX-1)
         end do
      end do
      !--->
      do k = 0, b%Hz%NZ-1
         do j = 0, b%Hz%NY-1
            n_block = int(b%Hz%NX / BLOCK_SIZE)
            ini = 0
            do i_block = 1, n_block
               fin = ini-1 + BLOCK_SIZE
               write(14,err=634) (Hz(i, j, k), i = ini, fin)
               ini = ini + BLOCK_SIZE
            end do
            write(14,err=634) (Hz(i, j, k), i = ini, b%Hz%NX-1)
         end do
      end do

      goto 635
634   call print11(0,SEPARADOR//separador//separador)
      call print11(0,'RESUMING STOREFIELDS: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING')
      call print11(0,SEPARADOR//separador//separador)          
635   return
   end subroutine


end module
