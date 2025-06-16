module MPIcomm

#ifdef CompileWithMPI
   !
   Use Report
   use fdetypes

   use wiresHolland_constants
   use HollandWires

   implicit none
   private

   type(Thinwires_t), pointer  ::  HwiresMPI


   Type buffer_t
      REAL (KIND=RKIND_wires), pointer, dimension( : )  ::  SendUP,SendDown,RecUp,RecDown !for the wire exchange
      integer (kind=4)  ::  SendSizeUp,SendSizeDown,RecSizeUp ,RecSizeDown  !,uno
   end type

   Type ibuffer_t
      integer (KIND=4), pointer, dimension( : )  ::  SendUP,SendDown,RecUp,RecDown !for the wire exchange
      integer (kind=4)  ::  SendSizeUp,SendSizeDown,RecSizeUp,RecSizeDown !,uno
   end type

   !!!LOCAL VARIABLES
   type (buffer_t) , save ::  buffer
   type (ibuffer_t), save  ::  ibuffer

   Logical, SAVE :: FlushExtraInfoDown,FlushExtraInfoUp
   integer (kind=4), SAVE  ::  sizeHx,sizeHy,HxXI,HxXE,HyXI,HyXE,HxYI,HxYE,HyYI,HyYE,comZ,finZ
   integer (kind=4), SAVE  ::  sizeEx,sizeEy,ExXI,ExXE,EyXI,EyXE,ExYI,ExYE,EyYI,EyYE
   integer (kind=4), SAVE  ::  sizeEz,sizeHz,EzXI,EzXE,HzXI,HzXE,EzYI,EzYE,HzYI,HzYE
   integer (kind=4), SAVE,pointer, dimension(:) :: mpiZcom ,mpiZfin


   public FlushMPI_E,FlushMPI_H,InitMPI,MPIupdateMin, InitGeneralMPI,MPIdivide
   public MPIupdateBloques, MPIinitSubcomm
   !public InitWiresMPI
   public newInitWiresMPI,NewFlushWiresMPI
   public InitExtraFlushMPI


   !jag bug Antares mas de 65295 steps
   type, public :: t_databuf
      integer :: ip_target
      integer :: sizex, sizey, sizez
      logical :: FlushExtraInfo
      real (kind = RKIND), dimension( :, :), pointer :: buf_x_rx, buf_y_rx, buf_z_rx
      real (kind = RKIND), dimension( :, :), pointer :: buf_x_tx, buf_y_tx, buf_z_tx
   endtype t_databuf

   type :: t_databuf_Set
      logical :: syncUp, pbcUp
      type( t_databuf) :: databuf_Up
      !--->
      logical :: syncDown, pbcDown
      type( t_databuf) :: databuf_Down
   endtype t_databuf_Set
   !
   type (t_databuf_Set),save, target :: databuf_SetH
   type (t_databuf_Set),save, target :: databuf_SetE
   public FlushMPI_E_Cray,FlushMPI_H_Cray,InitMPI_Cray

   public InitExtraFlushMPI_Cray

contains


   subroutine InitGeneralMPI(layoutnumber,size)
      integer (kind=4) :: layoutnumber,size  !ojo he quitado el , intent(in) pq ambas son de salida 031218
      character (len=MPI_MAX_PROCESSOR_NAME)  ::  name
      integer (kind=4) namelen, ierr
      call MPI_INIT (ierr)
      call MPI_COMM_SIZE (MPI_COMM_WORLD, size, ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD, layoutnumber, ierr)
      call MPI_GET_PROCESSOR_NAME (name, namelen, ierr)
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      return
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! DIVIDES COMPUTATIONAL REGION INTO LAYOUTS AND READS GEOM DATA
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine MPIdivide(sgg,fullsize,SINPML_FULLSIZE,layoutnumber,size,forcing,forced,slicesoriginales,resume,fatalerror)

      type (SGGFDTDINFO), intent(INOUT)         ::  sgg
      type (limit_t), dimension(1:6)  ::  fullsize,SINPML_fullsize
      integer(kind=4) size, layoutnumber,ilay,padding,index(1:1),j
      integer(kind=4) forced
      logical :: forcing
      integer(kind=4), dimension( : ), pointer ::  trancos
      REAL (kind=RKIND), dimension( : ), pointer :: cZI,cZE
      REAL (kind=RKIND) :: carga,guess,ZE(1:3),cargaZE(1:3)
      REAL (kind=RKIND) :: deltatrancos
      character (len=*), intent(in)  ::  slicesoriginales
      character (len=BUFSIZE_LONG)  ::  slices=' '
      character (LEN=BUFSIZE)     ::  dubuf
      logical :: resume,fatalerror,originalPML_up_or_down
      character(len=BUFSIZE) :: buff 
      character (len=BUFSIZE)  ::  whoami
      integer(kind=4), dimension( 1:2 ) :: sggPMLNumLayers_original
      
      
      sggPMLNumLayers_original(:)=sgg%PML%NumLayers(3,:) !bug 310124 slices justo en PML
      originalPML_up_or_down=sgg%Border%IsUpPML.or.sgg%Border%IsDownPML !bug 310124 slices justo en PML
      
      WRITE (whoami, '(a,i5,a,i5,a)') '(', layoutnumber + 1, '/', size, ') '
      cZE => null(); cZI=> null(); trancos=> null(); mpizcom=> null(); mpizfin=> null();
      !
      !clip the simulation region
      !Take into account the PML overhead factor plusCPU_PML (2= double overhead, 1=no overhead)
      allocate(trancos(0 : size-1),cZI(0 : size),cZE(0 : size-1))
      carga= 1.0_RKIND *(fullsize(iHz)%ZE        - fullsize(iHz)%ZI)/(1.0_RKIND * size) + &
      (plusCPU_PML-1.0_RKIND)*((SINPML_fullsize(iHz)%ZI - fullsize(iHz)%ZI)   + &
      (fullsize(iHz)%ZE - SINPML_fullsize(iHz)%ZE))/(1.0_RKIND * size)
      cZI(0)=fullsize(iHz)%ZI
      do ilay=0,size-1
         guess=carga+cZI(ilay)+(plusCPU_PML-1.0_RKIND)*(min(cZI(ilay),1.0_RKIND * sinpml_fullsize(iHz)%ZI) + &
         max(cZI(ilay),1.0_RKIND * sinpml_fullsize(iHz)%ZE))
         ZE(1)=(guess-(plusCPU_PML-1.0_RKIND)*(sinpml_fullsize(iHz)%ZI))/(1.0_RKIND+(plusCPU_PML-1.0_RKIND))
         ZE(2)=(guess-(plusCPU_PML-1.0_RKIND)*(sinpml_fullsize(iHz)%ZE))/(1.0_RKIND+(plusCPU_PML-1.0_RKIND))
         ZE(3)=(guess-(plusCPU_PML-1.0_RKIND)*(sinpml_fullsize(iHz)%ZE)-(plusCPU_PML-1.0_RKIND)*(sinpml_fullsize(iHz)%ZI))
         do j=1,3
            cargaZE(j)=ABS( (ZE(j)-cZI(ilay)) + &
            (plusCPU_PML-1.0_RKIND)*(min(1.0_RKIND * SINPML_fullsize(iHz)%ZI,ZE(j)) - Min(1.0_RKIND * SINPML_fullsize(iHz)%ZI,cZI(ilay))  + &
            max(1.0_RKIND * SINPML_fullsize(iHz)%ZE,ZE(j)) - Max(1.0_RKIND * SINPML_fullsize(iHz)%ZE,cZI(ilay)))-carga)
         end do
         !select the closest to 0 one
         index=minloc(cargaZE)
         cZE(ilay)=ZE(index(1))
         cZI(ilay+1)=cZE(ilay)
      end do
 
     if (forcing) then
        if (size==2) then
            write(dubuf,*) 'Forcing MPI cut at ',forced
            call print11(layoutnumber,dubuf)
            cZI=-1;cZE=-1; !voided
            ilay=0
            !
            cZI(ilay)= fullsize(iHz)%ZI
            cZE(ilay)=forced
            ilay=ilay+1
            !repeat as many times as wanted for more cuts in a future
            cZI(ilay)=cZE(ilay-1)
            cZE(ilay)= fullsize(iHz)%ZE  !14
        else
            write(dubuf,*) 'Cannot force for more than 1 cut in a size=2 MPI'
            call print11(layoutnumber,dubuf,.true.)
        endif
     endif
      
      do ilay=0,size-1
         cZE(ilay)=nint(cZE(ilay))
         cZI(ilay+1)=cZE(ilay)
         trancos(ilay)=int(cZE(ilay)-cZI(0))
      enddo
      !end PML CPU overhead tunning

      allocate (mpiZcom(0:size-1),mpiZfin(0:size-1))
      mpiZcom(0)=fullsize(iHz)%ZI
      mpiZfin(0)=fullsize(iHz)%ZI+trancos(0)
      do ilay=1,size-2
         mpiZcom(ilay)=fullsize(iHz)%ZI+trancos(ilay-1)
         mpiZfin(ilay)=fullsize(iHz)%ZI+trancos(ilay)
      end do
      mpiZcom(size-1)=fullsize(iHz)%ZI+trancos(size-2)
      mpiZfin(size-1)=fullsize(iHz)%ZE

      !asign the  limits
      if ((LAYOUTNUMBER>0).and.(LAYOUTNUMBER<size-1)) then
         sgg%Sweep(1:6)%ZI=fullsize(1:6)%ZI+trancos(LAYOUTNUMBER-1)
         sgg%Sweep(1:6)%ZE=fullsize(1:6)%ZI+trancos(LAYOUTNUMBER  )
      elseif ((layoutnumber == 0).and.(LAYOUTNUMBER/=size-1))  then
         sgg%Sweep(1:6)%ZI=fullsize(1:6)%ZI
         sgg%Sweep(1:6)%ZE=fullsize(1:6)%ZI+trancos(LAYOUTNUMBER)
      elseif ((LAYOUTNUMBER/=0).and.(LAYOUTNUMBER==size-1)) then
         sgg%Sweep(1:6)%ZI=fullsize(1:6)%ZI+trancos(LAYOUTNUMBER-1)
         sgg%Sweep(1:6)%ZE=fullsize(1:6)%ZE
      endif
      !adjust THE ENDINGS OF THE INTERMDIATE computational limits
      if ((LAYOUTNUMBER>0).and.(LAYOUTNUMBER<size-1)) then
         !adjust computational limits
         sgg%Sweep(iEz)%ZE=sgg%Sweep(iEz)%ZE-1
         sgg%Sweep(iHx)%ZE=sgg%Sweep(iHx)%ZE-1
         sgg%Sweep(iHy)%ZE=sgg%Sweep(iHy)%ZE-1
      elseif ((layoutnumber == 0).and.(LAYOUTNUMBER/=size-1))  then
         !adjust computational limits
         sgg%Sweep(iEz)%ZE=sgg%Sweep(iEz)%ZE-1
         sgg%Sweep(iHx)%ZE=sgg%Sweep(iHx)%ZE-1
         sgg%Sweep(iHy)%ZE=sgg%Sweep(iHy)%ZE-1
      elseif ((LAYOUTNUMBER/=0).and.(LAYOUTNUMBER==size-1)) then
         continue
         !adjustment not necessary since fullsize%ZE is already adjusted
      endif


      !set also the dimensions of the computational layout which will be different from those of the mediamatrix
      !because of the extra padding

      padding=1
      if (padding >= minval(trancos)) then
         buff='Number of cells per processor less than 2. Decrease the number of MPI processors'
!intento recuperarme de este error
         call stoponerror(layoutnumber,size,buff,.true.); 
         if (associated(cZe)) deallocate(cZe,cZi)
         if (associated(mpizcom)) deallocate(mpizcom,mpizfin)
         if (associated(trancos)) deallocate(trancos)
         fatalerror=.true.
         return
      endif
      if ( minval(trancos) <= 2) then
         buff='Number of cells per processor less than 2. Decrease the number of MPI processors'
         call stoponerror(layoutnumber,size,buff,.true.)
         if (associated(cZe)) deallocate(cZe,cZi)
         if (associated(mpizcom)) deallocate(mpizcom,mpizfin)
         if (associated(trancos)) deallocate(trancos)
         fatalerror=.true.
         return
      endif
      !must be in agreement with the timestepping padding stuff
      if ((layoutnumber>0).and.(layoutnumber<size-1)) then
         sgg%alloc(1:6)%ZI=sgg%Sweep(1:6)%ZI - padding !I read one more MM for routines (wires, e.g.) requiring it
         sgg%alloc(1:6)%ZE=sgg%Sweep(1:6)%ZE + padding !I read one more MM for routines (wires, e.g.) requiring it
      elseif ((layoutnumber == 0).and.(layoutnumber/=size-1)) then
         sgg%alloc(1:6)%ZI=sgg%Sweep(1:6)%ZI-1  !I use this extra length in the global boundaries
         sgg%alloc(1:6)%ZE=sgg%Sweep(1:6)%ZE + padding !I read one more MM for routines (wires, e.g.) requiring it
      elseif ((layoutnumber/=0).and.(layoutnumber==size-1)) then
         sgg%alloc(1:6)%ZI=sgg%Sweep(1:6)%ZI - padding
         sgg%alloc(1:6)%ZE=sgg%Sweep(1:6)%ZE+1   !I use this extra length in the global boundaries
      endif

      !ARRANGE PML BORDERS
      if (layoutnumber == 0) then
         sgg%Border%IsUpPML=.false.
         sgg%Border%IsUpmur=.false.
         sgg%Border%IsUpPMC=.false.
         sgg%Border%IsUpPEC=.false.   !no PML layers UP
      elseif (layoutnumber==size-1) then
         sgg%Border%IsDownPML=.false.
         sgg%Border%IsDownmur=.false.
         sgg%Border%IsDownPMC=.false.
         sgg%Border%IsDownPEC=.false. !no PML layers DOWN
      else
         sgg%Border%IsUpPMC=.false.
         sgg%Border%IsUpPEC=.false.
         sgg%Border%IsDownPMC=.false.
         sgg%Border%IsDownPEC=.false.   
         !ojoo  en un futuro con este < a secas  cuando la PML esta justo en la division mpi 1310124
         if ((sgg%Sweep(iEx)%ZI<SINPML_fullsize(iEx)%ZI)) then
            sgg%Border%IsDownPML=.true.
         else
            sgg%Border%IsDownPML=.false. !no PML layers DOWN
         endif
         !ojoo  en un futuro con este > a secas  cuando la PML esta justo en la division mpi 1310124
         if ((sgg%Sweep(iEx)%ZE>SINPML_fullsize(iEx)%ZE))   then 
            sgg%Border%IsUpPML=.true.
         else
            sgg%Border%IsUpPML=.false.  !no PML layers UP
         endif
         sgg%Border%IsUpMur=.false.
         sgg%Border%IsDownMur=.false.
      endif

      !I readjust this since mpicomm touches this variable.
      if (.not.(sgg%Border%IsUpPML  )) sgg%PML%NumLayers(3,2)=0 !no PML layers UP
      if (.not.(sgg%Border%IsDownPML)) sgg%PML%NumLayers(3,1)=0 !no PML layers DOWN

      !writing
      if (layoutnumber==0) then
         slices='!SLICES'
         do ilay=0,size-1
            write(buff,'(i7)') mpiZfin(ilay)-mpiZcom(ilay)
            Slices=trim(adjustl(slices))//'_'//trim(adjustl(buff))
         END DO
         if (resume.and.(Slices /= SlicesOriginales)) then
            buff='Different resumed/original MPI slices: ' &
            //trim(adjustl(Slices))//' '//trim(adjustl(SlicesOriginales))
            call StopOnError(layoutnumber,size,buff)
         endif
         call print11(layoutnumber,trim(adjustl(slices)))
         do ilay=0,size-1
            write(buff,'(a,i5,a,i5,a,i7,a,i7,a,i7)')  &
            '(',ilay+1,'/',size,') Spanning from z=',mpiZcom(ilay),' to ',mpiZfin(ilay), &
            ' = ',mpiZfin(ilay)-mpiZcom(ilay)
            call print11(layoutnumber,buff)
         end do
      endif
      !end writing

      !bug 310124 cuando PML coincide con la primera celda de la ultima particion
      !no me complico y fuerzo menos MPI size
      if ((originalPML_up_or_down).and. &
             (mpiZfin(layoutnumber)-mpiZcom(layoutnumber)<=minval(sggPMLNumLayers_original))) then
           write(buff,'(a,i3,i3)') trim(adjustl(whoami))//' Minimum slice sizes along MPI should be larger that PML number of layers ', &
              mpiZfin(layoutnumber)-mpiZcom(layoutnumber), minval(sggPMLNumLayers_original)
            call StopOnError(layoutnumber,size,buff)
      endif
      deallocate(cZE,cZI,trancos,mpizcom,mpizfin)

      return
   end subroutine mpidivide
   !

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! INITIALIZES THE BOUNDS AND SIZES TO COMMUNICATE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMPI(sggsweep,sggalloc)
      type (XYZlimit_t), dimension(1:6)  ::  sggalloc,sggsweep
      !not necessary at this moment since nothing is read at this mmoment
      FlushExtraInfoDown=.false.
      FlushExtraInfoUp=.false.

      ExXI=sggalloc(iEx)%XI
      ExXE=sggalloc(iEx)%XE
      EyXI=sggalloc(iEy)%XI
      EyXE=sggalloc(iEy)%XE
      EzXI=sggalloc(iEz)%XI
      EzXE=sggalloc(iEz)%XE

      ExYI=sggalloc(iEx)%YI
      ExYE=sggalloc(iEx)%YE
      EyYI=sggalloc(iEy)%YI
      EyYE=sggalloc(iEy)%YE
      EzYI=sggalloc(iEz)%YI
      EzYE=sggalloc(iEz)%YE

      HxXI=sggalloc(iHx)%XI
      HxXE=sggalloc(iHx)%XE
      HyXI=sggalloc(iHy)%XI
      HyXE=sggalloc(iHy)%XE
      HzXI=sggalloc(iHz)%XI
      HzXE=sggalloc(iHz)%XE

      HxYI=sggalloc(iHx)%YI
      HxYE=sggalloc(iHx)%YE
      HyYI=sggalloc(iHy)%YI
      HyYE=sggalloc(iHy)%YE
      HzYI=sggalloc(iHz)%YI
      HzYE=sggalloc(iHz)%YE

      sizeEx=(ExXE-ExXI+1)*(ExYE-ExYI+1)
      sizeEy=(EyXE-EyXI+1)*(EyYE-EyYI+1)
      sizeEz=(EzXE-EzXI+1)*(EzYE-EzYI+1)

      sizeHx=(HxXE-HxXI+1)*(HxYE-HxYI+1)
      sizeHy=(HyXE-HyXI+1)*(HyYE-HyYI+1)
      sizeHz=(HzXE-HzXI+1)*(HzYE-HzYI+1)
      !
      ComZ=sggsweep(iHx)%ZI !both Hx and Hy coincide in this
      FinZ=sggsweep(iHx)%ZE

      return
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FORCES THE TIME STEP TO BE EQUAL BY EVERY LAYER
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine MPIupdateMin(dtlay,dt)
      real (kind=RKIND), intent(in)  ::  dtlay
      real (kind=RKIND), intent(out)  ::  dt
      integer (kind=4)  ::  ierr
      call MPI_AllReduce(dtlay, dt, 1_4, REALSIZE, MPI_MIN, SUBCOMM_MPI, ierr)
      return
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Sync the Bloque current data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine MPIupdateBloques(layoutnumber,valores,newvalores,SubComm)
      integer (kind=4)  ::  ierr,sizeofvalores,SubComm
      real (kind=RKIND), intent(in), dimension( 0:BuffObse )  ::  valores
      real (kind=RKIND), intent(out), dimension( 0:BuffObse )  ::  newvalores
      integer :: layoutnumber
      
      sizeofvalores=BuffObse+1
      call MPI_AllReduce(valores, newvalores, sizeofvalores, REALSIZE, MPI_SUM, SubComm, ierr)
   end subroutine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Init the sync the Bloque current data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine MPIinitSubcomm(layoutnumber,size,SubComm,Root,group1)
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  count,i
      integer (kind=4)  ::  ierr,wgroup,GROUP1,SubComm,Root,NewRoot
      logical, pointer, dimension( : )  ::  newallranks ,allranks
      integer (kind=4), pointer, dimension( : )   ::  NGroup
      allocate (allranks(0 : size-1),newallranks(0 : size-1))
      !         
      allranks=.false.
      !print *,'---layoutnumber, subcomm',layoutnumber,subcomm
      if (Subcomm == 1) allranks(layoutnumber)=.true.
      call MPI_AllReduce(allranks, newallranks, size, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
      !choose the maximum layer to be root !ojo no poner el minimo pq es -1 (voided de entrada)
      call MPI_AllReduce(Root, newRoot, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)
      Root=newRoot
      count=-1
      do i=0,size-1
         if (newallranks(i)) count=count+1
      end do
      allocate (NGroup(0 : count))
      count=-1
      do i=0,size-1
         if (newallranks(i)) then
            count=count+1
            NGroup(count)=i
         endif
      end do
      !all must create the same subcomm group though only some will synchronize
      if (count >= 0) then
         CALL MPI_COMM_GROUP(SUBCOMM_MPI, WGROUP, IERR) ! get default group
         CALL MPI_GROUP_INCL(WGROUP, COUNT+1_4, NGroup, GROUP1, IERR) !create a group called group1
         CALL MPI_COMM_CREATE(SUBCOMM_MPI, GROUP1, SUBCOMM, IERR)!create a communicator for group1
      else
         SUBCOMM=-1
         group1=-1
         WGROUP=-1
      endif
      if (.not.newallranks(layoutnumber)) subcomm=-1 !override it since it will never synchronize
      !print *,'Init lay,count, subcomm',layoutnumber+1,count,subcomm
      deallocate (allranks,newallranks)
      deallocate (NGroup)
      return
   end subroutine



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! TRANSFERS THE FIELDS AMONG LAYERS
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine FlushMPI_H(sggalloc,layoutnumber,size, Hx,Hy,Hz)
      type (XYZlimit_t), dimension (1:6), intent(in)                      ::  sggAlloc
!!!! reutilizo esta subroutina para sincronizar matrices de medios. Las de campo las hacen los _Cray 210815
!!!! ojo hay que cambiar tambien si algun dia se vuelve a esta rutina para comunicar reales INTEGERSIZE por REALSIZE
!!!!      REAL (KIND=RKIND)   , intent(inout)      :: &
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) , intent(inout)  ::  &
      Hx(sggalloc(iHx)%XI : sggalloc(iHx)%XE,sggalloc(iHx)%YI : sggalloc(iHx)%YE,sggalloc(iHx)%ZI : sggalloc(iHx)%ZE), &
      Hy(sggalloc(iHy)%XI : sggalloc(iHy)%XE,sggalloc(iHy)%YI : sggalloc(iHy)%YE,sggalloc(iHy)%ZI : sggalloc(iHy)%ZE), &
      Hz(sggalloc(iHz)%XI : sggalloc(iHz)%XE,sggalloc(iHz)%YI : sggalloc(iHz)%YE,sggalloc(iHz)%ZI : sggalloc(iHz)%ZE)




      integer (kind=4)  ::  ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7,ierr8,ierr9,ierr10,ierr11,ierr12,ierr100,ierr100b
      integer (kind=4)  ::  jerr1,jerr2,jerr3,jerr4,jerr5,jerr6,jerr7,jerr8,jerr9,jerr10,jerr11,jerr12,jerr100,jerr100b
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  req1 (1 : 4),req2 (1 : 4),status1 (MPI_STATUS_SIZE,1 : 4),status2 (MPI_STATUS_SIZE,1 : 4)
      integer (kind=4)  ::  req1b(1 : 2),req2b(1 : 2),status1b(MPI_STATUS_SIZE,1 : 2),status2b(MPI_STATUS_SIZE,1 : 2)



      ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;ierr100=0;ierr100b=0;
      jerr1=0;jerr2=0;jerr3=0;jerr4=0;jerr5=0;jerr6=0;jerr7=0;jerr8=0;jerr9=0;jerr10=0;jerr11=0;jerr12=0;jerr100=0;jerr100b=0;
      if (layoutnumber/=size-1) then !syncUp
         !print *,'---fluHup>',layoutnumber
         call MPI_IRECV     (Hx(HxXI,HxYI,finZ+1), sizeHx, INTEGERSIZE, layoutnumber+1_4,  1_4, SUBCOMM_MPI, req1(1 ), ierr1 )
         call MPI_ISEND     (Hx(HxXI,HxYI,finZ  ), sizeHx, INTEGERSIZE, layoutnumber+1_4,  2_4, SUBCOMM_MPI, req1(2 ), ierr2 )
         call MPI_IRECV     (Hy(HyXI,HyYI,finZ+1), sizeHy, INTEGERSIZE, layoutnumber+1_4,  3_4, SUBCOMM_MPI, req1(3 ), ierr3 )
         call MPI_ISEND     (Hy(HyXI,HyYI,finZ  ), sizeHy, INTEGERSIZE, layoutnumber+1_4,  4_4, SUBCOMM_MPI, req1(4 ), ierr4 )
         if (FlushExtraInfoUp) then
            !print *,'---fluHextraup>',layoutnumber
            call MPI_IRECV (Hz(HzXI,HzYI,finZ+2), sizeHz, INTEGERSIZE, layoutnumber+1_4, 5_4, SUBCOMM_MPI, req1b(1), ierr11)
            call MPI_ISEND (Hz(HzXI,HzYI,finZ  ), sizeHz, INTEGERSIZE, layoutnumber+1_4, 6_4, SUBCOMM_MPI, req1b(2), ierr12)
         endif
      ELSE !only NEEDED BY THE PERIODIC BOUNDARY CONDITIONS !no real MPI burden since all layers communicate two sets
         !print *,'---fluHup>',layoutnumber
         call MPI_IRECV     (Hx(HxXI,HxYI,finZ+1), sizeHx, INTEGERSIZE, 0_4,  1_4, SUBCOMM_MPI, req1(1 ), ierr1 )
         call MPI_ISEND     (Hx(HxXI,HxYI,finZ  ), sizeHx, INTEGERSIZE, 0_4,  2_4, SUBCOMM_MPI, req1(2 ), ierr2 )
         call MPI_IRECV     (Hy(HyXI,HyYI,finZ+1), sizeHy, INTEGERSIZE, 0_4,  3_4, SUBCOMM_MPI, req1(3 ), ierr3 )
         call MPI_ISEND     (Hy(HyXI,HyYI,finZ  ), sizeHy, INTEGERSIZE, 0_4,  4_4, SUBCOMM_MPI, req1(4 ), ierr4 )
      endif
      if (layoutnumber/=0    ) then !syncDown
         !print *,'---fluHdown>',layoutnumber
         call MPI_ISEND     (Hx(HxXI,HxYI,comZ  ), sizeHx, INTEGERSIZE, layoutnumber-1_4,  1_4, SUBCOMM_MPI, req2(1 ), jerr1 )
         call MPI_IRECV     (Hx(HxXI,HxYI,comZ-1), sizeHx, INTEGERSIZE, layoutnumber-1_4,  2_4, SUBCOMM_MPI, req2(2 ), jerr2 )
         call MPI_ISEND     (Hy(HyXI,HyYI,comZ  ), sizeHy, INTEGERSIZE, layoutnumber-1_4,  3_4, SUBCOMM_MPI, req2(3 ), jerr3 )
         call MPI_IRECV     (Hy(HyXI,HyYI,comZ-1), sizeHy, INTEGERSIZE, layoutnumber-1_4,  4_4, SUBCOMM_MPI, req2(4 ), jerr4 )
         if (FlushExtraInfoDown) then
            !print *,'---fluHextradown>',layoutnumber
            call MPI_ISEND (Hz(HzXI,HzYI,comZ+1), sizeHz, INTEGERSIZE, layoutnumber-1_4, 5_4, SUBCOMM_MPI, req2b(1), jerr11)
            call MPI_IRECV (Hz(HzXI,HzYI,comZ-1), sizeHz, INTEGERSIZE, layoutnumber-1_4, 6_4, SUBCOMM_MPI, req2b(2), jerr12)
         endif
      ELSE !only NEEDED BY THE PERIODIC BOUNDARY CONDITIONS
         call MPI_ISEND     (Hx(HxXI,HxYI,comZ  ), sizeHx, INTEGERSIZE, size-1_4,  1_4, SUBCOMM_MPI, req2(1 ), jerr1 )
         call MPI_IRECV     (Hx(HxXI,HxYI,comZ-1), sizeHx, INTEGERSIZE, size-1_4,  2_4, SUBCOMM_MPI, req2(2 ), jerr2 )
         call MPI_ISEND     (Hy(HyXI,HyYI,comZ  ), sizeHy, INTEGERSIZE, size-1_4,  3_4, SUBCOMM_MPI, req2(3 ), jerr3 )
         call MPI_IRECV     (Hy(HyXI,HyYI,comZ-1), sizeHy, INTEGERSIZE, size-1_4,  4_4, SUBCOMM_MPI, req2(4 ), jerr4 )
      endif
      !
      if (layoutnumber/=0    )  then
         call MPI_WAITALL(4_4,req2,status2,ierr100)
         if (FlushExtraInfoDown) then
            call MPI_WAITALL(2_4,req2b,status2b,ierr100b)
         endif
      endif
      if (layoutnumber/=size-1) then
         call MPI_WAITALL(4_4,req1,status1,jerr100)
         if (FlushExtraInfoUp) then
            call MPI_WAITALL(2_4,req1b,status1b,jerr100b)
         endif
      endif

      !
      !call MPI_Barrier(SUBCOMM_MPI,ierr12)
      !

      if (ierr1+ierr2+ierr3+ierr4+ierr5+ierr6+ierr7+ierr8+ierr9+ierr10+ierr11+ierr12+ierr100+ierr100b+ &
      jerr1+jerr2+jerr3+jerr4+jerr5+jerr6+jerr7+jerr8+jerr9+jerr10+jerr11+jerr12+jerr100+jerr100b /= 0) &
      call StopOnError(layoutnumber,size,'FLUSHMPI')

      return
   end subroutine
   !
   !

   subroutine FlushMPI_E(sggalloc,layoutnumber,size, Ex,Ey,Ez)
      type (XYZlimit_t), dimension (1:6), intent(in)                      ::  sggAlloc
!!!! reutilizo esta subroutina para sincronizar matrices de medios. Las de campo las hacen los _Cray 210815
!!!! ojo hay que cambiar tambien si algun dia se vuelve a esta rutina para comunicar reales INTEGERSIZE por REALSIZE
!!!!      REAL (KIND=RKIND)   , intent(inout)      :: &
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES) , intent(inout)  ::  &
      Ex(sggalloc(iEx)%XI : sggalloc(iEx)%XE,sggalloc(iEx)%YI : sggalloc(iEx)%YE,sggalloc(iEx)%ZI : sggalloc(iEx)%ZE),&
      Ey(sggalloc(iEy)%XI : sggalloc(iEy)%XE,sggalloc(iEy)%YI : sggalloc(iEy)%YE,sggalloc(iEy)%ZI : sggalloc(iEy)%ZE),&
      Ez(sggalloc(iEz)%XI : sggalloc(iEz)%XE,sggalloc(iEz)%YI : sggalloc(iEz)%YE,sggalloc(iEz)%ZI : sggalloc(iEz)%ZE)




      integer (kind=4)  ::  ierr1,ierr2,ierr3,ierr4,ierr5,ierr6,ierr7,ierr8,ierr9,ierr10,ierr11,ierr12,ierr100,ierr100b
      integer (kind=4)  ::  jerr1,jerr2,jerr3,jerr4,jerr5,jerr6,jerr7,jerr8,jerr9,jerr10,jerr11,jerr12,jerr100,jerr100b
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  req1 (1 : 2),req2 (1 : 2),status1 (MPI_STATUS_SIZE,1 : 2),status2 (MPI_STATUS_SIZE,1 : 2)
      integer (kind=4)  ::  req1b(1 : 4),req2b(1 : 4),status1b(MPI_STATUS_SIZE,1 : 4),status2b(MPI_STATUS_SIZE,1 : 4)

      ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;ierr100=0;ierr100b=0;
      jerr1=0;jerr2=0;jerr3=0;jerr4=0;jerr5=0;jerr6=0;jerr7=0;jerr8=0;jerr9=0;jerr10=0;jerr11=0;jerr12=0;jerr100=0;jerr100b=0;
      if (layoutnumber/=size-1) then !syncUp
         if (FlushExtraInfoUp) then
            !print *,'---fluEextraup>',layoutnumber
            call MPI_IRECV (Ez(EzXI,EzYI,finZ+1), sizeEz, INTEGERSIZE, layoutnumber+1_4,  1_4, SUBCOMM_MPI, req1(1 ), ierr5 )
            call MPI_ISEND (Ez(EzXI,EzYI,finZ  ), sizeEz, INTEGERSIZE, layoutnumber+1_4,  2_4, SUBCOMM_MPI, req1(2 ), ierr6 )
            !
            call MPI_IRECV (Ex(ExXI,ExYI,finZ+2), sizeEx, INTEGERSIZE, layoutnumber+1_4,  3_4, SUBCOMM_MPI, req1b(1), ierr7 )
            call MPI_ISEND (Ex(ExXI,ExYI,finZ  ), sizeEx, INTEGERSIZE, layoutnumber+1_4,  4_4, SUBCOMM_MPI, req1b(2), ierr8 )
            call MPI_IRECV (Ey(EyXI,EyYI,finZ+2), sizeEy, INTEGERSIZE, layoutnumber+1_4,  5_4, SUBCOMM_MPI, req1b(3), ierr9 )
            call MPI_ISEND (Ey(EyXI,EyYI,finZ  ), sizeEy, INTEGERSIZE, layoutnumber+1_4,  6_4, SUBCOMM_MPI, req1b(4), ierr10)
         endif
      endif
      if (layoutnumber/=0    ) then !syncDown
         if (FlushExtraInfoDown) then
            !print *,'---fluEextradown>',layoutnumber
            call MPI_ISEND (Ez(EzXI,EzYI,comZ  ), sizeEz, INTEGERSIZE, layoutnumber-1_4,  1_4, SUBCOMM_MPI, req2(1 ), jerr5 )
            call MPI_IRECV (Ez(EzXI,EzYI,comZ-1), sizeEz, INTEGERSIZE, layoutnumber-1_4,  2_4, SUBCOMM_MPI, req2(2 ), jerr6 )
            !
            call MPI_ISEND (Ex(ExXI,ExYI,comZ+1), sizeEx, INTEGERSIZE, layoutnumber-1_4,  3_4, SUBCOMM_MPI, req2b(1), jerr7 )
            call MPI_IRECV (Ex(ExXI,ExYI,comZ-1), sizeEx, INTEGERSIZE, layoutnumber-1_4,  4_4, SUBCOMM_MPI, req2b(2), jerr8 )
            call MPI_ISEND (Ey(EyXI,EyYI,comZ+1), sizeEy, INTEGERSIZE, layoutnumber-1_4,  5_4, SUBCOMM_MPI, req2b(3), jerr9 )
            call MPI_IRECV (Ey(EyXI,EyYI,comZ-1), sizeEy, INTEGERSIZE, layoutnumber-1_4,  6_4, SUBCOMM_MPI, req2b(4), jerr10)
         endif
      endif
      !
      if (layoutnumber/=0    )  then
         if (FlushExtraInfoDown) then
            call MPI_WAITALL(2_4,req2,status2,ierr100)
            call MPI_WAITALL(4_4,req2b,status2b,ierr100b)
         endif
      endif
      if (layoutnumber/=size-1) then
         if (FlushExtraInfoUp) then
            call MPI_WAITALL(2_4,req1,status1,jerr100)
            call MPI_WAITALL(4_4,req1b,status1b,jerr100b)
         endif
      endif

      !
      !call MPI_Barrier(SUBCOMM_MPI,ierr12)
      !

      if (ierr1+ierr2+ierr3+ierr4+ierr5+ierr6+ierr7+ierr8+ierr9+ierr10+ierr11+ierr12+ierr100+ierr100b+ &
      jerr1+jerr2+jerr3+jerr4+jerr5+jerr6+jerr7+jerr8+jerr9+jerr10+jerr11+jerr12+jerr100+jerr100b /= 0) &
      call StopOnError(layoutnumber,size,'FLUSHMPI')

      return
   end subroutine
   !
   !



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! new routine: works without the MediaMatrix Info
   ! supports multiwires
   !ADJUSTS THE WIRE DATA NEEDED FOR TRANSVER
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   subroutine newInitWiresMPI(layoutnumber,therearewires,size,resume,c)
      type (XYZlimit_t), dimension(1:6), intent(in)  ::  c
      !
      logical, intent(in)  ::  resume,therearewires  !controls wether to read something for resuming
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  i1,i,j
      integer (kind=4)  ::  SharescontaMPIdown,SharescontaMPIup,NeedscontaMPIdown,NeedscontaMPIup


      integer (kind=4)  ::  ni,nj,nk,norigindex,idum
      type (CurrentSegments), pointer  :: segmento
      character (LEN=BUFSIZE)  ::  whoami
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '



      !Get info from wires
      IF (therearewires) then
         HwiresMPI => GetHwires()
      else
         allocate (HwiresMPI)
         HwiresMPI%NumChargeNodes=0
         HwiresMPI%NumCurrentSegments=0
         HwiresMPI%NumNeededCurrentUpMPI=0
         HwiresMPI%NumNeededCurrentDownMPI=0
      endif

      !chequea los segmentos que estan en el padding de 1 celda

      NeedscontaMPIdown=0
      NeedscontaMPIup=0
      !
      if (layoutnumber/=size-1) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZE+1).and.(segmento%tipofield==iEz)) NeedscontaMPIup = NeedscontaMPIup + 1
         end do
      endif
      !
      if (layoutnumber/=0) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZI-1).and.(segmento%tipofield==iEz)) NeedscontaMPIdown = NeedscontaMPIdown + 1
         end do
      endif

      !
      SharescontaMPIdown=0
      SharescontaMPIup=0
      if (layoutnumber/=size-1) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZE).and.(segmento%tipofield==iEz)) SharescontaMPIup = SharescontaMPIup + 1
         end do
      endif
      !
      if (layoutnumber/=0) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZI).and.(segmento%tipofield==iEz)) SharescontaMPIdown = SharescontaMPIdown + 1
         end do
      endif


      allocate (HwiresMPI%MPIUpSharedCurrentSegment(1 : SharescontaMPIup))
      allocate (HwiresMPI%MPIDownSharedCurrentSegment(1 : SharescontaMPIdown))
      HwiresMPI%NumSharedCurrentUpMPI   = SharescontaMPIup !solo lo defino para info mia
      HwiresMPI%NumSharedCurrentDownMPI = SharescontaMPIdown

      !create space for the new ghost MPI segments (only their actual current is needed)
      !already allocated by the thin-wires routine

      if (.not.resume) then
         allocate (HwiresMPI%MPIUpNeededCurrentSegment(1 : NeedscontaMPIup))
         allocate (HwiresMPI%MPIDownNeededCurrentSegment(1 : NeedscontaMPIdown))
         HwiresMPI%NumNeededCurrentUpMPI=NeedscontaMPIup
         HwiresMPI%NumNeededCurrentDownMPI=NeedscontaMPIdown
         HwiresMPI%MPIUpNeededCurrentSegment( : )%Current   =0.0_RKIND
         HwiresMPI%MPIDownNeededCurrentSegment( : )%Current =0.0_RKIND
      else
         ! otherwise should have already been read by the thin-wires routine
         !but override the needed number
         NeedscontaMPIup  =HwiresMPI%NumNeededCurrentUpMPI
         NeedscontaMPIdown=HwiresMPI%NumNeededCurrentDownMPI
      endif

      NeedscontaMPIdown=0
      NeedscontaMPIup=0
      !
      if (layoutnumber/=size-1) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZE+1).and.(segmento%tipofield==iEz)) then
               NeedscontaMPIup = NeedscontaMPIup + 1
               HwiresMPI%MPIUpNeededCurrentSegment (NeedscontaMPIup)%equivalentIndex = i1
            endif
         end do
      endif
      !
      if (layoutnumber/=0) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZI-1).and.(segmento%tipofield==iEz)) then
               NeedscontaMPIdown = NeedscontaMPIdown + 1
               HwiresMPI%MPIDownNeededCurrentSegment (NeedscontaMPIdown)%equivalentIndex = i1
            endif
         end do
      endif

      SharescontaMPIdown=0
      SharescontaMPIup=0
      if (layoutnumber/=size-1) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZE).and.(segmento%tipofield==iEz)) then
               SharescontaMPIup = SharescontaMPIup + 1
               HwiresMPI%MPIUpSharedCurrentSegment(SharescontaMPIup)%equivalentIndex = i1
            endif
         end do
      endif
      !
      if (layoutnumber/=0) then
         do i1=1,HwiresMPI%NumCurrentSegments
            segmento =>HwiresMPI%CurrentSegment(i1)
            if ((segmento%k==C(iEz)%ZI).and.(segmento%tipofield==iEz)) then
               SharescontaMPIdown = SharescontaMPIdown + 1
               HwiresMPI%MPIDownSharedCurrentSegment(SharescontaMPIdown)%equivalentIndex = i1
            endif
         end do
      endif



      !allocate buffers
      Buffer%SendSizeUp=SharescontaMPIup
      Buffer%SendSizeDown=SharescontaMPIdown
      Buffer%RecSizeUp=NeedscontaMPIup
      Buffer%RecSizeDown=NeedscontaMPIdown

      allocate (Buffer%SendUp(1 : SharescontaMPIup),buffer%SendDown(1 : SharescontaMPIdown), &
      Buffer%RecUp(1 : NeedscontaMPIup),buffer%RecDown(1 : NeedscontaMPIdown))


      !reorder correctly to match the order the data is sent and received

      !allocate ibuffers
      iBuffer%SendSizeUp=   4*SharescontaMPIup
      iBuffer%SendSizeDown= 4*SharescontaMPIdown
      iBuffer%RecSizeUp=    4*NeedscontaMPIup
      iBuffer%RecSizeDown=  4*NeedscontaMPIdown

      allocate (iBuffer%SendUp(1 : iBuffer%SendSizeUp),ibuffer%SendDown(1 : iBuffer%SendSizeDown), &
      iBuffer%RecUp(1 : iBuffer%RecSizeUp),ibuffer%RecDown(1 : iBuffer%RecSizeDown))

      do i=1,SharescontaMPIup
         iBuffer%SendUp(4*i-3)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIUpSharedCurrentSegment(i)%EquivalentIndex)%i
         iBuffer%SendUp(4*i-2)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIUpSharedCurrentSegment(i)%EquivalentIndex)%j
         iBuffer%SendUp(4*i-1)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIUpSharedCurrentSegment(i)%EquivalentIndex)%k
         iBuffer%SendUp(4*i  )  =HwiresMPI%CurrentSegment(HwiresMPI%MPIUpSharedCurrentSegment(i)%EquivalentIndex)%origindex
      end do
      do i=1,SharescontaMPIDown
         iBuffer%SendDown(4*i-3)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIDownSharedCurrentSegment(i)%EquivalentIndex)%i
         iBuffer%SendDown(4*i-2)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIDownSharedCurrentSegment(i)%EquivalentIndex)%j
         iBuffer%SendDown(4*i-1)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIDownSharedCurrentSegment(i)%EquivalentIndex)%k
         iBuffer%SendDown(4*i  )  =HwiresMPI%CurrentSegment(HwiresMPI%MPIDownSharedCurrentSegment(i)%EquivalentIndex)%origindex
      end do

      CALL newFlushWiresMPIorigindexInfo(layoutnumber,size)

      do j=1,NeedscontaMPIdown
         ni       =IBuffer%RecDown(4*j-3)
         nj       =IBuffer%RecDown(4*j-2)
         nk       =IBuffer%RecDown(4*j-1)
         norigindex=IBuffer%RecDown(4*j  )
         busca: do i=1,NeedscontaMPIdown
            if ((ni == HwiresMPI%CurrentSegment(HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex)%i) .and. &
            (nj == HwiresMPI%CurrentSegment(HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex)%j) .and. &
            (nk == HwiresMPI%CurrentSegment(HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex)%k) .and. &
            (norigindex == HwiresMPI%CurrentSegment(HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex)%origindex)) then
               idum=HwiresMPI%MPIDownNeededCurrentSegment(j)%EquivalentIndex !swap indexes
               HwiresMPI%MPIDownNeededCurrentSegment(j)%EquivalentIndex=HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex
               HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex=idum
               exit busca
            endif
         end do busca
      end do

      do j=1,NeedscontaMPIUp
         ni       =IBuffer%RecUp(4*j-3)
         nj       =IBuffer%RecUp(4*j-2)
         nk       =IBuffer%RecUp(4*j-1)
         norigindex=IBuffer%RecUp(4*j  )
         busca2: do i=1,NeedscontaMPIUp
            if ((ni == HwiresMPI%CurrentSegment(HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex)%i) .and. &
            (nj == HwiresMPI%CurrentSegment(HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex)%j) .and. &
            (nk == HwiresMPI%CurrentSegment(HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex)%k) .and. &
            (norigindex == HwiresMPI%CurrentSegment(HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex)%origindex)) then
               idum=HwiresMPI%MPIUpNeededCurrentSegment(j)%EquivalentIndex !swap indexes
               HwiresMPI%MPIUpNeededCurrentSegment(j)%EquivalentIndex=HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex
               HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex=idum
               exit busca2
            endif
         end do busca2
      end do
      deallocate (iBuffer%SendUp,ibuffer%SendDown,iBuffer%RecUp,ibuffer%RecDown)


      return
   end subroutine newInitWiresMPI


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FLUSH WIRE DATA
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine newFlushWiresMPIorigindexInfo(layoutnumber,size)
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  ierr1=0,ierr2=0,ierr3=0,ierr4=0,ierr5=0,ierr6=0,ierr7=0,ierr8=0,ierr9=0,ierr10=0,ierr11=0,ierr12=0
      integer  ::  status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)
      integer  ::  req1,req2,req11,req21
      character (LEN=BUFSIZE)  ::  whoami
      character(len=BUFSIZE) :: buff
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '


      ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;

      !
      if ((layoutnumber/=size-1).and.(ibuffer%RecSizeUp/=0)) & !syncUp
      call MPI_IRECV (ibuffer%RecUp(1),  &
      ibuffer%RecSizeUp,  MPI_INTEGER, layoutnumber+1_4, 1_4, SUBCOMM_MPI, req1,  ierr5)
      if ((layoutnumber/=size-1).and.(ibuffer%SendSizeUp/=0)) & !syncUp
      call MPI_ISEND (ibuffer%SendUp(1), &
      ibuffer%SendSizeUp, MPI_INTEGER, layoutnumber+1_4, 2_4, SUBCOMM_MPI, req11,  ierr6)
      if ((layoutnumber/=0    ).and.(ibuffer%SendSizeDown/=0)) & !syncDown
      call MPI_ISEND (ibuffer%SendDown(1), &
      ibuffer%SendSizeDown, MPI_INTEGER, layoutnumber-1_4, 1_4, SUBCOMM_MPI, req2,  ierr7)
      if ((layoutnumber/=0    ).and.(ibuffer%RecSizeDown/=0)) & !syncDown
      call MPI_IRECV (ibuffer%RecDown(1), &
      ibuffer%RecSizeDown, MPI_INTEGER, layoutnumber-1_4, 2_4, SUBCOMM_MPI, req21,  ierr8)

      !
      if ((layoutnumber/=size-1).and.(ibuffer%RecSizeUp/=0))    call MPI_WAIT(req1 ,status1,ierr9)
      if ((layoutnumber/=size-1).and.(ibuffer%SendSizeUp/=0))   call MPI_WAIT(req11,status1,ierr10)
      if ((layoutnumber/=0     ).and.(ibuffer%SendSizeDown/=0)) call MPI_WAIT(req2 ,status2,ierr10)
      if ((layoutnumber/=0     ).and.(ibuffer%RecSizeDown/=0))  call MPI_WAIT(req21,status2,ierr10)
      !

      !
      !call MPI_Barrier(SUBCOMM_MPI,ierr12)
      !ojo que aqui no entran todos y por tanto la barrera crea un deadlock
      !
      if ((layoutnumber/=0    ).and.(ierr1+ierr2+ierr3+ierr4 /= 0)) then
         write(buff,*) 'FLUSHMPI ierr1,ierr2,ierr3,ierr4',LAYOUTNUMBER+1_4,ierr1,ierr2,ierr3,ierr4
         call stoponerror (layoutnumber,size,buff)
      endif
      if ((layoutnumber/=size-1).and.(ierr5+ierr6+ierr7+ierr8 /= 0)) then
         write(buff,*) 'FLUSHMPI ierr5,ierr6,ierr7,ierr8',LAYOUTNUMBER+1_4,ierr5,ierr6,ierr7,ierr8
         call stoponerror (layoutnumber,size,buff)
      endif
      if (ierr9+ierr10+ierr11+ierr12 /= 0) then
         write(buff,*) 'FLUSHMPI ierr9,ierr10,ierr11,ierr12',LAYOUTNUMBER+1_4,ierr9,ierr10,ierr11,ierr12
         call stoponerror (layoutnumber,size,buff)
      endif
      return
   end subroutine newFlushWiresMPIorigindexInfo

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FLUSH WIRE DATA
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine newFlushWiresMPI(layoutnumber,size)
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  ierr1=0,ierr2=0,ierr3=0,ierr4=0,ierr5=0,ierr6=0,ierr7=0,ierr8=0,ierr9=0,ierr10=0,ierr11=0,ierr12=0,i
      integer ::  status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)
      integer ::  req1,req2,req11,req21
      character (LEN=BUFSIZE)  ::  whoami
      character(len=BUFSIZE) :: buff
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '


      ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;

      do i=1,Buffer%SendSizeUP
         Buffer%SendUp(i)  =HwiresMPI%CurrentSegment(HwiresMPI%MPIUpSharedCurrentSegment(i)%EquivalentIndex)%Current
      end do
      do i=1,Buffer%SendSizeDown
         Buffer%SendDown(i)=HwiresMPI%CurrentSegment(HwiresMPI%MPIDownSharedCurrentSegment(i)%EquivalentIndex)%Current
      end do

      !
      if ((layoutnumber/=size-1).and.(Buffer%RecSizeUp/=0)) & !syncUp
      call MPI_IRECV (Buffer%RecUp(1),  &
      Buffer%RecSizeUp, REALSIZE_wires, layoutnumber+1_4, 1_4, SUBCOMM_MPI, req1,  ierr5)
      if ((layoutnumber/=size-1).and.(Buffer%SendSizeUp/=0)) & !syncUp
      call MPI_ISEND (Buffer%SendUp(1), &
      Buffer%SendSizeUp, REALSIZE_wires, layoutnumber+1_4, 2_4, SUBCOMM_MPI, req11,  ierr6)
      if ((layoutnumber/=0    ).and.(Buffer%SendSizeDown/=0)) & !syncDown
      call MPI_ISEND (Buffer%SendDown(1), &
      Buffer%SendSizeDown, REALSIZE_wires, layoutnumber-1_4, 1_4, SUBCOMM_MPI, req2,  ierr7)
      if ((layoutnumber/=0    ).and.(Buffer%RecSizeDown/=0)) & !syncDown
      call MPI_IRECV (Buffer%RecDown(1), &
      Buffer%RecSizeDown, REALSIZE_wires, layoutnumber-1_4, 2_4, SUBCOMM_MPI, req21,  ierr8)

      !
      if ((layoutnumber/=size-1).and.(Buffer%RecSizeUp/=0))    call MPI_WAIT(req1 ,status1,ierr9)
      if ((layoutnumber/=size-1).and.(Buffer%SendSizeUp/=0))   call MPI_WAIT(req11,status1,ierr10)
      if ((layoutnumber/=0     ).and.(Buffer%SendSizeDown/=0)) call MPI_WAIT(req2 ,status2,ierr10)
      if ((layoutnumber/=0     ).and.(Buffer%RecSizeDown/=0))  call MPI_WAIT(req21,status2,ierr10)
      !

      do i=1,Buffer%RecSizeDown
         HwiresMPI%CurrentSegment(HwiresMPI%MPIDownNeededCurrentSegment(i)%EquivalentIndex)%Current=Buffer%RecDown(i)
      end do
      do i=1,Buffer%RecSizeUp
         HwiresMPI%CurrentSegment(HwiresMPI%MPIUpNeededCurrentSegment(i)%EquivalentIndex)%Current  =Buffer%RecUp(i)
      end do
      !
      !call MPI_Barrier(SUBCOMM_MPI,ierr12)
      !ojo que aqui no entran todos y por tanto la barrera crea un deadlock
      !
      if ((layoutnumber/=0    ).and.(ierr1+ierr2+ierr3+ierr4 /= 0)) then
         write(buff,*) 'FLUSHMPI ierr1,ierr2,ierr3,ierr4',LAYOUTNUMBER+1_4,ierr1,ierr2,ierr3,ierr4
         call stoponerror (layoutnumber,size,buff)
      endif
      if ((layoutnumber/=size-1).and.(ierr5+ierr6+ierr7+ierr8 /= 0)) then
         write(buff,*) 'FLUSHMPI ierr5,ierr6,ierr7,ierr8',LAYOUTNUMBER+1_4,ierr5,ierr6,ierr7,ierr8
         call stoponerror (layoutnumber,size,buff)
      endif
      if (ierr9+ierr10+ierr11+ierr12 /= 0) then
         write(buff,*) 'FLUSHMPI ierr9,ierr10,ierr11,ierr12',LAYOUTNUMBER+1_4,ierr9,ierr10,ierr11,ierr12
         call stoponerror (layoutnumber,size,buff)
      endif
      return
   end subroutine newFlushWiresMPI



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FLUSH WIRE additional info DATA
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine FlushWiresMPIorigindexInfo(layoutnumber,size)
      integer (kind=4), intent(in) :: layoutnumber,size
      integer (kind=4)  ::  ierr1=0,ierr2=0,ierr3=0,ierr4=0,ierr5=0,ierr6=0,ierr7=0,ierr8=0,ierr9=0,ierr10=0,ierr11=0,ierr12=0,i
      integer :: status1(MPI_STATUS_SIZE),status2(MPI_STATUS_SIZE)
      integer :: req1,req2,req11,req21
      character (LEN=BUFSIZE)  ::  whoami
      character(len=BUFSIZE) :: buff
      write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '


      ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;

      do i=1,Buffer%SendSizeUP
         Buffer%SendUp(i)  =HwiresMPI%MPIUpChargeNode(i)%MPIsharedCurrent%origindex * 1.0_RKIND_wires
      end do
      do i=1,Buffer%SendSizeDown
         Buffer%SendDown(i)=HwiresMPI%MPIDownChargeNode(i)%MPIsharedCurrent%origindex * 1.0_RKIND_wires
      end do

      !
      if ((layoutnumber/=size-1).and.(Buffer%RecSizeUp/=0)) & !syncUp
      call MPI_IRECV (Buffer%RecUp(1),  &
      Buffer%RecSizeUp, REALSIZE_wires, layoutnumber+1_4, 1_4, SUBCOMM_MPI, req1,  ierr5)
      if ((layoutnumber/=size-1).and.(Buffer%SendSizeUp/=0)) & !syncUp
      call MPI_ISEND (Buffer%SendUp(1), &
      Buffer%SendSizeUp, REALSIZE_wires, layoutnumber+1_4, 2_4, SUBCOMM_MPI, req11,  ierr6)
      if ((layoutnumber/=0    ).and.(Buffer%SendSizeDown/=0)) & !syncDown
      call MPI_ISEND (Buffer%SendDown(1), &
      Buffer%SendSizeDown, REALSIZE_wires, layoutnumber-1_4, 1_4, SUBCOMM_MPI, req2,  ierr7)
      if ((layoutnumber/=0    ).and.(Buffer%RecSizeDown/=0)) & !syncDown
      call MPI_IRECV (Buffer%RecDown(1), &
      Buffer%RecSizeDown, REALSIZE_wires, layoutnumber-1_4, 2_4, SUBCOMM_MPI, req21,  ierr8)

      !
      if ((layoutnumber/=size-1).and.(Buffer%RecSizeUp/=0))    call MPI_WAIT(req1 ,status1,ierr9)
      if ((layoutnumber/=size-1).and.(Buffer%SendSizeUp/=0))   call MPI_WAIT(req11,status1,ierr10)
      if ((layoutnumber/=0     ).and.(Buffer%SendSizeDown/=0)) call MPI_WAIT(req2 ,status2,ierr10)
      if ((layoutnumber/=0     ).and.(Buffer%RecSizeDown/=0))  call MPI_WAIT(req21,status2,ierr10)
      !

      do i=1,Buffer%RecSizeDown
         HwiresMPI%MPIDownNeededCurrentSegment(i)%origindex=nINT(Buffer%RecDown(i))
      end do
      do i=1,Buffer%RecSizeUp
         HwiresMPI%MPIUpNeededCurrentSegment(i)%origindex  =nINT(Buffer%RecUp(i))
      end do
      !
      !call MPI_Barrier(SUBCOMM_MPI,ierr12)
      !ojo que aqui no entran todos y por tanto la barrera crea un deadlock
      !
      if ((layoutnumber/=0    ).and.(ierr1+ierr2+ierr3+ierr4 /= 0)) then
         write(buff,*) 'FLUSHMPI ierr1,ierr2,ierr3,ierr4',LAYOUTNUMBER+1_4,ierr1,ierr2,ierr3,ierr4
         call stoponerror (layoutnumber,size,buff)
      endif
      if ((layoutnumber/=size-1).and.(ierr5+ierr6+ierr7+ierr8 /= 0)) then
         write(buff,*) 'FLUSHMPI ierr5,ierr6,ierr7,ierr8',LAYOUTNUMBER+1_4,ierr5,ierr6,ierr7,ierr8
         call stoponerror (layoutnumber,size,buff)
      endif
      if (ierr9+ierr10+ierr11+ierr12 /= 0) then
         write(buff,*) 'FLUSHMPI ierr9,ierr10,ierr11,ierr12',LAYOUTNUMBER+1_4,ierr9,ierr10,ierr11,ierr12
         call stoponerror (layoutnumber,size,buff)
      endif
      return
   end subroutine

   !!!!!!!!!


   !deprecated
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!! FLUSH WIRE DATA
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!subroutine FlushWiresMPI(layoutnumber,size)
   !!!integer (kind=4), intent(in) :: layoutnumber,size
   !!!integer (kind=4)  ::  ierr1=0,ierr2=0,ierr3=0,ierr4=0,ierr5=0,ierr6=0,ierr7=0,ierr8=0,ierr9=0,ierr10=0,ierr11=0,ierr12=0
   !!!integer (kind=4)  ::  status1(MPI_STATUS_SIZE,1 : 2),status2(MPI_STATUS_SIZE,1 : 2),i
   !!!integer (kind=4)  ::  req1,req2,req11,req21
   !!!character (LEN=BUFSIZE)  ::  whoami
   !!!write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
   !!!
   !!!
   !!!ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;
   !!!
   !!!do i=1,Buffer%SendSizeUP
   !!!   Buffer%SendUp(i)  =HwiresMPI%MPIUpChargeNode(i)%MPIsharedCurrent%Current
   !!!end do
   !!!do i=1,Buffer%SendSizeDown
   !!!    Buffer%SendDown(i)=HwiresMPI%MPIDownChargeNode(i)%MPIsharedCurrent%Current
   !!!end do
   !!!
   !!!!
   !!!if ((layoutnumber/=size-1).and.(Buffer%RecSizeUp/=0)) & !syncUp
   !!!    call MPI_IRECV (Buffer%RecUp(1),  &
   !!!                    Buffer%RecSizeUp, REALSIZE_wires, layoutnumber+1_4, 1_4, SUBCOMM_MPI, req1,  ierr5)
   !!!if ((layoutnumber/=size-1).and.(Buffer%SendSizeUp/=0)) & !syncUp
   !!!    call MPI_ISEND (Buffer%SendUp(1), &
   !!!                    Buffer%SendSizeUp, REALSIZE_wires, layoutnumber+1_4, 2_4, SUBCOMM_MPI, req11,  ierr6)
   !!!if ((layoutnumber/=0    ).and.(Buffer%SendSizeDown/=0)) & !syncDown
   !!!    call MPI_ISEND (Buffer%SendDown(1), &
   !!!                    Buffer%SendSizeDown, REALSIZE_wires, layoutnumber-1_4, 1_4, SUBCOMM_MPI, req2,  ierr7)
   !!!if ((layoutnumber/=0    ).and.(Buffer%RecSizeDown/=0)) & !syncDown
   !!!    call MPI_IRECV (Buffer%RecDown(1), &
   !!!                    Buffer%RecSizeDown, REALSIZE_wires, layoutnumber-1_4, 2_4, SUBCOMM_MPI, req21,  ierr8)
   !!!
   !!!!
   !!!if ((layoutnumber/=size-1).and.(Buffer%RecSizeUp/=0))    call MPI_WAIT(req1 ,status1,ierr9)
   !!!if ((layoutnumber/=size-1).and.(Buffer%SendSizeUp/=0))   call MPI_WAIT(req11,status1,ierr10)
   !!!if ((layoutnumber/=0     ).and.(Buffer%SendSizeDown/=0)) call MPI_WAIT(req2 ,status2,ierr10)
   !!!if ((layoutnumber/=0     ).and.(Buffer%RecSizeDown/=0))  call MPI_WAIT(req21,status2,ierr10)
   !!!!
   !!!
   !!!do i=1,Buffer%RecSizeDown
   !!!!   print *,whoami,'setting down',i,'ijk', &
   !!!!   HwiresMPI%MPIDownNeededCurrentSegment(i)%i, &
   !!!!   HwiresMPI%MPIDownNeededCurrentSegment(i)%j, &
   !!!!   HwiresMPI%MPIDownNeededCurrentSegment(i)%k, &
   !!!!   HwiresMPI%MPIDownNeededCurrentSegment(i)%tipofield, 'to',Buffer%RecDown(i)
   !!!    HwiresMPI%MPIDownNeededCurrentSegment(i)%Current=Buffer%RecDown(i)
   !!!end do
   !!!do i=1,Buffer%RecSizeUp
   !!!!   print *,whoami,'setting up',i,'ijk', &
   !!!!   HwiresMPI%MPIUpNeededCurrentSegment(i)%i, &
   !!!!   HwiresMPI%MPIUpNeededCurrentSegment(i)%j, &
   !!!!   HwiresMPI%MPIUpNeededCurrentSegment(i)%k, &
   !!!!   HwiresMPI%MPIUpNeededCurrentSegment(i)%tipofield,' to',Buffer%RecUp(i)
   !!!    HwiresMPI%MPIUpNeededCurrentSegment(i)%Current  =Buffer%RecUp(i)
   !!!end do
   !!!!
   !!!!call MPI_Barrier(SUBCOMM_MPI,ierr12)
   !!!!ojo que aqui no entran todos y por tanto la barrera crea un deadlock
   !!!!
   !!!if ((layoutnumber/=0    ).and.(ierr1+ierr2+ierr3+ierr4 /= 0)) &
   !!!      call StopOnError(layoutnumber,size,'FLUSHMPI ierr1,ierr2,ierr3,ierr4',LAYOUTNUMBER+1_4,ierr1,ierr2,ierr3,ierr4)
   !!!if ((layoutnumber/=size-1).and.(ierr5+ierr6+ierr7+ierr8 /= 0)) &
   !!!      call StopOnError(layoutnumber,size,'FLUSHMPI ierr5,ierr6,ierr7,ierr8',LAYOUTNUMBER+1_4,ierr5,ierr6,ierr7,ierr8)
   !!!if (ierr9+ierr10+ierr11+ierr12 /= 0) &
   !!!      call StopOnError(layoutnumber,size,'FLUSHMPI ierr9,ierr10,ierr11,ierr12',LAYOUTNUMBER+1_4,ierr9,ierr10,ierr11,ierr12)
   !!!return
   !!!end subroutine FlushWiresMPI
   !!!!
   !!!!


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! new routine: works without the MediaMatrix Info
   ! supports multiwires
   !ADJUSTS THE WIRE DATA NEEDED FOR TRANSVER
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! FLUSH WIRE DATA
   !! no se usa. detectado y comentado 260815 bug gra_simple.nfde en sync mpi de hilos que se doblan sobre si mismos y comparten tramos coincidentes
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!!subroutine newFlushWiresMPIindexmedInfo(layoutnumber,size)
   !!!   integer (kind=4)  ::  ierr1=0,ierr2=0,ierr3=0,ierr4=0,ierr5=0,ierr6=0,ierr7=0,ierr8=0,ierr9=0,ierr10=0,ierr11=0,ierr12=0
   !!!   integer (kind=4)  ::  status1(MPI_STATUS_SIZE,1 : 2),status2(MPI_STATUS_SIZE,1 : 2)
   !!!   integer (kind=4), intent(in) :: layoutnumber,size
   !!!   integer (kind=4)  ::  req1,req2,req11,req21
   !!!   character (LEN=BUFSIZE)  ::  whoami
   !!!   character(len=BUFSIZE) :: buff
   !!!   write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',size,') '
   !!!
   !!!
   !!!   ierr1=0;ierr2=0;ierr3=0;ierr4=0;ierr5=0;ierr6=0;ierr7=0;ierr8=0;ierr9=0;ierr10=0;ierr11=0;ierr12=0;
   !!!
   !!!   !
   !!!   if ((layoutnumber/=size-1).and.(ibuffer%RecSizeUp/=0)) & !syncUp
   !!!   call MPI_IRECV (ibuffer%RecUp(1),  &
   !!!   ibuffer%RecSizeUp,  MPI_INTEGER, layoutnumber+1_4, 1_4, SUBCOMM_MPI, req1,  ierr5)
   !!!   if ((layoutnumber/=size-1).and.(ibuffer%SendSizeUp/=0)) & !syncUp
   !!!   call MPI_ISEND (ibuffer%SendUp(1), &
   !!!   ibuffer%SendSizeUp, MPI_INTEGER, layoutnumber+1_4, 2_4, SUBCOMM_MPI, req11,  ierr6)
   !!!   if ((layoutnumber/=0    ).and.(ibuffer%SendSizeDown/=0)) & !syncDown
   !!!   call MPI_ISEND (ibuffer%SendDown(1), &
   !!!   ibuffer%SendSizeDown, MPI_INTEGER, layoutnumber-1_4, 1_4, SUBCOMM_MPI, req2,  ierr7)
   !!!   if ((layoutnumber/=0    ).and.(ibuffer%RecSizeDown/=0)) & !syncDown
   !!!   call MPI_IRECV (ibuffer%RecDown(1), &
   !!!   ibuffer%RecSizeDown, MPI_INTEGER, layoutnumber-1_4, 2_4, SUBCOMM_MPI, req21,  ierr8)
   !!!
   !!!   !
   !!!   if ((layoutnumber/=size-1).and.(ibuffer%RecSizeUp/=0))    call MPI_WAIT(req1 ,status1,ierr9)
   !!!   if ((layoutnumber/=size-1).and.(ibuffer%SendSizeUp/=0))   call MPI_WAIT(req11,status1,ierr10)
   !!!   if ((layoutnumber/=0     ).and.(ibuffer%SendSizeDown/=0)) call MPI_WAIT(req2 ,status2,ierr10)
   !!!   if ((layoutnumber/=0     ).and.(ibuffer%RecSizeDown/=0))  call MPI_WAIT(req21,status2,ierr10)
   !!!   !
   !!!
   !!!   !
   !!!   !call MPI_Barrier(SUBCOMM_MPI,ierr12)
   !!!   !ojo que aqui no entran todos y por tanto la barrera crea un deadlock
   !!!   !
   !!!   if ((layoutnumber/=0    ).and.(ierr1+ierr2+ierr3+ierr4 /= 0)) then
   !!!      write(buff,*) 'FLUSHMPI ierr1,ierr2,ierr3,ierr4',LAYOUTNUMBER+1_4,ierr1,ierr2,ierr3,ierr4
   !!!      call stoponerror (layoutnumber,size,buff)
   !!!   endif
   !!!   if ((layoutnumber/=size-1).and.(ierr5+ierr6+ierr7+ierr8 /= 0)) then
   !!!      write(buff,*) 'FLUSHMPI ierr5,ierr6,ierr7,ierr8',LAYOUTNUMBER+1_4,ierr5,ierr6,ierr7,ierr8
   !!!      call stoponerror (layoutnumber,size,buff)
   !!!   endif
   !!!   if (ierr9+ierr10+ierr11+ierr12 /= 0) then
   !!!      write(buff,*) 'FLUSHMPI ierr9,ierr10,ierr11,ierr12',LAYOUTNUMBER+1_4,ierr9,ierr10,ierr11,ierr12
   !!!      call stoponerror (layoutnumber,size,buff)
   !!!   endif
   !!!   return
   !!!end subroutine newFlushWiresMPIindexmedInfo

   subroutine InitExtraFlushMPI (layoutnumber,sggsweep,sggalloc,med,nummed,sggmiez,sggMiHz)
      type (XYZlimit_t), dimension(1:6)  ::  sggalloc,sggsweep
      integer (kind=4) :: layoutnumber
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES), intent(in)   :: &
      sggMiEz(sggalloc(iEz)%XI : sggalloc(iEz)%XE,sggalloc(iEz)%YI : sggalloc(iEz)%YE,sggalloc(iEz)%ZI : sggalloc(iEz)%ZE), &
      sggMiHz(sggalloc(iHz)%XI : sggalloc(iHz)%XE,sggalloc(iHz)%YI : sggalloc(iHz)%YE,sggalloc(iHz)%ZI : sggalloc(iHz)%ZE)

      integer (kind=4)  ::  j1,i1,jmed
      integer (kind=4) , intent(in) :: nummed
      type (MediaData_t), dimension ( 0:NumMed ), intent(in) :: med

      !thanks to overlapping of media matrix to detect media
      !an anistropic boundary will be detected both by the upper and the lower layouts
      !regardless of what happened with the wires (creo que esta formulacion no pasaria tampoco con
      !los wires, pero no lo he tocado)
      FlushExtraInfoDown=.false.
      FlushExtraInfoUp  =.false.

      !is enough to check Ez and Hz

      do jmed=1,NumMed
         if( Med( jmed)%Is%Anisotropic) then
            !!!Ez
            Do j1=sggsweep(iEz)%YI,sggsweep(iEz)%YE
               Do i1=sggsweep(iEz)%XI,sggsweep(iEz)%XE
                  if ((sggMiEz(i1,j1,   comZ)) == jmed)   then
                     FlushExtraInfoDown=.true.
                  endif
                  if ((sggMiEz(i1,j1,-1+comZ)) == jmed)   then
                     FlushExtraInfoDown=.true.
                  endif
                  if ((sggMiEz(i1,j1,   finZ)) == jmed)   then
                     FlushExtraInfoUp  =.true.
                  endif
                  if ((sggMiEz(i1,j1, 1+finZ)) == jmed)   then
                     FlushExtraInfoUp  =.true.
                  endif
               end do
            end do
            !!!Hz
            Do j1=sggsweep(iHz)%YI,sggsweep(iHz)%YE
               Do i1=sggsweep(iHz)%XI,sggsweep(iHz)%XE
                  if ((sggMiHz(i1,j1,   comZ)) == jmed)   then
                     FlushExtraInfoDown=.true.
                  endif
                  if ((sggMiHz(i1,j1, 1+comZ)) == jmed)   then
                     FlushExtraInfoDown=.true.
                  endif
                  if ((sggMiHz(i1,j1, 1+finZ)) == jmed)  then
                      FlushExtraInfoUp  =.true.
                  endif
                  if ((sggMiHz(i1,j1, 2+finZ)) == jmed)  then
                      FlushExtraInfoUp  =.true.
                  endif
               end do
            end do
         endif
         if(Med( jmed)%Is%SGBC .or. Med( jmed)%Is%Multiport .or. Med( jmed)%Is%AnisMultiport ) then
            !!!Hz
            Do j1=sggsweep(iHz)%YI,sggsweep(iHz)%YE
               Do i1=sggsweep(iHz)%XI,sggsweep(iHz)%XE
                  if ((sggMiHz(i1,j1,   comZ)) == jmed) then
                     FlushExtraInfoDown  = .true.
                  endif
                  if ((sggMiHz(i1,j1, 1+finZ)) == jmed) then
                       FlushExtraInfoUp    = .true.
                  endif
                  !creo que esto no es necesario para multiports de ss pero no creo que cargue mucho y no se si Ian lo necesita
                  !lo dejo por precaucion
                  if ((sggMiHz(i1,j1, 1+comZ)) == jmed) then
                       FlushExtraInfoDown  = .true.
                  endif
                  if ((sggMiHz(i1,j1, 2+finZ)) == jmed) then
                       FlushExtraInfoUp    = .true.
                  endif
               end do
            end do
         endif
      end do

!      print *,'----> layout,FlushExtraInfoDown,FlushExtraInfoUp',layoutnumber,FlushExtraInfoDown,FlushExtraInfoUp
      return

   end subroutine InitExtraFlushMPI



   !!!!!!!!Cray stuff



   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! INITIALIZES THE BOUNDS AND SIZES TO COMMUNICATE
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine InitMPI_Cray( layoutnumber, size, sggsweep, sggalloc, PBCDown, PBCUp, &
   Ex, Ey, Ez, Hx, Hy, Hz)
      !---------------- inputs -----------------------------------------------------------------------
      type( XYZlimit_t), dimension( 1: 6), intent( IN) :: sggalloc, sggsweep
      real (kind = RKIND), intent( IN), target :: &
      Hx( sggalloc( iHx)%XI: sggalloc( iHx)%XE, sggalloc( iHx)%YI: sggalloc( iHx)%YE, sggalloc( iHx)%ZI: sggalloc( iHx)%ZE), &
      Hy( sggalloc( iHy)%XI: sggalloc( iHy)%XE, sggalloc( iHy)%YI: sggalloc( iHy)%YE, sggalloc( iHy)%ZI: sggalloc( iHy)%ZE), &
      Hz( sggalloc( iHz)%XI: sggalloc( iHz)%XE, sggalloc( iHz)%YI: sggalloc( iHz)%YE, sggalloc( iHz)%ZI: sggalloc( iHz)%ZE)
      REAL (KIND=RKIND)   , intent( IN),target :: &
      Ex( sggalloc( iEx)%XI: sggalloc( iEx)%XE, sggalloc( iEx)%YI: sggalloc( iEx)%YE, sggalloc( iEx)%ZI: sggalloc( iEx)%ZE), &
      Ey( sggalloc( iEy)%XI: sggalloc( iEy)%XE, sggalloc( iEy)%YI: sggalloc( iEy)%YE, sggalloc( iEy)%ZI: sggalloc( iEy)%ZE), &
      Ez( sggalloc( iEz)%XI: sggalloc( iEz)%XE, sggalloc( iEz)%YI: sggalloc( iEz)%YE, sggalloc( iEz)%ZI: sggalloc( iEz)%ZE)
      integer( kind = 4), intent( IN) :: layoutnumber, size
      logical, intent( IN) :: PBCDown, PBCUp
      !---------------- variables locales ------------------------------------------------------------
      type( t_databuf), pointer :: databufH, databufE
      !---------------- empieza InitMPI_Cray ---------------------------------------------------------
      !not necessary at this moment since nothing is read at this mmoment
      ExXI=sggalloc(iEx)%XI
      ExXE=sggalloc(iEx)%XE
      EyXI=sggalloc(iEy)%XI
      EyXE=sggalloc(iEy)%XE
      EzXI=sggalloc(iEz)%XI
      EzXE=sggalloc(iEz)%XE
      !--->
      ExYI=sggalloc(iEx)%YI
      ExYE=sggalloc(iEx)%YE
      EyYI=sggalloc(iEy)%YI
      EyYE=sggalloc(iEy)%YE
      EzYI=sggalloc(iEz)%YI
      EzYE=sggalloc(iEz)%YE
      !--->
      HxXI=sggalloc(iHx)%XI
      HxXE=sggalloc(iHx)%XE
      HyXI=sggalloc(iHy)%XI
      HyXE=sggalloc(iHy)%XE
      HzXI=sggalloc(iHz)%XI
      HzXE=sggalloc(iHz)%XE
      !--->
      HxYI=sggalloc(iHx)%YI
      HxYE=sggalloc(iHx)%YE
      HyYI=sggalloc(iHy)%YI
      HyYE=sggalloc(iHy)%YE
      HzYI=sggalloc(iHz)%YI
      HzYE=sggalloc(iHz)%YE
      !--->
      sizeEx=(ExXE-ExXI+1)*(ExYE-ExYI+1)
      sizeEy=(EyXE-EyXI+1)*(EyYE-EyYI+1)
      sizeEz=(EzXE-EzXI+1)*(EzYE-EzYI+1)
      !--->
      sizeHx=(HxXE-HxXI+1)*(HxYE-HxYI+1)
      sizeHy=(HyXE-HyXI+1)*(HyYE-HyYI+1)
      sizeHz=(HzXE-HzXI+1)*(HzYE-HzYI+1)
      !--->
      ComZ=sggsweep(iHx)%ZI !both Hx and Hy coincide in this
      FinZ=sggsweep(iHx)%ZE
      !--->
      !jag: bug Antares mas de 65295 steps
      databuf_SetH%syncUp = layoutnumber /= (size-1)
      databuf_SetH%pbcUp = (layoutnumber == (size-1)) .and. PBCUp
      databufH => databuf_SetH%databuf_Up
      databufE => databuf_SetE%databuf_Up
      databufH%FlushExtraInfo = .FALSE. !later overriden by the anisotropic MPI init
      if( databuf_SetH%syncUp .OR. databuf_SetH%pbcUp) then
         databufH%sizex = sizeHx
         databufH%sizey = sizeHy
         databufH%sizez = -1
         !--->
         databufH%buf_x_rx => Hx( HxXI: HxXE, HxYI: HxYE, finZ+1)
         databufH%buf_x_tx => Hx( HxXI: HxXE, HxYI: HxYE, finZ  )
         databufH%buf_y_rx => Hy( HyXI: HyXE, HyYI: HyYE, finZ+1)
         databufH%buf_y_tx => Hy( HyXI: HyXE, HyYI: HyYE, finZ  )
         databufH%buf_z_rx => null( )
         databufH%buf_z_tx => null( )
         !--->
         databufE%sizex = -1
         databufE%sizey = -1
         databufE%sizez = -1
         !--->
         databufE%buf_z_rx => null( )
         databufE%buf_z_tx => null( )
         databufE%buf_x_rx => null( )
         databufE%buf_x_tx => null( )
         databufE%buf_y_rx => null( )
         databufE%buf_y_tx => null( )
         !--->
         if( databuf_SetH%pbcUp) then
            databufH%ip_target = 0
            databufE%ip_target = 0
         else
            databufH%ip_target = layoutnumber + 1
            databufE%ip_target = layoutnumber + 1
         endif
      else
         databufH%ip_target = -1
         !--->
         databufH%sizex = -1
         databufH%sizey = -1
         databufH%sizez = -1
         !--->
         databufH%buf_x_tx => null( )
         databufH%buf_y_tx => null( )
         databufH%buf_x_rx => null( )
         databufH%buf_y_rx => null( )
         databufH%buf_z_tx => null( )
         databufH%buf_z_rx => null( )
         !--->
         databufE%ip_target = -1
         !--->
         databufE%sizex = -1
         databufE%sizey = -1
         databufE%sizez = -1
         !--->
         databufE%buf_z_rx => null( )
         databufE%buf_z_tx => null( )
         databufE%buf_x_rx => null( )
         databufE%buf_x_tx => null( )
         databufE%buf_y_rx => null( )
         databufE%buf_y_tx => null( )
      endif
      !-----------------------------------------------------> DW
      databuf_SetH%syncDown = layoutnumber/=0
      databuf_SetH%pbcDown = (layoutnumber==0).and. PBCDown
      databufH => databuf_SetH%databuf_Down
      databufE => databuf_SetE%databuf_Down
      databufH%FlushExtraInfo = .FALSE. !later overriden by the anisotropic MPI init
      if( databuf_SetH%syncDown .OR. databuf_SetH%pbcDown) then
         databufH%sizex = sizeHx
         databufH%sizey = sizeHy
         databufH%sizez = -1
         !--->
         databufH%buf_x_tx => Hx( HxXI: HxXE, HxYI: HxYE, comZ)
         databufH%buf_x_rx => Hx( HxXI: HxXE, HxYI: HxYE, comZ-1)
         databufH%buf_y_tx => Hy( HyXI: HyXE, HyYI: HyYE, comZ)
         databufH%buf_y_rx => Hy( HyXI: HyXE, HyYI: HyYE, comZ-1)
         databufH%buf_z_rx => null( )
         databufH%buf_z_tx => null( )
         !--->
         databufE%sizex = -1
         databufE%sizey = -1
         databufE%sizez = -1
         !--->
         databufE%buf_z_rx => null( )
         databufE%buf_z_tx => null( )
         databufE%buf_x_rx => null( )
         databufE%buf_x_tx => null( )
         databufE%buf_y_rx => null( )
         databufE%buf_y_tx => null( )
         !--->
         if( databuf_SetH%pbcDown) then
            databufH%ip_target = size-1
            databufE%ip_target = size-1
         else
            databufH%ip_target = layoutnumber - 1
            databufE%ip_target = layoutnumber - 1
         endif
      else
         databufH%ip_target = -1
         !--->
         databufH%sizex = -1
         databufH%sizey = -1
         databufH%sizez = -1
         !--->
         databufH%buf_x_tx => null( )
         databufH%buf_y_tx => null( )
         databufH%buf_x_rx => null( )
         databufH%buf_y_rx => null( )
         databufH%buf_z_tx => null( )
         databufH%buf_z_rx => null( )
         !--->
         databufE%ip_target = -1
         !--->
         databufE%sizex = -1
         databufE%sizey = -1
         databufE%sizez = -1
         !--->
         databufE%buf_z_rx => null( )
         databufE%buf_z_tx => null( )
         databufE%buf_x_rx => null( )
         databufE%buf_x_tx => null( )
         databufE%buf_y_rx => null( )
         databufE%buf_y_tx => null( )
      endif
      !---------------- acaba InitMPI_Cray -----------------------------------------------------------
      return
   endsubroutine InitMPI_Cray
   !**************************************************************************************************
   subroutine FlushMPI_H_Cray
      !---------------- variables locales ------------------------------------------------------------
      type( t_databuf), pointer :: databuf_Up, databuf_Down
      integer :: ierr
      integer, dimension( 4) :: req1, req2
      integer, dimension( 2) :: req1b, req2b
      integer, dimension( MPI_STATUS_SIZE, 4) :: status1, status2
      integer, dimension( MPI_STATUS_SIZE, 2) :: status1b, status2b
      !---------------- empieza FlushMPI_H -----------------------------------------------------------
      databuf_Up => databuf_SetH%databuf_Up
      databuf_Down => databuf_SetH%databuf_Down
      !------------------------------------------->
      !print *,'flush H antes ',databuf_Down%ip_target,databuf_Up%ip_target
      !
      !syncUp AND PBC !only NEEDED BY THE PERIODIC BOUNDARY CONDITIONS !no real MPI burden since all layers communicate two sets
      if( databuf_SetH%syncUp .OR. databuf_SetH%pbcUp) then
         call MPI_VAMOS_ALLA_Hup( databuf_Up, req1, req1b)
      endif
      !------------------------------------------->
      !syncDown AND PBC !only NEEDED BY THE PERIODIC BOUNDARY CONDITIONS !no real MPI burden since all layers communicate two sets
      if( databuf_SetH%syncDown .OR. databuf_SetH%pbcDown) then
         call MPI_VAMOS_ALLA_Hdown( databuf_Down, req2, req2b)
      endif
      ! jag: yo compruebo que ha habido Rx de todos mis recepciones con MPI_TEST
      ! si asi continuo. no creo que haya mucha diferencia
      if( databuf_SetH%syncDown .OR. databuf_SetH%pbcDown) then
         call MPI_WAITALL( 4_4, req2, status2, ierr)
         if( databuf_Down%FlushExtraInfo) then
            call MPI_WAITALL( 2_4, req2b, status2b, ierr)
         endif
      endif
      if( databuf_SetH%syncUp .OR. databuf_SetH%pbcUp) then
         call MPI_WAITALL( 4_4, req1, status1, ierr)
         if( databuf_Up%FlushExtraInfo) then
            call MPI_WAITALL( 2_4, req1b, status1b, ierr)
         endif
      endif
      !   call MPI_Barrier(SUBCOMM_MPI,ierr)
      !-----------------------------------------------------------------------------------------------
      return
   contains
      subroutine MPI_VAMOS_ALLA_Hup( databufH, req, reqb)
         !---------------- inputs --------------------------------------------------------------------
         type( t_databuf), intent( IN) :: databufH
         !---------------- outputs -------------------------------------------------------------------
         integer, dimension( 4), intent( OUT) :: req
         integer, dimension( 2), intent( OUT) :: reqb
         integer :: ierr
         !---------------- empieza MPI_VAMOS_ALLA_Hup ------------------------------------------------
         call MPI_IRECV( databufH%buf_x_rx, databufH%sizex, REALSIZE, databufH%ip_target, 1_4, SUBCOMM_MPI, req( 1), ierr)
         call MPI_ISEND( databufH%buf_x_tx, databufH%sizex, REALSIZE, databufH%ip_target, 2_4, SUBCOMM_MPI, req( 2), ierr)
         call MPI_IRECV( databufH%buf_y_rx, databufH%sizey, REALSIZE, databufH%ip_target, 3_4, SUBCOMM_MPI, req( 3), ierr)
         call MPI_ISEND( databufH%buf_y_tx, databufH%sizey, REALSIZE, databufH%ip_target, 4_4, SUBCOMM_MPI, req( 4), ierr)
         if( databufH%FlushExtraInfo) then
            !         print *,'---fluHextraup>',layoutnumber
            call MPI_IRECV( databufH%buf_z_rx, databufH%sizez, REALSIZE, databufH%ip_target, 5_4, SUBCOMM_MPI, reqb( 1), ierr)
            call MPI_ISEND( databufH%buf_z_tx, databufH%sizez, REALSIZE, databufH%ip_target, 6_4, SUBCOMM_MPI, reqb( 2), ierr)
         endif
         !---------------- acaba MPI_VAMOS_ALLA_Hup --------------------------------------------------
         return
      endsubroutine MPI_VAMOS_ALLA_Hup
      !***********************************************************************************************
      subroutine MPI_VAMOS_ALLA_Hdown( databufH, req, reqb)
         !---------------- inputs --------------------------------------------------------------------
         type( t_databuf), intent( IN) :: databufH
         !---------------- outputs -------------------------------------------------------------------
         integer, dimension( 4), intent( OUT) :: req
         integer, dimension( 2), intent( OUT) :: reqb
         integer :: ierr
         !---------------- empieza MPI_VAMOS_ALLA_Hdown ----------------------------------------------
         call MPI_ISEND( databufH%buf_x_tx, databufH%sizex, REALSIZE, databufH%ip_target, 1_4, SUBCOMM_MPI, req( 1), ierr)
         call MPI_IRECV( databufH%buf_x_rx, databufH%sizex, REALSIZE, databufH%ip_target, 2_4, SUBCOMM_MPI, req( 2), ierr)
         call MPI_ISEND( databufH%buf_y_tx, databufH%sizey, REALSIZE, databufH%ip_target, 3_4, SUBCOMM_MPI, req( 3), ierr)
         call MPI_IRECV( databufH%buf_y_rx, databufH%sizey, REALSIZE, databufH%ip_target, 4_4, SUBCOMM_MPI, req( 4), ierr)
         if( databufH%FlushExtraInfo) then
            !         print *,'---fluHextraup>',layoutnumber
            call MPI_ISEND( databufH%buf_z_tx, databufH%sizez, REALSIZE, databufH%ip_target, 5_4, SUBCOMM_MPI, reqb( 1), ierr)
            call MPI_IRECV( databufH%buf_z_rx, databufH%sizez, REALSIZE, databufH%ip_target, 6_4, SUBCOMM_MPI, reqb( 2), ierr)
         endif
         !---------------- acaba MPI_VAMOS_ALLA_Hdown ------------------------------------------------
         return
      endsubroutine MPI_VAMOS_ALLA_Hdown
      !---------------- acaba FlushMPI_H -------------------------------------------------------------
   endsubroutine FlushMPI_H_Cray
   !**************************************************************************************************
   subroutine FlushMPI_E_Cray
      !---------------- variables locales ------------------------------------------------------------
      type( t_databuf), pointer :: databuf_Up, databuf_Down
      integer :: ierr
      integer, dimension( 2) :: req1, req2
      integer, dimension( 4) :: req1b, req2b
      integer, dimension( MPI_STATUS_SIZE, 2) :: status1, status2
      integer, dimension( MPI_STATUS_SIZE, 4) :: status1b, status2b
      !---------------- empieza FlushMPI_E -----------------------------------------------------------
      databuf_Up => databuf_SetE%databuf_Up
      databuf_Down => databuf_SetE%databuf_Down
      !------------------------------------------->
      !print *,'entrando en flush E antes up'
      if( databuf_SetE%syncUp .OR. databuf_SetE%pbcUp .or. databuf_up%FlushExtraInfo) then
         ! print *,'flush E antes up',databuf_Up%ip_target,databuf_Up%ip_target
         call MPI_VAMOS_ALLA_Eup( databuf_Up, req1, req1b)
      endif
      !------------------------------------------->
      if( databuf_SetE%syncDown .OR. databuf_SetE%pbcDown .or. databuf_down%FlushExtraInfo) then
         ! print *,'flush E antes down',databuf_Down%ip_target,databuf_Down%ip_target
         call MPI_VAMOS_ALLA_Edown( databuf_Down, req2, req2b)
      endif
      !------------------------------------------->
      if( databuf_SetE%syncDown .OR. databuf_SetE%pbcDown .or. databuf_down%FlushExtraInfo) then
         if( databuf_Down%FlushExtraInfo) then
            call MPI_WAITALL( 2_4, req2, status2, ierr)
            call MPI_WAITALL( 4_4, req2b, status2b, ierr)
         endif
      endif
      if( databuf_SetE%syncUp .OR. databuf_SetE%pbcUp .or. databuf_up%FlushExtraInfo) then
         if( databuf_Up%FlushExtraInfo) then
            call MPI_WAITALL( 2_4, req1, status1, ierr)
            call MPI_WAITALL( 4_4, req1b, status1b, ierr)
         endif
      endif
      !   call MPI_Barrier(SUBCOMM_MPI,ierr)
      !-----------------------------------------------------------------------------------------------
      return
   contains
      subroutine MPI_VAMOS_ALLA_Eup( databufE, req, reqb)
         !---------------- inputs --------------------------------------------------------------------
         type( t_databuf), intent( IN) :: databufE
         !---------------- outputs -------------------------------------------------------------------
         integer, dimension( 2), intent( OUT) :: req
         integer, dimension( 4), intent( OUT) :: reqb
         integer :: ierr
         !---------------- empieza MPI_VAMOS_ALLA_Eup ------------------------------------------------
         if( databufE%FlushExtraInfo) then
            !print *,'---fluEextraup>'
            call MPI_IRECV( databufE%buf_z_rx, databufE%sizez, REALSIZE, databufE%ip_target, 1_4, SUBCOMM_MPI, req( 1), ierr)
            call MPI_ISEND( databufE%buf_z_tx, databufE%sizez, REALSIZE, databufE%ip_target, 2_4, SUBCOMM_MPI, req( 2), ierr)
            !
            call MPI_IRECV( databufE%buf_x_rx, databufE%sizex, REALSIZE, databufE%ip_target, 3_4, SUBCOMM_MPI, reqb( 1), ierr)
            call MPI_ISEND( databufE%buf_x_tx, databufE%sizex, REALSIZE, databufE%ip_target, 4_4, SUBCOMM_MPI, reqb( 2), ierr)
            call MPI_IRECV( databufE%buf_y_rx, databufE%sizey, REALSIZE, databufE%ip_target, 5_4, SUBCOMM_MPI, reqb( 3), ierr)
            call MPI_ISEND( databufE%buf_y_tx, databufE%sizey, REALSIZE, databufE%ip_target, 6_4, SUBCOMM_MPI, reqb( 4), ierr)
         endif
         !---------------- acaba MPI_VAMOS_ALLA_Eup --------------------------------------------------
         return
      endsubroutine MPI_VAMOS_ALLA_Eup
      !***********************************************************************************************
      subroutine MPI_VAMOS_ALLA_Edown( databufE, req, reqb)
         !---------------- inputs --------------------------------------------------------------------
         type( t_databuf), intent( IN) :: databufE
         !---------------- outputs -------------------------------------------------------------------
         integer, dimension( 2), intent( OUT) :: req
         integer, dimension( 4), intent( OUT) :: reqb
         integer :: ierr
         !---------------- empieza MPI_VAMOS_ALLA_Edown ----------------------------------------------
         if( databufE%FlushExtraInfo) then
            !print *,'---fluEextradown>'
            call MPI_ISEND( databufE%buf_z_tx, databufE%sizez, REALSIZE, databufE%ip_target, 1_4, SUBCOMM_MPI, req( 1), ierr)
            call MPI_IRECV( databufE%buf_z_rx, databufE%sizez, REALSIZE, databufE%ip_target, 2_4, SUBCOMM_MPI, req( 2), ierr)
            !
            call MPI_ISEND( databufE%buf_x_tx, databufE%sizex, REALSIZE, databufE%ip_target, 3_4, SUBCOMM_MPI, reqb( 1), ierr)
            call MPI_IRECV( databufE%buf_x_rx, databufE%sizex, REALSIZE, databufE%ip_target, 4_4, SUBCOMM_MPI, reqb( 2), ierr)
            call MPI_ISEND( databufE%buf_y_tx, databufE%sizey, REALSIZE, databufE%ip_target, 5_4, SUBCOMM_MPI, reqb( 3), ierr)
            call MPI_IRECV( databufE%buf_y_rx, databufE%sizey, REALSIZE, databufE%ip_target, 6_4, SUBCOMM_MPI, reqb( 4), ierr)
         endif
         !---------------- acaba MPI_VAMOS_ALLA_Edown ------------------------------------------------
         return
      endsubroutine MPI_VAMOS_ALLA_Edown
      !---------------- acaba FlushMPI_E -------------------------------------------------------------
   endsubroutine FlushMPI_E_Cray
   !**************************************************************************************************
   !--->
   subroutine InitExtraFlushMPI_Cray( layoutnumber, sggsweep, sggalloc, med, nummed, sggmiez, sggMiHz, &
   Ex, Ey, Ez, Hx, Hy, Hz,therearemurborders)
      !---------------- inputs -----------------------------------------------------------------------
      logical, intent(in) :: therearemurborders
      type( XYZlimit_t), dimension( 1: 6), intent( IN) :: sggalloc, sggsweep
      integer( kind = 4), intent( IN) :: layoutnumber
      integer( kind = INTEGERSIZEOFMEDIAMATRICES), intent( IN) :: &
      sggMiEz( sggalloc( iEz)%XI: sggalloc( iEz)%XE, sggalloc( iEz)%YI: sggalloc( iEz)%YE, sggalloc( iEz)%ZI: sggalloc( iEz)%ZE), &
      sggMiHz( sggalloc( iHz)%XI: sggalloc( iHz)%XE, sggalloc( iHz)%YI: sggalloc( iHz)%YE, sggalloc( iHz)%ZI: sggalloc( iHz)%ZE)
      REAL (KIND=RKIND), intent( IN), target :: &
      Hx( sggalloc( iHx)%XI: sggalloc( iHx)%XE, sggalloc( iHx)%YI: sggalloc( iHx)%YE, sggalloc( iHx)%ZI: sggalloc( iHx)%ZE), &
      Hy( sggalloc( iHy)%XI: sggalloc( iHy)%XE, sggalloc( iHy)%YI: sggalloc( iHy)%YE, sggalloc( iHy)%ZI: sggalloc( iHy)%ZE), &
      Hz( sggalloc( iHz)%XI: sggalloc( iHz)%XE, sggalloc( iHz)%YI: sggalloc( iHz)%YE, sggalloc( iHz)%ZI: sggalloc( iHz)%ZE)
      REAL (KIND=RKIND), intent( IN), target :: &
      Ex( sggalloc( iEx)%XI: sggalloc( iEx)%XE, sggalloc( iEx)%YI: sggalloc( iEx)%YE, sggalloc( iEx)%ZI: sggalloc( iEx)%ZE), &
      Ey( sggalloc( iEy)%XI: sggalloc( iEy)%XE, sggalloc( iEy)%YI: sggalloc( iEy)%YE, sggalloc( iEy)%ZI: sggalloc( iEy)%ZE), &
      Ez( sggalloc( iEz)%XI: sggalloc( iEz)%XE, sggalloc( iEz)%YI: sggalloc( iEz)%YE, sggalloc( iEz)%ZI: sggalloc( iEz)%ZE)
      integer( kind = 4), intent( IN) :: nummed
      type( MediaData_t), dimension( 0: NumMed), intent( IN) :: med
      !---------------- variables locales ------------------------------------------------------------
      integer( kind = 4) :: j1, i1, jmed
      type( t_databuf), pointer :: databufH, databufE
      !---------------- empieza InitExtraFlushMPI_Cray ----------------------------------------------
      !thanks to overlapping of media matrix to detect media
      !an anistropic boundary will be detected both by the upper and the lower layouts
      !regardless of what happened with the wires (creo que esta formulacion no pasaria tampoco con
      !los wires, pero no lo he tocado)
      FlushExtraInfoDown = .false.
      FlushExtraInfoUp = .false.
      !
      if (therearemurborders) then
         FlushExtraInfoDown = .true.
         FlushExtraInfoUp = .true.
      endif
      !is enough to check Ez and Hz
      do jmed = 1, NumMed
         if( Med( jmed)%Is%Anisotropic) then
            !!!Ez
            do j1 = sggsweep( iEz)%YI, sggsweep( iEz)%YE
               do i1=sggsweep(iEz)%XI,sggsweep(iEz)%XE
                  if( (sggMiEz(i1, j1, comZ)) == jmed) then
                          FlushExtraInfoDown =.true.
                  endif
                  if( (sggMiEz(i1, j1, -1+comZ)) == jmed) then
                       FlushExtraInfoDown =.true.
                  endif
                  if( (sggMiEz(i1, j1, finZ)) == jmed)  then
                         FlushExtraInfoUp   =.true.
                  endif
                  if( (sggMiEz(i1, j1, 1+finZ)) == jmed) then
                        FlushExtraInfoUp   =.true.
                  endif
               enddo
            enddo
            !!!Hz
            do j1=sggsweep(iHz)%YI,sggsweep(iHz)%YE
               do i1=sggsweep(iHz)%XI,sggsweep(iHz)%XE
                  if ((sggMiHz(i1,j1,   comZ)) == jmed) then
                       FlushExtraInfoDown  = .true.
                  endif
                  if ((sggMiHz(i1,j1, 1+comZ)) == jmed) then
                       FlushExtraInfoDown  = .true.
                  endif
                  if ((sggMiHz(i1,j1, 1+finZ)) == jmed) then
                       FlushExtraInfoUp    = .true.
                  endif
                  if ((sggMiHz(i1,j1, 2+finZ)) == jmed) then
                       FlushExtraInfoUp    = .true.
                  endif
               enddo
            enddo
         endif
         if(Med( jmed)%Is%SGBC .or. Med( jmed)%Is%Multiport .or. Med( jmed)%Is%AnisMultiport) then
            !!!Hz
            do j1=sggsweep(iHz)%YI,sggsweep(iHz)%YE
               do i1=sggsweep(iHz)%XI,sggsweep(iHz)%XE
                  if ((sggMiHz(i1,j1,   comZ)) == jmed) then
                     FlushExtraInfoDown  = .true.
                  endif
                  if ((sggMiHz(i1,j1, 1+finZ)) == jmed) then
                       FlushExtraInfoUp    = .true.
                  endif
                  !creo que esto no es necesario para multiports de ss pero no creo que cargue mucho y no se si Ian lo necesita
                  !lo dejo por precaucion
                  if ((sggMiHz(i1,j1, 1+comZ)) == jmed) then
                       FlushExtraInfoDown  = .true.
                  endif
                  if ((sggMiHz(i1,j1, 2+finZ)) == jmed) then
                       FlushExtraInfoUp    = .true.
                  endif
               enddo
            enddo
         endif
      enddo
      !jag bug Antares mas de 65295 steps
      !print *,'------',FlushExtraInfoDown,FlushExtraInfoUp,comZ,finZ,sggMiHz(4,4,21)
      databufH => databuf_SetH%databuf_Up
      databufE => databuf_SetE%databuf_Up
      if(databuf_SetH%syncUp) then
         databufH%FlushExtraInfo = FlushExtraInfoUp
         databufE%FlushExtraInfo = FlushExtraInfoUp
         if(databufH%FlushExtraInfo) then
            databufE%sizex = sizeEx
            databufE%sizey = sizeEy
            databufE%sizez = sizeEz
            databufH%sizez = sizeHz
            databufH%buf_z_rx => Hz( HzXI: HzXE, HzYI: HzYE, finZ+2)
            databufH%buf_z_tx => Hz( HzXI: HzXE, HzYI: HzYE, finZ)
            databufE%buf_z_rx => Ez( EzXI: EzXE, EzYI: EzYE, finZ+1)
            databufE%buf_z_tx => Ez( EzXI: EzXE, EzYI: EzYE, finZ)
            databufE%buf_x_rx => Ex( ExXI: ExXE, ExYI: ExYE, finZ+2)
            databufE%buf_x_tx => Ex( ExXI: ExXE, ExYI: ExYE, finZ)
            databufE%buf_y_rx => Ey( EyXI: EyXE, EyYI: EyYE, finZ+2)
            databufE%buf_y_tx => Ey( EyXI: EyXE, EyYI: EyYE, finZ)
         endif
      endif
      !-----------------------------------------------------> DW
      databufH => databuf_SetH%databuf_Down
      databufE => databuf_SetE%databuf_Down
      if( databuf_SetH%syncDown) then
         databufH%FlushExtraInfo = FlushExtraInfoDown
         databufE%FlushExtraInfo = FlushExtraInfoDown
         if( databufH%FlushExtraInfo) then
            databufE%sizex = sizeEx
            databufE%sizey = sizeEy
            databufE%sizez = sizeEz
            databufH%sizez = sizeHz
            databufH%buf_z_tx => Hz( HzXI: HzXE, HzYI: HzYE, comZ+1)
            databufH%buf_z_rx => Hz( HzXI: HzXE, HzYI: HzYE, comZ-1)
            databufE%buf_z_tx => Ez( EzXI: EzXE, EzYI: EzYE, comZ)
            databufE%buf_z_rx => Ez( EzXI: EzXE, EzYI: EzYE, comZ-1)
            databufE%buf_x_tx => Ex( ExXI: ExXE, ExYI: ExYE, comZ+1)
            databufE%buf_x_rx => Ex( ExXI: ExXE, ExYI: ExYE, comZ-1)
            databufE%buf_y_tx => Ey( EyXI: EyXE, EyYI: EyYE, comZ+1)
            databufE%buf_y_rx => Ey( EyXI: EyXE, EyYI: EyYI, comZ-1)
         endif
      endif
      !---------------- acaba InitExtraFlushMPI_Cray ------------------------------------------------
      return
   endsubroutine InitExtraFlushMPI_Cray
  
    
#endif
   !del compilewithmpi

end module
 
MODULE build_t_linea_mpi

#ifdef CompileWithMPI
   !
   USE NFDETypes
   
contains


    subroutine build_derived_t_linea(mesg_mpi_t_linea)


    implicit none


    ! local
    integer(kind=4),parameter              :: number=2
    integer(kind=4)                        :: ierr, i
    integer(kind=4)                        :: block_lengths(1:number)
    integer(kind=MPI_ADDRESS_KIND) :: displacements(1:number)
    integer(kind=4)                        :: typelist(1:number)

    !output
    integer(kind=4),intent(out) :: mesg_mpi_t_linea

    !----------------------------------------

    !EL PRIMERO ES INTEGER
    typelist(1) = MPI_INTEGER4
    block_lengths(1) = 1
    displacements(1) = 0
    !EL SEGUNDO ES CHARACTER
    typelist(2) =  MPI_CHARACTER
    block_lengths(2) = BUFSIZE
    displacements(2) = 4 !el segundo se desplaza 4 porque el primero tiene 4 bytes

    ! build the derived data type
    call MPI_Type_create_struct(number,block_lengths,displacements,&
                    typelist,mesg_mpi_t_linea,ierr)
    if (ierr /= 0 ) then
        print *, 'got an error in type create: ', ierr
        call MPI_Abort(SUBCOMM_MPI, ierr, ierr)
    endif

    ! commit it to the system, so it knows we ll use it
    ! for communication
    call MPI_TYPE_COMMIT(mesg_mpi_t_linea,ierr)
    if (ierr /= 0 ) then
        print *, 'got an error in type commit: ', ierr
        call MPI_Abort(SUBCOMM_MPI, ierr, ierr)
    endif

    return

    end subroutine build_derived_t_linea

#endif
!------------- END SUBROUTINE----------------------------
end module build_t_linea_mpi