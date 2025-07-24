!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SEMBA_FDTD sOLVER MODULE
!  Creation date Date :  April, 8, 2010
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!__________________________________________________________________________________________________
!******************************** REVISAR PARA PGI (CRAY) *****************************************
!---> AdvanceMultiportE
!---> AdvanceAnisMultiportE
!---> AdvanceMultiportH
!---> AdvanceAnisMultiportH
!---> dfUpdateE
!---> dfUpdateH
!---> MinusCloneMagneticPMC
!________________________________________________________________________________________

module Solver_mod

   use fdetypes
   use report
   use PostProcessing
   use Ilumina
   use Observa
   use BORDERS_other
   use Borders_CPML
   use Borders_MUR
   use Resuming
   use nodalsources
   use Lumped
   use PMLbodies
   use xdmf
   use vtk
   use interpreta_switches_m, only: entrada_t
#ifdef CompileWithMPI
   use MPIcomm
#endif
#ifdef CompileWithStochastic
   use MPI_stochastic
#endif
#ifdef CompileWithNIBC
   use Multiports
#endif

#ifdef CompileWithStochastic
   use sgbc_stoch
#else
   use sgbc_NOstoch
#endif  
   use EDispersives
   use MDispersives
   use Anisotropic
   use HollandWires     

#ifdef CompileWithMTLN  
   use Wire_bundles_mtln_mod             
#endif       

#ifdef CompileWithBerengerWires
   use WiresBerenger
#ifdef CompileWithMPI
   use WiresBerenger_MPI
#endif
#endif

#ifdef CompileWithSlantedWires
   use WiresSlanted
   use estructura_slanted_m
#endif


#ifdef CompileWithConformal
   USE conformal_time_stepping_m
   USE CONFORMAL_MAPPED
#endif
   USE EpsMuTimeScale_m
   USE CALC_CONSTANTS
#ifdef CompileWithPrescale
   USE P_rescale
#endif              
#ifdef CompileWithMTLN
   ! use mtln_solver_mod, mtln_solver_t => mtln_t
   use mtln_types_mod, only: mtln_t
   use Wire_bundles_mtln_mod
#endif
!!
#ifdef CompileWithProfiling
   use nvtx
#endif
   implicit none


   type, public :: solver_t
      type(sim_control_t) :: control
      type(Logic_control) :: thereAre
      type(perform_t) :: perform, d_perform

      real(kind=rkind), pointer, dimension (:,:,:), contiguous :: Ex,Ey,Ez,Hx,Hy,Hz
      real(kind=rkind), pointer, dimension (:) :: Idxe, Idye, Idze, Idxh, Idyh, Idzh, dxe, dye, dze, dxh, dyh, dzh
      type(constants_t) :: g
      type(media_matrices_t) :: media
      real (kind=RKIND_tiempo) :: lastexecutedtime
      real (kind=RKIND) :: maxSourceValue

      integer (kind=4) :: initialtimestep, lastexecutedtimestep, ini_save, n_info, n

      type(bounds_t) :: bounds

      logical :: parar, everflushed = .false., still_planewave_time

#ifdef CompileWithMTLN
      type (mtln_t) :: mtln_parsed
#endif

   contains
      procedure :: init => solver_init
      procedure :: run => solver_run
      procedure :: end => solver_end
      procedure :: init_control => solver_init_control
      procedure, private :: init_fields
      procedure, private :: init_distances
      procedure :: launch_simulation
      procedure :: set_field_value
      procedure :: get_field_value
      procedure :: step
      procedure :: advanceE, advanceEx, advanceEy, advanceEz
      procedure :: advanceH, advanceHx, advanceHy, advanceHz
      procedure :: advancePlaneWaveE => solver_advancePlaneWaveE
      procedure :: advancePlaneWaveH => solver_advancePlaneWaveH
      procedure :: advanceWiresE => solver_advanceWiresE
      procedure :: advanceWiresH => solver_advanceWiresH
      procedure :: advancePMLE => solver_advancePMLE
      procedure :: advanceAnisotropicE => solver_advanceAnisotropicE
      procedure :: advanceAnisotropicH => solver_advanceAnisotropicH
      procedure :: advanceLumpedE => solver_advanceLumpedE
      procedure :: advanceNodalE => solver_advanceNodalE
      procedure :: advanceNodalH => solver_advanceNodalH
      procedure :: advancePMLbodyH => solver_advancePMLbodyH
      procedure :: advanceMagneticCPML => solver_advanceMagneticCPML
      procedure :: advanceSGBCE => solver_advanceSGBCE
      procedure :: advanceSGBCH => solver_advanceSGBCH
      procedure :: advanceEDispersiveE => solver_advanceEDispersiveE
      procedure :: advanceMDispersiveH => solver_advanceMDispersiveH
      procedure :: MinusCloneMagneticPMC => solver_MinusCloneMagneticPMC
      procedure :: CloneMagneticPeriodic => solver_CloneMagneticPeriodic
      procedure :: advanceMagneticMUR => solver_advanceMagneticMUR
      procedure :: destroy_and_deallocate
#ifdef CompileWithMTLN
      procedure :: launch_mtln_simulation
#endif
   end type

   contains

   subroutine solver_init_control(this, input, maxSourceValue, time_desdelanzamiento)
      class(solver_t) :: this
      type(entrada_t), intent(in) :: input
      real (kind=RKIND), intent(in) :: maxSourceValue
      REAL (kind=8), intent(in) :: time_desdelanzamiento


      this%control%maxSourceValue = maxSourceValue
      this%control%time_desdelanzamiento = time_desdelanzamiento

      this%control%simu_devia = input%simu_devia
      this%control%resume = input%resume
      this%control%saveall = input%saveall
      this%control%makeholes = input%makeholes
      this%control%connectendings = input%connectendings
      this%control%isolategroupgroups = input%isolategroupgroups
      this%control%createmap  = input%createmap
      this%control%groundwires = input%groundwires
      this%control%noSlantedcrecepelo = input%noSlantedcrecepelo
      this%control%SGBC  = input%SGBC
      this%control%SGBCDispersive = input%SGBCDispersive
      this%control%mibc = input%mibc
      this%control%ADE = input%ADE
      this%control%conformalskin  = input%conformalskin
      this%control%NOcompomur = input%NOcompomur
      this%control%strictOLD = input%strictOLD
      this%control%TAPARRABOS  = input%TAPARRABOS
      this%control%noconformalmapvtk = input%noconformalmapvtk
      this%control%hopf = input%hopf
      this%control%experimentalVideal = input%experimentalVideal
      this%control%forceresampled = input%forceresampled
      this%control%mur_second = input%mur_second
      this%control%MurAfterPML = input%MurAfterPML
      this%control%stableradholland = input%stableradholland
      this%control%singlefilewrite = input%singlefilewrite
      this%control%NF2FFDecim = input%NF2FFDecim
      this%control%sgbccrank = input%sgbccrank
      this%control%fieldtotl = input%fieldtotl
      this%control%permitscaling = input%permitscaling
      this%control%mtlnberenger = input%mtlnberenger
      this%control%niapapostprocess = input%niapapostprocess
      this%control%stochastic = input%stochastic
      this%control%verbose = input%verbose
      this%control%dontwritevtk = input%dontwritevtk
      this%control%use_mtln_wires = input%use_mtln_wires
      this%control%resume_fromold = input%resume_fromold
      this%control%vtkindex =  input%vtkindex
      this%control%createh5bin =  input%createh5bin
      this%control%wirecrank =  input%wirecrank
      this%control%fatalerror =  input%fatalerror

      this%control%cfl = input%cfl
      this%control%attfactorc = input%attfactorc
      this%control%attfactorw = input%attfactorw
      this%control%alphamaxpar = input%alphamaxpar
      this%control%alphaOrden = input%alphaOrden
      this%control%kappamaxpar = input%kappamaxpar
      this%control%mindistwires = input%mindistwires
      this%control%sgbcFreq = input%sgbcFreq
      this%control%sgbcresol = input%sgbcresol
      this%control%factorradius = input%factorradius
      this%control%factordelta = input%factordelta
      this%control%nEntradaRoot = trim(adjustl(input%nEntradaRoot))
      this%control%inductance_model = trim(adjustl(input%inductance_model))
      this%control%wiresflavor = trim(adjustl(input%wiresflavor))
      this%control%nresumeable2 = trim(adjustl(input%nresumeable2))
      this%control%opcionestotales = input%opcionestotales
      this%control%ficherohopf = input%ficherohopf
      this%control%finaltimestep = input%finaltimestep
      this%control%flushsecondsFields = input%flushsecondsFields
      this%control%flushsecondsData = input%flushsecondsData
      this%control%layoutnumber = input%layoutnumber
      this%control%mpidir = input%mpidir
      this%control%inductance_order = input%inductance_order
      this%control%wirethickness = input%wirethickness
      this%control%maxCPUtime = input%maxCPUtime
      this%control%SGBCDepth = input%SGBCDepth
      this%control%precision = input%precision
      this%control%size = input%size
      this%control%MEDIOEXTRA = input%MEDIOEXTRA
      this%control%facesNF2FF = input%facesNF2FF
      !this%control%EpsMuTimeScale_input_parameters = input%EpsMuTimeScale_input_parameters

#ifdef CompileWithConformal
      this%control%input_conformal_flag = input%input_conformal_flag
#endif

      call this%thereAre%reset()

   end subroutine


#ifdef CompileWithMTLN
   subroutine launch_mtln_simulation(this, mtln_parsed, nEntradaRoot, layoutnumber)
      class(solver_t) :: this
      type (mtln_t) :: mtln_parsed
      character (len=*), intent(in)  ::  nEntradaRoot
      integer (kind=4), intent(in) ::  layoutnumber

      call solveMTLNProblem(mtln_parsed)
      call reportSimulationEnd(layoutnumber)
      call FlushMTLNObservationFiles(nEntradaRoot, mtlnProblem = .true.)
   end subroutine
#endif

   subroutine init_fields(this, sgg)
      class(solver_t) :: this
      type (sggfdtdinfo), intent(in) :: sgg
      allocate ( &
      this%Ex(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE),&
      this%Ey(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE),&
      this%Ez(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE),&
      this%Hx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE),&
      this%Hy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE),&
      this%Hz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE))
      this%Ex = 0.0_RKIND; this%Ey = 0.0_RKIND; this%Ez = 0.0_RKIND; this%Hx = 0.0_RKIND; this%Hy = 0.0_RKIND; this%Hz = 0.0_RKIND
   end subroutine

   subroutine init_distances(this,sgg)
      class(solver_t) :: this
      type (sggfdtdinfo), intent(in) :: sgg
      integer :: i
      allocate ( & 
      this%dxe (sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      this%dye (sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      this%dze (sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      this%Idxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      this%Idye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      this%Idze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      this%dxh (sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      this%dyh (sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      this%dzh (sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE), &
      this%Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      this%Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      this%Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE))
      this%dxe=-1.0e10_RKIND
      this%dye=-1.0e10_RKIND
      this%dze=-1.0e10_RKIND
      this%dxh=-1.0e10_RKIND
      this%dyh=-1.0e10_RKIND
      this%dzh=-1.0e10_RKIND
      
      do i=sgg%ALLOC(iHx)%XI,sgg%ALLOC(iHx)%XE
         this%dxe(i)=sgg%DX(i)
      end do
      do i=sgg%ALLOC(iHy)%YI,sgg%ALLOC(iHy)%YE
         this%dye(i)=sgg%DY(i)
      end do
      do i=sgg%ALLOC(iHz)%ZI,sgg%ALLOC(iHz)%ZE
         this%dze(i)=sgg%DZ(i)
      end do
      do i=sgg%ALLOC(iEx)%XI,sgg%ALLOC(iEx)%XE
         this%dxh(i)=(sgg%DX(i)+sgg%DX(i-1))/2.0_RKIND
      end do
      do i=sgg%ALLOC(iEy)%YI,sgg%ALLOC(iEy)%YE
         this%dyh(i)=(sgg%DY(i)+sgg%DY(i-1))/2.0_RKIND
      end do
      do i=sgg%ALLOC(iEz)%ZI,sgg%ALLOC(iEz)%ZE
         this%dzh(i)=(sgg%DZ(i)+sgg%DZ(i-1))/2.0_RKIND
      end do

      this%Idxe=1.0_RKIND/this%dxe
      this%Idye=1.0_RKIND/this%dye
      this%Idze=1.0_RKIND/this%dze
      this%Idxh=1.0_RKIND/this%dxh
      this%Idyh=1.0_RKIND/this%dyh
      this%Idzh=1.0_RKIND/this%dzh
   end subroutine

   subroutine set_field_value(this, field_idx, i_range,j_range,k_range, field_value)
      class(solver_t) :: this
      integer (kind=4), intent(in) :: field_idx
      integer (kind=4), dimension(2), intent(in) :: i_range, j_range, k_range
      real (kind=rkind), intent(in) :: field_value
      
      real(kind=rkind), pointer, dimension (:,:,:) :: field
      integer(kind=4) :: i, j, k
      select case(field_idx)
      case(iEx)
         field => this%Ex
      case(iEy)
         field => this%Ey
      case(iEz)
         field => this%Ez
      case(iHx)
         field => this%Hx
      case(iHy)
         field => this%Hy
      case(iHz)
         field => this%Hz
      end select
      do i = i_range(1), i_range(2)
         do j = j_range(1), j_range(2)
            do k = k_range(1), k_range(2)
               field(i,j,k) = field_value
            end do
         end do
      end do
   end subroutine

   function get_field_value(this, field_idx, fi,fj,fk) result(res)
      class(solver_t) :: this
      integer (kind=4), intent(in) :: field_idx
      integer (kind=4), intent(in) :: fi, fj, fk
      real (kind=rkind) :: res
      
      real(kind=rkind), pointer, dimension (:,:,:) :: field
      select case(field_idx)
      case(iEx)
         field => this%Ex
      case(iEy)
         field => this%Ey
      case(iEz)
         field => this%Ez
      case(iHx)
         field => this%Hx
      case(iHy)
         field => this%Hy
      case(iHz)
         field => this%Hz
      end select
      res = field(fi,fj,fk)
   end function

   subroutine launch_simulation(this, sgg,media,tag_numbers,SINPML_Fullsize,fullsize,finishedwithsuccess,Eps0,Mu0,tagtype, &
                                input, maxSourceValue, time_desdelanzamiento)
      class(solver_t) :: this
      type (SGGFDTDINFO), intent(INOUT)   ::  sgg
      type(taglist_t), intent(in) :: tag_numbers
      type(media_matrices_t), intent(inout) :: media
      type (limit_t), dimension(1:6), intent(in) :: SINPML_fullsize,fullsize
      logical, intent(inout) :: finishedwithsuccess
      REAL (KIND=RKIND), intent(inout) :: eps0,mu0
      type (tagtype_t), intent(in) :: tagtype
      type(entrada_t), intent(in) :: input
      real (kind=RKIND), intent(in) :: maxSourceValue
      REAL (kind=8), intent(in) :: time_desdelanzamiento

      call this%init_control(input,maxSourceValue, time_desdelanzamiento)
      call this%init(sgg, eps0, mu0, media, SINPML_fullsize, fullsize, tag_numbers)
      call this%run(sgg, eps0, mu0, SINPML_fullsize, fullsize, tag_numbers, tagtype)
      call this%end(sgg, eps0, mu0, tagtype, finishedwithsuccess)

      
   end subroutine launch_simulation

   subroutine solver_init(this, sgg, eps0, mu0, media, sinPML_fullsize, fullsize, tag_numbers)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(inout) :: sgg
      real(kind=rkind), intent(inout) :: eps0,mu0
      type(media_matrices_t), intent(inout) :: media
      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize, fullsize
      type(taglist_t) :: tag_numbers

      integer(kind=4) :: i, j, k, field
      character (len=bufsize)  ::  whoami, chari, layoutcharID

      real(kind=rkind), pointer, dimension (:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=rkind), pointer, dimension (:) :: Idxe, Idye, Idze, Idxh, Idyh, Idzh, dxe, dye, dze, dxh, dyh, dzh

      real(kind=RKIND_tiempo) :: ultimodt
      
      character (len=bufsize) :: dubuf
      logical :: attinformado = .false.

! #ifdef compileWithMPI
      integer(kind=4) :: dummyMin,dummyMax, ierr
      real(kind=rkind) :: rdummy
! #endif
      this%media = media
      this%control%fatalerror=.false.

      this%parar=.false.
      call this%perform%reset()
      call this%d_perform%reset()
      call this%thereAre%reset()
      this%thereAre%MagneticMedia = sgg%thereareMagneticMedia
      this%thereAre%PMLMagneticMedia = sgg%therearePMLMagneticMedia

      !prechecking of no offsetting to prevent errors in case of modifications
      I=sgg%Alloc(iEx)%XI
      J=sgg%Alloc(iEx)%YI
      K=sgg%Alloc(iEx)%ZI
      do field=iEy,6
         if (sgg%Alloc(field)%XI /= I) call stoponerror(this%control%layoutnumber,this%control%size,'OFFSETS IN INITIAL COORD NOT ALLOWED')
         if (sgg%Alloc(field)%YI /= J) call stoponerror(this%control%layoutnumber,this%control%size,'OFFSETS IN INITIAL COORD NOT ALLOWED')
         if (sgg%Alloc(field)%ZI /= K) call stoponerror(this%control%layoutnumber,this%control%size,'OFFSETS IN INITIAL COORD NOT ALLOWED')
      END DO

      write(whoami,'(a,i5,a,i5,a)') '(',this%control%layoutnumber+1,'/',this%control%size,') '
      !file names
      write(chari,*) this%control%layoutnumber+1
      if ((this%control%layoutnumber == 0).and.this%control%verbose) call reportmedia(sgg)
      layoutcharID = trim(adjustl(this%control%nentradaroot))//'_'//trim(adjustl(chari))
      call findbounds(sgg,this%bounds)

      call this%init_distances(sgg)
      Idxe => this%Idxe; Idye => this%Idye; Idze => this%Idze; Idxh => this%Idxh; Idyh => this%Idyh; Idzh => this%Idzh; dxe => this%dxe; dye => this%dye; dze => this%dze; dxh => this%dxh; dyh => this%dyh; dzh => this%dzh
!!!lo cambio aqui permit scaling a 211118 por problemas con resuming: debe leer el eps0, mu0, antes de hacer numeros
      
      allocate (this%g%g1(0 : sgg%NumMedia),this%g%g2(0 : sgg%NumMedia),this%g%gm1(0 : sgg%NumMedia),this%g%gm2(0 : sgg%NumMedia))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Field matrices creation (an extra cell is padded at each limit and direction to deal with PMC imaging with no index errors)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !ojo las dimesniones deben ser giuales a las utlizadas en reallocate para las matrices sggmiEx, etc

      call this%init_fields(sgg)
      Ex => this%Ex; Ey => this%Ey; Ez => this%Ez; Hx => this%Hx; Hy => this%Hy; Hz => this%Hz
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Init the local variables and observation stuff needed by each module, taking into account resume status
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      dt0=sgg%dt !guardalo aqui para entrada pscale correcta si resume
      if (.not.this%control%resume) then
         Ex=0.0_RKIND; Ey=0.0_RKIND; Ez=0.0_RKIND; Hx=0.0_RKIND; Hy=0.0_RKIND; Hz=0.0_RKIND
         this%initialtimestep=0 
         this%lastexecutedtimestep=0
         this%lastexecutedtime=0.0_RKIND_tiempo
      else
         write(dubuf,*) 'Init processing resuming data'
         call print11(this%control%layoutnumber,dubuf)
         if (this%control%resume_fromold) then
            open (14,file=trim(adjustl(this%control%nresumeable2))//'.old',form='unformatted')
         else
            open (14,file=trim(adjustl(this%control%nresumeable2)),form='unformatted')
         endif
         call ReadFields(sgg%alloc,this%lastexecutedtimestep,this%lastexecutedtime,ultimodt,eps0,mu0,Ex,Ey,Ez,Hx,Hy,Hz)
         sgg%dt=ultimodt !para permit scaling
      !!!!!!!!!!!!No es preciso re-sincronizar pero lo hago !!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef CompileWithMPI
         rdummy=sgg%dt
         call MPIupdateMin(real(sgg%dt,RKIND),rdummy)
         rdummy=eps0
         call MPIupdateMin(eps0,rdummy)
         rdummy=mu0
         call MPIupdateMin(mu0,rdummy)
#endif
#ifdef CompileWithMPI
         call MPI_AllReduce( this%lastexecutedtimestep, dummyMin, 1_4, MPI_INTEGER, MPI_MIN, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( this%lastexecutedtimestep, dummyMax, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)
         if ((dummyMax /= this%lastexecutedtimestep).or.(dummyMin /= this%lastexecutedtimestep)) then
#ifdef CompileWithOldSaving
            if (this%control%resume_fromold) then
               close (14)
               write(dubuf,*) 'Incoherence between MPI saved steps for resuming.', dummyMin,dummyMax,this%lastexecutedtimesteP
               call stoponerror (this%control%layoutnumber,this%control%size,BUFF,.true.) !para que retorne
               call this%destroy_and_deallocate(sgg)
               return
            else
               write(dubuf,*) 'Incoherence between MPI saved steps for resuming. Retrying with -old....'
               call print11(this%control%layoutnumber,dubuf)
               this%control%resume_fromold=.true.
               close (14)
               open (14,file=trim(adjustl(this%control%nresumeable2))//'.old',form='unformatted')
               call ReadFields(sgg%alloc,this%lastexecutedtimestep,this%lastexecutedtime,ultimodt,eps0,mu0,Ex,Ey,Ez,Hx,Hy,Hz)
               sgg%dt=ultimodt !para permit scaling
               call MPI_AllReduce( this%lastexecutedtimestep, dummyMin, 1_4, MPI_INTEGER, MPI_MIN, SUBCOMM_MPI, ierr)
               call MPI_AllReduce( this%lastexecutedtimestep, dummyMax, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)
               if ((dummyMax /= this%lastexecutedtimestep).or.(dummyMin /= this%lastexecutedtimestep)) then
                  write(DUbuf,*) 'NO success. fields.old MPI are also incoherent for resuming.', dummyMin,dummyMax,this%lastexecutedtimestep
                  call stoponerror (this%control%layoutnumber,this%control%size,DUBUF,.true.) !para que retorne
                  call this%destroy_and_deallocate(sgg)
                  return
               else
                  write(dubuf,*) 'SUCCESS: Restarting from .fields.old instead. From n=',this%lastexecutedtimestep
                  call print11(this%control%layoutnumber,dubuf)
               endif
            endif
#else
            close (14)

            write(dubuf,*) 'Incoherence between MPI saved steps for resuming.',dummyMin,dummyMax,this%lastexecutedtimestep
            call stoponerror (this%control%layoutnumber,this%control%size,dubuf,.true.) !para que retorne
            call this%destroy_and_deallocate(sgg)
            return
#endif
         endif
#endif
         this%initialtimestep=this%lastexecutedtimestep+1
         write(dubuf,*) '[OK] processing resuming data. Last executed time step ',this%lastexecutedtimestep
         call print11(this%control%layoutnumber,dubuf)
      endif

      if (this%initialtimestep>this%control%finaltimestep) then
          call stoponerror (this%control%layoutnumber,this%control%size,'Initial time step greater than final one',.true.) !para que retorne
          call this%destroy_and_deallocate(sgg)
          return
      endif
!!!incializa el vector de tiempos para permit scaling 191118
      call crea_timevector(sgg,this%lastexecutedtimestep,this%control%finaltimestep,this%lastexecutedtime)
!!!!!!!!!!!!!!!!!!!!!

! !fin lo cambio aqui

      call updateSigmaM(attinformado)
      call updateThinWiresSigma(attinformado)
      call calc_G1G2Gm1Gm2(sgg,this%g,eps0,mu0)
      call revertThinWiresSigma()
 
      !
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) 'Init Reporting...';  call print11(this%control%layoutnumber,dubuf)
      call InitReporting(sgg,this%control)
      call reportSimulationOptions()

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
      !!!OJO SI SE CAMBIA EL ORDEN DE ESTAS INICIALIZACIONES HAY QUE CAMBIAR EL ORDEN DE STOREADO EN EL RESUMING
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      call initializeBorders()
      call initializeLumped()
      call initializeWires()
      call initializeAnisotropic()
      call initializeSGBC()
      call initializeMultiports()
      call initializeConformalElements()
      
      call initializeEDispersives()
      call initializeMDispersives()
      call initializePlanewave()
      call initializeNodalSources()

      call fillMtag(sgg, this%media%sggMiEx, this%media%sggMiEy, this%media%sggMiEz, this%media%sggMiHx, this%media%sggMiHy, this%media%sggMiHz,this%media%sggMtag, this%bounds, tag_numbers)
      call initializeObservation()

      !!!!voy a jugar con fuego !!!210815 sincronizo las matrices de medios porque a veces se precisan. Reutilizo rutinas viejas mias NO CRAY. Solo se usan aqui
      !MPI initialization
#ifdef CompileWithMPI
      call initializeMPI()
#endif
      
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      if (this%control%resume) close (14)
      !
      this%n=this%initialtimestep
      this%ini_save = this%initialtimestep
      this%n_info = 5 + this%initialtimestep

      write(dubuf,*) 'Init Timing...';  call print11(this%control%layoutnumber,dubuf)
      call InitTiming(sgg, this%control, this%control%time_desdelanzamiento, this%initialtimestep, this%control%maxSourceValue)


      CALL CLOSEWARNINGFILE(this%control%layoutnumber,this%control%size,this%control%fatalerror,.false.,this%control%simu_devia) !aqui ya esta dividido el stochastic y hay dos this%control%layoutnumber=0

      if (this%control%fatalerror) then
         dubuf='FATAL ERRORS. Revise *Warnings.txt file. ABORTING...'
         call stoponerror(this%control%layoutnumber,this%control%size,dubuf,.true.) !para que retorne
         call this%destroy_and_deallocate(sgg)
         return
      endif
#ifdef CompileWithMPI
      call flushMPIdata()
#endif

!!!no se si el orden wires - sgbcs del sync importa 150519
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
      if (this%control%stochastic)  then
         call syncstoch_mpi_sgbcs(this%control%simu_devia,this%control%layoutnumber,this%control%size)
         call syncstoch_mpi_lumped(this%control%simu_devia,this%control%layoutnumber,this%control%size)
      endif
#endif    
#endif    

      call printSimulationStart()
   
contains

      subroutine findbounds(sgg,b)
         !
         type (SGGFDTDINFO), intent(IN)     ::  sgg
         type (bounds_t), intent(out)  ::  b
         !

         !No tocar. Dejar como estan alocateados
         b%dxe%XI=sgg%alloc(iHx)%XI
         b%dxe%XE=sgg%alloc(iHx)%XE
         b%dye%YI=sgg%alloc(iHy)%YI
         b%dye%YE=sgg%alloc(iHy)%YE
         b%dze%ZI=sgg%alloc(iHz)%ZI
         b%dze%ZE=sgg%alloc(iHz)%ZE
         !
         b%dxh%XI=sgg%alloc(iEx)%XI
         b%dxh%XE=sgg%alloc(iEx)%XE
         b%dyh%YI=sgg%alloc(iEy)%YI
         b%dyh%YE=sgg%alloc(iEy)%YE
         b%dzh%ZI=sgg%alloc(iEz)%ZI
         b%dzh%ZE=sgg%alloc(iEz)%ZE

         !
         !No tocar. Dejar como estan alocateados
         b%Ex%XI=sgg%Alloc(iEx)%XI
         b%Ex%XE=sgg%Alloc(iEx)%XE
         b%Ey%XI=sgg%Alloc(iEy)%XI
         b%Ey%XE=sgg%Alloc(iEy)%XE
         b%Ez%XI=sgg%Alloc(iEz)%XI
         b%Ez%XE=sgg%Alloc(iEz)%XE
         !
         b%Hx%XI=sgg%Alloc(iHx)%XI
         b%Hx%XE=sgg%Alloc(iHx)%XE
         b%Hy%XI=sgg%Alloc(iHy)%XI
         b%Hy%XE=sgg%Alloc(iHy)%XE
         b%Hz%XI=sgg%Alloc(iHz)%XI
         b%Hz%XE=sgg%Alloc(iHz)%XE
         !
         b%Ex%YI=sgg%Alloc(iEx)%YI
         b%Ex%YE=sgg%Alloc(iEx)%YE
         b%Ey%YI=sgg%Alloc(iEy)%YI
         b%Ey%YE=sgg%Alloc(iEy)%YE
         b%Ez%YI=sgg%Alloc(iEz)%YI
         b%Ez%YE=sgg%Alloc(iEz)%YE
         !
         b%Hx%YI=sgg%Alloc(iHx)%YI
         b%Hx%YE=sgg%Alloc(iHx)%YE
         b%Hy%YI=sgg%Alloc(iHy)%YI
         b%Hy%YE=sgg%Alloc(iHy)%YE
         b%Hz%YI=sgg%Alloc(iHz)%YI
         b%Hz%YE=sgg%Alloc(iHz)%YE
         !
         b%Ex%ZI=sgg%Alloc(iEx)%ZI
         b%Ex%ZE=sgg%Alloc(iEx)%ZE
         b%Ey%ZI=sgg%Alloc(iEy)%ZI
         b%Ey%ZE=sgg%Alloc(iEy)%ZE
         b%Ez%ZI=sgg%Alloc(iEz)%ZI
         b%Ez%ZE=sgg%Alloc(iEz)%ZE
         !
         b%Hx%ZI=sgg%Alloc(iHx)%ZI
         b%Hx%ZE=sgg%Alloc(iHx)%ZE
         b%Hy%ZI=sgg%Alloc(iHy)%ZI
         b%Hy%ZE=sgg%Alloc(iHy)%ZE
         b%Hz%ZI=sgg%Alloc(iHz)%ZI
         b%Hz%ZE=sgg%Alloc(iHz)%ZE
         !
         !
         !

         !matrix indexes. Nothing to change. Asi estan alocateados
         b%sggMiEx%XI=sgg%Alloc(iEx)%XI
         b%sggMiEx%XE=sgg%Alloc(iEx)%XE
         b%sggMiEy%XI=sgg%Alloc(iEy)%XI
         b%sggMiEy%XE=sgg%Alloc(iEy)%XE
         b%sggMiEz%XI=sgg%Alloc(iEz)%XI
         b%sggMiEz%XE=sgg%Alloc(iEz)%XE
         !
         b%sggMiHx%XI=sgg%Alloc(iHx)%XI
         b%sggMiHx%XE=sgg%Alloc(iHx)%XE
         b%sggMiHy%XI=sgg%Alloc(iHy)%XI
         b%sggMiHy%XE=sgg%Alloc(iHy)%XE
         b%sggMiHz%XI=sgg%Alloc(iHz)%XI
         b%sggMiHz%XE=sgg%Alloc(iHz)%XE
         !
         b%sggMiEx%YI=sgg%Alloc(iEx)%YI
         b%sggMiEx%YE=sgg%Alloc(iEx)%YE
         b%sggMiEy%YI=sgg%Alloc(iEy)%YI
         b%sggMiEy%YE=sgg%Alloc(iEy)%YE
         b%sggMiEz%YI=sgg%Alloc(iEz)%YI
         b%sggMiEz%YE=sgg%Alloc(iEz)%YE
         !
         b%sggMiHx%YI=sgg%Alloc(iHx)%YI
         b%sggMiHx%YE=sgg%Alloc(iHx)%YE
         b%sggMiHy%YI=sgg%Alloc(iHy)%YI
         b%sggMiHy%YE=sgg%Alloc(iHy)%YE
         b%sggMiHz%YI=sgg%Alloc(iHz)%YI
         b%sggMiHz%YE=sgg%Alloc(iHz)%YE
         !
         b%sggMiEx%ZI=sgg%Alloc(iEx)%ZI
         b%sggMiEx%ZE=sgg%Alloc(iEx)%ZE
         b%sggMiEy%ZI=sgg%Alloc(iEy)%ZI
         b%sggMiEy%ZE=sgg%Alloc(iEy)%ZE
         b%sggMiEz%ZI=sgg%Alloc(iEz)%ZI
         b%sggMiEz%ZE=sgg%Alloc(iEz)%ZE
         !
         b%sggMiHx%ZI=sgg%Alloc(iHx)%ZI
         b%sggMiHx%ZE=sgg%Alloc(iHx)%ZE
         b%sggMiHy%ZI=sgg%Alloc(iHy)%ZI
         b%sggMiHy%ZE=sgg%Alloc(iHy)%ZE
         b%sggMiHz%ZI=sgg%Alloc(iHz)%ZI
         b%sggMiHz%ZE=sgg%Alloc(iHz)%ZE
         !
         !
         !
         b%sweepEx%XI=sgg%Sweep(iEx)%XI
         b%sweepEx%XE=sgg%Sweep(iEx)%XE
         b%sweepEy%XI=sgg%Sweep(iEy)%XI
         b%sweepEy%XE=sgg%Sweep(iEy)%XE
         b%sweepEz%XI=sgg%Sweep(iEz)%XI
         b%sweepEz%XE=sgg%Sweep(iEz)%XE
         !
         b%sweepHx%XI=sgg%Sweep(iHx)%XI
         b%sweepHx%XE=sgg%Sweep(iHx)%XE
         b%sweepHy%XI=sgg%Sweep(iHy)%XI
         b%sweepHy%XE=sgg%Sweep(iHy)%XE
         b%sweepHz%XI=sgg%Sweep(iHz)%XI
         b%sweepHz%XE=sgg%Sweep(iHz)%XE
         !
         !
         b%sweepEx%YI=sgg%Sweep(iEx)%YI
         b%sweepEx%YE=sgg%Sweep(iEx)%YE
         b%sweepEy%YI=sgg%Sweep(iEy)%YI
         b%sweepEy%YE=sgg%Sweep(iEy)%YE
         b%sweepEz%YI=sgg%Sweep(iEz)%YI
         b%sweepEz%YE=sgg%Sweep(iEz)%YE
         !
         b%sweepHx%YI=sgg%Sweep(iHx)%YI
         b%sweepHx%YE=sgg%Sweep(iHx)%YE
         b%sweepHy%YI=sgg%Sweep(iHy)%YI
         b%sweepHy%YE=sgg%Sweep(iHy)%YE
         b%sweepHz%YI=sgg%Sweep(iHz)%YI
         b%sweepHz%YE=sgg%Sweep(iHz)%YE
         !
         b%sweepEx%ZI=sgg%Sweep(iEx)%ZI
         b%sweepEx%ZE=sgg%Sweep(iEx)%ZE
         b%sweepEy%ZI=sgg%Sweep(iEy)%ZI
         b%sweepEy%ZE=sgg%Sweep(iEy)%ZE
         b%sweepEz%ZI=sgg%Sweep(iEz)%ZI
         b%sweepEz%ZE=sgg%Sweep(iEz)%ZE
         !
         b%sweepHx%ZI=sgg%Sweep(iHx)%ZI
         b%sweepHx%ZE=sgg%Sweep(iHx)%ZE
         b%sweepHy%ZI=sgg%Sweep(iHy)%ZI
         b%sweepHy%ZE=sgg%Sweep(iHy)%ZE
         b%sweepHz%ZI=sgg%Sweep(iHz)%ZI
         b%sweepHz%ZE=sgg%Sweep(iHz)%ZE
         !
         b%sweepSINPMLEx%XI=sgg%SINPMLSweep(iEx)%XI
         b%sweepSINPMLEy%XI=sgg%SINPMLSweep(iEy)%XI
         b%sweepSINPMLEz%XI=sgg%SINPMLSweep(iEz)%XI
         b%sweepSINPMLHx%XI=sgg%SINPMLSweep(iHx)%XI
         b%sweepSINPMLHy%XI=sgg%SINPMLSweep(iHy)%XI
         b%sweepSINPMLHz%XI=sgg%SINPMLSweep(iHz)%XI
         !
         b%sweepSINPMLEx%XE=sgg%SINPMLSweep(iEx)%XE
         b%sweepSINPMLEy%XE=sgg%SINPMLSweep(iEy)%XE
         b%sweepSINPMLEz%XE=sgg%SINPMLSweep(iEz)%XE
         b%sweepSINPMLHx%XE=sgg%SINPMLSweep(iHx)%XE
         b%sweepSINPMLHy%XE=sgg%SINPMLSweep(iHy)%XE
         b%sweepSINPMLHz%XE=sgg%SINPMLSweep(iHz)%XE
         !
         b%sweepSINPMLEx%YI=sgg%SINPMLSweep(iEx)%YI
         b%sweepSINPMLEy%YI=sgg%SINPMLSweep(iEy)%YI
         b%sweepSINPMLEz%YI=sgg%SINPMLSweep(iEz)%YI
         b%sweepSINPMLHx%YI=sgg%SINPMLSweep(iHx)%YI
         b%sweepSINPMLHy%YI=sgg%SINPMLSweep(iHy)%YI
         b%sweepSINPMLHz%YI=sgg%SINPMLSweep(iHz)%YI
         !
         b%sweepSINPMLEx%YE=sgg%SINPMLSweep(iEx)%YE
         b%sweepSINPMLEy%YE=sgg%SINPMLSweep(iEy)%YE
         b%sweepSINPMLEz%YE=sgg%SINPMLSweep(iEz)%YE
         b%sweepSINPMLHx%YE=sgg%SINPMLSweep(iHx)%YE
         b%sweepSINPMLHy%YE=sgg%SINPMLSweep(iHy)%YE
         b%sweepSINPMLHz%YE=sgg%SINPMLSweep(iHz)%YE
         !
         b%sweepSINPMLEx%ZI=sgg%SINPMLSweep(iEx)%ZI
         b%sweepSINPMLEy%ZI=sgg%SINPMLSweep(iEy)%ZI
         b%sweepSINPMLEz%ZI=sgg%SINPMLSweep(iEz)%ZI
         b%sweepSINPMLHx%ZI=sgg%SINPMLSweep(iHx)%ZI
         b%sweepSINPMLHy%ZI=sgg%SINPMLSweep(iHy)%ZI
         b%sweepSINPMLHz%ZI=sgg%SINPMLSweep(iHz)%ZI
         !
         b%sweepSINPMLEx%ZE=sgg%SINPMLSweep(iEx)%ZE
         b%sweepSINPMLEy%ZE=sgg%SINPMLSweep(iEy)%ZE
         b%sweepSINPMLEz%ZE=sgg%SINPMLSweep(iEz)%ZE
         b%sweepSINPMLHx%ZE=sgg%SINPMLSweep(iHx)%ZE
         b%sweepSINPMLHy%ZE=sgg%SINPMLSweep(iHy)%ZE
         b%sweepSINPMLHz%ZE=sgg%SINPMLSweep(iHz)%ZE

         !

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !find lenghts
         !this is automatic. Nothing to change
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         b%Ex%NX=b%Ex%XE-b%Ex%XI+1
         b%Ex%NY=b%Ex%YE-b%Ex%YI+1
         b%Ex%NZ=b%Ex%ZE-b%Ex%ZI+1

         b%Ey%NX=b%Ey%XE-b%Ey%XI+1
         b%Ey%NY=b%Ey%YE-b%Ey%YI+1
         b%Ey%NZ=b%Ey%ZE-b%Ey%ZI+1

         b%Ez%NX=b%Ez%XE-b%Ez%XI+1
         b%Ez%NY=b%Ez%YE-b%Ez%YI+1
         b%Ez%NZ=b%Ez%ZE-b%Ez%ZI+1
         !
         b%Hx%NX=b%Hx%XE-b%Hx%XI+1
         b%Hx%NY=b%Hx%YE-b%Hx%YI+1
         b%Hx%NZ=b%Hx%ZE-b%Hx%ZI+1
         !
         b%Hy%NX=b%Hy%XE-b%Hy%XI+1
         b%Hy%NY=b%Hy%YE-b%Hy%YI+1
         b%Hy%NZ=b%Hy%ZE-b%Hy%ZI+1
         !
         b%Hz%NX=b%Hz%XE-b%Hz%XI+1
         b%Hz%NY=b%Hz%YE-b%Hz%YI+1
         b%Hz%NZ=b%Hz%ZE-b%Hz%ZI+1
         !
         !
         b%sweepEx%NX=b%sweepEx%XE-b%sweepEx%XI+1
         b%sweepEx%NY=b%sweepEx%YE-b%sweepEx%YI+1
         b%sweepEx%NZ=b%sweepEx%ZE-b%sweepEx%ZI+1
         !
         b%sweepEy%NX=b%sweepEy%XE-b%sweepEy%XI+1
         b%sweepEy%NY=b%sweepEy%YE-b%sweepEy%YI+1
         b%sweepEy%NZ=b%sweepEy%ZE-b%sweepEy%ZI+1
         !
         b%sweepEz%NX=b%sweepEz%XE-b%sweepEz%XI+1
         b%sweepEz%NY=b%sweepEz%YE-b%sweepEz%YI+1
         b%sweepEz%NZ=b%sweepEz%ZE-b%sweepEz%ZI+1
         !
         b%sweepHx%NX=b%sweepHx%XE-b%sweepHx%XI+1
         b%sweepHx%NY=b%sweepHx%YE-b%sweepHx%YI+1
         b%sweepHx%NZ=b%sweepHx%ZE-b%sweepHx%ZI+1
         !
         b%sweepHy%NX=b%sweepHy%XE-b%sweepHy%XI+1
         b%sweepHy%NY=b%sweepHy%YE-b%sweepHy%YI+1
         b%sweepHy%NZ=b%sweepHy%ZE-b%sweepHy%ZI+1
         !
         b%sweepHz%NX=b%sweepHz%XE-b%sweepHz%XI+1
         b%sweepHz%NY=b%sweepHz%YE-b%sweepHz%YI+1
         b%sweepHz%NZ=b%sweepHz%ZE-b%sweepHz%ZI+1
         !
         !
         b%sggMiEx%NX=b%sggMiEx%XE-b%sggMiEx%XI+1
         b%sggMiEx%NY=b%sggMiEx%YE-b%sggMiEx%YI+1
         b%sggMiEx%NZ=b%sggMiEx%ZE-b%sggMiEx%ZI+1
         b%sggMiEy%NX=b%sggMiEy%XE-b%sggMiEy%XI+1
         b%sggMiEy%NY=b%sggMiEy%YE-b%sggMiEy%YI+1
         b%sggMiEy%NZ=b%sggMiEy%ZE-b%sggMiEy%ZI+1
         b%sggMiEz%NX=b%sggMiEz%XE-b%sggMiEz%XI+1
         b%sggMiEz%NY=b%sggMiEz%YE-b%sggMiEz%YI+1
         b%sggMiEz%NZ=b%sggMiEz%ZE-b%sggMiEz%ZI+1
         !
         b%sggMiHx%NX=b%sggMiHx%XE-b%sggMiHx%XI+1
         b%sggMiHx%NY=b%sggMiHx%YE-b%sggMiHx%YI+1
         b%sggMiHx%NZ=b%sggMiHx%ZE-b%sggMiHx%ZI+1
         b%sggMiHy%NX=b%sggMiHy%XE-b%sggMiHy%XI+1
         b%sggMiHy%NY=b%sggMiHy%YE-b%sggMiHy%YI+1
         b%sggMiHy%NZ=b%sggMiHy%ZE-b%sggMiHy%ZI+1
         b%sggMiHz%NX=b%sggMiHz%XE-b%sggMiHz%XI+1
         b%sggMiHz%NY=b%sggMiHz%YE-b%sggMiHz%YI+1
         b%sggMiHz%NZ=b%sggMiHz%ZE-b%sggMiHz%ZI+1
         !
         !
         !estas longitudes son relativas al layout !ojo
         b%dxe%NX=b%dxe%XE-b%dxe%XI+1
         b%dye%NY=b%dye%YE-b%dye%YI+1
         b%dze%NZ=b%dze%ZE-b%dze%ZI+1
         !
         b%dxh%NX=b%dxh%XE-b%dxh%XI+1
         b%dyh%NY=b%dyh%YE-b%dyh%YI+1
         b%dzh%NZ=b%dzh%ZE-b%dzh%ZI+1


      end subroutine

      subroutine updateSigmaM(att)
         logical, intent(inout) :: att
         real(kind=rkind) :: deltaespmax, fmax, skin_depth
         logical :: hayattmedia = .false.
         REAL (kind = rkind) :: mur,epr
         character(len=BUFSIZE) :: buff   
         integer :: i
         if (abs(this%control%attfactorc-1.0_RKIND) > 1.0e-12_RKIND) then
            att=.false.
            do i=1,sgg%nummedia
               if (sgg%Med(i)%Is%MultiportPadding) then
                  sgg%Med(i)%SigmaM =(-2.0_RKIND * (-1.0_RKIND + this%control%attfactorc)*mu0)/((1 + this%control%attfactorc)*sgg%dt)
                  hayattmedia=.true.
               endif
               deltaespmax=max(max(maxval(sgg%dx),maxval(sgg%dy)),maxval(sgg%dz))
               if (hayattmedia.and. .not. att) then
                  !!!!info on stabilization
                  epr   =1.0_RKIND
                  mur   =1.0_RKIND
                  !!
                  write(buff,'(a,2e10.2e3)') ' Composites stabilization att. factor=',this%control%attfactorc,sgg%Med(i)%SigmaM

                  call WarnErrReport(buff)
                  !!
                  fmax=1.0_RKIND / (10.0_RKIND * sgg%dt)
                  skin_depth=1.0_RKIND / (Sqrt(2.0_RKIND)*fmax*Pi*(epr*Eps0**2*(4*mur*mu0**2.0_RKIND + sgg%Med(i)%Sigmam**2/(fmax**2*Pi**2.0_RKIND )))**0.25_RKIND * &
                  Sin(atan2(2*Pi*epr*Eps0*mur*mu0, - (epr*eps0*sgg%Med(i)%Sigmam)/fmax)/2.0_RKIND))
                  write(buff,'(a,e9.2e2,a,e10.2e3)') ' At 10 samp/per f=',fmax,',Max Att(dB)=', &
                  -(0.0001295712360834271997*AIMAG(fmax*Sqrt((epr*((0,-2.825225e7) + &
                  8.8757061047382236e6*mur + this%control%attfactorc*((0,2.825225e7) + 8.8757061047382236e6*mur)))/ &
                  (1.124121310242e12 + 1.124121310242e12*this%control%attfactorc))*min(deltaespmax,skin_depth)))
                  if (this%control%layoutnumber == 0) call WarnErrReport(buff)
                  if (fmax > 3e9) then
                     fmax=3e9
                     write(buff,'(a,e9.2e2,a,e10.2e3)') '             At f=',fmax,',Max Att(dB)=', &
                     -(0.0001295712360834271997*AIMAG(fmax*Sqrt((epr*((0,-2.825225e7) + &
                     8.8757061047382236e6*mur + this%control%attfactorc*((0,2.825225e7) + 8.8757061047382236e6*mur)))/ &
                     (1.124121310242e12 + 1.124121310242e12*this%control%attfactorc))*min(deltaespmax,skin_depth)))
                     if (this%control%layoutnumber == 0) call WarnErrReport(buff)
                  endif
                  att=.true.
               endif
            end do
         endif
      end subroutine updateSigmaM

      subroutine updateThinWiresSigma(att)
         logical, intent(inout) :: att
         character(len=BUFSIZE) :: buff   
         integer :: i
         if (abs(this%control%attfactorw-1.0_RKIND) > 1.0e-12_RKIND) then
            att=.false.
            do i=1,sgg%nummedia
               if (sgg%Med(i)%Is%ThinWire) then
                  sgg%Med(i)%Sigma =(-2.0_RKIND * (-1.0_RKIND + this%control%attfactorw)*eps0)/((1 + this%control%attfactorw)*sgg%dt)
                  if (.not.att) then
                     write(buff,'(a,2e10.2e3)') ' WIREs stabilization att. factors=',this%control%attfactorw,sgg%Med(i)%Sigma
                     if (this%control%layoutnumber == 0) call WarnErrReport(buff)
                     att=.true.
                  endif
               endif
            end do
         endif
      end subroutine updateThinWiresSigma

      subroutine revertThinWiresSigma()
         integer :: i
         if (abs(this%control%attfactorw-1.0_RKIND) > 1.0e-12_RKIND) then
            do i=1,sgg%nummedia
               if (sgg%Med(i)%Is%ThinWire) then
                  sgg%Med(i)%Sigma = 0.0_RKIND !revert!!! !necesario para no lo tome como un lossy luego en wires !solo se toca el g1,g2
               endif
            end do
         endif
      end subroutine

      subroutine reportSimulationOptions()
         character(len=BUFSIZE) :: buff   
         if ((this%control%layoutnumber == 0).and.this%control%verbose) then
            write(buff,'(a,3e9.2e2)') 'CPML  alpha, alphaorder, kappa factors= ', this%control%alphamaxpar,this%control%alphaOrden,this%control%kappamaxpar
            call WarnErrReport(buff)
            if (this%control%medioextra%exists) then
               write(buff,'(a,i5,e9.2e2)') 'CPML correction size,factor to scale sigmamax = ', &
               this%control%medioextra%size,this%control%medioextra%sigma
               call WarnErrReport(buff)
            endif
            write(buff,*) 'saveall=',this%control%saveall,', flushsecondsFields=',this%control%flushsecondsFields,', flushsecondsData=',this%control%flushsecondsData,', maxCPUtime=',this%control%maxCPUtime,', singlefilewrite=',this%control%singlefilewrite
            call WarnErrReport(buff)
            write(buff,*) 'TAPARRABOS=',this%control%TAPARRABOS,', wiresflavor=',trim(adjustl(this%control%wiresflavor)),', mindistwires=',this%control%mindistwires,', wirecrank=',this%control%wirecrank , 'makeholes=',this%control%makeholes
            call WarnErrReport(buff)
            write(buff,*) 'use_mtln_wires=', this%control%use_mtln_wires
            write(buff,*) 'connectendings=',this%control%connectendings,', isolategroupgroups=',this%control%isolategroupgroups
            call WarnErrReport(buff)
            write(buff,*) 'wirethickness ', this%control%wirethickness, 'stableradholland=',this%control%stableradholland,'mtlnberenger=',this%control%mtlnberenger,' inductance_model=',trim(adjustl(this%control%inductance_model)), &
                        ', inductance_order=',this%control%inductance_order,', groundwires=',this%control%groundwires,' ,fieldtotl=',this%control%fieldtotl,' noSlantedcrecepelo =',this%control%noSlantedcrecepelo 
            call WarnErrReport(buff)
            write(buff,*) 'sgbc=',this%control%sgbc,', mibc=',this%control%mibc,', attfactorc=',this%control%attfactorc,', attfactorw=',this%control%attfactorw
            call WarnErrReport(buff)
            write(buff,*) 'NOcompomur=',this%control%NOcompomur,', ADE=',this%control%ADE,', conformalskin=',this%control%conformalskin,', sgbcFreq=',this%control%sgbcFreq,', sgbcresol=',this%control%sgbcresol,', sgbccrank=',this%control%sgbccrank,', sgbcDepth=',this%control%sgbcdepth
            call WarnErrReport(buff)
            write(buff,*) 'mur_second=',this%control%mur_second,', murafterpml=',this%control%murafterpml,', facesNF2FF%tr=',this%control%facesNF2FF%tr,', facesNF2FF%fr=',this%control%facesNF2FF%fr,', facesNF2FF%iz=',this%control%facesNF2FF%iz
            call WarnErrReport(buff)
            write(buff,*) 'facesNF2FF%de=',this%control%facesNF2FF%de,', facesNF2FF%ab=',this%control%facesNF2FF%ab,', facesNF2FF%ar=',this%control%facesNF2FF%ar,', NF2FFDecim=',this%control%NF2FFDecim
            call WarnErrReport(buff)
         endif
      end subroutine

      subroutine initializeBorders()
         character(len=BUFSIZE) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer (kind=4) :: ierr
#endif
         write(dubuf,*) 'Init Other Borders...';  call print11(this%control%layoutnumber,dubuf)
         call InitOtherBorders    (sgg,this%thereAre)
         l_auxinput=this%thereAre%PECBorders.or.this%thereAre%PMCBorders.or.this%thereAre%PeriodicBorders
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if ( l_auxoutput) then
               write (dubuf,*) '----> there are PEC, PMC or periodic Borders';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no PEC, PMC or periodic Borders found';  call print11(this%control%layoutnumber,dubuf)
            endif
            
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init CPML Borders...';  call print11(this%control%layoutnumber,dubuf)
         call InitCPMLBorders     (sgg,SINPML_Fullsize,this%thereAre%PMLBorders,this%control, &
                                 dxe,dye,dze,dxh,dyh,dzh,Idxe,Idye,Idze,Idxh,Idyh,Idzh,eps0,mu0)

         l_auxinput=this%thereAre%PMLBorders
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if (l_auxoutput ) then
               write (dubuf,*) '----> there are CPML Borders';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no CPML Borders found';  call print11(this%control%layoutnumber,dubuf)
         endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init PML Bodies...';  call print11(this%control%layoutnumber,dubuf)
         call InitPMLbodies(sgg,this%media,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,this%g%g2,this%g%gm2,this%thereAre%PMLbodies,this%control,eps0,mu0)
         l_auxinput=this%thereAre%PMLbodies
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if ( l_auxoutput) then
               write (dubuf,*) '----> there are PML Bodies';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no PML Bodies found';  call print11(this%control%layoutnumber,dubuf)
         endif
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Mur Borders...';  call print11(this%control%layoutnumber,dubuf)
         call InitMURBorders      (sgg,this%thereAre%MURBorders,this%control%resume,Idxh,Idyh,Idzh,eps0,mu0)
         l_auxinput= this%thereAre%MURBorders
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         if (l_auxoutput) then
               write (dubuf,*) '----> there are Mur Borders';  call print11(this%control%layoutnumber,dubuf)
         else
               write(dubuf,*) '----> no Mur Borders found';  call print11(this%control%layoutnumber,dubuf)
         endif

      end subroutine initializeBorders

      subroutine initializeLumped()
         character(len=BUFSIZE) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

         !init lumped debe ir antes de wires porque toca la conductividad del material !mmmm ojoooo 120123
         write(dubuf,*) 'Init Lumped Elements...';  call print11(this%control%layoutnumber,dubuf)
         CALL InitLumped(sgg,this%media,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,this%control,this%thereAre%Lumpeds,eps0,mu0)
         l_auxinput=this%thereAre%Lumpeds
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif   
         if (l_auxoutput ) then
             write (dubuf,*) '----> there are Structured lumped elements';  call print11(this%control%layoutnumber,dubuf)
         else
              write(dubuf,*) '----> no lumped Structured elements found';  call print11(this%control%layoutnumber,dubuf)
         endif
      end subroutine initializeLumped

      subroutine initializeWires()
         real (kind=rkind) :: dtcritico, newdtcritico
         character(len=BUFSIZE) :: dubuf, buff
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

         dtcritico=sgg%dt
         if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
            (trim(adjustl(this%control%wiresflavor))=='transition')) then
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init Holland Wires...';  call print11(this%control%layoutnumber,dubuf)
            call InitWires       (sgg,this%media%sggMiNo,this%media%sggMiEx,this%media%sggMiEy,this%media%sggMiEz,this%media%sggMiHx,this%media%sggMiHy,this%media%sggMiHz, & 
                                 this%thereAre%Wires, Ex,Ey,Ez,Hx,Hy,Hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh, &
                                 this%g%g2,SINPML_fullsize, fullsize,dtcritico,eps0,mu0,this%control)
            l_auxinput=this%thereAre%Wires
            l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if (l_auxoutput ) then
               write (dubuf,*) '----> there are Holland/transition wires';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no Holland/transition wires found';  call print11(this%control%layoutnumber,dubuf)
         endif
         endif

#ifdef CompileWithBerengerWires
         if (trim(adjustl(this%control%wiresflavor))=='berenger') then

#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init Multi-Wires...';  call print11(this%control%layoutnumber,dubuf)
            call InitWires_Berenger(sgg,this%media%sggMiNo,this%media%sggMiEx,this%media%sggMiEy,this%media%sggMiEz,this%media%sggMiHx,this%media%sggMiHy,this%media%sggMiHz,this%control%layoutnumber,this%control%size,this%thereAre%Wires,this%control%resume,this%control%makeholes, &
            this%control%isolategroupgroups,this%control%mtlnberenger,this%control%mindistwires, &
            this%control%groundwires,this%control%taparrabos,Ex,Ey,Ez, &
            Idxe,Idye,Idze,Idxh,Idyh,Idzh,this%control%inductance_model,this%g%g2,SINPML_fullsize,fullsize,dtcritico,eps0,mu0,this%control%verbose)
         l_auxinput= this%thereAre%Wires
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         
            if (l_auxoutput) then
               write (dubuf,*) '----> there are Multi-wires';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no Multi-wires found';  call print11(this%control%layoutnumber,dubuf)
            endif
         endif
#endif
#ifdef CompileWithSlantedWires
         if((trim(adjustl(this%control%wiresflavor))=='slanted').or.(trim(adjustl(this%control%wiresflavor))=='semistructured')) then

#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init Slanted Wires...';  call print11(this%control%layoutnumber,dubuf)
            if ((trim(adjustl(this%control%wiresflavor))=='semistructured')) then
               write(dubuf,*) '...',this%control%precision;  call print11(this%control%layoutnumber,dubuf)
               call estructura_slanted(sgg,this%control%precision)
            else
               continue
            endif
            call InitWires_Slanted(sgg, this%control%layoutnumber,this%control%size, Ex, Ey, Ez,   & 
                                    Idxe, Idye, Idze, Idxh, Idyh, Idzh,   &
                                    this%media%sggMiNo,                              &
                                    this%media%sggMiEx, this%media%sggMiEy, this%media%sggMiEz,            &
                                    this%media%sggMiHx, this%media%sggMiHy, this%media%sggMiHz,            &
                                    this%thereAre%Wires, this%control%resume,               &
                                    this%control%mindistwires, this%control%groundwires,this%control%noSlantedcrecepelo ,     &
                                    this%control%inductance_model, this%control%inductance_order,   &
                                    this%g%g2, SINPML_fullsize, dtcritico,eps0,mu0,this%control%verbose)
            l_auxinput=this%thereAre%Wires
            l_auxoutput=l_auxinput
!check for MUR1 nodes sgg 230124
            call init_murABC_slanted(sgg,SINPML_Fullsize,eps0,mu0)
!!!!!!         
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         
            if (l_auxoutput ) then
               write (dubuf,*) '----> there are Slanted wires';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no Slanted wires found';  call print11(this%control%layoutnumber,dubuf)
            endif
         endif
#endif
      !!!sincroniza el dtcritico
#ifdef CompileWithMPI
         call MPI_AllReduce( dtcritico, newdtcritico, 1_4, REALSIZE, MPI_MIN, SUBCOMM_MPI, ierr)
         dtcritico=newdtcritico
#endif
         if (sgg%dt <= dtcritico) then
            write(buff,'(a,e10.2e3)')  'WIR_INFO: deltat for stability OK: ',dtcritico
            if ((this%control%layoutnumber==0).and.this%control%verbose) call WarnErrReport(buff)
         else
            if (.not.(this%control%resume.and.this%control%permitscaling)) then !no abortasr solo advertir si permittivity scaling
               write(buff,'(a,e10.2e3)')  'WIR_ERROR: Possibly UNSTABLE dt, decrease wire radius, number of parallel WIREs, use -stableradholland or make dt < ',dtcritico
               if (this%control%layoutnumber==0) call WarnErrReport(buff,.true.)
            else
               write(buff,'(a,e10.2e3)')  'WIR_WARNING: Resume and Pscaling with wires. Possibly UNSTABLE dt, decrease wire radius, number of parallel WIREs: dt is over ',dtcritico
               if (this%control%layoutnumber==0) call WarnErrReport(buff,.false.)
            endif
         endif
      !!!
!!
         if (this%control%use_mtln_wires) then
#ifdef CompileWithMTLN
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init MTLN Wires...';  call print11(this%control%layoutnumber,dubuf)
            call InitWires_mtln(sgg,Ex,Ey,Ez,Idxh,Idyh,Idzh,eps0, mu0, this%mtln_parsed,this%thereAre%MTLNbundles)
#else
            write(buff,'(a)') 'WIR_ERROR: Executable was not compiled with MTLN modules.'
#endif
         endif

      end subroutine initializeWires

      subroutine initializeAnisotropic()
         character(len=BUFSIZE) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
         integer :: rank
#endif
      
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Anisotropic...';  call print11(this%control%layoutnumber,dubuf)
         call InitAnisotropic(sgg,this%media,this%thereAre%Anisotropic,this%thereAre%ThinSlot,eps0,mu0)
         l_auxinput=this%thereAre%Anisotropic.or.this%thereAre%ThinSlot
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_COMM_RANK(SUBCOMM_MPI, rank, ierr)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif   
         if (l_auxoutput) then
               write (dubuf,*) '----> there are Structured anisotropic elements';  call print11(this%control%layoutnumber,dubuf)
         else
               write(dubuf,*) '----> no Structured anisotropic elements found';  call print11(this%control%layoutnumber,dubuf)
         endif
      end subroutine initializeAnisotropic

      subroutine initializeSGBC()
         character(len=BUFSIZE) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

         IF (this%control%sgbc)  then
#ifdef CompileWithMPI
              call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
               write(dubuf,*) 'Init Multi sgbc...';  call print11(this%control%layoutnumber,dubuf)
               call Initsgbcs(sgg,this%media,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,this%control%layoutnumber,this%control%size, &
                  this%g,this%thereAre%sgbcs,this%control%resume,this%control%sgbccrank,this%control%sgbcFreq,this%control%sgbcresol,this%control%sgbcdepth,this%control%sgbcDispersive,eps0,mu0,this%control%simu_devia,this%control%stochastic)
               ! call Initsgbcs(sgg,this%media,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,this%control%layoutnumber,this%control%size, &
               !    this%g%G1,this%g%G2,this%g%GM1,this%g%GM2,this%thereAre%sgbcs,this%control%resume,this%control%sgbccrank,this%control%sgbcFreq,this%control%sgbcresol,this%control%sgbcdepth,this%control%sgbcDispersive,eps0,mu0,this%control%simu_devia,this%control%stochastic)

         l_auxinput= this%thereAre%sgbcs
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if (l_auxoutput) then
               write (dubuf,*) '----> there are Structured sgbc elements';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no Structured sgbc elements found';  call print11(this%control%layoutnumber,dubuf)
            endif
         endif
      end subroutine initializeSGBC
      
      subroutine initializeMultiports()
         character(len=BUFSIZE) :: dubuf
         logical :: l_auxinput, l_auxoutput

#ifdef CompileWithNIBC
         IF (this%control%mibc)  then
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init Multiports...';  call print11(this%control%layoutnumber,dubuf)
            call InitMultiports        (sgg,this%sggMiEx,this%sggMiEy,this%sggMiEz,this%sggMiHx ,this%sggMiHy ,this%sggMiHz,this%control%layoutnumber,this%control%size,this%thereAre%Multiports,this%control%resume, &
            Idxe,Idye,Idze,this%control%NOcompomur,this%control%AD,%this%control%cfl,eps0,mu0)
         l_auxinput= this%thereAre%Multiports
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if (l_auxoutput) then
               write (dubuf,*) '----> there are Structured  multiport elements';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no Structured multiport elements found';  call print11(this%control%layoutnumber,dubuf)
         endif
         endif
#endif
      end subroutine initializeMultiports

      subroutine initializeConformalElements()
         character(len=BUFSIZE) :: dubuf
         logical :: l_auxinput, l_auxoutput

#ifdef CompileWithConformal
         if(input_conformal_flag)then
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init Conformal Elements ...';  call print11(this%control%layoutnumber,dubuf)
!WIP
!DEBUG
            call initialize_memory_FDTD_conf_fields (sgg,this%sggMiEx, &
            & this%sggMiEy,this%sggMiEz,this%sggMiHx,this%sggMiHy,this%sggMiHz,Ex,Ey,Ez,Hx,Hy,Hz,&
            & this%control%layoutnumber,this%control%size, this%control%verbose);
            l_auxinput=input_conformal_flag
            l_auxoutput=l_auxinput
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         !       refactor JUN2015

         !!!!!!!sgg 051214 !rellena correctamente los campos magneticos. Necesario para construir los surfaces a partir del wireframe 
         !        call fillMagnetic(sgg, sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz, b)
         !!!!!!!ojo solo es valido para PEC!!!! cambiar luego !!?!?!?!?!?
            if (l_auxoutput ) then
               write (dubuf,*) '----> there are conformal elements';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no conformal elements found';  call print11(this%control%layoutnumber,dubuf)
            end if
      endif
#endif
      end subroutine initializeConformalElements

      subroutine initializeEDispersives()
         character (len=bufsize) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init EDispersives...';  call print11(this%control%layoutnumber,dubuf)
         call InitEDispersives(sgg,this%media,this%thereAre%EDispersives,this%control%resume,this%g%g1,this%g%g2,ex,ey,ez)
         l_auxinput=this%thereAre%EDispersives
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
            if (l_auxoutput ) then
               write (dubuf,*) '----> there are Structured Electric dispersive elements';  call print11(this%control%layoutnumber,dubuf)
            else
               write(dubuf,*) '----> no Structured Electric dispersive elements found';  call print11(this%control%layoutnumber,dubuf)
         endif
      end subroutine initializeEDispersives

      subroutine initializeMDispersives()
         character (len=bufsize) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init MDispersives...';  call print11(this%control%layoutnumber,dubuf)
         call InitMDispersives(sgg,this%media,this%thereAre%MDispersives,this%control%resume,this%g%gm1,this%g%gm2,hx,hy,hz)
         l_auxinput=this%thereAre%MDispersives
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         if ( l_auxoutput) then
             write (dubuf,*) '----> there are Structured Magnetic dispersive elements';  call print11(this%control%layoutnumber,dubuf)
         else
              write(dubuf,*) '----> no Structured Magnetic dispersive elements found';  call print11(this%control%layoutnumber,dubuf)
         endif
      end subroutine initializeMDispersives

      subroutine initializePlanewave()
         character (len=bufsize) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Multi Plane-Waves...';  call print11(this%control%layoutnumber,dubuf)
         call InitPlaneWave   (sgg,this%media,this%control%layoutnumber,this%control%size,SINPML_fullsize,this%thereAre%PlaneWaveBoxes,this%control%resume,eps0,mu0)
         l_auxinput=this%thereAre%PlaneWaveBoxes
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         if ( l_auxoutput) then
             write (dubuf,*) '----> there are Plane Wave';  call print11(this%control%layoutnumber,dubuf)
         else
              write(dubuf,*) '----> no Plane waves are found';  call print11(this%control%layoutnumber,dubuf)
         endif
      end subroutine initializePlanewave

      subroutine initializeNodalSources()
         character (len=bufsize) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Nodal Sources...';  call print11(this%control%layoutnumber,dubuf)
         if (.not.this%control%hopf) then
            call InitNodalSources(sgg,this%control%layoutnumber,sgg%NumNodalSources,sgg%NodalSource,sgg%Sweep,this%thereAre%NodalE,this%thereAre%NodalH)
         else
            call InitHopf(sgg,sgg%NumNodalSources,sgg%NodalSource,sgg%Sweep,this%control%ficherohopf) !lo manejara antonio con las entradas que precise
            this%thereAre%NodalE=.false. !no habra mas nodales excepto la de Hopf
            this%thereAre%NodalH=.false. 
         endif
         
         l_auxinput=this%thereAre%NodalH.or.this%thereAre%NodalE
         l_auxoutput=l_auxinput
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         if ( l_auxoutput) then
             write (dubuf,*) '----> there are Structured Nodal sources';  call print11(this%control%layoutnumber,dubuf)
         else
              write(dubuf,*) '----> no Structured Nodal sources are found';  call print11(this%control%layoutnumber,dubuf)
         endif

      end subroutine initializeNodalSources

      subroutine initializeObservation()
         character(len=bufsize) :: dubuf
         logical :: l_auxinput, l_auxoutput
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Observation...';  call print11(this%control%layoutnumber,dubuf)
         call InitObservation (sgg,this%media,tag_numbers, &
                              this%thereAre%Observation,this%thereAre%wires,this%thereAre%FarFields,this%control%resume,this%initialtimestep,this%control%finaltimestep,this%lastexecutedtime, &
                              this%control%nentradaroot,this%control%layoutnumber,this%control%size,this%control%saveall,this%control%singlefilewrite,this%control%wiresflavor,&
                              SINPML_FULLSIZE,this%control%facesNF2FF,this%control%NF2FFDecim,eps0,mu0,this%control%simu_devia,this%control%mpidir,this%control%niapapostprocess,this%bounds)

         l_auxinput=this%thereAre%Observation.or.this%thereAre%FarFields
         l_auxoutput=l_auxinput

#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( l_auxinput, l_auxoutput, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
#endif
         if (l_auxoutput ) then
               write (dubuf,*) '----> there are observation requests';  call print11(this%control%layoutnumber,dubuf)
         else
               write(dubuf,*) '----> no observation requests are found';  call print11(this%control%layoutnumber,dubuf)
         endif
      end subroutine initializeObservation

#ifdef CompileWithMPI
      subroutine initializeMPI()
         character(len=bufsize) :: dubuf      
         integer(kind=4) :: ierr
         if (this%control%size>1) then
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            write(dubuf,*) 'Init MPI MediaMatrix flush...';  call print11(this%control%layoutnumber,dubuf)
            call InitMPI(sgg%sweep,sgg%alloc)
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call InitExtraFlushMPI(this %control%layoutnumber,sgg%sweep,sgg%alloc,sgg%med,sgg%nummedia,this%media%sggMiEz,this%media%sggMiHz)
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call FlushMPI_H(sgg%alloc,this%control%layoutnumber,this%control%size, this%media%sggMiHx,this%media%sggMiHy,this%media%sggMiHz)
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call FlushMPI_E(sgg%alloc,this%control%layoutnumber,this%control%size, this%media%sggMiEx,this%media%sggMiEy,this%media%sggMiEz)
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
         endif

!!!!!!!!!!!!!!!!!!!!!fin juego con fuego 210815

      !MPI initialization
         if (this%control%size>1) then
            write(dubuf,*) 'Init MPI Cray...';  call print11(this%control%layoutnumber,dubuf)
            call InitMPI_Cray(this%control%layoutnumber,this%control%size,sgg%sweep,sgg%alloc, &
            sgg%Border%IsDownPeriodic,sgg%Border%IsUpPeriodic, &
            Ex,Ey,Ez,Hx,Hy,Hz)
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)

         !this modifies the initwires stuff and must be called after initwires (typically at the end)
         !llamalo siempre aunque no HAYA WIRES!!! para que no se quede colgado en hilos terminales
            if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
               (trim(adjustl(this%control%wiresflavor))=='transition') .or. & 
               this%control%use_mtln_wires) then
               write(dubuf,*) 'Init MPI Holland Wires...';  call print11(this%control%layoutnumber,dubuf)
               call newInitWiresMPI(this%control%layoutnumber,this%thereAre%wires,this%control%size,this%control%resume,sgg%sweep)
               call MPI_Barrier(SUBCOMM_MPI,ierr)
               write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
            endif

#ifdef CompileWithBerengerWires
            if (trim(adjustl(this%control%wiresflavor))=='berenger') then
               write(dubuf,*) 'Init MPI Multi-Wires...';  call print11(this%control%layoutnumber,dubuf)
               call InitWiresMPI_Berenger(this%control%layoutnumber,this%thereAre%wires,this%control%size,this%control%resume,sgg%sweep)
               call MPI_Barrier(SUBCOMM_MPI,ierr)
               write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
            endif
#endif
         !llamalo siempre para forzar los flush extra en caso de materiales anisotropos o multiport
            write(dubuf,*) 'Init Extra Flush MPI...';  call print11(this%control%layoutnumber,dubuf)
            call InitExtraFlushMPI_Cray(this%control%layoutnumber,sgg%sweep,sgg%alloc,sgg%Med,sgg%NumMedia,this%media%sggMiez,this%media%sggMiHz, &
            Ex,Ey,Ez,Hx,Hy,Hz,this%thereAre%MURBorders)
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
         endif

      
      !must be called now in case the MPI has changed the connectivity info
         if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
            (trim(adjustl(this%control%wiresflavor))=='transition')) then
            call ReportWireJunctions(this%control%layoutnumber,this%control%size,this%thereAre%wires,sgg%Sweep(iHz)%ZI, sgg%Sweep(iHz)%ZE,this%control%groundwires,this%control%strictOLD,this%control%verbose)
         endif

#ifdef CompileWithBerengerWires
      if (trim(adjustl(this%control%wiresflavor))=='berenger') then
               call ReportWireJunctionsBerenger(this%control%layoutnumber,this%control%size,this%thereAre%wires,sgg%Sweep(iHz)%ZI, sgg%Sweep(iHz)%ZE,this%control%groundwires,this%control%strictOLD,this%control%verbose)
                  !dama no tenia el equivalente 050416
      endif
#endif
#ifdef CompileWithSlantedWires
      if ((trim(adjustl(this%control%wiresflavor))=='slanted').or.(trim(adjustl(this%control%wiresflavor))=='semistructured')) then
         continue
      endif
#endif
      
      end subroutine initializeMPI
#endif

#ifdef CompileWithMPI
      subroutine flushMPIdata()
         integer(kind=4) :: ierr
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         !!Flush all the MPI data (needed a initial flush for correct resuming)
         if (this%control%size>1) then
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call   FlushMPI_H_Cray
         endif
         if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
            (trim(adjustl(this%control%wiresflavor))=='transition')) then
            if ((this%control%size>1).and.(this%thereAre%wires))   then
               call newFlushWiresMPI(this%control%layoutnumber,this%control%size)
            endif
#ifdef CompileWithStochastic
            if (this%control%stochastic) then
               call syncstoch_mpi_wires(this%control%simu_devia,this%control%layoutnumber,this%control%size)
            endif
#endif
         endif

#ifdef CompileWithBerengerWires
         if (trim(adjustl(this%control%wiresflavor))=='berenger') then
            if ((this%control%size>1).and.(this%thereAre%wires))   call FlushWiresMPI_Berenger(this%control%layoutnumber,this%control%size)
         endif
#endif
      end subroutine flushMPIdata
#endif

      subroutine printSimulationStart()
         character(len=bufsize) :: dubuf
         TYPE (tiempo_t) :: time_out2
#ifdef CompileWithMPI
         integer (kind=4) :: ierr
#endif

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         IF (this%control%resume) then
            write(dubuf,*)'END PREPROCESSING. RESUMING simulation from n=',this%n
            call print11(this%control%layoutnumber,dubuf)
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%control%layoutnumber,dubuf)
         else
            write(dubuf,*)'END PREPROCESSING. STARTING simulation.'
            call print11(this%control%layoutnumber,dubuf)
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%control%layoutnumber,dubuf)
#ifdef CompileWithMPI
            call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            CALL get_secnds (time_out2)
            write(dubuf,*)  'Start Date/time ', time_out2%fecha( 7: 8),'/',&
               &time_out2%fecha( 5: 6),'   ',time_out2%hora( 1: 2), ':',&
               &time_out2%hora( 3: 4),':',time_out2%hora( 5: 6)
            call print11(this%control%layoutnumber,dubuf)
            write(dubuf,*) SEPARADOR//separador//separador
            call print11(this%control%layoutnumber,dubuf)
         endif  
      end subroutine printSimulationStart

      subroutine fillMtag(sgg,sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz,sggMtag, b, tag_numbers)

         !------------------------>
         type (SGGFDTDINFO), intent(IN)    ::  sgg
         type (bounds_t), intent( IN)  ::  b
         INTEGER(KIND = IKINDMTAG), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( INOUT)     ::  sggMtag
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN   )     ::  sggMiHx
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN   )     ::  sggMiHy
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN   )     ::  sggMiHz
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN   )     ::  sggMiEx
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( IN   )     ::  sggMiEy
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN   )     ::  sggMiEz
         type (taglist_t) :: tag_numbers
         !------------------------> Variables locales
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio1,medio2,medio3,medio4,medio5
         logical  ::  mediois1,mediois2,mediois3,mediois4
         integer, dimension(3) :: lbx, lby, lbz
         lbx = lbound(tag_numbers%face%x)
         lby = lbound(tag_numbers%face%y)   
         lbz = lbound(tag_numbers%face%z)

         mediois3=.true.; mediois4=.true.
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4,medio5,mediois1,mediois2,mediois3,mediois4)
#endif
         Do k=1,b%sweepHx%NZ
            Do j=1,b%sweepHx%NY
               Do i=1,b%sweepHx%NX
                  medio1 =sggMiEy(i,j,k)
                  medio2 =sggMiEy(i,j,k+1)
                  medio3 =sggMiEz(i,j,k)
                  medio4 =sggMiEz(i,j+1,k)
                  medio5 =sggMiHx(i,j,k)
                  mediois1= (medio5==1).and.(medio1/=1).and.(medio2/=1).and.(medio3==1).and.(medio4==1)
                  mediois2= (medio5==1).and.(medio3/=1).and.(medio4/=1).and.(medio1==1).and.(medio2==1)
                  mediois3= .true. !.not.((medio5==1).and.(((sggMiHx(i-1,j,k)/=1).or.(sggMiHx(i+1,j,k)/=1)))) !esta condicion en realidad no detecta alabeos de una celda que siendo slots son acoples de un agujerito solo en el peor de los casos
                  if ((mediois1.or.mediois2).and.(mediois3))  then
                      !solo lo hace con celdas de vacio porque en particular el mismo medio sgbc con diferentes orientaciones tiene distintos indices de medio y lo activaria erroneamente si lo hago para todos los medios
                      tag_numbers%face%x(i+lbx(1)-1,j+lbx(2)-1,k+lbx(3)-1)=-ibset(iabs(tag_numbers%face%x(i+lbx(1)-1,j+lbx(2)-1,k+lbx(3)-1)),3) 
                      !ojo no cambiar: interacciona con observation tags 141020 !151020 a efectos de mapvtk el signo importa
                  endif
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4,medio5,mediois1,mediois2,mediois3,mediois4)
#endif
         Do k=1,b%sweepHy%NZ
            Do j=1,b%sweepHy%NY
               Do i=1,b%sweepHy%NX
                  medio1 =sggMiEz(i,j,k)
                  medio2 =sggMiEz(i+1,j,k)
                  medio3 =sggMiEx(i,j,k)
                  medio4 =sggMiEx(i,j,k+1)
                  medio5 =sggMiHy(i,j,k)
                  mediois1= (medio5==1).and.(medio1/=1).and.(medio2/=1).and.(medio3==1).and.(medio4==1)
                  mediois2= (medio5==1).and.(medio3/=1).and.(medio4/=1).and.(medio1==1).and.(medio2==1)
                  mediois3= .true. !.not.((medio5==1).and.(((sggMiHy(i,j-1,k)/=1).or.(sggMiHy(i,j+1,k)/=1))))
                  if ((mediois1.or.mediois2).and.(mediois3))  then
                     tag_numbers%face%y(i+lby(1)-1,j+lby(2)-1,k+lby(3)-1)=-ibset(iabs(tag_numbers%face%y(i+lby(1)-1,j+lby(2)-1,k+lby(3)-1)),4) 
                  endif
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4,medio5,mediois1,mediois2,mediois3,mediois4)
#endif
         Do k=1,b%sweepHz%NZ
            Do j=1,b%sweepHz%NY
               Do i=1,b%sweepHz%NX
                  medio1 =sggMiEx(i,j,k)
                  medio2 =sggMiEx(i,j+1,k)
                  medio3 =sggMiEy(i,j,k)
                  medio4 =sggMiEy(i+1,j,k)
                  medio5 =sggMiHz(i,j,k)
                  mediois1= (medio5==1).and.(medio1/=1).and.(medio2/=1).and.(medio3==1).and.(medio4==1)
                  mediois2= (medio5==1).and.(medio3/=1).and.(medio4/=1).and.(medio1==1).and.(medio2==1)
                  mediois3= .true. !.not.((medio5==1).and.(((sggMiHz(i,j,k-1)/=1).or.(sggMiHz(i,j,k+1)/=1))))
                  if ((mediois1.or.mediois2).and.(mediois3))  then
                     tag_numbers%face%z(i+lbz(1)-1,j+lbz(2)-1,k+lbz(3)-1)=-ibset(iabs(tag_numbers%face%z(i+lbz(1)-1,j+lbz(2)-1,k+lbz(3)-1)),5) 
                  endif
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine fillMtag

      subroutine crea_timevector(sgg,lastexecutedtimestep,finaltimestep,lastexecutedtime)
         integer (kind=4) :: lastexecutedtimestep,finaltimestep,i
         real (kind=RKIND_tiempo) :: lastexecutedtime
         type (SGGFDTDINFO), intent(INOUT)   ::  sgg
         allocate (sgg%tiempo(lastexecutedtimestep:finaltimestep+2))
         sgg%tiempo(lastexecutedtimestep)=lastexecutedtime
         do i=lastexecutedtimestep+1,finaltimestep+2
               sgg%tiempo(i)=sgg%tiempo(i-1)+sgg%dt !equiespaciados por defecto !luego los modifica prescale
         end do
         return
      end subroutine

   end subroutine solver_init

   subroutine solver_run(this, sgg, eps0, mu0, sinPML_fullsize, fullsize, tag_numbers, tagtype)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg

      real(kind=rkind), intent(inout) :: eps0,mu0

      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize, fullsize
      type(taglist_t), intent(in) :: tag_numbers
      type (tagtype_t), intent(in) :: tagtype

      real(kind=rkind), pointer, dimension (:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=rkind), pointer, dimension (:) :: Idxe, Idye, Idze, Idxh, Idyh, Idzh, dxe, dye, dze, dxh, dyh, dzh

      logical :: call_timing, l_aux, flushFF, somethingdone, newsomethingdone
      integer :: i
      real (kind=rkind) :: pscale_alpha
      REAL (kind=rkind_tiempo) :: at
      character(len=bufsize) :: dubuf
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
#ifdef CompileWithProfiling
      call nvtxStartRange("Antes del bucle N")
#endif
!240424 sgg creo el comunicador mpi de las sondas conformal aqui. debe irse con el nuevo conformal
#ifdef CompileWithConformal                
#ifdef CompileWithMPI
      call initMPIConformalProbes()
#endif  
#endif
      this%still_planewave_time=.true. !inicializacion de la variable 
      flushFF = .false.
      pscale_alpha=1.0 !se le entra con 1.0 

      Ex => this%Ex; Ey => this%Ey; Ez => this%Ez; Hx => this%Hx; Hy => this%Hy; Hz => this%Hz
      
      Idxe => this%Idxe; Idye => this%Idye; Idze => this%Idze; Idxh => this%Idxh; Idyh => this%Idyh; Idzh => this%Idzh; dxe => this%dxe; dye => this%dye; dze => this%dze; dxh => this%dxh; dyh => this%dyh; dzh => this%dzh


      ciclo_temporal :  DO while (this%n <= this%control%finaltimestep)
      
         call this%step(sgg, eps0, mu0, sinPML_fullsize, tag_numbers)
         call updateAndFlush()

         if(this%n >= this%n_info) then
             call_timing=.true.
         else
             call_timing=.false.
         endif
#ifdef CompileWithMPI
         l_aux=call_timing
         call MPI_AllReduce( l_aux, call_timing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr) !050619 incluido problemas stochastic stopflusing
#endif
         
         if (call_timing) then
            call Timing(sgg,this%bounds,this%n,this%n_info,this%control%layoutnumber,this%control%size, this%control%maxCPUtime,this%control%flushsecondsFields,this%control%flushsecondsData,this%initialtimestep, &
            this%control%finaltimestep,this%perform,this%parar,.FALSE., &
            Ex,Ey,Ez,this%everflushed,this%control%nentradaroot,this%control%maxSourceValue,this%control%opcionestotales,this%control%simu_devia,this%control%dontwritevtk,this%control%permitscaling)

            if (.not.this%parar) then !!! si es por parada se gestiona al final
!!!!! si esta hecho lo flushea todo pero poniendo de acuerdo a todos los mpi
                do i=1,sgg%NumberRequest
                   if  (sgg%Observation(i)%done.and.(.not.sgg%Observation(i)%flushed)) then
                      this%perform%flushXdmf=.true.
                      this%perform%flushVTK=.true.
                   endif
                end do
#ifdef CompileWithMPI
                l_aux=this%perform%flushVTK
                call MPI_AllReduce( l_aux, this%perform%flushVTK, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=this%perform%flushXdmf
                call MPI_AllReduce( l_aux, this%perform%flushXdmf, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=this%perform%flushDATA
                call MPI_AllReduce( l_aux, this%perform%flushDATA, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=this%perform%flushFIELDS
                call MPI_AllReduce( l_aux, this%perform%flushFIELDS, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=this%perform%postprocess
                call MPI_AllReduce( l_aux, this%perform%postprocess, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
#endif
!!!!!!!!!!!!
                if (this%perform%flushFIELDS) then
                   write(dubuf,*)  SEPARADOR,trim(adjustl(this%control%nentradaroot)),separador
                   call print11(this%control%layoutnumber,dubuf)
                   write(dubuf,*)  'INIT FLUSHING OF RESTARTING FIELDS n=',this%n
                   call print11(this%control%layoutnumber,dubuf)
                   call flush_and_save_resume(sgg, this%bounds, this%control%layoutnumber, this%control%size, this%control%nentradaroot, this%control%nresumeable2, this%thereare, this%n,eps0,mu0, this%everflushed,  &
                   Ex, Ey, Ez, Hx, Hy, Hz,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic)
#ifdef CompileWithMPI
                   call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
                   write(dubuf,*) SEPARADOR//separador//separador
                   call print11(this%control%layoutnumber,dubuf)
                   write(dubuf,*) 'DONE FLUSHING OF RESTARTING FIELDS n=',this%n
                   call print11(this%control%layoutnumber,dubuf)
                   write(dubuf,*) SEPARADOR//separador//separador
                   call print11(this%control%layoutnumber,dubuf)
                endif
                if (this%perform%isFlush()) then
                      !
                      flushFF=this%perform%postprocess
                      if (this%thereAre%FarFields.and.flushFF) then
                          write(dubuf,'(a,i9)')  ' INIT OBSERVATION DATA FLUSHING and Near-to-Far field n= ',this%n
                      else
                          write(dubuf,'(a,i9)')  ' INIT OBSERVATION DATA FLUSHING n= ',this%n
                      endif
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                      call print11(this%control%layoutnumber,dubuf)
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
    !!
                      if (this%thereAre%Observation) call FlushObservationFiles(sgg,this%ini_save, this%n,this%control%layoutnumber, this%control%size, dxe, dye, dze, dxh, dyh, dzh,this%bounds,this%control%singlefilewrite,this%control%facesNF2FF,flushFF)
                      !!
#ifdef CompileWithMPI
                      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
                      if (this%thereAre%FarFields.and.flushFF) then
                          write(dubuf,'(a,i9)')  ' Done OBSERVATION DATA FLUSHED and Near-to-Far field n= ',this%n
                      else
                          write(dubuf,'(a,i9)')  ' Done OBSERVATION DATA FLUSHED n= ',this%n
                      endif
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                      call print11(this%control%layoutnumber,dubuf)
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
    !
                      if (this%perform%postprocess) then
                         write(dubuf,'(a,i9)') 'Postprocessing frequency domain probes, if any, at n= ',this%n
                         call print11(this%control%layoutnumber,dubuf)
                         write(dubuf,*) SEPARADOR//separador//separador
                         call print11(this%control%layoutnumber,dubuf)
                         somethingdone=.false.
                         at=this%n*sgg%dt
                         if (this%thereAre%Observation) call PostProcessOnthefly(this%control%layoutnumber,this%control%size,sgg,this%control%nentradaroot,at,somethingdone,this%control%niapapostprocess,this%control%forceresampled)
#ifdef CompileWithMPI
                         call MPI_Barrier(SUBCOMM_MPI,ierr)
                         call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                         somethingdone=newsomethingdone
#endif
                         if (somethingdone) then
                           write(dubuf,*) 'End Postprocessing frequency domain probes.'
                           call print11(this%control%layoutnumber,dubuf)
                           write(dubuf,*) SEPARADOR//separador//separador
                           call print11(this%control%layoutnumber,dubuf)
                         else
                           write(dubuf,*) 'No frequency domain probes snapshots found to be postrocessed'
                           call print11(this%control%layoutnumber,dubuf)
                           write(dubuf,*) SEPARADOR//separador//separador
                           call print11(this%control%layoutnumber,dubuf)
                         endif
                      endif
                  !!       
                      if (this%perform%flushvtk) then   
                         write(dubuf,'(a,i9)')  ' Post-processing .vtk files n= ',this%n
                         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                         call print11(this%control%layoutnumber,dubuf)
                         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                         somethingdone=.false.
                         if (this%thereAre%Observation) call createvtkOnTheFly(this%control%layoutnumber,this%control%size,sgg,this%control%vtkindex,somethingdone,this%control%mpidir,tagtype,this%media%sggMtag,this%control%dontwritevtk)
#ifdef CompileWithMPI
                         call MPI_Barrier(SUBCOMM_MPI,ierr)
                         call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                         somethingdone=newsomethingdone
#endif
                          if (somethingdone) then
                                write(dubuf,*) 'End flushing .vtk snapshots'
                                call print11(this%control%layoutnumber,dubuf)
                                write(dubuf,*) SEPARADOR//separador//separador
                                call print11(this%control%layoutnumber,dubuf)
                          else
                                write(dubuf,*) 'No .vtk snapshots found to be flushed'
                                call print11(this%control%layoutnumber,dubuf)
                                write(dubuf,*) SEPARADOR//separador//separador
                                call print11(this%control%layoutnumber,dubuf)
                          endif
                      endif  
                         if (this%perform%flushXdmf) then
                            write(dubuf,'(a,i9)')  ' Post-processing .xdmf files n= ',this%n
                            call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                            call print11(this%control%layoutnumber,dubuf)
                            call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                            somethingdone=.false.

                            if (this%thereAre%Observation) call createxdmfOnTheFly(sgg,this%control%layoutnumber,this%control%size,this%control%vtkindex,this%control%createh5bin,somethingdone,this%control%mpidir)                          
                            if (this%control%createh5bin) call createh5bintxt(sgg,this%control%layoutnumber,this%control%size) !lo deben llamar todos haya on on this%thereAre%observation

#ifdef CompileWithMPI
                        call MPI_Barrier(SUBCOMM_MPI,ierr)
                        call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                        somethingdone=newsomethingdone
#endif
                            if (somethingdone) then
                                      write(dubuf,*) 'End flushing .xdmf snapshots'
                                      call print11(this%control%layoutnumber,dubuf)
                                      write(dubuf,*) SEPARADOR//separador//separador
                                      call print11(this%control%layoutnumber,dubuf)
                             else
                                      write(dubuf,*) 'No .xdmf snapshots found to be flushed'
                                      call print11(this%control%layoutnumber,dubuf)
                                      write(dubuf,*) SEPARADOR//separador//separador
                                      call print11(this%control%layoutnumber,dubuf)
                            endif
                      endif

#ifdef CompileWithMPI
                     call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
                 endif !del if (this%performflushDATA.or....
    !


                  if (this%control%singlefilewrite.and.this%perform%Unpack) call singleUnpack()
                  if ((this%control%singlefilewrite.and.this%perform%Unpack).or.this%perform%isFlush()) then
                     write(dubuf,'(a,i9)')  ' Continuing simulation at n= ',this%n
                     call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                     call print11(this%control%layoutnumber,dubuf)
                     call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                  endif

                endif !!!del if (.not.this%parar)
             endif !!!del if(n >= n_info
!          !!!!!!!!all the previous must be together
              
         this%control%fatalerror=.false.
         if (this%parar) then
             this%control%fatalerror=.true.
             exit ciclo_temporal
         endif
#ifdef CompileWithPrescale
         if (this%control%permitscaling) then
#ifndef miguelPscaleStandAlone
            if ((sgg%tiempo(this%n)>=EpsMuTimeScale_input_parameters%tini).and.&
                &(sgg%tiempo(this%n)<=EpsMuTimeScale_input_parameters%tend)) then
#endif
             call updateconstants(sgg,this%n,this%thereare,this%g%g1,this%g%g2,this%g%gM1,this%g%gM2, & 
                               Idxe,Idye,Idze,Idxh,Idyh,Idzh, &  !needed by  CPML to be updated
                               this%control%sgbc,this%control%mibc,input_conformal_flag, &
                               this%control%wiresflavor, this%control%wirecrank, this%control%fieldtotl,&
                               this%control%sgbcDispersive,this%control%finaltimestep, &
                               eps0,mu0, &
                               this%control%simu_devia, &
                               EpsMuTimeScale_input_parameters,pscale_alpha,this%still_planewave_time &
#ifdef CompileWithMPI
                               ,this%control%layoutnumber,this%control%size &
#endif
                               ,this%control%stochastic,this%control%verbose)
#ifndef miguelPscaleStandAlone
         endif
#endif
      endif
#endif
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !!!  Increase time step
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ! write(*write(*,*) 'timestepping: ', n
         this%n=this%n+1 !sube de iteracion
      end do ciclo_temporal ! End of the time-stepping loop

contains
      subroutine updateAndFlush()
         integer(kind=4) :: mindum
         IF (this%thereAre%Observation) then
            call UpdateObservation(sgg,this%media,tag_numbers, this%n,this%ini_save, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh,this%control%wiresflavor,SINPML_FULLSIZE,this%control%wirecrank, this%control%noconformalmapvtk,this%bounds)
            if (this%n>=this%ini_save+BuffObse)  then
               mindum=min(this%control%finaltimestep,this%ini_save+BuffObse)
               call FlushObservationFiles(sgg,this%ini_save,mindum,this%control%layoutnumber,this%control%size, dxe, dye, dze, dxh, dyh, dzh,this%bounds,this%control%singlefilewrite,this%control%facesNF2FF,.FALSE.) !no se flushean los farfields ahora
            endif
         endif
      end subroutine

      subroutine singleUnpack()
         character (LEN=BUFSIZE) :: dubuf
         logical :: somethingdone
         real (kind=rkind_tiempo) :: at
#ifdef CompileWithMPI
         integer(kind=4) :: ierr
#endif
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         write(dubuf,'(a,i9)')  ' Unpacking .bin files and prostprocessing them at n= ',this%n
         call print11(this%control%layoutnumber,dubuf)
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         if (this%thereAre%Observation) call unpacksinglefiles(sgg,this%control%layoutnumber,this%control%size,this%control%singlefilewrite,this%initialtimestep,this%control%resume) !dump the remaining to disk
         somethingdone=.false.
         if (this%control%singlefilewrite.and.this%perform%Unpack) then
            at=this%n*sgg%dt
            if (this%thereAre%Observation) call PostProcessOnthefly(this%control%layoutnumber,this%control%size,sgg,this%control%nentradaroot,at,somethingdone,this%control%niapapostprocess,this%control%forceresampled)
         endif
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
         somethingdone=newsomethingdone
#endif
         write(dubuf,'(a,i9)')  ' Done Unpacking .bin files and prostprocessing them at n= ',this%n
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         call print11(this%control%layoutnumber,dubuf)
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)

      end subroutine singleUnpack


   end subroutine solver_run



   subroutine step(this, sgg, eps0, mu0, sinPML_fullsize, tag_numbers)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      real(kind=rkind), intent(inout) :: eps0,mu0
      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize
      type(taglist_t), intent(in) :: tag_numbers

      logical :: planewave_switched_off = .false., thereareplanewave

#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif

      call flushPlanewaveOff(planewave_switched_off, this%still_planewave_time, thereareplanewave)
      call this%AdvanceAnisotropicE(sgg%alloc)
      call this%advanceE()
#ifdef CompileWithConformal
      if(this%control%input_conformal_flag) call conformal_advance_E()
#endif
      call this%advanceWiresE(sgg, eps0, mu0)
      call this%advancePMLE(sgg%NumMedia)
#ifdef CompileWithNIBC
      IF (this%thereAre%Multiports.and.(this%control%mibc)) call AdvanceMultiportE(sgg%alloc, this%Ex, this%Ey, this%Ez)
#endif
      call this%AdvancesgbcE(real(sgg%dt,RKIND))
      call this%advanceLumpedE(sgg)
      call this%advanceEDispersiveE(sgg)
      call this%advancePlaneWaveE(sgg)
      call this%advanceNodalE(sgg)

#ifdef CompileWithMPI
      if (this%control%size>1) then
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call FlushMPI_E_Cray
      endif
#endif

      call this%advanceAnisotropicH(sgg%alloc)
      call this%advanceH()
      call this%advancePMLbodyH()
      call this%AdvanceMagneticCPML(sgg%NumMedia)
      call this%MinusCloneMagneticPMC(sgg%alloc, sgg%border, sgg%sweep)
      call this%CloneMagneticPeriodic(sgg%alloc,sgg%Border, sgg%sweep)
      call this%AdvancesgbcH()
      call this%AdvanceMDispersiveH(sgg)
#ifdef CompileWithNIBC
      IF (this%thereAre%Multiports .and.(this%control%mibc))  &
         call AdvanceMultiportH (sgg%alloc,this%Hx,this%Hy,this%Hz, & 
                                 this%Ex,this%Ey,this%Ez,& 
                                 this%Idxe,this%Idye,this%Idze, & 
                                 this%sggMiHx,this%sggMiHy,this%sggMiHz, & 
                                 this%g%gm2,sgg%nummedia,this%control%conformalskin)
#endif
      call this%advancePlaneWaveH(sgg)
      call this%advanceNodalH(sgg)
      call this%advanceWiresH(sgg, eps0, mu0)
      call this%MinusCloneMagneticPMC(sgg%alloc, sgg%border, sgg%sweep)
      call this%CloneMagneticPeriodic(sgg%alloc, sgg%Border, sgg%sweep)

#ifdef CompileWithConformal                      
      if(this%control%input_conformal_flag) call conformal_advance_H()
#endif

#ifdef CompileWithMPI
      !!Flush all the MPI (esto estaba justo al principo del bucle temporal diciendo que era necesario para correcto resuming)
      !lo he movido aqui a 16/10/2012 porque el farfield necesita tener los campos magneticos correctos
      !e intuyo que el Bloque current tambien a tenor del comentario siguiente
      !Incluyo un flush inicial antes de entrar al bucle para que el resuming sea correcto
      if (this%control%size>1) then
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         call FlushMPI_H_Cray
      endif
      if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
            (trim(adjustl(this%control%wiresflavor))=='transition')) then
         if ((this%control%size>1).and.(this%thereAre%wires)) call newFlushWiresMPI(this%control%layoutnumber,this%control%size)
#ifdef CompileWithStochastic
         if (this%control%stochastic) call syncstoch_mpi_wires(this%control%simu_devia,this%control%layoutnumber,this%control%size)
#endif
      endif
#ifdef CompileWithBerengerWires
      if (trim(adjustl(this%control%wiresflavor))=='berenger') then
         if ((this%control%size>1).and.(this%thereAre%wires))   call FlushWiresMPI_Berenger(this%control%layoutnumber,this%control%size)
      endif
#endif
#endif

!!!no se si el orden wires - sgbcs del sync importa 150519
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
         if (this%control%stochastic) call syncstoch_mpi_sgbcs(this%control%simu_devia,this%control%layoutnumber,this%control%size)
#endif    
#endif

#ifdef CompileWithMPI
#ifdef CompileWithStochastic
         if (this%control%stochastic) call syncstoch_mpi_lumped(this%control%simu_devia,this%control%layoutnumber,this%control%size)
#endif    
#endif 
      call this%advanceMagneticMUR(sgg)
contains

   subroutine flushPlanewaveOff(pw_switched_off, pw_still_time, pw_thereAre)
      logical, intent(inout) :: pw_switched_off, pw_still_time, pw_thereAre
      logical :: pw_still_time_aux, pw_thereAre_aux
      integer (kind=4) :: ierr
      character(len=bufsize) :: dubuf
      if (.not.pw_switched_off) then
         pw_still_time = pw_still_time.and.this%thereAre%PlaneWaveBoxes
         pw_thereAre = this%thereAre%PlaneWaveBoxes
#ifdef CompileWithMPI
         if (this%control%size>1) then
            pw_still_time_aux = pw_still_time
            call MPI_AllReduce(pw_still_time_aux, pw_still_time, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
            pw_thereAre_aux = pw_thereAre
            call MPI_AllReduce(pw_thereAre_aux, pw_thereAre, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
         endif
#endif
         if (.not.pw_still_time)  then
            pw_switched_off=.true.
            write(dubuf,*) 'Switching plane-wave off at n=', this%n
            if (pw_thereAre) call print11(this%control%layoutnumber,dubuf)
         endif
      endif
   end subroutine 


!       !PML E-field advancing (IT IS IMPORTANT TO FIRST CALL THE PML ADVANCING ROUTINES, SINCE THE DISPERSIVE
!       !ROUTINES INJECT THE POLARIZATION CURRENTS EVERYWHERE (PML INCLUDED)
!       !SO THAT DISPERSIVE MATERIALS CAN ALSO BE TRUNCATED BY CPML)

#ifdef CompileWithMPI
   subroutine initMPIConformalProbes()
      integer (kind=4) :: group_conformalprobes_dummy, ierr
!!!!sgg250424 niapa para que funcionen sondas conformal mpi
!todos deben crear el subcomunicador mpi una sola vez   
      if (input_conformal_flag) then
         SUBCOMM_MPI_conformal_probes=1   
         MPI_conformal_probes_root=this%control%layoutnumber
      else  
         SUBCOMM_MPI_conformal_probes=0 
         MPI_conformal_probes_root=-1
      endif
      call MPIinitSubcomm(this%control%layoutnumber,this%control%size,SUBCOMM_MPI_conformal_probes,&
                           MPI_conformal_probes_root,group_conformalprobes_dummy)
      ! print *,'-----creating--->',this%control%layoutnumber,SIZE,SUBCOMM_MPI_conformal_probes,MPI_conformal_probes_root
      call MPI_BARRIER(SUBCOMM_MPI, ierr)
      !!!no lo hago pero al salir deberia luego destruir el grupo call MPI_Group_free(output(ii)%item(i)%MPIgroupindex,ierr)                   
   end subroutine initMPIConformalProbes
#endif

   end subroutine step

   subroutine advanceE(this)
      class(solver_t) :: this
#ifdef CompileWithProfiling
      call nvtxStartRange("Antes del bucle EX")
#endif
      call this%advanceEx(this%media%sggMiEx)
#ifdef CompileWithProfiling
      call nvtxEndRange

      call nvtxStartRange("Antes del bucle EY")
#endif
      call this%advanceEy(this%media%sggMiEy)
      
#ifdef CompileWithProfiling    
      call nvtxEndRange

      call nvtxStartRange("Antes del bucle EZ")
#endif
      call this%advanceEz(this%media%sggMiEz)
#ifdef CompileWithProfiling    
      call nvtxEndRange
#endif
   end subroutine

   subroutine advanceEx(this, sggMiEx)
      class(solver_t) :: this
      integer(kind=integersizeofmediamatrices), dimension(0:this%bounds%sggMiEx%NX-1,0:this%bounds%sggMiEx%NY-1,0:this%bounds%sggMiEx%NZ-1), intent(in) :: sggMiEx

      real(kind=rkind), dimension(:,:,:), pointer, contiguous ::  Ex
      real(kind=rkind), dimension(:,:,:), pointer, contiguous ::  Hy
      real(kind=rkind), dimension(:,:,:), pointer, contiguous ::  Hz
      real(kind=rkind), dimension(:), pointer :: Idyh
      real(kind=rkind), dimension(:), pointer :: Idzh

      real(kind=rkind) :: Idzhk, Idyhj
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Ex(0:this%bounds%Ex%NX-1,0:this%bounds%Ex%NY-1,0:this%bounds%Ex%NZ-1) => this%Ex
      Hy(0:this%bounds%Hy%NX-1,0:this%bounds%Hy%NY-1,0:this%bounds%Hy%NZ-1) => this%Hy
      Hz(0:this%bounds%Hz%NX-1,0:this%bounds%Hz%NY-1,0:this%bounds%Hz%NZ-1) => this%Hz

      Idyh(0:this%bounds%dyh%NY-1) => this%Idyh
      Idzh(0:this%bounds%dzh%NZ-1) => this%Idzh

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzhk,Idyhj) 
#endif
#ifdef CompileWithACC   
!$ACC parallel loop DEFAULT(present) collapse (2) private (i,j,k,medio,Idzhk,Idyhj)  copyin(Ex,sggMiEx,Hy,Hz,Idyh,Idzh,b,G1,G2) copyout(Ex) 
#endif
      Do k=1,this%bounds%sweepEx%NZ
         Do j=1,this%bounds%sweepEx%NY
            Do i=1,this%bounds%sweepEx%NX
               Idzhk=Idzh(k)
               Idyhj=Idyh(j)
               medio =sggMiEx(i,j,k)
               Ex(i,j,k)=this%g%g1(MEDIO)*Ex(i,j,k)+this%g%g2(MEDIO)* &
               ((Hz(i,j,k)-Hz(i,j-1,k))*Idyhj-(Hy(i,j,k)-Hy(i,j,k-1))*Idzhk)
            End do
         End do
      End do
#ifdef CompileWithOpenMP   
!$OMP  END PARALLEL DO
#endif
   end subroutine advanceEx

   subroutine advanceEy(this,sggMiEy)
      class(solver_t) :: this
      integer(kind=integersizeofmediamatrices), dimension (0:this%bounds%sggMiEy%NX-1,0:this%bounds%sggMiEy%NY-1,0:this%bounds%sggMiEy%NZ-1), intent(in) :: sggMiEy

      real(kind=rkind), dimension(:,:,:), pointer, contiguous ::  Ey
      real(kind=rkind), dimension(:,:,:), pointer, contiguous ::  Hz
      real(kind=rkind), dimension(:,:,:), pointer, contiguous ::  Hx
      real(kind=rkind), dimension(:), pointer :: Idzh
      real(kind=rkind), dimension(:), pointer :: Idxh

      real (kind=rkind) :: Idzhk
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Ey(0:this%bounds%Ey%NX-1,0:this%bounds%Ey%NY-1,0:this%bounds%Ey%NZ-1) => this%Ey
      Hz(0:this%bounds%Hz%NX-1,0:this%bounds%Hz%NY-1,0:this%bounds%Hz%NZ-1) => this%Hz
      Hx(0:this%bounds%Hx%NX-1,0:this%bounds%Hx%NY-1,0:this%bounds%Hx%NZ-1) => this%Hx

      Idzh(0:this%bounds%dzh%NZ-1) => this%Idzh
      Idxh(0:this%bounds%dxh%NX-1) => this%Idxh

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzhk)  
#endif
#ifdef CompileWithACC   
!$ACC parallel loop  DEFAULT(present) collapse (2) private (i,j,k,medio,Idzhk)     copyin(Ey,sggMiEy,Hz,Hx,Idzh,Idxh,b,G1,G2) copyout(Ey) 
#endif
      Do k=1,this%bounds%sweepEy%NZ
         Do j=1,this%bounds%sweepEy%NY
            Do i=1,this%bounds%sweepEy%NX
               Idzhk=Idzh(k)
               medio =sggMiEy(i,j,k)
               Ey(i,j,k)=this%g%g1(MEDIO)*Ey(i,j,k)+this%g%g2(MEDIO)*((Hx(i,j,k)-Hx(i,j,k-1))*Idzhk-(Hz(i,j,k)-Hz(i-1,j,k))*Idxh(i))
            End do
         End do
      End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif



      return
   end subroutine advanceEy

   subroutine advanceEz(this,sggMiEz)
      class(solver_t) :: this
      integer(kind=integersizeofmediamatrices), dimension(0:this%bounds%sggMiEz%NX-1,0:this%bounds%sggMiEz%NY-1,0:this%bounds%sggMiEz%NZ-1), intent(in) ::  sggMiEz

      real (kind=rkind), dimension(:,:,:), pointer, contiguous :: Ez
      real (kind=rkind), dimension(:,:,:), pointer, contiguous :: Hx
      real (kind=rkind), dimension(:,:,:), pointer, contiguous :: Hy
      real (kind=rkind), dimension(:), pointer ::  Idyh
      real (kind=rkind), dimension(:), pointer ::  Idxh
      !------------------------> Variables locales
      real (kind = RKIND)  ::   Idyhj
      integer(kind = 4)  ::  i, j, k
      integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio


      Ez(0:this%bounds%Ez%NX-1,0:this%bounds%Ez%NY-1,0:this%bounds%Ez%NZ-1) => this%Ez
      Hx(0:this%bounds%HX%NX-1,0:this%bounds%HX%NY-1,0:this%bounds%HX%NZ-1) => this%Hx
      Hy(0:this%bounds%Hy%NX-1,0:this%bounds%Hy%NY-1,0:this%bounds%Hy%NZ-1) => this%Hy

      Idyh(0:this%bounds%dyh%NY-1) => this%Idyh
      Idxh(0:this%bounds%dxh%NX-1) => this%Idxh

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idyhj)    
#endif
#ifdef CompileWithACC   
!$ACC parallel loop   DEFAULT(present) collapse (2) private (i,j,k,medio,Idyhj)        copyin(Ez,sggMiEz,Hx,Hy,Idxh,Idyh,b,G1,G2) copyout(Ez) 
#endif
      Do k=1,this%bounds%sweepEz%NZ
         Do j=1,this%bounds%sweepEz%NY
            Do i=1,this%bounds%sweepEz%NX
               Idyhj=Idyh(j)
               medio =sggMiEz(i,j,k)
               Ez(i,j,k)=this%g%g1(MEDIO)*Ez(i,j,k)+this%g%g2(MEDIO)*((Hy(i,j,k)-Hy(i-1,j,k))*Idxh(i)-(Hx(i,j,k)-Hx(i,j-1,k))*Idyhj)
            End do
         End do
      End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
      return
   end subroutine advanceEz

   subroutine advanceH(this)
      class(solver_t) :: this
#ifdef CompileWithProfiling    
      call nvtxStartRange("Antes del bucle HX")
#endif
      call this%advanceHx(this%media%sggMiHx)
#ifdef CompileWithProfiling    
      call nvtxEndRange
      call nvtxStartRange("Antes del bucle HY")
#endif
      call this%advanceHy(this%media%sggMiHy)
#ifdef CompileWithProfiling    
      call nvtxEndRange
      call nvtxStartRange("Antes del bucle HZ")
#endif
      call this%advanceHz(this%media%sggMiHz)  
#ifdef CompileWithProfiling    
      call nvtxEndRange
#endif
   end subroutine advanceH

   subroutine advanceHx(this, sggMiHx)
      class(solver_t) :: this
      integer(kind=integersizeofmediamatrices), dimension(0:this%bounds%sggMiHx%NX-1,0:this%bounds%sggMiHx%NY-1,0:this%bounds%sggMiHx%NZ-1), intent(in) :: sggMiHx

      real (kind=rkind), dimension(:,:,:), pointer, contiguous  ::  Hx
      real (kind=rkind), dimension(:,:,:), pointer, contiguous  ::  Ey
      real (kind=rkind), dimension(:,:,:), pointer, contiguous  ::  Ez
      real (kind=rkind), dimension(:), pointer:: IdyE
      real (kind=rkind), dimension(:), pointer:: IdzE

      real (kind=rkind) :: Idzek, Idyej
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Hx(0:this%bounds%Hx%NX-1,0:this%bounds%Hx%NY-1,0:this%bounds%Hx%NZ-1) => this%Hx
      Ey(0:this%bounds%Ey%NX-1,0:this%bounds%Ey%NY-1,0:this%bounds%Ey%NZ-1) => this%Ey
      Ez(0:this%bounds%Ez%NX-1,0:this%bounds%Ez%NY-1,0:this%bounds%Ez%NZ-1) => this%Ez

      IdyE(0:this%bounds%dyE%NY-1) => this%IdyE
      IdzE(0:this%bounds%dzE%NZ-1) => this%IdzE

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzek,Idyej)     
#endif
#ifdef CompileWithACC   
!$ACC parallel loop  DEFAULT(present) collapse (2) private (i,j,k,medio,Idzek,Idyej)       copyin(Hx,sggMiHx,Ey,Ez,Idye,Idze,b,GM1,GM2) copyout(Hx) 
#endif
      Do k=1,this%bounds%sweepHx%NZ
         Do j=1,this%bounds%sweepHx%NY
            Do i=1,this%bounds%sweepHx%NX
            Idzek=Idze(k)
            Idyej=Idye(j)
               medio =sggMiHx(i,j,k)
               Hx(i,j,k)=this%g%gm1(medio)*Hx(i,j,k)+this%g%gm2(medio)*((Ey(i,j,k+1)-Ey(i,j,k))*Idzek-(Ez(i,j+1,k)-Ez(i,j,k))*Idyej)
            End do
         End do
      End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
      return
   end subroutine advanceHx

   subroutine advanceHy(this, sggMiHy)
      class(solver_t) :: this
      integer(kind=integersizeofmediamatrices), dimension(0:this%bounds%sggMiHy%NX-1,0:this%bounds%sggMiHy%NY-1,0:this%bounds%sggMiHy%NZ-1), intent(in) :: sggMiHy

      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hy
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ez
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ex
      real(kind=rkind), dimension(:), pointer :: IdzE
      real(kind=rkind), dimension(:), pointer :: IdxE

      real (kind=rkind) :: Idzek
      integer(kind=4) :: i, j, k
      integer(kind=integersizeofmediamatrices) :: medio

      Hy(0:this%bounds%Hy%NX-1,0:this%bounds%Hy%NY-1,0:this%bounds%Hy%NZ-1) => this%Hy
      Ez(0:this%bounds%Ez%NX-1,0:this%bounds%Ez%NY-1,0:this%bounds%Ez%NZ-1) => this%Ez
      Ex(0:this%bounds%Ex%NX-1,0:this%bounds%Ex%NY-1,0:this%bounds%Ex%NZ-1) => this%Ex

      IdzE(0:this%bounds%dzE%NZ-1) => this%IdzE
      IdxE(0:this%bounds%dxE%NX-1) => this%IdxE

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzek)     
#endif
#ifdef CompileWithACC   
!$ACC parallel loop DEFAULT(present) collapse (2) private (i,j,k,medio,Idzek)         copyin(Hy,sggMiHy,Ez,Ex,Idze,Idxe,b,GM1,GM2) copyout(Hy) 
#endif
      Do k=1,this%bounds%sweepHy%NZ
         Do j=1,this%bounds%sweepHy%NY
            Do i=1,this%bounds%sweepHy%NX
               Idzek=Idze(k)
               medio =sggMiHy(i,j,k)
               Hy(i,j,k)=this%g%gm1(medio)*Hy(i,j,k)+this%g%gm2(medio)*((Ez(i+1,j,k)-Ez(i,j,k))*Idxe(i)-(Ex(i,j,k+1)-Ex(i,j,k))*Idzek)
            End do
         End do
      End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
      return
   end subroutine advanceHy

   subroutine advanceHz(this, sggMiHz)
      class(solver_t) :: this
      integer(kind=integersizeofmediamatrices), dimension(0:this%bounds%sggMiHz%NX-1,0:this%bounds%sggMiHz%NY-1,0:this%bounds%sggMiHz%NZ-1), intent(in) :: sggMiHz
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Hz
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ex
      real(kind=rkind), dimension(:,:,:), pointer, contiguous :: Ey
      real(kind=rkind), dimension(:), pointer :: IdyE
      real(kind=rkind), dimension(:), pointer :: IdxE

      real (kind = RKIND)  ::  Idyej
      integer(kind = 4)  ::  i, j, k
      integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
      Hz(0:this%bounds%Hz%NX-1,0:this%bounds%Hz%NY-1,0:this%bounds%Hz%NZ-1) => this%Hz
      Ex(0:this%bounds%EX%NX-1,0:this%bounds%EX%NY-1,0:this%bounds%EX%NZ-1) => this%Ex
      Ey(0:this%bounds%Ey%NX-1,0:this%bounds%Ey%NY-1,0:this%bounds%Ey%NZ-1) => this%Ey
      IdyE(0:this%bounds%dyE%NY-1) => this%IdyE
      IdxE(0:this%bounds%dxE%NX-1) => this%IdxE

#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idyej)  
#endif
#ifdef CompileWithACC   
!$ACC parallel loop  DEFAULT(present) collapse (2) private (i,j,k,medio,Idyej)       copyin(Hz,sggMiHz,Ex,Ey,Idxe,Idye,b,GM1,GM2) copyout(Hz)
#endif
      Do k=1,this%bounds%sweepHz%NZ
         Do j=1,this%bounds%sweepHz%NY
            Do i=1,this%bounds%sweepHz%NX
               Idyej=Idye(j)
               medio =sggMiHz(i,j,k)
               Hz(i,j,k)=this%g%gm1(medio)*Hz(i,j,k)+this%g%gm2(medio)*((Ex(i,j+1,k)-Ex(i,j,k))*Idyej-(Ey(i+1,j,k)-Ey(i,j,k))*Idxe(i))
            End do
         End do
      End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
      return
   end subroutine advanceHz


   subroutine solver_advanceEDispersiveE(this, sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      if (this%thereAre%Edispersives) call AdvanceEDispersiveE(sgg)
   end subroutine

   subroutine solver_advanceMDispersiveH(this, sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      if (this%thereAre%Mdispersives) call AdvanceMDispersiveH(sgg)
   end subroutine

   subroutine solver_advanceLumpedE(this, sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      if (this%thereAre%Lumpeds) call AdvanceLumpedE(sgg,this%n,this%control%simu_devia,this%control%stochastic)
   end subroutine


   subroutine solver_advancePlaneWaveE(this,sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      If (this%thereAre%PlaneWaveBoxes.and.this%still_planewave_time) then 
         if(.not.this%control%simu_devia) call AdvancePlaneWaveE(sgg,this%n, this%bounds,this%g%G2, &
                                                                 this%Idxh,this%Idyh,this%Idzh, & 
                                                                 this%Ex,this%Ey,this%Ez, & 
                                                                 this%still_planewave_time)
      end if
   end subroutine

   subroutine solver_advancePlaneWaveH(this,sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      If (this%thereAre%PlaneWaveBoxes.and.this%still_planewave_time)  then
         if (.not.this%control%simu_devia) call AdvancePlaneWaveH(sgg,this%n, this%bounds, this%g%GM2, & 
                                                                  this%Idxe, this%Idye, this%Idze, & 
                                                                  this%Hx, this%Hy, this%Hz, & 
                                                                  this%still_planewave_time)
      endif
   end subroutine

   subroutine solver_advanceNodalE(this, sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
         if (this%thereAre%NodalE) then 
            call advanceNodalE(sgg,this%media%sggMiEx,this%media%sggMiEy,this%media%sggMiEz,& 
                               sgg%NumMedia,this%n, this%bounds, this%g%G2,& 
                               this%Idxh,this%Idyh,this%Idzh,&
                               this%Ex,this%Ey,this%Ez,&
                               this%control%simu_devia)
         end if
   end subroutine

   subroutine solver_advanceNodalH(this, sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      if (this%thereAre%NodalH) then 
         call AdvanceNodalH(sgg,this%media%sggMiHx,this%media%sggMiHy,this%media%sggMiHz,&
                            sgg%NumMedia,this%n, this%bounds,this%g%GM2, & 
                            this%Idxe,this%Idye,this%Idze, & 
                            this%Hx,this%Hy,this%Hz,&
                            this%control%simu_devia)
      end if
   end subroutine

   subroutine solver_advanceAnisotropicE(this, alloc)
      class(solver_t) :: this
      type(XYZlimit_t), dimension (1:6), intent(in) :: alloc
      if (this%thereAre%Anisotropic) call AdvanceAnisotropicE(alloc,this%ex,this%ey,this%ez, & 
                                                             this%hx, this%hy, this%hz, & 
                                                             this%Idxe, this%Idye, this%Idze, & 
                                                             this%Idxh, this%Idyh, this%Idzh)
   end subroutine

   subroutine solver_advanceAnisotropicH(this, alloc)
      class(solver_t) :: this
      type(XYZlimit_t), dimension (1:6), intent(in) :: alloc
      IF (this%thereAre%Anisotropic) call AdvanceAnisotropicH(alloc, this%ex, this%ey, this%ez, & 
                                                              this%hx, this%hy, this%hz, & 
                                                              this%Idxe, this%Idye, this%Idze, &
                                                              this%Idxh, this%Idyh, this%Idzh)
   end subroutine

   subroutine solver_advancePMLbodyH(this)
      class(solver_t) :: this
      if (this%thereAre%PMLbodies) call AdvancePMLbodyH()
   end subroutine

   subroutine solver_advanceMagneticCPML(this, numMedia)
      class(solver_t) :: this
      integer(kind=4), intent(in) :: numMedia
      If (this%thereAre%PMLBorders) call advanceMagneticCPML(numMedia, this%bounds, & 
                                                             this%media%sggMiHx, this%media%sggMiHy, this%media%sggMiHz, & 
                                                             this%g%gm2, this%Hx, this%Hy, this%Hz, & 
                                                             this%Ex, this%Ey, this%Ez)
   end subroutine

   subroutine solver_MinusCloneMagneticPMC(this, alloc, border, sweep)
      class(solver_t) :: this
      type(XYZlimit_t), dimension (1:6), intent(in) :: alloc
      type (Border_t), intent(in) :: border
      type (XYZlimit_t), dimension(1:6) :: sweep
      If (this%thereAre%PMCBorders) call MinusCloneMagneticPMC(alloc,border,this%Hx,this%Hy,this%Hz,sweep, & 
                                                               this%control%layoutnumber,this%control%size)
   end subroutine

   subroutine solver_CloneMagneticPeriodic(this, alloc, border, sweep)
      class(solver_t) :: this
      type(XYZlimit_t), dimension (1:6), intent(in) :: alloc
      type (Border_t), intent(in) :: border
      type (XYZlimit_t), dimension(1:6) :: sweep
      If (this%thereAre%PeriodicBorders) call CloneMagneticPeriodic(alloc,border,this%Hx,this%Hy,this%Hz,sweep,& 
                                                                    this%control%layoutnumber,this%control%size)
   end subroutine


   subroutine solver_advancePMLE(this, numMedia)
      class (solver_t) :: this
      integer (kind=4) :: numMedia
      If (this%thereAre%PMLbodies) then !waveport absorbers
         call AdvancePMLbodyE()
      endif
      If (this%thereAre%PMLBorders) then
         call AdvanceelectricCPML(numMedia, this%bounds,this%media%sggMiEx,this%media%sggMiEy,this%media%sggMiEz, & 
                                  this%g%G2, this%Ex, this%Ey, this%Ez, this%Hx, this%Hy, this%Hz)
      endif
   end subroutine

   subroutine solver_advancesgbcE(this, dt)
      class(solver_t) :: this
      real(kind=rkind), intent(in) :: dt

      if (this%thereAre%sgbcs.and.(this%control%sgbc)) then 
         call AdvancesgbcE(dt,this%control%sgbcDispersive, & 
                                this%control%simu_devia,this%control%stochastic)
      end if
   end subroutine

   subroutine solver_advancesgbcH(this)
      class(solver_t) :: this
      if (this%thereAre%sgbcs.and.(this%control%sgbc)) call AdvancesgbcH()
   end subroutine

   subroutine solver_advanceWiresE(this, sgg, eps0, mu0)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      real(kind=rkind), intent(inout) :: eps0,mu0
      character(len=bufsize) :: buff

      if (( (trim(adjustl(this%control%wiresflavor))=='holland') .or. &
            (trim(adjustl(this%control%wiresflavor))=='transition')) .and. .not. this%control%use_mtln_wires) then
         IF (this%thereAre%Wires) then
            if (this%control%wirecrank) then
               call AdvanceWiresEcrank(sgg, this%n, this%control%layoutnumber,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic)
            else
#ifdef CompileWithMTLN
               if (this%mtln_parsed%has_multiwires) then
                  write(buff, *) 'ERROR: Multiwires in simulation but -mtlnwires flag has not been selected'
                  call WarnErrReport(buff)
               end if
#endif
               call AdvanceWiresE(sgg,this%n, this%control%layoutnumber,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic,this%control%experimentalVideal,this%control%wirethickness,eps0,mu0)
            endif
         endif
      endif
#ifdef CompileWithBerengerWires
      if (trim(adjustl(this%control%wiresflavor))=='berenger') then
         IF (this%thereAre%Wires) call AdvanceWiresE_Berenger(sgg,n)
      endif
#endif
#ifdef CompileWithSlantedWires
      if((trim(adjustl(this%control%wiresflavor))=='slanted').or.(trim(adjustl(this%control%wiresflavor))=='semistructured')) then
         call AdvanceWiresE_Slanted(sgg,n) 
      endif
#endif
      if (this%control%use_mtln_wires) then
#ifdef CompileWithMTLN
         call AdvanceWiresE_mtln(sgg,this%Idxh,this%Idyh,this%Idzh,eps0,mu0)
#else
         write(buff,'(a)') 'WIR_ERROR: Executable was not compiled with MTLN modules.'
#endif   
      end if

   end subroutine

   subroutine solver_advancewiresH(this, sgg, eps0, mu0)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      real(kind=rkind), intent(inout) :: eps0,mu0
      if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
            (trim(adjustl(this%control%wiresflavor))=='transition')) then
         IF (this%thereAre%Wires) then
            if (this%control%wirecrank) then
               continue
            else
               call AdvanceWiresH(sgg,this%n, this%control%layoutnumber,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic,this%control%experimentalVideal,this%control%wirethickness,eps0,mu0)
            endif
         endif
      endif

   end subroutine

   subroutine solver_advanceMagneticMUR(this, sgg)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
#ifdef CompileWithMPI
      integer(kind=4) :: ierr
#endif
      If (this%thereAre%MURBorders) then
         call AdvanceMagneticMUR(this%bounds, sgg, & 
                                 this%media%sggMiHx, this%media%sggMiHy, this%media%sggMiHz, &
                                 this%Hx, this%Hy, this%Hz, & 
                                 this%control%mur_second)
#ifdef CompileWithMPI
         if (this%control%mur_second) then
            if (this%control%size>1) then
               call MPI_Barrier(SUBCOMM_MPI,ierr)
               call FlushMPI_H_Cray
            endif
         endif
#endif
      endif
   end subroutine


   subroutine solver_end(this, sgg, eps0, mu0, tagtype, finishedwithsuccess)
      class(solver_t) :: this
      type(sggfdtdinfo), intent(in) :: sgg
      real(kind=rkind), intent(in) :: eps0,mu0
      type (tagtype_t), intent(in) :: tagtype
      logical, intent(inout) :: finishedwithsuccess

      real(kind=rkind), pointer, dimension (:,:,:) :: Ex, Ey, Ez, Hx, Hy, Hz
      real(kind=rkind), pointer, dimension (:) :: dxe, dye, dze, dxh, dyh, dzh
      real(kind=rkind_tiempo) :: at
      integer (kind=4) :: ndummy
      logical :: dummylog, somethingdone, newsomethingdone
      character(len=bufsize) :: dubuf

#ifdef CompileWithMPI
      integer (kind=4) :: ierr
#endif
      Ex => this%Ex; Ey => this%Ey; Ez => this%Ez; Hx => this%Hx; Hy => this%Hy; Hz => this%Hz;
      dxe => this%dxe; dye => this%dye; dze => this%dze; dxh => this%dxh; dyh => this%dyh; dzh => this%dzh

#ifdef CompileWithProfiling
      call nvtxEndRange
#endif      
      
#ifdef CompileWithConformal
      if(this%control%input_conformal_flag) call conformal_final_simulation  (conf_timeSteps, this%n)
#endif

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      if (this%n>this%control%finaltimestep) this%n = this%control%finaltimestep !readjust n since after finishing it is increased
      this%control%finaltimestep = this%n
      this%lastexecutedtime=sgg%tiempo(this%control%finaltimestep)
      !se llama con dummylog para no perder los flags de parada
      call Timing(sgg,this%bounds,this%n,ndummy,this%control%layoutnumber, this%control%size, & 
                  this%control%maxCPUtime,this%control%flushsecondsFields, this%control%flushsecondsData, &
                  this%initialtimestep, this%control%finaltimestep,this%d_perform,dummylog,.FALSE., &
                  Ex,Ey,Ez,this%everflushed,this%control%nentradaroot,this%control%maxSourceValue, & 
                  this%control%opcionestotales,this%control%simu_devia,this%control%dontwritevtk,this%control%permitscaling)

      write(dubuf,*)'END FDTD time stepping. Beginning posprocessing at n= ',this%n
      call print11(this%control%layoutnumber,dubuf)

      if ((this%control%flushsecondsFields/=0).or.this%perform%flushFIELDS) then
         write(dubuf,'(a,i9)')  ' INIT FINAL FLUSHING OF RESTARTING FIELDS n= ',this%n
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         call flush_and_save_resume(sgg, this%bounds, this%control%layoutnumber, this%control%size, this%control%nentradaroot, this%control%nresumeable2, this%thereare, this%n,eps0,mu0, this%everflushed,  &
         Ex, Ey, Ez, Hx, Hy, Hz,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic)
         write(dubuf,'(a,i9)')  ' DONE FINAL FLUSHING OF RESTARTING FIELDS N=',this%n
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         call print11(this%control%layoutnumber,dubuf)
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      endif

      if (this%thereAre%FarFields) then
          write(dubuf,'(a,i9)')  ' INIT FINAL OBSERVATION DATA FLUSHING and Near-to-Far field  n= ',this%n
      else
          write(dubuf,'(a,i9)')  ' INIT FINAL OBSERVATION DATA FLUSHING n= ',this%n
      endif
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      call print11(this%control%layoutnumber,dubuf)
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      if (this%thereAre%Observation) THEN
         call FlushObservationFiles(sgg,this%ini_save, this%n,this%control%layoutnumber, this%control%size, dxe, dye, dze, dxh, dyh, dzh,this%bounds,this%control%singlefilewrite,this%control%facesNF2FF,.TRUE.)
         call CloseObservationFiles(sgg,this%control%layoutnumber,this%control%size,this%control%singlefilewrite,this%initialtimestep,this%lastexecutedtime,this%control%resume) !dump the remaining to disk
#ifdef CompileWithMTLN      
         if (this%control%use_mtln_wires) then
            call FlushMTLNObservationFiles(this%control%nentradaroot, mtlnProblem = .false.)
         end if
#endif
      endif
      
      if (this%thereAre%FarFields) then
         write(dubuf,'(a,i9)')   ' DONE FINAL OBSERVATION DATA FLUSHED and Near-to-Far field  n= ',this%n
      else
         write(dubuf,'(a,i9)')    ' DONE FINAL OBSERVATION  DATA FLUSHED n= ',this%n
      endif
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      call print11(this%control%layoutnumber,dubuf)
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      write(dubuf,'(a,i9)') 'INIT FINAL Postprocessing frequency domain probes, if any, at n= ',this%n
      call print11(this%control%layoutnumber,dubuf)
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(this%control%layoutnumber,dubuf)
      somethingdone=.false.
      at=this%n*sgg%dt
      if (this%thereAre%Observation) call PostProcess(this%control%layoutnumber,this%control%size,sgg,this%control%nentradaroot,at,somethingdone,this%control%niapapostprocess,this%control%forceresampled)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
      call MPI_AllReduce(somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
      somethingdone=newsomethingdone
#endif

      if (somethingdone) then
         write(dubuf,*) 'DONE FINAL Postprocessing frequency domain probes.'
         call print11(this%control%layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%control%layoutnumber,dubuf)
      else
         write(dubuf,*) 'No FINAL frequency domain probes snapshots found to be postrocessed'
         call print11(this%control%layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%control%layoutnumber,dubuf)
      endif

      write(dubuf,*)'INIT FINAL FLUSHING .vtk if any.'
      call print11(this%control%layoutnumber,dubuf)
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(this%control%layoutnumber,dubuf)
      somethingdone=.false.

      if (this%thereAre%Observation) call createvtk(this%control%layoutnumber,this%control%size,sgg,this%control%vtkindex,somethingdone,this%control%mpidir,tagtype,this%media%sggMtag,this%control%dontwritevtk)

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
      call MPI_AllReduce(somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
      somethingdone=newsomethingdone
#endif
      if (somethingdone) then
        write(dubuf,*) 'DONE FINAL FLUSHING .vtk snapshots'
        call print11(this%control%layoutnumber,dubuf)
        write(dubuf,*) SEPARADOR//separador//separador
        call print11(this%control%layoutnumber,dubuf)
      else
        write(dubuf,*) 'No FINAL .vtk snapshots found to be flushed'
        call print11(this%control%layoutnumber,dubuf)
        write(dubuf,*) SEPARADOR//separador//separador
        call print11(this%control%layoutnumber,dubuf)
      endif

      write(dubuf,*)'INIT FINAL FLUSHING .xdmf if any.'
      call print11(this%control%layoutnumber,dubuf)
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(this%control%layoutnumber,dubuf)
      somethingdone=.false.
      if (this%thereAre%Observation) call createxdmf(sgg,this%control%layoutnumber,this%control%size,this%control%vtkindex,this%control%createh5bin,somethingdone,this%control%mpidir)
      if (this%control%createh5bin) call createh5bintxt(sgg,this%control%layoutnumber,this%control%size) !lo deben llamar todos haya o no this%thereAre%observation
!        call create_interpreted_mesh(sgg)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
      call MPI_AllReduce(somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
      somethingdone=newsomethingdone
#endif
      if (somethingdone) then
         write(dubuf,*) 'DONE FINAL FLUSHING .xdmf snapshots'
         call print11(this%control%layoutnumber,dubuf)  
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%control%layoutnumber,dubuf)
      else
         write(dubuf,*) 'No FINAL .xdmf snapshots found to be flushed'
         call print11(this%control%layoutnumber,dubuf)
         write(dubuf,*) SEPARADOR//separador//separador
         call print11(this%control%layoutnumber,dubuf)
      endif

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      call Timing(sgg,this%bounds,this%n,ndummy,this%control%layoutnumber, &
                  this%control%size, this%control%maxCPUtime,this%control%flushsecondsFields, &
                  this%control%flushsecondsData,this%initialtimestep, &
                  this%control%finaltimestep,this%perform,this%parar,.FALSE., &
                  Ex,Ey,Ez,this%everflushed,this%control%nentradaroot,this%control%maxSourceValue,this%control%opcionestotales, & 
                  this%control%simu_devia,this%control%dontwritevtk,this%control%permitscaling)
      write(dubuf,*)'END FINAL POSTPROCESSING at n= ',this%n
      call print11(this%control%layoutnumber,dubuf)
      finishedwithsuccess=.true.

      return

   end subroutine

   !las sggmixx se desctruyen el en main pq se alocatean alli
   subroutine Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,thereare,wiresflavor )
      character (len=*) , intent(in)    ::  wiresflavor
      type (Logic_control), intent(IN)  ::  thereare
      type (SGGFDTDINFO), intent(INOUT)     ::  sgg
      REAL (KIND=RKIND), intent(INOUT)     , pointer, dimension ( : , : , : )  ::  Ex,Ey,Ez,Hx,Hy,Hz
      REAL (KIND=RKIND), intent(INOUT)     , pointer, dimension ( : )  ::  G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh

      call DestroyObservation(sgg)
      Call DestroyNodal(sgg)
      call DestroyIlumina(sgg)
#ifdef CompileWithNIBC
      call DestroyMultiports(sgg)
#endif

      call destroysgbcs(sgg) !!todos deben destruir pq alocatean en funcion de sgg no de si contienen estos materiales que lo controla therearesgbcs. Lo que habia era IF ((this%thereAre%sgbcs).and.(sgbc))
      call destroyLumped(sgg)
      call DestroyEDispersives(sgg)
      call DestroyMDispersives(sgg)
      if ((trim(adjustl(wiresflavor))=='holland') .or. &
          (trim(adjustl(wiresflavor))=='transition')) then
         call DestroyWires(sgg)
      endif
#ifdef CompileWithBerengerWires
      if (trim(adjustl(wiresflavor))=='berenger') then
         call DestroyWires_Berenger(sgg)
      endif
#endif
#ifdef CompileWithSlantedWires
      if((trim(adjustl(wiresflavor))=='slanted').or.(trim(adjustl(wiresflavor))=='semistructured')) then
         call DestroyWires_Slanted(sgg)
      endif
#endif      

      call DestroyCPMLBorders
      call DestroyPMLbodies(sgg)
      call DestroyMURBorders
      !Destroy the remaining
      deallocate (sgg%Med,sgg%LineX,sgg%LineY,sgg%LineZ,sgg%DX,sgg%DY,sgg%DZ,sgg%tiempo)
      deallocate (G1,G2,GM1,GM2)
      deallocate (Ex, Ey, Ez, Hx, Hy, Hz)
      deallocate (dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh )
      return
   end subroutine Destroy_All_exceptSGGMxx

   subroutine destroy_and_deallocate(this,sgg)
      class(solver_t) :: this
      type (SGGFDTDINFO), intent(INOUT)     ::  sgg

      call DestroyObservation(sgg)
      Call DestroyNodal(sgg)
      call DestroyIlumina(sgg)
#ifdef CompileWithNIBC
      call DestroyMultiports(sgg)
#endif

      call destroysgbcs(sgg) !!todos deben destruir pq alocatean en funcion de sgg no de si contienen estos materiales que lo controla therearesgbcs. Lo que habia era IF ((this%thereAre%sgbcs).and.(sgbc))
      call destroyLumped(sgg)
      call DestroyEDispersives(sgg)
      call DestroyMDispersives(sgg)
      if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
          (trim(adjustl(this%control%wiresflavor))=='transition')) then
         call DestroyWires(sgg)
      endif
#ifdef CompileWithBerengerWires
      if (trim(adjustl(this%control%wiresflavor))=='berenger') then
         call DestroyWires_Berenger(sgg)
      endif
#endif
#ifdef CompileWithSlantedWires
      if((trim(adjustl(this%control%wiresflavor))=='slanted').or.(trim(adjustl(this%control%wiresflavor))=='semistructured')) then
         call DestroyWires_Slanted(sgg)
      endif
#endif      

      call DestroyCPMLBorders
      call DestroyPMLbodies(sgg)
      call DestroyMURBorders
      !Destroy the remaining
      deallocate (sgg%Med,sgg%LineX,sgg%LineY,sgg%LineZ,sgg%DX,sgg%DY,sgg%DZ,sgg%tiempo)
      call this%g%destroy()
      deallocate (this%Ex, this%Ey, this%Ez, this%Hx, this%Hy, this%Hz)
      deallocate (this%dxe, this%dye, this%dze, this%Idxe, this%Idye, this%Idze, this%dxh, this%dyh, this%dzh, this%Idxh, this%Idyh, this%Idzh)
      return
   end subroutine destroy_and_deallocate

   
   

end module
