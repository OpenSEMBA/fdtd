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
      REAL (kind=rkind), pointer, dimension ( : , : , : )  ::  Ex,Ey,Ez,Hx,Hy,Hz
   contains
      procedure :: init => solver_init
      procedure :: allocate_fields => solver_allocate_fields
      procedure :: launch_simulation
#ifdef CompileWithMTLN
      procedure :: launch_mtln_simulation
#endif
   end type

!    private

!    public launch_simulation
! #ifdef CompileWithMTLN
!    public launch_mtln_simulation
! #endif

   contains

   subroutine solver_init(this, input)
      class(solver_t) :: this
      type(entrada_t) :: input
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
      ! this%control%finishedwithsuccess = 
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

      ! time_desdelanzamiento = input%time_desdelanzamiento
      this%control%cfl = input%cfl
      this%control%attfactorc = input%attfactorc
      this%control%attfactorw = input%attfactorw
!       maxSourceValue = input%maxSourceValue
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

   subroutine solver_allocate_fields(this, sgg)
      class(solver_t) :: this
      type (sggfdtdinfo), intent(in)   ::  sgg
      allocate ( &
      this%Ex(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE),&
      this%Ey(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE),&
      this%Ez(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE),&
      this%Hx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE),&
      this%Hy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE),&
      this%Hz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE))
   end subroutine

#ifdef CompileWithMTLN
   subroutine launch_simulation(this, sgg,sggMtag,tag_numbers,sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, &
   SINPML_Fullsize,fullsize,finishedwithsuccess,Eps0,Mu0,tagtype,  &
   time_desdelanzamiento, maxSourceValue, EpsMuTimeScale_input_parameters, mtln_parsed)
#else
   subroutine launch_simulation(this, sgg,sggMtag,tag_numbers,sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, &
   SINPML_Fullsize,fullsize,finishedwithsuccess,Eps0,Mu0,tagtype,  &
   time_desdelanzamiento,  maxSourceValue, EpsMuTimeScale_input_parameters)
#endif
   !!!
      class(solver_t) :: this
#ifdef CompileWithMTLN
      type (mtln_t) :: mtln_parsed
#endif


      logical :: dummylog
      type (tagtype_t) :: tagtype

      !!for tuning
      !real (kind=rkind) :: time_elec=0.0_RKIND,time_magnet=0.0_RKIND
      !type (tiempo_t)  ::  time_MagnetInit,time_ElecInit,time_MagnetFin,time_ElecFin
      !!for tuning

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! SIMULATION VARIABLES
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      real (kind=rkind) :: maxSourceValue
      REAL (kind=8) :: time_desdelanzamiento
      type (EpsMuTimeScale_input_parameters_t) :: EpsMuTimeScale_input_parameters

      REAL (KIND=RKIND), intent(inout)              :: eps0,mu0
      real (kind=RKIND_tiempo) :: tiempoinicial,lastexecutedtime,ultimodt
      type (SGGFDTDINFO), intent(INOUT)   ::  sgg
      REAL (KIND=RKIND)                                   ::  dtcritico,newdtcritico
      REAL (KIND=RKIND)     , pointer, dimension ( : , : , : )  ::  Ex,Ey,Ez,Hx,Hy,Hz
      !!!! 
      integer (KIND=IKINDMTAG)   ::  &
      sggMtag(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      type(taglist_t) :: tag_numbers
      integer (KIND=INTEGERSIZEOFMEDIAMATRICES)   ::  &
      sggMiNo(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE), &
      sggMiEx(sgg%alloc(iEx)%XI : sgg%alloc(iEx)%XE,sgg%alloc(iEx)%YI : sgg%alloc(iEx)%YE,sgg%alloc(iEx)%ZI : sgg%alloc(iEx)%ZE), &
      sggMiEy(sgg%alloc(iEy)%XI : sgg%alloc(iEy)%XE,sgg%alloc(iEy)%YI : sgg%alloc(iEy)%YE,sgg%alloc(iEy)%ZI : sgg%alloc(iEy)%ZE), &
      sggMiEz(sgg%alloc(iEz)%XI : sgg%alloc(iEz)%XE,sgg%alloc(iEz)%YI : sgg%alloc(iEz)%YE,sgg%alloc(iEz)%ZI : sgg%alloc(iEz)%ZE), &
      sggMiHx(sgg%alloc(iHx)%XI : sgg%alloc(iHx)%XE,sgg%alloc(iHx)%YI : sgg%alloc(iHx)%YE,sgg%alloc(iHx)%ZI : sgg%alloc(iHx)%ZE), &
      sggMiHy(sgg%alloc(iHy)%XI : sgg%alloc(iHy)%XE,sgg%alloc(iHy)%YI : sgg%alloc(iHy)%YE,sgg%alloc(iHy)%ZI : sgg%alloc(iHy)%ZE), &
      sggMiHz(sgg%alloc(iHz)%XI : sgg%alloc(iHz)%XE,sgg%alloc(iHz)%YI : sgg%alloc(iHz)%YE,sgg%alloc(iHz)%ZI : sgg%alloc(iHz)%ZE)
      REAL (KIND=RKIND)     , pointer, dimension ( : )      ::  Idxe,Idye,Idze,Idxh,Idyh,Idzh,dxe,dye,dze,dxh,dyh,dzh
      !!!REAL (KIND=RKIND)     , pointer, dimension ( : )      ::  dxe_orig,dye_orig,dze_orig !deprecado 28/04/2014
      REAL (KIND=RKIND)     , pointer, dimension ( : )      ::  g1,g2,gM1,gM2
      !for lossy paddings
      REAL (KIND=RKIND)     :: Mur,epr,fmax,deltaespmax,Sigma,Epsilon,Mu,Sigmam,skin_depth,averagefactor,width,sigmatemp,eprtemp,murtemp,rdummy
      REAL (KIND=RKIND_tiempo)     :: at,rdummydt
      logical :: hayattmedia = .false.,attinformado = .false., somethingdone,newsomethingdone,call_timing,l_auxoutput,l_auxinput
      character(len=BUFSIZE) :: buff   
      integer (kind=4)      :: group_conformalprobes_dummy
      !
      !!!!!!!PML params!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical :: finishedwithsuccess
      !!!!!!!
      !Input
      type (bounds_t)  ::  b

      type (limit_t), dimension(1:6), intent(in)  ::  SINPML_fullsize,fullsize
      !
      character (LEN=BUFSIZE)     ::  chari,layoutcharID,dubuf
      integer (kind=4)   ::  ini_save,mindum
      !Generic
      type (Logic_control)  ::  thereare
      integer (kind=4) :: ierr,ndummy
      type (perform_t) :: perform, d_perform
      
      Logical  ::  parar,flushFF, &
                   everflushed,still_planewave_time,planewave_switched_off,thereareplanewave,l_aux

      integer (kind=4)  ::  i,J,K,r,n,initialtimestep,lastexecutedtimestep,n_info,FIELD,dummyMin,dummyMax
      !
      character (LEN=BUFSIZE)  ::  whoami
      !
      TYPE (tiempo_t) :: time_out2
       real (kind=RKIND) :: pscale_alpha
       integer :: rank
      !*******************************************************************************
      !*******************************************************************************
      !*******************************************************************************

      planewave_switched_off=.false.
      this%control%fatalerror=.false.
      parar=.false.
      call perform%reset()
      call d_perform%reset()
      flushFF=.false.
      everflushed=.false.
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
      !!!!!!!!!!!!!!!!!!!!!!!!END PRECHECKING

      write(whoami,'(a,i5,a,i5,a)') '(',this%control%layoutnumber+1,'/',this%control%size,') '

      !file names
      write(chari,*) this%control%layoutnumber+1
      !

      !!!!!!!write the material data in the Warnings file
      if ((this%control%layoutnumber == 0).and.this%control%verbose) call reportmedia(sgg)
      !
      layoutcharID = trim(adjustl(this%control%nentradaroot))//'_'//trim(adjustl(chari))

      !
      call findbounds(sgg,b)



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Space steps matrices creation (an extra cell is padded to deal with PMC imaging with no index errors
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !I need the whole increments to properly find the same time step in each layout

      allocate (dxe (sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      dye (sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      dze (sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      Idxe(sgg%ALLOC(iHx)%XI : sgg%ALLOC(iHx)%XE), &
      Idye(sgg%ALLOC(iHy)%YI : sgg%ALLOC(iHy)%YE), &
      Idze(sgg%ALLOC(iHz)%ZI : sgg%ALLOC(iHz)%ZE), &
      dxh (sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      dyh (sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      dzh (sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE), &
      Idxh(sgg%ALLOC(iEx)%XI : sgg%ALLOC(iEx)%XE), &
      Idyh(sgg%ALLOC(iEy)%YI : sgg%ALLOC(iEy)%YE), &
      Idzh(sgg%ALLOC(iEz)%ZI : sgg%ALLOC(iEz)%ZE))
      !
      dxe=-1.0e10_RKIND; dye=-1.0e10_RKIND; dze=-1.0e10_RKIND; dxh=-1.0e10_RKIND; dyh=-1.0e10_RKIND; dzh=-1.0e10_RKIND ;  !default values (ABSURD TO PREVEN ERRORS)


      !original
      do i=sgg%ALLOC(iHx)%XI,sgg%ALLOC(iHx)%XE
         dxe(i)=sgg%DX(i)
      end do

      do J=sgg%ALLOC(iHy)%YI,sgg%ALLOC(iHy)%YE
         dye(J)=sgg%DY(J)
      end do

      do K=sgg%ALLOC(iHz)%ZI,sgg%ALLOC(iHz)%ZE
         dze(K)=sgg%DZ(K)
      end do

      do i=sgg%ALLOC(iEx)%XI,sgg%ALLOC(iEx)%XE
         dxh(i)=(sgg%DX(i)+sgg%DX(i-1))/2.0_RKIND
      end do
      do J=sgg%ALLOC(iEy)%YI,sgg%ALLOC(iEy)%YE
         dyh(J)=(sgg%DY(J)+sgg%DY(J-1))/2.0_RKIND
      end do
      do K=sgg%ALLOC(iEz)%ZI,sgg%ALLOC(iEz)%ZE
         dzh(K)=(sgg%DZ(K)+sgg%DZ(K-1))/2.0_RKIND
      end do

      !!!lo he deprecado 28/04/2014 por incoherencia global con los deltas usados por todos lados
      !!!!mittra libro used in the stepping
      !!!dxe_orig=-1.0e10_RKIND; dye_orig=-1.0e10_RKIND; dze_orig=-1.0e10_RKIND;
      !!!do i=sgg%ALLOC(iHx)%XI,sgg%ALLOC(iHx)%XE
      !!!    dxe_mittra(i)=1.0_RKIND / 8.0_RKIND * (6.0_RKIND * sgg%DX(i)+sgg%DX(i-1)+sgg%DX(i+1))
      !!!end do
      !!!
      !!!do J=sgg%ALLOC(iHy)%YI,sgg%ALLOC(iHy)%YE
      !!!    dyee_mittra(J)=1.0_RKIND / 8.0_RKIND * (6.0_RKIND * sgg%DY(J)+sgg%DY(J-1)+sgg%DY(J+1))
      !!!end do
      !!!
      !!!do K=sgg%ALLOC(iHz)%ZI,sgg%ALLOC(iHz)%ZE
      !!!    dze_mittra(K)=1.0_RKIND / 8.0_RKIND * (6.0_RKIND * sgg%DZ(K)+sgg%DZ(K-1)+sgg%DZ(K+1))
      !!!end do
      !!!Idxe_mittra=1.0_RKIND/dxe_mittra ; Idye=1.0_RKIND/dye_mittra; Idze=1.0_RKIND/dze_mittra;
      !fin mitrra solo usado en time-stepping
!!!ojo que cpml toca los idxe, etc. para stretchaarlos con kappa (es 1 por lo general). Pero cuidado 251018
      Idxe=1.0_RKIND/dxe ; Idye=1.0_RKIND/dye; Idze=1.0_RKIND/dze; Idxh=1.0_RKIND/dxh; Idyh=1.0_RKIND/dyh; Idzh=1.0_RKIND/dzh;


!!!lo cambio aqui permit scaling a 211118 por problemas con resuming: debe leer el eps0, mu0, antes de hacer numeros

      allocate (G1(0 : sgg%NumMedia),G2(0 : sgg%NumMedia),GM1(0 : sgg%NumMedia),GM2(0 : sgg%NumMedia))
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Field matrices creation (an extra cell is padded at each limit and direction to deal with PMC imaging with no index errors)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !ojo las dimesniones deben ser giuales a las utlizadas en reallocate para las matrices sggmiEx, etc

      call this%allocate_fields(sgg)
      Ex => this%Ex
      Ey => this%Ey
      Ez => this%Ez
      Hx => this%Hx
      Hy => this%Hy
      Hz => this%Hz
      ! ALLOCATE ( &
      ! Ex(sgg%Alloc(iEx)%XI : sgg%Alloc(iEx)%XE,sgg%Alloc(iEx)%YI : sgg%Alloc(iEx)%YE,sgg%Alloc(iEx)%ZI : sgg%Alloc(iEx)%ZE),&
      ! Ey(sgg%Alloc(iEy)%XI : sgg%Alloc(iEy)%XE,sgg%Alloc(iEy)%YI : sgg%Alloc(iEy)%YE,sgg%Alloc(iEy)%ZI : sgg%Alloc(iEy)%ZE),&
      ! Ez(sgg%Alloc(iEz)%XI : sgg%Alloc(iEz)%XE,sgg%Alloc(iEz)%YI : sgg%Alloc(iEz)%YE,sgg%Alloc(iEz)%ZI : sgg%Alloc(iEz)%ZE),&
      ! Hx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE,sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE,sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE),&
      ! Hy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE,sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE,sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE),&
      ! Hz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE,sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE,sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE))

      !!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Init the local variables and observation stuff needed by each module, taking into account resume status
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      dt0=sgg%dt !guardalo aqui para entrada pscale correcta si resume
      if (.not.this%control%resume) then
         Ex=0.0_RKIND; Ey=0.0_RKIND; Ez=0.0_RKIND; Hx=0.0_RKIND; Hy=0.0_RKIND; Hz=0.0_RKIND
      !!!!
  !!!!!!!!!!!!!!!!?!?!?!?!?!?!?        Ex=1.0_RKIND; Ey=2.0_RKIND; Ez=3.0_RKIND; Hx=4.0_RKIND; Hy=5.0_RKIND; Hz=6.0_RKIND
         initialtimestep=0 !vamos a empezar en 0 para escribir el tiempo 0 !sgg sept'16 !?
         tiempoinicial = 0.0_RKIND_tiempo
         lastexecutedtimestep=0
         lastexecutedtime=0.0_RKIND_tiempo
      else
         write(dubuf,*) 'Init processing resuming data'
         call print11(this%control%layoutnumber,dubuf)
         !In case of resuming, the fields are read from disk
         if (this%control%resume_fromold) then
            open (14,file=trim(adjustl(this%control%nresumeable2))//'.old',form='unformatted')
         else
            open (14,file=trim(adjustl(this%control%nresumeable2)),form='unformatted')
         endif
         call ReadFields(sgg%alloc,lastexecutedtimestep,lastexecutedtime,ultimodt,eps0,mu0,Ex,Ey,Ez,Hx,Hy,Hz)
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
         call MPI_AllReduce( lastexecutedtimestep, dummyMin, 1_4, MPI_INTEGER, MPI_MIN, SUBCOMM_MPI, ierr)
         call MPI_AllReduce( lastexecutedtimestep, dummyMax, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)
         if ((dummyMax /= lastexecutedtimestep).or.(dummyMin /= lastexecutedtimestep)) then
#ifdef CompileWithOldSaving
            if (this%control%resume_fromold) then
               close (14)
               write(DUbuf,*) 'Incoherence between MPI saved steps for resuming.', dummyMin,dummyMax,lastexecutedtimesteP
               call stoponerror (this%control%layoutnumber,this%control%size,BUFF,.true.) !para que retorne
               call Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,this%thereare,this%control%wiresflavor )
               return
            else
               write(dubuf,*) 'Incoherence between MPI saved steps for resuming. Retrying with -old....'
               call print11(this%control%layoutnumber,dubuf)
               this%control%resume_fromold=.true.
               close (14)
               open (14,file=trim(adjustl(this%control%nresumeable2))//'.old',form='unformatted')
               call ReadFields(sgg%alloc,lastexecutedtimestep,lastexecutedtime,ultimodt,eps0,mu0,Ex,Ey,Ez,Hx,Hy,Hz)
               sgg%dt=ultimodt !para permit scaling
               call MPI_AllReduce( lastexecutedtimestep, dummyMin, 1_4, MPI_INTEGER, MPI_MIN, SUBCOMM_MPI, ierr)
               call MPI_AllReduce( lastexecutedtimestep, dummyMax, 1_4, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr)
               if ((dummyMax /= lastexecutedtimestep).or.(dummyMin /= lastexecutedtimestep)) then
                  write(DUbuf,*) 'NO success. fields.old MPI are also incoherent for resuming.', dummyMin,dummyMax,lastexecutedtimestep
                  call stoponerror (this%control%layoutnumber,this%control%size,DUBUF,.true.) !para que retorne
                  call Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,this%thereare,this%control%wiresflavor )
                  return
               else
                  write(dubuf,*) 'SUCCESS: Restarting from .fields.old instead. From n=',lastexecutedtimestep
                  call print11(this%control%layoutnumber,dubuf)
               endif
            endif
#else
            close (14)

            write(dubuf,*) 'Incoherence between MPI saved steps for resuming.',dummyMin,dummyMax,lastexecutedtimestep
            call stoponerror (this%control%layoutnumber,this%control%size,dubuf,.true.) !para que retorne
            call Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,this%thereare,this%control%wiresflavor )
            return
#endif
         endif
#endif
         initialtimestep=lastexecutedtimestep+1
         tiempoinicial = lastexecutedtime
         write(dubuf,*) '[OK] processing resuming data. Last executed time step ',lastexecutedtimestep
         call print11(this%control%layoutnumber,dubuf)
      endif
      if (initialtimestep>this%control%finaltimestep) then
          call stoponerror (this%control%layoutnumber,this%control%size,'Initial time step greater than final one',.true.) !para que retorne
          call Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,this%thereare,this%control%wiresflavor )
          return
      endif
!!!incializa el vector de tiempos para permit scaling 191118
      call crea_timevector(sgg,lastexecutedtimestep,this%control%finaltimestep,lastexecutedtime)
!!!!!!!!!!!!!!!!!!!!!

!fin lo cambio aqui

      call updateSigmaM()
      call updateThinWiresSigma()
      call calc_G1G2Gm1Gm2(sgg,G1,G2,Gm1,Gm2,eps0,mu0)
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
      call InitPMLbodies(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,g2,Gm2,this%thereAre%PMLbodies,this%control,eps0,mu0)
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

      !Initborders must be called in first place . Double check that Idxe,h are not needed by other initialization modules
      !llamalo antes de sgbc y composites, porque overrideo nodos sgbc conectados a hilos

      !init lumped debe ir antes de wires porque toca la conductividad del material !mmmm ojoooo 120123
      write(dubuf,*) 'Init Lumped Elements...';  call print11(this%control%layoutnumber,dubuf)
      CALL InitLumped(sgg,sggMiEx,sggMiEy,sggMiEz,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,this%control,this%thereAre%Lumpeds,eps0,mu0)
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
         
      ! one more MM for right adjancencies
      dtcritico=sgg%dt
      if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
          (trim(adjustl(this%control%wiresflavor))=='transition')) then
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Holland Wires...';  call print11(this%control%layoutnumber,dubuf)
         call InitWires       (sgg,sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz, & 
                               this%thereAre%Wires, Ex,Ey,Ez,Hx,Hy,Hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh, &
                               g2,SINPML_fullsize, fullsize,dtcritico,eps0,mu0,this%control)
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
         call InitWires_Berenger(sgg,sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,this%control%layoutnumber,this%control%size,this%thereAre%Wires,this%control%resume,this%control%makeholes, &
         this%control%isolategroupgroups,this%control%mtlnberenger,this%control%mindistwires, &
         this%control%groundwires,this%control%taparrabos,Ex,Ey,Ez, &
         Idxe,Idye,Idze,Idxh,Idyh,Idzh,this%control%inductance_model,g2,SINPML_fullsize,fullsize,dtcritico,eps0,mu0,this%control%verbose)
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
!!!does not work             
!!!!             precision=0
!!!!             call estructura_slanted(sgg,precision)   
             continue
         endif
         call InitWires_Slanted(sgg, this%control%layoutnumber,this%control%size, Ex, Ey, Ez,   & 
                                 Idxe, Idye, Idze, Idxh, Idyh, Idzh,   &
                                 sggMiNo,                              &
                                 sggMiEx, sggMiEy, sggMiEz,            &
                                 sggMiHx, sggMiHy, sggMiHz,            &
                                 this%thereAre%Wires, this%control%resume,               &
                                 this%control%mindistwires, this%control%groundwires,this%control%noSlantedcrecepelo ,     &
                                 this%control%inductance_model, this%control%inductance_order,   &
                                 g2, SINPML_fullsize, dtcritico,eps0,mu0,this%control%verbose)
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
         call InitWires_mtln(sgg,Ex,Ey,Ez,Idxh,Idyh,Idzh,eps0, mu0, mtln_parsed,this%thereAre%MTLNbundles)
#else
         write(buff,'(a)') 'WIR_ERROR: Executable was not compiled with MTLN modules.'
#endif
      endif


      !Anisotropic
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) 'Init Anisotropic...';  call print11(this%control%layoutnumber,dubuf)
      call InitAnisotropic(sgg,sggmiex,sggmiey,sggmiez,sggMiHx ,sggMiHy ,sggMiHz,this%thereAre%Anisotropic,this%thereAre%ThinSlot,eps0,mu0)
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

      IF (this%control%sgbc)  then
#ifdef CompileWithMPI
           call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
            write(dubuf,*) 'Init Multi sgbc...';  call print11(this%control%layoutnumber,dubuf)
            call Initsgbcs(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,Ex,Ey,Ez,Hx,Hy,Hz,IDxe,IDye,IDze,IDxh,IDyh,IDzh,this%control%layoutnumber,this%control%size, &
                 G1,G2,GM1,GM2,this%thereAre%sgbcs,this%control%resume,this%control%sgbccrank,this%control%sgbcFreq,this%control%sgbcresol,this%control%sgbcdepth,this%control%sgbcDispersive,eps0,mu0,this%control%simu_devia,this%control%stochastic)
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

!!!!
#ifdef CompileWithNIBC
      IF (this%control%mibc)  then
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Multiports...';  call print11(this%control%layoutnumber,dubuf)
         call InitMultiports        (sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx ,sggMiHy ,sggMiHz,this%control%layoutnumber,this%control%size,this%thereAre%Multiports,this%control%resume, &
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

      ![conformal] ##ref##
      !poner aqui mi inicializador de campos...
      !todo aquello que dependa sgg%dt
      !**************************************************************************************************
      !***[conformal]  *******************************************************************
      !**************************************************************************************************
      !--->[conformal](initiaizate memory EM fields)-----------------------------------------------------
      !ref: ##conf_timestepping_memory_ini##
#ifdef CompileWithConformal
      if(input_conformal_flag)then
#ifdef CompileWithMPI
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init Conformal Elements ...';  call print11(this%control%layoutnumber,dubuf)
!WIP
!DEBUG
         call initialize_memory_FDTD_conf_fields (sgg,sggMiEx, &
         & sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,Ex,Ey,Ez,Hx,Hy,Hz,&
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

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) 'Init EDispersives...';  call print11(this%control%layoutnumber,dubuf)
      call InitEDispersives(sgg,sggMiEx,sggMiEy,sggMiEz,this%thereAre%EDispersives,this%control%resume,g1,g2,ex,ey,ez)
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
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) 'Init MDispersives...';  call print11(this%control%layoutnumber,dubuf)
      call InitMDispersives(sgg,sggMiHx,sggMiHy,sggMiHz,this%thereAre%MDispersives,this%control%resume,gm1,gm2,hx,hy,hz)
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

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) 'Init Multi Plane-Waves...';  call print11(this%control%layoutnumber,dubuf)
      call InitPlaneWave   (sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,this%control%layoutnumber,this%control%size,SINPML_fullsize,this%thereAre%PlaneWaveBoxes,this%control%resume,eps0,mu0)

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

         !!!!!!!sgg 121020 !rellena la matriz Mtag con los slots de una celda
                 call fillMtag(sgg, sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz,sggMtag, b, tag_numbers)
         !!!!!!!fin
                 
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      write(dubuf,*) 'Init Observation...';  call print11(this%control%layoutnumber,dubuf)
      call InitObservation (sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag,tag_numbers, &
                            this%thereAre%Observation,this%thereAre%wires,this%thereAre%FarFields,this%control%resume,initialtimestep,this%control%finaltimestep,lastexecutedtime, &
                            this%control%nentradaroot,this%control%layoutnumber,this%control%size,this%control%saveall,this%control%singlefilewrite,this%control%wiresflavor,&
                            SINPML_FULLSIZE,this%control%facesNF2FF,this%control%NF2FFDecim,eps0,mu0,this%control%simu_devia,this%control%mpidir,this%control%niapapostprocess,b)
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
      !observation must be the last one to initialize
      
!!!!voy a jugar con fuego !!!210815 sincronizo las matrices de medios porque a veces se precisan. Reutilizo rutinas viejas mias NO CRAY. Solo se usan aqui
#ifdef CompileWithMPI
      !MPI initialization
      if (this%control%size>1) then
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) 'Init MPI MediaMatrix flush...';  call print11(this%control%layoutnumber,dubuf)
         call InitMPI(sgg%sweep,sgg%alloc)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         !write(dubuf,*) trim(adjustl(whoami))//' [OK]';  call print11(this%control%layoutnumber,dubuf,.true.)
         !write(dubuf,*) trim(adjustl(whoami))//' Init MPI Extra flushings...';  call print11(this%control%layoutnumber,dubuf,.true.)
         call InitExtraFlushMPI(this%control%layoutnumber,sgg%sweep,sgg%alloc,sgg%med,sgg%nummedia,sggmiEz,sggMiHz)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         !write(dubuf,*) trim(adjustl(whoami))//' [OK]';  call print11(this%control%layoutnumber,dubuf,.true.)
         !write(dubuf,*) trim(adjustl(whoami))//' First MPI H flushing...';  call print11(this%control%layoutnumber,dubuf,.true.)
         call FlushMPI_H(sgg%alloc,this%control%layoutnumber,this%control%size, sggmiHx,sggmiHy,sggmiHz)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         !write(dubuf,*) trim(adjustl(whoami))//' [OK]';  call print11(this%control%layoutnumber,dubuf,.true.)
         !write(dubuf,*) trim(adjustl(whoami))//' First MPI E flushing...';  call print11(this%control%layoutnumber,dubuf,.true.)
         call FlushMPI_E(sgg%alloc,this%control%layoutnumber,this%control%size, sggmiEx,sggmiEy,sggmiEz)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
         write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
      endif
#endif
!!!!!!!!!!!!!!!se supone que la inicializacion de Cray machacara luego a esta que solo uso para flushear medios (lo preciso en sgbcs de momento, pero es bueno tener esta info)
!!!!!!!!!!!!!!!!!!!!!fin juego con fuego 210815

#ifdef CompileWithMPI
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
         call InitExtraFlushMPI_Cray(this%control%layoutnumber,sgg%sweep,sgg%alloc,sgg%Med,sgg%NumMedia,sggMiez,sggMiHz, &
         Ex,Ey,Ez,Hx,Hy,Hz,this%thereAre%MURBorders)
         call MPI_Barrier(SUBCOMM_MPI,ierr)
         write(dubuf,*) '[OK]';  call print11(this%control%layoutnumber,dubuf)
      endif
#endif


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

 
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

      if (this%control%resume) close (14)
      !
      n=initialtimestep
      ini_save = initialtimestep
      n_info = 5 + initialtimestep
      !
!      if (verbose) call ReportExistence(sgg,this%control%layoutnumber,size, thereare,mur_second,MurAfterPML)


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! For Timing
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(dubuf,*) 'Init Timing...';  call print11(this%control%layoutnumber,dubuf)
      call InitTiming(sgg, this%control, time_desdelanzamiento, Initialtimestep,maxSourceValue)


      !!!if (createmap) then
      !!!    call writemmdxf(this%control%layoutnumber,sgg,sggMiHx,sggMiHy,sggMiHz)
      !!!endif
      !!!CALL CLOSEdxfFILE(this%control%layoutnumber,size)
      !!!!NO MORE WARNINGS SHOULD BE PRODUCED

      CALL CLOSEWARNINGFILE(this%control%layoutnumber,this%control%size,this%control%fatalerror,.false.,this%control%simu_devia) !aqui ya esta dividido el stochastic y hay dos this%control%layoutnumber=0

      if (this%control%fatalerror) then
         dubuf='FATAL ERRORS. Revise *Warnings.txt file. ABORTING...'
         call stoponerror(this%control%layoutnumber,this%control%size,dubuf,.true.) !para que retorne
         call Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,this%thereare,this%control%wiresflavor )
         return
      endif


#ifdef CompileWithMPI
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
#endif
!!!no se si el orden wires - sgbcs del sync importa 150519
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
          if (this%control%stochastic)  then
             call syncstoch_mpi_sgbcs(this%control%simu_devia,this%control%layoutnumber,this%control%size)
          endif
#endif    
#endif    

#ifdef CompileWithMPI
#ifdef CompileWithStochastic
          if (this%control%stochastic)  then
             call syncstoch_mpi_lumped(this%control%simu_devia,this%control%layoutnumber,this%control%size)
          endif
#endif    
#endif    
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF (this%control%resume) then
         write(dubuf,*)'END PREPROCESSING. RESUMING simulation from n=',n
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
      still_planewave_time=.true. !inicializacion de la variable 
     !!!aqui no. bug resume pscale 131020      ! dt0=sgg%dt !entrada pscale
      pscale_alpha=1.0 !se le entra con 1.0 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  TIME STEPPING
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                           
      
#ifdef CompileWithProfiling
      call nvtxStartRange("Antes del bucle N")
#endif
!240424 sgg creo el comunicador mpi de las sondas conformal aqui. debe irse con el nuevo conformal
#ifdef CompileWithConformal                
#ifdef CompileWithMPI
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
#endif  
#endif

      ciclo_temporal :  DO while (N <= this%control%finaltimestep)
      
         call step()

         IF (this%thereAre%Observation) then
            call UpdateObservation(sgg,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz,sggMtag,tag_numbers, n,ini_save, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh,this%control%wiresflavor,SINPML_FULLSIZE,this%control%wirecrank, this%control%noconformalmapvtk,b)
            if (n>=ini_save+BuffObse)  then
               mindum=min(this%control%finaltimestep,ini_save+BuffObse)
               call FlushObservationFiles(sgg,ini_save,mindum,this%control%layoutnumber,this%control%size, dxe, dye, dze, dxh, dyh, dzh,b,this%control%singlefilewrite,this%control%facesNF2FF,.FALSE.) !no se flushean los farfields ahora
            endif
         endif

         if(n >= n_info) then
             call_timing=.true.
         else
             call_timing=.false.
         endif
#ifdef CompileWithMPI
         l_aux=call_timing
         call MPI_AllReduce( l_aux, call_timing, 1_4, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
         call MPI_Barrier(MPI_COMM_WORLD,ierr) !050619 incluido problemas stochastic stopflusing
#endif
         
!!!
         if (call_timing) then
            call Timing(sgg,b,        n,n_info,this%control%layoutnumber,this%control%size, this%control%maxCPUtime,this%control%flushsecondsFields,this%control%flushsecondsData,initialtimestep, &
            this%control%finaltimestep,perform,parar,.FALSE., &
            Ex,Ey,Ez,everflushed,this%control%nentradaroot,maxSourceValue,this%control%opcionestotales,this%control%simu_devia,this%control%dontwritevtk,this%control%permitscaling)

            if (.not.parar) then !!! si es por parada se gestiona al final
!!!!! si esta hecho lo flushea todo pero poniendo de acuerdo a todos los mpi
                do i=1,sgg%NumberRequest
                   if  (sgg%Observation(i)%done.and.(.not.sgg%Observation(i)%flushed)) then
                      perform%flushXdmf=.true.
                      perform%flushVTK=.true.
                   endif
                end do
#ifdef CompileWithMPI
                l_aux=perform%flushVTK
                call MPI_AllReduce( l_aux, perform%flushVTK, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=perform%flushXdmf
                call MPI_AllReduce( l_aux, perform%flushXdmf, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=perform%flushDATA
                call MPI_AllReduce( l_aux, perform%flushDATA, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=perform%flushFIELDS
                call MPI_AllReduce( l_aux, perform%flushFIELDS, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                !
                l_aux=perform%postprocess
                call MPI_AllReduce( l_aux, perform%postprocess, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
#endif
!!!!!!!!!!!!
                if (perform%flushFIELDS) then
                   write(dubuf,*)  SEPARADOR,trim(adjustl(this%control%nentradaroot)),separador
                   call print11(this%control%layoutnumber,dubuf)
                   write(dubuf,*)  'INIT FLUSHING OF RESTARTING FIELDS n=',N
                   call print11(this%control%layoutnumber,dubuf)
                   call flush_and_save_resume(sgg, b, this%control%layoutnumber, this%control%size, this%control%nentradaroot, this%control%nresumeable2, this%thereare, n,eps0,mu0, everflushed,  &
                   Ex, Ey, Ez, Hx, Hy, Hz,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic)
#ifdef CompileWithMPI
                   call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
                   write(dubuf,*) SEPARADOR//separador//separador
                   call print11(this%control%layoutnumber,dubuf)
                   write(dubuf,*) 'DONE FLUSHING OF RESTARTING FIELDS n=',N
                   call print11(this%control%layoutnumber,dubuf)
                   write(dubuf,*) SEPARADOR//separador//separador
                   call print11(this%control%layoutnumber,dubuf)
                endif
                if (perform%isFlush()) then
                      !
                      flushFF=perform%postprocess
                      if (this%thereAre%FarFields.and.flushFF) then
                          write(dubuf,'(a,i9)')  ' INIT OBSERVATION DATA FLUSHING and Near-to-Far field n= ',n
                      else
                          write(dubuf,'(a,i9)')  ' INIT OBSERVATION DATA FLUSHING n= ',n
                      endif
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                      call print11(this%control%layoutnumber,dubuf)
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
    !!
                      if (this%thereAre%Observation) call FlushObservationFiles(sgg,ini_save, n,this%control%layoutnumber, this%control%size, dxe, dye, dze, dxh, dyh, dzh,b,this%control%singlefilewrite,this%control%facesNF2FF,flushFF)
                      !!
#ifdef CompileWithMPI
                      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
                      if (this%thereAre%FarFields.and.flushFF) then
                          write(dubuf,'(a,i9)')  ' Done OBSERVATION DATA FLUSHED and Near-to-Far field n= ',n
                      else
                          write(dubuf,'(a,i9)')  ' Done OBSERVATION DATA FLUSHED n= ',n
                      endif
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                      call print11(this%control%layoutnumber,dubuf)
                      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
    !
                      if (perform%postprocess) then
                         write(dubuf,'(a,i9)') 'Postprocessing frequency domain probes, if any, at n= ',n
                         call print11(this%control%layoutnumber,dubuf)
                         write(dubuf,*) SEPARADOR//separador//separador
                         call print11(this%control%layoutnumber,dubuf)
                         somethingdone=.false.
                         at=n*sgg%dt
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
                      if (perform%flushvtk) then   
                         write(dubuf,'(a,i9)')  ' Post-processing .vtk files n= ',n
                         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                         call print11(this%control%layoutnumber,dubuf)
                         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                         somethingdone=.false.
                         if (this%thereAre%Observation) call createvtkOnTheFly(this%control%layoutnumber,this%control%size,sgg,this%control%vtkindex,somethingdone,this%control%mpidir,tagtype,sggMtag,this%control%dontwritevtk)
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
                         if (perform%flushXdmf) then
                            write(dubuf,'(a,i9)')  ' Post-processing .xdmf files n= ',n
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
                 endif !del if (performflushDATA.or....
    !
                 if (this%control%singlefilewrite.and.perform%Unpack) then
                       call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                       write(dubuf,'(a,i9)')  ' Unpacking .bin files and prostprocessing them at n= ',n
                       call print11(this%control%layoutnumber,dubuf)
                       call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                       if (this%thereAre%Observation) call unpacksinglefiles(sgg,this%control%layoutnumber,this%control%size,this%control%singlefilewrite,initialtimestep,this%control%resume) !dump the remaining to disk
                       somethingdone=.false.
                       if (this%control%singlefilewrite.and.perform%Unpack) then
                           at=n*sgg%dt
                           if (this%thereAre%Observation) call PostProcessOnthefly(this%control%layoutnumber,this%control%size,sgg,this%control%nentradaroot,at,somethingdone,this%control%niapapostprocess,this%control%forceresampled)
                       endif
#ifdef CompileWithMPI
                       call MPI_Barrier(SUBCOMM_MPI,ierr)
                       call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
                       somethingdone=newsomethingdone
#endif
                       write(dubuf,'(a,i9)')  ' Done Unpacking .bin files and prostprocessing them at n= ',n
                       call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                       call print11(this%control%layoutnumber,dubuf)
                       call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                    endif !del if (singlefilewrite....
!!!si ha hecho algo reporta que va a continuar          
                    if ((this%control%singlefilewrite.and.perform%Unpack).or.perform%isFlush()) then
                        write(dubuf,'(a,i9)')  ' Continuing simulation at n= ',n
                        call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                        call print11(this%control%layoutnumber,dubuf)
                        call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
                    endif
                endif !!!del if (.not.parar)
             endif !!!del if(n >= n_info
         !!!!!!!!all the previous must be together
              
         this%control%fatalerror=.false.
         if (parar) then
             this%control%fatalerror=.true.
             exit ciclo_temporal
         endif
#ifdef CompileWithPrescale
         if (this%control%permitscaling) then
#ifndef miguelPscaleStandAlone
            if ((sgg%tiempo(n)>=EpsMuTimeScale_input_parameters%tini).and.&
                &(sgg%tiempo(n)<=EpsMuTimeScale_input_parameters%tend)) then
#endif
             call updateconstants(sgg,n,this%thereare,g1,g2,gM1,gM2, & 
                               Idxe,Idye,Idze,Idxh,Idyh,Idzh, &  !needed by  CPML to be updated
                               this%control%sgbc,this%control%mibc,input_conformal_flag, &
                               this%control%wiresflavor, this%control%wirecrank, this%control%fieldtotl,&
                               this%control%sgbcDispersive,this%control%finaltimestep, &
                               eps0,mu0, &
                               this%control%simu_devia, &
                               EpsMuTimeScale_input_parameters,pscale_alpha,still_planewave_time &
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
         n=n+1 !sube de iteracion
      end do ciclo_temporal ! End of the time-stepping loop
      
                        
      
#ifdef CompileWithProfiling
      call nvtxEndRange
#endif      
      
#ifdef CompileWithConformal
      if(input_conformal_flag)then
            call conformal_final_simulation  (conf_timeSteps, n)
      endif
#endif

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (n>this%control%finaltimestep) n=this%control%finaltimestep !readjust n since after finishing it is increased
      this%control%finaltimestep=n
      lastexecutedtime=sgg%tiempo(this%control%finaltimestep)
      !se llama con dummylog para no perder los flags de parada
      call Timing(sgg,b,        n,ndummy,this%control%layoutnumber,this%control%size, this%control%maxCPUtime,this%control%flushsecondsFields,this%control%flushsecondsData,initialtimestep, &
            this%control%finaltimestep,d_perform,dummylog,.FALSE., &
            Ex,Ey,Ez,everflushed,this%control%nentradaroot,maxSourceValue,this%control%opcionestotales,this%control%simu_devia,this%control%dontwritevtk,this%control%permitscaling)

      write(dubuf,*)'END FDTD time stepping. Beginning posprocessing at n= ',n
      call print11(this%control%layoutnumber,dubuf)

      if ((this%control%flushsecondsFields/=0).or.perform%flushFIELDS) then
         write(dubuf,'(a,i9)')  ' INIT FINAL FLUSHING OF RESTARTING FIELDS n= ',n
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         call flush_and_save_resume(sgg, b, this%control%layoutnumber, this%control%size, this%control%nentradaroot, this%control%nresumeable2, this%thereare, n,eps0,mu0, everflushed,  &
         Ex, Ey, Ez, Hx, Hy, Hz,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic)
         write(dubuf,'(a,i9)')  ' DONE FINAL FLUSHING OF RESTARTING FIELDS N=',n
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
         call print11(this%control%layoutnumber,dubuf)
         call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      endif
!
      if (this%thereAre%FarFields) then
          write(dubuf,'(a,i9)')  ' INIT FINAL OBSERVATION DATA FLUSHING and Near-to-Far field  n= ',n
      else
          write(dubuf,'(a,i9)')  ' INIT FINAL OBSERVATION DATA FLUSHING n= ',n
      endif
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      call print11(this%control%layoutnumber,dubuf)
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      if (this%thereAre%Observation) THEN
         !dump the remaining to disk
         call FlushObservationFiles(sgg,ini_save, n,this%control%layoutnumber, this%control%size, dxe, dye, dze, dxh, dyh, dzh,b,this%control%singlefilewrite,this%control%facesNF2FF,.TRUE.)
         call CloseObservationFiles(sgg,this%control%layoutnumber,this%control%size,this%control%singlefilewrite,initialtimestep,lastexecutedtime,this%control%resume) !dump the remaining to disk
#ifdef CompileWithMTLN      
         if (this%control%use_mtln_wires) then
            ! call GatherMPI_MTL()
            call FlushMTLNObservationFiles(this%control%nentradaroot, mtlnProblem = .false.)
         end if
#endif
      endif
      
      if (this%thereAre%FarFields) then
          write(dubuf,'(a,i9)')   ' DONE FINAL OBSERVATION DATA FLUSHED and Near-to-Far field  n= ',n
      else
         write(dubuf,'(a,i9)')    ' DONE FINAL OBSERVATION  DATA FLUSHED n= ',n
      endif
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)
      call print11(this%control%layoutnumber,dubuf)
      call print11(this%control%layoutnumber,SEPARADOR//separador//separador)

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif

       write(dubuf,'(a,i9)') 'INIT FINAL Postprocessing frequency domain probes, if any, at n= ',n
       call print11(this%control%layoutnumber,dubuf)
       write(dubuf,*) SEPARADOR//separador//separador
       call print11(this%control%layoutnumber,dubuf)
       somethingdone=.false.
        at=n*sgg%dt
       if (this%thereAre%Observation) call PostProcess(this%control%layoutnumber,this%control%size,sgg,this%control%nentradaroot,at,somethingdone,this%control%niapapostprocess,this%control%forceresampled)
#ifdef CompileWithMPI
       call MPI_Barrier(SUBCOMM_MPI,ierr)
       call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
       somethingdone=newsomethingdone
#endif
      !!!!!!!!!!
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
!
      write(dubuf,*)'INIT FINAL FLUSHING .vtk if any.'
      call print11(this%control%layoutnumber,dubuf)
      write(dubuf,*) SEPARADOR//separador//separador
      call print11(this%control%layoutnumber,dubuf)
      somethingdone=.false.

      if (this%thereAre%Observation) call createvtk(this%control%layoutnumber,this%control%size,sgg,this%control%vtkindex,somethingdone,this%control%mpidir,tagtype,sggMtag,this%control%dontwritevtk)

#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
      call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
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
!
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
      call MPI_AllReduce( somethingdone, newsomethingdone, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr)
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
      call Timing(sgg,b,        n,ndummy,this%control%layoutnumber,this%control%size, this%control%maxCPUtime,this%control%flushsecondsFields,this%control%flushsecondsData,initialtimestep, &
            this%control%finaltimestep,perform,parar,.FALSE., &
            Ex,Ey,Ez,everflushed,this%control%nentradaroot,maxSourceValue,this%control%opcionestotales,this%control%simu_devia,this%control%dontwritevtk,this%control%permitscaling)
      write(dubuf,*)'END FINAL POSTPROCESSING at n= ',n
      call print11(this%control%layoutnumber,dubuf)
      finishedwithsuccess=.true.
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
!!!me he dado cuenta de que nunca entra aqui hoy 120617 pero no me he atrevido a borrar las lineas que siguen
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!  Tell each module to free-up memory
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call Destroy_All_exceptSGGMxx(sgg,Ex, Ey, Ez, Hx, Hy, Hz,G1,G2,GM1,GM2,dxe  ,dye  ,dze  ,Idxe ,Idye ,Idze ,dxh  ,dyh  ,dzh  ,Idxh ,Idyh ,Idzh,this%thereare,this%control%wiresflavor )
      !
#ifdef CompileWithMPI
      call MPI_Barrier(SUBCOMM_MPI,ierr)
#endif
      !---------------------------------------------------->

   contains

      subroutine updateSigmaM()
         if (abs(this%control%attfactorc-1.0_RKIND) > 1.0e-12_RKIND) then
            hayattmedia=.false.
            attinformado=.false.
            do i=1,sgg%nummedia
               if (sgg%Med(i)%Is%MultiportPadding) then
                  sgg%Med(i)%SigmaM =(-2.0_RKIND * (-1.0_RKIND + this%control%attfactorc)*mu0)/((1 + this%control%attfactorc)*sgg%dt)
                  hayattmedia=.true.
               endif
               deltaespmax=max(max(maxval(sgg%dx),maxval(sgg%dy)),maxval(sgg%dz))
               if (hayattmedia.and. .not. attinformado) then
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
                  attinformado=.true.
               endif
            end do
         endif
      end subroutine updateSigmaM

      subroutine updateThinWiresSigma
         !thin wires !
         if (abs(this%control%attfactorw-1.0_RKIND) > 1.0e-12_RKIND) then
            attinformado=.false.
            do i=1,sgg%nummedia
               if (sgg%Med(i)%Is%ThinWire) then
                  sgg%Med(i)%Sigma =(-2.0_RKIND * (-1.0_RKIND + this%control%attfactorw)*eps0)/((1 + this%control%attfactorw)*sgg%dt)
                  if (.not.attinformado) then
                     write(buff,'(a,2e10.2e3)') ' WIREs stabilization att. factors=',this%control%attfactorw,sgg%Med(i)%Sigma
                     if (this%control%layoutnumber == 0) call WarnErrReport(buff)
                     attinformado=.true.
                  endif
               endif
            end do
         endif
      end subroutine updateThinWiresSigma

      subroutine revertThinWiresSigma()
         if (abs(this%control%attfactorw-1.0_RKIND) > 1.0e-12_RKIND) then
            do i=1,sgg%nummedia
               !thin wires
               if (sgg%Med(i)%Is%ThinWire) then
                  sgg%Med(i)%Sigma = 0.0_RKIND !revert!!! !necesario para no lo tome como un lossy luego en wires !solo se toca el g1,g2
               endif
            end do
         endif
      end subroutine

      subroutine reportSimulationOptions()
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


      subroutine flushPlanewaveOff(pw_switched_off, pw_still_time, pw_thereAre)
         logical, intent(inout) :: pw_switched_off, pw_still_time, pw_thereAre
         logical :: pw_still_time_aux, pw_thereAre_aux
         integer (kind=4) :: ierr
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
               write(dubuf,*) 'Switching plane-wave off at n=', N
               if (pw_thereAre) call print11(this%control%layoutnumber,dubuf)
            endif
         endif
      end subroutine 

      subroutine step()
         call flushPlanewaveOff(planewave_switched_off, still_planewave_time, thereareplanewave)
         IF (this%thereAre%Anisotropic) call AdvanceAnisotropicE(sgg%alloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)
         call advanceE()
#ifdef CompileWithConformal
         call advanceConformalE()
#endif
         call advanceWires()
         call advancePMLE()

#ifdef CompileWithNIBC
         IF (this%thereAre%Multiports.and.(this%control%mibc))      call AdvanceMultiportE(sgg%alloc,Ex, Ey, Ez)
#endif
         IF (this%thereAre%sgbcs.and.(this%control%sgbc))  then
            call AdvancesgbcE(real(sgg%dt,RKIND),this%control%sgbcDispersive,this%control%simu_devia,this%control%stochastic)
         endif
         if (this%thereAre%Lumpeds) call AdvanceLumpedE(sgg,n,this%control%simu_devia,this%control%stochastic)
         IF (this%thereAre%Edispersives)     call AdvanceEDispersiveE(sgg)
         If (this%thereAre%PlaneWaveBoxes.and.still_planewave_time) then
              if (.not.this%control%simu_devia) then
                 call AdvancePlaneWaveE(sgg,n, b       ,G2,Idxh,Idyh,Idzh,Ex,Ey,Ez,still_planewave_time)
              endif
          endif
         If (this%thereAre%NodalE) then
            call AdvanceNodalE(sgg,sggMiEx,sggMiEy,sggMiEz,sgg%NumMedia,n, b,G2,Idxh,Idyh,Idzh,Ex,Ey,Ez,this%control%simu_devia)
         endif

#ifdef CompileWithMPI
         if (this%control%size>1) then
            call MPI_Barrier(SUBCOMM_MPI,ierr)
            call   FlushMPI_E_Cray
         endif
#endif
         IF (this%thereAre%Anisotropic) call AdvanceAnisotropicH(sgg%alloc,ex,ey,ez,hx,hy,hz,Idxe,Idye,Idze,Idxh,Idyh,Idzh)
#ifdef CompileWithProfiling    
         call nvtxStartRange("Antes del bucle HX")
#endif
         call Advance_Hx           (Hx, Ey, Ez, Idye, Idze, sggMiHx, b,gm1,gm2)        
#ifdef CompileWithProfiling    
         call nvtxEndRange
         call nvtxStartRange("Antes del bucle HY")
#endif
         call Advance_Hy           (Hy, Ez, Ex, Idze, Idxe, sggMiHy, b,gm1,gm2)     
#ifdef CompileWithProfiling    
         call nvtxEndRange
         call nvtxStartRange("Antes del bucle HZ")
#endif
         call Advance_Hz           (Hz, Ex, Ey, Idxe, Idye, sggMiHz, b,gm1,gm2)  
#ifdef CompileWithProfiling    
         call nvtxEndRange
#endif
         
         If (this%thereAre%PMLbodies) then !waveport absorbers
            call AdvancePMLbodyH
         endif
         If (this%thereAre%PMLBorders) then
               call AdvanceMagneticCPML          ( sgg%NumMedia, b, sggMiHx, sggMiHy, sggMiHz, gm2, Hx, Hy, Hz, Ex, Ey, Ez)
         endif
         
         If (this%thereAre%PMCBorders)     then
            call MinusCloneMagneticPMC(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,this%control%layoutnumber,this%control%size)
         endif
         If (this%thereAre%PeriodicBorders)     then
            call CloneMagneticPeriodic(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,this%control%layoutnumber,this%control%size)
         endif
         IF (this%thereAre%sgbcs.and.(this%control%sgbc))  then
            call AdvancesgbcH()
         endif
         IF (this%thereAre%Mdispersives)     call AdvanceMDispersiveH(sgg)
#ifdef CompileWithNIBC
         IF (this%thereAre%Multiports    .and.(this%control%mibc))  &
         call AdvanceMultiportH    (sgg%alloc,Hx,Hy,Hz,Ex,Ey,Ez,Idxe,Idye,Idze,sggMiHx,sggMiHy,sggMiHz,gm2,sgg%nummedia,this%control%conformalskin)
#endif
         If (this%thereAre%PlaneWaveBoxes.and.still_planewave_time)  then
              if (.not.this%control%simu_devia) then
                 call AdvancePlaneWaveH(sgg,n, b        , GM2, Idxe,Idye, Idze, Hx, Hy, Hz,still_planewave_time)
              endif
         endif
         If (this%thereAre%NodalH) then
                 call AdvanceNodalH(sgg,sggMiHx,sggMiHy,sggMiHz,sgg%NumMedia,n, b       ,GM2,Idxe,Idye,Idze,Hx,Hy,Hz,this%control%simu_devia)
         endif

         if ((trim(adjustl(this%control%wiresflavor))=='holland') .or. &
             (trim(adjustl(this%control%wiresflavor))=='transition')) then
            IF (this%thereAre%Wires) then
               if (this%control%wirecrank) then
                  continue
               else
                  call AdvanceWiresH(sgg,n, this%control%layoutnumber,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic,this%control%experimentalVideal,this%control%wirethickness,eps0,mu0)
               endif
            endif
         endif
         If (this%thereAre%PMCBorders)     call MinusCloneMagneticPMC(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,this%control%layoutnumber,this%control%size)
         If (this%thereAre%PeriodicBorders)     then
            call CloneMagneticPeriodic(sgg%alloc,sgg%Border,Hx,Hy,Hz,sgg%sweep,this%control%layoutnumber,this%control%size)
         endif

#ifdef CompileWithConformal                      
         if(input_conformal_flag)then     
            call conformal_advance_H()
         endif
#endif

#ifdef CompileWithMPI
         !!Flush all the MPI (esto estaba justo al principo del bucle temporal diciendo que era necesario para correcto resuming)
         !lo he movido aqui a 16/10/2012 porque el farfield necesita tener los campos magneticos correctos
         !e intuyo que el Bloque current tambien a tenor del comentario siguiente
         !Incluyo un flush inicial antes de entrar al bucle para que el resuming sea correcto
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
#endif

!!!no se si el orden wires - sgbcs del sync importa 150519
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
          if (this%control%stochastic)  then
             call syncstoch_mpi_sgbcs(this%control%simu_devia,this%control%layoutnumber,this%control%size)
          endif
#endif    
#endif

#ifdef CompileWithMPI
#ifdef CompileWithStochastic
          if (this%control%stochastic)  then
             call syncstoch_mpi_lumped(this%control%simu_devia,this%control%layoutnumber,this%control%size)
          endif
#endif    
#endif 
         If (this%thereAre%MURBorders) then
            call AdvanceMagneticMUR              (b, sgg,sggMiHx, sggMiHy, sggMiHz, Hx, Hy, Hz,this%control%mur_second)
#ifdef CompileWithMPI
            if (this%control%mur_second) then
               if (this%control%size>1) then
                  call MPI_Barrier(SUBCOMM_MPI,ierr)
                  call   FlushMPI_H_Cray
               endif
            endif
#endif
         ENDIF


      end subroutine



      subroutine advanceE()
#ifdef CompileWithProfiling
         call nvtxStartRange("Antes del bucle EX")
#endif
         call Advance_Ex          (Ex, Hy, Hz, Idyh, Idzh, sggMiEx, b,g1,g2)    
#ifdef CompileWithProfiling
         call nvtxEndRange

         call nvtxStartRange("Antes del bucle EY")
#endif
         call Advance_Ey          (Ey, Hz, Hx, Idzh, Idxh, sggMiEy, b,g1,g2)
         
#ifdef CompileWithProfiling    
         call nvtxEndRange

         call nvtxStartRange("Antes del bucle EZ")
#endif
         call Advance_Ez          (Ez, Hx, Hy, Idxh, Idyh, sggMiEz, b,g1,g2)
#ifdef CompileWithProfiling    
         call nvtxEndRange
#endif
      end subroutine

      subroutine Advance_Ex(Ex,Hy,Hz,Idyh,Idzh,sggMiEx,b,g1,g2)

         !------------------------>
         type (bounds_t), intent( IN)  ::  b
         REAL (KIND=RKIND)     , pointer, dimension ( : )  ::  g1, g2
         !
         real (kind = RKIND), dimension    ( 0 :     b%dyh%NY-1     )  , intent( IN)  ::  Idyh
         real (kind = RKIND), dimension    ( 0 :     b%dzh%NZ-1     )  , intent( IN)  ::  Idzh
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN)     ::  sggMiEx
         real (kind = RKIND), dimension    ( 0 :      b%Ex%NX-1 , 0 :      b%Ex%NY-1 , 0 :      b%Ex%NZ-1 )  , intent( INOUT)  ::  Ex
         real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 )  , intent( IN)  ::  HY
         real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 )  , intent( IN)  ::  HZ
         !------------------------> Variables locales
         real (kind = RKIND)  ::  Idzhk, Idyhj
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzhk,Idyhj) 
#endif
#ifdef CompileWithACC   
!$ACC parallel loop DEFAULT(present) collapse (2) private (i,j,k,medio,Idzhk,Idyhj)  copyin(Ex,sggMiEx,Hy,Hz,Idyh,Idzh,b,G1,G2) copyout(Ex) 
#endif
         Do k=1,b%sweepEx%NZ
            Do j=1,b%sweepEx%NY
               Do i=1,b%sweepEx%NX
                  Idzhk=Idzh(k)
                  Idyhj=Idyh(j)
                  medio =sggMiEx(i,j,k)
                  Ex(i,j,k)=G1(MEDIO)*Ex(i,j,k)+G2(MEDIO)* &
                  ((Hz(i,j,k)-Hz(i,j-1,k))*Idyhj-(Hy(i,j,k)-Hy(i,j,k-1))*Idzhk)
               End do
            End do
         End do
#ifdef CompileWithOpenMP   
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine Advance_Ex
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Advance_Ey(Ey,Hz,Hx,Idzh,Idxh,sggMiEy,b,g1,g2)

         !------------------------>
         type (bounds_t), intent( IN)  ::  b
         REAL (KIND=RKIND)     , pointer, dimension ( : )   ::  g1, g2
         !
         real (kind = RKIND), dimension    ( 0 :     b%dzh%NZ-1     )  , intent( IN)  ::  Idzh
         real (kind = RKIND), dimension    ( 0 :     b%dxh%NX-1     )  , intent( IN)  ::  Idxh
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1   )  , intent( IN)     ::  sggMiEy
         real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 )  , intent( INOUT)  ::  EY
         real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 )  , intent( IN)  ::  HZ
         real (kind = RKIND), dimension    ( 0 :      b%Hx%NX-1 , 0 :      b%Hx%NY-1 , 0 :      b%Hx%NZ-1 )  , intent( IN)  ::  HX
         !------------------------> Variables locales
         real (kind = RKIND)  ::  Idzhk
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzhk)  
#endif
#ifdef CompileWithACC   
!$ACC parallel loop  DEFAULT(present) collapse (2) private (i,j,k,medio,Idzhk)     copyin(Ey,sggMiEy,Hz,Hx,Idzh,Idxh,b,G1,G2) copyout(Ey) 
#endif
         Do k=1,b%sweepEy%NZ
            Do j=1,b%sweepEy%NY
               Do i=1,b%sweepEy%NX
                  Idzhk=Idzh(k)
                  medio =sggMiEy(i,j,k)
                  Ey(i,j,k)=G1(MEDIO)*Ey(i,j,k)+G2(MEDIO)*((Hx(i,j,k)-Hx(i,j,k-1))*Idzhk-(Hz(i,j,k)-Hz(i-1,j,k))*Idxh(i))
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif



         return
      end subroutine Advance_Ey

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Advance_Ez(Ez,Hx,Hy,Idxh,Idyh,sggMiEz,b,g1,g2)

         !------------------------>
         type (bounds_t), intent( IN)  ::  b
         REAL (KIND=RKIND)     , pointer, dimension ( : )   ::  g1, g2
         !
         real (kind = RKIND), dimension    ( 0 :     b%dyh%NY-1     )  , intent( IN)  ::  Idyh
         real (kind = RKIND), dimension    ( 0 :     b%dxh%NX-1     )  , intent( IN)  ::  Idxh
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN)     ::  sggMiEz
         real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 )  , intent( INOUT)  ::  Ez
         real (kind = RKIND), dimension    ( 0 :      b%HX%NX-1 , 0 :      b%HX%NY-1 , 0 :      b%HX%NZ-1 )  , intent( IN)  ::  HX
         real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 )  , intent( IN)  ::  HY
         !------------------------> Variables locales
         real (kind = RKIND)  ::   Idyhj
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idyhj)    
#endif
#ifdef CompileWithACC   
!$ACC parallel loop   DEFAULT(present) collapse (2) private (i,j,k,medio,Idyhj)        copyin(Ez,sggMiEz,Hx,Hy,Idxh,Idyh,b,G1,G2) copyout(Ez) 
#endif
         Do k=1,b%sweepEz%NZ
            Do j=1,b%sweepEz%NY
               Do i=1,b%sweepEz%NX
                  Idyhj=Idyh(j)
                  medio =sggMiEz(i,j,k)
                  Ez(i,j,k)=G1(MEDIO)*Ez(i,j,k)+G2(MEDIO)*((Hy(i,j,k)-Hy(i-1,j,k))*Idxh(i)-(Hx(i,j,k)-Hx(i,j-1,k))*Idyhj)
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine Advance_Ez

      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Advance_Hx(Hx,Ey,Ez,IdyE,IdzE,sggMiHx,b,gm1,gm2)

         !------------------------>
         type (bounds_t), intent( IN)  ::  b
         REAL (KIND=RKIND)     , pointer, dimension ( : )   ::  gm1 ,gm2
         !!
         real (kind = RKIND), dimension    ( 0 :     b%dyE%NY-1    )  , intent( IN)  ::  IdyE
         real (kind = RKIND), dimension    ( 0 :     b%dzE%NZ-1    )  , intent( IN)  ::  IdzE
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN)     ::  sggMiHx
         real (kind = RKIND), dimension    ( 0 :      b%Hx%NX-1 , 0 :      b%Hx%NY-1 , 0 :      b%Hx%NZ-1 )  , intent( INOUT)  ::  Hx
         real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 )  , intent( IN)  ::  EY
         real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 )  , intent( IN)  ::  EZ
         !------------------------> Variables locales
         real (kind = RKIND)  ::  Idzek, Idyej
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzek,Idyej)     
#endif
#ifdef CompileWithACC   
!$ACC parallel loop  DEFAULT(present) collapse (2) private (i,j,k,medio,Idzek,Idyej)       copyin(Hx,sggMiHx,Ey,Ez,Idye,Idze,b,GM1,GM2) copyout(Hx) 
#endif
         Do k=1,b%sweepHx%NZ
            Do j=1,b%sweepHx%NY
               Do i=1,b%sweepHx%NX
               Idzek=Idze(k)
               Idyej=Idye(j)
                  medio =sggMiHx(i,j,k)
                  Hx(i,j,k)=GM1(MEDIO)*Hx(i,j,k)+GM2(MEDIO)*((Ey(i,j,k+1)-Ey(i,j,k))*Idzek-(Ez(i,j+1,k)-Ez(i,j,k))*Idyej)
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine Advance_Hx

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Advance_Hy(Hy,Ez,Ex,IdzE,IdxE,sggMiHy,b,gm1,gm2)

         !------------------------>
         type (bounds_t), intent( IN)  ::  b
         REAL (KIND=RKIND)     , pointer, dimension ( : )   ::  gm1 ,gm2
         !
         real (kind = RKIND), dimension    ( 0 :     b%dzE%NZ-1     )  , intent( IN)  ::  IdzE
         real (kind = RKIND), dimension    ( 0 :     b%dxE%NX-1     )  , intent( IN)  ::  IdxE
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN)     ::  sggMiHy
         real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 )  , intent( INOUT)  ::  HY
         real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 )  , intent( IN)  ::  EZ
         real (kind = RKIND), dimension    ( 0 :      b%Ex%NX-1 , 0 :      b%Ex%NY-1 , 0 :      b%Ex%NZ-1 )  , intent( IN)  ::  EX
         !------------------------> Variables locales
         real (kind = RKIND)  ::  Idzek
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idzek)     
#endif
#ifdef CompileWithACC   
!$ACC parallel loop DEFAULT(present) collapse (2) private (i,j,k,medio,Idzek)         copyin(Hy,sggMiHy,Ez,Ex,Idze,Idxe,b,GM1,GM2) copyout(Hy) 
#endif
         Do k=1,b%sweepHy%NZ
            Do j=1,b%sweepHy%NY
               Do i=1,b%sweepHy%NX
                  Idzek=Idze(k)
                  medio =sggMiHy(i,j,k)
                  Hy(i,j,k)=GM1(MEDIO)*Hy(i,j,k)+GM2(MEDIO)*((Ez(i+1,j,k)-Ez(i,j,k))*Idxe(i)-(Ex(i,j,k+1)-Ex(i,j,k))*Idzek)
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine Advance_Hy

 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine Advance_Hz(Hz,Ex,Ey,IdxE,IdyE,sggMiHz,b,gm1,gm2)

         !------------------------>
         type (bounds_t), intent( IN)  ::  b
         REAL (KIND=RKIND)     , pointer, dimension ( : )   ::  gm1 ,gm2
         !
         real (kind = RKIND), dimension    ( 0 :     b%dyE%NY-1     )  , intent( IN)  ::  IdyE
         real (kind = RKIND), dimension    ( 0 :     b%dxE%NX-1     )  , intent( IN)  ::  IdxE
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN)     ::  sggMiHz
         real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 )  , intent( INOUT)  ::  Hz
         real (kind = RKIND), dimension    ( 0 :      b%EX%NX-1 , 0 :      b%EX%NY-1 , 0 :      b%EX%NZ-1 )  , intent( IN)  ::  EX
         real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 )  , intent( IN)  ::  EY
         !------------------------> Variables locales
         real (kind = RKIND)  ::  Idyej
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO DEFAULT(SHARED) collapse (2) private (i,j,k,medio,Idyej)  
#endif
#ifdef CompileWithACC   
!$ACC parallel loop  DEFAULT(present) collapse (2) private (i,j,k,medio,Idyej)       copyin(Hz,sggMiHz,Ex,Ey,Idxe,Idye,b,GM1,GM2) copyout(Hz)
#endif
         Do k=1,b%sweepHz%NZ
            Do j=1,b%sweepHz%NY
               Do i=1,b%sweepHz%NX
               Idyej=Idye(j)
                  medio =sggMiHz(i,j,k)
                  Hz(i,j,k)=GM1(MEDIO)*Hz(i,j,k)+GM2(MEDIO)*((Ex(i,j+1,k)-Ex(i,j,k))*Idyej-(Ey(i+1,j,k)-Ey(i,j,k))*Idxe(i))
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine Advance_Hz

      subroutine advanceConformalE()
#ifdef CompileWithConformal
         if(input_conformal_flag)then
            call conformal_advance_E()
         endif
#endif
      end subroutine


      subroutine advanceWires()
         if (( (trim(adjustl(this%control%wiresflavor))=='holland') .or. &
               (trim(adjustl(this%control%wiresflavor))=='transition')) .and. .not. this%control%use_mtln_wires) then
            IF (this%thereAre%Wires) then
               if (this%control%wirecrank) then
                  call AdvanceWiresEcrank(sgg, N, this%control%layoutnumber,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic)
               else
#ifdef CompileWithMTLN
                  if (mtln_parsed%has_multiwires) then
                     write(buff, *) 'ERROR: Multiwires in simulation but -mtlnwires flag has not been selected'
                     call WarnErrReport(buff)
                  end if
#endif
                  call AdvanceWiresE(sgg,N, this%control%layoutnumber,this%control%wiresflavor,this%control%simu_devia,this%control%stochastic,this%control%experimentalVideal,this%control%wirethickness,eps0,mu0)
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
            call AdvanceWiresE_mtln(sgg,Idxh,Idyh,Idzh,eps0,mu0)
#else
            write(buff,'(a)') 'WIR_ERROR: Executable was not compiled with MTLN modules.'
#endif   
         end if

      end subroutine

      !PML E-field advancing (IT IS IMPORTANT TO FIRST CALL THE PML ADVANCING ROUTINES, SINCE THE DISPERSIVE
      !ROUTINES INJECT THE POLARIZATION CURRENTS EVERYWHERE (PML INCLUDED)
      !SO THAT DISPERSIVE MATERIALS CAN ALSO BE TRUNCATED BY CPML)
      subroutine advancePMLE()
         If (this%thereAre%PMLbodies) then !waveport absorbers
            call AdvancePMLbodyE
         endif
         If (this%thereAre%PMLBorders) then
            call AdvanceelectricCPML(sgg%NumMedia, b,sggMiEx,sggMiEy,sggMiEz,G2,Ex,Ey,Ez,Hx,Hy,Hz)
         endif
      end subroutine



      !!!!!!!!!sgg 051214 fill in the magnetic walls after the wireframe info


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine XXXXfillMagnetic(sgg,sggMiEx, sggMiEy, sggMiEz, sggMiHx, sggMiHy, sggMiHz, b)

         !------------------------>
         type (SGGFDTDINFO), intent(IN)    ::  sgg
         type (bounds_t), intent( IN)  ::  b
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( INOUT)     ::  sggMiHx
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( INOUT)     ::  sggMiHy
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( INOUT)     ::  sggMiHz
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( INOUT)     ::  sggMiEx
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( INOUT)     ::  sggMiEy
         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( INOUT)     ::  sggMiEz
         !------------------------> Variables locales
         integer(kind = 4)  ::  i, j, k
         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio1,medio2,medio3,medio4
         logical  ::  mediois1,mediois2,mediois3,mediois4
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4)
#endif
         Do k=1,b%sweepHx%NZ
            Do j=1,b%sweepHx%NY
               Do i=1,b%sweepHx%NX
                  medio1 =sggMiEy(i,j,k)
                  medio2 =sggMiEy(i,j,k+1)
                  medio3 =sggMiEz(i,j,k)
                  medio4 =sggMiEz(i,j+1,k)
                  !mediois1= sgg%med(medio1)%is%already_YEEadvanced_byconformal .or. sgg%med(medio1)%is%split_and_useless .or. (medio1==0)   !!!errror mio de concepto 061214
                  !mediois2= sgg%med(medio2)%is%already_YEEadvanced_byconformal .or. sgg%med(medio2)%is%split_and_useless .or. (medio2==0)
                  !mediois3= sgg%med(medio3)%is%already_YEEadvanced_byconformal .or. sgg%med(medio3)%is%split_and_useless .or. (medio3==0)
                  !mediois4= sgg%med(medio4)%is%already_YEEadvanced_byconformal .or. sgg%med(medio4)%is%split_and_useless .or. (medio4==0)
                  mediois1= (medio1==0)
                  mediois2= (medio2==0)
                  mediois3= (medio3==0)
                  mediois4= (medio4==0)
                  if (mediois1.and.mediois2.and.mediois3.and.mediois4)  sggMiHx(i,j,k)=0
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4)
#endif
         Do k=1,b%sweepHy%NZ
            Do j=1,b%sweepHy%NY
               Do i=1,b%sweepHy%NX
                  medio1 =sggMiEz(i,j,k)
                  medio2 =sggMiEz(i+1,j,k)
                  medio3 =sggMiEx(i,j,k)
                  medio4 =sggMiEx(i,j,k+1)
                  mediois1= (medio1==0)
                  mediois2= (medio2==0)
                  mediois3= (medio3==0)
                  mediois4= (medio4==0)
                  if (mediois1.and.mediois2.and.mediois3.and.mediois4)  sggMiHy(i,j,k)=0
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4)
#endif
         Do k=1,b%sweepHz%NZ
            Do j=1,b%sweepHz%NY
               Do i=1,b%sweepHz%NX
                  medio1 =sggMiEx(i,j,k)
                  medio2 =sggMiEx(i,j+1,k)
                  medio3 =sggMiEy(i,j,k)
                  medio4 =sggMiEy(i+1,j,k)
                  mediois1= (medio1==0)
                  mediois2= (medio2==0)
                  mediois3= (medio3==0)
                  mediois4= (medio4==0)
                  if (mediois1.and.mediois2.and.mediois3.and.mediois4)  sggMiHz(i,j,k)=0
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         !
#ifdef CompileWithOpenMP
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4)
#endif
         Do k=1,b%sweepHx%NZ
            Do j=1,b%sweepHx%NY
               Do i=1,b%sweepHx%NX
                  if ((sggMiHx(i,j,k)==0).or.(sgg%med(sggMiHx(i,j,k))%is%pec))  THEN
                     sggMiEy(i,j,k)   =0
                     sggMiEy(i,j,k+1) =0
                     sggMiEz(i,j,k)   =0
                     sggMiEz(i,j+1,k) =0
                  ENDIF
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4)
#endif
         Do k=1,b%sweepHy%NZ
            Do j=1,b%sweepHy%NY
               Do i=1,b%sweepHy%NX
                  if ((sggMiHy(i,j,k)==0).or.(sgg%med(sggMiHy(i,j,k))%is%pec)) THEN
                     sggMiEz(i,j,k)   =0
                     sggMiEz(i+1,j,k) =0
                     sggMiEx(i,j,k)   =0
                     sggMiEx(i,j,k+1) =0
                  ENDIF
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio1,medio2,medio3,medio4)
#endif
         Do k=1,b%sweepHz%NZ
            Do j=1,b%sweepHz%NY
               Do i=1,b%sweepHz%NX
                  if ((sggMiHz(i,j,k)==0).or.(sgg%med(sggMiHz(i,j,k))%is%pec)) THEN
                     sggMiEx(i,j,k)    =0
                     sggMiEx(i,j+1,k)  =0
                     sggMiEy(i,j,k)    =0
                     sggMiEy(i+1,j,k)  =0
                  ENDIF
               End do
            End do
         End do
#ifdef CompileWithOpenMP
!$OMP  END PARALLEL DO
#endif
         return
      end subroutine XXXXfillMagnetic

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

      
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
      
      
      

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Find the bounds and store everything under the bounds_t variable b
      ! There is redundancy which should be corrected to leave everything in terms of b
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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


   end subroutine launch_simulation




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
   
   !!!!pruebas mergeando bucles (06/09/2016) !no se gana nada (ademas no he validado resultados, solo testeado velocidad)
   
   
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      subroutine Advance_ExEyEz(Ex,Ey,Ez,Hx,Hy,Hz,Idxh,Idyh,Idzh,sggMiEx,sggMiEy,sggMiEz,b,g1,g2)
!!
!!         !------------------------>
!!         type (bounds_t), intent( IN)  ::  b
!!         REAL (KIND=RKIND)     , pointer, dimension ( : )  ::  g1, g2
!!         !
!!         real (kind = RKIND), dimension    ( 0 :     b%dxh%NX-1     )  , intent( IN)  ::  Idxh
!!         real (kind = RKIND), dimension    ( 0 :     b%dyh%NY-1     )  , intent( IN)  ::  Idyh
!!         real (kind = RKIND), dimension    ( 0 :     b%dzh%NZ-1     )  , intent( IN)  ::  Idzh
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN)     ::  sggMiEx
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( IN)     ::  sggMiEy
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN)     ::  sggMiEz
!!         real (kind = RKIND), dimension    ( 0 :      b%Ex%NX-1 , 0 :      b%Ex%NY-1 , 0 :      b%Ex%NZ-1 )  , intent( INOUT)  ::  Ex
!!         real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 )  , intent( INOUT)  ::  EY
!!         real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 )  , intent( INOUT)  ::  Ez
!!         real (kind = RKIND), dimension    ( 0 :      b%Hx%NX-1 , 0 :      b%Hx%NY-1 , 0 :      b%Hx%NZ-1 )  , intent( IN)  ::  HX
!!         real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 )  , intent( IN)  ::  HY
!!         real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 )  , intent( IN)  ::  HZ
!!         !------------------------> Variables locales
!!         real (kind = RKIND)  ::  Idzhk, Idyhj
!!         integer(kind = 4)  ::  i, j, k
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
!!         
!!         
!!#ifdef CompileWithOpenMP
!!!$OMP  PARALLEL DO DEFAULT(SHARED) private (i,j,k,medio,Idzhk,Idyhj)
!!#endif
!!         Do k=1,max(b%sweepEx%NZ,b%sweepEy%NZ,b%sweepEz%NZ)
!!            Idzhk=Idzh(k)
!!            Do j=1,max(b%sweepEx%NY,b%sweepEy%NY,b%sweepEz%NY)
!!               Idyhj=Idyh(j)
!!               Do i=1,max(b%sweepEx%NX,b%sweepEy%NX,b%sweepEz%NX)
!!                  medio =sggMiEx(i,j,k)
!!                  Ex(i,j,k)=G1(MEDIO)*Ex(i,j,k)+G2(MEDIO)*((Hz(i,j,k)-Hz(i,j-1,k))*Idyhj-(Hy(i,j,k)-Hy(i,j,k-1))*Idzhk)
!!                  medio =sggMiEy(i,j,k)
!!                  Ey(i,j,k)=G1(MEDIO)*Ey(i,j,k)+G2(MEDIO)*((Hx(i,j,k)-Hx(i,j,k-1))*Idzhk-(Hz(i,j,k)-Hz(i-1,j,k))*Idxh(i))
!!                  medio =sggMiEz(i,j,k)
!!                  Ez(i,j,k)=G1(MEDIO)*Ez(i,j,k)+G2(MEDIO)*((Hy(i,j,k)-Hy(i-1,j,k))*Idxh(i)-(Hx(i,j,k)-Hx(i,j-1,k))*Idyhj)
!!               End do
!!            End do
!!         End do
!!#ifdef CompileWithOpenMP
!!!$OMP  END PARALLEL DO
!!#endif
!!         return
!!      end subroutine Advance_ExEyEz
!!
!!
!!
!!
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      subroutine Advance_HxHyHz(Hx,Hy,Hz,Ex,Ey,Ez,IdxE,IdyE,IdzE,sggMiHx,sggMiHy,sggMiHz,b,gm1,gm2)
!!
!!         !------------------------>
!!         type (bounds_t), intent( IN)  ::  b
!!         REAL (KIND=RKIND)     , pointer, dimension ( : )   ::  gm1 ,gm2
!!         !!
!!         real (kind = RKIND), dimension    ( 0 :     b%dxE%NX-1    )  , intent( IN)  ::  IdxE
!!         real (kind = RKIND), dimension    ( 0 :     b%dyE%NY-1    )  , intent( IN)  ::  IdyE
!!         real (kind = RKIND), dimension    ( 0 :     b%dzE%NZ-1    )  , intent( IN)  ::  IdzE
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN)     ::  sggMiHx
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN)     ::  sggMiHy
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension ( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN)     ::  sggMiHz
!!         real (kind = RKIND), dimension    ( 0 :      b%Hx%NX-1 , 0 :      b%Hx%NY-1 , 0 :      b%Hx%NZ-1 )  , intent( INOUT)  ::  Hx
!!         real (kind = RKIND), dimension    ( 0 :      b%Hy%NX-1 , 0 :      b%Hy%NY-1 , 0 :      b%Hy%NZ-1 )  , intent( INOUT)  ::  HY
!!         real (kind = RKIND), dimension    ( 0 :      b%Hz%NX-1 , 0 :      b%Hz%NY-1 , 0 :      b%Hz%NZ-1 )  , intent( INOUT)  ::  Hz
!!         real (kind = RKIND), dimension    ( 0 :      b%Ex%NX-1 , 0 :      b%Ex%NY-1 , 0 :      b%Ex%NZ-1 )  , intent( IN)  ::  EX
!!         real (kind = RKIND), dimension    ( 0 :      b%Ey%NX-1 , 0 :      b%Ey%NY-1 , 0 :      b%Ey%NZ-1 )  , intent( IN)  ::  EY
!!         real (kind = RKIND), dimension    ( 0 :      b%Ez%NX-1 , 0 :      b%Ez%NY-1 , 0 :      b%Ez%NZ-1 )  , intent( IN)  ::  EZ
!!         !------------------------> Variables locales
!!         real (kind = RKIND)  ::  Idzek, Idyej
!!         integer(kind = 4)  ::  i, j, k
!!         integer(kind = INTEGERSIZEOFMEDIAMATRICES)  ::  medio
!!#ifdef CompileWithOpenMP
!!!$OMP  PARALLEL DO  DEFAULT(SHARED) private (i,j,k,medio,Idzek,Idyej)
!!#endif
!!         Do k=1,max(b%sweepHx%NZ,b%sweepHy%NZ,b%sweepHz%NZ)
!!            Idzek=Idze(k)
!!            Do j=1,max(b%sweepHx%NY,b%sweepHy%NY,b%sweepHz%NY)
!!               Idyej=Idye(j)
!!               Do i=1,max(b%sweepHx%NX,b%sweepHy%NX,b%sweepHz%NX)
!!                  medio =sggMiHx(i,j,k)
!!                  Hx(i,j,k)=GM1(MEDIO)*Hx(i,j,k)+GM2(MEDIO)*((Ey(i,j,k+1)-Ey(i,j,k))*Idzek-(Ez(i,j+1,k)-Ez(i,j,k))*Idyej)
!!                  medio =sggMiHy(i,j,k)
!!                  Hy(i,j,k)=GM1(MEDIO)*Hy(i,j,k)+GM2(MEDIO)*((Ez(i+1,j,k)-Ez(i,j,k))*Idxe(i)-(Ex(i,j,k+1)-Ex(i,j,k))*Idzek)
!!                  medio =sggMiHz(i,j,k)
!!                  Hz(i,j,k)=GM1(MEDIO)*Hz(i,j,k)+GM2(MEDIO)*((Ex(i,j+1,k)-Ex(i,j,k))*Idyej-(Ey(i+1,j,k)-Ey(i,j,k))*Idxe(i))
!!               End do
!!            End do
!!         End do
!!#ifdef CompileWithOpenMP
!!!$OMP  END PARALLEL DO
!!#endif
!!         return
!!      end subroutine Advance_HxHyHz


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
   

end module
