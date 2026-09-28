
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains the types and parameters shared by all the rest of the modules
! No public variables are defined. Only types and parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module  FDETYPES_m


#ifdef CompileWithOpenMP
   use omp_lib
#endif

 
#ifdef CompileWithReal16 
#undef CompileWithReal8
#undef CompileWithReal4
#endif

#ifdef CompileWithReal8 
#undef CompileWithReal16
#undef CompileWithReal4
#endif

#ifndef CompileWithReal16
#ifndef CompileWithReal8
#ifndef CompileWithReal4
#define CompileWithReal4
#endif
#endif
#endif


#ifndef CompileWithInt4
#ifndef CompileWithInt2
#ifndef CompileWithInt1
#define CompileWithInt4
#endif
#endif
#endif


#ifdef CompileWithMPI
   use MPI
   implicit none
#endif

   !Every type and parameter is public
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   public
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Tunable Parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer(kind=4) :: quienmpi,tamaniompi
   integer(kind=4) :: SUBCOMM_MPI
!240424 para que funcionen las sondas slice de conformal lo pongo como general. niapaa. algun dia hay que reahacer el conformal 
   !y esto debe desaparecer
   integer(kind=4) :: SUBCOMM_MPI_conformal_probes,MPI_conformal_probes_root
!!!
   integer(kind=8),  parameter  :: MAXMPIBYTES = 2**27
   integer(kind=4),  parameter  :: BUFFOBSE=2**10 !Steps of the temporal buffer to store evolution data
   integer(kind=8),  parameter  :: MAXMEMORYPROBES=2_8**37_8 !128 Gb Maximum bytes of the buffer to store evolution data
   integer(kind=8),  parameter  :: MAXPROBES=150000 !Maximum number of probes (a limit of 200000 is set with ulimit in Linux)
   !
   !
   integer, parameter :: TOPCPUTIME=10000000 !maximum cpu time in minutes 
   !size of character strings 
   integer, parameter :: BUFSIZE=1024
   integer, parameter :: BUFSIZE_LONG=16384
   !!!integer :: maxmessages=20000 !numero maximo mensajes para alocatear en MPI overrideable con -maxmessages y quitado como parameter fijo !deprecated 07/03/15
   !dxf output stuff
   !!!integer, parameter :: maxdxf= 20000,dxflinesize=14
   !!!character(len=dxflinesize) :: dxfbuff
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !Rest of Parameters
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef CompileWithInt1
#undef CompileWithInt2
#undef CompileWithInt4
   integer(kind=2), parameter  :: INTEGERSIZEOFMEDIAMATRICES=1
#ifdef CompileWithMPI
   integer(kind=4), parameter  :: INTEGERSIZE=MPI_INTEGER1
#endif
#endif
#ifdef CompileWithInt2
#undef CompileWithInt1
#undef CompileWithInt4
   integer(kind=2), parameter  :: INTEGERSIZEOFMEDIAMATRICES=2
#ifdef CompileWithMPI
   integer(kind=4), parameter  :: INTEGERSIZE=MPI_INTEGER2
#endif
#endif
#ifdef CompileWithInt4
#undef CompileWithInt1
#undef CompileWithInt2
   integer(kind=2), parameter  :: INTEGERSIZEOFMEDIAMATRICES=4
#ifdef CompileWithMPI
   integer(kind=4), parameter  :: INTEGERSIZE=MPI_INTEGER4
#endif
#endif
   integer(kind=4), parameter  :: IKINDMTAG=4 !PARA SGGMTAG 151020 !dejarlo en 4 bytes. No tocar

   integer(kind=2), parameter  :: SINGLE=4
   integer(kind=2), parameter  :: DOUBLE_KIND=8
   integer(kind=2), parameter  :: LONG_DOUBLE=16
#ifdef CompileWithReal8
   integer(kind=2), parameter  :: RKIND=DOUBLE_KIND
   integer(kind=2), parameter  :: RKIND_WIRES=DOUBLE_KIND
   integer(kind=2), parameter  :: RKIND_TIEMPO=DOUBLE_KIND
   integer(kind=2), parameter  :: CKIND=DOUBLE_KIND
#else
#ifdef CompileWithReal16
   integer(kind=2), parameter  :: RKIND=LONG_DOUBLE
   integer(kind=2), parameter  :: RKIND_WIRES=LONG_DOUBLE
   integer(kind=2), parameter  :: RKIND_TIEMPO=LONG_DOUBLE
   integer(kind=2), parameter  :: CKIND=LONG_DOUBLE
#else
   !default
   integer(kind=2), parameter  :: RKIND=SINGLE
   integer(kind=2), parameter  :: RKIND_WIRES=DOUBLE_KIND !020719 a peticion 
   integer(kind=2), parameter  :: RKIND_TIEMPO=DOUBLE_KIND
   !! integer(kind=2), parameter  :: CKIND=SINGLE
   integer(kind=2), parameter  :: CKIND=DOUBLE_KIND  !LOS COMPLEJOS LOS VOY A MANEJAR SIEMPRE EN DOBLE PRECISION como minimo
#endif
#endif

   !

#ifdef CompileWithMPI
   real(kind=RKIND), parameter  :: PLUSCPU_PML=2.0_RKIND !heuristic (1=No CPU overhead, 2=double CPU overhead)
#endif
#ifdef CompileWithMPI
#ifdef CompileWithReal8
   integer(kind=4), parameter  :: REALSIZE=MPI_DOUBLE_PRECISION
   integer(kind=4), parameter  :: REALSIZE_WIRES=MPI_DOUBLE_PRECISION
   integer(kind=4), parameter  :: COMPLEXSIZE=MPI_DOUBLE_COMPLEX
   integer(kind=4), parameter  :: REALSIZE_TIEMPO=MPI_DOUBLE_PRECISION
#else
#ifdef CompileWithReal16
   integer(kind=4), parameter  :: REALSIZE=MPI_REAL16
   integer(kind=4), parameter  :: COMPLEXSIZE=MPI_COMPLEX32
   integer(kind=4), parameter  :: REALSIZE_TIEMPO=MPI_REAL_16
#else
   integer(kind=4), parameter  :: REALSIZE=MPI_REAL
   integer(kind=4), parameter  :: REALSIZE_WIRES=MPI_DOUBLE_PRECISION
   integer(kind=4), parameter  :: REALSIZE_TIEMPO=MPI_DOUBLE_PRECISION

!!!   integer(kind=4), parameter  :: COMPLEXSIZE=MPI_COMPLEX
   integer(kind=4), parameter  :: COMPLEXSIZE=MPI_DOUBLE_COMPLEX  !LOS COMPLEJOS LOS VOY A MANEJAR SIEMPRE EN DOBLE PRECISION como minimo !esto debe ir ligado a la definicion de ckind
#endif
#endif
#endif
   real(kind=RKIND) , parameter  :: HEURCFL=0.8_RKIND
   real(kind=RKIND) , parameter  :: &
   pi=3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148, &
   unmedio = 0.5_RKIND
   complex(kind=CKIND), parameter :: MCPI2 = - (0.0_RKIND, 1.0_RKIND) * 2.0_RKIND * pi;
   !
   integer(kind=4), parameter  :: DOWN=1, UP=2,  LEFT=3, RIGHT=4, BACK=5, FRONT=6
   !
   integer(kind=4),  parameter  :: IEX=1,IEY=2,IEZ=3,IHX=4,IHY=5,IHZ=6,CENTROIDE=8,NOTHING=666
   !
   integer(kind=4),  parameter  :: IMEC=51 !modulus, TANGENTIAL, NORMAL fields in cuts for Volumic probes
   integer(kind=4),  parameter  :: IMHC=52
   integer(kind=4),  parameter  :: ICUR=53 !Bloque currents along edges in thin wires, PEC and surface edges
   integer(kind=4),  parameter  :: ICURX=54 !Bloque currents along edges in surface with normal X
   integer(kind=4),  parameter  :: ICURY=55 !Bloque currents along edges in surface with normal Y
   integer(kind=4),  parameter  :: ICURZ=56 !Bloque currents along edges in surface with normal Z
   integer(kind=4),  parameter  :: MAPVTK=57 !Bloque currents along edges in surface with normal Z
   integer(kind=4),  parameter  :: IEXC=61 !components in cuts for Volumic probes
   integer(kind=4),  parameter  :: IEYC=62
   integer(kind=4),  parameter  :: IEZC=63
   integer(kind=4),  parameter  :: IHXC=64
   integer(kind=4),  parameter  :: IHYC=65
   integer(kind=4),  parameter  :: IHZC=66
   integer(kind=4),  parameter  :: FARFIELD=67
   integer(kind=4),  parameter  :: LINEINTEGRAL=68
   ! do not change
   integer(kind=4),  parameter  :: IJX=10*iEx,IJY=10*iEy,IJZ=10*IEZ
   integer(kind=4),  parameter  :: IQX=10000*iEx,IQY=10000*iEy,IQZ=10000*IEZ
   integer(kind=4),  parameter  :: IVX=1000*iEx,IVY=1000*iEy,IVZ=1000*IEZ
   integer(kind=4),  parameter  :: IBLOQUEJX=100*iEx,IBLOQUEJY=100*iEy,IBLOQUEJZ=100*IEZ
   integer(kind=4),  parameter  :: IBLOQUEMX=100*IHX,IBLOQUEMY=100*IHY,IBLOQUEMZ=100*IHZ
   !
   integer(kind=4), parameter :: VOLUMIC_M_MEASURE(3) = [ICUR, IMEC, IMHC]
   integer(kind=4), parameter :: VOLUMIC_X_MEASURE(3) = [ICURX, IEXC, IHXC]
   integer(kind=4), parameter :: VOLUMIC_Y_MEASURE(3) = [ICURY, IEYC, IHYC]
   integer(kind=4), parameter :: VOLUMIC_Z_MEASURE(3) = [ICURZ, IEZC, IHZC]

   integer(kind=4), parameter :: ELECTRIC_FIELD_DIRECTION(3) = [iEx, iEy, IEZ]
   integer(kind=4), parameter :: MAGNETIC_FIELD_DIRECTION(3) = [IHX, IHY, IHZ]
   integer(kind=4), parameter :: CURRENT_MEASURE(4) = [ICUR, ICURX, ICURY, ICURZ]
   integer(kind=4), parameter :: ELECTRIC_FIELD_MEASURE(4) = [IMEC, IEXC, IEYC, IEZC]
   integer(kind=4), parameter :: MAGNETIC_FIELD_MEASURE(4) = [IMHC, IHXC, IHYC, IHZC]
   !
   character(len=*), parameter  :: SEPARADOR='______________'
   integer(kind=4), parameter  :: COMI=1,FINE=2, ICOORD=1,JCOORD=2,KCOORD=3

   real(kind=RKIND), parameter :: EPSILON_VACUUM   =   &
   8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12
   real(kind=RKIND), parameter :: MU_VACUUM        =   &
   1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6
   real(kind=rkind), parameter :: C_VACUUM = 1.0_RKIND/sqrt(EPSILON_VACUUM*MU_VACUUM)
   
   real(kind=RKIND_TIEMPO) :: dt0 !aqui para OLDrlo accesible en resume pscale
   
   integer(kind=4), parameter :: FACE_X = 1
   integer(kind=4), parameter :: FACE_Y = 2
   integer(kind=4), parameter :: FACE_Z = 3
   
   integer(kind=4), parameter :: EDGE_X = 1
   integer(kind=4), parameter :: EDGE_Y = 2
   integer(kind=4), parameter :: EDGE_Z = 3

   !source types
   character(len=*), parameter :: F_SOURCE_VOLTAGE = 'VOLT'
   character(len=*), parameter :: F_SOURCE_CURRENT = 'CURR'

   
#ifdef CompileWithReal4
   character(len=*), parameter  :: FMT='(e27.17e3,11(e19.9e3))'  !IEEE 754 single-precision 6 to 9 decimals -1.123456789E-001
#else
#ifdef CompileWithReal8
   character(len=*), parameter  :: FMT='(12(e27.17e3))' !IEEE 754 single-precision 15 to 17 decimals 
#else   
#ifdef CompileWithReal16
   character(len=*), parameter  :: FMT='(12(e46.36e3))'  !IEEE 754 single-precision 33 to 36 decimals  
#else !default
   character(len=*), parameter  :: FMT='(e27.17e3,11(e19.9e3))'  !IEEE 754 single-precision 6 to 9 decimals -1.123456789E-001
#endif
#endif
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !solo tipos
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    type tagtype_t
        character(len=BUFSIZE), allocatable, dimension(:) :: tag
        integer(kind=4) :: numertags
    end type

   type coorsxyz_t
      real(kind=RKIND), pointer, dimension(:) :: x,y,z
   end type coorsxyz_t
   !
   type coorsxyzP_t
      type(coorsxyz_t), dimension(1:6) :: PhysCoor
   end type coorsxyzP_t

   type ExtraMedium_t
      integer(kind=4) :: pml_size,elementIndex
      real(kind=rkind) :: sigma,sigmam
      logical :: exists
   end type
   type logic_control_t
      logical  :: Wires  , &
      PMLbodies  , &
      MultiportS  , &
      AnisMultiportS  , &
      SGBCs , &
      Lumpeds , &
      EDispersives  , &
      MDispersives  , &
      PlaneWaveBoxes  , &
      Observation  , &
      FarFields  , &
      PMCBorders  , &
      PMLBorders  , &
      MurBorders  , &
      PECBorders  , &
      PeriodicBorders, &
      Anisotropic  , &
      ThinSlot  , &
      NodalE  , &
      NodalH  , &
      MagneticMedia, PMLMagneticMedia, &
      MTLNbundles
   contains 
      procedure :: reset => logic_reset
   end type



   !computational limits
   type Xlimit_t
      integer(kind=4) :: XI,XE,NX
   end type
   type Ylimit_t
      integer(kind=4) :: YI,YE,NY
   end type
   type Zlimit_t
      integer(kind=4) :: ZI,ZE,NZ
   end type
   type limit_t
      integer(kind=4) :: XI,XE,YI,YE,ZI,ZE,NX,NY,NZ
   end type
   type XYZlimit_t
      integer(kind=4) :: XI,XE,YI,YE,ZI,ZE
   end type
   type xyzlimit_scaled_t
      integer(kind=4) :: XI,XE,YI,YE,ZI,ZE
      real(kind=RKIND) :: xc,yc,zc
      integer(kind=4) :: Or   !to include possible orientations (nodal sources 180915)
   end type

   type tagnumber_t
      integer(kind=IKINDMTAG) , allocatable , dimension(:,:,:) :: x, y, z
   end type

   type taglist_t
      type(tagnumber_t) :: edge, face
   contains
      private
      procedure, public :: getFaceTag => taglist_getFaceTag
      procedure, public :: getEdgeTag => taglist_getEdgeTag
   end type

   !
   type bounds_t
      type(limit_t) :: sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz
      type(limit_t) :: Ex,Ey,Ez,Hx,Hy,Hz
      type(limit_t) :: sweepEx,sweepEy,sweepEz,sweepHx,sweepHy,sweepHz
      type(limit_t) :: sweepSINPMLEx,sweepSINPMLEy,sweepSINPMLEz,sweepSINPMLHx,sweepSINPMLHy,sweepSINPMLHz
      type(Xlimit_t) :: dxe,dxh
      type(Ylimit_t) :: dye,dyh
      type(Zlimit_t) :: dze,dzh
   end type
   type  :: PML_t
      real(kind=RKIND) :: orden(3,2)
      real(kind=RKIND) :: CoeffReflPML(3,2) !(icor : jcor : kcor,start : ende)
      integer(kind=4) :: NumLayers(3,2)
   end type
   !

   type  :: fichevol_t
      character(len=BUFSIZE) :: Name
      integer(kind=4) :: NumSamples
      real(kind=RKIND) :: DeltaSamples
      real(kind=RKIND), dimension(:), pointer  :: Samples
   end type
   !

   !wires
   type  :: fichevol_wires_t
      character(len=BUFSIZE) :: Name
      integer(kind=4) :: NumSamples
      real(kind=RKIND_WIRES) :: DeltaSamples
      real(kind=RKIND_WIRES), dimension(:), pointer  :: Samples
   end type
   type  :: source_t
      type(fichevol_wires_t) :: sourceFile
      real(kind=RKIND_WIRES) :: Resistance
      real(kind=RKIND_WIRES) :: Multiplier
      integer(kind=4) :: i,j,k
   end type

   type  :: NodalSource_t
      type(fichevol_t) :: sourceFile
      type(xyzlimit_scaled_t), pointer, dimension(:) :: punto
      integer(kind=4) :: numpuntos
      logical :: IsInitialValue
      logical :: IsHard
      logical :: IsElec
   end type NodalSource_t
   !
   type  :: WireDispersiveParams_t
      integer(kind=4)                            :: numPoles
      complex(kind=CKIND), pointer, dimension(:) :: res, p
      complex(kind=CKIND)                        :: d, e
   end type

   type  :: oriented_point_t
      integer(kind=4) :: ori
      integer(kind=4) :: i,j,k,origIndex,ilibre,jlibre,klibre,multiraboDE !si es multirabo de que indice lo es
      logical :: Is_LeftEnd,Is_RightEnd,IsEnd_norLeft_norRight
      logical :: repetido,multirabo !marca segmentos que aparecen repetidos en un mismo thin wire!los bundles deberan estar thin-wires distintos
      logical :: orientadoalreves
   end type oriented_point_t

#ifdef CompileWithMTLN   
   type  :: Multiwires_t
   end type
#endif

   type  :: Wires_t
      real(kind=RKIND_WIRES) :: Radius,R,L,C,P_R,P_L,P_C
      real(kind=RKIND_WIRES) :: Radius_devia,R_devia,L_devia,C_devia
      type(WireDispersiveParams_t), allocatable, dimension(:) :: disp
      integer(kind=4) :: numsegmentos,NUMVOLTAGESOURCES,NUMCURRENTSOURCES
      type(oriented_point_t), pointer, dimension(:) :: segm
      type(source_t), pointer, dimension(:) :: Vsource
      type(source_t), pointer, dimension(:) :: Isource
      logical  :: VsourceExists ,IsourceExists
      logical  :: HasParallel_LeftEnd ,HasParallel_RightEnd ,&
                   HasSeries_LeftEnd ,HasSeries_RightEnd,HasAbsorbing_LeftEnd,HasAbsorbing_RightEnd
      real(kind=RKIND_WIRES) :: Parallel_R_RightEnd,Parallel_R_LeftEnd
      real(kind=RKIND_WIRES) :: Series_R_RightEnd,Series_R_LeftEnd
      real(kind=RKIND_WIRES) :: Parallel_L_RightEnd,Parallel_L_LeftEnd
      real(kind=RKIND_WIRES) :: Series_L_RightEnd,Series_L_LeftEnd
      real(kind=RKIND_WIRES) :: Parallel_C_RightEnd,Parallel_C_LeftEnd
      real(kind=RKIND_WIRES) :: Series_C_RightEnd,Series_C_LeftEnd
!
      real(kind=RKIND_WIRES) :: Parallel_R_RightEnd_devia ,Parallel_R_LeftEnd_devia
      real(kind=RKIND_WIRES) :: Series_R_RightEnd_devia ,  Series_R_LeftEnd_devia
      real(kind=RKIND_WIRES) :: Parallel_L_RightEnd_devia ,Parallel_L_LeftEnd_devia
      real(kind=RKIND_WIRES) :: Series_L_RightEnd_devia ,  Series_L_LeftEnd_devia
      real(kind=RKIND_WIRES) :: Parallel_C_RightEnd_devia ,Parallel_C_LeftEnd_devia
      real(kind=RKIND_WIRES) :: Series_C_RightEnd_devia ,  Series_C_LeftEnd_devia
      type(WireDispersiveParams_t), allocatable, dimension(:) :: disp_LeftEnd, disp_RightEnd
      ! integer(kind=4) :: LextremoI,LextremoJ,LextremoK,RextremoI,RextremoJ,RextremoK !no ncesario: yo luego calculo bien los extremos
      integer(kind=4) :: LeftEnd,RightEnd
   end type Wires_t
   
   type  :: SlantedNode_t
      integer(kind=4) :: elementIndex
      real(kind=RKIND_WIRES) :: x, y, z
      logical                 :: VsourceExists, IsourceExists
      type(source_t), pointer  :: Vsource, Isource
   end type SlantedNode_t
   
   type  :: SlantedWires_t
      real(kind=RKIND_WIRES) :: radius,R,L,C,P_R,P_L,P_C
      type(WireDispersiveParams_t), allocatable, dimension(:) :: disp
      integer(kind=4) :: LeftEnd, RightEnd
      integer(kind=4) :: NumNodes
      type(SlantedNode_t), pointer, dimension(:) :: nodes
      logical           :: HasParallel_LeftEnd
      real(kind=RKIND_WIRES) :: Parallel_R_LeftEnd, Parallel_L_LeftEnd, Parallel_C_LeftEnd
      logical           :: HasParallel_RightEnd
      real(kind=RKIND_WIRES) :: Parallel_R_RightEnd, Parallel_L_RightEnd, Parallel_C_RightEnd
      logical           :: HasSeries_LeftEnd
      real(kind=RKIND_WIRES) :: Series_R_LeftEnd, Series_L_LeftEnd, Series_C_LeftEnd
      logical           :: HasSeries_RightEnd
      real(kind=RKIND_WIRES) :: Series_R_RightEnd, Series_L_RightEnd, Series_C_RightEnd
      type(WireDispersiveParams_t), allocatable, dimension(:) :: disp_LeftEnd, disp_RightEnd
   end type SlantedWires_t
   !
   type  :: Lumped_t
      integer(kind=4) :: Orient = 0 !orientation +iEx, -iEx,+iEy.......
!deprecado 201222      real(kind=RKIND_wires) :: epr,mur,sigma,sigmam
      real(kind=RKIND_WIRES) :: R,L,C,DiodB,DiodIsat,Rtime_on,Rtime_off
      logical :: resistor , inductor , capacitor , diode 
      real(kind=RKIND_WIRES) ::R_devia,L_devia,C_devia
   end type Lumped_t
!!!
   !end wires
   type  :: PMLbody_t
      integer(kind=4) :: orient = 0 !orientation +iEx, -iEx,+iEy.......el signo aqui es intranscendente
   end type PMLbody_t
!!!
   type  :: Multiport_t
      integer(kind=4) :: Multiportdir = 0 !orientation +iEx, -iEx,+iEy.......
      character(len=BUFSIZE)                            :: multiportFileZ11,multiportFileZ22,multiportFileZ12,multiportFileZ21
      real(kind=rkind), dimension(:), pointer :: epr,mur,sigma,sigmam,width   
                  !_for_devia 090519
      real(kind=rkind), dimension(:), pointer :: epr_devia,mur_devia,sigma_devia,sigmam_devia,width_devia
                  !!!
!!old pre 17/08/115: no es valido para mallados NO uniformes. Hay que hacerlo punto a punto
!!!                     real(kind=rkind) :: transversalSpaceDelta
      integer(kind=4) :: numcapas
   end type Multiport_t
   !
   type  :: AnisMultiport_t
      integer(kind=4) :: Multiportdir = 0 !orientation +iEx, -iEx,+iEy.......
      character(len=BUFSIZE)                            :: MultiportFileZ11,MultiportFileZ22, &
      MultiportFileZ12,MultiportFileZ21
      real(kind=rkind), pointer, dimension(:) :: epr,mur,sigma,sigmam,width
   end type AnisMultiport_t
   !
   type planeonde_t
      real(kind=RKIND) :: INCERTMAX
      real(kind=RKIND), allocatable, dimension(:) :: px,py,pz,ex,ey,ez,incert
      integer(kind=4) :: esqx1,esqy1,esqz1,esqx2,esqy2,esqz2
      type(fichevol_t) :: sourceFile
      integer(kind=4) :: nummodes
      logical :: isRC 
   end type planeonde_t
   !
   type  :: Border_t
      logical  :: IsBackPEC , &
      IsFrontPEC , &
      IsLeftPEC , &
      IsRightPEC , &
      IsUpPEC , &
      IsDownPEC , &
      IsBackPMC , &
      IsFrontPMC , &
      IsLeftPMC , &
      IsRightPMC , &
      IsUpPMC , &
      IsDownPMC , &
      IsBackPML , &
      IsFrontPML , &
      IsLeftPML , &
      IsRightPML , &
      IsUpPML , &
      IsDownPML , &
      IsBackPeriodic , &
      IsFrontPeriodic , &
      IsLeftPeriodic , &
      IsRightPeriodic , &
      IsUpPeriodic , &
      IsDownPeriodic, &
      IsBackMUR , &
      IsFrontMUR , &
      IsLeftMUR , &
      IsRightMUR , &
      IsUpMUR , &
      IsDownMUR
   end type
   !

   type, public :: direction_t
      integer(kind=4) :: x,y,z, orientation
   contains
      private
      procedure :: direction_eq
      generic, public :: operator(==) => direction_eq
   end type

   type  :: observable_t
      integer(kind=4) :: XI,YI,ZI,XE,YE,ZE,What,Node  !los valores finales XE,YE,ZE solo se precisan para las CurrentProbes
      integer(kind=4) :: Xtrancos,Ytrancos,Ztrancos
      type(direction_t), dimension(:), allocatable :: line
      
   end type observable_t
   !
   type  :: Obses_t
      integer(kind=4) :: nP
      type(observable_t), pointer, dimension(:) :: P
      real(kind=RKIND) :: InitialTime,FinalTime,TimeStep
      real(kind=RKIND) :: InitialFreq,FinalFreq,FreqStep

      real(kind=RKIND) :: thetaStart,thetaStop,thetaStep
      real(kind=RKIND) :: phiStart,phiStop,phiStep

      character(len=BUFSIZE) :: outputrequest
      character(len=BUFSIZE) :: FileNormalize
      logical :: FreqDomain ,TimeDomain , Saveall,  &
      transferFlag, Volumic,Done,Begun,Flushed
   end type

   type SharedElement_t
      integer(kind=4) :: i,j,k,field,PropMed,SharedMed,times !field(i,j,k)=PropMed shares ShareMed
   end type
   type Shared_t
      integer(kind=4) :: Conta = 0, MaxConta = 10
      type(SharedElement_t), pointer, dimension(:) :: elem
   end type


   type  :: DispersiveParams_t
      integer(kind=4) :: NumPolRes11,NumPolRes12,NumPolRes13,NumPolRes22,NumPolRes23,NumPolRes33
      complex(kind=CKIND), pointer, dimension(:) :: C11,A11,C12,A12,C13,A13,C22,A22,C23,A23,C33,A33
      real(kind=RKIND) :: eps11,MU11,SIGMA11,SIGMAM11
      real(kind=RKIND) :: eps12,MU12,SIGMA12,SIGMAM12
      real(kind=RKIND) :: EPs13,MU13,SIGMA13,SIGMAM13
      real(kind=RKIND) :: EPs22,MU22,SIGMA22,SIGMAM22
      real(kind=RKIND) :: EPs23,MU23,SIGMA23,SIGMAM23
      real(kind=RKIND) :: EPs33,MU33,SIGMA33,SIGMAM33
   end type

   type :: Anisotropic_t
      real(kind=RKIND),  dimension(3,3) :: sigma,epr,mur,sigmaM
   end type


   type Exists_t
      logical                    :: &
      PML , &
      PEC , &
      ConformalPEC , &
      PMC , &
      ThinWire , &
      Multiwire , &
      SlantedWire, &
      EDispersive , &
      MDispersive , &
      EDispersiveAnis , &
      MDispersiveAnis , &
      ThinSlot , &
      PMLbody , &
      SGBC , &
      SGBCDispersive , &
      Lumped , &
      Lossy, &
      AnisMultiport , &
      Multiport , &
      MultiportPadding , &
      DIELECTRIC , &
      Anisotropic , &
      Volume , &
      Line , &
      Surface , &
      Needed , &
      Interfase,&
      already_YEEadvanced_byconformal,  &
      split_and_useless
   end type



   type  :: MediaData_t
      integer(kind=SINGLE) :: Id
      real(kind=RKIND) :: Priority,Epr,Sigma,Mur,SigmaM
      logical :: sigmareasignado !solo afecta a un chequeo de errores en lumped 120123
      type(exists_t)            :: Is
      type(Wires_t)           , dimension(:), pointer  :: Wire
      type(SlantedWires_t)    , dimension(:), pointer  :: SlantedWire
      type(PMLbody_t)         , dimension(:), pointer  :: PMLbody
      type(Multiport_t)       , dimension(:), pointer  :: Multiport
      type(AnisMultiport_t)   , dimension(:), pointer  :: AnisMultiport
      type(DispersiveParams_t), dimension(:), pointer  :: EDispersive
      type(DispersiveParams_t), dimension(:), pointer  :: MDispersive
      type(Anisotropic_t)     , dimension(:), pointer  :: Anisotropic
      type(Lumped_t)          , dimension(:), pointer  :: Lumped
#ifdef CompileWithMTLN
      type(Multiwires_t)      , dimension(:), pointer  :: Multiwire
#endif
   end type

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! This is the  class which stores all the simulation data
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   type  :: SGGFDTDINFO_t
      real(kind=RKIND_TIEMPO)     , pointer, dimension(:) :: time !para permit scaling
      real(kind=RKIND_TIEMPO) :: dt
      character(len=BUFSIZE) :: extraSwitches
      !!
      integer(kind=4) :: NumMedia,AllocMed
      integer(kind=4) :: IniPMLMedia,EndPMLMedia
      integer(kind=4) :: NumPlaneWaves,TimeSteps,InitialTimeStep
      integer(kind=4) :: NumNodalSources
      integer(kind=4) :: NumberRequest
      !!!
      real(kind=RKIND)     , pointer, dimension(:) :: LineX,LineY,LineZ
      real(kind=RKIND)     , pointer, dimension(:) :: DX,DY,DZ
      integer(kind=4)                                        :: AllocDxI,AllocDyI,AllocDzI,AllocDxE,AllocDyE,AllocDzE
      type(planeonde_t), pointer, dimension(:)            :: PlaneWave
      type(Border_t)                                         :: Border
      type(PML_t)                                            :: PML
      !    !
      type(Shared_t)                                        :: Eshared !etangetial info
      !only needed by Slots and processed by anisotropic
      type(Shared_t)                                        :: Hshared !hnormal info
      type(XYZlimit_t), dimension(1:6)                      :: Alloc,Sweep,SINPMLSweep
      type(MediaData_t), pointer, dimension(:)            :: Med
      type(NodalSource_t), dimension(:), pointer           :: NodalSource
      type(obses_t)  , pointer, dimension(:)              :: Observation
      !
      logical  :: thereAreMagneticMedia
      logical  :: thereArePMLMagneticMedia
      character(len=BUFSIZE) :: nEntradaRoot
      type(coorsxyzP_t) :: Punto
   end type

   type media_matrices_t
      integer(kind=INTEGERSIZEOFMEDIAMATRICES) , allocatable , dimension(:,:,:) :: sggMiNo,sggMiEx,sggMiEy,sggMiEz,sggMiHx,sggMiHy,sggMiHz
      integer(kind=IKINDMTAG) , allocatable , dimension(:,:,:) :: sggMtag
   end type
      

   type :: constants_t
      real(kind=rkind), pointer, dimension(:) :: g1,g2,gM1,gM2
   contains
      procedure :: destroy => constants_destroy 
   end type


   type nf2ff_t
      logical :: tr,fr,iz,de,ab,ar
   end type

   type :: perform_t
      logical :: flushFields = .false.
      logical :: flushData = .false.
      logical :: unpackFlag = .false.
      logical :: postprocess = .false.
      logical :: flushXdmf = .false.
      logical :: flushVTK = .false.
   contains
      procedure :: isFlush
      procedure :: reset => perform_reset
   end type

   ! variables for timestepping solver control
   type :: sim_control_t
      logical :: simu_devia, resume,saveall,makeholes,& 
                 connectendings,isolategroupgroups,createmap, & 
                 groundwires,noSlantedcrecepelo, & 
                 mibc,ADE,conformalskin,sgbc, sgbcDispersive, sgbccrank, & 
                 NOcompomur,strictOLD,TAPARRABOS, & 
                 noconformalmapvtk, experimentalVideal, &
                 forceresampled, mur_second,MurAfterPML, &
                 stableradholland,singlefilewrite,NF2FFDecim, &
                 fieldtotl,finishedwithsuccess, &
                 permitscaling,mtlnberenger,niapapostprocess, &
                 stochastic, verbose, dontwritevtk, &
                 resume_fromold, vtkindex,createh5bin,wirecrank,fatalerror
      real(kind=8) :: time_desdelanzamiento
      real(kind=RKIND) :: cfl, attfactorc,attfactorw, alphamaxpar, &
                           alphaOrden, kappamaxpar, mindistwires,sgbcFreq,sgbcresol, maxSourceValue
      real(kind=RKIND_WIRES) :: factorradius,factordelta
      
      character(len=BUFSIZE) :: nEntradaRoot, inductance_model,wiresflavor, nresumeable2
      character(len=BUFSIZE) :: opcionestotales
      
      integer(kind=4) :: finaltimestep, flushsecondsFields,flushsecondsData, layoutnumber,& 
                          mpidir, inductance_order, wirethickness, maxCPUtime, SGBCDepth, precision, num_procs
      
      type(ExtraMedium_t) :: extraMedium
      type(nf2ff_T) :: facesNF2FF

   end type sim_control_t

   !!!!!!!!VARIABLES GLOBALES
   integer(kind=4), save, public :: prior_BV     , &
   prior_IB     , &
   prior_pmlbody, &
   prior_AB     , &
   prior_FDB    , &
   prior_IS     , &
   prior_AS     , &
   prior_FDS    , &
   prior_IL     , &
   prior_AL     , &
   prior_FDL    , &
   prior_IP     , &
   prior_AP     , &
   prior_FDP    , &
   prior_PEC    , &
   prior_PMC    , &
   prior_TG     , &
   prior_CS     , &
   prior_TW

   !**************************************************************************************************
   !**************************************************************************************************
   !conformal existence flags   ref: ##Confflag##
   logical, save, public  :: input_conformal_flag
   !**************************************************************************************************
   !**************************************************************************************************

contains

   subroutine constants_destroy(this)
      class(constants_t) :: this
      deallocate(this%g1,this%g2,this%gm1,this%gm2)
   end subroutine

   logical function isFlush(this)
      class(perform_t) :: this
      isFlush = this%flushDATA.or.this%flushFIELDS.or.this%postprocess.or.this%flushXdmf.or.this%flushVTK
   end function

   subroutine perform_reset(this)
      class(perform_t) :: this
      this%flushFields = .false.
      this%flushData = .false.
      this%unpackFlag = .false.
      this%postprocess = .false.
      this%flushXdmf = .false.
      this%flushVTK = .false.
   end subroutine 

   subroutine logic_reset(this)
      class(logic_control_t) :: this
      this%Wires = .false.
      this%PMLbodies = .false.
      this%MultiportS = .false.
      this%AnisMultiportS = .false.
      this%SGBCs= .false.
      this%Lumpeds= .false.
      this%EDispersives = .false.
      this%MDispersives = .false.
      this%PlaneWaveBoxes = .false.
      this%Observation = .false.
      this%FarFields = .false.
      this%PMCBorders = .false.
      this%PMLBorders = .false.
      this%MurBorders = .false.
      this%PECBorders = .false.
      this%Anisotropic = .false.
      this%ThinSlot = .false.
      this%NodalE = .false.
      this%NodalH = .false.
      this%PeriodicBorders = .false.
      this%MagneticMedia = .false.
      this%PMLMagneticMedia= .false.
      this%MTLNbundles = .false.
   end subroutine 

   subroutine setglobal(iu1,iu2)
       integer(kind=4) :: iu1,iu2
       quienmpi=iu1
       tamaniompi=iu2
       return
   end subroutine

   subroutine set_priorities(prioritizeCOMPOoverPEC,prioritizeISOTROPICBODYoverall,prioritizeTHINWIRE)
      logical :: prioritizeCOMPOoverPEC,prioritizeISOTROPICBODYoverall,prioritizeTHINWIRE
      !!movido aqui el sistema de prioridades para poder controlarlos con switches. util para siva 070815 (bug de PEC con prioridad sobre compo del siva
      prior_BV      =10 !background volume
      prior_AB      =30 !anisotropic body
      prior_FDB     =40 !Frequency dependent body
      prior_IS      =50 !Isotropic surface
      prior_AS      =60   !Anisotropic surface
      prior_FDS     =70   !Frequency dependent surface
      prior_IL      =90   !Isotropic line
      prior_AL      =100  !Anisotropic line
      prior_FDL     =110  !Frequency dependent line
      prior_IP      =120  !Isotropic point
      prior_AP      =130  !Anisotropic point
      prior_FDP     =140  !Frequency dependent point
      prior_PEC     =150  !Perfectly electric conducting body, surface, line, or point
      prior_PMC     =160  !Perfectly magnetic conducting body, surface, line, or point
      prior_TG      =155       !thin Slot has more priority than PEC
      !!!!!se aniade la opcion -prioritizeCOMPOoverPEC para subir su prioridad y poder simular SIVA (sgg 070815)
      if (prioritizeTHINWIRE) then
        prior_TW   = 1500   !cambiado a 231024 y puesto con maxima prioridad. es solo experimental y por visualizacion    
      else !opcion correcta. lo anterior es solo experimental y por visualizacion      
        prior_TW   = 15   !prioridad del thin-wire por debajo de todos (excepto del background)  
      end if  
!      prior_pmlbody = prior_TW-1 !el hilo tiene prioridad sobre el pmlbody (prueba HOLD coax sgg 251019)
      prior_pmlbody = prior_BV+1 !el pml body puede ser penetrado por todo 311019 sgg
      !!!!
      if (prioritizeCOMPOoverPEC) then  !Composite surface
         prior_CS=prior_PEC+2
      else
         prior_CS=prior_PEC-2       !composites has lower than PEC to properly handle junctions PEC-composite !(ss's 210312 mail)
      end if
      if (prioritizeISOTROPICBODYoverall) then  ! Isotropic body
         prior_IB      = 200   !SOLO PARA EL CASO DEL SIVA SACAR BOCADOS DE vacio 
      else
         prior_IB      =   20 !EL SUSUAL
      end if 
      return
      

   end subroutine set_priorities
   
   function taglist_getFaceTag(this, field, i, j, k) result(res)
      class(taglist_t) :: this
      integer(kind=IKINDMTAG) :: res 
      integer(kind = 4) :: field, i, j, k
      select case(field)
      case(IHX)
         res = this%face%x(i, j, k)
      case(IHY)
         res = this%face%y(i, j, k)
      case(IHZ)
         res = this%face%z(i, j, k)
      end select
   end function

   function taglist_getEdgeTag(this, field, i, j, k) result(res)
      class(taglist_t) :: this
      integer(kind=IKINDMTAG) :: res 
      integer(kind = 4) :: field, i, j, k
      select case(field)
      case(iEx)
         res = this%edge%x(i, j, k)
      case(iEy)
         res = this%edge%y(i, j, k)
      case(IEZ)
         res = this%edge%z(i, j, k)
      end select
   end function

   logical function direction_eq(a,b)
      class(direction_t), intent(in) :: a,b 
      direction_eq = .true.
      direction_eq = direction_eq .and. (a%x == b%x)
      direction_eq = direction_eq .and. (a%y == b%y)
      direction_eq = direction_eq .and. (a%z == b%z)
      direction_eq = direction_eq .and. (a%orientation == b%orientation)

   end function
end module FDETYPES_m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         STRUCTURE OF SGG
! NumMedia                       : Number of different media
! NumPlaneWaves                  : Number of Plane Wave Boxes (only one currently allowed)
! TimeSteps                      : Number of simulation time steps
! InitialTimeStep                : Initial time step (1 to start, otherwise read from file to resume a previous simulation)
!
! LineX( : ),LineY( : ),LineZ( : )     : Positions of the lines containing Ex, Ey, Ez ('discretization lines')
! PlaneWave                      : Plane Wave info
!       px,py,pz                              : components of the incidence vector
!       ex,ey,ez                              : amplitudes of the electric field (must be perpendicular to de incident vector)
!       esqx1,esqy1,esqz1,esqx2,esqy2,esqz2   : discretization lines bounding the Huygens surface
!       fichero                               : name of the field with the time profile of the transinet excitation
!                                                     (must be well sampled).
! Border                        : Limits of the compuational domain info
!       IsBackPEC,IsFrontPEC,IsLeftPEC,IsRightPEC,IsUpPEC,IsDownPEC     : Whether each limit is PEC, PMC or PML
!       IsBackPMC,IsFrontPMC,IsLeftPMC,IsRightPMC,IsUpPMC,IsDownPMC
!       IsBackPML,IsFrontPML,IsLeftPML,IsRightPML,IsUpPML,IsDownPML
! PML                           : PML info (meningless if Is...PEC or Is...PMC are set
!       CoeffReflPML(3,2)                           : Refflection coeffients at the end of the PML at each
!                                                     termination ({1-x,2-y,3-z} : {1-start,2-end})
!       NumLayers(3,2)                              : Number of PML layers ({1-x,2-y,3-z} : {1-start,2-end})
!
! M(1:6)                        : six matrices (one per field component) with the index of the medium present at each Yee location
!                               : plus 1 for the centroid of the cell
!       XI,XE,YI,YE,ZI,ZE                   : Mediamatrix dimensions
!       Mediamatrix( : , : , :                   : Index of the medium present at each Yee location
!
! Med                           : Info of each medium
!      Epr(:),Sigma(:),Mur(:),SigmaM(:)   : Relative permittivity, electric conductivity,
!                                           Relative permeability, Magnetic Conductivity
!      Priority(:)                   : To decide overlapping of media (meningless in the simulation, only needed during PREPROCESS)
!      IsPML( : )                      : If the medium is a PML (needed to calculate the especific PML updating coefficients)
!      Wire( : )                       : If the medium is a wire, this type contains its parameters
!
!           TipoWire     : Info on the wire parameters
!                      radius,R,L       : radius, resistance per unit length, inductance per unit length
!                      Vsource,Isource  : Info with the voltage/current source on the wire
!                              Exists          : Wheter this wire is a source (a single wire, with a single segment)
!                                                is needed for the source
!                              Fichero         : name of the field with the time profile of the transinet excitation
!                                                (must be well sampled).
!      Multiport( : )                          : If the medium is a Multiport, this type contains its parameters
!           multiportFileZ11,multiportFileZ22,multiportFileZ12 : Files with the pole/residues info
! Observation : Observation info
!           Size     : How many observation points
!           XI,YI,ZI : index of the voxel to be observed (a voxel is limited by 8 discretization lines. An average to find
!                      the magnitude at the center is used)
!           What     : What to observe, fields (Ex, Ey, Ez, Hx, Hy, Hz) or currents at wires (Jx,Jy,Jz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
