module NFDETypes
   !
   use FDETYPES
#ifdef CompileWithMTLN   
   use mtln_types_mod
#endif
   use conformal_types_mod
   !
   implicit none
   integer(kind=4), parameter :: RK = RKIND
   !------------------------------------------------------------------------------
   ! CONSTANTS FOR THE PARSER
   !------------------------------------------------------------------------------
   ! global variable stochastic
   ! MATERIALS
   real(kind=RK), parameter :: SIGMA_PEC = 1e19_RK
   real(kind=RK), parameter :: SIGMA_PMC = 1e19_RK 
   ! PROBES
   !!!!
   integer(kind=4), parameter :: NP_T1_PLAIN = 0
   integer(kind=4), parameter :: NP_T1_AMBOS = 2
   integer(kind=4), parameter :: NP_T2_TIME = 0
   integer(kind=4), parameter :: NP_T2_FREQ = 1
   integer(kind=4), parameter :: NP_T2_TRANSFER = 2
   integer(kind=4), parameter :: NP_T2_TIMEFREQ  = 3
   integer(kind=4), parameter :: NP_T2_TIMETRANSF = 4
   integer(kind=4), parameter :: NP_T2_FREQTRANSF = 5
   integer(kind=4), parameter :: NP_T2_TIMEFRECTRANSF = 6
   integer(kind=4), parameter :: NP_COR_EX = 0
   integer(kind=4), parameter :: NP_COR_EY = 1
   integer(kind=4), parameter :: NP_COR_EZ = 2
   integer(kind=4), parameter :: NP_COR_HX = 3
   integer(kind=4), parameter :: NP_COR_HY = 4
   integer(kind=4), parameter :: NP_COR_HZ = 5
   integer(kind=4), parameter :: NP_COR_WIRECURRENT = 6
   integer(kind=4), parameter :: NP_COR_DDP = 7
   integer(kind=4), parameter :: NP_COR_LINE = 8
   integer(kind=4), parameter :: NP_COR_CHARGE = 9
   LOGICAL, parameter :: BcELECT = .TRUE.
   LOGICAL, parameter :: BcMAGNE = .FALSE.
   ! THIN WIRES
   integer(kind=4), parameter :: MATERIAL_CONS = 0
   integer(kind=4), parameter :: MATERIAL_absorbing = 100
   integer(kind=4), parameter :: Parallel_CONS = 1
   integer(kind=4), parameter :: SERIES_CONS = 2
   integer(kind=4), parameter :: DISPERSIVE_CONS = 3
   ! BORDERS
   integer(kind=4), parameter :: F_PEC = 1
   integer(kind=4), parameter :: F_PMC = 2
   integer(kind=4), parameter :: F_PER = 4
   integer(kind=4), parameter :: F_MUR = 7
   integer(kind=4), parameter :: F_PML = 9
   integer(kind=4), parameter :: F_XL = 1
   integer(kind=4), parameter :: F_XU = 2
   integer(kind=4), parameter :: F_YL = 3
   integer(kind=4), parameter :: F_YU = 4
   integer(kind=4), parameter :: F_ZL = 5
   integer(kind=4), parameter :: F_ZU = 6
   integer(kind=4), parameter :: F_TIMEFRECTRANSF = 0
   ! rlc y diodos
   integer(kind=4), parameter :: inductor = 20
   integer(kind=4), parameter :: capacitor = 21
   integer(kind=4), parameter :: resistor = 22
   integer(kind=4), parameter :: diodo = 23
   integer(kind=4), parameter :: Dielectric = 24
   integer(kind=4), parameter :: PMLbody = 25

   !------------------------------------------------------------------------------
   ! TYPES
   !------------------------------------------------------------------------------
   !-----------------> Cordinate Types
   !------------------------------------------------------------------------------
   ! Basic cordinate type for two points and orientation
   !------------------------------------------------------------------------------
   type, public :: coords
      integer(kind=4) :: Xi = - 1
      integer(kind=4) :: Xe = - 1
      integer(kind=4) :: Yi = - 1
      integer(kind=4) :: Ye = - 1
      integer(kind=4) :: Zi = - 1
      integer(kind=4) :: Ze = - 1
      integer(kind=4) :: Xtrancos = 1
      integer(kind=4) :: Ytrancos = 1
      integer(kind=4) :: Ztrancos = 1
      integer(kind=4) :: Or = 0 !f1eld orientation
      character(len=BUFSIZE) :: tag
   END type coords
   type, public :: coords_scaled
      integer(kind=4) :: Xi = - 1
      integer(kind=4) :: Xe = - 1
      integer(kind=4) :: Yi = - 1
      integer(kind=4) :: Ye = - 1
      integer(kind=4) :: Zi = - 1
      integer(kind=4) :: Ze = - 1
      real(kind=RK) :: xc = 0.0_RKIND
      real(kind=RK) :: yc = 0.0_RKIND
      real(kind=RK) :: zc = 0.0_RKIND
      integer(kind=4) :: Or = 0 !field orientation nuevo 2015
      character(len=BUFSIZE) :: tag
   END type coords_scaled
   !-----------------> Material Types
   !------------------------------------------------------------------------------
   ! Basic constants for materials
   !------------------------------------------------------------------------------
   type, public :: Material
      real(kind=RK) :: eps = 0.0_RKIND
      real(kind=RK) :: mu = 0.0_RKIND
      real(kind=RK) :: sigma = 0.0_RKIND
      real(kind=RK) :: sigmam = 0.0_RKIND
      integer(kind=4) :: id = 0
   END type Material
   !------------------------------------------------------------------------------
   ! New Class which is a collection of different materials
   !------------------------------------------------------------------------------
   type, public :: Materials
      integer(kind=4) :: n_Mats = 0
      integer(kind=4) :: n_Mats_max = 0
      type(Material), dimension(:), pointer :: Mats => NULL ()
   END type Materials

   !------------------------------------------------------------------------------
   ! Identifies conformal PEC "media"
   !------------------------------------------------------------------------------

   type, public :: ConformalPECElements
      type(triangle_t), dimension(:), allocatable :: triangles
      type(interval_t), dimension(:), allocatable :: intervals
      character(len=bufsize) :: tag
   end type 

   type, public :: ConformalPECRegions
      type(ConformalPECElements), dimension(:), pointer :: volumes => null()
      type(ConformalPECElements), dimension(:), pointer :: surfaces => null()
   end type


   type, public :: edge_t 
      integer(kind=4), dimension(3) :: cell
      integer(kind=4) :: direction = -1
      real(kind=rkind) :: ratio = -1
      real(kind=rkind), dimension(2) :: material_coords
   end type 
   type, public :: face_t 
      integer(kind=4), dimension(3) :: cell
      integer(kind=4) :: direction = -1
      real(kind=rkind) :: ratio = -1
   end type 

   type, public :: conformal_edge_media_t
      type(edge_t), dimension(:), allocatable :: edges
      real(kind=rkind) :: ratio
      integer(kind=4) :: size
   end type
   type, public :: conformal_face_media_t
      type(face_t), dimension(:), allocatable :: faces
      real(kind=rkind) :: ratio
      integer(kind=4) :: size
   end type

   type, public :: ConformalMedia_t
      integer(kind=4) :: n_edges_media = 0
      integer(kind=4) :: n_faces_media = 0
      type(conformal_face_media_t), dimension(:), pointer :: face_media => NULL ()
      type(conformal_edge_media_t), dimension(:), pointer :: edge_media => NULL ()
      real(kind=rkind) :: time_step_scale_factor = 1.0
      character(len=bufsize) :: tag
   END type ConformalMedia_t


   !------------------------------------------------------------------------------
   ! Locates all the different PEC media found
   !------------------------------------------------------------------------------

   type, public :: PECRegions
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nSurfs = 0
      integer(kind=4) :: nLins = 0
      integer(kind=4) :: nVols_max = 0
      integer(kind=4) :: nSurfs_max = 0
      integer(kind=4) :: nLins_max = 0
      type(coords), dimension(:), pointer :: Vols => NULL ()
      type(coords), dimension(:), pointer :: Surfs => NULL ()
      type(coords), dimension(:), pointer :: Lins => NULL ()
   END type PECRegions
   !------------------------------------------------------------------------------
   ! Defines a Non Metal Body
   !------------------------------------------------------------------------------
   type, public :: Dielectric_t
      type(coords), dimension(:), pointer :: c1P => NULL ()
      type(coords), dimension(:), pointer :: c2P => NULL ()
      real(kind=RK) :: sigma = 0.0_RKIND
      real(kind=RK) :: eps = 0.0_RKIND
      real(kind=RK) :: mu = 0.0_RKIND
      real(kind=RK) :: sigmam = 0.0_RKIND
      integer(kind=4) :: n_C1P = 0
      integer(kind=4) :: n_C2P = 0
      !!! Lumped vars.
      real(kind=RK) :: Rtime_on = 0.0_RKIND, Rtime_off = 0.0_RKIND
      real(kind=RK) :: R = 0.0_RKIND
      real(kind=RK) :: L = 0.0_RKIND
      real(kind=RK) :: C = 0.0_RKIND
      !!! Stoch vars.
      real(kind=RK) :: R_devia = 0.0_RKIND
      real(kind=RK) :: L_devia = 0.0_RKIND
      real(kind=RK) :: C_devia = 0.0_RKIND
      !
      real(kind=RK) :: DiodB = 0.0_RKIND
      real(kind=RK) :: DiodIsat = 0.0_RKIND
      integer(kind=4) :: DiodOri = 0
      !!! Berenger's waveports
      integer(kind=4) :: orient = 0
!!!!!!!!!
      logical :: resistor=.false. , inductor=.false. , capacitor=.false. , diodo=.false. , plain=.false. , PMLbody=.false.
   END type Dielectric_t
   !------------------------------------------------------------------------------
   ! Locates all the different Non Metal Media found
   !------------------------------------------------------------------------------
   type, public :: DielectricRegions
      type(Dielectric_t), dimension(:), pointer :: Vols => NULL ()
      type(Dielectric_t), dimension(:), pointer :: Surfs => NULL ()
      type(Dielectric_t), dimension(:), pointer :: Lins => NULL ()
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nSurfs = 0
      integer(kind=4) :: nLins = 0
      integer(kind=4) :: nVols_max = 0
      integer(kind=4) :: nSurfs_max = 0
      integer(kind=4) :: nLins_max = 0
      integer(kind=4) :: n_C1P_max = 0
      integer(kind=4) :: n_C2P_max = 0
   END type DielectricRegions
   !------------------------------------------------------------------------------
   ! type that defines the information of a frequency depENDent material,
   ! it inherits from the material class and it adds the possible values needed
   ! in the frequency depENDent section of the
   !------------------------------------------------------------------------------
   type, public :: FreqDepenMaterial
      complex, dimension(:), pointer :: a11 => NULL ()
      complex, dimension(:), pointer :: b11 => NULL ()
      complex, dimension(:), pointer :: am11 => NULL ()
      complex, dimension(:), pointer :: bm11 => NULL ()
      complex, dimension(:), pointer :: a12 => NULL ()
      complex, dimension(:), pointer :: b12 => NULL ()
      complex, dimension(:), pointer :: am12 => NULL ()
      complex, dimension(:), pointer :: bm12 => NULL ()
      complex, dimension(:), pointer :: a13 => NULL ()
      complex, dimension(:), pointer :: b13 => NULL ()
      complex, dimension(:), pointer :: am13 => NULL ()
      complex, dimension(:), pointer :: bm13 => NULL ()
      complex, dimension(:), pointer :: a22 => NULL ()
      complex, dimension(:), pointer :: b22 => NULL ()
      complex, dimension(:), pointer :: am22 => NULL ()
      complex, dimension(:), pointer :: bm22 => NULL ()
      complex, dimension(:), pointer :: a23 => NULL ()
      complex, dimension(:), pointer :: b23 => NULL ()
      complex, dimension(:), pointer :: am23 => NULL ()
      complex, dimension(:), pointer :: bm23 => NULL ()
      complex, dimension(:), pointer :: a33 => NULL ()
      complex, dimension(:), pointer :: b33 => NULL ()
      complex, dimension(:), pointer :: am33 => NULL ()
      complex, dimension(:), pointer :: bm33 => NULL ()
      real(kind=RK), dimension(:), pointer :: alpha => NULL ()
      real(kind=RK), dimension(:), pointer :: beta => NULL ()
      real(kind=RK), dimension(:), pointer :: gamma => NULL ()
      real(kind=RK), dimension(:), pointer :: alpham => NULL ()
      real(kind=RK), dimension(:), pointer :: betam => NULL ()
      real(kind=RK), dimension(:), pointer :: gammam => NULL ()
      type(coords), dimension(:), pointer :: c => NULL ()
      real(kind=RK) ::  eps11 = 0.0_RKIND ,    eps12 = 0.0_RKIND ,    eps13 = 0.0_RKIND ,    eps22 = 0.0_RKIND ,    eps23 = 0.0_RKIND ,    eps33 = 0.0_RKIND
      real(kind=RK) ::   mu11 = 0.0_RKIND ,     mu12 = 0.0_RKIND ,     mu13 = 0.0_RKIND ,     mu22 = 0.0_RKIND ,     mu23 = 0.0_RKIND ,     mu33 = 0.0_RKIND
      real(kind=RK) :: sigma11 = 0.0_RKIND ,  sigma12 = 0.0_RKIND ,  sigma13 = 0.0_RKIND ,  sigma22 = 0.0_RKIND ,  sigma23 = 0.0_RKIND ,  sigma33 = 0.0_RKIND
      real(kind=RK) :: sigmam11 = 0.0_RKIND , sigmam12 = 0.0_RKIND , sigmam13 = 0.0_RKIND , sigmam22 = 0.0_RKIND , sigmam23 = 0.0_RKIND , sigmam33 = 0.0_RKIND
      integer(kind=4) :: K11 = 0
      integer(kind=4) :: Km11 = 0
      integer(kind=4) :: K12 = 0
      integer(kind=4) :: Km12 = 0
      integer(kind=4) :: K13 = 0
      integer(kind=4) :: Km13 = 0
      integer(kind=4) :: K22 = 0
      integer(kind=4) :: Km22 = 0
      integer(kind=4) :: K23 = 0
      integer(kind=4) :: Km23 = 0
      integer(kind=4) :: K33 = 0
      integer(kind=4) :: Km33 = 0
      integer(kind=4) :: L = 0
      integer(kind=4) :: Lm = 0
      integer(kind=4) :: n_c = 0
      character(len=BUFSIZE) :: files = ' ' !2015 si esta presente lee los polos/residuos desde fichero
   END type FreqDepenMaterial
   !------------------------------------------------------------------------------
   ! type that defines the list of frequency depedent materials
   !------------------------------------------------------------------------------
   type, public :: FreqDepenMaterials
      type(FreqDepenMaterial), dimension(:), pointer :: Vols => NULL ()
      type(FreqDepenMaterial), dimension(:), pointer :: Surfs => NULL ()
      type(FreqDepenMaterial), dimension(:), pointer :: Lins => NULL ()
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nSurfs = 0
      integer(kind=4) :: nLins = 0
      integer(kind=4) :: nVols_max = 0
      integer(kind=4) :: nSurfs_max = 0
      integer(kind=4) :: nLins_max = 0
      integer(kind=4) :: n_c_max = 0 !cota superior
   END type FreqDepenMaterials
   !------------------------------------------------------------------------------
   ! Type for the ANISOTROPIC body, surface and lines since they will contain
   ! the same information
   !------------------------------------------------------------------------------
   type, public :: ANISOTROPICbody_t
      type(coords), dimension(:), pointer :: c1P => NULL ()
      type(coords), dimension(:), pointer :: c2P => NULL ()
      real(kind=RK), dimension(3, 3) :: sigma, eps, mu, sigmam
      integer(kind=4) :: n_C1P = 0
      integer(kind=4) :: n_C2P = 0
   END type ANISOTROPICbody_t
   !------------------------------------------------------------------------------
   ! Type that contains the elements found in the nfde File
   !------------------------------------------------------------------------------
   type, public :: ANISOTROPICelements_t
      type(ANISOTROPICbody_t), dimension(:), pointer :: Vols => NULL ()
      type(ANISOTROPICbody_t), dimension(:), pointer :: Surfs => NULL ()
      type(ANISOTROPICbody_t), dimension(:), pointer :: Lins => NULL ()
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nSurfs = 0
      integer(kind=4) :: nLins = 0
      integer(kind=4) :: nVols_max = 0
      integer(kind=4) :: nSurfs_max = 0
      integer(kind=4) :: nLins_max = 0
      integer(kind=4) :: n_C1P_max = 0 !cota superior de c1p y c2p en vols,sufs,lins
      integer(kind=4) :: n_C2P_max = 0
   END type ANISOTROPICelements_t
   !------------------------------------------------------------------------------
   ! Defines a Comp Surface
   !------------------------------------------------------------------------------
   type, public :: LossyThinSurface
      type(coords), dimension(:), pointer :: c => NULL ()
      real(kind=RK), dimension(:), pointer :: sigma
      real(kind=RK), dimension(:), pointer :: eps
      real(kind=RK), dimension(:), pointer :: mu
      real(kind=RK), dimension(:), pointer :: sigmam
      real(kind=RK), dimension(:), pointer :: thk
      !for_devia
      real(kind=RK), dimension(:), pointer :: sigma_devia
      real(kind=RK), dimension(:), pointer :: eps_devia
      real(kind=RK), dimension(:), pointer :: mu_devia
      real(kind=RK), dimension(:), pointer :: sigmam_devia
      real(kind=RK), dimension(:), pointer :: thk_devia
      !
      integer(kind=4) :: nc = 0
      character(len=BUFSIZE) :: files = ' ' 
      integer(kind=4) :: numcapas  
   END type LossyThinSurface
   !------------------------------------------------------------------------------
   ! Locates all the different Comp media found
   !------------------------------------------------------------------------------
   type, public :: LossyThinSurfaces
      type(LossyThinSurface), dimension(:), pointer :: cs => NULL ()
      integer(kind=4) :: length = 0
      integer(kind=4) :: length_max = 0
      integer(kind=4) :: nC_max = 0 !cota de todos los nc de LossyThinSurface
   END type LossyThinSurfaces
   !------------------------------------------------------------------------------
   ! Component for Thin Wires there is a list of this inside the component
   ! that defines the whole Thin Wire Reference
   !------------------------------------------------------------------------------
   type, public :: ThinWireComp
      character(len=BUFSIZE) :: srctype, srcfile
      integer(kind=4) :: i = - 1
      integer(kind=4) :: j = - 1
      integer(kind=4) :: K = - 1
      integer(kind=4) :: nd = - 1
      integer(kind=4) :: d = - 1
      real(kind=RK) :: m = 0.0_RKIND
      character(len=BUFSIZE) :: tag
   END type ThinWireComp
   !------------------------------------------------------------------------------
   ! ThinWire component that defines the overall properties of the definition
   ! of ThinWires
   !------------------------------------------------------------------------------
   type, public :: ThinWire
      type(ThinWireComp), dimension(:), pointer :: twc => NULL ()
      real(kind=RK) :: rad = 0 , rad_devia = 0
      LOGICAL :: disp = .false.
      character(len=BUFSIZE) :: dispfile
      real(kind=RK) :: res = 0 , res_devia = 0
      real(kind=RK) :: ind = 0 , ind_devia = 0
      real(kind=RK) :: cap = 0 , cap_devia = 0
      real(kind=RK) :: P_res = 0
      real(kind=RK) :: P_ind = 0
      real(kind=RK) :: P_cap = 0
      character(len=BUFSIZE) :: dispfile_LeftEnd
      real(kind=RK) :: R_LeftEnd = 0 , R_LeftEnd_devia = 0
      real(kind=RK) :: L_LeftEnd = 0 , L_LeftEnd_devia = 0
      real(kind=RK) :: C_LeftEnd = 0 , C_LeftEnd_devia = 0
      character(len=BUFSIZE) :: dispfile_RightEnd
      real(kind=RK) :: R_RightEnd = 0 , R_RightEnd_devia = 0
      real(kind=RK) :: L_RightEnd = 0 , L_RightEnd_devia = 0
      real(kind=RK) :: C_RightEnd = 0 , C_RightEnd_devia = 0
      integer(kind=4) :: LeftEnd = 0
      integer(kind=4) :: RightEnd = 0
      ! Components
      integer(kind=4) :: tl = 0
      integer(kind=4) :: tr = 0
      integer(kind=4) :: n_twc = 0
      integer(kind=4) :: n_twc_max = 0
   END type ThinWire
   !------------------------------------------------------------------------------
   ! List of the different thin wires that were found in the file
   !------------------------------------------------------------------------------
   type, public :: ThinWires
      type(ThinWire), dimension(:), pointer :: tw => NULL ()
      integer(kind=4) :: n_tw = 0
      integer(kind=4) :: n_tw_max = 0
   END type ThinWires
   !------------------------------------------------------------------------------
   ! Component for Slanted Wires there is a list of this inside the component
   ! that defines the whole Slanted Wire Reference
   !------------------------------------------------------------------------------
   type, public :: SlantedWireComp
      character(len=BUFSIZE) :: srctype, srcfile
      real(kind=RK) :: x = - 1.0_RKIND
      real(kind=RK) :: y = - 1.0_RKIND
      real(kind=RK) :: z = - 1.0_RKIND
      integer(kind=4) :: nd = - 1
      real(kind=RK) :: m = 0.0_RKIND
      character(len=BUFSIZE) :: tag
   END type SlantedWireComp
   !------------------------------------------------------------------------------
   ! ThinWire component that defines the overall properties of the definition
   ! of ThinWires
   !------------------------------------------------------------------------------
   type, public :: SlantedWire
      type(SlantedWireComp), dimension(:), pointer :: swc => NULL ()
      real(kind=RK) :: rad = 0
      LOGICAL :: disp = .false.
      character(len=BUFSIZE) :: dispfile
      real(kind=RK) :: res = 0
      real(kind=RK) :: ind = 0
      real(kind=RK) :: cap = 0
      real(kind=RK) :: P_res = 0
      real(kind=RK) :: P_ind = 0
      real(kind=RK) :: P_cap = 0
      character(len=BUFSIZE) :: dispfile_LeftEnd
      real(kind=RK) :: R_LeftEnd = 0
      real(kind=RK) :: L_LeftEnd = 0
      real(kind=RK) :: C_LeftEnd = 0
      character(len=BUFSIZE) :: dispfile_RightEnd
      real(kind=RK) :: R_RightEnd = 0
      real(kind=RK) :: L_RightEnd = 0
      real(kind=RK) :: C_RightEnd = 0
      integer(kind=4) :: LeftEnd = 0
      integer(kind=4) :: RightEnd = 0
      ! Components
      integer(kind=4) :: tl = 0
      integer(kind=4) :: tr = 0
      integer(kind=4) :: n_swc = 0
      integer(kind=4) :: n_swc_max = 0
   END type SlantedWire
   !------------------------------------------------------------------------------
   ! List of the different thin wires that were found in the file
   !------------------------------------------------------------------------------
   type, public :: SlantedWires
      type(SlantedWire), dimension(:), pointer :: sw => NULL ()
      integer(kind=4) :: n_sw = 0
      integer(kind=4) :: n_sw_max = 0
   END type SlantedWires
   !--------------------------------------------------------------------------
   ! Component for Thin Slots there is a list of this inside the component
   ! that defines the whole Thin Slot Reference
   !--------------------------------------------------------------------------
   type, public :: ThinSlotComp
      integer(kind=4) :: i = 0
      integer(kind=4) :: j = 0
      integer(kind=4) :: K = 0
      integer(kind=4) :: node = 0
      integer(kind=4) :: dir = - 1
      integer(kind=4) :: Or = - 1
      character(len=BUFSIZE) :: tag
   END type ThinSlotComp
   !--------------------------------------------------------------------------
   ! ThinSlot component that defines the overall properties of the definition
   ! of ThinSlots in ORIGINAL
   !--------------------------------------------------------------------------
   type, public :: ThinSlot
      type(ThinSlotComp), dimension(:), pointer :: tgc => NULL ()
      real(kind=RK) :: width = 0
      integer(kind=4) :: n_tgc = 0
      integer(kind=4) :: n_tgc_max = 0
   END type ThinSlot
   !--------------------------------------------------------------------------
   ! List of the different thin Slots that were found in the file
   !--------------------------------------------------------------------------
   type, public :: ThinSlots
      type(ThinSlot), dimension(:), pointer :: tg => NULL ()
      integer(kind=4) :: n_tg = 0
      integer(kind=4) :: n_tg_max = 0
   END type ThinSlots
   !-----------------> Border Types
   !------------------------------------------------------------------------------
   ! PML Border Type
   !------------------------------------------------------------------------------
   type, public :: FronteraPML
      real(kind=RK) :: orden = 2.0_RK
      real(kind=RK) :: refl = 1e-3_RK
      integer(kind=4) :: numCapas = 8
   END type FronteraPML
   !------------------------------------------------------------------------------
   ! Tipo de la frontera
   !------------------------------------------------------------------------------
   type, public :: Frontera
      integer(kind=4), dimension(6) :: tipoFrontera
      type(FronteraPML), dimension(6) :: propiedadesPML
   END type Frontera
   !-----------------> Probe Types
   !------------------------------------------------------------------------------
   ! type to define the new probe object which contains the type of calculation
   ! the type of analysis, time and frequency step and the filename where
   ! it should be saved
   !------------------------------------------------------------------------------
   type, public :: MasSonda
      character(len=BUFSIZE) :: filename
      type(coords), dimension(:), pointer :: cordinates => NULL ()
      real(kind=RK) :: tstart, tstop, tstep
      real(kind=RK) :: fstart, fstop, fstep
      integer(kind=4) :: type1, type2
      integer(kind=4) :: len_cor = 0
      character(len=BUFSIZE) :: outputrequest
   END type MasSonda
   !------------------------------------------------------------------------------
   ! type that defines a list of probes to be appended and accesed
   !------------------------------------------------------------------------------
   type, public :: MasSondas
      type(MasSonda), dimension(:), pointer :: collection => NULL ()
      integer(kind=4) :: length = 0
      integer(kind=4) :: length_max = 0
      integer(kind=4) :: len_cor_max = 0 !cota
   END type MasSondas
   !------------------------------------------------------------------------------
   ! This type contains the basic information in nearly all the different PROBES
   !------------------------------------------------------------------------------
   type, public :: Sonda
      character(len=BUFSIZE) :: grname
      integer(kind=4), dimension(:), pointer :: i => NULL ()
      integer(kind=4), dimension(:), pointer :: j => NULL ()
      integer(kind=4), dimension(:), pointer :: K => NULL ()
      integer(kind=4), dimension(:), pointer :: node => NULL ()
      integer(kind=4) :: n_cord = 0
      integer(kind=4) :: n_cord_max = 0
      real(kind=RK) :: tstart, tstop, tstep
      character(len=BUFSIZE) :: outputrequest
      !por si se precisa para el Far Field
      real(kind=RK) :: fstart, fstop, fstep
      real(kind=RK) :: phistart, phistop, phistep
      real(kind=RK) :: thetastart, thetastop, thetastep
      character(len=BUFSIZE) :: FileNormalize
   END type Sonda
   !------------------------------------------------------------------------------
   ! type for the electric far field
   !------------------------------------------------------------------------------
   type, public :: FarField_Sonda
      type(Sonda) :: probe
   END type FarField_Sonda
   !------------------------------------------------------------------------------
   ! type for the electric field
   !------------------------------------------------------------------------------
   type, public :: Electric_Sonda
      type(Sonda) :: probe
   END type Electric_Sonda
   !------------------------------------------------------------------------------
   ! type for the magnetic field
   !------------------------------------------------------------------------------
   type, public :: Magnetic_Sonda
      type(Sonda) :: probe
   END type Magnetic_Sonda
   !------------------------------------------------------------------------------
   ! type for the normal electric field
   !------------------------------------------------------------------------------
   type, public :: NormalElectric_Sonda
      type(Sonda) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   END type NormalElectric_Sonda
   !------------------------------------------------------------------------------
   ! type for the normal magnetic field
   !------------------------------------------------------------------------------
   type, public :: NormalMagnetic_Sonda
      type(Sonda) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   END type NormalMagnetic_Sonda
   !------------------------------------------------------------------------------
   ! type for the electric surface current density
   !------------------------------------------------------------------------------
   type, public :: SurfaceElectricCurrent_Sonda
      type(Sonda) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   END type SurfaceElectricCurrent_Sonda
   !------------------------------------------------------------------------------
   ! type for the magnetic surface current density
   !------------------------------------------------------------------------------
   type, public :: SurfaceMagneticCurrent_Sonda
      type(Sonda) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   END type SurfaceMagneticCurrent_Sonda
   !------------------------------------------------------------------------------
   ! Abstract class which performs the dynamic dispatching
   !------------------------------------------------------------------------------
   type, public :: abstractSonda
      integer(kind=4) :: n_FarField = 0 
      integer(kind=4) :: n_Electric = 0
      integer(kind=4) :: n_Magnetic = 0
      integer(kind=4) :: n_NormalElectric = 0
      integer(kind=4) :: n_NormalMagnetic = 0
      integer(kind=4) :: n_SurfaceElectricCurrent = 0
      integer(kind=4) :: n_SurfaceMagneticCurrent = 0
      !
      integer(kind=4) :: n_FarField_max = 0
      integer(kind=4) :: n_Electric_max = 0
      integer(kind=4) :: n_Magnetic_max = 0
      integer(kind=4) :: n_NormalElectric_max = 0
      integer(kind=4) :: n_NormalMagnetic_max = 0
      integer(kind=4) :: n_SurfaceElectricCurrent_max = 0
      integer(kind=4) :: n_SurfaceMagneticCurrent_max = 0
      type(FarField_Sonda), dimension(:), pointer :: FarField => NULL ()
      type(Electric_Sonda), dimension(:), pointer :: Electric => NULL ()
      type(Magnetic_Sonda), dimension(:), pointer :: Magnetic => NULL ()
      type(NormalElectric_Sonda), dimension(:), pointer :: NormalElectric => NULL ()
      type(NormalMagnetic_Sonda), dimension(:), pointer :: NormalMagnetic => NULL ()
      type(SurfaceElectricCurrent_Sonda), dimension(:), pointer :: SurfaceElectricCurrent => NULL ()
      type(SurfaceMagneticCurrent_Sonda), dimension(:), pointer :: SurfaceMagneticCurrent => NULL ()
   END type abstractSonda
   !------------------------------------------------------------------------------
   ! Class to account as a list for all the probes
   ! that might be required during the parsing process
   !------------------------------------------------------------------------------
   type, public :: Sondas
      type(abstractSonda), dimension(:), pointer :: probes => NULL ()
      integer(kind=4) :: n_probes = 0
      integer(kind=4) :: n_probes_max = 0
   END type Sondas
   !------------------------------------------------------------------------------
   ! Object type defined for the Bloque current probe
   !------------------------------------------------------------------------------
   type, public :: BloqueProbe
      real(kind=RK) :: tstart, tstop, tstep
      real(kind=RK) :: fstart, fstop, fstep
      character(len=BUFSIZE) :: FileNormalize
      integer(kind=4) :: type2
      integer(kind=4) :: i1, i2, j1, j2, k1, k2, skip
      integer(kind=4) :: nml
      LOGICAL :: t
      character(len=BUFSIZE) :: outputrequest
      character(len=BUFSIZE) :: tag
   END type BloqueProbe
   ! Object made for the collection of defined Bloque probes
   type, public :: BloqueProbes
      type(BloqueProbe), dimension(:), pointer :: bp => NULL ()
      integer(kind=4) :: n_bp = 0
      integer(kind=4) :: n_bp_max = 0
   END type BloqueProbes

   !------------------------------------------------------------------------------
   ! Object type defined for the Volumic probes
   !------------------------------------------------------------------------------
   type, public :: VolProbe
      type(coords), dimension(:), pointer :: cordinates => NULL ()
      real(kind=RK) :: tstart, tstop, tstep
      character(len=BUFSIZE) :: outputrequest
      integer(kind=4) :: len_cor = 0
      !para freq domain
      real(kind=RK) :: fstart, fstop, fstep
      integer(kind=4) :: type2
      character(len=BUFSIZE) :: filename
   END type VolProbe
   ! Object made for the collection of defined Volumic probes
   type, public :: VolProbes
      type(VolProbe), dimension(:), pointer :: collection => NULL ()
      integer(kind=4) :: length = 0
      integer(kind=4) :: length_max = 0
      integer(kind=4) :: len_cor_max = 0 !cota
   END type VolProbes

   !-----------------> Source Types
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: Box
      character(len=BUFSIZE) :: nombre_fichero
      integer(kind=4), dimension(3) :: coor1, coor2
   END type Box
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: Boxes
      type(Box), dimension(:), pointer :: Vols => NULL ()
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nVols_max = 0
   END type Boxes
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: PlaneWave
      character(len=BUFSIZE) :: nombre_fichero
      character(len=BUFSIZE) :: atributo
      integer(kind=4), dimension(3) :: coor1, coor2
      real(kind=RK) :: theta, phi, alpha, beta
      logical :: isRC !for reververation chambers
      real(kind=RK) :: INCERTMAX
      integer(kind=4) :: numModes !for reververation chambers
   END type PlaneWave
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: PlaneWaves
      type(PlaneWave), dimension(:), pointer :: collection => NULL ()
      integer(kind=4) :: nc = 0
      integer(kind=4) :: nC_max = 0
   END type PlaneWaves
   !------------------------------------------------------------------------------
   ! Definicin de los tipos current density que existirn en el ficero
   ! nfde
   !------------------------------------------------------------------------------
   type, public :: Curr_Field_Src
      type(coords_scaled), dimension(:), pointer :: c1P => NULL ()
      type(coords_scaled), dimension(:), pointer :: c2P => NULL ()
      character(len=BUFSIZE) :: nombre
      integer(kind=4) :: n_C1P = 0
      integer(kind=4) :: n_C2P = 0
      LOGICAL :: isElec, isHard, isInitialValue
   END type Curr_Field_Src
   !------------------------------------------------------------------------------
   ! Definicin de las Nodal Source global
   !------------------------------------------------------------------------------
   type, public :: NodSource
      type(Curr_Field_Src), dimension(:), pointer :: NodalSource => NULL ()
      integer(kind=4) :: n_nodSrc = 0
      integer(kind=4) :: n_nodSrc_max = 0
      integer(kind=4) :: n_C1P_max = 0
      integer(kind=4) :: n_C2P_max = 0
   END type NodSource
   !-----------------> General Types
   !------------------------------------------------------------------------------
   ! Matrix attributes.
   ! Total[XYZ] -> Is the cell number for each axis.
   !------------------------------------------------------------------------------
   type, public :: MatrizMedios
      integer(kind=4) :: totalX, totalY, totalZ
   END type MatrizMedios
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type NFDEGeneral
      real(kind=RK) :: dt
      integer(kind=4) :: nmax
      LOGICAL :: mtlnProblem
   END type NFDEGeneral
   !------------------------------------------------------------------------------
   ! Definition of the type. Three vectors are defined, for each axis X,Y,Z. If
   ! their size is equal to 1 then there is a constant increment. If it is not
   ! then it will be one for each cell position.
   !! WARNING
   !! Even though, Type MatrizMedios defines the total number of cell for each
   !! axis. n[XYZ] defines it too partially. Meaning that if there is a
   !! constant increment, the pointer will be a scalar. However, when it is
   !! variable the pointer will have the same size as total[XYZ] in the
   !! MatrizMedios type and for each vector position the increment for those
   !! Cells
   !------------------------------------------------------------------------------
   type Desplazamiento
      real(kind=RK), dimension(:), pointer :: desX => NULL ()
      real(kind=RK), dimension(:), pointer :: desY => NULL ()
      real(kind=RK), dimension(:), pointer :: desZ => NULL ()
      integer(kind=4) :: nX = 0, mx1 = 0, mx2 = 0 !2012
      integer(kind=4) :: nY = 0, my1 = 0, my2 = 0 !2012
      integer(kind=4) :: nZ = 0, mz1 = 0, mz2 = 0 !2012
      real(kind=RK) ::originx= 0.0_RKIND  !2012
      real(kind=RK) ::originy= 0.0_RKIND  !2012
      real(kind=RK) ::originz= 0.0_RKIND  !2012
   END type Desplazamiento
   !-----------------> Program Types
   !------------------------------------------------------------------------------
   ! Parameters needed for the parser
   !------------------------------------------------------------------------------
   type, public :: Parseador
      character(len=BUFSIZE) :: switches=' '  
      ! Basics
      type(NFDEGeneral), pointer ::         general => NULL ()
      type(MatrizMedios), pointer ::        matriz => NULL ()
      type(Desplazamiento), pointer ::      despl => NULL ()
      type(Frontera), pointer ::            front => NULL ()
      ! Materials
      type(Materials), pointer ::           Mats => NULL ()
      type(PECRegions), pointer ::          pecRegs => NULL ()
      type(PECRegions), pointer ::          pmcRegs => NULL ()
      type(DielectricRegions), pointer ::   DielRegs => NULL ()
      type(LossyThinSurfaces), pointer ::   LossyThinSurfs => NULL ()
      type(FreqDepenMaterials), pointer ::  frqDepMats => NULL ()
      type(ANISOTROPICelements_t), pointer :: aniMats => NULL ()
      ! Sources
      type(Boxes), pointer ::               boxSrc => NULL ()
      type(PlaneWaves), pointer ::          plnSrc => NULL ()
      type(NodSource), pointer ::           nodSrc => NULL ()
      ! Probes
      type(Sondas), pointer ::              oldSONDA => NULL ()
      type(MasSondas), pointer ::           Sonda => NULL ()
      type(BloqueProbes), pointer ::        BloquePrb => NULL ()
      type(VolProbes), pointer ::           VolPrb => NULL ()
      ! Thin Elements                         
      type(ThinWires), pointer ::           tWires => NULL ()
      type(SlantedWires), pointer ::        sWires => NULL ()
      type(ThinSlots), pointer ::           tSlots => NULL ()
      ! Conformal
      type(ConformalPECRegions), pointer ::  conformalRegs => NULL()
#ifdef CompileWithMTLN
      type(mtln_t), pointer ::              mtln => NULL () 
#endif
   END type Parseador
   
   !---> definicion de tipos
   type, public :: t_linea
      integer(kind=4) :: LEN
      character(len=BUFSIZE) :: dato
   END type t_linea
   !--->
   type, public :: t_NFDE_FILE
      integer(kind=4) mpidir !x=1,y=2,z=3
      integer(kind=8) :: targ
      !--->
      integer(kind=8) :: numero
      type(t_linea), dimension(:), pointer :: lineas
      logical :: thereare_stoch
   END type t_NFDE_FILE
!--->

contains


end module NFDETypes

    
