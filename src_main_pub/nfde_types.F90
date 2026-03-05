module NFDETypes
   !
   use FDETYPES
#ifdef CompileWithMTLN   
   use mtln_types_mod
#endif
   use conformal_types_mod
   !
   implicit none
   integer (KIND=4), PARAMETER :: RK = RKIND
   !------------------------------------------------------------------------------
   ! CONSTANTS FOR THE PARSER
   !------------------------------------------------------------------------------
   ! global variable stochastic
   ! MATERIALS
   real(KIND=RK), PARAMETER :: SIGMA_PEC = 1e19_RK
   real(KIND=RK), PARAMETER :: SIGMA_PMC = 1e19_RK 
   ! PROBES
   !!!!
   integer (KIND=4), PARAMETER :: NP_T1_PLAIN = 0
   integer (KIND=4), PARAMETER :: NP_T1_AMBOS = 2
   integer (KIND=4), PARAMETER :: NP_T2_TIME = 0
   integer (KIND=4), PARAMETER :: NP_T2_FREQ = 1
   integer (KIND=4), PARAMETER :: NP_T2_TRANSFER = 2
   integer (KIND=4), PARAMETER :: NP_T2_TIMEFREQ  = 3
   integer (KIND=4), PARAMETER :: NP_T2_TIMETRANSF = 4
   integer (KIND=4), PARAMETER :: NP_T2_FREQTRANSF = 5
   integer (KIND=4), PARAMETER :: NP_T2_TIMEFRECTRANSF = 6
   integer (KIND=4), PARAMETER :: NP_COR_EX = 0
   integer (KIND=4), PARAMETER :: NP_COR_EY = 1
   integer (KIND=4), PARAMETER :: NP_COR_EZ = 2
   integer (KIND=4), PARAMETER :: NP_COR_HX = 3
   integer (KIND=4), PARAMETER :: NP_COR_HY = 4
   integer (KIND=4), PARAMETER :: NP_COR_HZ = 5
   integer (KIND=4), PARAMETER :: NP_COR_WIRECURRENT = 6
   integer (KIND=4), PARAMETER :: NP_COR_DDP = 7
   integer (KIND=4), PARAMETER :: NP_COR_LINE = 8
   integer (KIND=4), PARAMETER :: NP_COR_CHARGE = 9
   LOGICAL, PARAMETER :: BcELECT = .TRUE.
   LOGICAL, PARAMETER :: BcMAGNE = .FALSE.
   ! THIN WIRES
   integer (KIND=4), PARAMETER :: MATERIAL_CONS = 0
   integer (KIND=4), PARAMETER :: MATERIAL_absorbing = 100
   integer (KIND=4), PARAMETER :: Parallel_CONS = 1
   integer (KIND=4), PARAMETER :: SERIES_CONS = 2
   integer (KIND=4), PARAMETER :: DISPERSIVE_CONS = 3
   ! BORDERS
   integer (KIND=4), PARAMETER :: F_PEC = 1
   integer (KIND=4), PARAMETER :: F_PMC = 2
   integer (KIND=4), PARAMETER :: F_PER = 4
   integer (KIND=4), PARAMETER :: F_MUR = 7
   integer (KIND=4), PARAMETER :: F_PML = 9
   integer (KIND=4), PARAMETER :: F_XL = 1
   integer (KIND=4), PARAMETER :: F_XU = 2
   integer (KIND=4), PARAMETER :: F_YL = 3
   integer (KIND=4), PARAMETER :: F_YU = 4
   integer (KIND=4), PARAMETER :: F_ZL = 5
   integer (KIND=4), PARAMETER :: F_ZU = 6
   integer (KIND=4), PARAMETER :: F_TIMEFRECTRANSF = 0
   ! rlc y diodos
   integer (KIND=4), PARAMETER :: inductor = 20
   integer (KIND=4), PARAMETER :: capacitor = 21
   integer (KIND=4), PARAMETER :: resistor = 22
   integer (KIND=4), PARAMETER :: diodo = 23
   integer (KIND=4), PARAMETER :: Dielectric = 24
   integer (KIND=4), PARAMETER :: PMLbody = 25

   !------------------------------------------------------------------------------
   ! TYPES
   !------------------------------------------------------------------------------
   !-----------------> Cordinate Types
   !------------------------------------------------------------------------------
   ! Basic cordinate type for two points and orientation
   !------------------------------------------------------------------------------
   type, public :: coords
      integer (KIND=4) :: Xi = - 1
      integer (KIND=4) :: Xe = - 1
      integer (KIND=4) :: Yi = - 1
      integer (KIND=4) :: Ye = - 1
      integer (KIND=4) :: Zi = - 1
      integer (KIND=4) :: Ze = - 1
      integer (KIND=4) :: Xtrancos = 1
      integer (KIND=4) :: Ytrancos = 1
      integer (KIND=4) :: Ztrancos = 1
      integer (KIND=4) :: Or = 0 !f1eld orientation
      character (LEN=BUFSIZE) :: tag
   END type coords
   type, public :: coords_scaled
      integer (KIND=4) :: Xi = - 1
      integer (KIND=4) :: Xe = - 1
      integer (KIND=4) :: Yi = - 1
      integer (KIND=4) :: Ye = - 1
      integer (KIND=4) :: Zi = - 1
      integer (KIND=4) :: Ze = - 1
      real(KIND=RK) :: xc = 0.0_RKIND
      real(KIND=RK) :: yc = 0.0_RKIND
      real(KIND=RK) :: zc = 0.0_RKIND
      integer (KIND=4) :: Or = 0 !field orientation nuevo 2015
      character (LEN=BUFSIZE) :: tag
   END type coords_scaled
   !-----------------> Material Types
   !------------------------------------------------------------------------------
   ! Basic constants for materials
   !------------------------------------------------------------------------------
   type, public :: Material
      real(KIND=RK) :: eps = 0.0_RKIND
      real(KIND=RK) :: mu = 0.0_RKIND
      real(KIND=RK) :: sigma = 0.0_RKIND
      real(KIND=RK) :: sigmam = 0.0_RKIND
      integer (KIND=4) :: id = 0
   END type Material
   !------------------------------------------------------------------------------
   ! New Class which is a collection of different materials
   !------------------------------------------------------------------------------
   type, public :: Materials
      integer (KIND=4) :: n_Mats = 0
      integer (KIND=4) :: n_Mats_max = 0
      type (Material), dimension(:), POINTER :: Mats => NULL ()
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
      integer (kind=4), dimension(3) :: cell
      integer(kind=4) :: direction = -1
      real(kind=rkind) :: ratio = -1
      real(kind=rkind), dimension(2) :: material_coords
   end type 
   type, public :: face_t 
      integer (kind=4), dimension(3) :: cell
      integer(kind=4) :: direction = -1
      real(kind=rkind) :: ratio = -1
   end type 

   type, public :: conformal_edge_media_t
      type(edge_t), dimension(:), allocatable :: edges
      real(kind=rkind) :: ratio
      integer (kind=4) :: size
   end type
   type, public :: conformal_face_media_t
      type(face_t), dimension(:), allocatable :: faces
      real(kind=rkind) :: ratio
      integer (kind=4) :: size
   end type

   type, public :: ConformalMedia_t
      integer (KIND=4) :: n_edges_media = 0
      integer (KIND=4) :: n_faces_media = 0
      type (conformal_face_media_t), dimension(:), POINTER :: face_media => NULL ()
      type (conformal_edge_media_t), dimension(:), POINTER :: edge_media => NULL ()
      real(kind=rkind) :: time_step_scale_factor = 1.0
      character(len=bufsize) :: tag
   END type ConformalMedia_t


   !------------------------------------------------------------------------------
   ! Locates all the different PEC media found
   !------------------------------------------------------------------------------

   type, public :: PECRegions
      integer (KIND=4) :: nVols = 0
      integer (KIND=4) :: nSurfs = 0
      integer (KIND=4) :: nLins = 0
      integer (KIND=4) :: nVols_max = 0
      integer (KIND=4) :: nSurfs_max = 0
      integer (KIND=4) :: nLins_max = 0
      type (coords), dimension(:), POINTER :: Vols => NULL ()
      type (coords), dimension(:), POINTER :: Surfs => NULL ()
      type (coords), dimension(:), POINTER :: Lins => NULL ()
   END type PECRegions
   !------------------------------------------------------------------------------
   ! Defines a Non Metal Body
   !------------------------------------------------------------------------------
   type, public :: Dielectric_t
      type (coords), dimension(:), POINTER :: c1P => NULL ()
      type (coords), dimension(:), POINTER :: c2P => NULL ()
      real(KIND=RK) :: sigma = 0.0_RKIND
      real(KIND=RK) :: eps = 0.0_RKIND
      real(KIND=RK) :: mu = 0.0_RKIND
      real(KIND=RK) :: sigmam = 0.0_RKIND
      integer (KIND=4) :: n_C1P = 0
      integer (KIND=4) :: n_C2P = 0
      !!! Lumped vars.
      real(KIND=RK)   :: Rtime_on = 0.0_RKIND, Rtime_off = 0.0_RKIND
      real(KIND=RK) :: R = 0.0_RKIND
      real(KIND=RK) :: L = 0.0_RKIND
      real(KIND=RK) :: C = 0.0_RKIND
      !!! Stoch vars.
      real(KIND=RK) :: R_devia = 0.0_RKIND
      real(KIND=RK) :: L_devia = 0.0_RKIND
      real(KIND=RK) :: C_devia = 0.0_RKIND
      !
      real(KIND=RK) :: DiodB = 0.0_RKIND
      real(KIND=RK)   :: DiodIsat = 0.0_RKIND
      integer (KIND=4) :: DiodOri = 0
      !!! Berenger's waveports
      integer (KIND=4) :: orient = 0
!!!!!!!!!
      logical :: resistor=.false. , inductor=.false. , capacitor=.false. , diodo=.false. , plain=.false. , PMLbody=.false.
   END type Dielectric_t
   !------------------------------------------------------------------------------
   ! Locates all the different Non Metal Media found
   !------------------------------------------------------------------------------
   type, public :: DielectricRegions
      type (Dielectric_t), dimension(:), POINTER :: Vols => NULL ()
      type (Dielectric_t), dimension(:), POINTER :: Surfs => NULL ()
      type (Dielectric_t), dimension(:), POINTER :: Lins => NULL ()
      integer (KIND=4) :: nVols = 0
      integer (KIND=4) :: nSurfs = 0
      integer (KIND=4) :: nLins = 0
      integer (KIND=4) :: nVols_max = 0
      integer (KIND=4) :: nSurfs_max = 0
      integer (KIND=4) :: nLins_max = 0
      integer (KIND=4) :: n_C1P_max = 0
      integer (KIND=4) :: n_C2P_max = 0
   END type DielectricRegions
   !------------------------------------------------------------------------------
   ! type that defines the information of a frequency depENDent material,
   ! it inherits from the material class and it adds the possible values needed
   ! in the frequency depENDent section of the
   !------------------------------------------------------------------------------
   type, public :: FreqDepenMaterial
      COMPLEX, dimension(:), POINTER ::  a11 => NULL ()
      COMPLEX, dimension(:), POINTER ::  b11 => NULL ()
      COMPLEX, dimension(:), POINTER :: am11 => NULL ()
      COMPLEX, dimension(:), POINTER :: bm11 => NULL ()
      COMPLEX, dimension(:), POINTER ::  a12 => NULL ()
      COMPLEX, dimension(:), POINTER ::  b12 => NULL ()
      COMPLEX, dimension(:), POINTER :: am12 => NULL ()
      COMPLEX, dimension(:), POINTER :: bm12 => NULL ()
      COMPLEX, dimension(:), POINTER ::  a13 => NULL ()
      COMPLEX, dimension(:), POINTER ::  b13 => NULL ()
      COMPLEX, dimension(:), POINTER :: am13 => NULL ()
      COMPLEX, dimension(:), POINTER :: bm13 => NULL ()
      COMPLEX, dimension(:), POINTER ::  a22 => NULL ()
      COMPLEX, dimension(:), POINTER ::  b22 => NULL ()
      COMPLEX, dimension(:), POINTER :: am22 => NULL ()
      COMPLEX, dimension(:), POINTER :: bm22 => NULL ()
      COMPLEX, dimension(:), POINTER ::  a23 => NULL ()
      COMPLEX, dimension(:), POINTER ::  b23 => NULL ()
      COMPLEX, dimension(:), POINTER :: am23 => NULL ()
      COMPLEX, dimension(:), POINTER :: bm23 => NULL ()
      COMPLEX, dimension(:), POINTER ::  a33 => NULL ()
      COMPLEX, dimension(:), POINTER ::  b33 => NULL ()
      COMPLEX, dimension(:), POINTER :: am33 => NULL ()
      COMPLEX, dimension(:), POINTER :: bm33 => NULL ()
      real(KIND=RK), dimension(:), POINTER :: alpha => NULL ()
      real(KIND=RK), dimension(:), POINTER :: beta => NULL ()
      real(KIND=RK), dimension(:), POINTER :: gamma => NULL ()
      real(KIND=RK), dimension(:), POINTER :: alpham => NULL ()
      real(KIND=RK), dimension(:), POINTER :: betam => NULL ()
      real(KIND=RK), dimension(:), POINTER :: gammam => NULL ()
      type (coords), dimension(:), POINTER :: c => NULL ()
      real(KIND=RK) ::    eps11 = 0.0_RKIND ,    eps12 = 0.0_RKIND ,    eps13 = 0.0_RKIND ,    eps22 = 0.0_RKIND ,    eps23 = 0.0_RKIND ,    eps33 = 0.0_RKIND
      real(KIND=RK) ::     mu11 = 0.0_RKIND ,     mu12 = 0.0_RKIND ,     mu13 = 0.0_RKIND ,     mu22 = 0.0_RKIND ,     mu23 = 0.0_RKIND ,     mu33 = 0.0_RKIND
      real(KIND=RK) ::  sigma11 = 0.0_RKIND ,  sigma12 = 0.0_RKIND ,  sigma13 = 0.0_RKIND ,  sigma22 = 0.0_RKIND ,  sigma23 = 0.0_RKIND ,  sigma33 = 0.0_RKIND
      real(KIND=RK) :: sigmam11 = 0.0_RKIND , sigmam12 = 0.0_RKIND , sigmam13 = 0.0_RKIND , sigmam22 = 0.0_RKIND , sigmam23 = 0.0_RKIND , sigmam33 = 0.0_RKIND
      integer (KIND=4) ::  K11 = 0
      integer (KIND=4) :: Km11 = 0
      integer (KIND=4) ::  K12 = 0
      integer (KIND=4) :: Km12 = 0
      integer (KIND=4) ::  K13 = 0
      integer (KIND=4) :: Km13 = 0
      integer (KIND=4) ::  K22 = 0
      integer (KIND=4) :: Km22 = 0
      integer (KIND=4) ::  K23 = 0
      integer (KIND=4) :: Km23 = 0
      integer (KIND=4) ::  K33 = 0
      integer (KIND=4) :: Km33 = 0
      integer (KIND=4) :: L = 0
      integer (KIND=4) :: Lm = 0
      integer (KIND=4) :: n_c = 0
      character (LEN=BUFSIZE) :: files = ' ' !2015 si esta presente lee los polos/residuos desde fichero
   END type FreqDepenMaterial
   !------------------------------------------------------------------------------
   ! type that defines the list of frequency depedent materials
   !------------------------------------------------------------------------------
   type, public :: FreqDepenMaterials
      type (FreqDepenMaterial), dimension(:), POINTER :: Vols => NULL ()
      type (FreqDepenMaterial), dimension(:), POINTER :: Surfs => NULL ()
      type (FreqDepenMaterial), dimension(:), POINTER :: Lins => NULL ()
      integer (KIND=4) :: nVols = 0
      integer (KIND=4) :: nSurfs = 0
      integer (KIND=4) :: nLins = 0
      integer (KIND=4) :: nVols_max = 0
      integer (KIND=4) :: nSurfs_max = 0
      integer (KIND=4) :: nLins_max = 0
      integer (KIND=4) :: n_c_max = 0 !cota superior
   END type FreqDepenMaterials
   !------------------------------------------------------------------------------
   ! Type for the ANISOTROPIC body, surface and lines since they will contain
   ! the same information
   !------------------------------------------------------------------------------
   type, public :: ANISOTROPICbody_t
      type (coords), dimension(:), POINTER :: c1P => NULL ()
      type (coords), dimension(:), POINTER :: c2P => NULL ()
      real(KIND=RK), dimension(3, 3) :: sigma, eps, mu, sigmam
      integer (KIND=4) :: n_C1P = 0
      integer (KIND=4) :: n_C2P = 0
   END type ANISOTROPICbody_t
   !------------------------------------------------------------------------------
   ! Type that contains the elements found in the nfde File
   !------------------------------------------------------------------------------
   type, public :: ANISOTROPICelements_t
      type (ANISOTROPICbody_t), dimension(:), POINTER :: Vols => NULL ()
      type (ANISOTROPICbody_t), dimension(:), POINTER :: Surfs => NULL ()
      type (ANISOTROPICbody_t), dimension(:), POINTER :: Lins => NULL ()
      integer (KIND=4) :: nVols = 0
      integer (KIND=4) :: nSurfs = 0
      integer (KIND=4) :: nLins = 0
      integer (KIND=4) :: nVols_max = 0
      integer (KIND=4) :: nSurfs_max = 0
      integer (KIND=4) :: nLins_max = 0
      integer (KIND=4) :: n_C1P_max = 0 !cota superior de c1p y c2p en vols,sufs,lins
      integer (KIND=4) :: n_C2P_max = 0
   END type ANISOTROPICelements_t
   !------------------------------------------------------------------------------
   ! Defines a Comp Surface
   !------------------------------------------------------------------------------
   type, public :: LossyThinSurface
      type (coords), dimension(:), POINTER :: c => NULL ()
      real(KIND=RK), dimension(:), POINTER :: sigma
      real(KIND=RK), dimension(:), POINTER :: eps
      real(KIND=RK), dimension(:), POINTER :: mu
      real(KIND=RK), dimension(:), POINTER :: sigmam
      real(KIND=RK), dimension(:), POINTER :: thk
      !for_devia
      real(KIND=RK), dimension(:), POINTER :: sigma_devia
      real(KIND=RK), dimension(:), POINTER :: eps_devia
      real(KIND=RK), dimension(:), POINTER :: mu_devia
      real(KIND=RK), dimension(:), POINTER :: sigmam_devia
      real(KIND=RK), dimension(:), POINTER :: thk_devia
      !
      integer (KIND=4) :: nc = 0
      character (LEN=BUFSIZE) :: files = ' ' 
      integer (KIND=4)  :: numcapas  
   END type LossyThinSurface
   !------------------------------------------------------------------------------
   ! Locates all the different Comp media found
   !------------------------------------------------------------------------------
   type, public :: LossyThinSurfaces
      type (LossyThinSurface), dimension(:), POINTER :: cs => NULL ()
      integer (KIND=4) :: length = 0
      integer (KIND=4) :: length_max = 0
      integer (KIND=4) :: nC_max = 0 !cota de todos los nc de LossyThinSurface
   END type LossyThinSurfaces
   !------------------------------------------------------------------------------
   ! Component for Thin Wires there is a list of this inside the component
   ! that defines the whole Thin Wire Reference
   !------------------------------------------------------------------------------
   type, public :: ThinWireComp
      character (LEN=BUFSIZE) :: srctype, srcfile
      integer (KIND=4) :: i = - 1
      integer (KIND=4) :: j = - 1
      integer (KIND=4) :: K = - 1
      integer (KIND=4) :: nd = - 1
      integer (KIND=4) :: d = - 1
      real(KIND=RK) :: m = 0.0_RKIND
      character (LEN=BUFSIZE) :: tag
   END type ThinWireComp
   !------------------------------------------------------------------------------
   ! ThinWire component that defines the overall properties of the definition
   ! of ThinWires
   !------------------------------------------------------------------------------
   type, public :: ThinWire
      type (ThinWireComp), dimension(:), POINTER :: twc => NULL ()
      real(KIND=RK) :: rad = 0 , rad_devia = 0
      LOGICAL :: disp = .false.
      character (LEN=BUFSIZE) :: dispfile
      real(KIND=RK) :: res = 0 , res_devia = 0
      real(KIND=RK) :: ind = 0 , ind_devia = 0
      real(KIND=RK) :: cap = 0 , cap_devia = 0
      real(KIND=RK) :: P_res = 0
      real(KIND=RK) :: P_ind = 0
      real(KIND=RK) :: P_cap = 0
      character (LEN=BUFSIZE) :: dispfile_LeftEnd
      real(KIND=RK) :: R_LeftEnd = 0 , R_LeftEnd_devia = 0
      real(KIND=RK) :: L_LeftEnd = 0 , L_LeftEnd_devia = 0
      real(KIND=RK) :: C_LeftEnd = 0 , C_LeftEnd_devia = 0
      character (LEN=BUFSIZE) :: dispfile_RightEnd
      real(KIND=RK) :: R_RightEnd = 0 , R_RightEnd_devia = 0
      real(KIND=RK) :: L_RightEnd = 0 , L_RightEnd_devia = 0
      real(KIND=RK) :: C_RightEnd = 0 , C_RightEnd_devia = 0
      integer (KIND=4) :: LeftEnd = 0
      integer (KIND=4) :: RightEnd = 0
      ! Components
      integer (KIND=4) :: tl = 0
      integer (KIND=4) :: tr = 0
      integer (KIND=4) :: n_twc = 0
      integer (KIND=4) :: n_twc_max = 0
   END type ThinWire
   !------------------------------------------------------------------------------
   ! List of the different thin wires that were found in the file
   !------------------------------------------------------------------------------
   type, public :: ThinWires
      type (ThinWire), dimension(:), POINTER :: tw => NULL ()
      integer (KIND=4) :: n_tw = 0
      integer (KIND=4) :: n_tw_max = 0
   END type ThinWires
   !------------------------------------------------------------------------------
   ! Component for Slanted Wires there is a list of this inside the component
   ! that defines the whole Slanted Wire Reference
   !------------------------------------------------------------------------------
   type, public :: SlantedWireComp
      character (LEN=BUFSIZE) :: srctype, srcfile
      real(KIND=RK) :: x = - 1.0_RKIND
      real(KIND=RK) :: y = - 1.0_RKIND
      real(KIND=RK) :: z = - 1.0_RKIND
      integer (KIND=4) :: nd = - 1
      real(KIND=RK) :: m = 0.0_RKIND
      character (LEN=BUFSIZE) :: tag
   END type SlantedWireComp
   !------------------------------------------------------------------------------
   ! ThinWire component that defines the overall properties of the definition
   ! of ThinWires
   !------------------------------------------------------------------------------
   type, public :: SlantedWire
      type (SlantedWireComp), dimension(:), POINTER :: swc => NULL ()
      real(KIND=RK) :: rad = 0
      LOGICAL :: disp = .false.
      character (LEN=BUFSIZE) :: dispfile
      real(KIND=RK) :: res = 0
      real(KIND=RK) :: ind = 0
      real(KIND=RK) :: cap = 0
      real(KIND=RK) :: P_res = 0
      real(KIND=RK) :: P_ind = 0
      real(KIND=RK) :: P_cap = 0
      character (LEN=BUFSIZE) :: dispfile_LeftEnd
      real(KIND=RK) :: R_LeftEnd = 0
      real(KIND=RK) :: L_LeftEnd = 0
      real(KIND=RK) :: C_LeftEnd = 0
      character (LEN=BUFSIZE) :: dispfile_RightEnd
      real(KIND=RK) :: R_RightEnd = 0
      real(KIND=RK) :: L_RightEnd = 0
      real(KIND=RK) :: C_RightEnd = 0
      integer (KIND=4) :: LeftEnd = 0
      integer (KIND=4) :: RightEnd = 0
      ! Components
      integer (KIND=4) :: tl = 0
      integer (KIND=4) :: tr = 0
      integer (KIND=4) :: n_swc = 0
      integer (KIND=4) :: n_swc_max = 0
   END type SlantedWire
   !------------------------------------------------------------------------------
   ! List of the different thin wires that were found in the file
   !------------------------------------------------------------------------------
   type, public :: SlantedWires
      type (SlantedWire), dimension(:), POINTER :: sw => NULL ()
      integer (KIND=4) :: n_sw = 0
      integer (KIND=4) :: n_sw_max = 0
   END type SlantedWires
   !--------------------------------------------------------------------------
   ! Component for Thin Slots there is a list of this inside the component
   ! that defines the whole Thin Slot Reference
   !--------------------------------------------------------------------------
   type, public :: ThinSlotComp
      integer (KIND=4) :: i = 0
      integer (KIND=4) :: j = 0
      integer (KIND=4) :: K = 0
      integer (KIND=4) :: node = 0
      integer (KIND=4) :: dir = - 1
      integer (KIND=4) :: Or = - 1
      character (LEN=BUFSIZE) :: tag
   END type ThinSlotComp
   !--------------------------------------------------------------------------
   ! ThinSlot component that defines the overall properties of the definition
   ! of ThinSlots in ORIGINAL
   !--------------------------------------------------------------------------
   type, public :: ThinSlot
      type (ThinSlotComp), dimension(:), POINTER :: tgc => NULL ()
      real(KIND=RK) :: width = 0
      integer (KIND=4) :: n_tgc = 0
      integer (KIND=4) :: n_tgc_max = 0
   END type ThinSlot
   !--------------------------------------------------------------------------
   ! List of the different thin Slots that were found in the file
   !--------------------------------------------------------------------------
   type, public :: ThinSlots
      type (ThinSlot), dimension(:), POINTER :: tg => NULL ()
      integer (KIND=4) :: n_tg = 0
      integer (KIND=4) :: n_tg_max = 0
   END type ThinSlots
   !-----------------> Border Types
   !------------------------------------------------------------------------------
   ! PML Border Type
   !------------------------------------------------------------------------------
   type, public :: FronteraPML
      real(KIND=RK) :: orden = 2.0_RK
      real(KIND=RK) :: refl = 1e-3_RK
      integer (KIND=4) :: numCapas = 8
   END type FronteraPML
   !------------------------------------------------------------------------------
   ! Tipo de la frontera
   !------------------------------------------------------------------------------
   type, public :: Frontera
      integer (KIND=4), dimension(6) :: tipoFrontera
      type (FronteraPML), dimension(6) :: propiedadesPML
   END type Frontera
   !-----------------> Probe Types
   !------------------------------------------------------------------------------
   ! type to define the new probe object which contains the type of calculation
   ! the type of analysis, time and frequency step and the filename where
   ! it should be saved
   !------------------------------------------------------------------------------
   type, public :: MasSonda
      character (LEN=BUFSIZE) :: filename
      type (coords), dimension(:), POINTER :: cordinates => NULL ()
      real(KIND=RK) :: tstart, tstop, tstep
      real(KIND=RK) :: fstart, fstop, fstep
      integer (KIND=4) :: type1, type2
      integer (KIND=4) :: len_cor = 0
      character (LEN=BUFSIZE) :: outputrequest
   END type MasSonda
   !------------------------------------------------------------------------------
   ! type that defines a list of probes to be appended and accesed
   !------------------------------------------------------------------------------
   type, public :: MasSondas
      type (MasSonda), dimension(:), POINTER :: collection => NULL ()
      integer (KIND=4) :: length = 0
      integer (KIND=4) :: length_max = 0
      integer (KIND=4) :: len_cor_max = 0 !cota
   END type MasSondas
   !------------------------------------------------------------------------------
   ! This type contains the basic information in nearly all the different PROBES
   !------------------------------------------------------------------------------
   type, public :: Sonda
      character (LEN=BUFSIZE) :: grname
      integer (KIND=4), dimension(:), POINTER :: i => NULL ()
      integer (KIND=4), dimension(:), POINTER :: j => NULL ()
      integer (KIND=4), dimension(:), POINTER :: K => NULL ()
      integer (KIND=4), dimension(:), POINTER :: node => NULL ()
      integer (KIND=4) :: n_cord = 0
      integer (KIND=4) :: n_cord_max = 0
      real(KIND=RK) :: tstart, tstop, tstep
      character (LEN=BUFSIZE) :: outputrequest
      !por si se precisa para el Far Field
      real(KIND=RK) :: fstart, fstop, fstep
      real(KIND=RK) :: phistart, phistop, phistep
      real(KIND=RK) :: thetastart, thetastop, thetastep
      character (LEN=BUFSIZE) :: FileNormalize
   END type Sonda
   !------------------------------------------------------------------------------
   ! type for the electric far field
   !------------------------------------------------------------------------------
   type, public :: FarField_Sonda
      type (Sonda) :: probe
   END type FarField_Sonda
   !------------------------------------------------------------------------------
   ! type for the electric field
   !------------------------------------------------------------------------------
   type, public :: Electric_Sonda
      type (Sonda) :: probe
   END type Electric_Sonda
   !------------------------------------------------------------------------------
   ! type for the magnetic field
   !------------------------------------------------------------------------------
   type, public :: Magnetic_Sonda
      type (Sonda) :: probe
   END type Magnetic_Sonda
   !------------------------------------------------------------------------------
   ! type for the normal electric field
   !------------------------------------------------------------------------------
   type, public :: NormalElectric_Sonda
      type (Sonda) :: probe
      integer (KIND=4), dimension(:), POINTER :: nml => NULL ()
      integer (KIND=4) :: n_nml = 0
      integer (KIND=4) :: n_nml_max = 0
   END type NormalElectric_Sonda
   !------------------------------------------------------------------------------
   ! type for the normal magnetic field
   !------------------------------------------------------------------------------
   type, public :: NormalMagnetic_Sonda
      type (Sonda) :: probe
      integer (KIND=4), dimension(:), POINTER :: nml => NULL ()
      integer (KIND=4) :: n_nml = 0
      integer (KIND=4) :: n_nml_max = 0
   END type NormalMagnetic_Sonda
   !------------------------------------------------------------------------------
   ! type for the electric surface current density
   !------------------------------------------------------------------------------
   type, public :: SurfaceElectricCurrent_Sonda
      type (Sonda) :: probe
      integer (KIND=4), dimension(:), POINTER :: nml => NULL ()
      integer (KIND=4) :: n_nml = 0
      integer (KIND=4) :: n_nml_max = 0
   END type SurfaceElectricCurrent_Sonda
   !------------------------------------------------------------------------------
   ! type for the magnetic surface current density
   !------------------------------------------------------------------------------
   type, public :: SurfaceMagneticCurrent_Sonda
      type (Sonda) :: probe
      integer (KIND=4), dimension(:), POINTER :: nml => NULL ()
      integer (KIND=4) :: n_nml = 0
      integer (KIND=4) :: n_nml_max = 0
   END type SurfaceMagneticCurrent_Sonda
   !------------------------------------------------------------------------------
   ! Abstract class which performs the dynamic dispatching
   !------------------------------------------------------------------------------
   type, public :: abstractSonda
      integer (KIND=4) :: n_FarField = 0 
      integer (KIND=4) :: n_Electric = 0
      integer (KIND=4) :: n_Magnetic = 0
      integer (KIND=4) :: n_NormalElectric = 0
      integer (KIND=4) :: n_NormalMagnetic = 0
      integer (KIND=4) :: n_SurfaceElectricCurrent = 0
      integer (KIND=4) :: n_SurfaceMagneticCurrent = 0
      !
      integer (KIND=4) :: n_FarField_max = 0
      integer (KIND=4) :: n_Electric_max = 0
      integer (KIND=4) :: n_Magnetic_max = 0
      integer (KIND=4) :: n_NormalElectric_max = 0
      integer (KIND=4) :: n_NormalMagnetic_max = 0
      integer (KIND=4) :: n_SurfaceElectricCurrent_max = 0
      integer (KIND=4) :: n_SurfaceMagneticCurrent_max = 0
      type (FarField_Sonda), dimension(:), POINTER :: FarField => NULL ()
      type (Electric_Sonda), dimension(:), POINTER :: Electric => NULL ()
      type (Magnetic_Sonda), dimension(:), POINTER :: Magnetic => NULL ()
      type (NormalElectric_Sonda), dimension(:), POINTER :: NormalElectric => NULL ()
      type (NormalMagnetic_Sonda), dimension(:), POINTER :: NormalMagnetic => NULL ()
      type (SurfaceElectricCurrent_Sonda), dimension(:), POINTER :: SurfaceElectricCurrent => NULL ()
      type (SurfaceMagneticCurrent_Sonda), dimension(:), POINTER :: SurfaceMagneticCurrent => NULL ()
   END type abstractSonda
   !------------------------------------------------------------------------------
   ! Class to account as a list for all the probes
   ! that might be required during the parsing process
   !------------------------------------------------------------------------------
   type, public :: Sondas
      type (abstractSonda), dimension(:), POINTER :: probes => NULL ()
      integer (KIND=4) :: n_probes = 0
      integer (KIND=4) :: n_probes_max = 0
   END type Sondas
   !------------------------------------------------------------------------------
   ! Object type defined for the Bloque current probe
   !------------------------------------------------------------------------------
   type, public :: BloqueProbe
      real(KIND=RK) :: tstart, tstop, tstep
      real(KIND=RK) :: fstart, fstop, fstep
      character (LEN=BUFSIZE) :: FileNormalize
      integer (KIND=4) :: type2
      integer (KIND=4) :: i1, i2, j1, j2, k1, k2, skip
      integer (KIND=4) :: nml
      LOGICAL :: t
      character (LEN=BUFSIZE) :: outputrequest
      character (LEN=BUFSIZE) :: tag
   END type BloqueProbe
   ! Object made for the collection of defined Bloque probes
   type, public :: BloqueProbes
      type (BloqueProbe), dimension(:), POINTER :: bp => NULL ()
      integer (KIND=4) :: n_bp = 0
      integer (KIND=4) :: n_bp_max = 0
   END type BloqueProbes

   !------------------------------------------------------------------------------
   ! Object type defined for the Volumic probes
   !------------------------------------------------------------------------------
   type, public :: VolProbe
      type (coords), dimension(:), POINTER :: cordinates => NULL ()
      real(KIND=RK) :: tstart, tstop, tstep
      character (LEN=BUFSIZE) :: outputrequest
      integer (KIND=4) :: len_cor = 0
      !para freq domain
      real(KIND=RK) :: fstart, fstop, fstep
      integer (KIND=4) ::  type2
      character (LEN=BUFSIZE) :: filename
   END type VolProbe
   ! Object made for the collection of defined Volumic probes
   type, public :: VolProbes
      type (VolProbe), dimension(:), POINTER :: collection => NULL ()
      integer (KIND=4) :: length = 0
      integer (KIND=4) :: length_max = 0
      integer (KIND=4) :: len_cor_max = 0 !cota
   END type VolProbes

   !-----------------> Source Types
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: Box
      character (LEN=BUFSIZE) :: nombre_fichero
      integer (KIND=4), dimension(3) :: coor1, coor2
   END type Box
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: Boxes
      type (Box), dimension(:), POINTER :: Vols => NULL ()
      integer (KIND=4) :: nVols = 0
      integer (KIND=4) :: nVols_max = 0
   END type Boxes
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: PlaneWave
      character (LEN=BUFSIZE) :: nombre_fichero
      character (LEN=BUFSIZE) :: atributo
      integer (KIND=4), dimension(3) :: coor1, coor2
      real(KIND=RK) :: theta, phi, alpha, beta
      logical :: isRC !for reververation chambers
      real(KIND=RK) :: INCERTMAX
      integer (KIND=4) :: numModes !for reververation chambers
   END type PlaneWave
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: PlaneWaves
      type (PlaneWave), dimension(:), POINTER :: collection => NULL ()
      integer (KIND=4) :: nc = 0
      integer (KIND=4) :: nC_max = 0
   END type PlaneWaves
   !------------------------------------------------------------------------------
   ! Definicin de los tipos current density que existirn en el ficero
   ! nfde
   !------------------------------------------------------------------------------
   type, public :: Curr_Field_Src
      type (coords_scaled), dimension(:), POINTER :: c1P => NULL ()
      type (coords_scaled), dimension(:), POINTER :: c2P => NULL ()
      character (LEN=BUFSIZE) :: nombre
      integer (KIND=4) :: n_C1P = 0
      integer (KIND=4) :: n_C2P = 0
      LOGICAL :: isElec, isHard, isInitialValue
   END type Curr_Field_Src
   !------------------------------------------------------------------------------
   ! Definicin de las Nodal Source global
   !------------------------------------------------------------------------------
   type, public :: NodSource
      type (Curr_Field_Src), dimension(:), POINTER :: NodalSource => NULL ()
      integer (KIND=4) :: n_nodSrc = 0
      integer (KIND=4) :: n_nodSrc_max = 0
      integer (KIND=4) :: n_C1P_max = 0
      integer (KIND=4) :: n_C2P_max = 0
   END type NodSource
   !-----------------> General Types
   !------------------------------------------------------------------------------
   ! Matrix attributes.
   ! Total[XYZ] -> Is the cell number for each axis.
   !------------------------------------------------------------------------------
   type, public :: MatrizMedios
      integer (KIND=4) :: totalX, totalY, totalZ
   END type MatrizMedios
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type NFDEGeneral
      real(KIND=RK) :: dt
      integer (KIND=4) :: nmax
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
      real(KIND=RK), dimension(:), POINTER :: desX => NULL ()
      real(KIND=RK), dimension(:), POINTER :: desY => NULL ()
      real(KIND=RK), dimension(:), POINTER :: desZ => NULL ()
      integer (KIND=4) :: nX = 0, mx1 = 0, mx2 = 0 !2012
      integer (KIND=4) :: nY = 0, my1 = 0, my2 = 0 !2012
      integer (KIND=4) :: nZ = 0, mz1 = 0, mz2 = 0 !2012
      real(KIND=RK) ::originx= 0.0_RKIND  !2012
      real(KIND=RK) ::originy= 0.0_RKIND  !2012
      real(KIND=RK) ::originz= 0.0_RKIND  !2012
   END type Desplazamiento
   !-----------------> Program Types
   !------------------------------------------------------------------------------
   ! Parameters needed for the parser
   !------------------------------------------------------------------------------
   type, public :: Parseador
      character (len=BUFSIZE) :: switches=' '  
      ! Basics
      type (NFDEGeneral), POINTER ::           general => NULL ()
      type (MatrizMedios), POINTER ::          matriz => NULL ()
      type (Desplazamiento), POINTER ::        despl => NULL ()
      type (Frontera), POINTER ::              front => NULL ()
      ! Materials
      type (Materials), POINTER ::             Mats => NULL ()
      type (PECRegions), POINTER ::            pecRegs => NULL ()
      type (PECRegions), POINTER ::            pmcRegs => NULL ()
      type (DielectricRegions), POINTER ::     DielRegs => NULL ()
      type (LossyThinSurfaces), POINTER ::     LossyThinSurfs => NULL ()
      type (FreqDepenMaterials), POINTER ::    frqDepMats => NULL ()
      type (ANISOTROPICelements_t), POINTER :: aniMats => NULL ()
      ! Sources
      type (Boxes), POINTER ::                 boxSrc => NULL ()
      type (PlaneWaves), POINTER ::            plnSrc => NULL ()
      type (NodSource), POINTER ::             nodSrc => NULL ()
      ! Probes
      type (Sondas), POINTER ::                oldSONDA => NULL ()
      type (MasSondas), POINTER ::             Sonda => NULL ()
      type (BloqueProbes), POINTER ::          BloquePrb => NULL ()
      type (VolProbes), POINTER ::             VolPrb => NULL ()
      ! Thin Elements                         
      type (ThinWires), POINTER ::             tWires => NULL ()
      type (SlantedWires), POINTER ::          sWires => NULL ()
      type (ThinSlots), POINTER ::             tSlots => NULL ()
      ! Conformal
      type(ConformalPECRegions), pointer ::    conformalRegs => NULL()
#ifdef CompileWithMTLN
      type (mtln_t), POINTER ::                mtln => NULL () 
#endif
   END type Parseador
   
   !---> definicion de tipos
   type, public :: t_linea
      integer (KIND=4) :: LEN
      character (LEN=BUFSIZE) :: dato
   END type t_linea
   !--->
   type, public :: t_NFDE_FILE
      integer (KIND=4) mpidir !x=1,y=2,z=3
      integer (KIND=8) :: targ
      !--->
      integer (KIND=8) :: numero
      type (t_linea), dimension(:), POINTER :: lineas
      logical :: thereare_stoch
   END type t_NFDE_FILE
!--->

contains


end module NFDETypes

    
