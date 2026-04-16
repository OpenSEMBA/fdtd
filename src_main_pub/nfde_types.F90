module NFDETypes_m
   !
   use FDETYPES_m
#ifdef CompileWithMTLN   
   use mtln_types_m
#endif
   use conformal_types_m
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
   type, public :: coords_t
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
   end type coords_t
   type, public :: coords_scaled_t
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
   end type coords_scaled_t
   !-----------------> Material Types
   !------------------------------------------------------------------------------
   ! Basic constants for materials
   !------------------------------------------------------------------------------
   type, public :: Material_t
      real(kind=RK) :: eps = 0.0_RKIND
      real(kind=RK) :: mu = 0.0_RKIND
      real(kind=RK) :: sigma = 0.0_RKIND
      real(kind=RK) :: sigmam = 0.0_RKIND
      integer(kind=4) :: id = 0
   end type Material_t
   !------------------------------------------------------------------------------
   ! New Class which is a collection of different materials
   !------------------------------------------------------------------------------
   type, public :: Materials_t
      integer(kind=4) :: n_Mats = 0
      integer(kind=4) :: n_Mats_max = 0
      type(Material_t), dimension(:), pointer :: Mats => NULL ()
   end type Materials_t

   !------------------------------------------------------------------------------
   ! Identifies conformal PEC "media"
   !------------------------------------------------------------------------------

   type, public :: ConformalPECElements_t
      type(triangle_t), dimension(:), allocatable :: triangles
      type(interval_t), dimension(:), allocatable :: intervals
      character(len=bufsize) :: tag
   end type 

   type, public :: ConformalPECRegions_t
      type(ConformalPECElements_t), dimension(:), pointer :: volumes => null()
      type(ConformalPECElements_t), dimension(:), pointer :: surfaces => null()
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
      integer(kind=4) :: n_elements
   end type
   type, public :: conformal_face_media_t
      type(face_t), dimension(:), allocatable :: faces
      real(kind=rkind) :: ratio
      integer(kind=4) :: n_elements
   end type

   type, public :: ConformalMedia_t
      integer(kind=4) :: n_edges_media = 0
      integer(kind=4) :: n_faces_media = 0
      type(conformal_face_media_t), dimension(:), pointer :: face_media => NULL ()
      type(conformal_edge_media_t), dimension(:), pointer :: edge_media => NULL ()
      real(kind=rkind) :: time_step_scale_factor = 1.0
      character(len=bufsize) :: tag
   end type ConformalMedia_t


   !------------------------------------------------------------------------------
   ! Locates all the different PEC media found
   !------------------------------------------------------------------------------

   type, public :: PECRegions_t
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nSurfs = 0
      integer(kind=4) :: nLins = 0
      integer(kind=4) :: nVols_max = 0
      integer(kind=4) :: nSurfs_max = 0
      integer(kind=4) :: nLins_max = 0
      type(coords_t), dimension(:), pointer :: Vols => NULL ()
      type(coords_t), dimension(:), pointer :: Surfs => NULL ()
      type(coords_t), dimension(:), pointer :: Lins => NULL ()
   end type PECRegions_t
   !------------------------------------------------------------------------------
   ! Defines a Non Metal Body
   !------------------------------------------------------------------------------
   type, public :: Dielectric_t
      type(coords_t), dimension(:), pointer :: c1P => NULL ()
      type(coords_t), dimension(:), pointer :: c2P => NULL ()
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
   end type Dielectric_t
   !------------------------------------------------------------------------------
   ! Locates all the different Non Metal Media found
   !------------------------------------------------------------------------------
   type, public :: DielectricRegions_t
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
   end type DielectricRegions_t
   !------------------------------------------------------------------------------
   ! type that defines the information of a frequency depENDent material,
   ! it inherits from the material class and it adds the possible values needed
   ! in the frequency depENDent section of the
   !------------------------------------------------------------------------------
   type, public :: FreqDepenMaterial_t
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
      type(coords_t), dimension(:), pointer :: c => NULL ()
      real(kind=RK) :: eps11 = 0.0_RKIND ,    eps12 = 0.0_RKIND ,    eps13 = 0.0_RKIND ,    eps22 = 0.0_RKIND ,    eps23 = 0.0_RKIND ,    eps33 = 0.0_RKIND
      real(kind=RK) :: mu11 = 0.0_RKIND ,     mu12 = 0.0_RKIND ,     mu13 = 0.0_RKIND ,     mu22 = 0.0_RKIND ,     mu23 = 0.0_RKIND ,     mu33 = 0.0_RKIND
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
   end type FreqDepenMaterial_t
   !------------------------------------------------------------------------------
   ! type that defines the list of frequency depedent materials
   !------------------------------------------------------------------------------
   type, public :: FreqDepenMaterials_t
      type(FreqDepenMaterial_t), dimension(:), pointer :: Vols => NULL ()
      type(FreqDepenMaterial_t), dimension(:), pointer :: Surfs => NULL ()
      type(FreqDepenMaterial_t), dimension(:), pointer :: Lins => NULL ()
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nSurfs = 0
      integer(kind=4) :: nLins = 0
      integer(kind=4) :: nVols_max = 0
      integer(kind=4) :: nSurfs_max = 0
      integer(kind=4) :: nLins_max = 0
      integer(kind=4) :: n_c_max = 0 !cota superior
   end type FreqDepenMaterials_t
   !------------------------------------------------------------------------------
   ! Type for the ANISOTROPIC body, surface and lines since they will contain
   ! the same information
   !------------------------------------------------------------------------------
   type, public :: ANISOTROPICbody_t
      type(coords_t), dimension(:), pointer :: c1P => NULL ()
      type(coords_t), dimension(:), pointer :: c2P => NULL ()
      real(kind=RK), dimension(3, 3) :: sigma, eps, mu, sigmam
      integer(kind=4) :: n_C1P = 0
      integer(kind=4) :: n_C2P = 0
   end type ANISOTROPICbody_t
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
   end type ANISOTROPICelements_t
   !------------------------------------------------------------------------------
   ! Defines a Comp Surface
   !------------------------------------------------------------------------------
   type, public :: LossyThinSurface_t
      type(coords_t), dimension(:), pointer :: c => NULL ()
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
   end type LossyThinSurface_t
   !------------------------------------------------------------------------------
   ! Locates all the different Comp media found
   !------------------------------------------------------------------------------
   type, public :: LossyThinSurfaces_t
      type(LossyThinSurface_t), dimension(:), pointer :: cs => NULL ()
      integer(kind=4) :: length = 0
      integer(kind=4) :: length_max = 0
      integer(kind=4) :: nC_max = 0 !cota de todos los nc de LossyThinSurface
   end type LossyThinSurfaces_t
   !------------------------------------------------------------------------------
   ! Component for Thin Wires there is a list of this inside the component
   ! that defines the whole Thin Wire Reference
   !------------------------------------------------------------------------------
   type, public :: ThinWireComp_t
      character(len=BUFSIZE) :: srctype, srcfile
      integer(kind=4) :: i = - 1
      integer(kind=4) :: j = - 1
      integer(kind=4) :: K = - 1
      integer(kind=4) :: nd = - 1
      integer(kind=4) :: d = - 1
      real(kind=RK) :: m = 0.0_RKIND
      character(len=BUFSIZE) :: tag
   end type ThinWireComp_t
   !------------------------------------------------------------------------------
   ! ThinWire component that defines the overall properties of the definition
   ! of ThinWires
   !------------------------------------------------------------------------------
   type, public :: ThinWire_t
      type(ThinWireComp_t), dimension(:), pointer :: twc => NULL ()
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
   end type ThinWire_t
   !------------------------------------------------------------------------------
   ! List of the different thin wires that were found in the file
   !------------------------------------------------------------------------------
   type, public :: ThinWires_t
      type(ThinWire_t), dimension(:), pointer :: tw => NULL ()
      integer(kind=4) :: n_tw = 0
      integer(kind=4) :: n_tw_max = 0
   end type ThinWires_t
   !------------------------------------------------------------------------------
   ! Component for Slanted Wires there is a list of this inside the component
   ! that defines the whole Slanted Wire Reference
   !------------------------------------------------------------------------------
   type, public :: SlantedWireComp_t
      character(len=BUFSIZE) :: srctype, srcfile
      real(kind=RK) :: x = - 1.0_RKIND
      real(kind=RK) :: y = - 1.0_RKIND
      real(kind=RK) :: z = - 1.0_RKIND
      integer(kind=4) :: nd = - 1
      real(kind=RK) :: m = 0.0_RKIND
      character(len=BUFSIZE) :: tag
   end type SlantedWireComp_t
   !------------------------------------------------------------------------------
   ! ThinWire component that defines the overall properties of the definition
   ! of ThinWires
   !------------------------------------------------------------------------------
   type, public :: SlantedWire_t
      type(SlantedWireComp_t), dimension(:), pointer :: swc => NULL ()
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
   end type SlantedWire_t
   !------------------------------------------------------------------------------
   ! List of the different thin wires that were found in the file
   !------------------------------------------------------------------------------
   type, public :: SlantedWiresInfo_t
      type(SlantedWire_t), dimension(:), pointer :: sw => NULL ()
      integer(kind=4) :: n_sw = 0
      integer(kind=4) :: n_sw_max = 0
   end type
   !--------------------------------------------------------------------------
   ! Component for Thin Slots there is a list of this inside the component
   ! that defines the whole Thin Slot Reference
   !--------------------------------------------------------------------------
   type, public :: ThinSlotComp_t
      integer(kind=4) :: i = 0
      integer(kind=4) :: j = 0
      integer(kind=4) :: K = 0
      integer(kind=4) :: node = 0
      integer(kind=4) :: dir = - 1
      integer(kind=4) :: Or = - 1
      character(len=BUFSIZE) :: tag
   end type ThinSlotComp_t
   !--------------------------------------------------------------------------
   ! ThinSlot component that defines the overall properties of the definition
   ! of ThinSlots in ORIGINAL
   !--------------------------------------------------------------------------
   type, public :: ThinSlot_t
      type(ThinSlotComp_t), dimension(:), pointer :: tgc => NULL ()
      real(kind=RK) :: width = 0
      integer(kind=4) :: n_tgc = 0
      integer(kind=4) :: n_tgc_max = 0
   end type ThinSlot_t
   !--------------------------------------------------------------------------
   ! List of the different thin Slots that were found in the file
   !--------------------------------------------------------------------------
   type, public :: ThinSlots_t
      type(ThinSlot_t), dimension(:), pointer :: tg => NULL ()
      integer(kind=4) :: n_tg = 0
      integer(kind=4) :: n_tg_max = 0
   end type ThinSlots_t
   !-----------------> Border Types
   !------------------------------------------------------------------------------
   ! PML Border Type
   !------------------------------------------------------------------------------
   type, public :: FronteraPML_t
      real(kind=RK) :: orden = 2.0_RK
      real(kind=RK) :: refl = 1e-3_RK
      integer(kind=4) :: numCapas = 8
   end type FronteraPML_t
   !------------------------------------------------------------------------------
   ! Tipo de la frontera
   !------------------------------------------------------------------------------
   type, public :: Frontera_t
      integer(kind=4), dimension(6) :: tipoFrontera
      type(FronteraPML_t), dimension(6) :: propiedadesPML
   end type Frontera_t
   !-----------------> Probe Types
   !------------------------------------------------------------------------------
   ! type to define the new probe object which contains the type of calculation
   ! the type of analysis, time and frequency step and the filename where
   ! it should be saved
   !------------------------------------------------------------------------------
   type, public :: MasSonda_t
      character(len=BUFSIZE) :: filename
      type(coords_t), dimension(:), pointer :: cordinates => NULL ()
      real(kind=RK) :: tstart, tstop, tstep
      real(kind=RK) :: fstart, fstop, fstep
      integer(kind=4) :: type1, type2
      integer(kind=4) :: len_cor = 0
      character(len=BUFSIZE) :: outputrequest
   end type MasSonda_t
   !------------------------------------------------------------------------------
   ! type that defines a list of probes to be appended and accesed
   !------------------------------------------------------------------------------
   type, public :: MasSondas_t
      type(MasSonda_t), dimension(:), pointer :: collection => NULL ()
      integer(kind=4) :: length = 0
      integer(kind=4) :: length_max = 0
      integer(kind=4) :: len_cor_max = 0 !cota
   end type MasSondas_t
   !------------------------------------------------------------------------------
   ! This type contains the basic information in nearly all the different PROBES
   !------------------------------------------------------------------------------
   type, public :: Sonda_t
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
   end type Sonda_t
   !------------------------------------------------------------------------------
   ! type for the electric far field
   !------------------------------------------------------------------------------
   type, public :: FarField_Sonda_t
      type(Sonda_t) :: probe
   end type FarField_Sonda_t
   !------------------------------------------------------------------------------
   ! type for the electric field
   !------------------------------------------------------------------------------
   type, public :: Electric_Sonda_t
      type(Sonda_t) :: probe
   end type Electric_Sonda_t
   !------------------------------------------------------------------------------
   ! type for the magnetic field
   !------------------------------------------------------------------------------
   type, public :: Magnetic_Sonda_t
      type(Sonda_t) :: probe
   end type Magnetic_Sonda_t
   !------------------------------------------------------------------------------
   ! type for the normal electric field
   !------------------------------------------------------------------------------
   type, public :: NormalElectric_Sonda_t
      type(Sonda_t) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   end type NormalElectric_Sonda_t
   !------------------------------------------------------------------------------
   ! type for the normal magnetic field
   !------------------------------------------------------------------------------
   type, public :: NormalMagnetic_Sonda_t
      type(Sonda_t) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   end type NormalMagnetic_Sonda_t
   !------------------------------------------------------------------------------
   ! type for the electric surface current density
   !------------------------------------------------------------------------------
   type, public :: SurfaceElectricCurrent_Sonda_t
      type(Sonda_t) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   end type SurfaceElectricCurrent_Sonda_t
   !------------------------------------------------------------------------------
   ! type for the magnetic surface current density
   !------------------------------------------------------------------------------
   type, public :: SurfaceMagneticCurrent_Sonda_t
      type(Sonda_t) :: probe
      integer(kind=4), dimension(:), pointer :: nml => NULL ()
      integer(kind=4) :: n_nml = 0
      integer(kind=4) :: n_nml_max = 0
   end type SurfaceMagneticCurrent_Sonda_t
   !------------------------------------------------------------------------------
   ! Abstract class which performs the dynamic dispatching
   !------------------------------------------------------------------------------
   type, public :: abstractSonda_t
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
      type(FarField_Sonda_t), dimension(:), pointer :: FarField => NULL ()
      type(Electric_Sonda_t), dimension(:), pointer :: Electric => NULL ()
      type(Magnetic_Sonda_t), dimension(:), pointer :: Magnetic => NULL ()
      type(NormalElectric_Sonda_t), dimension(:), pointer :: NormalElectric => NULL ()
      type(NormalMagnetic_Sonda_t), dimension(:), pointer :: NormalMagnetic => NULL ()
      type(SurfaceElectricCurrent_Sonda_t), dimension(:), pointer :: SurfaceElectricCurrent => NULL ()
      type(SurfaceMagneticCurrent_Sonda_t), dimension(:), pointer :: SurfaceMagneticCurrent => NULL ()
   end type abstractSonda_t
   !------------------------------------------------------------------------------
   ! Class to account as a list for all the probes
   ! that might be required during the parsing process
   !------------------------------------------------------------------------------
   type, public :: Sondas_t
      type(abstractSonda_t), dimension(:), pointer :: probes => NULL ()
      integer(kind=4) :: n_probes = 0
      integer(kind=4) :: n_probes_max = 0
   end type Sondas_t
   !------------------------------------------------------------------------------
   ! Object type defined for the Bloque current probe
   !------------------------------------------------------------------------------
   type, public :: BloqueProbe_t
      real(kind=RK) :: tstart, tstop, tstep
      real(kind=RK) :: fstart, fstop, fstep
      character(len=BUFSIZE) :: FileNormalize
      integer(kind=4) :: type2
      integer(kind=4) :: i1, i2, j1, j2, k1, k2, skip
      integer(kind=4) :: nml
      LOGICAL :: t
      character(len=BUFSIZE) :: outputrequest
      character(len=BUFSIZE) :: tag
   end type BloqueProbe_t
   ! Object made for the collection of defined Bloque probes
   type, public :: BloqueProbes_t
      type(BloqueProbe_t), dimension(:), pointer :: bp => NULL ()
      integer(kind=4) :: n_bp = 0
      integer(kind=4) :: n_bp_max = 0
   end type BloqueProbes_t

   !------------------------------------------------------------------------------
   ! Object type defined for the Volumic probes
   !------------------------------------------------------------------------------
   type, public :: VolProbe_t
      type(coords_t), dimension(:), pointer :: cordinates => NULL ()
      real(kind=RK) :: tstart, tstop, tstep
      character(len=BUFSIZE) :: outputrequest
      integer(kind=4) :: len_cor = 0
      !para freq domain
      real(kind=RK) :: fstart, fstop, fstep
      integer(kind=4) :: type2
      character(len=BUFSIZE) :: filename
   end type VolProbe_t
   ! Object made for the collection of defined Volumic probes
   type, public :: VolProbes_t
      type(VolProbe_t), dimension(:), pointer :: collection => NULL ()
      integer(kind=4) :: length = 0
      integer(kind=4) :: length_max = 0
      integer(kind=4) :: len_cor_max = 0 !cota
   end type VolProbes_t

   !-----------------> Source Types
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: Box_t
      character(len=BUFSIZE) :: nombre_fichero
      integer(kind=4), dimension(3) :: coor1, coor2
   end type Box_t
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: Boxes_t
      type(Box_t), dimension(:), pointer :: Vols => NULL ()
      integer(kind=4) :: nVols = 0
      integer(kind=4) :: nVols_max = 0
   end type Boxes_t
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: PlaneWave_t
      character(len=BUFSIZE) :: nombre_fichero
      character(len=BUFSIZE) :: atributo
      integer(kind=4), dimension(3) :: coor1, coor2
      real(kind=RK) :: theta, phi, alpha, beta
      logical :: isRC !for reververation chambers
      real(kind=RK) :: INCERTMAX
      integer(kind=4) :: numModes !for reververation chambers
   end type PlaneWave_t
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type, public :: PlaneWaves_t
      type(PlaneWave_t), dimension(:), pointer :: collection => NULL ()
      integer(kind=4) :: nc = 0
      integer(kind=4) :: nC_max = 0
   end type PlaneWaves_t
   !------------------------------------------------------------------------------
   ! Definicin de los tipos current density que existirn en el ficero
   ! nfde
   !------------------------------------------------------------------------------
   type, public :: Curr_Field_Src_t
      type(coords_scaled_t), dimension(:), pointer :: c1P => NULL ()
      type(coords_scaled_t), dimension(:), pointer :: c2P => NULL ()
      character(len=BUFSIZE) :: nombre
      integer(kind=4) :: n_C1P = 0
      integer(kind=4) :: n_C2P = 0
      LOGICAL :: isElec, isHard, isInitialValue
   end type Curr_Field_Src_t
   !------------------------------------------------------------------------------
   ! Definicin de las Nodal Source global
   !------------------------------------------------------------------------------
   type, public :: NodSource_t
      type(Curr_Field_Src_t), dimension(:), pointer :: NodalSource => NULL ()
      integer(kind=4) :: n_nodSrc = 0
      integer(kind=4) :: n_nodSrc_max = 0
      integer(kind=4) :: n_C1P_max = 0
      integer(kind=4) :: n_C2P_max = 0
   end type NodSource_t
   !-----------------> General Types
   !------------------------------------------------------------------------------
   ! Matrix attributes.
   ! Total[XYZ] -> Is the cell number for each axis.
   !------------------------------------------------------------------------------
   type, public :: MatrizMedios_t
      integer(kind=4) :: totalX, totalY, totalZ
   end type MatrizMedios_t
   !------------------------------------------------------------------------------
   !------------------------------------------------------------------------------
   type NFDEGeneral_t
      real(kind=RK) :: dt
      integer(kind=4) :: nmax
      LOGICAL :: mtlnProblem
   end type NFDEGeneral_t
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
   type Desplazamiento_t
      real(kind=RK), dimension(:), pointer :: desX => NULL ()
      real(kind=RK), dimension(:), pointer :: desY => NULL ()
      real(kind=RK), dimension(:), pointer :: desZ => NULL ()
      integer(kind=4) :: nX = 0, mx1 = 0, mx2 = 0 !2012
      integer(kind=4) :: nY = 0, my1 = 0, my2 = 0 !2012
      integer(kind=4) :: nZ = 0, mz1 = 0, mz2 = 0 !2012
      real(kind=RK) ::originx= 0.0_RKIND  !2012
      real(kind=RK) ::originy= 0.0_RKIND  !2012
      real(kind=RK) ::originz= 0.0_RKIND  !2012
   end type Desplazamiento_t
   !-----------------> Program Types
   !------------------------------------------------------------------------------
   ! Parameters needed for the parser
   !------------------------------------------------------------------------------
   type, public :: Parseador_t
      character(len=BUFSIZE) :: switches=' '  
      ! Basics
      type(NFDEGeneral_t), pointer :: general => NULL ()
      type(MatrizMedios_t), pointer :: matriz => NULL ()
      type(Desplazamiento_t), pointer :: despl => NULL ()
      type(Frontera_t), pointer :: front => NULL ()
      ! Materials
      type(Materials_t), pointer :: Mats => NULL ()
      type(PECRegions_t), pointer :: pecRegs => NULL ()
      type(PECRegions_t), pointer :: pmcRegs => NULL ()
      type(DielectricRegions_t), pointer :: DielRegs => NULL ()
      type(LossyThinSurfaces_t), pointer :: LossyThinSurfs => NULL ()
      type(FreqDepenMaterials_t), pointer :: frqDepMats => NULL ()
      type(ANISOTROPICelements_t), pointer :: aniMats => NULL ()
      ! Sources
      type(Boxes_t), pointer :: boxSrc => NULL ()
      type(PlaneWaves_t), pointer :: plnSrc => NULL ()
      type(NodSource_t), pointer :: nodSrc => NULL ()
      ! Probes
      type(Sondas_t), pointer :: oldSONDA => NULL ()
      type(MasSondas_t), pointer :: Sonda => NULL ()
      type(BloqueProbes_t), pointer :: BloquePrb => NULL ()
      type(VolProbes_t), pointer :: VolPrb => NULL ()
      ! Thin Elements                         
      type(ThinWires_t), pointer :: tWires => NULL ()
      type(SlantedWiresInfo_t), pointer :: sWires => NULL ()
      type(ThinSlots_t), pointer :: tSlots => NULL ()
      ! Conformal
      type(ConformalPECRegions_t), pointer :: conformalRegs => NULL()
#ifdef CompileWithMTLN
      type(mtln_t), pointer :: mtln => NULL () 
#endif
   end type Parseador_t
   
   !---> definicion de tipos
   type, public :: t_linea_t
      integer(kind=4) :: LEN
      character(len=BUFSIZE) :: dato
   end type t_linea_t
   !--->
   type, public :: t_NFDE_FILE_t
      integer(kind=8) :: targ
      !--->
      integer(kind=8) :: numero
      type(t_linea_t), dimension(:), pointer :: lineas
      logical :: thereare_stoch
   end type t_NFDE_FILE_t
!--->

contains


end module NFDETypes_m

    
