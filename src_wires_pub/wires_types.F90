
 
module wiresHolland_constants
   use fdetypes
   
   !Types definitions
   
   integer(kind=4), parameter             :: MaxNumCurrentMinusPlus=9
   type, public  :: ChargeNodes_t
      integer(kind=4)                        :: IndexNode
      type(CurrentSegments_t), pointer         :: CurrentPlus_1,CurrentMinus_1  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_2,CurrentMinus_2  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_3,CurrentMinus_3  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_4,CurrentMinus_4  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_5,CurrentMinus_5  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_6,CurrentMinus_6  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_7,CurrentMinus_7  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_8,CurrentMinus_8  !neighbours in the plus and Minus direction
      type(CurrentSegments_t), pointer         :: CurrentPlus_9,CurrentMinus_9  !neighbours in the plus and Minus direction
      logical                                 :: IsMur,IsPeriodic,IsAttachedtoVoltage,IsPEC,HasIsource,Exists,&
                                                  Is_LeftEnd,Is_RightEnd,IsLossy
      logical                                 :: IsBackDownLeftMur,IsFrontUpRightMur
      logical                                         :: proc !dama
      logical                                 :: IsHeterogeneousJunction,IsInSingleRLCsegment   
      real(kind=RKIND_wires)                           :: cteMur,ctePlain,origctePlain,cteprop
      !to apply Mur. Needs extra storage everywhere but it is only 1D
      real(kind=RKIND_wires)                           :: ChargePresent,ChargePast
      type(ChargeNodes_t), pointer             :: NodeInside
      integer(kind=4)                        :: NumCurrentMinus,NumCurrentPlus
      integer(kind=4)                        :: i,j,k
      type(source_t), pointer                  :: Isource
      integer(kind=4), dimension(1:2*MaxNumCurrentMinusPlus) :: YESsegment
      
      real(kind=RKIND) , pointer             :: already_YEEadvanced_byconformal_changedtoPECfield1 => null()
      real(kind=RKIND) , pointer             :: already_YEEadvanced_byconformal_changedtoPECfield2 => null()
      real(kind=RKIND) , pointer             :: already_YEEadvanced_byconformal_changedtoPECfield3 => null()
      real(kind=RKIND) , pointer             :: already_YEEadvanced_byconformal_changedtoPECfield4 => null()
      real(kind=RKIND) , pointer             :: already_YEEadvanced_byconformal_changedtoPECfield5 => null()
      real(kind=RKIND) , pointer             :: already_YEEadvanced_byconformal_changedtoPECfield6 => null()
#ifdef CompileWithMPI
      !For MPI purposes !only handled and initialized in MPIcomm
      type(CurrentSegments_t), pointer         :: MPISharedCurrent
#endif

   end type ChargeNodes_t

#ifdef CompileWithThickWires    
   type container_t                
        real(kind=RKIND), pointer :: punt
        real(kind=RKIND) :: retardo      
        real(kind=RKIND), dimension(:), allocatable :: field_retard  
   end type container_t
   type :: thick_t    
      integer(kind=4)                        :: Enumero,Hnumero
      type(container_t), dimension(:), allocatable :: Efield_wire2main
      type(container_t), dimension(:), allocatable :: Hfield_wire2main, H_Efield_wire2main
      real(kind=RKIND_wires), dimension(:), allocatable :: EArea,rEArea,HArea,rHArea,rEfractionArea,Hsigno,Hcte  
      integer, dimension(:), allocatable :: i, j, k, field     
      logical :: Hplus         
      real(kind=RKIND), dimension(:), allocatable :: Current_ret    
      integer          :: maxretardo
   end type thick_t
#endif       
   type, public  :: CurrentSegments_t
      integer(kind=4)                        :: IndexSegment,NumParallel,OrigIndex
      type(wires_t), pointer              :: TipoWire
      real(kind=RKIND_wires)                    :: Lind,inv_Lind_acum,HEUR_safety,Lind_acum
      real(kind=RKIND_wires)                    :: delta,deltaTransv1,deltaTransv2
      real(kind=RKIND_wires)                    :: givenautoin, resist
      real(kind=RKIND_wires)                    :: givenautoin_devia, resist_devia
      type(ChargeNodes_t), pointer          :: ChargePlus ,ChargeMinus !neighbours in the plus and Minus direction
      logical                              :: IsPMC,HasVsource,IsShielded,HasParallel_RightEnd,HasParallel_LeftEnd, &
                                               HasSeries_RightEnd,HasSeries_LeftEnd,HasAbsorbing_RightEnd,HasAbsorbing_LeftEnd
      logical                              :: Is_LeftEnd, Is_RightEnd,IsEnd_norLeft_norRight,proc,IsConformal
      real(kind=RKIND_wires)                           :: cte1,cte2,cte3,cte5,FractionPlus,FractionMinus
      real(kind=RKIND_wires)                           :: Current,qplus_qminus
      real(kind=RKIND_wires)                           :: CurrentPast !added just for right observation
      !at the desired time step in observation.f90       
      real(kind=RKIND) , pointer                 :: Efield_wire2main,Efield_main2wire
#ifdef CompileWithThickWires      
      type(thick_t) :: thick
#endif
!      real(kind=RKIND_wires)                           :: Efield_wire2main_past  !no sirve para nada 171216
      integer(kind=4) :: i,j,k,indexmed,ILIBRE,JLIBRE,KLIBRE
      !dama
      integer(kind=4) :: ie,je,ke
      real(kind=RKIND_wires) :: x,y
      real    (kind=RKIND_wires)                            :: L, C, R
      real    (kind=RKIND_wires)                            :: L_devia, C_devia, R_devia
      real    (kind=RKIND_wires)                            :: cI
      real    (kind=RKIND_wires)                            :: bI
      real    (kind=RKIND_wires)                            :: Lintrinsic
      !fin dama
      integer(kind=4) :: tipofield !iEx,iEy o iEz
      logical :: orientadoalreves
      type(source_t), pointer                  :: Vsource
#ifdef CompileWithMPI
      !only required by the new MPI wires routines march'12 2012 bug multiwires MPI
      integer(kind=4) :: equivalentIndex
#endif

      !!!crank-nicolson coefficients
      real    (kind=RKIND_wires) :: upperdiag, diag, lowerdiag, rightCHminus, rightCHplus,rightCU,rightCUminus,rightCUplus
      !!!!!!!!!!end crank-nicolson
!!!se aniade siempre aunque solo lo use stochastic
      real(kind=RKIND_wires) :: qplus_qminus_for_devia,current_for_devia,Efield_main2wire_for_devia ,Lind_devia
      real(kind=RKIND_wires)                           :: cte1_for_devia ,cte2_for_devia ,cte3_for_devia  
   end type CurrentSegments_t
   !

   !dama
   type, public    :: TSegmentPtr_t
      type    (CurrentSegments_t)   , pointer                  :: ptr
   end type            TSegmentPtr_t

   type, public    :: TMultiline_t
      integer(kind=4)                                :: NumParallel
      type    (TSegmentPtr_t), pointer, dimension(:) :: Segments
      real    (kind=RKIND_wires) , pointer, dimension(:,:) :: R, L, C
      real    (kind=RKIND_wires) , pointer, dimension(:,:) :: b1I, b2I, b3I
   end type            TMultiline_t
   !!!!!!!!!!!!fin dama

   type, public   :: ThinWires_t
      integer(kind=4)                                :: NumMultilines !dama
      type    (TMultiline_t) , pointer, dimension(:) :: Multilines    !dama
      integer(kind=4) :: NumDifferentWires,NumCurrentSegments,NumChargeNodes
      integer(kind=4), pointer, dimension( : ) :: WireTipoMedio
      type(CurrentSegments_t) :: NullSegment !contiene informacion nula precisada por segmentos voided pero observados en la rutina de observacion 12/09/13
      type(ChargeNodes_t) :: NullNode
      type(CurrentSegments_t), pointer, dimension( : ) :: CurrentSegment
      type(ChargeNodes_t), pointer, dimension( : ) :: ChargeNode
#ifdef CompileWithMPI
      !For MPI purposes !only handled and initialized in MPIcomm
      type(CurrentSegments_t), pointer, dimension( : ) :: MPIUpNeededCurrentSegment,MPIDownNeededCurrentSegment
      integer(kind=4)                                 :: NumNeededCurrentUpMPI,NumNeededCurrentDownMPI
      type(ChargeNodes_t), pointer, dimension( : ) :: MPIUpChargeNode,MPIDownChargeNode
      !only required by the new MPI wires routines march'12 2012 bug multiwires MPI
      type(CurrentSegments_t), pointer, dimension( : ) :: MPIUpSharedCurrentSegment,MPIDownSharedCurrentSegment
      integer(kind=4)                                 :: NumSharedCurrentUpMPI,NumSharedCurrentDownMPI
#endif
      real(kind=RKIND)                   :: null_field !en los segmentos embeddeds y en los paralelos no hay acople entre thin-wire y medio
      real(kind=RKIND_wires)                   :: olddt !para permit scaling 141118
      ! apunto  a null_field el pointer field anterior en vez de al campo fdtd y lo obligo a ser cero
   end type ThinWires_t
   !
   type, public:: adyacc_t
      logical  :: Is, Parallel,IsHeterogeneousJunction,BothEndingsConnected
      integer(kind=4) :: i,j,k
      integer(kind=4), dimension(1:2) :: YESsegment
   end type
   !end type definitions

end module wiresHolland_constants