
    
module lumped_vars
    !structures needed by the Lumped
   use fdetypes   
   
   type, public  :: Nodes_t
      real(kind=RKIND) :: EfieldPrev,EfieldPrevPrev,Jcur,sigmaEffResistInduct,alignedDeltaE ,transversalDeltaHa ,transversalDeltaHb, currentCoeff 
      real(kind=RKIND) :: diodepreA,diodeB
      real(kind=RKIND), pointer :: Efield,Ha_Plus,Ha_Minu,Hb_Plus,Hb_Minu
      real(kind=RKIND) :: g1
      real(kind=RKIND) :: g2a,g2b,GJ
      real(kind=RKIND) :: g1_usual
      real(kind=RKIND) :: g2a_usual,g2b_usual
      integer(kind=4) :: jmed, Orient !!! positivo o negativo...... 
!!!!for_devia 151222
#ifdef CompileWithStochastic
      real(kind=RKIND) :: EfieldPrev_for_devia,EfieldPrevPrev_for_devia,Jcur_for_devia,sigmaEffResistInduct_devia
      real(kind=RKIND) :: Efield_for_devia,Ha_Plus_for_devia,Ha_Minu_for_devia,Hb_Plus_for_devia,Hb_Minu_for_devia !no son punteros para stochastic sino valores que recibe desde mpi
      real(kind=RKIND) :: g1_devia
      real(kind=RKIND) :: g2a_devia,g2b_devia,GJ_devia
      real(kind=RKIND) :: g1_usual_devia
      real(kind=RKIND) :: g2a_usual_devia,g2b_usual_devia
#endif 
   END type Nodes_t


   type, public  :: LumpedElem_t
      integer(kind=4) :: NumNodes
      type(Nodes_t), allocatable, dimension(:) :: nodes
   end type LumpedElem_t

   

end module lumped_vars