
    
module lumped_vars
    !structures needed by the Lumped
   use fdetypes   
   
   TYPE, public  ::  Nodes_t
      REAL (KIND=RKIND) :: EfieldPrev,EfieldPrevPrev,Jcur,sigmaEffResistInduct,alignedDeltaE ,transversalDeltaHa ,transversalDeltaHb 
      REAL (KIND=RKIND) :: diodepreA,diodeB
      REAL (KIND=RKIND), pointer ::  Efield,Ha_Plus,Ha_Minu,Hb_Plus,Hb_Minu
      REAL (KIND=RKIND)          ::  g1
      REAL (KIND=RKIND)          ::  g2a,g2b,GJ
      REAL (KIND=RKIND)          ::  g1_usual
      REAL (KIND=RKIND)          ::  g2a_usual,g2b_usual
      integer (kind=4) :: jmed, Orient !!! positivo o negativo...... 
!!!!for_devia 151222
#ifdef CompileWithStochastic
      REAL (KIND=RKIND) :: EfieldPrev_for_devia,EfieldPrevPrev_for_devia,Jcur_for_devia,sigmaEffResistInduct_devia
      REAL (KIND=RKIND) ::  Efield_for_devia,Ha_Plus_for_devia,Ha_Minu_for_devia,Hb_Plus_for_devia,Hb_Minu_for_devia !no son punteros para stochastic sino valores que recibe desde mpi
      REAL (KIND=RKIND)          ::  g1_devia
      REAL (KIND=RKIND)          ::  g2a_devia,g2b_devia,GJ_devia
      REAL (KIND=RKIND)          ::  g1_usual_devia
      REAL (KIND=RKIND)          ::  g2a_usual_devia,g2b_usual_devia
#endif 
   END TYPE Nodes_t


   TYPE, public  ::  LumpedElem_t
      integer (kind=4)   ::   NumNodes
      type (Nodes_t), allocatable, dimension (:) :: nodes
   end type LumpedElem_t

   

end module lumped_vars