#include <vector>
#include <cstdint>

// Assuming FDETYPES_m provides RKIND definition. 
// In a real translation, this would be an include or a typedef.
// For this translation, we assume RKIND is defined as double.
#ifndef RKIND
#define RKIND double
#endif

namespace lumped_vars_m {

    struct Nodes_t {
        double EfieldPrev;
        double EfieldPrevPrev;
        double Jcur;
        double sigmaEffResistInduct;
        double alignedDeltaE;
        double transversalDeltaHa;
        double transversalDeltaHb;
        double currentCoeff;
        double diodepreA;
        double diodeB;
        
        // Pointers in Fortran become raw pointers in C++. 
        // Note: Memory management for these pointers must be handled externally.
        double* Efield;
        double* Ha_Plus;
        double* Ha_Minu;
        double* Hb_Plus;
        double* Hb_Minu;
        
        double g1;
        double g2a;
        double g2b;
        double GJ;
        double g1_usual;
        double g2a_usual;
        double g2b_usual;
        
        int32_t jmed;
        int32_t Orient; // positivo o negativo......

#ifdef CompileWithStochastic
        double EfieldPrev_for_devia;
        double EfieldPrevPrev_for_devia;
        double Jcur_for_devia;
        double sigmaEffResistInduct_devia;
        double Efield_for_devia;
        double Ha_Plus_for_devia;
        double Ha_Minu_for_devia;
        double Hb_Plus_for_devia;
        double Hb_Minu_for_devia; // no son punteros para stochastic sino valores que recibe desde mpi
        double g1_devia;
        double g2a_devia;
        double g2b_devia;
        double GJ_devia;
        double g1_usual_devia;
        double g2a_usual_devia;
        double g2b_usual_devia;
#endif
    };

    struct LumpedElem_t {
        int32_t NumNodes;
        std::vector<Nodes_t> nodes; // allocatable array becomes std::vector
    };

} // namespace lumped_vars_m