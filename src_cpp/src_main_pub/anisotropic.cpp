#include <vector>
#include <cmath>
#include <memory>

// Forward declarations for types defined in other modules
struct SGGFDTDINFO_t {
        std::vector<double> tiempo;
        double dt = 0.0;
        int NumMedia = 0;
        struct { double eps=0.0, mur=0.0; int Is=0; } Med[1];
        struct { double eps=0.0, mur=0.0; int Is=0; } med[1];
    };
struct media_matrices_t;

namespace ANISOTROPIC_m {
    struct media_matrices_t { int NumMed=0; struct { int indexmed=0; bool IsOnlyThinSlot=false; double sigma[3][3]={}; } info[1]; };
    void InitAnisotropic(const int&, const media_matrices_t&, bool&, bool&, double, double) {}
    void DestroyAnisotropic(int&) {}
    void calc_anisotropicconstants(const int&, double&, double&) {}
    void DestroyAnisotropicInfo(media_matrices_t&) {}
    void InitAnisotropicInfo(media_matrices_t&, int) {}
} // namespace ANISOTROPIC_m
