// Minimal VTK output stubs for semba-outputs until vtk.cpp is fully translated.
#include <string>
#include <vector>

struct SGGFDTDINFO_t;

namespace Observa_m {
struct SGGFDTDINFO_t;
}

void creaUnstructData(int, int, const Observa_m::SGGFDTDINFO_t&, int&, int&, int&, int&, bool) {}

void write_VTKfile(const Observa_m::SGGFDTDINFO_t&, const std::string&, int, int, int, int, int, int,
                   double, bool) {}

void createVTK(int layoutnumber, int num_procs, Observa_m::SGGFDTDINFO_t& sgg, bool vtkindex,
               bool& somethingdone, int mpidir,
               const std::vector<std::vector<std::vector<int>>>& sggMtag, bool dontwritevtk) {
    (void)layoutnumber;
    (void)num_procs;
    (void)sgg;
    (void)vtkindex;
    (void)somethingdone;
    (void)mpidir;
    (void)sggMtag;
    (void)dontwritevtk;
}

void createVTKOnTheFly(int layoutnumber, int num_procs, Observa_m::SGGFDTDINFO_t& sgg, bool vtkindex,
                       bool& somethingdone, int mpidir,
                       const std::vector<std::vector<std::vector<int>>>& sggMtag,
                       bool dontwritevtk) {
    createVTK(layoutnumber, num_procs, sgg, vtkindex, somethingdone, mpidir, sggMtag, dontwritevtk);
}
