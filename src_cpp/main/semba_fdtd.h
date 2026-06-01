#ifndef SEMBA_FDTD_H
#define SEMBA_FDTD_H

#include <memory>
#include <string>
#include <vector>

namespace SEMBA_FDTD_m {

std::string extractCaseNameFromInput(const std::string& input_file);
std::string resolveInputFileFromFlags(const std::string& input_flags);

struct semba_fdtd_t {
    struct Impl;
    std::unique_ptr<Impl> impl_;
    bool finishedwithsuccess = false;

    semba_fdtd_t();
    ~semba_fdtd_t();

    void init(const std::string& input_flags = "");
    void launch();
    void end(const std::string& case_name);
};

namespace SEMBA_FDTD_test {

struct MurAbsorptionResult {
    double max_ex_initial = 0.0;
    double max_ex_final = 0.0;
    double probe_ex_initial = 0.0;
    double probe_ex_final = 0.0;
    double energy_initial = 0.0;
    double energy_final = 0.0;
};

struct BoundaryModeInfo {
    bool useMur = false;
    bool usePml = false;
    bool murBack = false, murFront = false, murLeft = false;
    bool murRight = false, murDown = false, murUp = false;
    bool pmlBack = false, pmlFront = false, pmlLeft = false;
    bool pmlRight = false, pmlDown = false, pmlUp = false;
    int pmlElectricCalls = 0;
    int pmlBodyHCalls = 0;
    int pmlMagneticCpmlCalls = 0;
};

struct MpiSliceInfo {
    int rank = 0;
    int ranks = 1;
    int axis = 3;
    int com = 0;
    int fin = 0;
    int sweepZI = 0;
    int sweepZE = 0;
    int allocZI = 0;
    int allocZE = 0;
    bool physicalDown = true;
    bool physicalUp = true;
    bool pmlDown = false;
    bool pmlUp = false;
};

struct PlaneWaveInitInfo {
    double px = 0.0, py = 0.0, pz = 0.0;
    double ex = 0.0, ey = 0.0, ez = 0.0;
    double hx = 0.0, hy = 0.0, hz = 0.0;
    double distanciaInicial = 0.0;
    double dt = 0.0;
    int numSteps = 0;
    int esqx1 = 0, esqx2 = 0, esqy1 = 0, esqy2 = 0, esqz1 = 0, esqz2 = 0;
    bool iluminaAb = false, iluminaAr = false;
    double murCx = 0.0, murCy = 0.0, murCz = 0.0;
    double deltaevol = 0.0;
    int numSamples = 0;
};

int run_init_solver_test(const std::string& json_path);
double test_evolucion(const std::string& json_path, int pw_idx, double t_delay);
double test_compute_incid(const std::string& json_path, int pw_idx, int nfield,
                          double time, int i, int j, int k);
double test_grid_inverse_z(const std::string& json_path, int k);
PlaneWaveInitInfo test_plane_wave_init(const std::string& json_path, int pw_idx);
BoundaryModeInfo test_boundary_mode(const std::string& json_path,
                                    bool step_once = false);
int test_mpi_axis_from_flags(const std::string& flags);
std::vector<MpiSliceInfo> test_mpi_one_axis_slices(int cells,
                                                   int ranks,
                                                   int pml_down_layers = 0,
                                                   int pml_up_layers = 0,
                                                   int forced_cut = -1,
                                                   int axis = 3);
int test_mpi_exchange_electric_ghost_planes(int axis);
int test_mpi_exchange_magnetic_ghost_planes(int axis);
double test_mur_apply_back_hy(const std::string& json_path);
MurAbsorptionResult test_mur_pulse_absorption(const std::string& json_path,
                                              int num_steps,
                                              int pulse_i, int pulse_j, int pulse_k,
                                              double amplitude,
                                              bool apply_mur = true);
double test_field_after_tfsf_e_step(const std::string& json_path, int component, int i, int j, int k);
int test_run_pw_in_box_probes(const std::string& json_path,
                              const std::string& ref_before,
                              const std::string& ref_inbox,
                              const std::string& ref_after,
                              int max_steps = -1);
int test_run_pw_in_box_probe_files_exact(const std::string& json_path,
                                         const std::string& ref_before,
                                         const std::string& ref_inbox,
                                         const std::string& ref_after,
                                         int max_steps = -1);
int test_run_holland_probe_output(const std::string& json_path,
                                  int max_steps = 10);
int test_run_bulk_current_probe_output(const std::string& json_path,
                                       const std::string& expected_name,
                                       int max_steps = 5);

} // namespace SEMBA_FDTD_test
} // namespace SEMBA_FDTD_m

#endif
