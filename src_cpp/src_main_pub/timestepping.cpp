#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>

// Forward declarations and includes for external modules/types would go here
// Assuming these types exist in the translated headers
struct sim_control_t;
struct logic_control_t;
struct perform_t;
struct constants_t;
struct bounds_t;
struct EpsMuTimeScale_input_parameters_t;
struct SGGFDTDINFO_t;
struct media_matrices_t;
struct taglist_t;
struct limit_t;
struct tagtype_t;
struct entrada_t;
struct mtln_t;

// Constants and types assumed from FDETYPES_m and others
using rkind = double;
using RKIND_tiempo = double;
using RKIND = double;
using bufsize = 256; // Example size

// Placeholder for external functions/types that are not fully defined in the snippet
// In a real translation, these would be proper classes/structs from their respective modules
namespace ExternalModules {
    void stoponerror(int layout, int procs, const std::string& msg, bool fatal = false);
    void reportmedia(const SGGFDTDINFO_t& sgg);
    void findbounds(bounds_t& bounds);
    void print11(int layout, const std::string& msg);
    void ReadFields(const SGGFDTDINFO_t& alloc, int& laststep, double& lasttime, double& ultimodt, double& eps0, double& mu0, 
                    std::vector<std::vector<std::vector<rkind>>>& Ex, std::vector<std::vector<std::vector<rkind>>>& Ey, 
                    std::vector<std::vector<std::vector<rkind>>>& Ez, std::vector<std::vector<std::vector<rkind>>>& Hx, 
                    std::vector<std::vector<std::vector<rkind>>>& Hy, std::vector<std::vector<std::vector<rkind>>>& Hz);
    void MPIupdateMin(double val, double& min_val);
    void MPI_AllReduce(int* val, int* result, int count, int type, int op, int comm, int& ierr);
    void MPIupdateMin(double val, double& min_val); // Overload for double
    void MPIupdateMin(int val, int& min_val); // Overload for int
    void MPIupdateMin(double val, double& min_val); // Overload for double
    void crea_timevector(const SGGFDTDINFO_t& sgg, int laststep, int finalstep, double lasttime);
    void updateSigmaM(bool& attinformado);
    void updateThinWiresSigma(bool& attinformado);
    void calc_G1G2Gm1Gm2(const SGGFDTDINFO_t& sgg, constants_t& g, double eps0, double mu0);
    void revertThinWiresSigma();
    void InitReporting(const SGGFDTDINFO_t& sgg, sim_control_t* control);
    void reportSimulationOptions();
    void initializeBorders();
    void initializeLumped();
    void initializeWires();
    void initializeAnisotropic();
    void initializeSGBC();
    void initializeMultiports();
    void initializeEDispersives();
    void initializeMDispersives();
    void initializePlanewave();
    void initializeNodalSources();
    void fillMtag(const SGGFDTDINFO_t& sgg, const std::vector<std::vector<std::vector<int>>>& MiEx, 
                  const std::vector<std::vector<std::vector<int>>>& MiEy, const std::vector<std::vector<std::vector<int>>>& MiEz,
                  const std::vector<std::vector<std::vector<int>>>& MiHx, const std::vector<std::vector<std::vector<int>>>& MiHy, 
                  const std::vector<std::vector<std::vector<int>>>& MiHz, std::vector<std::vector<std::vector<int>>>& Mtag, 
                  const bounds_t& bounds, const taglist_t& tag_numbers);
    void initializeObservation();
    void initializeMPI();
    void InitTiming(const SGGFDTDINFO_t& sgg, sim_control_t* control, double time_start, int initial_step, double max_source);
    void CLOSEWARNINGFILE(int layout, int procs, bool& fatal_error, bool write, bool simu_devia);
    void solveMTLNProblem(const mtln_t& mtln_parsed, const std::string& nEntradaRoot);
    void reportSimulationEnd(int layoutnumber);
}

// Placeholder for MPI constants
#ifdef CompileWithMPI
    int SUBCOMM_MPI = 0;
    int MPI_INTEGER = 0;
    int MPI_MIN = 0;
    int MPI_MAX = 0;
#endif

namespace Solver_m {

    class solver_t {
    public:
        sim_control_t control;
        logic_control_t thereAre;
        perform_t perform;
        perform_t d_perform;

        // Using 3D vectors for Ex, Ey, Ez, Hx, Hy, Hz
        std::vector<std::vector<std::vector<rkind>>> Ex, Ey, Ez;
        std::vector<std::vector<std::vector<rkind>>> Hx, Hy, Hz;

        // Using 1D vectors for Idxe, Idye, etc.
        std::vector<rkind> Idxe, Idye, Idze;
        std::vector<rkind> Idxh, Idyh, Idzh;
        std::vector<rkind> dxe, dye, dze;
        std::vector<rkind> dxh, dyh, dzh;

        constants_t g;
        RKIND_tiempo lastexecutedtime;
        rkind maxSourceValue;

        int initialtimestep;
        int lastexecutedtimestep;
        int ini_save;
        int n_info;
        int n;

        bounds_t bounds;
        EpsMuTimeScale_input_parameters_t EpsMuTimeScale_input_parameters;

        bool parar;
        bool everflushed;
        bool still_planewave_time;

        // semba variables
        SGGFDTDINFO_t sgg;
        media_matrices_t media;
        taglist_t tag_numbers;
        std::vector<limit_t> SINPML_fullsize;
        std::vector<limit_t> fullsize;
        bool finishedwithsuccess;
        rkind eps0;
        rkind mu0;
        tagtype_t tagtype;

#ifdef CompileWithMTLN
        mtln_t mtln_parsed;
#endif

        // Constructor
        static solver_t create(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, const taglist_t& tag_numbers,
                               const std::vector<limit_t>& SINPML_Fullsize, const std::vector<limit_t>& fullsize,
                               bool finishedwithsuccess, rkind eps0, rkind mu0, const tagtype_t& tagtype,
                               const entrada_t& input, rkind maxSourceValue, double time_desdelanzamiento) {
            solver_t res;
            res.init_control(input, maxSourceValue, time_desdelanzamiento);
            res.sgg = sgg;
            res.media = media;
            res.tag_numbers = tag_numbers;
            res.SINPML_fullsize = SINPML_Fullsize;
            res.fullsize = fullsize;
            res.eps0 = eps0;
            res.mu0 = mu0;
            res.tagtype = tagtype;
            return res;
        }

        void init_control(const entrada_t& input, rkind maxSourceValue, double time_desdelanzamiento) {
            control.maxSourceValue = maxSourceValue;
            control.time_desdelanzamiento = time_desdelanzamiento;

            control.simu_devia = input.simu_devia;
            control.resume = input.resume;
            control.saveall = input.saveall;
            control.makeholes = input.makeholes;
            control.connectendings = input.connectendings;
            control.isolategroupgroups = input.isolategroupgroups;
            control.createmap = input.createmap;
            control.groundwires = input.groundwires;
            control.noSlantedcrecepelo = input.noSlantedcrecepelo;
            control.SGBC = input.SGBC;
            control.SGBCDispersive = input.SGBCDispersive;
            control.mibc = input.mibc;
            control.ADE = input.ADE;
            control.conformalskin = input.conformalskin;
            control.NOcompomur = input.NOcompomur;
            control.strictOLD = input.strictOLD;
            control.TAPARRABOS = input.TAPARRABOS;
            control.noconformalmapvtk = input.noconformalmapvtk;
            control.experimentalVideal = input.experimentalVideal;
            control.forceresampled = input.forceresampled;
            control.mur_second = input.mur_second;
            control.MurAfterPML = input.MurAfterPML;
            control.stableradholland = input.stableradholland;
            control.singlefilewrite = input.singlefilewrite;
            control.NF2FFDecim = input.NF2FFDecim;
            control.sgbccrank = input.sgbccrank;
            control.fieldtotl = input.fieldtotl;
            control.permitscaling = input.permitscaling;
            control.mtlnberenger = input.mtlnberenger;
            control.niapapostprocess = input.niapapostprocess;
            control.stochastic = input.stochastic;
            control.verbose = input.verbose;
            control.dontwritevtk = input.dontwritevtk;
            control.resume_fromold = input.resume_fromold;
            control.vtkindex = input.vtkindex;
            control.createh5bin = input.createh5bin;
            control.wirecrank = input.wirecrank;
            control.fatalerror = input.fatalerror;

            control.cfl = input.cfl;
            control.attfactorc = input.attfactorc;
            control.attfactorw = input.attfactorw;
            control.alphamaxpar = input.alphamaxpar;
            control.alphaOrden = input.alphaOrden;
            control.kappamaxpar = input.kappamaxpar;
            control.mindistwires = input.mindistwires;
            control.sgbcFreq = input.sgbcFreq;
            control.sgbcresol = input.sgbcresol;
            control.factorradius = input.factorradius;
            control.factordelta = input.factordelta;
            control.nEntradaRoot = input.nEntradaRoot;
            std::transform(control.nEntradaRoot.begin(), control.nEntradaRoot.end(), control.nEntradaRoot.begin(), ::toupper); // trim(adjustl()) equivalent roughly
            control.inductance_model = input.inductance_model;
            std::transform(control.inductance_model.begin(), control.inductance_model.end(), control.inductance_model.begin(), ::toupper);
            control.wiresflavor = input.wiresflavor;
            std::transform(control.wiresflavor.begin(), control.wiresflavor.end(), control.wiresflavor.begin(), ::toupper);
            control.nresumeable2 = input.nresumeable2;
            std::transform(control.nresumeable2.begin(), control.nresumeable2.end(), control.nresumeable2.begin(), ::toupper);
            control.opcionestotales = input.opcionestotales;
            control.finaltimestep = input.finaltimestep;
            control.flushsecondsFields = input.flushsecondsFields;
            control.flushsecondsData = input.flushsecondsData;
            control.layoutnumber = input.layoutnumber;
            control.mpidir = input.mpidir;
            control.inductance_order = input.inductance_order;
            control.wirethickness = input.wirethickness;
            control.maxCPUtime = input.maxCPUtime;
            control.SGBCDepth = input.SGBCDepth;
            control.precision = input.precision;
            control.num_procs = input.num_procs;
            control.MEDIOEXTRA = input.MEDIOEXTRA;
            control.facesNF2FF = input.facesNF2FF;
            control.EpsMuTimeScale_input_parameters = input.EpsMuTimeScale_input_parameters;

            thereAre.reset();
        }

#ifdef CompileWithMTLN
        void launch_mtln_simulation(const mtln_t& mtln_parsed, const std::string& nEntradaRoot, int layoutnumber) {
            ExternalModules::solveMTLNProblem(mtln_parsed, nEntradaRoot);
            ExternalModules::reportSimulationEnd(layoutnumber);
        }
#endif

        void init_fields() {
            // Allocate Ex, Ey, Ez, Hx, Hy, Hz based on sgg dimensions
            // Note: In Fortran, this uses lower bounds. In C++, we use 0-based indexing.
            // We assume sgg.Alloc(iEx)%XI etc. provide the size or offset.
            // For simplicity, we assume the allocated size is (XE-XI+1) etc.
            
            int xi_ex = sgg.Alloc(0).XI; // Assuming iEx=0, iEy=1, etc. Adjust indices as needed
            int xe_ex = sgg.Alloc(0).XE;
            int yi_ex = sgg.Alloc(0).YI;
            int ye_ex = sgg.Alloc(0).YE;
            int zi_ex = sgg.Alloc(0).ZI;
            int ze_ex = sgg.Alloc(0).ZE;

            int size_x = xe_ex - xi_ex + 1;
            int size_y = ye_ex - yi_ex + 1;
            int size_z = ze_ex - zi_ex + 1;

            Ex.assign(size_x, std::vector<std::vector<rkind>>(size_y, std::vector<rkind>(size_z, 0.0)));
            Ey.assign(size_x, std::vector<std::vector<rkind>>(size_y, std::vector<rkind>(size_z, 0.0)));
            Ez.assign(size_x, std::vector<std::vector<rkind>>(size_y, std::vector<rkind>(size_z, 0.0)));
            Hx.assign(size_x, std::vector<std::vector<rkind>>(size_y, std::vector<rkind>(size_z, 0.0)));
            Hy.assign(size_x, std::vector<std::vector<rkind>>(size_y, std::vector<rkind>(size_z, 0.0)));
            Hz.assign(size_x, std::vector<std::vector<rkind>>(size_y, std::vector<rkind>(size_z, 0.0)));
        }

        void init_distances() {
            // Allocate 1D vectors based on sgg dimensions
            // Similar logic to init_fields but for 1D arrays
            
            // dxe, Idxe associated with iHx
            int xi_hx = sgg.Alloc(3).XI; // Assuming iHx=3
            int xe_hx = sgg.Alloc(3).XE;
            int size_x_hx = xe_hx - xi_hx + 1;
            dxe.resize(size_x_hx);
            Idxe.resize(size_x_hx);

            // dye, Idye associated with iHy
            int yi_hy = sgg.Alloc(4).YI; // Assuming iHy=4
            int ye_hy = sgg.Alloc(4).YE;
            int size_y_hy = ye_hy - yi_hy + 1;
            dye.resize(size_y_hy);
            Idye.resize(size_y_hy);

            // dze, Idze associated with iHz
            int zi_hz = sgg.Alloc(5).ZI; // Assuming iHz=5
            int ze_hz = sgg.Alloc(5).ZE;
            int size_z_hz = ze_hz - zi_hz + 1;
            dze.resize(size_z_hz);
            Idze.resize(size_z_hz);

            // dxh, Idxh associated with iEx
            int xi_ex = sgg.Alloc(0).XI;
            int xe_ex = sgg.Alloc(0).XE;
            int size_x_ex = xe_ex - xi_ex + 1;
            dxh.resize(size_x_ex);
            Idxh.resize(size_x_ex);

            // dyh, Idyh associated with iEy
            int yi_ey = sgg.Alloc(1).YI;
            int ye_ey = sgg.Alloc(1).YE;
            int size_y_ey = ye_ey - yi_ey + 1;
            dyh.resize(size_y_ey);
            Idyh.resize(size_y_ey);

            // dzh, Idzh associated with iEz
            int zi_ez = sgg.Alloc(2).ZI;
            int ze_ez = sgg.Alloc(2).ZE;
            int size_z_ez = ze_ez - zi_ez + 1;
            dzh.resize(size_z_ez);
            Idzh.resize(size_z_ez);

            // Initialize with -1.0e10
            std::fill(dxe.begin(), dxe.end(), -1.0e10);
            std::fill(dye.begin(), dye.end(), -1.0e10);
            std::fill(dze.begin(), dze.end(), -1.0e10);
            std::fill(dxh.begin(), dxh.end(), -1.0e10);
            std::fill(dyh.begin(), dyh.end(), -1.0e10);
            std::fill(dzh.begin(), dzh.end(), -1.0e10);

            // Fill dxe, dye, dze from sgg.DX, DY, DZ
            for (int i = 0; i < size_x_hx; ++i) {
                dxe[i] = sgg.DX(xi_hx + i);
            }
            for (int i = 0; i < size_y_hy; ++i) {
                dye[i] = sgg.DY(yi_hy + i);
            }
            for (int i = 0; i < size_z_hz; ++i) {
                dze[i] = sgg.DZ(zi_hz + i);
            }

            // Fill dxh, dyh, dzh as averages
            for (int i = 0; i < size_x_ex; ++i) {
                dxh[i] = (sgg.DX(xi_ex + i) + sgg.DX(xi_ex + i - 1)) / 2.0;
            }
            for (int i = 0; i < size_y_ey; ++i) {
                dyh[i] = (sgg.DY(yi_ey + i) + sgg.DY(yi_ey + i - 1)) / 2.0;
            }
            for (int i = 0; i < size_z_ez; ++i) {
                dzh[i] = (sgg.DZ(zi_ez + i) + sgg.DZ(zi_ez + i - 1)) / 2.0;
            }

            // Calculate inverses
            for (int i = 0; i < size_x_hx; ++i) Idxe[i] = 1.0 / dxe[i];
            for (int i = 0; i < size_y_hy; ++i) Idye[i] = 1.0 / dye[i];
            for (int i = 0; i < size_z_hz; ++i) Idze[i] = 1.0 / dze[i];
            for (int i = 0; i < size_x_ex; ++i) Idxh[i] = 1.0 / dxh[i];
            for (int i = 0; i < size_y_ey; ++i) Idyh[i] = 1.0 / dyh[i];
            for (int i = 0; i < size_z_ez; ++i) Idzh[i] = 1.0 / dzh[i];
        }

        void set_field_value(int field_idx, const std::vector<int>& i_range, const std::vector<int>& j_range, const std::vector<int>& k_range, rkind field_value) {
            std::vector<std::vector<std::vector<rkind>>>* field_ptr = nullptr;
            switch (field_idx) {
                case 0: field_ptr = &Ex; break; // iEx
                case 1: field_ptr = &Ey; break; // iEy
                case 2: field_ptr = &Ez; break; // iEz
                case 3: field_ptr = &Hx; break; // iHx
                case 4: field_ptr = &Hy; break; // iHy
                case 5: field_ptr = &Hz; break; // iHz
                default: return;
            }

            for (int i = i_range[0]; i <= i_range[1]; ++i) {
                for (int j = j_range[0]; j <= j_range[1]; ++j) {
                    for (int k = k_range[0]; k <= k_range[1]; ++k) {
                        (*field_ptr)[i][j][k] = field_value;
                    }
                }
            }
        }

        rkind get_field_value(int field_idx, int fi, int fj, int fk) {
            std::vector<std::vector<std::vector<rkind>>>* field_ptr = nullptr;
            switch (field_idx) {
                case 0: field_ptr = &Ex; break;
                case 1: field_ptr = &Ey; break;
                case 2: field_ptr = &Ez; break;
                case 3: field_ptr = &Hx; break;
                case 4: field_ptr = &Hy; break;
                case 5: field_ptr = &Hz; break;
                default: return 0.0;
            }
            return (*field_ptr)[fi][fj][fk];
        }

        void launch_simulation() {
            init();
            run();
            end();
        }

        void init() {
            // This is a simplified version of solver_init due to length and external dependencies
            // In a real translation, all external function calls would be present
            
            control.fatalerror = false;
            parar = false;
            perform.reset();
            d_perform.reset();
            thereAre.reset();
            thereAre.MagneticMedia = sgg.thereareMagneticMedia;
            thereAre.PMLMagneticMedia = sgg.therearePMLMagneticMedia;

            // Prechecking offsets
            int I = sgg.Alloc(0).XI;
            int J = sgg.Alloc(0).YI;
            int K = sgg.Alloc(0).ZI;
            for (int field = 1; field <= 5; ++field) {
                if (sgg.Alloc(field).XI != I) ExternalModules::stoponerror(control.layoutnumber, control.num_procs, "OFFSETS IN INITIAL COORD NOT ALLOWED");
                if (sgg.Alloc(field).YI != J) ExternalModules::stoponerror(control.layoutnumber, control.num_procs, "OFFSETS IN INITIAL COORD NOT ALLOWED");
                if (sgg.Alloc(field).ZI != K) ExternalModules::stoponerror(control.layoutnumber, control.num_procs, "OFFSETS IN INITIAL COORD NOT ALLOWED");
            }

            // File names and bounds
            std::string whoami = "(" + std::to_string(control.layoutnumber + 1) + "/" + std::to_string(control.num_procs) + ") ";
            std::string chari = std::to_string(control.layoutnumber + 1);
            if ((control.layoutnumber == 0) && control.verbose) ExternalModules::reportmedia(sgg);
            std::string layoutcharID = control.nentradaroot + "_" + chari;
            ExternalModules::findbounds(bounds);

            init_distances();

            // Allocate g matrices
            g.g1.assign(sgg.NumMedia + 1, 0.0);
            g.g2.assign(sgg.NumMedia + 1, 0.0);
            g.gm1.assign(sgg.NumMedia + 1, 0.0);
            g.gm2.assign(sgg.NumMedia + 1, 0.0);

            init_fields();

            // Initialize local variables and observation stuff
            // dt0 = sgg.dt; // Assuming dt0 is a member or global, not shown in struct
            if (!control.resume) {
                std::fill(Ex[0][0].begin(), Ex[0][0].end(), 0.0);
                std::fill(Ey[0][0].begin(), Ey[0][0].end(), 0.0);
                std::fill(Ez[0][0].begin(), Ez[0][0].end(), 0.0);
                std::fill(Hx[0][0].begin(), Hx[0][0].end(), 0.0);
                std::fill(Hy[0][0].begin(), Hy[0][0].end(), 0.0);
                std::fill(Hz[0][0].begin(), Hz[0][0].end(), 0.0);
                initialtimestep = 0;
                lastexecutedtimestep = 0;
                lastexecutedtime = 0.0;
            } else {
                std::string dubuf = "Init processing resuming data";
                ExternalModules::print11(control.layoutnumber, dubuf);
                
                std::string filename = control.nresumeable2;
                if (control.resume_fromold) {
                    filename += ".old";
                }
                
                // Open file and read fields (simplified)
                // In real code, this would use std::fstream with binary mode
                // For now, we assume ExternalModules::ReadFields handles the file I/O
                double ultimodt = 0.0;
                ExternalModules::ReadFields(sgg.alloc, lastexecutedtimestep, lastexecutedtime, ultimodt, eps0, mu0, Ex, Ey, Ez, Hx, Hy, Hz);
                sgg.dt = ultimodt;

#ifdef CompileWithMPI
                double rdummy = sgg.dt;
                ExternalModules::MPIupdateMin(sgg.dt, rdummy);
                rdummy = eps0;
                ExternalModules::MPIupdateMin(eps0, rdummy);
                rdummy = mu0;
                ExternalModules::MPIupdateMin(mu0, rdummy);
#endif

#ifdef CompileWithMPI
                int dummyMin, dummyMax, ierr;
                ExternalModules::MPI_AllReduce(&lastexecutedtimestep, &dummyMin, 1, MPI_INTEGER, MPI_MIN, SUBCOMM_MPI, ierr);
                ExternalModules::MPI_AllReduce(&lastexecutedtimestep, &dummyMax, 1, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr);
                
                if ((dummyMax != lastexecutedtimestep) || (dummyMin != lastexecutedtimestep)) {
#ifdef CompileWithOldSaving
                    if (control.resume_fromold) {
                        ExternalModules::stoponerror(control.layoutnumber, control.num_procs, "Incoherence between MPI saved steps for resuming.", true);
                        destroy_and_deallocate();
                        return;
                    } else {
                        std::string dubuf = "Incoherence between MPI saved steps for resuming. Retrying with -old....";
                        ExternalModules::print11(control.layoutnumber, dubuf);
                        control.resume_fromold = true;
                        filename = control.nresumeable2 + ".old";
                        ExternalModules::ReadFields(sgg.alloc, lastexecutedtimestep, lastexecutedtime, ultimodt, eps0, mu0, Ex, Ey, Ez, Hx, Hy, Hz);
                        sgg.dt = ultimodt;
                        ExternalModules::MPI_AllReduce(&lastexecutedtimestep, &dummyMin, 1, MPI_INTEGER, MPI_MIN, SUBCOMM_MPI, ierr);
                        ExternalModules::MPI_AllReduce(&lastexecutedtimestep, &dummyMax, 1, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr);
                        if ((dummyMax != lastexecutedtimestep) || (dummyMin != lastexecutedtimestep)) {
                            std::string DUbuf = "NO success. fields.old MPI are also incoherent for resuming.";
                            ExternalModules::stoponerror(control.layoutnumber, control.num_procs, DUbuf, true);
                            destroy_and_deallocate();
                            return;
                        } else {
                            std::string dubuf = "SUCCESS: Restarting from .fields.old instead. From n=" + std::to_string(lastexecutedtimestep);
                            ExternalModules::print11(control.layoutnumber, dubuf);
                        }
                    }
#else
                    std::string dubuf = "Incoherence between MPI saved steps for resuming.";
                    ExternalModules::stoponerror(control.layoutnumber, control.num_procs, dubuf, true);
                    destroy_and_deallocate();
                    return;
#endif
                }
#endif
                initialtimestep = lastexecutedtimestep + 1;
                std::string dubuf = "[OK] processing resuming data. Last executed time step " + std::to_string(lastexecutedtimestep);
                ExternalModules::print11(control.layoutnumber, dubuf);
            }

            if (initialtimestep > control.finaltimestep) {
                ExternalModules::stoponerror(control.layoutnumber, control.num_procs, "Initial time step greater than final one", true);
                destroy_and_deallocate();
                return;
            }

            ExternalModules::crea_timevector(sgg, lastexecutedtimestep, control.finaltimestep, lastexecutedtime);

            bool attinformado = false;
            ExternalModules::updateSigmaM(attinformado);
            ExternalModules::updateThinWiresSigma(attinformado);
            ExternalModules::calc_G1G2Gm1Gm2(sgg, g, eps0, mu0);
            ExternalModules::revertThinWiresSigma();

#ifdef CompileWithMPI
            int ierr;
            ExternalModules::MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
            std::string dubuf = "Init Reporting...";
            ExternalModules::print11(control.layoutnumber, dubuf);
            ExternalModules::InitReporting(sgg, &control);
            ExternalModules::reportSimulationOptions();

#ifdef CompileWithMPI
            ExternalModules::MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
            dubuf = "[OK]";
            ExternalModules::print11(control.layoutnumber, dubuf);

#ifdef CompileWithMPI
            ExternalModules::MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

            ExternalModules::initializeBorders();
            ExternalModules::initializeLumped();
            ExternalModules::initializeWires();
            ExternalModules::initializeAnisotropic();
            ExternalModules::initializeSGBC();
            ExternalModules::initializeMultiports();
            ExternalModules::initializeEDispersives();
            ExternalModules::initializeMDispersives();
            ExternalModules::initializePlanewave();
            ExternalModules::initializeNodalSources();

            // fillMtag call would go here with proper matrix references
            // ExternalModules::fillMtag(...);

            ExternalModules::initializeObservation();

#ifdef CompileWithMPI
            ExternalModules::initializeMPI();
#endif

#ifdef CompileWithMPI
            ExternalModules::MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

            if (control.resume) {
                // Close file handle if opened
            }

            n = initialtimestep;
            ini_save = initialtimestep;
            n_info = 5 + initialtimestep;

            dubuf = "Init Timing...";
            ExternalModules::print11(control.layoutnumber, dubuf);
            ExternalModules::InitTiming(sgg, &control, control.time_desdelanzamiento, initialtimestep, control.maxSourceValue);

            ExternalModules::CLOSEWARNINGFILE(control.layoutnumber, control.num_procs, control.fatalerror, false, control.simu_devia);

            if (control.fatalerror) {
                dubuf = "FATAL ERRORS. Revise *Warnings.txt file. ABORTING...";
                ExternalModules::print11(control.layoutnumber, dubuf);
                // Abort logic would go here
            }
        }

        void run() {
            // Placeholder for run logic
        }

        void end() {
            // Placeholder for end logic
        }

        void destroy_and_deallocate() {
            // Deallocate vectors
            Ex.clear();
            Ey.clear();
            Ez.clear();
            Hx.clear();
            Hy.clear();
            Hz.clear();
            Idxe.clear();
            Idye.clear();
            Idze.clear();
            Idxh.clear();
            Idyh.clear();
            Idzh.clear();
            dxe.clear();
            dye.clear();
            dze.clear();
            dxh.clear();
            dyh.clear();
            dzh.clear();
            g.g1.clear();
            g.g2.clear();
            g.gm1.clear();
            g.gm2.clear();
            SINPML_fullsize.clear();
            fullsize.clear();
        }

        // Forward declarations for other methods that are not implemented in this snippet
        void advanceE();
        void advanceEx();
        void advanceEy();
        void advanceEz();
        void advanceH();
        void advanceHx();
        void advanceHy();
        void advanceHz();
        void advancePlaneWaveE();
        void advancePlaneWaveH();
        void advanceWiresE();
        void advanceWiresH();
        void advancePMLE();
        void advanceAnisotropicE();
        void advanceAnisotropicH();
        void advanceLumpedE();
        void advanceNodalE();
        void advanceNodalH();
        void advancePMLbodyH();
        void advanceMagneticCPML();
        void advanceSGBCE();
        void advanceSGBCH();
        void advanceEDispersiveE();
        void advanceMDispersiveH();
        void MinusCloneMagneticPMC();
        void CloneMagneticPeriodic();
        void advanceMagneticMUR();

#ifdef CompileWithMPI
        void init_MPIConformalProbes();
#endif
    };

} // namespace Solver_m

stoponerror(this->control.layoutnumber, this->control.num_procs, dubuf, true); // para que retorne
            this->destroy_and_deallocate();
            return;
        }
#ifdef CompileWithMPI
        flushMPIdata();
#endif

        //!!!no se si el orden wires - sgbcs del sync importa 150519
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
        if (this->control.stochastic) {
            syncstoch_mpi_sgbcs(this->control.simu_devia, this->control.layoutnumber, this->control.num_procs);
            syncstoch_mpi_lumped(this->control.simu_devia, this->control.layoutnumber, this->control.num_procs);
        }
#endif
#endif

        printSimulationStart();

    } // end of public methods block, entering contains/private methods

    void findbounds(bounds_t& b) {
        //
        //No tocar. Dejar como estan alocateados
        b.dxe.XI = this->sgg.alloc(iHx).XI;
        b.dxe.XE = this->sgg.alloc(iHx).XE;
        b.dye.YI = this->sgg.alloc(iHy).YI;
        b.dye.YE = this->sgg.alloc(iHy).YE;
        b.dze.ZI = this->sgg.alloc(iHz).ZI;
        b.dze.ZE = this->sgg.alloc(iHz).ZE;
        //
        b.dxh.XI = this->sgg.alloc(iEx).XI;
        b.dxh.XE = this->sgg.alloc(iEx).XE;
        b.dyh.YI = this->sgg.alloc(iEy).YI;
        b.dyh.YE = this->sgg.alloc(iEy).YE;
        b.dzh.ZI = this->sgg.alloc(iEz).ZI;
        b.dzh.ZE = this->sgg.alloc(iEz).ZE;

        //
        //No tocar. Dejar como estan alocateados
        b.Ex.XI = this->sgg.Alloc(iEx).XI;
        b.Ex.XE = this->sgg.Alloc(iEx).XE;
        b.Ey.XI = this->sgg.Alloc(iEy).XI;
        b.Ey.XE = this->sgg.Alloc(iEy).XE;
        b.Ez.XI = this->sgg.Alloc(iEz).XI;
        b.Ez.XE = this->sgg.Alloc(iEz).XE;
        //
        b.Hx.XI = this->sgg.Alloc(iHx).XI;
        b.Hx.XE = this->sgg.Alloc(iHx).XE;
        b.Hy.XI = this->sgg.Alloc(iHy).XI;
        b.Hy.XE = this->sgg.Alloc(iHy).XE;
        b.Hz.XI = this->sgg.Alloc(iHz).XI;
        b.Hz.XE = this->sgg.Alloc(iHz).XE;
        //
        b.Ex.YI = this->sgg.Alloc(iEx).YI;
        b.Ex.YE = this->sgg.Alloc(iEx).YE;
        b.Ey.YI = this->sgg.Alloc(iEy).YI;
        b.Ey.YE = this->sgg.Alloc(iEy).YE;
        b.Ez.YI = this->sgg.Alloc(iEz).YI;
        b.Ez.YE = this->sgg.Alloc(iEz).YE;
        //
        b.Hx.YI = this->sgg.Alloc(iHx).YI;
        b.Hx.YE = this->sgg.Alloc(iHx).YE;
        b.Hy.YI = this->sgg.Alloc(iHy).YI;
        b.Hy.YE = this->sgg.Alloc(iHy).YE;
        b.Hz.YI = this->sgg.Alloc(iHz).YI;
        b.Hz.YE = this->sgg.Alloc(iHz).YE;
        //
        b.Ex.ZI = this->sgg.Alloc(iEx).ZI;
        b.Ex.ZE = this->sgg.Alloc(iEx).ZE;
        b.Ey.ZI = this->sgg.Alloc(iEy).ZI;
        b.Ey.ZE = this->sgg.Alloc(iEy).ZE;
        b.Ez.ZI = this->sgg.Alloc(iEz).ZI;
        b.Ez.ZE = this->sgg.Alloc(iEz).ZE;
        //
        b.Hx.ZI = this->sgg.Alloc(iHx).ZI;
        b.Hx.ZE = this->sgg.Alloc(iHx).ZE;
        b.Hy.ZI = this->sgg.Alloc(iHy).ZI;
        b.Hy.ZE = this->sgg.Alloc(iHy).ZE;
        b.Hz.ZI = this->sgg.Alloc(iHz).ZI;
        b.Hz.ZE = this->sgg.Alloc(iHz).ZE;
        //
        //
        //

        //matrix indexes. Nothing to change. Asi estan alocateados
        b.sggMiEx.XI = this->sgg.Alloc(iEx).XI;
        b.sggMiEx.XE = this->sgg.Alloc(iEx).XE;
        b.sggMiEy.XI = this->sgg.Alloc(iEy).XI;
        b.sggMiEy.XE = this->sgg.Alloc(iEy).XE;
        b.sggMiEz.XI = this->sgg.Alloc(iEz).XI;
        b.sggMiEz.XE = this->sgg.Alloc(iEz).XE;
        //
        b.sggMiHx.XI = this->sgg.Alloc(iHx).XI;
        b.sggMiHx.XE = this->sgg.Alloc(iHx).XE;
        b.sggMiHy.XI = this->sgg.Alloc(iHy).XI;
        b.sggMiHy.XE = this->sgg.Alloc(iHy).XE;
        b.sggMiHz.XI = this->sgg.Alloc(iHz).XI;
        b.sggMiHz.XE = this->sgg.Alloc(iHz).XE;
        //
        b.sggMiEx.YI = this->sgg.Alloc(iEx).YI;
        b.sggMiEx.YE = this->sgg.Alloc(iEx).YE;
        b.sggMiEy.YI = this->sgg.Alloc(iEy).YI;
        b.sggMiEy.YE = this->sgg.Alloc(iEy).YE;
        b.sggMiEz.YI = this->sgg.Alloc(iEz).YI;
        b.sggMiEz.YE = this->sgg.Alloc(iEz).YE;
        //
        b.sggMiHx.YI = this->sgg.Alloc(iHx).YI;
        b.sggMiHx.YE = this->sgg.Alloc(iHx).YE;
        b.sggMiHy.YI = this->sgg.Alloc(iHy).YI;
        b.sggMiHy.YE = this->sgg.Alloc(iHy).YE;
        b.sggMiHz.YI = this->sgg.Alloc(iHz).YI;
        b.sggMiHz.YE = this->sgg.Alloc(iHz).YE;
        //
        b.sggMiEx.ZI = this->sgg.Alloc(iEx).ZI;
        b.sggMiEx.ZE = this->sgg.Alloc(iEx).ZE;
        b.sggMiEy.ZI = this->sgg.Alloc(iEy).ZI;
        b.sggMiEy.ZE = this->sgg.Alloc(iEy).ZE;
        b.sggMiEz.ZI = this->sgg.Alloc(iEz).ZI;
        b.sggMiEz.ZE = this->sgg.Alloc(iEz).ZE;
        //
        b.sggMiHx.ZI = this->sgg.Alloc(iHx).ZI;
        b.sggMiHx.ZE = this->sgg.Alloc(iHx).ZE;
        b.sggMiHy.ZI = this->sgg.Alloc(iHy).ZI;
        b.sggMiHy.ZE = this->sgg.Alloc(iHy).ZE;
        b.sggMiHz.ZI = this->sgg.Alloc(iHz).ZI;
        b.sggMiHz.ZE = this->sgg.Alloc(iHz).ZE;
        //
        //
        //
        b.sweepEx.XI = this->sgg.Sweep(iEx).XI;
        b.sweepEx.XE = this->sgg.Sweep(iEx).XE;
        b.sweepEy.XI = this->sgg.Sweep(iEy).XI;
        b.sweepEy.XE = this->sgg.Sweep(iEy).XE;
        b.sweepEz.XI = this->sgg.Sweep(iEz).XI;
        b.sweepEz.XE = this->sgg.Sweep(iEz).XE;
        //
        b.sweepHx.XI = this->sgg.Sweep(iHx).XI;
        b.sweepHx.XE = this->sgg.Sweep(iHx).XE;
        b.sweepHy.XI = this->sgg.Sweep(iHy).XI;
        b.sweepHy.XE = this->sgg.Sweep(iHy).XE;
        b.sweepHz.XI = this->sgg.Sweep(iHz).XI;
        b.sweepHz.XE = this->sgg.Sweep(iHz).XE;
        //
        //
        b.sweepEx.YI = this->sgg.Sweep(iEx).YI;
        b.sweepEx.YE = this->sgg.Sweep(iEx).YE;
        b.sweepEy.YI = this->sgg.Sweep(iEy).YI;
        b.sweepEy.YE = this->sgg.Sweep(iEy).YE;
        b.sweepEz.YI = this->sgg.Sweep(iEz).YI;
        b.sweepEz.YE = this->sgg.Sweep(iEz).YE;
        //
        b.sweepHx.YI = this->sgg.Sweep(iHx).YI;
        b.sweepHx.YE = this->sgg.Sweep(iHx).YE;
        b.sweepHy.YI = this->sgg.Sweep(iHy).YI;
        b.sweepHy.YE = this->sgg.Sweep(iHy).YE;
        b.sweepHz.YI = this->sgg.Sweep(iHz).YI;
        b.sweepHz.YE = this->sgg.Sweep(iHz).YE;
        //
        b.sweepEx.ZI = this->sgg.Sweep(iEx).ZI;
        b.sweepEx.ZE = this->sgg.Sweep(iEx).ZE;
        b.sweepEy.ZI = this->sgg.Sweep(iEy).ZI;
        b.sweepEy.ZE = this->sgg.Sweep(iEy).ZE;
        b.sweepEz.ZI = this->sgg.Sweep(iEz).ZI;
        b.sweepEz.ZE = this->sgg.Sweep(iEz).ZE;
        //
        b.sweepHx.ZI = this->sgg.Sweep(iHx).ZI;
        b.sweepHx.ZE = this->sgg.Sweep(iHx).ZE;
        b.sweepHy.ZI = this->sgg.Sweep(iHy).ZI;
        b.sweepHy.ZE = this->sgg.Sweep(iHy).ZE;
        b.sweepHz.ZI = this->sgg.Sweep(iHz).ZI;
        b.sweepHz.ZE = this->sgg.Sweep(iHz).ZE;
        //
        b.sweepSINPMLEx.XI = this->sgg.SINPMLSweep(iEx).XI;
        b.sweepSINPMLEy.XI = this->sgg.SINPMLSweep(iEy).XI;
        b.sweepSINPMLEz.XI = this->sgg.SINPMLSweep(iEz).XI;
        b.sweepSINPMLHx.XI = this->sgg.SINPMLSweep(iHx).XI;
        b.sweepSINPMLHy.XI = this->sgg.SINPMLSweep(iHy).XI;
        b.sweepSINPMLHz.XI = this->sgg.SINPMLSweep(iHz).XI;
        //
        b.sweepSINPMLEx.XE = this->sgg.SINPMLSweep(iEx).XE;
        b.sweepSINPMLEy.XE = this->sgg.SINPMLSweep(iEy).XE;
        b.sweepSINPMLEz.XE = this->sgg.SINPMLSweep(iEz).XE;
        b.sweepSINPMLHx.XE = this->sgg.SINPMLSweep(iHx).XE;
        b.sweepSINPMLHy.XE = this->sgg.SINPMLSweep(iHy).XE;
        b.sweepSINPMLHz.XE = this->sgg.SINPMLSweep(iHz).XE;
        //
        b.sweepSINPMLEx.YI = this->sgg.SINPMLSweep(iEx).YI;
        b.sweepSINPMLEy.YI = this->sgg.SINPMLSweep(iEy).YI;
        b.sweepSINPMLEz.YI = this->sgg.SINPMLSweep(iEz).YI;
        b.sweepSINPMLHx.YI = this->sgg.SINPMLSweep(iHx).YI;
        b.sweepSINPMLHy.YI = this->sgg.SINPMLSweep(iHy).YI;
        b.sweepSINPMLHz.YI = this->sgg.SINPMLSweep(iHz).YI;
        //
        b.sweepSINPMLEx.YE = this->sgg.SINPMLSweep(iEx).YE;
        b.sweepSINPMLEy.YE = this->sgg.SINPMLSweep(iEy).YE;
        b.sweepSINPMLEz.YE = this->sgg.SINPMLSweep(iEz).YE;
        b.sweepSINPMLHx.YE = this->sgg.SINPMLSweep(iHx).YE;
        b.sweepSINPMLHy.YE = this->sgg.SINPMLSweep(iHy).YE;
        b.sweepSINPMLHz.YE = this->sgg.SINPMLSweep(iHz).YE;
        //
        b.sweepSINPMLEx.ZI = this->sgg.SINPMLSweep(iEx).ZI;
        b.sweepSINPMLEy.ZI = this->sgg.SINPMLSweep(iEy).ZI;
        b.sweepSINPMLEz.ZI = this->sgg.SINPMLSweep(iEz).ZI;
        b.sweepSINPMLHx.ZI = this->sgg.SINPMLSweep(iHx).ZI;
        b.sweepSINPMLHy.ZI = this->sgg.SINPMLSweep(iHy).ZI;
        b.sweepSINPMLHz.ZI = this->sgg.SINPMLSweep(iHz).ZI;
        //
        b.sweepSINPMLEx.ZE = this->sgg.SINPMLSweep(iEx).ZE;
        b.sweepSINPMLEy.ZE = this->sgg.SINPMLSweep(iEy).ZE;
        b.sweepSINPMLEz.ZE = this->sgg.SINPMLSweep(iEz).ZE;
        b.sweepSINPMLHx.ZE = this->sgg.SINPMLSweep(iHx).ZE;
        b.sweepSINPMLHy.ZE = this->sgg.SINPMLSweep(iHy).ZE;
        b.sweepSINPMLHz.ZE = this->sgg.SINPMLSweep(iHz).ZE;

        //

        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //find lenghts
        //this is automatic. Nothing to change
        //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //
        b.Ex.NX = b.Ex.XE - b.Ex.XI + 1;
        b.Ex.NY = b.Ex.YE - b.Ex.YI + 1;
        b.Ex.NZ = b.Ex.ZE - b.Ex.ZI + 1;

        b.Ey.NX = b.Ey.XE - b.Ey.XI + 1;
        b.Ey.NY = b.Ey.YE - b.Ey.YI + 1;
        b.Ey.NZ = b.Ey.ZE - b.Ey.ZI + 1;

        b.Ez.NX = b.Ez.XE - b.Ez.XI + 1;
        b.Ez.NY = b.Ez.YE - b.Ez.YI + 1;
        b.Ez.NZ = b.Ez.ZE - b.Ez.ZI + 1;
        //
        b.Hx.NX = b.Hx.XE - b.Hx.XI + 1;
        b.Hx.NY = b.Hx.YE - b.Hx.YI + 1;
        b.Hx.NZ = b.Hx.ZE - b.Hx.ZI + 1;
        //
        b.Hy.NX = b.Hy.XE - b.Hy.XI + 1;
        b.Hy.NY = b.Hy.YE - b.Hy.YI + 1;
        b.Hy.NZ = b.Hy.ZE - b.Hy.ZI + 1;
        //
        b.Hz.NX = b.Hz.XE - b.Hz.XI + 1;
        b.Hz.NY = b.Hz.YE - b.Hz.YI + 1;
        b.Hz.NZ = b.Hz.ZE - b.Hz.ZI + 1;
        //
        //
        b.sweepEx.NX = b.sweepEx.XE - b.sweepEx.XI + 1;
        b.sweepEx.NY = b.sweepEx.YE - b.sweepEx.YI + 1;
        b.sweepEx.NZ = b.sweepEx.ZE - b.sweepEx.ZI + 1;
        //
        b.sweepEy.NX = b.sweepEy.XE - b.sweepEy.XI + 1;
        b.sweepEy.NY = b.sweepEy.YE - b.sweepEy.YI + 1;
        b.sweepEy.NZ = b.sweepEy.ZE - b.sweepEy.ZI + 1;
        //
        b.sweepEz.NX = b.sweepEz.XE - b.sweepEz.XI + 1;
        b.sweepEz.NY = b.sweepEz.YE - b.sweepEz.YI + 1;
        b.sweepEz.NZ = b.sweepEz.ZE - b.sweepEz.ZI + 1;
        //
        b.sweepHx.NX = b.sweepHx.XE - b.sweepHx.XI + 1;
        b.sweepHx.NY = b.sweepHx.YE - b.sweepHx.YI + 1;
        b.sweepHx.NZ = b.sweepHx.ZE - b.sweepHx.ZI + 1;
        //
        b.sweepHy.NX = b.sweepHy.XE - b.sweepHy.XI + 1;
        b.sweepHy.NY = b.sweepHy.YE - b.sweepHy.YI + 1;
        b.sweepHy.NZ = b.sweepHy.ZE - b.sweepHy.ZI + 1;
        //
        b.sweepHz.NX = b.sweepHz.XE - b.sweepHz.XI + 1;
        b.sweepHz.NY = b.sweepHz.YE - b.sweepHz.YI + 1;
        b.sweepHz.NZ = b.sweepHz.ZE - b.sweepHz.ZI + 1;
        //
        //
        b.sggMiEx.NX = b.sggMiEx.XE - b.sggMiEx.XI + 1;
        b.sggMiEx.NY = b.sggMiEx.YE - b.sggMiEx.YI + 1;
        b.sggMiEx.NZ = b.sggMiEx.ZE - b.sggMiEx.ZI + 1;
        b.sggMiEy.NX = b.sggMiEy.XE - b.sggMiEy.XI + 1;
        b.sggMiEy.NY = b.sggMiEy.YE - b.sggMiEy.YI + 1;
        b.sggMiEy.NZ = b.sggMiEy.ZE - b.sggMiEy.ZI + 1;
        b.sggMiEz.NX = b.sggMiEz.XE - b.sggMiEz.XI + 1;
        b.sggMiEz.NY = b.sggMiEz.YE - b.sggMiEz.YI + 1;
        b.sggMiEz.NZ = b.sggMiEz.ZE - b.sggMiEz.ZI + 1;
        //
        b.sggMiHx.NX = b.sggMiHx.XE - b.sggMiHx.XI + 1;
        b.sggMiHx.NY = b.sggMiHx.YE - b.sggMiHx.YI + 1;
        b.sggMiHx.NZ = b.sggMiHx.ZE - b.sggMiHx.ZI + 1;
        b.sggMiHy.NX = b.sggMiHy.XE - b.sggMiHy.XI + 1;
        b.sggMiHy.NY = b.sggMiHy.YE - b.sggMiHy.YI + 1;
        b.sggMiHy.NZ = b.sggMiHy.ZE - b.sggMiHy.ZI + 1;
        b.sggMiHz.NX = b.sggMiHz.XE - b.sggMiHz.XI + 1;
        b.sggMiHz.NY = b.sggMiHz.YE - b.sggMiHz.YI + 1;
        b.sggMiHz.NZ = b.sggMiHz.ZE - b.sggMiHz.ZI + 1;
        //
        //
        //estas longitudes son relativas al layout !ojo
        b.dxe.NX = b.dxe.XE - b.dxe.XI + 1;
        b.dye.NY = b.dye.YE - b.dye.YI + 1;
        b.dze.NZ = b.dze.ZE - b.dze.ZI + 1;
        //
        b.dxh.NX = b.dxh.XE - b.dxh.XI + 1;
        b.dyh.NY = b.dyh.YE - b.dyh.YI + 1;
        b.dzh.NZ = b.dzh.ZE - b.dzh.ZI + 1;

    }

    void updateSigmaM(bool& att) {
        double deltaespmax, fmax, skin_depth;
        bool hayattmedia = false;
        double mur, epr;
        char buff[BUFSIZE];
        int i;
        if (std::abs(this->control.attfactorc - 1.0_RKIND) > 1.0e-12_RKIND) {
            att = false;
            for (i = 1; i <= this->sgg.nummedia; ++i) {
                if (this->sgg.Med(i).Is.MultiportPadding) {
                    this->sgg.Med(i).SigmaM = (-2.0_RKIND * (-1.0_RKIND + this->control.attfactorc) * this->mu0) / ((1 + this->control.attfactorc) * this->sgg.dt);
                    hayattmedia = true;
                }
                deltaespmax = std::max(std::max(std::maxval(this->sgg.dx), std::maxval(this->sgg.dy)), std::maxval(this->sgg.dz));
                if (hayattmedia && !att) {
                    !!!!info on stabilization
                    epr = 1.0_RKIND;
                    mur = 1.0_RKIND;
                    !!
                    sprintf(buff, " Composites stabilization att. factor=%e%e", this->control.attfactorc, this->sgg.Med(i).SigmaM);

                    WarnErrReport(buff);
                    !!
                    fmax = 1.0_RKIND / (10.0_RKIND * this->sgg.dt);
                    skin_depth = 1.0_RKIND / (std::sqrt(2.0_RKIND) * fmax * Pi * std::pow(epr * this->eps0 * std::pow(4 * mur * std::pow(this->mu0, 2.0_RKIND) + std::pow(this->sgg.Med(i).Sigmam, 2) / (std::pow(fmax, 2) * std::pow(Pi, 2.0_RKIND))), 0.25_RKIND) * std::sin(std::atan2(2 * Pi * epr * this->eps0 * mur * this->mu0, -(epr * this->eps0 * this->sgg.Med(i).Sigmam) / fmax) / 2.0_RKIND));
                    sprintf(buff, " At 10 samp/per f=%e,Max Att(dB)=%e", fmax,
                            -(0.0001295712360834271997 * std::imag(fmax * std::sqrt((epr * ((0, -2.825225e7) +
                            8.8757061047382236e6 * mur + this->control.attfactorc * ((0, 2.825225e7) + 8.8757061047382236e6 * mur))) /
                            (1.124121310242e12 + 1.124121310242e12 * this->control.attfactorc)) * std::min(deltaespmax, skin_depth))));
                    if (this->control.layoutnumber == 0) WarnErrReport(buff);
                    if (fmax > 3e9) {
                        fmax = 3e9;
                        sprintf(buff, "             At f=%e,Max Att(dB)=%e", fmax,
                                -(0.0001295712360834271997 * std::imag(fmax * std::sqrt((epr * ((0, -2.825225e7) +
                                8.8757061047382236e6 * mur + this->control.attfactorc * ((0, 2.825225e7) + 8.8757061047382236e6 * mur))) /
                                (1.124121310242e12 + 1.124121310242e12 * this->control.attfactorc)) * std::min(deltaespmax, skin_depth))));
                        if (this->control.layoutnumber == 0) WarnErrReport(buff);
                    }
                    att = true;
                }
            }
        }
    }

    void updateThinWiresSigma(bool& att) {
        char buff[BUFSIZE];
        int i;
        if (std::abs(this->control.attfactorw - 1.0_RKIND) > 1.0e-12_RKIND) {
            att = false;
            for (i = 1; i <= this->sgg.nummedia; ++i) {
                if (this->sgg.Med(i).Is.ThinWire) {
                    this->sgg.Med(i).Sigma = (-2.0_RKIND * (-1.0_RKIND + this->control.attfactorw) * this->eps0) / ((1 + this->control.attfactorw) * this->sgg.dt);
                    if (!att) {
                        sprintf(buff, " WIREs stabilization att. factors=%e%e", this->control.attfactorw, this->sgg.Med(i).Sigma);
                        if (this->control.layoutnumber == 0) WarnErrReport(buff);
                        att = true;
                    }
                }
            }
        }
    }

    void revertThinWiresSigma() {
        int i;
        if (std::abs(this->control.attfactorw - 1.0_RKIND) > 1.0e-12_RKIND) {
            for (i = 1; i <= this->sgg.nummedia; ++i) {
                if (this->sgg.Med(i).Is.ThinWire) {
                    this->sgg.Med(i).Sigma = 0.0_RKIND; //revert!!! !necesario para no lo tome como un lossy luego en wires !solo se toca el g1,g2
                }
            }
        }
    }

    void reportSimulationOptions() {
        char buff[BUFSIZE];
        if ((this->control.layoutnumber == 0) && this->control.verbose) {
            sprintf(buff, "CPML  alpha, alphaorder, kappa factors= %e%e%e", this->control.alphamaxpar, this->control.alphaOrden, this->control.kappamaxpar);
            WarnErrReport(buff);
            if (this->control.medioextra.exists) {
                sprintf(buff, "CPML correction size,factor to scale sigmamax = %i%e", this->control.medioextra.pml_size, this->control.medioextra.sigma);
                WarnErrReport(buff);
            }
            sprintf(buff, "saveall=%i, flushsecondsFields=%i, flushsecondsData=%i, maxCPUtime=%i, singlefilewrite=%i", this->control.saveall, this->control.flushsecondsFields, this->control.flushsecondsData, this->control.maxCPUtime, this->control.singlefilewrite);
            WarnErrReport(buff);
            sprintf(buff, "TAPARRABOS=%i, wiresflavor=%s, mindistwires=%i, wirecrank=%i makeholes=%i", this->control.TAPARRABOS, trim(adjustl(this->control.wiresflavor)), this->control.mindistwires, this->control.wirecrank, this->control.makeholes);
            WarnErrReport(buff);
            sprintf(buff, "connectendings=%i, isolategroupgroups=%i", this->control.connectendings, this->control.isolategroupgroups);
            WarnErrReport(buff);
            sprintf(buff, "wirethickness %i stableradholland=%i mtlnberenger=%i inductance_model=%s, inductance_order=%i, groundwires=%i ,fieldtotl=%i noSlantedcrecepelo =%i", this->control.wirethickness, this->control.stableradholland, this->control.mtlnberenger, trim(adjustl(this->control.inductance_model)), this->control.inductance_order, this->control.groundwires, this->control.fieldtotl, this->control.noSlantedcrecepelo);
            WarnErrReport(buff);
            sprintf(buff, "sgbc=%i, mibc=%i, attfactorc=%i, attfactorw=%i", this->control.sgbc, this->control.mibc, this->control.attfactorc, this->control.attfactorw);
            WarnErrReport(buff);
            sprintf(buff, "NOcompomur=%i, ADE=%i, conformalskin=%i, sgbcFreq=%i, sgbcresol=%i, sgbccrank=%i, sgbcDepth=%i", this->control.NOcompomur, this->control.ADE, this->control.conformalskin, this->control.sgbcFreq, this->control.sgbcresol, this->control.sgbccrank, this->control.sgbcdepth);
            WarnErrReport(buff);
            sprintf(buff, "mur_second=%i, murafterpml=%i, facesNF2FF%%tr=%i, facesNF2FF%%fr=%i, facesNF2FF%%iz=%i", this->control.mur_second, this->control.murafterpml, this->control.facesNF2FF.tr, this->control.facesNF2FF.fr, this->control.facesNF2FF.iz);
            WarnErrReport(buff);
            sprintf(buff, "facesNF2FF%%de=%i, facesNF2FF%%ab=%i, facesNF2FF%%ar=%i, NF2FFDecim=%i", this->control.facesNF2FF.de, this->control.facesNF2FF.ab, this->control.facesNF2FF.ar, this->control.NF2FFDecim);
            WarnErrReport(buff);
        }
    }

    void initializeBorders() {
        char dubuf[BUFSIZE];
        bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
        int ierr;
#endif
        sprintf(dubuf, "Init Other Borders...");
        print11(this->control.layoutnumber, dubuf);
        InitOtherBorders(this->sgg, this->thereAre);
        l_auxinput = this->thereAre.PECBorders || this->thereAre.PMCBorders || this->thereAre.PeriodicBorders;
        l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
        if (l_auxoutput) {
            sprintf(dubuf, "----> there are PEC, PMC or periodic Borders");
            print11(this->control.layoutnumber, dubuf);
        } else {
            sprintf(dubuf, "----> no PEC, PMC or periodic Borders found");
            print11(this->control.layoutnumber, dubuf);
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
        sprintf(dubuf, "Init CPML Borders...");
        print11(this->control.layoutnumber, dubuf);
        InitCPMLBorders(this->sgg, this->sinPML_fullsize, this->thereAre.PMLBorders, this->control,
                        dxe, dye, dze, dxh, dyh, dzh, Idxe, Idye, Idze, Idxh, Idyh, Idzh, this->eps0, this->mu0);

        l_auxinput = this->thereAre.PMLBorders;
        l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
        if (l_auxoutput) {
            sprintf(dubuf, "----> there are CPML Borders");
            print11(this->control.layoutnumber, dubuf);
        } else {
            sprintf(dubuf, "----> no CPML Borders found");
            print11(this->control.layoutnumber, dubuf);
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
        sprintf(dubuf, "Init PML Bodies...");
        print11(this->control.layoutnumber, dubuf);
        InitPMLbodies(this->sgg, this->media, Ex, Ey, Ez, Hx, Hy, Hz, IDxe, IDye, IDze, IDxh, IDyh, IDzh, this->g.g2, this->g.gm2, this->thereAre.PMLbodies, this->control, this->eps0, this->mu0);
        l_auxinput = this->thereAre.PMLbodies;
        l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
        if (l_auxoutput) {
            sprintf(dubuf, "----> there are PML Bodies");
            print11(this->control.layoutnumber, dubuf);
        } else {
            sprintf(dubuf, "----> no PML Bodies found");
            print11(this->control.layoutnumber, dubuf);
        }
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
        sprintf(dubuf, "Init Mur Borders...");
        print11(this->control.layoutnumber, dubuf);
        InitMURBorders(this->sgg, this->thereAre.MURBorders, this->control.resume, Idxh, Idyh, Idzh, this->eps0, this->mu0);
        l_auxinput = this->thereAre.MURBorders;
        l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
        if (l_auxoutput) {
            sprintf(dubuf, "----> there are Mur Borders");
            print11(this->control.layoutnumber, dubuf);
        } else {
            sprintf(dubuf, "----> no Mur Borders found");
            print11(this->control.layoutnumber, dubuf);
        }

    }

    void initializeLumped() {
        char dubuf[BUFSIZE];
        bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
        int ierr;
#endif

        //init lumped debe ir antes de wires porque toca la conductividad del material !mmmm ojoooo 120123
        sprintf(dubuf, "Init Lumped Elements...");
        print11(this->control.layoutnumber, dubuf);
        InitLumped(this->sgg, this->media, Ex, Ey, Ez, Hx, Hy, Hz, IDxe, IDye, IDze, IDxh, IDyh, IDzh, this->control, this->thereAre.Lumpeds, this->eps0, this->mu0);
        l_auxinput = this->thereAre.Lumpeds;
        l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
        if (l_auxoutput) {
            sprintf(dubuf, "----> there are Structured lumped elements");
            print11(this->control.layoutnumber, dubuf);
        } else {
            sprintf(dubuf, "----> no lumped Structured elements found");
            print11(this->control.layoutnumber, dubuf);
        }
    }

    void initializeWires() {
        double dtcritico, newdtcritico;
        char dubuf[BUFSIZE], buff[BUFSIZE];
        bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
        int ierr;
#endif

        dtcritico = this->sgg.dt;
#ifndef CompileWithMTLN
        if ((trim(adjustl(this->control.wiresflavor)) == "holland") ||
            (trim(adjustl(this->control.wiresflavor)) == "transition")) {
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
            sprintf(dubuf, "Init Holland Wires...");
            print11(this->control.layoutnumber, dubuf);
            InitWires(this->sgg, this->media.sggMiNo, this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz, this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz,
                      this->thereAre.Wires, Ex, Ey, Ez, Hx, Hy, Hz, Idxe, Idye, Idze, Idxh, Idyh, Idzh,
                      this->g.g2, this->sinPML_fullsize, this->fullsize, dtcritico, this->eps0, this->mu0, this->control);
            l_auxinput = this->thereAre.Wires;
            l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, &ierr);
            MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
            if (l_auxoutput) {
                sprintf(dubuf, "----> there are Holland/transition wires");
                print11(this->control.layoutnumber, dubuf);
            } else {
                sprintf(dubuf, "----> no Holland/transition wires found");
                print11(this->control.layoutnumber, dubuf);
            }
        }

#ifdef CompileWithBerengerWires
        if (trim(adjustl(this->control.wiresflavor)) == "berenger") {

#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
            sprintf(dubuf, "Init Multi-Wires...");
            print11(this->control.layoutnumber, dubuf);
            InitWires_Berenger(
#endif

this->sgg, this->media.sggMiNo, this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz, this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz, this->control.layoutnumber, this->control.num_procs, this->thereAre.Wires, this->control.resume, this->control.makeholes,
            this->control.isolategroupgroups, this->control.mtlnberenger, this->control.mindistwires,
            this->control.groundwires, this->control.taparrabos, Ex, Ey, Ez,
            Idxe, Idye, Idze, Idxh, Idyh, Idzh, this->control.inductance_model, this->g.g2, this->sinPML_fullsize, this->fullsize, dtcritico, this->eps0, this->mu0, this->control.verbose);
         l_auxinput = this->thereAre.Wires;
         l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif

         if (l_auxoutput) {
            sprintf(dubuf, "----> there are Multi-wires");
            print11(this->control.layoutnumber, dubuf);
         } else {
            sprintf(dubuf, "----> no Multi-wires found");
            print11(this->control.layoutnumber, dubuf);
         }
      }
#endif
#ifdef CompileWithSlantedWires
      if ((std::string(this->control.wiresflavor).find("slanted") != std::string::npos) || (std::string(this->control.wiresflavor).find("semistructured") != std::string::npos)) {

#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
         sprintf(dubuf, "Init Slanted Wires...");
         print11(this->control.layoutnumber, dubuf);
         if (std::string(this->control.wiresflavor).find("semistructured") != std::string::npos) {
            sprintf(dubuf, "...%d", this->control.precision);
            print11(this->control.layoutnumber, dubuf);
            estructura_slanted(this->sgg, this->control.precision);
         } else {
            // continue
         }
         InitWires_Slanted(this->sgg, this->control.layoutnumber, this->control.num_procs, Ex, Ey, Ez,
                           Idxe, Idye, Idze, Idxh, Idyh, Idzh,
                           this->media.sggMiNo,
                           this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz,
                           this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz,
                           this->thereAre.Wires, this->control.resume,
                           this->control.mindistwires, this->control.groundwires, this->control.noSlantedcrecepelo,
                           this->control.inductance_model, this->control.inductance_order,
                           this->g.g2, this->sinPML_fullsize, dtcritico, this->eps0, this->mu0, this->control.verbose);
         l_auxinput = this->thereAre.Wires;
         l_auxoutput = l_auxinput;
         // check for MUR1 nodes sgg 230124
         init_murABC_slanted(this->sgg, this->sinPML_fullsize, this->eps0, this->mu0);
         // !!!!!

#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif

         if (l_auxoutput) {
            sprintf(dubuf, "----> there are Slanted wires");
            print11(this->control.layoutnumber, dubuf);
         } else {
            sprintf(dubuf, "----> no Slanted wires found");
            print11(this->control.layoutnumber, dubuf);
         }
      }
#endif


#else // else of #ifndef CompileWithMTLN
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init MTLN Wires...");
      print11(this->control.layoutnumber, dubuf);
      InitWires_mtln(this->sgg, Ex, Ey, Ez,
                     this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz,
                     this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz,
                     this->eps0, this->mu0, this->mtln_parsed, this->thereAre.MTLNbundles, dtcritico);
#endif


      // !!!sincroniza el dtcritico
#ifdef CompileWithMPI
      newdtcritico = 0.0_RKIND_tiempo;
      MPI_AllReduce(&dtcritico, &newdtcritico, 1, REALSIZE_tiempo, MPI_MIN, SUBCOMM_MPI, &ierr);
      dtcritico = newdtcritico;
#endif
      if (this->sgg.dt <= dtcritico) {
         sprintf(buff, "WIR_INFO: deltat for stability OK: %e", dtcritico);
         if ((this->control.layoutnumber == 0) && this->control.verbose) WarnErrReport(buff);
      } else {
         if (!(this->control.resume && this->control.permitscaling)) { // no abortasr solo advertir si permittivity scaling
#ifdef CompileWithMTLN
            sprintf(buff, "WIR_ERROR: Possibly UNSTABLE dt, make dt < %e", dtcritico);
#else
            sprintf(buff, "WIR_ERROR: Possibly UNSTABLE dt, decrease wire radius, number of parallel WIREs, use -stableradholland or make dt < %e", dtcritico);
#endif
            if (this->control.layoutnumber == 0) WarnErrReport(buff, true);
         } else {
            sprintf(buff, "WIR_WARNING: Resume and Pscaling with wires. Possibly UNSTABLE dt, decrease wire radius, number of parallel WIREs: dt is over %e", dtcritico);
            if (this->control.layoutnumber == 0) WarnErrReport(buff, false);
         }
      }
      // !!!
      // !!

   } // end subroutine initializeWires

   void initializeAnisotropic() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
      int rank;
#endif

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init Anisotropic...");
      print11(this->control.layoutnumber, dubuf);
      InitAnisotropic(this->sgg, this->media, this->thereAre.Anisotropic, this->thereAre.ThinSlot, this->eps0, this->mu0);
      l_auxinput = this->thereAre.Anisotropic || this->thereAre.ThinSlot;
      l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
      MPI_COMM_RANK(SUBCOMM_MPI, &rank, &ierr);
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
      if (l_auxoutput) {
         sprintf(dubuf, "----> there are Structured anisotropic elements");
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "----> no Structured anisotropic elements found");
         print11(this->control.layoutnumber, dubuf);
      }
   }

   void initializeSGBC() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
#endif

      if (this->control.sgbc) {
#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
         sprintf(dubuf, "Init Multi sgbc...");
         print11(this->control.layoutnumber, dubuf);
         Initsgbcs(this->sgg, this->media, Ex, Ey, Ez, Hx, Hy, Hz, IDxe, IDye, IDze, IDxh, Idyh, Idzh,
                   this->control.layoutnumber, this->control.num_procs, this->g, this->thereAre.sgbcs, this->control.resume,
                   this->control.sgbccrank, this->control.sgbcFreq, this->control.sgbcresol, this->control.sgbcdepth, this->control.sgbcDispersive,
                   this->eps0, this->mu0, this->control.simu_devia, this->control.stochastic);

         l_auxinput = this->thereAre.sgbcs;
         l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
         if (l_auxoutput) {
            sprintf(dubuf, "----> there are Structured sgbc elements");
            print11(this->control.layoutnumber, dubuf);
         } else {
            sprintf(dubuf, "----> no Structured sgbc elements found");
            print11(this->control.layoutnumber, dubuf);
         }
      }
   }

   void initializeMultiports() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;

#ifdef CompileWithNIBC
      if (this->control.mibc) {
#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
         sprintf(dubuf, "Init Multiports...");
         print11(this->control.layoutnumber, dubuf);
         InitMultiports(this->sgg, this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz, this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz, this->control.layoutnumber, this->control.num_procs, this->thereAre.Multiports, this->control.resume,
                        Idxe, Idye, Idze, this->control.NOcompomur, this->control.ADE, this->control.cfl, this->eps0, this->mu0);
         l_auxinput = this->thereAre.Multiports;
         l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
         if (l_auxoutput) {
            sprintf(dubuf, "----> there are Structured  multiport elements");
            print11(this->control.layoutnumber, dubuf);
         } else {
            sprintf(dubuf, "----> no Structured multiport elements found");
            print11(this->control.layoutnumber, dubuf);
         }
      }
#endif
   }

   void initializeEDispersives() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
#endif

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init EDispersives...");
      print11(this->control.layoutnumber, dubuf);
      InitEDispersives(this->sgg, this->media, this->thereAre.EDispersives, this->control.resume, this->g.g1, this->g.g2, ex, ey, ez);
      l_auxinput = this->thereAre.EDispersives;
      l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
      if (l_auxoutput) {
         sprintf(dubuf, "----> there are Structured Electric dispersive elements");
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "----> no Structured Electric dispersive elements found");
         print11(this->control.layoutnumber, dubuf);
      }
   }

   void initializeMDispersives() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
#endif

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init MDispersives...");
      print11(this->control.layoutnumber, dubuf);
      InitMDispersives(this->sgg, this->media, this->thereAre.MDispersives, this->control.resume, this->g.gm1, this->g.gm2, hx, hy, hz);
      l_auxinput = this->thereAre.MDispersives;
      l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
      if (l_auxoutput) {
         sprintf(dubuf, "----> there are Structured Magnetic dispersive elements");
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "----> no Structured Magnetic dispersive elements found");
         print11(this->control.layoutnumber, dubuf);
      }
   }

   void initializePlanewave() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
#endif

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init Multi Plane-Waves...");
      print11(this->control.layoutnumber, dubuf);
      InitPlaneWave(this->sgg, this->media, this->control.layoutnumber, this->control.num_procs, this->sinPML_fullsize, this->thereAre.PlaneWaveBoxes, this->control.resume, this->eps0, this->mu0);
      l_auxinput = this->thereAre.PlaneWaveBoxes;
      l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
      if (l_auxoutput) {
         sprintf(dubuf, "----> there are Plane Wave");
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "----> no Plane waves are found");
         print11(this->control.layoutnumber, dubuf);
      }
   }

   void initializeNodalSources() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
#endif

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init Nodal Sources...");
      print11(this->control.layoutnumber, dubuf);
      InitNodalSources(this->sgg, this->control.layoutnumber, this->sgg.NumNodalSources, this->sgg.NodalSource, this->sgg.Sweep, this->thereAre.NodalE, this->thereAre.NodalH);
      l_auxinput = this->thereAre.NodalH || this->thereAre.NodalE;
      l_auxoutput = l_auxinput;
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
      if (l_auxoutput) {
         sprintf(dubuf, "----> there are Structured Nodal sources");
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "----> no Structured Nodal sources are found");
         print11(this->control.layoutnumber, dubuf);
      }
   }

   void initializeObservation() {
      char dubuf[BUFSIZE];
      bool l_auxinput, l_auxoutput;
#ifdef CompileWithMPI
      int ierr;
#endif

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      sprintf(dubuf, "Init Observation...");
      print11(this->control.layoutnumber, dubuf);
      InitObservation(this->sgg, this->media, this->tag_numbers,
                      this->thereAre.Observation, this->thereAre.wires, this->thereAre.FarFields, this->initialtimestep, this->lastexecutedtime,
                      this->sinPML_fullsize, this->eps0, this->mu0, this->bounds, this->control);
      l_auxinput = this->thereAre.Observation || this->thereAre.FarFields;
      l_auxoutput = l_auxinput;

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&l_auxinput, &l_auxoutput, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#endif
      if (l_auxoutput) {
         sprintf(dubuf, "----> there are observation requests");
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "----> no observation requests are found");
         print11(this->control.layoutnumber, dubuf);
      }
   }

#ifdef CompileWithMPI
   void initializeMPI() {
      char dubuf[BUFSIZE];
      int ierr;
      if (this->control.num_procs > 1) {
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         sprintf(dubuf, "Init MPI MediaMatrix flush...");
         print11(this->control.layoutnumber, dubuf);
         InitMPI(this->sgg.sweep, this->sgg.alloc);
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         InitExtraFlushMPI(this->control.layoutnumber, this->sgg.sweep, this->sgg.alloc, this->sgg.med, this->sgg.nummedia, this->media.sggMiEz, this->media.sggMiHz);
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         FlushMPI_H(this->sgg.alloc, this->control.layoutnumber, this->control.num_procs, this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz);
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         FlushMPI_E(this->sgg.alloc, this->control.layoutnumber, this->control.num_procs, this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz);
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         sprintf(dubuf, "[OK]");
         print11(this->control.layoutnumber, dubuf);
      }

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! fin juego con fuego 210815

      // MPI initialization
      if (this->control.num_procs > 1) {
         sprintf(dubuf, "Init MPI Cray...");
         print11(this->control.layoutnumber, dubuf);
         InitMPI_Cray(this->control.layoutnumber, this->control.num_procs, this->sgg.sweep, this->sgg.alloc,
                      this->sgg.Border.IsDownPeriodic, this->sgg.Border.IsUpPeriodic,
                      Ex, Ey, Ez, Hx, Hy, Hz);
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         sprintf(dubuf, "[OK]");
         print11(this->control.layoutnumber, dubuf);

         // this modifies the initwires stuff and must be called after initwires (typically at the end)
         // llamalo siempre aunque no HAYA WIRES!!! para que no se quede colgado en hilos terminales
         if ((std::string(this->control.wiresflavor).find("holland") != std::string::npos) ||
             (std::string(this->control.wiresflavor).find("transition") != std::string::npos)) {
            sprintf(dubuf, "Init MPI Holland Wires...");
            print11(this->control.layoutnumber, dubuf);
            newInitWiresMPI(this->control.layoutnumber, this->thereAre.wires, this->control.num_procs, this->control.resume, this->sgg.sweep);
            MPI_Barrier(SUBCOMM_MPI, &ierr);
            sprintf(dubuf, "[OK]");
            print11(this->control.layoutnumber, dubuf);
         }

#ifdef CompileWithBerengerWires
         if (std::string(this->control.wiresflavor).find("berenger") != std::string::npos) {
            sprintf(dubuf, "Init MPI Multi-Wires...");
            print11(this->control.layoutnumber, dubuf);
            InitWiresMPI_Berenger(this->control.layoutnumber, this->thereAre.wires, this->control.num_procs, this->control.resume, this->sgg.sweep);
            MPI_Barrier(SUBCOMM_MPI, &ierr);
            sprintf(dubuf, "[OK]");
            print11(this->control.layoutnumber, dubuf);
         }
#endif
         // llamalo siempre para forzar los flush extra en caso de materiales anisotropos o multiport
         sprintf(dubuf, "Init Extra Flush MPI...");
         print11(this->control.layoutnumber, dubuf);
         InitExtraFlushMPI_Cray(this->control.layoutnumber, this->sgg.sweep, this->sgg.alloc, this->sgg.Med, this->sgg.NumMedia, this->media.sggMiez, this->media.sggMiHz,
                                Ex, Ey, Ez, Hx, Hy, Hz, this->thereAre.MURBorders);
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         sprintf(dubuf, "[OK]");
         print11(this->control.layoutnumber, dubuf);
      }


      // must be called now in case the MPI has changed the connectivity info
      if ((std::string(this->control.wiresflavor).find("holland") != std::string::npos) ||
          (std::string(this->control.wiresflavor).find("transition") != std::string::npos)) {
         ReportWireJunctions(this->control.layoutnumber, this->control.num_procs, this->thereAre.wires, this->sgg.Sweep[iHz]->ZI, this->sgg.Sweep[iHz]->ZE, this->control.groundwires, this->control.strictOLD, this->control.verbose);
      }

#ifdef CompileWithBerengerWires
      if (std::string(this->control.wiresflavor).find("berenger") != std::string::npos) {
         ReportWireJunctionsBerenger(this->control.layoutnumber, this->control.num_procs, this->thereAre.wires, this->sgg.Sweep[iHz]->ZI, this->sgg.Sweep[iHz]->ZE, this->control.groundwires, this->control.strictOLD, this->control.verbose);
         // dama no tenia el equivalente 050416
      }
#endif
#ifdef CompileWithSlantedWires
      if ((std::string(this->control.wiresflavor).find("slanted") != std::string::npos) || (std::string(this->control.wiresflavor).find("semistructured") != std::string::npos)) {
         // continue
      }
#endif

   }
#endif

#ifdef CompileWithMPI
   void flushMPIdata() {
      int ierr;
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      // !Flush all the MPI data (needed a initial flush for correct resuming)
      if (this->control.num_procs > 1) {
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         FlushMPI_H_Cray();
      }
      if ((std::string(this->control.wiresflavor).find("holland") != std::string::npos) ||
          (std::string(this->control.wiresflavor).find("transition") != std::string::npos)) {
         if ((this->control.num_procs > 1) && (this->thereAre.wires)) {
            newFlushWiresMPI(this->control.layoutnumber, this->control.num_procs);
         }
#ifdef CompileWithStochastic
         if (this->control.stochastic) {
            syncstoch_mpi_wires(this->control.simu_devia, this->control.layoutnumber, this->control.num_procs);
         }
#endif
      }

#ifdef CompileWithBerengerWires
      if (std::string(this->control.wiresflavor).find("berenger") != std::string::npos) {
         if ((this->control.num_procs > 1) && (this->thereAre.wires)) FlushWiresMPI_Berenger(this->control.layoutnumber, this->control.num_procs);
      }
#endif
   }
#endif

   void printSimulationStart() {
      char dubuf[BUFSIZE];
      tiempo_t time_out2;
#ifdef CompileWithMPI
      int ierr;
#endif

      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (this->control.resume) {
         sprintf(dubuf, "END PREPROCESSING. RESUMING simulation from n=%d", this->n);
         print11(this->control.layoutnumber, dubuf);
         sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
         print11(this->control.layoutnumber, dubuf);
      } else {
         sprintf(dubuf, "END PREPROCESSING. STARTING simulation.");
         print11(this->control.layoutnumber, dubuf);
         sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
         print11(this->control.layoutnumber, dubuf);
#ifdef CompileWithMPI
         MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
         get_secnds(time_out2);
         sprintf(dubuf, "Start Date/time %s/%s   %s:%s:%s",
                 time_out2.fecha.substr(6, 2).c_str(),
                 time_out2.fecha.substr(4, 2).c_str(),
                 time_out2.hora.substr(0, 2).c_str(),
                 time_out2.hora.substr(2, 2).c_str(),
                 time_out2.hora.substr(4, 2).c_str());
         print11(this->control.layoutnumber, dubuf);
         sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
         print11(this->control.layoutnumber, dubuf);
      }
   }

   void fillMtag(SGGFDTDINFO_t& sgg, bounds_t& b,
                 std::vector<std::vector<std::vector<int>>>& sggMtag,
                 std::vector<std::vector<std::vector<int>>>& sggMiHx,
                 std::vector<std::vector<std::vector<int>>>& sggMiHy,
                 std::vector<std::vector<std::vector<int>>>& sggMiHz,
                 std::vector<std::vector<std::vector<int>>>& sggMiEx,
                 std::vector<std::vector<std::vector<int>>>& sggMiEy,
                 std::vector<std::vector<std::vector<int>>>& sggMiEz,
                 taglist_t& tag_numbers) {

      // ------------------------>
      // type(SGGFDTDINFO_t), intent(in) :: sgg
      // type(bounds_t), intent( IN) :: b
      // integer(KIND = IKINDMTAG), dimension( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( INOUT) :: sggMtag
      // integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 : b%sggMiHx%NX-1 , 0 : b%sggMiHx%NY-1 , 0 : b%sggMiHx%NZ-1 )  , intent( IN   ) :: sggMiHx
      // integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 : b%sggMiHy%NX-1 , 0 : b%sggMiHy%NY-1 , 0 : b%sggMiHy%NZ-1 )  , intent( IN   ) :: sggMiHy
      // integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 : b%sggMiHz%NX-1 , 0 : b%sggMiHz%NY-1 , 0 : b%sggMiHz%NZ-1 )  , intent( IN   ) :: sggMiHz
      // integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 : b%sggMiEx%NX-1 , 0 : b%sggMiEx%NY-1 , 0 : b%sggMiEx%NZ-1 )  , intent( IN   ) :: sggMiEx
      // integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 : b%sggMiEy%NX-1 , 0 : b%sggMiEy%NY-1 , 0 : b%sggMiEy%NZ-1 )  , intent( IN   ) :: sggMiEy
      // integer(kind = INTEGERSIZEOFMEDIAMATRICES), dimension( 0 : b%sggMiEz%NX-1 , 0 : b%sggMiEz%NY-1 , 0 : b%sggMiEz%NZ-1 )  , intent( IN   ) :: sggMiEz
      // type(taglist_t) :: tag_numbers
      // ------------------------> Variables locales
      int i, j, k;
      int medio1, medio2, medio3, medio4, medio5;
      bool mediois1, mediois2, mediois3, mediois4;
      std::vector<int> lbx(3), lby(3), lbz(3);
      
      // Assuming tag_numbers has x, y, z members which are vectors or arrays
      // Fortran lbound usually returns the lower index. In C++ std::vector, lower is 0.
      // However, the code uses 0-based indexing in the loops below (1 to NX-1 implies 0 is boundary or unused, or NX is size)
      // Fortran: dimension( 0 : b%sggMiHx%NX-1 ) means size is NX. Indices 0..NX-1.
      // Loops: Do k=1,b%sweepHx%NZ ... Do i=1,b%sweepHx%NX
      // This suggests 1-based indexing for the loops, but array is 0-based?
      // Let's assume standard 1-based logic mapped to 0-based vector access if needed, 
      // but here the loops go 1 to N. If vector is 0-indexed, we access [i-1].
      // However, looking at `medio1 =sggMiEy(i,j,k)`, if i starts at 1, and vector is 0-indexed, this is out of bounds unless vector is 1-indexed (padded).
      // Given the complexity, I will assume the vectors are accessed directly with 1-based indices if they were resized to N+1, 
      // OR the loops are 0-based in C++. 
      // Fortran: `dimension( 0 : ... )` usually means index 0 is valid.
      // Loop `Do i=1, ...` skips index 0.
      // I will translate loops to 1-based if the Fortran code implies it, but C++ vectors are 0-based.
      // To preserve logic exactly: if Fortran uses 1..N, C++ should use 0..N-1.
      // But Fortran array is 0..N-1. So Fortran loop 1..N accesses indices 1..N. Index N is out of bounds for size N (0..N-1).
      // Wait, `dimension( 0 : b%sggMiHx%NX-1 )` has size `NX`. Indices `0` to `NX-1`.
      // Loop `Do i=1,b%sweepHx%NX`. If `b%sweepHx%NX` equals `NX`, then `i` goes up to `NX`.
      // Accessing `sggMiEy(i,j,k)` where `i=NX` is out of bounds for a vector of size `NX` (indices 0..NX-1).
      // This suggests either `NX` in loop is different from array dimension size, or there's padding.
      // I will assume the loops are meant to iterate over valid indices. 
      // Let's assume 1-based indexing in Fortran maps to 0-based in C++ by subtracting 1.
      
      // Re-reading: `dimension( 0 : b%sggMiHx%NX-1 )`. Size is `NX`.
      // Loop `Do i=1,b%sweepHx%NX`.
      // If `b%sweepHx%NX` is the dimension size, then loop goes 1..Size.
      // Access `sggMiEy(i,j,k)`.
      // This is likely a bug in my interpretation or the Fortran code relies on 1-based arrays (common in older Fortran if not declared 0-based explicitly, but here it IS declared 0-based).
      // Actually, if declared `0:N-1`, valid indices are `0..N-1`.
      // If loop is `1:N`, it accesses `1..N`. `N` is invalid.
      // Perhaps `b%sweepHx%NX` is not the same as the array dimension?
      // Or perhaps the array is 1-based in reality despite the declaration?
      // I will stick to the loop bounds provided: `1` to `N`.
      // And access `vec[i-1]` to map to 0-based C++ vector if the vector is 0-based.
      // BUT, if the Fortran code accesses `i` directly, and `i` starts at 1, it expects 1-based indexing.
      // I will assume the C++ vectors are padded or accessed with 1-based logic (e.g. `vec[i]`).
      // To be safe and "preserve names/logic", I will use 1-based indexing for loops and direct access, assuming the vectors are sized appropriately (N+1).

      lbx = tag_numbers.face.x.lbound();
      lby = tag_numbers.face.y.lbound();
      lbz = tag_numbers.face.z.lbound();

      mediois3 = true;
      mediois4 = true;
#ifdef CompileWithOpenMP
#pragma omp parallel for default(shared) private(i,j,k,medio1,medio2,medio3,medio4,medio5,mediois1,mediois2,mediois3,mediois4)
#endif
      for (k = 1; k <= b.sweepHx.NZ; ++k) {
         for (j = 1; j <= b.sweepHx.NY; ++j) {
            for (i = 1; i <= b.sweepHx.NX; ++i) {
               medio1 = sggMiEy[i][j][k];
               medio2 = sggMiEy[i][j][k + 1];
               medio3 = sggMiEz[i][j][k];
               medio4 = sggMiEz[i][j + 1][k];
               medio5 = sggMiHx[i][j][k];
               mediois1 = (medio5 == 1) && (medio1 != 1) && (medio2 != 1) && (medio3 == 1) && (medio4 == 1);
               mediois2 = (medio5 == 1) && (medio3 != 1) && (medio4 != 1) && (medio1 == 1) && (medio2 == 1);
               mediois3 = true; // .not.((medio5==1).and.(((sggMiHx(i-1,j,k)/=1).or.(sggMiHx(i+1,j,k)/=1)))) !esta condicion en realidad no detecta alabeos de una celda que siendo slots son acoples de un agujerito solo en el peor de los casos
               if ((mediois1 || mediois2) && mediois3) {
                  // ... rest of the block not provided in chunk
               }
            }
         }
      }
   }

// Continuation chunk

                      // solo lo hace con celdas de vacio porque en particular el mismo medio sgbc con diferentes orientaciones tiene distintos indices de medio y lo activaria erroneamente si lo hago para todos los medios
                      tag_numbers.face.x(i + lbx(0) - 1, j + lbx(1) - 1, k + lbx(2) - 1) = -ibset(std::abs(tag_numbers.face.x(i + lbx(0) - 1, j + lbx(1) - 1, k + lbx(2) - 1)), 3);
                      // ojo no cambiar: interacciona con observation tags 141020 !151020 a efectos de mapvtk el signo importa
                  }
               }
            }
         }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#pragma omp parallel do default(shared) private(i, j, k, medio1, medio2, medio3, medio4, medio5, mediois1, mediois2, mediois3, mediois4)
#endif
         for (int k = 1; k <= b.sweepHy.NZ; ++k) {
            for (int j = 1; j <= b.sweepHy.NY; ++j) {
               for (int i = 1; i <= b.sweepHy.NX; ++i) {
                  medio1 = sggMiEz(i, j, k);
                  medio2 = sggMiEz(i + 1, j, k);
                  medio3 = sggMiEx(i, j, k);
                  medio4 = sggMiEx(i, j, k + 1);
                  medio5 = sggMiHy(i, j, k);
                  mediois1 = (medio5 == 1) && (medio1 != 1) && (medio2 != 1) && (medio3 == 1) && (medio4 == 1);
                  mediois2 = (medio5 == 1) && (medio3 != 1) && (medio4 != 1) && (medio1 == 1) && (medio2 == 1);
                  mediois3 = true; //.not.((medio5==1).and.(((sggMiHy(i,j-1,k)/=1).or.(sggMiHy(i,j+1,k)/=1))))
                  if ((mediois1 || mediois2) && mediois3) {
                     tag_numbers.face.y(i + lby(0) - 1, j + lby(1) - 1, k + lby(2) - 1) = -ibset(std::abs(tag_numbers.face.y(i + lby(0) - 1, j + lby(1) - 1, k + lby(2) - 1)), 4);
                  }
               }
            }
         }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#pragma omp parallel do default(shared) private(i, j, k, medio1, medio2, medio3, medio4, medio5, mediois1, mediois2, mediois3, mediois4)
#endif
         for (int k = 1; k <= b.sweepHz.NZ; ++k) {
            for (int j = 1; j <= b.sweepHz.NY; ++j) {
               for (int i = 1; i <= b.sweepHz.NX; ++i) {
                  medio1 = sggMiEx(i, j, k);
                  medio2 = sggMiEx(i, j + 1, k);
                  medio3 = sggMiEy(i, j, k);
                  medio4 = sggMiEy(i + 1, j, k);
                  medio5 = sggMiHz(i, j, k);
                  mediois1 = (medio5 == 1) && (medio1 != 1) && (medio2 != 1) && (medio3 == 1) && (medio4 == 1);
                  mediois2 = (medio5 == 1) && (medio3 != 1) && (medio4 != 1) && (medio1 == 1) && (medio2 == 1);
                  mediois3 = true; //.not.((medio5==1).and.(((sggMiHz(i,j,k-1)/=1).or.(sggMiHz(i,j,k+1)/=1))))
                  if ((mediois1 || mediois2) && mediois3) {
                     tag_numbers.face.z(i + lbz(0) - 1, j + lbz(1) - 1, k + lbz(2) - 1) = -ibset(std::abs(tag_numbers.face.z(i + lbz(0) - 1, j + lbz(1) - 1, k + lbz(2) - 1)), 5);
                  }
               }
            }
         }
#ifdef CompileWithOpenMP
#pragma omp end parallel do
#endif
         return;
      } // end subroutine fillMtag

      void crea_timevector(SGGFDTDINFO_t& sgg, int lastexecutedtimestep, int finaltimestep, double lastexecutedtime) {
         sgg.tiempo.resize(finaltimestep + 3); // Assuming 1-based indexing in Fortran, mapped to 0-based or adjusted in C++ logic. 
         // Note: Fortran allocate (sgg%tiempo(lastexecutedtimestep:finaltimestep+2)) implies size is finaltimestep+2 - lastexecutedtimestep + 1.
         // In C++, we usually use 0-based or resize to fit. Let's assume the vector is resized to accommodate the indices.
         // If the vector is 0-indexed in C++, we need to adjust. However, keeping names, let's assume the vector is large enough or resized.
         // Standard practice: resize to finaltimestep + 3 to safely access up to finaltimestep+2 if 0-based, or just resize to size.
         // Given the complexity of index mapping, we'll resize to max index + 1.
         if (sgg.tiempo.size() < finaltimestep + 3) {
             sgg.tiempo.resize(finaltimestep + 3);
         }
         
         sgg.tiempo[lastexecutedtimestep] = lastexecutedtime;
         for (int i = lastexecutedtimestep + 1; i <= finaltimestep + 2; ++i) {
            sgg.tiempo[i] = sgg.tiempo[i - 1] + sgg.dt; // equiespaciados por defecto !luego los modifica prescale
         }
         return;
      }

   } // end subroutine solver_init (contextual closing brace from previous chunk likely, but here we just output code)

   void solver_run(solver_t& this_obj) {

      double* Ex = this_obj.Ex.data();
      double* Ey = this_obj.Ey.data();
      double* Ez = this_obj.Ez.data();
      double* Hx = this_obj.Hx.data();
      double* Hy = this_obj.Hy.data();
      double* Hz = this_obj.Hz.data();
      
      int* Idxe = this_obj.Idxe.data();
      int* Idye = this_obj.Idye.data();
      int* Idze = this_obj.Idze.data();
      int* Idxh = this_obj.Idxh.data();
      int* Idyh = this_obj.Idyh.data();
      int* Idzh = this_obj.Idzh.data();
      double* dxe = this_obj.dxe.data();
      double* dye = this_obj.dye.data();
      double* dze = this_obj.dze.data();
      double* dxh = this_obj.dxh.data();
      double* dyh = this_obj.dyh.data();
      double* dzh = this_obj.dzh.data();

      bool call_timing, l_aux, flushFF, somethingdone, newsomethingdone;
      int i;
      double pscale_alpha;
      double at;
      char dubuf[bufsize];
#ifdef CompileWithMPI
      int ierr;
#endif
#ifdef CompileWithProfiling
      nvtxStartRange("Antes del bucle N");
#endif

      this_obj.still_planewave_time = true; // inicializacion de la variable 
      flushFF = false;
      pscale_alpha = 1.0; // se le entra con 1.0 

      // Pointers are already assigned above via .data()

      while (this_obj.n <= this_obj.control.finaltimestep) {
      
         this_obj.step();
         updateAndFlush(this_obj, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh);

         if (this_obj.n >= this_obj.n_info) {
             call_timing = true;
         } else {
             call_timing = false;
         }
#ifdef CompileWithMPI
         l_aux = call_timing;
         MPI_AllReduce(&l_aux, &call_timing, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
         MPI_Barrier(MPI_COMM_WORLD, &ierr); // 050619 incluido problemas stochastic stopflusing
#endif
         
         if (call_timing) {
            Timing(this_obj.sgg, this_obj.bounds, this_obj.n, this_obj.n_info, this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.control.maxCPUtime, this_obj.control.flushsecondsFields, this_obj.control.flushsecondsData, this_obj.initialtimestep,
            this_obj.control.finaltimestep, this_obj.perform, this_obj.parar, false,
            Ex, Ey, Ez, this_obj.everflushed, this_obj.control.nentradaroot, this_obj.control.maxSourceValue, this_obj.control.opcionestotales, this_obj.control.simu_devia, this_obj.control.dontwritevtk, this_obj.control.permitscaling);

            if (!this_obj.parar) { // !!! si es por parada se gestiona al final
!!!!! si esta hecho lo flushea todo pero poniendo de acuerdo a todos los mpi
                for (i = 1; i <= this_obj.sgg.NumberRequest; ++i) {
                   if (this_obj.sgg.Observation[i].done && (!this_obj.sgg.Observation[i].flushed)) {
                      this_obj.perform.flushXdmf = true;
                      this_obj.perform.flushVTK = true;
                   }
                }
#ifdef CompileWithMPI
                l_aux = this_obj.perform.flushVTK;
                MPI_AllReduce(&l_aux, &this_obj.perform.flushVTK, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                //
                l_aux = this_obj.perform.flushXdmf;
                MPI_AllReduce(&l_aux, &this_obj.perform.flushXdmf, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                //
                l_aux = this_obj.perform.flushDATA;
                MPI_AllReduce(&l_aux, &this_obj.perform.flushDATA, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                //
                l_aux = this_obj.perform.flushFIELDS;
                MPI_AllReduce(&l_aux, &this_obj.perform.flushFIELDS, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                //
                l_aux = this_obj.perform.postprocess;
                MPI_AllReduce(&l_aux, &this_obj.perform.postprocess, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
#endif
!!!!!!!!!!!!
                if (this_obj.perform.flushFIELDS) {
                   sprintf(dubuf, "%s%s%s", SEPARADOR, trim(adjustl(this_obj.control.nentradaroot)), separador);
                   print11(this_obj.control.layoutnumber, dubuf);
                   sprintf(dubuf, "INIT FLUSHING OF RESTARTING FIELDS n=%d", this_obj.n);
                   print11(this_obj.control.layoutnumber, dubuf);
                   flush_and_save_resume(this_obj.sgg, this_obj.bounds, this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.control.nentradaroot, this_obj.control.nresumeable2, this_obj.thereare, this_obj.n, this_obj.eps0, this_obj.mu0, this_obj.everflushed,
                   Ex, Ey, Ez, Hx, Hy, Hz, this_obj.control.wiresflavor, this_obj.control.simu_devia, this_obj.control.stochastic);
#ifdef CompileWithMPI
                   MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
                   sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                   print11(this_obj.control.layoutnumber, dubuf);
                   sprintf(dubuf, "DONE FLUSHING OF RESTARTING FIELDS n=%d", this_obj.n);
                   print11(this_obj.control.layoutnumber, dubuf);
                   sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                   print11(this_obj.control.layoutnumber, dubuf);
                }
                if (this_obj.perform.isFlush()) {
                      //
                      flushFF = this_obj.perform.postprocess;
                      if (this_obj.thereAre.FarFields && flushFF) {
                          sprintf(dubuf, " INIT OBSERVATION DATA FLUSHING and Near-to-Far field n= %d", this_obj.n);
                      } else {
                          sprintf(dubuf, " INIT OBSERVATION DATA FLUSHING n= %d", this_obj.n);
                      }
                      print11(this_obj.control.layoutnumber, SEPARADOR);
                      print11(this_obj.control.layoutnumber, separador);
                      print11(this_obj.control.layoutnumber, separador);
                      print11(this_obj.control.layoutnumber, dubuf);
                      print11(this_obj.control.layoutnumber, SEPARADOR);
                      print11(this_obj.control.layoutnumber, separador);
                      print11(this_obj.control.layoutnumber, separador);
    !!
                      if (this_obj.thereAre.Observation) FlushObservationFiles(this_obj.sgg, this_obj.ini_save, this_obj.n, this_obj.control.layoutnumber, this_obj.control.num_procs, dxe, dye, dze, dxh, dyh, dzh, this_obj.bounds, this_obj.control.singlefilewrite, this_obj.control.facesNF2FF, flushFF);
                      !!
#ifdef CompileWithMPI
                      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
                      if (this_obj.thereAre.FarFields && flushFF) {
                          sprintf(dubuf, " Done OBSERVATION DATA FLUSHED and Near-to-Far field n= %d", this_obj.n);
                      } else {
                          sprintf(dubuf, " Done OBSERVATION DATA FLUSHED n= %d", this_obj.n);
                      }
                      print11(this_obj.control.layoutnumber, SEPARADOR);
                      print11(this_obj.control.layoutnumber, separador);
                      print11(this_obj.control.layoutnumber, separador);
                      print11(this_obj.control.layoutnumber, dubuf);
                      print11(this_obj.control.layoutnumber, SEPARADOR);
                      print11(this_obj.control.layoutnumber, separador);
                      print11(this_obj.control.layoutnumber, separador);
    !
                      if (this_obj.perform.postprocess) {
                         sprintf(dubuf, "Postprocessing frequency domain probes, if any, at n= %d", this_obj.n);
                         print11(this_obj.control.layoutnumber, dubuf);
                         sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                         print11(this_obj.control.layoutnumber, dubuf);
                         somethingdone = false;
                         at = this_obj.n * this_obj.sgg.dt;
                         if (this_obj.thereAre.Observation) PostProcessOnthefly(this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.sgg, this_obj.control.nentradaroot, at, somethingdone, this_obj.control.niapapostprocess, this_obj.control.forceresampled);
#ifdef CompileWithMPI
                         MPI_Barrier(SUBCOMM_MPI, &ierr);
                         MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                         somethingdone = newsomethingdone;
#endif
                         if (somethingdone) {
                           sprintf(dubuf, "End Postprocessing frequency domain probes.");
                           print11(this_obj.control.layoutnumber, dubuf);
                           sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                           print11(this_obj.control.layoutnumber, dubuf);
                         } else {
                           sprintf(dubuf, "No frequency domain probes snapshots found to be postrocessed");
                           print11(this_obj.control.layoutnumber, dubuf);
                           sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                           print11(this_obj.control.layoutnumber, dubuf);
                         }
                      }
                  !!       
                      if (this_obj.perform.flushvtk) {   
                         sprintf(dubuf, " Post-processing .vtk files n= %d", this_obj.n);
                         print11(this_obj.control.layoutnumber, SEPARADOR);
                         print11(this_obj.control.layoutnumber, separador);
                         print11(this_obj.control.layoutnumber, separador);
                         print11(this_obj.control.layoutnumber, dubuf);
                         print11(this_obj.control.layoutnumber, SEPARADOR);
                         print11(this_obj.control.layoutnumber, separador);
                         print11(this_obj.control.layoutnumber, separador);
                         somethingdone = false;
                         if (this_obj.thereAre.Observation) createvtkOnTheFly(this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.sgg, this_obj.control.vtkindex, somethingdone, this_obj.control.mpidir, this_obj.media.sggMtag, this_obj.control.dontwritevtk);
#ifdef CompileWithMPI
                         MPI_Barrier(SUBCOMM_MPI, &ierr);
                         MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                         somethingdone = newsomethingdone;
#endif
                          if (somethingdone) {
                                sprintf(dubuf, "End flushing .vtk snapshots");
                                print11(this_obj.control.layoutnumber, dubuf);
                                sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                                print11(this_obj.control.layoutnumber, dubuf);
                          } else {
                                sprintf(dubuf, "No .vtk snapshots found to be flushed");
                                print11(this_obj.control.layoutnumber, dubuf);
                                sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                                print11(this_obj.control.layoutnumber, dubuf);
                          }
                      }  
                         if (this_obj.perform.flushXdmf) {
                            sprintf(dubuf, " Post-processing .xdmf files n= %d", this_obj.n);
                            print11(this_obj.control.layoutnumber, SEPARADOR);
                            print11(this_obj.control.layoutnumber, separador);
                            print11(this_obj.control.layoutnumber, separador);
                            print11(this_obj.control.layoutnumber, dubuf);
                            print11(this_obj.control.layoutnumber, SEPARADOR);
                            print11(this_obj.control.layoutnumber, separador);
                            print11(this_obj.control.layoutnumber, separador);
                            somethingdone = false;

                            if (this_obj.thereAre.Observation) createxdmfOnTheFly(this_obj.sgg, this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.control.vtkindex, this_obj.control.createh5bin, somethingdone, this_obj.control.mpidir);                          
                            if (this_obj.control.createh5bin) createh5bintxt(this_obj.sgg, this_obj.control.layoutnumber, this_obj.control.num_procs); // lo deben llamar todos haya on on this%thereAre%observation

#ifdef CompileWithMPI
                        MPI_Barrier(SUBCOMM_MPI, &ierr);
                        MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
                        somethingdone = newsomethingdone;
#endif
                            if (somethingdone) {
                                      sprintf(dubuf, "End flushing .xdmf snapshots");
                                      print11(this_obj.control.layoutnumber, dubuf);
                                      sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                                      print11(this_obj.control.layoutnumber, dubuf);
                             } else {
                                      sprintf(dubuf, "No .xdmf snapshots found to be flushed");
                                      print11(this_obj.control.layoutnumber, dubuf);
                                      sprintf(dubuf, "%s%s%s", SEPARADOR, separador, separador);
                                      print11(this_obj.control.layoutnumber, dubuf);
                            }
                      end if

#ifdef CompileWithMPI
                     MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
                 end if // del if (this%performflushDATA.or....
    !


                  if (this_obj.control.singlefilewrite && this_obj.perform.Unpack) singleUnpack(this_obj);
                  if ((this_obj.control.singlefilewrite && this_obj.perform.Unpack) || this_obj.perform.isFlush()) {
                     sprintf(dubuf, " Continuing simulation at n= %d", this_obj.n);
                     print11(this_obj.control.layoutnumber, SEPARADOR);
                     print11(this_obj.control.layoutnumber, separador);
                     print11(this_obj.control.layoutnumber, separador);
                     print11(this_obj.control.layoutnumber, dubuf);
                     print11(this_obj.control.layoutnumber, SEPARADOR);
                     print11(this_obj.control.layoutnumber, separador);
                     print11(this_obj.control.layoutnumber, separador);
                  end if

                end if // !!!del if (.not.this%parar)
             end if // !!!del if(n >= n_info
!          !!!!!!!!all the previous must be together
              
         this_obj.control.fatalerror = false;
         if (this_obj.parar) {
             this_obj.control.fatalerror = true;
             break; // exit ciclo_temporal
         }
#ifdef CompileWithPrescale
         if (this_obj.control.permitscaling) {
#ifndef miguelPscaleStandAlone
            if ((this_obj.sgg.tiempo[this_obj.n] >= this_obj.EpsMuTimeScale_input_parameters.tini) &&
                (this_obj.sgg.tiempo[this_obj.n] <= this_obj.EpsMuTimeScale_input_parameters.tend)) {
#endif
             updateconstants(this_obj.sgg, this_obj.n, this_obj.thereare, this_obj.g,
                               Idxe, Idye, Idze, Idxh, Idyh, Idzh, // needed by CPML to be updated
                               this_obj.control.sgbc, this_obj.control.mibc, input_conformal_flag,
                               this_obj.control.wiresflavor, this_obj.control.wirecrank, this_obj.control.fieldtotl,
                               this_obj.control.sgbcDispersive, this_obj.control.finaltimestep,
                               this_obj.eps0, this_obj.mu0,
                               this_obj.control.simu_devia,
                               this_obj.EpsMuTimeScale_input_parameters, pscale_alpha, this_obj.still_planewave_time
#ifdef CompileWithMPI
                               , this_obj.control.layoutnumber, this_obj.control.num_procs
#endif
                               , this_obj.control.stochastic, this_obj.control.verbose);
#ifndef miguelPscaleStandAlone
            }
#endif
         }
#endif
         // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         // !!!  Increase time step
         // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         // write(*write(*,*) 'timestepping: ', n
         this_obj.n = this_obj.n + 1; // sube de iteracion
      } // end do ciclo_temporal ! End of the time-stepping loop

   } // end subroutine solver_run

   void updateAndFlush(solver_t& this_obj, double* Ex, double* Ey, double* Ez, double* Hx, double* Hy, double* Hz, double* dxe, double* dye, double* dze, double* dxh, double* dyh, double* dzh) {
      int mindum;
      if (this_obj.thereAre.Observation) {
         UpdateObservation(this_obj.sgg, this_obj.media, this_obj.tag_numbers, this_obj.n, this_obj.ini_save, Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh, this_obj.control.wiresflavor, this_obj.sinPML_fullsize, this_obj.control.wirecrank, this_obj.control.noconformalmapvtk, this_obj.bounds);
         if (this_obj.n >= this_obj.ini_save + BuffObse) {
            mindum = std::min(this_obj.control.finaltimestep, this_obj.ini_save + BuffObse);
            FlushObservationFiles(this_obj.sgg, this_obj.ini_save, mindum, this_obj.control.layoutnumber, this_obj.control.num_procs, dxe, dye, dze, dxh, dyh, dzh, this_obj.bounds, this_obj.control.singlefilewrite, this_obj.control.facesNF2FF, false); // no se flushean los farfields ahora
         }
      }
   }

   void singleUnpack(solver_t& this_obj) {
      char dubuf[BUFSIZE];
      bool somethingdone;
      double at;
#ifdef CompileWithMPI
      int ierr;
#endif
      print11(this_obj.control.layoutnumber, SEPARADOR);
      print11(this_obj.control.layoutnumber, separador);
      print11(this_obj.control.layoutnumber, separador);
      sprintf(dubuf, " Unpacking .bin files and prostprocessing them at n= %d", this_obj.n);
      print11(this_obj.control.layoutnumber, dubuf);
      print11(this_obj.control.layoutnumber, SEPARADOR);
      print11(this_obj.control.layoutnumber, separador);
      print11(this_obj.control.layoutnumber, separador);
      if (this_obj.thereAre.Observation) unpacksinglefiles(this_obj.sgg, this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.control.singlefilewrite, this_obj.initialtimestep, this_obj.control.resume); // dump the remaining to disk
      somethingdone = false;
      if (this_obj.control.singlefilewrite && this_obj.perform.Unpack) {
         at = this_obj.n * this_obj.sgg.dt;
         if (this_obj.thereAre.Observation) PostProcessOnthefly(this_obj.control.layoutnumber, this_obj.control.num_procs, this_obj.sgg, this_obj.control.nentradaroot, at, somethingdone, this_obj.control.niapapostprocess, this_obj.control.forceresampled);
      }
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
      MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
      somethingdone = newsomethingdone;
#endif
      sprintf(dubuf, " Done Unpacking .bin files and prostprocessing them at n= %d", this_obj.n);
      print11(this_obj.control.layoutnumber, SEPARADOR);
      print11(this_obj.control.layoutnumber, separador);
      print11(this_obj.control.layoutnumber, separador);
      print11(this_obj.control.layoutnumber, dubuf);
      print11(this_obj.control.layoutnumber, SEPARADOR);
      print11(this_obj.control.layoutnumber, separador);
      print11(this_obj.control.layoutnumber, separador);

   } // end subroutine singleUnpack

   } // end subroutine solver_run

   void step(solver_t& this_obj) {
      bool planewave_switched_off = false, thereareplanewave;

#ifdef CompileWithMPI
      int ierr;
#endif

      flushPlanewaveOff(planewave_switched_off, this_obj.still_planewave_time, thereareplanewave, this_obj);
      this_obj.AdvanceAnisotropicE();
      this_obj.advanceE();

      this_obj.advanceWiresE();
      this_obj.advancePMLE();
#ifdef CompileWithNIBC
      if (this_obj.thereAre.Multiports && (this_obj.control.mibc)) AdvanceMultiportE(this_obj.sgg.alloc, this_obj.Ex, this_obj.Ey, this_obj.Ez);
#endif
      this_obj.AdvancesgbcE();
      this_obj.advanceLumpedE();
      this_obj.advanceEDispersiveE();
      this_obj.advancePlaneWaveE();
      this_obj.advanceNodalE();

#ifdef CompileWithMPI
      if (this_obj.control.num_procs > 1) {
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         FlushMPI_E_Cray();
      }
#endif

      this_obj.advanceAnisotropicH();
      this_obj.advanceH();
      this_obj.advancePMLbodyH();
      this_obj.AdvanceMagneticCPML();
      this_obj.MinusCloneMagneticPMC();
      this_obj.CloneMagneticPeriodic();
      this_obj.AdvancesgbcH();
      this_obj.AdvanceMDispersiveH();
#ifdef CompileWithNIBC
      if (this_obj.thereAre.Multiports && (this_obj.control.mibc))
         AdvanceMultiportH(this_obj.sgg.alloc, this_obj.Hx, this_obj.Hy, this_obj.Hz,
                           this_obj.Ex, this_obj.Ey, this_obj.Ez,
                           this_obj.Idxe, this_obj.Idye, this_obj.Idze,
                           this_obj.media.sggMiHx, this_obj.media.sggMiHy, this_obj.media.sggMiHz,
                           this_obj.g.gm2, this_obj.sgg.nummedia, this_obj.control.conformalskin);
#endif
      this_obj.advancePlaneWaveH();
      this_obj.advanceNodalH();
      this_obj.advanceWiresH();
      this_obj.MinusCloneMagneticPMC();
      this_obj.CloneMagneticPeriodic();

#ifdef CompileWithMPI
      // Flush all the MPI (esto estaba justo al principo del bucle temporal diciendo que era necesario para correcto resuming)
      // lo he movido aqui a 16/10/2012 porque el farfield necesita tener los campos magneticos correctos
      // e intuyo que el Bloque current tambien a tenor del comentario siguiente
      // Incluyo un flush inicial antes de entrar al bucle para que el resuming sea correcto
      if (this_obj.control.num_procs > 1) {
         MPI_Barrier(SUBCOMM_MPI, &ierr);
         FlushMPI_H_Cray();
      }
      if ((trim(adjustl(this_obj.control.wiresflavor)) == "holland") ||
            (trim(adjustl(this_obj.control.wiresflavor)) == "transition")) {
         if ((this_obj.control.num_procs > 1) && (this_obj.thereAre.wires)) newFlushWiresMPI(this_obj.control.layoutnumber, this_obj.control.num_procs);
#ifdef CompileWithStochastic
         if (this_obj.control.stochastic) syncstoch_mpi_wires(this_obj.control.simu_devia, this_obj.control.layoutnumber, this_obj.control.num_procs);
#endif
      }
#ifdef CompileWithBerengerWires
      if (trim(adjustl(this_obj.control.wiresflavor)) == "berenger") {
         if ((this_obj.control.num_procs > 1) && (this_obj.thereAre.wires)) FlushWiresMPI_Berenger(this_obj.control.layoutnumber, this_obj.control.num_procs);
      }
#endif
#endif

      // no se si el orden wires - sgbcs del sync importa 150519
#ifdef CompileWithMPI
#ifdef CompileWithStochastic
         if (this_obj.control.stochastic) syncstoch_mpi_sgbcs(this_obj.control.simu_devia, this_obj.control.layoutnumber, this_obj.control.num_procs);
#endif    
#endif

#ifdef CompileWithMPI
#ifdef CompileWithStochastic
         if (this_obj.control.stochastic) syncstoch_mpi_lumped(this_obj.control.simu_devia, this_obj.control.layoutnumber, this_obj.control.num_procs);
#endif    
#endif 
      this_obj.advanceMagneticMUR();

   } // end subroutine step

   void flushPlanewaveOff(bool& pw_switched_off, bool& pw_still_time, bool& pw_thereAre, solver_t& this_obj) {
      bool pw_still_time_aux, pw_thereAre_aux;
      int ierr;
      char dubuf[bufsize];
      if (!pw_switched_off) {
         pw_still_time = pw_still_time && this_obj.thereAre.PlaneWaveBoxes;
         pw_thereAre = this_obj.thereAre.PlaneWaveBoxes;
#ifdef CompileWithMPI
         if (this_obj.control.num_procs > 1) {
            pw_still_time_aux = pw_still_time;
            MPI_AllReduce(&pw_still_time_aux, &pw_still_time, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
            pw_thereAre_aux = pw_thereAre;
            MPI_AllReduce(&pw_thereAre_aux, &pw_thereAre, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
         }
#endif
         if (!pw_still_time) {
            pw_switched_off = true;
            sprintf(dubuf, "Switching plane-wave off at n=%d", this_obj.n);
            if (pw_thereAre) print11(this_obj.control.layoutnumber, dubuf);
         }
      }
   } // end subroutine flushPlanewaveOff

#ifdef CompileWithMPI
void solver_t::init_MPIConformalProbes() {
    int group_conformalprobes_dummy = 0;
    int ierr = 0;
    // sgg250424 niapa para que funcionen sondas conformal mpi
    // todos deben crear el subcomunicador mpi una sola vez   
    if (input_conformal_flag) {
        SUBCOMM_MPI_conformal_probes = 1;
        MPI_conformal_probes_root = this->control.layoutnumber;
    } else {
        SUBCOMM_MPI_conformal_probes = 0;
        MPI_conformal_probes_root = -1;
    }
    MPIinitSubcomm(this->control.layoutnumber, this->control.num_procs, SUBCOMM_MPI_conformal_probes,
                   MPI_conformal_probes_root, group_conformalprobes_dummy);
    // print *,'-----creating--->',this%control%layoutnumber,SIZE,SUBCOMM_MPI_conformal_probes,MPI_conformal_probes_root
    MPI_BARRIER(SUBCOMM_MPI, ierr);
    // !!!no lo hago pero al salir deberia luego destruir el grupo call MPI_Group_free(output(ii)%item(i)%MPIgroupindex,ierr)                   
}
#endif

void solver_t::advanceE() {
#ifdef CompileWithProfiling
    nvtxStartRange("Antes del bucle EX");
#endif
    advanceEx(this->media.sggMiEx);
#ifdef CompileWithProfiling
    nvtxEndRange();

    nvtxStartRange("Antes del bucle EY");
#endif
    advanceEy(this->media.sggMiEy);

#ifdef CompileWithProfiling
    nvtxEndRange();

    nvtxStartRange("Antes del bucle EZ");
#endif
    advanceEz(this->media.sggMiEz);
#ifdef CompileWithProfiling
    nvtxEndRange();
#endif
}

void solver_t::advanceEx(const std::vector<std::vector<std::vector<integersizeofmediamatrices>>>& sggMiEx) {
    // Map pointers to local references/vectors for easier access
    
    for (k = 0; k < this->bounds.sweepEx.NZ; ++k) {
        for (j = 0; j < this->bounds.sweepEx.NY; ++j) {
            for (i = 0; i < this->bounds.sweepEx.NX; ++i) {
                Idzhk = this->Idzh[k]; // Idzh(0:dzh%NZ-1) => this%Idzh. Idzh is 1D.
                Idyhj = this->Idyh[j]; // Idyh(0:dyh%NY-1) => this%Idyh.
                
                medio = sggMiEx[i][j][k]; // Assuming sggMiEx is passed as 0-indexed 3D vector
                
                // Ex(i,j,k) = ...
                // Using 0-based indexing for this->Ex
                this->Ex[i][j][k] = this->g.g1(medio) * this->Ex[i][j][k] + 
                                    this->g.g2(medio) * 
                                    ((this->Hz[i][j][k] - this->Hz[i][j-1][k]) * Idyhj - 
                                     (this->Hy[i][j][k] - this->Hy[i][j][k-1]) * Idzhk);
            }
        }
    }
}

void solver_t::advanceEy(const std::vector<std::vector<std::vector<integersizeofmediamatrices>>>& sggMiEy) {
    double Idzhk;
    int i, j, k;
    integersizeofmediamatrices medio;

    for (k = 0; k < this->bounds.sweepEy.NZ; ++k) {
        for (j = 0; j < this->bounds.sweepEy.NY; ++j) {
            for (i = 0; i < this->bounds.sweepEy.NX; ++i) {
                Idzhk = this->Idzh[k];
                medio = sggMiEy[i][j][k];
                
                this->Ey[i][j][k] = this->g.g1(medio) * this->Ey[i][j][k] + 
                                    this->g.g2(medio) * 
                                    ((this->Hx[i][j][k] - this->Hx[i][j][k-1]) * Idzhk - 
                                     (this->Hz[i][j][k] - this->Hz[i-1][j][k]) * this->Idxh[i]);
            }
        }
    }
}

void solver_t::advanceEz(const std::vector<std::vector<std::vector<integersizeofmediamatrices>>>& sggMiEz) {
    double Idyhj;
    int i, j, k;
    integersizeofmediamatrices medio;

    for (k = 0; k < this->bounds.sweepEz.NZ; ++k) {
        for (j = 0; j < this->bounds.sweepEz.NY; ++j) {
            for (i = 0; i < this->bounds.sweepEz.NX; ++i) {
                Idyhj = this->Idyh[j];
                medio = sggMiEz[i][j][k];
                
                this->Ez[i][j][k] = this->g.g1(medio) * this->Ez[i][j][k] + 
                                    this->g.g2(medio) * 
                                    ((this->Hy[i][j][k] - this->Hy[i-1][j][k]) * this->Idxh[i] - 
                                     (this->Hx[i][j][k] - this->Hx[i][j-1][k]) * Idyhj);
            }
        }
    }
}

void solver_t::advanceH() {
#ifdef CompileWithProfiling
    nvtxStartRange("Antes del bucle HX");
#endif
    advanceHx(this->media.sggMiHx);
#ifdef CompileWithProfiling
    nvtxEndRange();
    nvtxStartRange("Antes del bucle HY");
#endif
    advanceHy(this->media.sggMiHy);
#ifdef CompileWithProfiling
    nvtxEndRange();
    nvtxStartRange("Antes del bucle HZ");
#endif
    advanceHz(this->media.sggMiHz);
#ifdef CompileWithProfiling
    nvtxEndRange();
#endif
}

void solver_t::advanceHx(const std::vector<std::vector<std::vector<integersizeofmediamatrices>>>& sggMiHx) {
    double Idzek, Idyej;
    int i, j, k;
    integersizeofmediamatrices medio;

    for (k = 0; k < this->bounds.sweepHx.NZ; ++k) {
        for (j = 0; j < this->bounds.sweepHx.NY; ++j) {
            for (i = 0; i < this->bounds.sweepHx.NX; ++i) {
                Idzek = this->IdzE[k];
                Idyej = this->IdyE[j];
                medio = sggMiHx[i][j][k];
                
                this->Hx[i][j][k] = this->g.gm1(medio) * this->Hx[i][j][k] + 
                                    this->g.gm2(medio) * 
                                    ((this->Ey[i][j][k+1] - this->Ey[i][j][k]) * Idzek - 
                                     (this->Ez[i][j+1][k] - this->Ez[i][j][k]) * Idyej);
            }
        }
    }
}

void solver_t::advanceHy(const std::vector<std::vector<std::vector<integersizeofmediamatrices>>>& sggMiHy) {
    double Idzek;
    int i, j, k;
    integersizeofmediamatrices medio;

    for (k = 0; k < this->bounds.sweepHy.NZ; ++k) {
        for (j = 0; j < this->bounds.sweepHy.NY; ++j) {
            for (i = 0; i < this->bounds.sweepHy.NX; ++i) {
                Idzek = this->IdzE[k];
                medio = sggMiHy[i][j][k];
                
                this->Hy[i][j][k] = this->g.gm1(medio) * this->Hy[i][j][k] + 
                                    this->g.gm2(medio) * 
                                    ((this->Ez[i+1][j][k] - this->Ez[i][j][k]) * this->IdxE[i] - 
                                     (this->Ex[i][j][k+1] - this->Ex[i][j][k]) * Idzek);
            }
        }
    }
}

void solver_t::advanceHz(const std::vector<std::vector<std::vector<integersizeofmediamatrices>>>& sggMiHz) {
    double Idyej;
    int i, j, k;
    integersizeofmediamatrices medio;

    for (k = 0; k < this->bounds.sweepHz.NZ; ++k) {
        for (j = 0; j < this->bounds.sweepHz.NY; ++j) {
            for (i = 0; i < this->bounds.sweepHz.NX; ++i) {
                Idyej = this->IdyE[j];
                medio = sggMiHz[i][j][k];
                
                this->Hz[i][j][k] = this->g.gm1(medio) * this->Hz[i][j][k] + 
                                    this->g.gm2(medio) * 
                                    ((this->Ex[i][j+1][k] - this->Ex[i][j][k]) * Idyej - 
                                     (this->Ey[i+1][j][k] - this->Ey[i][j][k]) * this->IdxE[i]);
            }
        }
    }
}

void solver_t::solver_advanceEDispersiveE() {
    if (this->thereAre.Edispersives) {
        AdvanceEDispersiveE(this->sgg);
    }
}

void solver_t::solver_advanceMDispersiveH() {
    if (this->thereAre.Mdispersives) {
        AdvanceMDispersiveH(this->sgg);
    }
}

void solver_t::solver_advanceLumpedE() {
    if (this->thereAre.Lumpeds) {
        AdvanceLumpedE(this->sgg, this->n, this->control.simu_devia, this->control.stochastic);
    }
}

void solver_t::solver_advancePlaneWaveE() {
    if (this->thereAre.PlaneWaveBoxes && this->still_planewave_time) {
        if (!this->control.simu_devia) {
            AdvancePlaneWaveE(this->sgg, this->n, this->bounds, this->g.G2,
                              this->Idxh, this->Idyh, this->Idzh,
                              this->Ex, this->Ey, this->Ez,
                              this->still_planewave_time);
        }
    }
}

void solver_t::solver_advancePlaneWaveH() {
    if (this->thereAre.PlaneWaveBoxes && this->still_planewave_time) {
        if (!this->control.simu_devia) {
            AdvancePlaneWaveH(this->sgg, this->n, this->bounds, this->g.GM2,
                              this->Idxe, this->Idye, this->Idze,
                              this->Hx, this->Hy, this->Hz,
                              this->still_planewave_time);
        }
    }
}

void solver_t::solver_advanceNodalE() {
    if (this->thereAre.NodalE) {
        advanceNodalE(this->sgg, this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz,
                      this->sgg.NumMedia, this->n, this->bounds, this->g.G2,
                      this->Idxh, this->Idyh, this->Idzh,
                      this->Ex, this->Ey, this->Ez,
                      this->control.simu_devia);
    }
}

void solver_t::solver_advanceNodalH() {
    if (this->thereAre.NodalH) {
        AdvanceNodalH(this->sgg, this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz,
                      this->sgg.NumMedia, this->n, this->bounds, this->g.GM2,
                      this->Idxe, this->Idye, this->Idze,
                      this->Hx, this->Hy, this->Hz,
                      this->control.simu_devia);
    }
}

void solver_t::solver_advanceAnisotropicE() {
    if (this->thereAre.Anisotropic) {
        AdvanceAnisotropicE(this->sgg.alloc, this->ex, this->ey, this->ez,
                            this->hx, this->hy, this->hz,
                            this->Idxe, this->Idye, this->Idze,
                            this->Idxh, this->Idyh, this->Idzh);
    }
}

void solver_t::solver_advanceAnisotropicH() {
    if (this->thereAre.Anisotropic) {
        AdvanceAnisotropicH(this->sgg.alloc, this->ex, this->ey, this->ez,
                            this->hx, this->hy, this->hz,
                            this->Idxe, this->Idye, this->Idze,
                            this->Idxh, this->Idyh, this->Idzh);
    }
}

void solver_t::solver_advancePMLbodyH() {
    if (this->thereAre.PMLbodies) {
        AdvancePMLbodyH();
    }
}

void solver_t::solver_advanceMagneticCPML() {
    if (this->thereAre.PMLBorders) {
        advanceMagneticCPML(this->sgg.numMedia, this->bounds,
                            this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz,
                            this->g.gm2, this->Hx, this->Hy, this->Hz,
                            this->Ex, this->Ey, this->Ez);
    }
}

void solver_t::solver_MinusCloneMagneticPMC() {
    if (this->thereAre.PMCBorders) {
        MinusCloneMagneticPMC(this->sgg.alloc, this->sgg.border, this->Hx, this->Hy, this->Hz, this->sgg.sweep,
                              this->control.layoutnumber, this->control.num_procs);
    }
}

void solver_t::solver_CloneMagneticPeriodic() {
    if (this->thereAre.PeriodicBorders) {
        CloneMagneticPeriodic(this->sgg.alloc, this->sgg.border, this->Hx, this->Hy, this->Hz, this->sgg.sweep,
                              this->control.layoutnumber, this->control.num_procs);
    }
}

void solver_t::solver_advancePMLE() {
    if (this->thereAre.PMLbodies) { // waveport absorbers
        AdvancePMLbodyE();
    }
    if (this->thereAre.PMLBorders) {
        AdvanceelectricCPML(this->sgg.numMedia, this->bounds, this->media.sggMiEx, this->media.sggMiEy, this->media.sggMiEz,
                            this->g.G2, this->Ex, this->Ey, this->Ez, this->Hx, this->Hy, this->Hz);
    }
}

void solver_t::solver_advancesgbcE() {
    if (this->thereAre.sgbcs && (this->control.sgbc)) {
        AdvancesgbcE(static_cast<double>(this->sgg.dt), this->control.sgbcDispersive,
                     this->control.simu_devia, this->control.stochastic);
    }
}

void solver_t::solver_advancesgbcH() {
    if (this->thereAre.sgbcs && (this->control.sgbc)) {
        AdvancesgbcH();
    }
}

void solver_t::solver_advanceWiresE() {
#ifdef CompileWithMTLN
    AdvanceWiresE_mtln(this->sgg, this->Idxh, this->Idyh, this->Idzh, this->eps0, this->mu0);
#else
    if ((trim(adjustl(this->control.wiresflavor)) == "holland") ||
        (trim(adjustl(this->control.wiresflavor)) == "transition")) {
        if (this->thereAre.Wires) {
            if (this->control.wirecrank) {
                AdvanceWiresEcrank(this->sgg, this->n, this->control.layoutnumber, this->control.wiresflavor, this->control.simu_devia, this->control.stochastic);
            } else {
                AdvanceWiresE(this->sgg, this->n, this->control.layoutnumber, this->control.wiresflavor, this->control.simu_devia, this->control.stochastic, this->control.experimentalVideal, this->control.wirethickness, this->eps0, this->mu0);
            }
        }
    }
#ifdef CompileWithBerengerWires
    if (trim(adjustl(this->control.wiresflavor)) == "berenger") {
        if (this->thereAre.Wires) {
            AdvanceWiresE_Berenger(this->sgg, this->n);
        }
    }
#endif
#ifdef CompileWithSlantedWires
    if ((trim(adjustl(this->control.wiresflavor)) == "slanted") || (trim(adjustl(this->control.wiresflavor)) == "semistructured")) {
        AdvanceWiresE_Slanted(this->sgg, this->n);
    }
#endif
#endif
}

void solver_t::solver_advancewiresH() {
    if ((trim(adjustl(this->control.wiresflavor)) == "holland") ||
        (trim(adjustl(this->control.wiresflavor)) == "transition")) {
        if (this->thereAre.Wires) {
            if (this->control.wirecrank) {
                // continue
            } else {
                AdvanceWiresH(this->sgg, this->n, this->control.layoutnumber,
                              this->control.wiresflavor, this->control.simu_devia, this->control.stochastic,
                              this->control.experimentalVideal, this->control.wirethickness, this->eps0, this->mu0);
            }
        }
    }
}

void solver_t::solver_advanceMagneticMUR() {
#ifdef CompileWithMPI
    int ierr = 0;
#endif
    if (this->thereAre.MURBorders) {
        AdvanceMagneticMUR(this->bounds, this->sgg,
                           this->media.sggMiHx, this->media.sggMiHy, this->media.sggMiHz,
                           this->Hx, this->Hy, this->Hz,
                           this->control.mur_second);
#ifdef CompileWithMPI
        if (this->control.mur_second) {
            if (this->control.num_procs > 1) {
                MPI_Barrier(SUBCOMM_MPI, ierr);
                FlushMPI_H_Cray();
            }
        }
#endif
    }
}

void solver_t::solver_end() {
    // Ex, Ey, Ez, Hx, Hy, Hz, dxe, dye, dze, dxh, dyh, dzh are pointers in Fortran
    // In C++, we just use the members directly.
    
    real_kind_tiempo at; // Unused in the visible code snippet end
    int ndummy = 0;
    bool dummylog = false;
    bool somethingdone = false; // Unused
    bool newsomethingdone = false; // Unused
    char dubuf[bufsize];

#ifdef CompileWithMPI
    int ierr = 0;
#endif

#ifdef CompileWithProfiling
    nvtxEndRange();
#endif

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

    if (this->n > this->control.finaltimestep) {
        this->n = this->control.finaltimestep;
    }
    this->control.finaltimestep = this->n;
    this->lastexecutedtime = this->sgg.tiempo(this->control.finaltimestep);

    Timing(this->sgg, this->bounds, this->n, ndummy, this->control.layoutnumber, this->control.num_procs,
           this->control.maxCPUtime, this->control.flushsecondsFields, this->control.flushsecondsData,
           this->initialtimestep, this->control.finaltimestep, this->d_perform, dummylog, false,
           this->Ex, this->Ey, this->Ez, this->everflushed, this->control.nentradaroot, this->control.maxSourceValue,
           this->control.opcionestotales, this->control.simu_devia, this->control.dontwritevtk, this->control.permitscaling);

    snprintf(dubuf, bufsize, "END FDTD time stepping. Beginning posprocessing at n= %d", this->n);
    print11(this->control.layoutnumber, dubuf);

    if ((this->control.flushsecondsFields != 0) || this->perform.flushFIELDS) {
        snprintf(dubuf, bufsize, " INIT FINAL FLUSHING OF RESTARTING FIELDS n= %d", this->n);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, SEPARADOR);
        flush_and_save_resume(this->sgg, this->bounds, this->control.layoutnumber, this->control.num_procs, this->control.nentradaroot, this->control.nresumeable2, this->thereare, this->n, this->eps0, this->mu0, this->everflushed,
                              this->Ex, this->Ey, this->Ez, this->Hx, this->Hy, this->Hz, this->control.wiresflavor, this->control.simu_devia, this->control.stochastic);
        snprintf(dubuf, bufsize, " DONE FINAL FLUSHING OF RESTARTING FIELDS N=%d", this->n);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, dubuf);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, SEPARADOR);
        print11(this->control.layoutnumber, SEPARADOR);
    }

    if (this->thereAre.FarFields) {
        snprintf(dubuf, bufsize, " INIT FINAL OBSERVATION DATA FLUSHING and Near-to-Far field  n= %d", this->n);
    } else {
        snprintf(dubuf, bufsize, " INIT FINAL OBSERVATION DATA FLUSHING n= %d", this->n);
    }
    print11(this->control.layoutnumber, SEPARADOR);
    print11(this->control.layoutnumber, SEPARADOR);
    print11(this->control.layoutnumber, SEPARADOR);
    print11(this->control.layoutnumber, dubuf);
    print11(this->control.layoutnumber, SEPARADOR);
    print11(this->control.layoutnumber, SEPARADOR);
    print11(this->control.layoutnumber, SEPARADOR);
    
    if (this->thereAre.Observation) {
        FlushObservationFiles(this->sgg, this->ini_save, this->n, this->control.layoutnumber, this->control.num_procs, this->dxe, this->dye, this->dze, this->dxh, this->dyh, this->dzh, this->bounds, this->control.singlefilewrite, this->control.facesNF2FF, true);
    }
}

CloseObservationFiles(this->sgg, this->control.layoutnumber, this->control.num_procs, this->control.singlefilewrite, this->initialtimestep, this->lastexecutedtime, this->control.resume); //dump the remaining to disk
        }
        
        if (this->thereAre.FarFields) {
            snprintf(dubuf, sizeof(dubuf), " DONE FINAL OBSERVATION DATA FLUSHED and Near-to-Far field  n= %9d", this->n);
        } else {
            snprintf(dubuf, sizeof(dubuf), " DONE FINAL OBSERVATION  DATA FLUSHED n= %9d", this->n);
        }
        print11(this->control.layoutnumber, std::string(SEPARADOR) + separador + separador);
        print11(this->control.layoutnumber, dubuf);
        print11(this->control.layoutnumber, std::string(SEPARADOR) + separador + separador);

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

        snprintf(dubuf, sizeof(dubuf), "INIT FINAL Postprocessing frequency domain probes, if any, at n= %9d", this->n);
        print11(this->control.layoutnumber, dubuf);
        snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
        print11(this->control.layoutnumber, dubuf);
        somethingdone = false;
        at = this->n * this->sgg.dt;
        if (this->thereAre.Observation) {
            PostProcess(this->control.layoutnumber, this->control.num_procs, this->sgg, this->control.nentradaroot, at, somethingdone, this->control.niapapostprocess, this->control.forceresampled);
        }
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
        somethingdone = newsomethingdone;
#endif

        if (somethingdone) {
            snprintf(dubuf, sizeof(dubuf), "DONE FINAL Postprocessing frequency domain probes.");
            print11(this->control.layoutnumber, dubuf);
            snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
            print11(this->control.layoutnumber, dubuf);
        } else {
            snprintf(dubuf, sizeof(dubuf), "No FINAL frequency domain probes snapshots found to be postrocessed");
            print11(this->control.layoutnumber, dubuf);
            snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
            print11(this->control.layoutnumber, dubuf);
        }

        snprintf(dubuf, sizeof(dubuf), "INIT FINAL FLUSHING .vtk if any.");
        print11(this->control.layoutnumber, dubuf);
        snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
        print11(this->control.layoutnumber, dubuf);
        somethingdone = false;

        if (this->thereAre.Observation) {
            createvtk(this->control.layoutnumber, this->control.num_procs, this->sgg, this->control.vtkindex, somethingdone, this->control.mpidir, this->media.sggMtag, this->control.dontwritevtk);
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
        somethingdone = newsomethingdone;
#endif
        if (somethingdone) {
            snprintf(dubuf, sizeof(dubuf), "DONE FINAL FLUSHING .vtk snapshots");
            print11(this->control.layoutnumber, dubuf);
            snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
            print11(this->control.layoutnumber, dubuf);
        } else {
            snprintf(dubuf, sizeof(dubuf), "No FINAL .vtk snapshots found to be flushed");
            print11(this->control.layoutnumber, dubuf);
            snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
            print11(this->control.layoutnumber, dubuf);
        }

        snprintf(dubuf, sizeof(dubuf), "INIT FINAL FLUSHING .xdmf if any.");
        print11(this->control.layoutnumber, dubuf);
        snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
        print11(this->control.layoutnumber, dubuf);
        somethingdone = false;
        if (this->thereAre.Observation) {
            createxdmf(this->sgg, this->control.layoutnumber, this->control.num_procs, this->control.vtkindex, this->control.createh5bin, somethingdone, this->control.mpidir);
        }
        if (this->control.createh5bin) {
            createh5bintxt(this->sgg, this->control.layoutnumber, this->control.num_procs); //lo deben llamar todos haya o no this%thereAre%observation
        }
        //        call create_interpreted_mesh(sgg)
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
        MPI_AllReduce(&somethingdone, &newsomethingdone, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, &ierr);
        somethingdone = newsomethingdone;
#endif
        if (somethingdone) {
            snprintf(dubuf, sizeof(dubuf), "DONE FINAL FLUSHING .xdmf snapshots");
            print11(this->control.layoutnumber, dubuf);
            snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
            print11(this->control.layoutnumber, dubuf);
        } else {
            snprintf(dubuf, sizeof(dubuf), "No FINAL .xdmf snapshots found to be flushed");
            print11(this->control.layoutnumber, dubuf);
            snprintf(dubuf, sizeof(dubuf), "%s%s%s", SEPARADOR, separador, separador);
            print11(this->control.layoutnumber, dubuf);
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
        Timing(this->sgg, this->bounds, this->n, ndummy, this->control.layoutnumber,
               this->control.num_procs, this->control.maxCPUtime, this->control.flushsecondsFields,
               this->control.flushsecondsData, this->initialtimestep,
               this->control.finaltimestep, this->perform, parar, .FALSE.,
               Ex, Ey, Ez, this->everflushed, this->control.nentradaroot, this->control.maxSourceValue, this->control.opcionestotales,
               this->control.simu_devia, this->control.dontwritevtk, this->control.permitscaling);
        snprintf(dubuf, sizeof(dubuf), "END FINAL POSTPROCESSING at n= %9d", this->n);
        print11(this->control.layoutnumber, dubuf);
        this->finishedwithsuccess = true;
        return;

    } // end subroutine

    //las sggmixx se desctruyen el en main pq se alocatean alli
    void Destroy_All_exceptSGGMxx(SGGFDTDINFO_t& sgg, std::vector<std::vector<std::vector<double>>>& Ex, std::vector<std::vector<std::vector<double>>>& Ey, std::vector<std::vector<std::vector<double>>>& Ez,
                                  std::vector<std::vector<std::vector<double>>>& Hx, std::vector<std::vector<std::vector<double>>>& Hy, std::vector<std::vector<std::vector<double>>>& Hz,
                                  std::vector<double>& G1, std::vector<double>& G2, std::vector<double>& GM1, std::vector<double>& GM2,
                                  std::vector<double>& dxe, std::vector<double>& dye, std::vector<double>& dze, std::vector<double>& Idxe, std::vector<double>& Idye, std::vector<double>& Idze,
                                  std::vector<double>& dxh, std::vector<double>& dyh, std::vector<double>& dzh, std::vector<double>& Idxh, std::vector<double>& Idyh, std::vector<double>& Idzh,
                                  const logic_control_t& thereare, const std::string& wiresflavor) {
        DestroyObservation(sgg);
        DestroyNodal(sgg);
        DestroyIlumina(sgg);
#ifdef CompileWithNIBC
        DestroyMultiports(sgg);
#endif

        destroysgbcs(sgg); //todos deben destruir pq alocatean en funcion de sgg no de si contienen estos materiales que lo controla therearesgbcs. Lo que habia era if ((this%thereAre%sgbcs).and.(sgbc))
        destroyLumped(sgg);
        DestroyEDispersives(sgg);
        DestroyMDispersives(sgg);
        if ((wiresflavor == "holland") || (wiresflavor == "transition")) {
            DestroyWires(sgg);
        }
#ifdef CompileWithBerengerWires
        if (wiresflavor == "berenger") {
            DestroyWires_Berenger(sgg);
        }
#endif
#ifdef CompileWithSlantedWires
        if ((wiresflavor == "slanted") || (wiresflavor == "semistructured")) {
            DestroyWires_Slanted(sgg);
        }
#endif

        DestroyCPMLBorders();
        DestroyPMLbodies(sgg);
        DestroyMURBorders();
        //Destroy the remaining
        sgg.Med.clear();
        sgg.LineX.clear();
        sgg.LineY.clear();
        sgg.LineZ.clear();
        sgg.DX.clear();
        sgg.DY.clear();
        sgg.DZ.clear();
        sgg.tiempo.clear();
        G1.clear();
        G2.clear();
        GM1.clear();
        GM2.clear();
        Ex.clear();
        Ey.clear();
        Ez.clear();
        Hx.clear();
        Hy.clear();
        Hz.clear();
        dxe.clear();
        dye.clear();
        dze.clear();
        Idxe.clear();
        Idye.clear();
        Idze.clear();
        dxh.clear();
        dyh.clear();
        dzh.clear();
        Idxh.clear();
        Idyh.clear();
        Idzh.clear();
    }

    void destroy_and_deallocate(solver_t& this_obj) {
        DestroyObservation(this_obj.sgg);
        DestroyNodal(this_obj.sgg);
        DestroyIlumina(this_obj.sgg);
#ifdef CompileWithNIBC
        DestroyMultiports(this_obj.sgg);
#endif

        destroysgbcs(this_obj.sgg); //todos deben destruir pq alocatean en funcion de this%sgg no de si contienen estos materiales que lo controla therearesgbcs. Lo que habia era if ((this%thereAre%sgbcs).and.(sgbc))
        destroyLumped(this_obj.sgg);
        DestroyEDispersives(this_obj.sgg);
        DestroyMDispersives(this_obj.sgg);
        if ((this_obj.control->wiresflavor == "holland") || (this_obj.control->wiresflavor == "transition")) {
            DestroyWires(this_obj.sgg);
        }
#ifdef CompileWithBerengerWires
        if (this_obj.control->wiresflavor == "berenger") {
            DestroyWires_Berenger(this_obj.sgg);
        }
#endif
#ifdef CompileWithSlantedWires
        if ((this_obj.control->wiresflavor == "slanted") || (this_obj.control->wiresflavor == "semistructured")) {
            DestroyWires_Slanted(this_obj.sgg);
        }
#endif

        DestroyCPMLBorders();
        DestroyPMLbodies(this_obj.sgg);
        DestroyMURBorders();
        //Destroy the remaining
        this_obj.sgg->Med.clear();
        this_obj.sgg->LineX.clear();
        this_obj.sgg->LineY.clear();
        this_obj.sgg->LineZ.clear();
        this_obj.sgg->DX.clear();
        this_obj.sgg->DY.clear();
        this_obj.sgg->DZ.clear();
        this_obj.sgg->tiempo.clear();
        this_obj.g->destroy();
        this_obj.Ex.clear();
        this_obj.Ey.clear();
        this_obj.Ez.clear();
        this_obj.Hx.clear();
        this_obj.Hy.clear();
        this_obj.Hz.clear();
        this_obj.dxe.clear();
        this_obj.dye.clear();
        this_obj.dze.clear();
        this_obj.Idxe.clear();
        this_obj.Idye.clear();
        this_obj.Idze.clear();
        this_obj.dxh.clear();
        this_obj.dyh.clear();
        this_obj.dzh.clear();
        this_obj.Idxh.clear();
        this_obj.Idyh.clear();
        this_obj.Idzh.clear();
    }

} // end module

