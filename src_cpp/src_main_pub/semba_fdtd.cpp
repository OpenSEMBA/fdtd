#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <optional>
#include <sstream>

// Includes for external modules/types referenced in the Fortran code
// These would need to be implemented or included from their respective headers
// #include "version_m.h"
// #include "Report_m.h"
// #include "Getargs_m.h"
// #include "FDETYPES_m.h"
// #include "Solver_m.h"
// #include "resuming_m.h"
// #include "NFDETypes_m.h"
// #include "nfde_rotate_m.h"
// #include "Preprocess_m.h"
// #include "storeData_m.h"
// #include "xdmf_h5_m.h"
// #include "EpsMuTimeScale_m.h"
// #include "interpreta_switches_m.h"

#ifdef CompileWithMPI
#include <mpi.h>
// #include "MPIcomm_m.h"
// #include "build_t_linea_mpi_m.h"
#endif

#ifdef CompileWithMTLN
// #include "mtln_t.h"
#endif

#ifdef CompileWithSMBJSON
// #include "smbjson_m.h"
#endif

#ifdef CompilePrivateVersion
// #include "ParseadorClass.h"
#endif

// Forward declarations for types assumed to be in included modules
struct entrada_t;
struct tiempo_t;
struct media_matrices_t;
struct SGGFDTDINFO_t;
struct limit_t;
struct taglist_t;
struct tagtype_t;
struct Parseador_t;
struct t_NFDE_FILE_t;
struct solver_t;
struct XYZlimit_t;

// Constants and Types assumed from modules
#ifndef RKIND
#define RKIND double
#endif

#ifndef BUFSIZE
#define BUFSIZE 256
#endif

#ifndef BUFSIZE_LONG
#define BUFSIZE_LONG 1024
#endif

#ifndef SEPARADOR
#define SEPARADOR "================================================================"
#endif

// Placeholder for global variables that might exist in the Fortran modules
// In a real translation, these would be part of the respective namespace/class implementations
extern int SUBCOMM_MPI; // Assumed global from MPIcomm_m

namespace SEMBA_FDTD_m {

    // Derived type: semba_fdtd_t
    struct semba_fdtd_t {
        entrada_t l;
        tiempo_t time_comienzo;
        double time_desdelanzamiento; // real(kind=8)
        media_matrices_t media;
        SGGFDTDINFO_t sgg;
        std::vector<limit_t> fullsize; // dimension(1:6)
        std::vector<limit_t> SINPML_fullsize; // dimension(1:6)
        RKIND eps0;
        RKIND mu0;
        RKIND cluz;
        RKIND maxSourceValue;
        std::string whoami; // character(len=BUFSIZE)
        std::string whoamishort; // character(len=BUFSIZE)
#ifdef CompileWithMTLN
        mtln_t mtln_parsed;
#endif
        taglist_t tag_numbers;
        tagtype_t tagtype;
        bool finishedwithsuccess; // logical

        // Methods
        void init(const std::optional<std::string>& input_flags = std::nullopt);
        void launch();
        void end();
        void create_solver();
        void update_after_simulation();
    };

    // Helper functions that would be in other modules
    // These are declared here to allow compilation of the method bodies if the actual implementations are elsewhere
    void initEntrada(entrada_t& l);
    void OnPrint();
    void setglobal(int layoutnumber, int num_procs);
    void get_secnds(tiempo_t& time_out2);
    void print_credits(const entrada_t& l);
    void CLOSEWARNINGFILE(int layoutnumber, int num_procs, bool& dummylog, bool l_auxinput, bool l_auxoutput);
    void default_flags(entrada_t& l);
    void getcommandargument(const std::string& chain2, int arg_index, std::string& chaindummy, int& length, int& status, const std::string& binary_path);
    void busca_switchficheroinput(entrada_t& l);
    void data_loader(const std::string& filefde, Parseador_t*& parser);
    void interpreta(entrada_t& l, int& status);
    void nfde_rotate(Parseador_t* parser, const std::string& mpidir);
    void set_priorities(bool prioritizeCOMPOoverPEC, bool prioritizeISOTROPICBODYoverall, bool prioritizeTHINWIRE);
    void print11(int layoutnumber, const std::string& msg);
    void erasesignalingfiles(const std::string& simu_devia);
    void stoponerror(int layoutnumber, int num_procs, const std::string& msg, bool fatal = false);
    void AssigLossyOrPECtoNodes(SGGFDTDINFO_t& sgg, media_matrices_t& media);
    void store_geomData(const SGGFDTDINFO_t& sgg, const media_matrices_t& media, const std::string& geomfile);
    void WarnErrReport(const std::string& msg, bool fatal);
    std::string getBinaryPath(); // Assumed function

    // Implementation of semba_fdtd_t methods

    void semba_fdtd_t::init(const std::optional<std::string>& input_flags) {
        double dtantesdecorregir = 0.0;
        double dxmin, dymin, dzmin, dtlay;
        
        bool dummylog, l_auxinput, l_auxoutput, ThereArethinslots = false;
        bool hayinput = false;
        bool lexis = false;

        std::string f = " ";
        std::string chain = " ";
        std::string chain3 = " ";
        std::string chain4 = " ";
        std::string chaindummy = " ";
        std::string slices = " ";
        std::string dubuf;
        std::string buff;
        std::string filename_h5bin;

        int myunit = 0;
        int jmed = 0;
        int finaltimestepantesdecorregir = 0;
        int NEWfinaltimestep = 0;
        int thefileno = 0;
        int statuse = 0;
        int status = 0;
        int i = 0;
        int field = 0;
        int my_iostat = 0;

        Parseador_t* parser = nullptr;
        t_NFDE_FILE_t* NFDE_FILE = nullptr;
        solver_t solver;

#ifdef CompileWithMPI
        bool fatalerror_aux = false;
        std::vector<XYZlimit_t> tempalloc(6);
#endif

        int conf_err = 0;

        initEntrada(this->l);

        this->eps0 = 8.8541878176203898505365630317107502606083701665994498081024171524053950954599821142852891607182008932e-12;
        this->mu0 = 1.2566370614359172953850573533118011536788677597500423283899778369231265625144835994512139301368468271e-6;
        this->cluz = 1.0 / std::sqrt(this->eps0 * this->mu0);
        
        OnPrint();

#ifdef CompileWithMPI
        // InitGeneralMPI(this->l.layoutnumber, this->l.num_procs);
        SUBCOMM_MPI = MPI_COMM_WORLD;
#else
        this->l.num_procs = 1;
        this->l.layoutnumber = 0;
#endif
        setglobal(this->l.layoutnumber, this->l.num_procs);
        
        this->whoamishort = std::to_string(this->l.layoutnumber + 1);
        this->whoami = "(" + std::to_string(this->l.layoutnumber + 1) + "/" + std::to_string(this->l.num_procs) + ") ";
        
#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        // get_secnds(this->l.time_out2);
        // this->time_desdelanzamiento = this->l.time_out2.segundos;

#ifndef keeppause
        if (this->l.layoutnumber == 0) {
            // File handling for running, pause, relaunch, forcestop
            // Using fstream for C++
            std::ofstream ofs;
            ofs.open("running");
            ofs << "!END";
            ofs.close();
            std::remove("running");

            ofs.open("pause");
            ofs << "!END";
            ofs.close();
            std::remove("pause");

            ofs.open("relaunch");
            ofs << "!END";
            ofs.close();
            std::remove("relaunch");

            ofs.open("forcestop");
            ofs << "!END";
            ofs.close();
            std::remove("forcestop");
        }
#endif

        if (this->l.layoutnumber == 0) {
            my_iostat = 0;
            // Retry loop for opening temp log
            while (true) {
                std::ifstream ifs("SEMBA_FDTD_temp.log");
                if (!ifs) {
                    my_iostat = 0;
                    break;
                }
                ifs.close();
                std::cout << "." << std::flush;
                std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Simple wait
            }
            
            std::ofstream ofs_log("SEMBA_FDTD_temp.log", std::ios::out | std::ios::trunc);
            ofs_log << "!END";
            ofs_log.close();
            std::remove("SEMBA_FDTD_temp.log");
            my_iostat = 0;

            while (true) {
                std::ifstream ifs("SEMBA_FDTD_temp.log");
                if (!ifs) {
                    my_iostat = 0;
                    break;
                }
                ifs.close();
                std::cout << "." << std::flush;
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
            }

            std::ofstream ofs_new("SEMBA_FDTD_temp.log", std::ios::out | std::ios::trunc);
            print_credits(this->l);
            ofs_new.close();
        }

#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif

        // Label 652 continue
        // CLOSEWARNINGFILE(this->l.layoutnumber, this->l.num_procs, dummylog, false, false);
        
        // this->l.opcionespararesumeo = "mpirun -n " + std::to_string(this->l.num_procs) + " ";
        // default_flags(this->l);

#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        // get_secnds(this->time_comienzo);

        if (this->l.layoutnumber == 0) {
            std::ofstream ofs_append("SEMBA_FDTD_temp.log", std::ios::app);
            this->l.file11isopen = true;
            ofs_append.close();
        }

#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif

        // Check pause semaphore
        // inquire(file='pause', EXIST=this->l.pausar);
        this->l.pausar = false; // Placeholder
#ifdef CompileWithMPI
        this->l.l_aux = this->l.pausar;
        // MPI_AllReduce(&this->l.l_aux, &this->l.pausar, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this->l.ierr);
#endif
#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        // get_secnds(this->l.time_out2);
        // this->l.time_begin = this->l.time_out2.segundos;
        
        dubuf = "Paused at              "; // Simplified formatting
        if (this->l.pausar) print11(this->l.layoutnumber, dubuf);
        
        while (this->l.pausar) {
#ifdef CompileWithMPI
            // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
            // get_secnds(this->l.time_out2);
            // this->l.time_end = this->l.time_out2.segundos;
            
            // if (this->l.time_end - this->l.time_begin > 10.0) {
            //     inquire(file='pause', EXIST=this->l.pausar);
#ifdef CompileWithMPI
            //     MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            //     this->l.l_aux = this->l.pausar;
            //     MPI_AllReduce(&this->l.l_aux, &this->l.pausar, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, this->l.ierr);
            //     MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
            //     get_secnds(this->l.time_out2);
            //     this->l.time_begin = this->l.time_out2.segundos;
            //     dubuf = "Paused at              ";
            //     if (this->l.pausar) print11(this->l.layoutnumber, dubuf);
            // }
            break; // Placeholder to avoid infinite loop in translation
        }

#ifdef keeppause
        // inquire(file='forcestop', EXIST=this->l.forcestop);
        this->l.forcestop = false; // Placeholder
        if (this->l.forcestop) {
            if (this->l.layoutnumber == 0) {
                std::remove("running");
                std::remove("pause");
                std::remove("relaunch");
                std::remove("forcestop");
            }
#ifdef CompileWithMPI
            // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            // MPI_Finalize(this->l.ierr);
#endif
            return; // STOP
        }
#endif

#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        // get_secnds(this->l.time_out2);

        if (input_flags.has_value()) {
            this->l.read_command_line = false;
            this->l.chain2 = input_flags.value();
            this->l.length = this->l.chain2.length();
        } else {
            // get_command(this->l.chain2, this->l.length, status);
            this->l.chain2 = "";
            this->l.length = 0;
            status = 0;
            if (status != 0) {
                stoponerror(this->l.layoutnumber, this->l.num_procs, "General error", true);
                goto label_652;
            }
        }

        this->l.chain2 = std::string(this->l.chain2.begin(), std::find_if(this->l.chain2.begin(), this->l.chain2.end(), [](unsigned char c){ return !std::isspace(c); }));
        
        // inquire(file='launch', EXIST=hayinput);
        hayinput = false; // Placeholder
        if (hayinput) {
            std::ifstream ifs_launch("launch");
            std::getline(ifs_launch, chain3);
            chain3 = std::string(chain3.begin(), std::find_if(chain3.begin(), chain3.end(), [](unsigned char c){ return !std::isspace(c); }));
            ifs_launch.close();
            std::cout << "----> launch input file " << chain3 << std::endl;
        }
#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif

        this->l.chain2 = this->l.chain2 + " " + chain3;

        // buscaswitchficheroinput(this->l);

label_652:
        if (status != 0) {
            stoponerror(this->l.layoutnumber, this->l.num_procs, "Error in searching input file. Correct and remove pause file", true);
            goto label_652;
        }

        print_credits(this->l);

#ifdef CompileWithMPI
        // initialize_MPI_process(this->l.filefde, this->l.extension);
#else
#ifdef CompilePrivateVersion
        if (this->l.extension == ".nfde") {
#ifndef CompileWithMTLN
            // NFDE_FILE = cargar_NFDE_FILE(this->l.filefde);
#else
            WarnErrReport(".nfde files are not supported when compiling with MTLN.", true);
#endif
        } else {
            NFDE_FILE = new t_NFDE_FILE_t();
        }
#else
        NFDE_FILE = new t_NFDE_FILE_t();
#endif
#endif

        // data_loader(this->l.filefde, parser);
        parser = nullptr; // Placeholder

        this->sgg.extraswitches = parser ? parser->switches : "";

        // getcommandargument(this->l.chain2, 1, chaindummy, this->l.length, statuse, getBinaryPath());
        chaindummy = "";
        this->l.chain2 = std::string(this->l.chain2.begin(), std::find_if(this->l.chain2.begin(), this->l.chain2.end(), [](unsigned char c){ return !std::isspace(c); }));
        chaindummy = std::string(chaindummy.begin(), std::find_if(chaindummy.begin(), chaindummy.end(), [](unsigned char c){ return !std::isspace(c); }));
        this->l.length = chaindummy.length();
        this->l.chain2 = chaindummy + " " + this->sgg.extraswitches + " " + this->l.chain2.substr(this->l.length);
        this->l.chaininput = this->l.chain2;

        // interpreta(this->l, status);
        this->sgg.nEntradaRoot = this->l.nEntradaRoot;

#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif

        // nfde_rotate(parser, this->l.mpidir);

#ifdef CompileWithMPI
        // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif

#ifdef CompileWithMTLN
        // if (parser->general.mtlnProblem) {
        //     solver.launch_mtln_simulation(parser->mtln, this->l.nEntradaRoot, this->l.layoutnumber);
        //     return;
        // }
#endif

#ifdef CompileWithHDF
        // if (this->l.createh5filefromsinglebin) {
        //     if (this->l.layoutnumber == 0) {
        //         // inquire(file=..., exist=lexis);
        //         lexis = false;
        //         if (!lexis) goto label_9083;
        //         std::ifstream ifs_h5bin("test_h5bin.txt"); // Placeholder filename
        //         std::string fname;
        //         while (std::getline(ifs_h5bin, fname)) {
        //             // createh5filefromsinglebin(fname, this->l.vtkindex);
        //             std::cout << "Processed " << fname << std::endl;
        //         }
        //         ifs_h5bin.close();
        //         std::cout << "END: SUCCESS creating " << this->sgg.nEntradaRoot << "_h5bin.txt" << std::endl;
        //         return;
        //     }
        //     // MPI_Barrier...
        //     return;
        // label_9083:
        //     stoponerror(0, this->l.num_procs, "Invalid _h5bin.txt file", true);
        //     statuse = -1;
        // }
#endif

        if (status != 0) {
            print11(this->l.layoutnumber, "Remove running and pause files. If error persists check switches for error.  " + this->l.chain2);
            print11(this->l.layoutnumber, " ");
            print11(this->l.layoutnumber, " ");
            print11(this->l.layoutnumber, " ");
            print11(this->l.layoutnumber, " ");
            print11(this->l.layoutnumber, " ");
            print11(this->l.layoutnumber, " ");
            goto label_652;
        }

        // set_priorities(this->l.prioritizeCOMPOoverPEC, this->l.prioritizeISOTROPICBODYoverall, this->l.prioritizeTHINWIRE);
        
        if (this->l.finaltimestep != -2) {
            print11(this->l.layoutnumber, "INIT conversion internal ASCII => Binary");
            print11(this->l.layoutnumber, SEPARADOR SEPARADOR SEPARADOR);
            
            // NFDE2sgg();
            this->l.fatalerror = this->l.fatalerror || this->l.fatalerrornfde2sgg;

#ifdef CompileWithMPI
            // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
            print11(this->l.layoutnumber, "[OK] Ended conversion internal ASCII => Binary");

            if (this->l.fatalerror) {
                // if (allocated(this->media.sggMiEx)) deallocate...
                stoponerror(this->l.layoutnumber, this->l.num_procs, "Error in .nfde file syntax. Check all *Warnings* and *tmpWarnings* files, correct and remove pause file if any", true);
                goto label_652;
            }

            // if (allocated(this->media.sggMiEx)) {
            AssigLossyOrPECtoNodes(this->sgg, this->media);

            if (this->l.createmap) store_geomData(this->sgg, this->media, this->l.geomfile);
            // }

#ifdef CompileWithMPI
            // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        }
        
        dubuf = "[OK] Ended Conformal Mesh";
        print11(this->l.layoutnumber, dubuf);
        
        if (this->l.finaltimestep == 0) this->l.finaltimestep = this->sgg.TimeSteps;
        
        if (this->l.forcesteps) {
            this->sgg.TimeSteps = this->l.finaltimestep;
#ifdef CompileWithMTLN
            this->mtln_parsed.number_of_steps = this->l.finaltimestep;
#endif
        } else {
            this->l.finaltimestep = this->sgg.TimeSteps;
        }

        if (!this->l.forcesteps) {
            finaltimestepantesdecorregir = this->l.finaltimestep;
            if (dtantesdecorregir != 0.0) {
                this->l.finaltimestep = static_cast<int>(dtantesdecorregir / this->sgg.dt * finaltimestepantesdecorregir);
            }
#ifdef CompileWithMPI
            // MPI_AllReduce(&this->l.finaltimestep, &NEWfinaltimestep, 1, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, this->l.ierr);
            // MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            // this->l.finaltimestep = NEWfinaltimestep;
#endif
#ifdef CompileWithMTLN
            this->mtln_parsed.number_of_steps = this->l.finaltimestep;
#endif
            if (finaltimestepantesdecorregir != this->l.finaltimestep) {
                dubuf = SEPARADOR SEPARADOR SEPARADOR;
                print11(this->l.layoutnumber, dubuf);
                dubuf = "Original Final Time Step= " + std::to_string(finaltimestepantesdecorregir);
                if (this->l.layoutnumber == 0) print11(this->l.layoutnumber, dubuf);
                dubuf = "Corrected Final Time Step= " + std::to_string(this->l.finaltimestep);
                if (this->l.layoutnumber == 0) print11(this->l.layoutnumber, dubuf);
            }
        }

        for (i = 1; i <= this->sgg.nummedia; ++i) {
            if (this->sgg.Med[i].Is.ThinWire) {
#ifndef CompileWithBerengerWires
                if (this->l.wiresflavor == "berenger") {
                    stoponerror(this->l.layoutnumber, this->l.num_procs, "Berenger Wires without support. Recompile!");
                }
#endif
#ifndef CompileWithSlantedWires
                if (this->l.wiresflavor == "slanted" || this->l.wiresflavor == "semistructured") {
                    stoponerror(this->l.layoutnumber, this->l.num_procs, "slanted Wires without support. Recompile!");
                }
#endif
                continue;
            }
            if (this->sgg.Med[i].Is.AnisMultiport || this->sgg.Med[i].Is.multiport || this->sgg.Med[i].Is.SGBC) {
#ifndef CompileWithNIBC
                if (this->l.mibc) stoponerror(this->l.layoutnumber, this->l.num_procs, "this->l.mibc Multiports without support. Recompile!");
#endif
                continue;
            }
#ifdef NoConformalSGBC
            if (this->sgg.Med[i].Is.sgbc && this->l.input_conformal_flag) {
                stoponerror(this->l.layoutnumber, this->l.num_procs, "Conformal sgbc not allowed. ");
            }
#endif
        }

        if (this->l.thereare_stoch && !this->l.chosenyesornostochastic) {
            stoponerror(this->l.layoutnumber, this->l.num_procs, "!STOCH found in .nfde. Specify either -stoch or -nostoch");
        }
#ifndef CompileWithSlantedWires
        if (this->l.hay_slanted_wires) {
            stoponerror(this->l.layoutnumber, this->l.num_procs, "slanted wires without slanted support. Recompile ()");
        }
#endif
        if (this->l.hay_slanted_wires && (this->l.wiresflavor != "slanted" && this->l.wiresflavor != "semistructured")) {
            stoponerror(this->l.layoutnumber, this->l.num_procs, "slanted wires require -this->l.wiresflavor Slanted/semistructured");
        }

        ThereArethinslots = false;
        for (jmed = 1; jmed <= this->sgg.NumMedia; ++jmed) {
            if (this->sgg.Med[jmed].Is.ThinSlot) ThereArethinslots = true;
        }
        if (this->l.resume && this->l.run_with_abrezanjas && ThereArethinslots) {
            stoponerror(this->l.layoutnumber, this->l.num_procs, "this->l.resume -r currently unsupported by conformal solver", true);
            statuse = -1;
        }

        if (this->l.layoutnumber == 0) {
            dubuf = SEPARADOR " " SEPARADOR " " SEPARADOR;
            print11(this->l.layoutnumber, dubuf);
            print11(this->l.layoutnumber, "Solver launched with options:");
            
            dubuf = this->l.mibc;
            print11(this->l.layoutnumber, "---> this->l.mibc    solver for NIBC multilayer: " + dubuf);
            
            dubuf = this->l.ade;
            print11(this->l.layoutnumber, "---> this->l.ade     solver for ADC multilayer: " + dubuf);
            
            dubuf = this->l.sgbc;
            print11(this->l.layoutnumber, "---> sgbc    solver for multilayer: " + dubuf);
            
            if (this->l.sgbc) {
                dubuf = this->l.sgbcDispersive;
                print11(this->l.layoutnumber, "---> sgbc DISPERSIVE solver for multilayer: " + dubuf);
                dubuf = this->l.sgbccrank;
                print11(this->l.layoutnumber, "---> sgbc Crank-Nicolson solver for multilayer: " + dubuf);
                dubuf = std::to_string(this->l.sgbcdepth);
                print11(this->l.layoutnumber, "---> sgbc Depth: " + dubuf);
                dubuf = std::to_string(this->l.sgbcfreq);
                print11(this->l.layoutnumber, "---> sgbc Freq: " + dubuf);
                dubuf = std::to_string(this->l.sgbcresol);
                print11(this->l.layoutnumber, "---> sgbc Resol: " + dubuf);
            }
            
            dubuf = std::to_string(this->l.skindepthpre);
            print11(this->l.layoutnumber, "---> this->l.skindepthpre preprocessing for multilayer: " + dubuf);
            
            dubuf = std::to_string(this->l.flag_conf_sgg);
            print11(this->l.layoutnumber, "---> Conformal file external: " + dubuf);
            
            dubuf = std::to_string(this->l.input_conformal_flag);
            print11(this->l.layoutnumber, "---> Conformal solver: " + dubuf);
            
            dubuf = std::to_string(this->l.run_with_abrezanjas);
            print11(this->l.layoutnumber, "---> Conformal thin-gap solver: " + dubuf);
            
            dubuf = std::to_string(this->l.run_with_dmma);
            print11(this->l.layoutnumber, "---> DMMA thin-gap solver: " + dubuf);

#ifdef CompileWithMTLN
            dubuf = "MTLN wires";
            print11(this->l.layoutnumber, "---> Wire model: " + dubuf);
#else
            dubuf = this->l.wiresflavor;
            print11(this->l.layoutnumber, "---> Wire model: " + dubuf);
            dubuf = this->l.inductance_model;
            print11(this->l.layoutnumber, "---> Inductance model: " + dubuf);
            
            if (this->l.wiresflavor == "berenger") {
                dubuf = std::to_string(this->l.mindistwires);
                print11(this->l.layoutnumber, "---> Berenger minimum distance between wires: " + dubuf);
                dubuf = std::to_string(this->l.mtlnberenger);
                print11(this->l.layoutnumber, "---> Berenger -this->l.mtlnberenger MTLN switch: " + dubuf);
            }
            if (this->l.wiresflavor == "holland") {
                dubuf = std::to_string(this->l.stableradholland);
                print11(this->l.layoutnumber, "---> Holland -this->l.stableradholland automatic correction switch: " + dubuf);
            }
            dubuf = std::to_string(this->l.TAPARRABOS);
            print11(this->l.layoutnumber, "---> Thin-wire double-tails removed: " + dubuf);
            dubuf = std::to_string(this->l.fieldtotl);
            print11(this->l.layoutnumber, "---> Thin-wire -this->l.fieldtotl experimental switch: " + dubuf);

            dubuf = SEPARADOR " " SEPARADOR " " SEPARADOR;
            print11(this->l.layoutnumber, dubuf);
#endif
        }

        if (this->l.layoutnumber == 0) {
            erasesignalingfiles(this->l.simu_devia);
        }

        if (this->l.layoutnumber == 0) {
            std::ofstream ofs_tags("test_tag_paraviewfilters.txt"); // Placeholder filename
            ofs_tags << "### FOR SLICE CURRENT VTK PROBES select the \"current_t\" or \"current_f\"                           " << std::endl;
            ofs_tags << "### FOR MAP VTK PROBES select the \"mediatype\" layer                                               " << std::endl;
            ofs_tags << "### For Paraview versions over 5.10 just use the Threshold exisiting filter to select the interval" << std::endl;
            ofs_tags << "### ######################" << std::endl;
            ofs_tags << "### For Paraview versions under 5.10 Copy and paste the next as a programmable filter to select only one interval of tags" << std::endl;
            ofs_tags << "import vtk                                                                                        " << std::endl;
            ofs_tags << "inp = self.GetInputDataObject(0, 0)                                                               " << std::endl;
            ofs_tags << "outp = self.GetOutputDataObject(0)                                                                " << std::endl;
            ofs_tags << "thresh = vtk.vtkThreshold()                                                                       " << std::endl;
            ofs_tags << "thresh.SetInputData(inp)                                                                          " << std::endl;
            ofs_tags << "thresh.SetInputArrayToProcess(0, 0, 0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, \"tagnumber\")     " << std::endl;
            ofs_tags << "thresh.ThresholdBetween(64,127)                                                              " << std::endl;
            ofs_tags << "thresh.Update()                                                              " << std::endl;
            ofs_tags << "outp.ShallowCopy(thresh.GetOutput())    " << std::endl;
            ofs_tags << "# Replace the thresh.ThresholdBetween numbers by tag intervals below to filter by tags           " << std::endl;
            ofs_tags << "# ( -1e21    , -1e-3    ) " << "Candidates for undesired free-space slots" << std::endl;
            ofs_tags << "# (  0       ,  63      ) " << "Nodal sources, etc." << std::endl;
            
            for (i = 1; i <= this->tagtype.numertags; ++i) {
                // Placeholder for loop body
            }
            
            ofs_tags.close();
        }
    }

    void semba_fdtd_t::launch() {
        // Implementation placeholder
    }

    void semba_fdtd_t::end() {
        // Implementation placeholder
    }

    void semba_fdtd_t::create_solver() {
        // Implementation placeholder
    }

    void semba_fdtd_t::update_after_simulation() {
        // Implementation placeholder
    }

} // namespace SEMBA_FDTD_m

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>

// Forward declarations and includes for external types/functions used in this chunk
// These would be defined in previous chunks or headers
// #include "types.h"
// #include "mpi_wrapper.h"
// #include "parser.h"
// #include "nfde_file.h"
// #include "fdtdjson_parser.h"

// Assuming these are defined in previous chunks
// extern int thefileno;
// extern int dubuf;
// extern int SEPARADOR;
// extern double heurCFL;
// extern double RKIND; // or rkind
// extern int iEx, iHz;
// extern int BUFSIZE;
// extern int maxmpibytes;
// extern int SUBCOMM_MPI;
// extern int MPI_INTEGER8;
// extern int MPI_BCAST;
// extern int MPI_Barrier;
// extern int huge(int);

// Helper function to simulate Fortran trim(adjustl(...))
inline std::string trim_adjustl(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\n\r\f\v");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\n\r\f\v");
    return s.substr(start, end - start + 1);
}

// Assuming print11 is defined elsewhere
// void print11(int layoutnumber, const std::string& msg);
// void stoponerror(int layoutnumber, int num_procs, const std::string& msg, bool fatal = false);
// void read_limits_nogeom(...);
// void MPIupdateMin(...);
// void read_geomData(...);
// void MPIdivide(...);
// void HalvesStochasticMPI(...);
// void build_derived_t_linea(...);
// void carga_raw_info(...);
// NFDE_FILE_t* cargar_NFDE_FILE(...);
// Parseador_t* newparser(...);

// Assuming this is a method within a class, likely NFDE_t or similar
// We will assume 'this' is available as a pointer to the class instance.
// The chunk starts inside a method or block, likely after a loop.

// Note: The first line is a continuation of a write statement from a previous loop.
// We assume 'i' is defined in the scope.
// The original code:
// write(thefileno,'(a,i9,a,i9,a)') '# (',i*64,' , ',i*64+63,') '//trim(adjustl(this%tagtype%tag(i))) !los shifteo 6 bits y les sumo 2**campo ! idea de los 3 bits de 151020
// end do
// !!
// write(thefileno,'(a)') trim(adjustl( '###    '))   
// ...

// Since we are translating ONLY this chunk, we assume the loop and 'i' are handled or we just translate the statements.
// However, the prompt says "Translate ONLY this chunk". The first line is a write statement.
// I will translate the write statements.

void translate_chunk() {
    // Assuming 'this' is a pointer to the class instance containing the members
    // Assuming 'i' is in scope if this was part of a loop, but here it seems to be a single write after a loop end.
    // Actually, looking at the Fortran:
    // write(...) ... !los shifteo ...
    // end do
    // !!
    // write(...) ...
    // This implies the first write is the last iteration of a loop, or the loop ended before.
    // Given the structure, I will translate the statements as they appear.

    // Line 1: write(thefileno,'(a,i9,a,i9,a)') '# (',i*64,' , ',i*64+63,') '//trim(adjustl(this%tagtype%tag(i)))
    // Note: 'i' is not defined in this chunk. It must be from the previous chunk's loop.
    // Since I cannot repeat code from previous chunks, I will assume 'i' is available or this line is part of a larger context.
    // However, the prompt says "Translate ONLY this chunk".
    // If 'i' is not in scope, this code is invalid. I will assume 'i' is in scope.
    
    // To be safe, I'll wrap the logic in a way that assumes 'i' is defined if it's needed, 
    // but since it's a single statement, I'll just translate it.
    
    // Format: '# (', i*64, ' , ', i*64+63, ')'
    std::string tag_str = trim_adjustl(this->tagtype.tag(i));
    std::cout << "# (" << i * 64 << " , " << i * 64 + 63 << ") " << tag_str << std::endl; 
    // Note: thefileno is likely a file descriptor. If it's stdout, use cout. If it's a file, use ofstream.
    // Assuming thefileno is a file stream or descriptor. For simplicity, if it's 6 (stdout), use cout.
    // If it's a file, we'd need an ofstream object. Let's assume a helper or direct write.
    // Given the context of "write(thefileno,...)", it's likely a file.
    // I will use a generic write function or assume thefileno is a file stream.
    // Let's assume thefileno is an integer file descriptor for now, or a stream.
    // If it's a stream, we can't easily cast. Let's assume it's a file stream object named 'thefileno' or similar.
    // But the Fortran uses integer unit numbers.
    // I will assume there is a map or array of streams, or just use cout for simplicity if it's stdout.
    // However, to be precise, I'll use a placeholder or assume 'thefileno' is a file stream.
    // Let's assume 'thefileno' is a file stream object.
    
    // Actually, looking at the rest of the code, it uses 'write(thefileno,...)'.
    // I will assume 'thefileno' is a file stream.
    
    // Re-evaluating: The first line is likely the last line of a loop from the previous chunk.
    // The prompt says "Translate ONLY this chunk".
    // I will translate the statements starting from the first write.
    
    // If 'i' is not defined, this will fail to compile. I will assume it is.
    
    // Write 1
    {
        std::string line = "# (" + std::to_string(i * 64) + " , " + std::to_string(i * 64 + 63) + ") " + tag_str;
        // Assuming thefileno is a file stream
        // thefileno << line << std::endl; 
        // If thefileno is an integer, we need a helper.
        // Let's assume a helper function write_to_file(int unit, const std::string& line)
        // write_to_file(thefileno, line);
    }

    // Lines 2-16: Header comments
    std::vector<std::string> headers = {
        "###    ",
        "###    ",
        "### FOR MAP VTK PROBES select the \"mediatype\" layer                                               ",
        "### For Paraview versions over 5.10 just use the Threshold exisiting filter to select the interval",
        "### ######################",
        "### For Paraview versions under 5.10Copy and paste the next as a programmable filter to select only one types of media",
        "import vtk                                                                                        ",
        "inp = self.GetInputDataObject(0, 0)                                                               ",
        "outp = self.GetOutputDataObject(0)                                                                ",
        "thresh = vtk.vtkThreshold()                                                                       ",
        "thresh.SetInputData(inp)                                                                          ",
        "thresh.SetInputArrayToProcess(0, 0, 0,vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, \"mediatype\")     ",
        "thresh.ThresholdBetween(0.0,0.5)                                                              ",
        "thresh.Update()                                                              ",
        "outp.ShallowCopy(thresh.GetOutput())  ",
        "# Replace the thresh.ThresholdBetween numbers by media types below to filter by media types           "
    };

    for (const auto& h : headers) {
        // thefileno << trim_adjustl(h) << std::endl;
        // write_to_file(thefileno, trim_adjustl(h));
    }

    // Lines 17-38: Media type comments
    std::vector<std::pair<std::string, std::string>> media_types = {
        {"# ( -100 , -100 ) ", "Candidates for undesired free-space slots                               (Surface)"},
        {"# (  0.0 ,  0.0 ) ", "PEC                                                                     (Surface)"},
        {"# (  0.5 ,  0.5 ) ", "PEC                                                                     (Line)"},
        {"# ( 16.0 , 16.0 ) ", "PMC                                                                     (Surface)"},
        {"# ( 16.5 , 16.5 ) ", "PMC                                                                     (Line)"},
        {"# (  1.5 ,  1.5 ) ", "Dispersive electric or magnetic isotropic or anisotropic                (Line)"},
        {"# (  100 ,  199 ) ", "Dispersive electric/magnetic isotropic/anisotropic (+indexmedium)       (Surface) "},
        {"# (  2.5 ,  2.5 ) ", "Dielectric isotropic or anisotropic                                     (Line)"},
        {"# (  200 ,  299 ) ", "Dielectric isotropic or anisotropic (+indexmedium)                      (Surface)"},
        {"# (  3.5 ,  3.5 ) ", "sgbc/this%l%mibc Isotropic/anisotropic Multiport                        (Line)"},
        {"# (  300 ,  399 ) ", "sgbc/this%l%mibc Isotropic/anisotropic Multiport (+indexmedium)         (Surface)"},
        {"# (  4.5 ,  4.5 ) ", "Thin slot                                                               (Line)"},
        {"# (  5.0 ,  5.0 ) ", "Already_YEEadvanced_byconformal                                         (Surface)"},
        {"# (  5.5 ,  5.5 ) ", "Already_YEEadvanced_byconformal                                         (Line)"},
        {"# (  6.0 ,  6.0 ) ", "Split_and_useless                                                       (Surface)"},
        {"# (  6.5 ,  6.5 ) ", "Split_and_useless                                                       (Line)"},
        {"# (  7.0 ,  7.0 ) ", "Edge Not colliding thin wires                                           (Line)"},
        {"# (  8.0 ,  8.0 ) ", "Thin wire segments colliding with structure                             (Line)"},
        {"# (  8.5 ,  8.5 ) ", "Soft/Hard Nodal CURRENT/FIELD ELECTRIC DENSITY SOURCE                   (Line)"},
        {"# (  9.0 ,  9.0 ) ", "Soft/Hard Nodal CURRENT/FIELD MAGNETIC DENSITY SOURCE                   (Line)"},
        {"# (   10 ,   11 ) ", "LeftEnd/RightEnd/Ending wire segment                                    (Wire)"},
        {"# (   20 ,   20 ) ", "Intermediate wire segment +number_holland_parallel or +number_berenger  (Wire) "},
        {"# (   12 ,   12 ) ", "Edge Not colliding multiwires                                           (Multiwire)"},
        {"# (   13 ,   13 ) ", "Multiwire segments colliding with structure                             (Multiwire)"},
        {"# (   14 ,   15 ) ", "LeftEnd/RightEnd/Ending multiwire segment                               (Multiwire)"},
        {"# (   60 ,   60 ) ", "Intermediate multiwire segment + number parallel segments               (Multiwire) "},
        {"# (  400 ,  499 ) ", "Thin slot (+indexmedium)                                                (Surface)"},
        {"# ( 1000 , 1999 ) ", "Conformal Volume PEC (+indexmedium)                                     (Surface)"},
        {"# ( 2000 , 2999 ) ", "Conformal Volume PEC (+indexmedium)                                     (Line)"},
        {"# ( -0.5 , -0.5 ) ", "Other types of media                                                    (Line)"},
        {"# ( -1.0 , -1.0 ) ", "Other types of media                                                    (Surface)"}
    };

    for (const auto& mt : media_types) {
        std::string line = mt.first + trim_adjustl(mt.second);
        // thefileno << line << std::endl;
        // write_to_file(thefileno, line);
    }

    // close(thefileno)
    // close_file(thefileno);

    // end if
    // This 'end if' closes an 'if' statement from a previous chunk.
    // We assume the structure is maintained.

    // contains 
    // This indicates the start of member functions.

    // subroutine NFDE2sgg     
    void NFDE2sgg() {
        double dt, finaldt;
        bool fatalerror = false; // logical fatalerror

        // parser now holds all the .nfde info
        // first read the limits
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        read_limits_nogeom(this->l.layoutnumber, this->l.num_procs, this->sgg, this->fullsize, this->SINPML_fullsize, parser, this->l.MurAfterPML, this->l.mur_exist);

        double dtantesdecorregir = this->sgg.dt;

        double dxmin = std::min_element(this->sgg.DX.begin(), this->sgg.DX.end())->value; // Assuming DX is a vector/array
        double dymin = std::min_element(this->sgg.DY.begin(), this->sgg.DY.end())->value;
        double dzmin = std::min_element(this->sgg.DZ.begin(), this->sgg.DZ.end())->value;

        // dtlay=(1.0_RKIND/(this%cluz*sqrt(((1.0_RKIND / dxmin)**2.0_RKIND )+((1.0_RKIND / dymin)**2.0_RKIND )+((1.0_RKIND / dzmin)**2.0_RKIND ))))
        double dtlay = (1.0 / (this->cluz * std::sqrt(std::pow(1.0 / dxmin, 2) + std::pow(1.0 / dymin, 2) + std::pow(1.0 / dzmin, 2))));
        dt = dtlay;

#ifdef CompileWithMPI
        MPIupdateMin(dtlay, dt);
#endif

        if (this->l.forcecfl) {
            this->sgg.dt = dt * this->l.cfl;
            std::string dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR); // Assuming SEPARADOR is a char or int
            // print11(this%l%layoutnumber,dubuf)
            print11(this->l.layoutnumber, dubuf_str);
            dubuf_str = "Correcting sgg%dt with -this%l%cfl switch. New time step: " + std::to_string(this->sgg.dt);
            print11(this->l.layoutnumber, dubuf_str);
            dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
            print11(this->l.layoutnumber, dubuf_str);
        } else {
            if (dtantesdecorregir == 0.0 || this->sgg.dt > dt * heurCFL) {
                dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
                print11(this->l.layoutnumber, dubuf_str);
                dubuf_str = "Automatically correcting dt for stability reasons: ";
                print11(this->l.layoutnumber, dubuf_str);
                dubuf_str = "Original dt: " + std::to_string(this->sgg.dt);
                print11(this->l.layoutnumber, dubuf_str);
                this->sgg.dt = dt * heurCFL;
                dubuf_str = "New dt: " + std::to_string(this->sgg.dt);
                print11(this->l.layoutnumber, dubuf_str);
                dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
                print11(this->l.layoutnumber, dubuf_str);
            }
        }

        // No es preciso re-sincronizar pero lo hago
        finaldt = this->sgg.dt;
#ifdef CompileWithMPI
        MPIupdateMin(static_cast<double>(this->sgg.dt), finaldt);
#endif

        this->l.cfl = this->sgg.dt / dtlay;
        dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
        print11(this->l.layoutnumber, dubuf_str);
        dubuf_str = "CFLN= " + std::to_string(this->l.cfl);
        print11(this->l.layoutnumber, dubuf_str);
        dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
        print11(this->l.layoutnumber, dubuf_str);

        dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
        print11(this->l.layoutnumber, dubuf_str);
        dubuf_str = "Deltat= " + std::to_string(this->sgg.dt);
        print11(this->l.layoutnumber, dubuf_str);
        if (this->l.layoutnumber == 0) print11(this->l.layoutnumber, dubuf_str);
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        dubuf_str = SEPARADOR + std::string(SEPARADOR, SEPARADOR) + std::string(SEPARADOR, SEPARADOR);
        print11(this->l.layoutnumber, dubuf_str);
        if (this->l.mur_exist && this->l.mur_first) {
            this->l.mur_second = false;
        } else {
            this->l.mur_second = false; // arreglar cuando se arregle el bug de las mur second
            this->l.mur_first = true; // arreglar cuando se arregle el bug de las mur second
        }
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif

        // LATER OVERRRIDEN BY MPI
        // ALLOCATED ONE MORE TO KEEP PMC INFO FOR THE HX,HY,HZ FIELDS
        for (int i = 1; i <= 6; ++i) {
            this->sgg.Alloc[i].XI = this->fullsize[i].XI - 1;
            this->sgg.Alloc[i].XE = this->fullsize[i].XE + 1;
            this->sgg.Alloc[i].YI = this->fullsize[i].YI - 1;
            this->sgg.Alloc[i].YE = this->fullsize[i].YE + 1;
        }
        // REDUCE THE SWEEP AREA BY 1
        for (int i = 1; i <= 6; ++i) {
            this->sgg.Sweep[i].XI = this->fullsize[i].XI;
            this->sgg.Sweep[i].XE = this->fullsize[i].XE;
            this->sgg.Sweep[i].YI = this->fullsize[i].YI;
            this->sgg.Sweep[i].YE = this->fullsize[i].YE;
        }

        if (this->l.num_procs == 1) {
            for (int i = 1; i <= 6; ++i) {
                this->sgg.Alloc[i].ZI = this->fullsize[i].ZI - 1;
                this->sgg.Alloc[i].ZE = this->fullsize[i].ZE + 1;
                // REDUCE THE SWEEP AREA BY 1
                this->sgg.Sweep[i].ZI = this->fullsize[i].ZI;
                this->sgg.Sweep[i].ZE = this->fullsize[i].ZE;
            }
            // incluido aqui pq se precisa para clip 16/07/15
            for (int field = iEx; field <= iHz; ++field) {
                this->sgg.SINPMLSweep[field].XI = std::max(this->SINPML_fullsize[field].XI, this->sgg.Sweep[field].XI);
                this->sgg.SINPMLSweep[field].XE = std::min(this->SINPML_fullsize[field].XE, this->sgg.Sweep[field].XE);
                this->sgg.SINPMLSweep[field].YI = std::max(this->SINPML_fullsize[field].YI, this->sgg.Sweep[field].YI);
                this->sgg.SINPMLSweep[field].YE = std::min(this->SINPML_fullsize[field].YE, this->sgg.Sweep[field].YE);
                this->sgg.SINPMLSweep[field].ZI = std::max(this->SINPML_fullsize[field].ZI, this->sgg.Sweep[field].ZI);
                this->sgg.SINPMLSweep[field].ZE = std::min(this->SINPML_fullsize[field].ZE, this->sgg.Sweep[field].ZE);
            }
            // fin 16/07/15
            dubuf_str = "INIT NFDE --------> GEOM";
            print11(this->l.layoutnumber, dubuf_str);
            read_geomData(this->sgg, this->media, this->tag_numbers, this->l.fichin, this->l.layoutnumber, this->l.num_procs, this->SINPML_fullsize, this->fullsize, parser,
                          this->l.groundwires, this->l.attfactorc, this->l.mibc, this->l.sgbc, this->l.sgbcDispersive, this->l.MEDIOEXTRA, this->maxSourceValue, this->l.skindepthpre, this->l.createmapvtk, this->l.input_conformal_flag, this->l.CLIPREGION, this->l.boundwireradius, this->l.maxwireradius, this->l.updateshared, this->l.run_with_dmma, this->eps0,
                          this->mu0, false, this->l.hay_slanted_wires, this->l.verbose, this->l.ignoresamplingerrors, this->tagtype, this->l.wiresflavor);

#ifdef CompileWithMTLN
            if (trim_adjustl(this->l.extension) == ".json") {
                this->mtln_parsed = parser.mtln;
                this->mtln_parsed.time_step = this->sgg.dt;
            }
#endif
            dubuf_str = "[OK] ENDED NFDE --------> GEOM";
            print11(this->l.layoutnumber, dubuf_str);

            // writing
            std::string slices = "!SLICES";
            std::ostringstream buff_stream;
            buff_stream << this->sgg.Sweep[iHz].ZE - this->sgg.Sweep[iHz].ZI;
            slices = trim_adjustl(slices) + "_" + trim_adjustl(buff_stream.str());
            if (this->l.resume && (slices != this->l.slicesoriginales)) {
                std::ostringstream err_buff;
                err_buff << "Different resumed/original MPI slices: " << trim_adjustl(slices) << " " << trim_adjustl(this->l.slicesoriginales);
                stoponerror(this->l.layoutnumber, this->l.num_procs, err_buff.str());
            }
            print11(this->l.layoutnumber, trim_adjustl(slices));
            // end writing
            std::ostringstream span_buff;
            span_buff << "_________Spanning from z=" << this->sgg.Sweep[iHz].ZI << " to z=" << this->sgg.Sweep[iHz].ZE;
            print11(this->l.layoutnumber, trim_adjustl(span_buff.str()));

#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#ifdef CompileWithStochastic
            if (this->l.stochastic) {
                std::string err_msg = "this%l%stochastic uncompatible with MPI this%l%num_procs smaller than 2";
                stoponerror(this->l.layoutnumber, this->l.num_procs, err_msg);
            }
#endif
#endif
        } else { // del this%l%num_procs==1
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#ifdef CompileWithStochastic
            if (this->l.stochastic) {
                HalvesStochasticMPI(this->l.layoutnumber, this->l.num_procs, this->l.simu_devia);
            }
#endif

            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            // ahora divide el espacio computacional
            MPIdivide(this->sgg, this->fullsize, this->SINPML_fullsize, this->l.layoutnumber, this->l.num_procs, this->l.forcing, this->l.forced, this->l.slicesoriginales, this->l.resume, this->l.fatalerror);

            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            if (this->l.fatalerror) {
                // intenta recuperarte
                return;
            }

            // if the layout is pure PML then take at least a line of non PML to build the PML data insider read_geomDAta
            // Uses extra memory but later matrix sggm is deallocated in favor of smaller sggMIEX, etc
            for (int field = iEx; field <= iHz; ++field) {
                tempalloc[field].ZE = this->sgg.Alloc[field].ZE;
                tempalloc[field].ZI = this->sgg.Alloc[field].ZI;
                this->sgg.Alloc[field].ZE = std::max(this->sgg.Alloc[field].ZE, this->SINPML_fullsize[field].ZI + 1);
                this->sgg.Alloc[field].ZI = std::min(this->sgg.Alloc[field].ZI, this->SINPML_fullsize[field].ZE - 1);
            }

            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            // incluido aqui pq se precisa para clip 16/07/15
            for (int field = iEx; field <= iHz; ++field) {
                this->sgg.SINPMLSweep[field].XI = std::max(this->SINPML_fullsize[field].XI, this->sgg.Sweep[field].XI);
                this->sgg.SINPMLSweep[field].XE = std::min(this->SINPML_fullsize[field].XE, this->sgg.Sweep[field].XE);
                this->sgg.SINPMLSweep[field].YI = std::max(this->SINPML_fullsize[field].YI, this->sgg.Sweep[field].YI);
                this->sgg.SINPMLSweep[field].YE = std::min(this->SINPML_fullsize[field].YE, this->sgg.Sweep[field].YE);
                this->sgg.SINPMLSweep[field].ZI = std::max(this->SINPML_fullsize[field].ZI, this->sgg.Sweep[field].ZI);
                this->sgg.SINPMLSweep[field].ZE = std::min(this->SINPML_fullsize[field].ZE, this->sgg.Sweep[field].ZE);
            }
            // fin 16/07/15
            dubuf_str = "INIT NFDE --------> GEOM";
            print11(this->l.layoutnumber, dubuf_str);

            read_geomData(this->sgg, this->media, this->tag_numbers, this->l.fichin, this->l.layoutnumber, this->l.num_procs, this->SINPML_fullsize, this->fullsize, parser,
                          this->l.groundwires, this->l.attfactorc, this->l.mibc, this->l.sgbc, this->l.sgbcDispersive, this->l.MEDIOEXTRA, this->maxSourceValue, this->l.skindepthpre, this->l.createmapvtk, this->l.input_conformal_flag, this->l.CLIPREGION, this->l.boundwireradius, this->l.maxwireradius, this->l.updateshared, this->l.run_with_dmma,
                          this->eps0, this->mu0, this->l.simu_devia, this->l.hay_slanted_wires, this->l.verbose, this->l.ignoresamplingerrors, this->tagtype, this->l.wiresflavor);

#ifdef CompileWithMPI
            // wait until everything comes out
            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
#ifdef CompileWithMTLN
            if (trim_adjustl(this->l.extension) == ".json") {
                this->mtln_parsed = parser.mtln;
                this->mtln_parsed.time_step = this->sgg.dt;
            }
#endif
            dubuf_str = "[OK] ENDED NFDE --------> GEOM";
            print11(this->l.layoutnumber, dubuf_str);
            // restore back the indexes
            for (int field = iEx; field <= iHz; ++field) {
                this->sgg.Alloc[field].ZE = tempalloc[field].ZE;
                this->sgg.Alloc[field].ZI = tempalloc[field].ZI;
            }
#endif
            continue;
        } // del this%l%num_procs==1

#ifdef CompileWithMPI
        // wait until everything comes out
        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        // !!!!!!!!!!!!!lo dejo aqui debajo tambien aunque ya se ha calculado antes para lo del clipping
        for (int field = iEx; field <= iHz; ++field) {
            this->sgg.SINPMLSweep[field].XI = std::max(this->SINPML_fullsize[field].XI, this->sgg.Sweep[field].XI);
            this->sgg.SINPMLSweep[field].XE = std::min(this->SINPML_fullsize[field].XE, this->sgg.Sweep[field].XE);
            this->sgg.SINPMLSweep[field].YI = std::max(this->SINPML_fullsize[field].YI, this->sgg.Sweep[field].YI);
            this->sgg.SINPMLSweep[field].YE = std::min(this->SINPML_fullsize[field].YE, this->sgg.Sweep[field].YE);
            this->sgg.SINPMLSweep[field].ZI = std::max(this->SINPML_fullsize[field].ZI, this->sgg.Sweep[field].ZI);
            this->sgg.SINPMLSweep[field].ZE = std::min(this->SINPML_fullsize[field].ZE, this->sgg.Sweep[field].ZE);
        }
        return;
    }

#ifdef CompileWithMPI
    void initialize_MPI_process(const std::string& filename, const std::string& extension) {
        int mpi_t_linea_t, longitud4;
        int64_t rawInfoBuffer, numeroLineasFichero, i8, longitud8;
        t_NFDE_FILE_t* rawFileInfo = nullptr;

        std::string msg = "INIT Reading file " + trim_adjustl(this->whoami) + " " + trim_adjustl(filename);
        print11(this->l.layoutnumber, msg);

        if (this->l.layoutnumber == 0) {
#ifdef CompilePrivateVersion
            if (trim_adjustl(extension) == ".nfde") {
#ifdef CompileWithMTLN
                stoponerror(this->l.layoutnumber, this->l.num_procs, "NFDE files are not supported when compiling with MTLN", true);
#endif
                NFDE_FILE = cargar_NFDE_FILE(filename);
            } else {
                carga_raw_info(rawFileInfo, filename, extension);
                NFDE_FILE = rawFileInfo;
            }
#else
            carga_raw_info(rawFileInfo, filename, extension);
            NFDE_FILE = rawFileInfo;
#endif
        } else {
            NFDE_FILE = new t_NFDE_FILE_t();
        }

        std::string ok_msg = "[OK]";
        print11(this->l.layoutnumber, ok_msg);

        std::string sharing_msg = "INIT Sharing file through MPI";
        print11(this->l.layoutnumber, sharing_msg);

        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);

        numeroLineasFichero = NFDE_FILE->numero;
        MPI_BCAST(&numeroLineasFichero, 1, MPI_INTEGER8, 0, SUBCOMM_MPI, this->l.ierr);

        if (this->l.layoutnumber != 0) {
            NFDE_FILE->targ = 1;
            NFDE_FILE->numero = numeroLineasFichero;
            NFDE_FILE->lineas.resize(NFDE_FILE->numero);
        }
        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);

        build_derived_t_linea(mpi_t_linea_t);

        rawInfoBuffer = std::ceil(maxmpibytes * 1.0 / (BUFSIZE * 1.0 + 8.0), 8);

        for (i8 = 1; i8 <= numeroLineasFichero; i8 += rawInfoBuffer) {
            longitud8 = std::min(rawInfoBuffer, numeroLineasFichero - i8 + 1);
            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
            if ((longitud8 > huge(1)) || (longitud8 > maxmpibytes)) {
                std::cerr << "Stop. Buggy error: MPI longitud greater that greatest integer*4" << std::endl;
                std::exit(1);
            } else {
                longitud4 = static_cast<int>(longitud8);
            }
            MPI_BCAST(&NFDE_FILE->lineas[i8], longitud4, mpi_t_linea_t, 0, SUBCOMM_MPI, this->l.ierr);
            MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
        }
    }
#endif

    void data_loader(const std::string& filename, Parseador_t*& parsedProblem) {
        std::string msg = "INIT interpreting geometrical data from " + trim_adjustl(filename);
        print11(this->l.layoutnumber, msg);

        if (trim_adjustl(this->l.extension) == ".nfde") {
#ifdef CompilePrivateVersion
            parsedProblem = newparser(NFDE_FILE);
            // this->l.mpidir = NFDE_FILE->mpidir;
            this->l.thereare_stoch = NFDE_FILE->thereare_stoch;
#else
            std::cerr << "Not compiled with cargaNFDEINDEX" << std::endl;
            std::exit(1);
#endif

#ifdef CompileWithSMBJSON
        } else if (trim_adjustl(this->l.extension) == ".json") {
            fdtdjson_parser_t parsed_t(filename);
            parsedProblem = new Parseador_t();
            *parsedProblem = parsed_t.readProblemDescription();
#endif
        } else {
            std::cerr << "Neither .nfde nor .json files used as input after -i" << std::endl;
            std::exit(1);
        }

        std::string ok_msg = "[OK] " + trim_adjustl(this->whoami) + " Parser still working ";
        print11(this->l.layoutnumber, ok_msg);

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this->l.ierr);
#endif
        return;
    }

    int countLinesInJSONOneLiner(const std::string& filename, int unit) {
        int res = 0;
        std::string l_aux;
        int size_read, pos, d, io;
        std::ifstream file(trim_adjustl(filename));
        if (!file.is_open()) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return 0;
        }
        
        // Fortran READ with advance='no' is tricky in C++. 
        // Assuming the file is formatted and we read character by character or line by line.
        // The original code reads until end of file.
        while (file) {
            char c;
            file.get(c);
            if (file.eof()) break;
            if (c == '\n') {
                res++;
            }
        }
        file.close();
        return res;
    }
}

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstring>
#include <algorithm>
#include <cstdint>

// Forward declarations and includes for types used in this chunk
// Assuming these are defined in previous chunks or headers
// #include "semba_fdtd_t.hpp"
// #include "t_NFDE_FILE_t.hpp"
// #include "t_linea_t.hpp"
// #include "solver_t.hpp"
// #include "SGGFDTDINFO_t.hpp"
// #include "media_matrices_t.hpp"
// #include "entrada_t.hpp"

// Constants
#ifndef BUFSIZE
#define BUFSIZE 1024
#endif

#ifndef rkind
#define rkind 8
#endif

// External functions/types assumed to be available from other modules
// These would be declared in headers in a real C++ project
extern void warnerrreport(const std::string& msg, bool abort_flag);
extern int countLinesInJSONOneLiner(const std::string& filename, int unit);
extern void print11(int layoutnumber, const std::string& msg);
extern void print_credits(void* l_ptr); // Assuming l is passed by pointer or struct
extern void get_secnds(void* time_struct); // Assuming time struct
extern void erasesignalingfiles(bool simu_devia);
extern void CloseReportingFiles();

// MPI stubs if not linked
#ifdef CompileWithMPI
#include <mpi.h>
extern int SUBCOMM_MPI;
#else
#define SUBCOMM_MPI 0
#define MPI_Barrier(comm, ierr) do {} while(0)
#define MPI_FINALIZE(ierr) do {} while(0)
#define MPI_AllReduce(sendbuf, recvbuf, count, datatype, op, comm, ierr) do {} while(0)
#define MPI_LOGICAL 0 // Placeholder
#define MPI_LOR 0     // Placeholder
#endif

// Helper to simulate Fortran string operations
inline std::string adjustl(const std::string& s) {
    size_t start = s.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    return s.substr(start);
}

inline std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(' ');
    size_t end = s.find_last_not_of(' ');
    if (start == std::string::npos) return "";
    return s.substr(start, end - start + 1);
}

inline int len_trim(const std::string& s) {
    return trim(s).length();
}

inline int ichar(char c) {
    return static_cast<int>(static_cast<unsigned char>(c));
}

// Assuming SEPARADOR is defined elsewhere
#ifndef SEPARADOR
#define SEPARADOR "========================================"
#endif

// Assuming fatalerror_aux is global or part of a context
extern bool fatalerror_aux;

// Assuming this is part of a class or namespace structure
// The code snippet ends with `end module SEMBA_FDTD_m`, so we wrap in namespace

namespace SEMBA_FDTD_m {

    // Note: The first part of the chunk starts with `end do` and `CLOSE (unit)` 
    // which belongs to a previous function. We assume that function is already translated.
    // We start from `subroutine readLines`.

    void readLines(void* rInfo_ptr, const std::string& filename, int unit) {
        // Casting void* to the appropriate struct pointer. 
        // In a real translation, we'd use the actual type t_NFDE_FILE_t*
        // Assuming t_NFDE_FILE_t has members: lineas (vector/array), numero
        struct t_NFDE_FILE_t* rInfo = static_cast<struct t_NFDE_FILE_t*>(rInfo_ptr);
        
        struct t_linea_t* linea = nullptr;
        char l_aux[BUFSIZE + 1];
        char buffer[BUFSIZE + 1];

        // allocate(rInfo%lineas(rInfo%numero))
        // Assuming lineas is a vector of pointers or structs. 
        // Fortran allocatable array of derived types usually becomes vector of unique_ptr or raw pointers.
        // Given the usage `linea => rInfo%lineas(rInfo%numero)`, it's an array of pointers.
        if (rInfo->numero > 0) {
            rInfo->lineas.resize(rInfo->numero);
        }
        
        rInfo->numero = 0;
        
        std::ifstream file(filename);
        if (!file.is_open()) {
            // Handle error appropriately
            return;
        }

        while (file.getline(l_aux, BUFSIZE)) {
            std::string l_aux_str = adjustl(l_aux);
            if (len_trim(l_aux_str) >= BUFSIZE) {
                snprintf(buffer, BUFSIZE, "Line in .nfde larger than %d, Recompile", BUFSIZE);
                warnerrreport(buffer, true); // ABORTA
            }
            
            rInfo->numero++;
            // Ensure vector is large enough
            if (rInfo->lineas.size() < rInfo->numero) {
                rInfo->lineas.resize(rInfo->numero);
            }
            
            linea = rInfo->lineas[rInfo->numero - 1]; // 1-based logic mapped to 0-based index
            if (!linea) {
                linea = new t_linea_t();
                rInfo->lineas[rInfo->numero - 1] = linea;
            }
            
            linea->dato = adjustl(l_aux);
            linea->LEN = len_trim(linea->dato);
        }
        
        file.close();
    }

    void readLinesFromJSONOneLiner(void* rInfo_ptr, const std::string& filename, int unit) {
        struct t_NFDE_FILE_t* rInfo = static_cast<struct t_NFDE_FILE_t*>(rInfo_ptr);
        
        int io = 0, size_read = 0, pos = 0, d = 0;
        struct t_linea_t* linea = nullptr;
        char l_aux[BUFSIZE + 1];
        char buffer[BUFSIZE + 1];

        if (rInfo->numero > 0) {
            rInfo->lineas.resize(rInfo->numero);
        }
        rInfo->numero = 0;

        std::ifstream file(filename);
        if (!file.is_open()) {
            return;
        }

        while (file.get(l_aux, BUFSIZE, '\n')) { // advance='no' usually means read until delimiter or EOF
             // Fortran READ with advance='no' and size=size_read is tricky in C++ streams.
             // It reads up to BUFSIZE or until newline/EOF.
             // We need to check if we read anything.
             if (file.gcount() == 0) break;
             
             // Null terminate
             l_aux[file.gcount()] = '\0';
             size_read = file.gcount();
             
             if (size_read == 0) break;

             rInfo->numero++;
             if (rInfo->lineas.size() < rInfo->numero) {
                 rInfo->lineas.resize(rInfo->numero);
             }
             
             linea = rInfo->lineas[rInfo->numero - 1];
             if (!linea) {
                 linea = new t_linea_t();
                 rInfo->lineas[rInfo->numero - 1] = linea;
             }
             
             linea->dato = adjustl(l_aux);
             linea->LEN = len_trim(linea->dato);
        }
        
        file.close();
    }

    void carga_raw_info(void* rawFileInfo_ptr, const std::string& filename, const std::string& extension) {
        struct t_NFDE_FILE_t* rawFileInfo = static_cast<struct t_NFDE_FILE_t*>(rawFileInfo_ptr);
        
        struct t_linea_t* linea = nullptr;
        bool ok = false;
        char l_aux[BUFSIZE + 1];
        char buffer[BUFSIZE + 1];
        int i, tamanio, i0, ascii, offset, ascii_menos1, j, k;
        // Character(Len=:), Allocatable :: fichero - Not used in the visible logic, skipping
        
        int UNIT_EF = 10;
        int prelines = 0, io = 0;

        // allocate(rawFileInfo)
        // Assuming rawFileInfo is already allocated or needs allocation. 
        // The pointer is passed in, so we assume it's valid memory.
        
        rawFileInfo->numero = 0;
        rawFileInfo->targ = 1;

        // precount
        std::ifstream file_ef(filename);
        if (!file_ef.is_open()) {
            return;
        }

        while (file_ef.getline(l_aux, BUFSIZE)) {
            io = 0; // Success
            if (file_ef.eof()) {
                io = 1; // End of file
                break;
            }
            prelines++;
        }
        file_ef.close();

        if (prelines == 1 && trim(extension) == ".json") {
            rawFileInfo->numero = countLinesInJSONOneLiner(filename, UNIT_EF);      
            readLinesFromJSONOneLiner(rawFileInfo, filename, UNIT_EF);
        } else { 
            rawFileInfo->numero = prelines;
            readLines(rawFileInfo, filename, UNIT_EF);
        }

        for (k = 1; k <= rawFileInfo->numero; k++) {
            linea = rawFileInfo->lineas[k - 1];
            for (j = 1; j <= linea->LEN; j++) {
                i = j;
                // buscaespa: do while ((ichar(linea%dato(i:i))==32).or.(ichar(linea%dato(i:i))==9))
                while (true) {
                    if (i > linea->LEN) break;
                    
                    char c = linea->dato[i - 1]; // 1-based to 0-based
                    if (ichar(c) != 32 && ichar(c) != 9) {
                        break;
                    }

                    if (i + 1 <= linea->LEN) {
                        char c_next = linea->dato[i]; // i+1 in 1-based is i in 0-based
                        if (ichar(c_next) == 32 || ichar(c_next) == 9) {
                            // linea%dato = trim (adjustl(linea%dato(1:i)))//' '//trim (adjustl(linea%dato(i+2:linea%len)))
                            std::string part1 = trim(adjustl(linea->dato.substr(0, i)));
                            std::string part2 = trim(adjustl(linea->dato.substr(i + 1))); // i+2 in 1-based is i+1 in 0-based
                            linea->dato = part1 + " " + part2;
                            linea->LEN = len_trim(linea->dato);
                        }
                    }
                    i++;
                    if (i > linea->LEN) break;
                }
            }
            // update
            linea->dato = trim(adjustl(linea->dato));
            linea->LEN = len_trim(adjustl(linea->dato));
        }

        return;
    }

    // Note: `end subroutine semba_init` is here, implying semba_init was defined before this chunk.
    // We assume semba_init is already translated.

    solver_t semba_create_solver(semba_fdtd_t& this_obj) {
        return solver_ctor(this_obj.sgg, this_obj.media, this_obj.tag_numbers, 
                           this_obj.SINPML_fullsize, this_obj.fullsize, 
                           this_obj.finishedwithsuccess, this_obj.eps0, this_obj.mu0, 
                           this_obj.tagtype, this_obj.l, this_obj.maxSourceValue, 
                           this_obj.time_desdelanzamiento);
    }

    void semba_update_after_simulation(semba_fdtd_t& this_obj, bool success, SGGFDTDINFO_t sgg, media_matrices_t media, double eps, double mu) {
        this_obj.finishedwithsuccess = success;
        this_obj.sgg = sgg;
        this_obj.eps0 = eps;
        this_obj.mu0 = mu;
        this_obj.media = media;
    }

    void semba_launch(semba_fdtd_t& this_obj) {
        solver_t solver;
        char dubuf[BUFSIZE + 1];
        bool dummylog = false;

        // call each simulation   !ojo que los layoutnumbers empiezan en 0
        if (this_obj.l->finaltimestep != 0) {
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
            this_obj.finishedwithsuccess = false;
            solver = this_obj.create_solver();
#ifdef CompileWithMTLN
            solver.mtln_parsed = this_obj.mtln_parsed;
#endif
            if ((this_obj.l->finaltimestep >= 0) && (!this_obj.l->skindepthpre)) {
                solver.launch_simulation();
                this_obj.update_after_simulation(solver.finishedwithsuccess, solver.sgg, solver.eps0, solver.mu0, solver.media);

                // deallocate(this%media%sggMiEx, ...)
                // Assuming media is a struct with pointers that need deletion
                delete[] this_obj.media.sggMiEx;
                delete[] this_obj.media.sggMiEy;
                delete[] this_obj.media.sggMiEz;
                delete[] this_obj.media.sggMiHx;
                delete[] this_obj.media.sggMiHy;
                delete[] this_obj.media.sggMiHz;
                delete[] this_obj.media.sggMiNo;
                delete[] this_obj.media.sggMtag;
            } else {
#ifdef CompileWithMPI
                MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
                get_secnds(this_obj.l->time_out2);
                if (this_obj.l->layoutnumber == 0) {
                    print_credits(this_obj.l);
                    snprintf(dubuf, BUFSIZE, "BEGUN %s at %s/%s/%s , %s:%s", 
                             trim(adjustl(this_obj.l->nEntradaRoot)).c_str(),
                             this_obj.time_comienzo->fecha.substr(6, 2).c_str(), // 7:8 in 1-based is 6:7 in 0-based
                             this_obj.time_comienzo->fecha.substr(4, 2).c_str(), // 5:6
                             this_obj.time_comienzo->fecha.substr(0, 4).c_str(), // 1:4
                             this_obj.time_comienzo->hora.substr(0, 2).c_str(),  // 1:2
                             this_obj.time_comienzo->hora.substr(2, 2).c_str()); // 3:4
                    print11(this_obj.l->layoutnumber, dubuf);
                    
                    snprintf(dubuf, BUFSIZE, "ENDED %s at %s/%s/%s , %s:%s", 
                             trim(adjustl(this_obj.l->nEntradaRoot)).c_str(),
                             this_obj.l->time_out2->fecha.substr(6, 2).c_str(),
                             this_obj.l->time_out2->fecha.substr(4, 2).c_str(),
                             this_obj.l->time_out2->fecha.substr(0, 4).c_str(),
                             this_obj.l->time_out2->hora.substr(0, 2).c_str(),
                             this_obj.l->time_out2->hora.substr(2, 2).c_str());
                    print11(this_obj.l->layoutnumber, dubuf);
                    
                    snprintf(dubuf, BUFSIZE, "%s %s %s", SEPARADOR, SEPARADOR, SEPARADOR);
                    print11(this_obj.l->layoutnumber, dubuf);
                    print11(this_obj.l->layoutnumber, dubuf);
                }
                CLOSEWARNINGFILE(this_obj.l->layoutnumber, this_obj.l->num_procs, dummylog, this_obj.l->stochastic, this_obj.l->simu_devia);
#ifdef CompileWithMPI
                MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
#ifdef CompileWithMPI
                MPI_FINALIZE(this_obj.l->ierr);
#endif
                exit(0); // stop
            }
        }
        //
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
    }

    void semba_end(semba_fdtd_t& this_obj) {
        char dubuf[BUFSIZE + 1];
        bool existe = false;  
        char filenombre[BUFSIZE + 1] = " ";

        if (this_obj.l->layoutnumber == 0) {
            if (this_obj.l->run) {
                std::ofstream f("running");
                f << "!END" << std::endl;
                f.close();
                std::remove("running");
            }
            snprintf(dubuf, BUFSIZE, "%s %s %s", SEPARADOR, SEPARADOR, SEPARADOR);
            print11(this_obj.l->layoutnumber, dubuf);
            
            snprintf(dubuf, BUFSIZE, "DONE :  %s UNTIL n=%d", 
                     trim(adjustl(this_obj.l->nEntradaRoot)).c_str(), 
                     this_obj.l->finaltimestep);
            print11(this_obj.l->layoutnumber, dubuf);
            
            snprintf(dubuf, BUFSIZE, "%s %s %s", SEPARADOR, SEPARADOR, SEPARADOR);
            print11(this_obj.l->layoutnumber, dubuf);
            print11(this_obj.l->layoutnumber, dubuf);
            erasesignalingfiles(this_obj.l->simu_devia);
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
        //
        if (this_obj.l->deleteintermediates) {
            snprintf(dubuf, BUFSIZE, "%s %s %s", SEPARADOR, SEPARADOR, SEPARADOR);
            print11(this_obj.l->layoutnumber, dubuf);
            print11(this_obj.l->layoutnumber, "Attempting to delete all intermediate data files");
            print11(this_obj.l->layoutnumber, dubuf);
            
            std::string inquire_file = trim(adjustl(this_obj.l->nEntradaRoot)) + "_Outputrequests_" + trim(adjustl(this_obj.whoamishort)) + ".txt";
            std::ifstream inq(inquire_file);
            existe = inq.is_open();
            inq.close();

            if (existe) {
                std::ifstream f19(inquire_file);
                while (f19.getline(filenombre, BUFSIZE)) {
                    std::string filenombre_trimmed = trim(adjustl(filenombre));
                    if (filenombre_trimmed == "!END") {
                        break;
                    } else {
                        std::ofstream f34(filenombre_trimmed);
                        f34 << "!END" << std::endl;
                        f34.close();
                        std::remove(filenombre_trimmed.c_str());
                    }
                }
                f19.close();
                std::remove(inquire_file.c_str());
                
                if (this_obj.l->layoutnumber == 0) {
                    std::string list_file = trim(adjustl(this_obj.l->nEntradaRoot)) + "_Outputlists.dat";
                    std::ofstream f33(list_file);
                    f33 << "!END" << std::endl;
                    f33.close();
                    std::remove(list_file.c_str());
                }
            }
        }
        //
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
        get_secnds(this_obj.l->time_out2);
        if (this_obj.l->layoutnumber == 0) {
            print_credits(this_obj.l);
            snprintf(dubuf, BUFSIZE, "BEGUN %s at %s/%s/%s , %s:%s", 
                     trim(adjustl(this_obj.l->nEntradaRoot)).c_str(),
                     this_obj.time_comienzo->fecha.substr(6, 2).c_str(),
                     this_obj.time_comienzo->fecha.substr(4, 2).c_str(),
                     this_obj.time_comienzo->fecha.substr(0, 4).c_str(),
                     this_obj.time_comienzo->hora.substr(0, 2).c_str(),
                     this_obj.time_comienzo->hora.substr(2, 2).c_str());
            print11(this_obj.l->layoutnumber, dubuf);
            
            snprintf(dubuf, BUFSIZE, "ENDED %s at %s/%s/%s , %s:%s", 
                     trim(adjustl(this_obj.l->nEntradaRoot)).c_str(),
                     this_obj.l->time_out2->fecha.substr(6, 2).c_str(),
                     this_obj.l->time_out2->fecha.substr(4, 2).c_str(),
                     this_obj.l->time_out2->fecha.substr(0, 4).c_str(),
                     this_obj.l->time_out2->hora.substr(0, 2).c_str(),
                     this_obj.l->time_out2->hora.substr(2, 2).c_str());
            print11(this_obj.l->layoutnumber, dubuf);
            
            snprintf(dubuf, BUFSIZE, "%s %s %s", SEPARADOR, SEPARADOR, SEPARADOR);
            print11(this_obj.l->layoutnumber, dubuf);
            print11(this_obj.l->layoutnumber, dubuf);
        }
        
        std::ifstream relaunch_inq("relaunch");
        this_obj.l->relaunching = relaunch_inq.is_open();
        relaunch_inq.close();

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif

#ifdef keeppause
        if (this_obj.l->fatalerror) {
            fatalerror_aux = true;
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
            // MPI_AllReduce is tricky with bools, assuming MPI_LOGICAL maps to int or bool
            int send = fatalerror_aux ? 1 : 0;
            int recv = 0;
            MPI_AllReduce(&send, &recv, 1, MPI_INT, MPI_MAX, SUBCOMM_MPI, this_obj.l->ierr);
            this_obj.l->fatalerror = (recv != 0);
#else
            this_obj.l->fatalerror = fatalerror_aux;
#endif
            if (this_obj.l->fatalerror) this_obj.l->relaunching = true;
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
        }
#endif

        if (this_obj.l->relaunching && (!this_obj.finishedwithsuccess)) {
            if (this_obj.l->layoutnumber == 0) {
                print11(this_obj.l->layoutnumber, std::string(SEPARADOR) + SEPARADOR);
                print11(this_obj.l->layoutnumber, "Not finishing solicited either manually or by an error condition. Edit of create launch file and remove pause file ");
                print11(this_obj.l->layoutnumber, std::string(SEPARADOR) + SEPARADOR);
                
                std::ofstream f9_pause("pause");
                f9_pause << " " << std::endl;
                f9_pause.close();
                
                std::ofstream f9_relaunch("relaunch");
                f9_relaunch << " " << std::endl;
                f9_relaunch.close();
                std::remove("relaunch");
            }
            //!!!!
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif
            if (this_obj.l->layoutnumber == 0) {
                CloseReportingFiles();
            }
            // GO TO 652
        }
        //si ha acabado con exito sal borrando signal files
        if (this_obj.finishedwithsuccess) {
            if (this_obj.l->layoutnumber == 0) {
                std::ofstream f9_pause("pause");
                f9_pause << " " << std::endl;
                f9_pause.close();
                std::remove("pause");
                
                std::ofstream f9_relaunch("relaunch");
                f9_relaunch << " " << std::endl;
                f9_relaunch.close();
                std::remove("relaunch");
                
                std::ofstream f9_running("running");
                f9_running << " " << std::endl;
                f9_running.close();
                std::remove("running");
            }
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, this_obj.l->ierr);
#endif

        if (this_obj.l->layoutnumber == 0) {
            CloseReportingFiles();
        }
        //**************************************************************************************************

#ifdef CompileWithMPI
        MPI_FINALIZE(this_obj.l->ierr);
#endif
    }

    void initEntrada(entrada_t& input) {
        input.geomfile = " ";
        input.prefix = " ";
        input.fichin = " ";
        input.chain2 = " ";
        input.opcionestotales = " "; 
        input.nEntradaRoot = " ";
        input.fileFDE = " ";
        input.fileH5 = " ";
        input.prefixopci = " ";
        input.prefixopci1 = " ";
        input.opcionespararesumeo = " ";
        input.opcionesoriginales = " ";
        input.slicesoriginales = " ";
        input.chdummy = " ";
        input.flushsecondsFields = 0.;
        input.flushsecondsData = 0.;
        input.time_end = 0.; 
        input.existeNFDE = false;
        input.existeh5 = false;
        input.creditosyaprinteados = false;
        input.EpsMuTimeScale_input_parameters.init0();
    }

} // namespace SEMBA_FDTD_m