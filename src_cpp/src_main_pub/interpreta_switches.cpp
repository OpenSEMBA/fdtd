#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <chrono>
#include <cmath>

// Forward declarations and stubs for external modules/types to make the code compile
// In a real translation, these would be implemented in their respective files.

namespace FDETYPES_m {
    constexpr int RKIND = 8; // Assuming double precision
    constexpr int RKIND_wires = 8;
    constexpr int BUFSIZE = 1024;
    constexpr int BUFSIZE_LONG = 4096;
}

namespace Getargs_m {
    std::string getBinaryPath();
    int commandargumentcount(const std::string& input, const std::string& binaryPath);
    void getcommandargument(const std::string& input, int i, std::string& chain, int& length, int& statuse, const std::string& binaryPath);
}

namespace EpsMuTimeScale_m {
    struct EpsMuTimeScale_input_parameters_t {
        // Placeholder fields
    };
}

namespace Report_m {
    void print11(int layoutnumber, const std::string& msg);
    void stoponerror(int layoutnumber, int num_procs, const std::string& msg, bool fatal);
    void print_basic_help(void* l); // Simplified signature
    void print_help(void* l);
    void print_credits(void* l);
    void removeintraspaces(std::string& s);
    void buscaswitchficheroinput(void* l);
    void default_flags(void* l);
    void get_secnds(void* time_out2);
}

namespace version_m {
    int input_conformal_flag = 0;
}

// External MPI stubs
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#else
#define MPI_Barrier(comm, err) do {} while(0)
#endif

// Struct definitions based on Fortran types
struct nf2ff_T {
    // Placeholder
};

struct MedioExtra_t {
    // Placeholder
};

struct tiempo_t {
    double segundos;
};

struct entrada_t {
    // Logicals
    bool forcing = false;
    bool singlefilewrite = false;
    bool ignoresamplingerrors = false;
    bool ignoreerrors = false;
    bool updateshared = false;
    bool prioritizeISOTROPICBODYoverall = false;
    bool wirecrank = false;
    bool CLIPREGION = false;
    bool verbose = false;
    bool resume = false;
    bool forcesteps = false;
    bool resume_fromold = false;
    bool freshstart = false;
    bool run = false;
    bool createmap = false;
    bool dontwritevtk = false;
    bool vtkindex = false;
    bool createmapvtk = false;
    bool run_with_dmma = false;
    bool run_with_abrezanjas = false;
    bool input_conformal_flag = false;
    bool pausar = false;
    bool relaunching = false;
    bool forcestop = false;
    bool l_aux = false;
    bool flag_conf_sgg = false;
    bool takeintcripte = false;
    bool skindepthpre = false;
    bool sgbc = false;
    bool conformalskin = false;
    bool ade = false;
    bool mibc = false;
    bool NOcompomur = false;
    bool MurAfterPML = false;
    bool sgbccrank = false;
    bool sgbcDispersive = false;
    bool saveall = false;
    bool boundwireradius = false;
    bool hay_slanted_wires = false;
    bool makeholes = false;
    bool mur_first = false;
    bool mur_second = false;
    bool mur_exist = false;
    bool connectendings = false;
    bool strictOLD = false;
    bool mtlnberenger = false;
    bool stableradholland = false;
    bool TAPARRABOS = false;
    bool fieldtotl = false;
    bool forceresampled = false;
    bool isolategroupgroups = false;
    bool groundwires = false;
    bool noSlantedcrecepelo = false;
    bool forcecfl = false;
    bool niapapostprocess = false;
    bool permitscaling = false;
    bool stochastic = false;
    bool chosenyesornostochastic = false;
    bool prioritizeCOMPOoverPEC = false;
    bool prioritizeTHINWIRE = false;
    bool createh5bin = false;
    bool deleteintermediates = false;
    bool existeNFDE = false;
    bool file11isopen = false;
    bool NF2FFDecim = false;
    bool existeh5 = false;
    bool fatalerror = false;
    bool fatalerrornfde2sgg = false;
    bool thereare_stoch = false;
    bool experimentalVideal = false;
    bool simu_devia = false;
    bool noconformalmapvtk = false;
    bool createh5filefromsinglebin = false;
    bool creditosyaprinteados = false;
    bool read_command_line = false;

    // Integers
    int wirethickness = 0;
    int inductance_order = 0;
    int finaltimestep = 0;
    int ierr = 0;
    int layoutnumber = 0;
    int num_procs = 0;
    int length = 0;
    int mpidir = 0;
    int flushminutesFields = 0;
    int flushminutesData = 0;
    int flushsecondsFields = 0;
    int flushsecondsData = 0;
    int forced = 0;
    int maxCPUtime = 0;
    int sgbcdepth = 0;
    int precision = 0;
    int statuse = 0;

    // Reals (RKIND usually 8 for double in modern Fortran codes, mapped to double)
    double maxwireradius = 0.0;
    double mindistwires = 0.0;
    double attfactorc = 0.0;
    double attfactorw = 0.0;
    double cfltemp = 0.0;
    double cfl = 0.0;
    double sgbcfreq = 0.0;
    double sgbcresol = 0.0;
    double alphamaxpar = 0.0;
    double kappamaxpar = 0.0;
    double alphaOrden = 0.0;

    // Reals (kind=8)
    double time_begin = 0.0;
    double time_end = 0.0;

    // Reals (RKIND_wires)
    double factorradius = 0.0;
    double factordelta = 0.0;

    // Types
    nf2ff_T facesNF2FF;
    MedioExtra_t MEDIOEXTRA;
    EpsMuTimeScale_m::EpsMuTimeScale_input_parameters_t EpsMuTimeScale_input_parameters;
    tiempo_t time_out2;

    // Characters
    std::string prefix;
    std::string prefixopci;
    std::string prefixopci1;
    std::string opcionespararesumeo;
    std::string opcionesoriginales;
    std::string slicesoriginales;
    std::string chdummy;
    std::string chaininput;
    std::string chain;
    std::string chain2;
    std::string fichin;
    std::string extension;
    std::string nresumeable2;
    std::string fileFDE;
    std::string fileH5;
    std::string inductance_model;
    std::string wiresflavor;
    std::string nEntradaRoot;
    std::string opcionestotales;
    std::string conformal_file_input_name;
    std::string geomfile;
};

namespace interpreta_switches_m {

    void interpreta(entrada_t& l, int& statuse) {
        std::string chari, f, dubuf, buff, binaryPath;
        bool existiarunningigual = false;
        bool mpidirset = false;
        bool resume3 = false;
        int i, j, donde, n, newmpidir;
        double pausetime;
        int iostatus = 0;

        l.input_conformal_flag = (version_m::input_conformal_flag != 0);

        mpidirset = false;
        existiarunningigual = false;
        statuse = 0;

        binaryPath = Getargs_m::getBinaryPath();
        n = Getargs_m::commandargumentcount(l.chaininput, binaryPath);
        if (n == 0) {
            Report_m::print_basic_help(&l);
            Report_m::stoponerror(l.layoutnumber, l.num_procs, "Error: NO arguments neither command line nor in launch file. Correct and remove pause...", true);
            statuse = -1;
            return;
        }
        l.opcionestotales = "";
        for (i = 1; i <= n; ++i) {
            Getargs_m::getcommandargument(l.chaininput, i, l.chain, l.length, statuse, binaryPath);
            if (statuse != 0) {
                Report_m::stoponerror(l.layoutnumber, l.num_procs, "Reading input", true);
                statuse = -1;
                return;
            }
            // trim and adjustl equivalent
            std::string trimmed_opts = l.opcionestotales;
            trimmed_opts.erase(0, trimmed_opts.find_first_not_of(" \t"));
            trimmed_opts.erase(trimmed_opts.find_last_not_of(" \t") + 1);
            
            std::string trimmed_chain = l.chain;
            trimmed_chain.erase(0, trimmed_chain.find_first_not_of(" \t"));
            trimmed_chain.erase(trimmed_chain.find_last_not_of(" \t") + 1);

            l.opcionestotales = trimmed_opts + " " + trimmed_chain;
        }
        Report_m::print11(l.layoutnumber, "Switches " + l.opcionestotales);

        if (n > 0) {
            i = 2; // Start at 2 because first arg is executable name
            while (i <= n) {
                Getargs_m::getcommandargument(l.chaininput, i, l.chain, l.length, statuse, binaryPath);
                if (statuse != 0) {
                    Report_m::stoponerror(l.layoutnumber, l.num_procs, "Reading input", true);
                    statuse = -1;
                    return;
                }
                
                std::string chain_trimmed = l.chain;
                chain_trimmed.erase(0, chain_trimmed.find_first_not_of(" \t"));
                chain_trimmed.erase(chain_trimmed.find_last_not_of(" \t") + 1);

                if (chain_trimmed == "-i") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    continue; // already interpreted
                } else if (chain_trimmed == "-a") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    continue; // already interpreted
                } else if (chain_trimmed == "-mpidir") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    
                    std::string f_trimmed = f;
                    f_trimmed.erase(0, f_trimmed.find_first_not_of(" \t"));
                    f_trimmed.erase(f_trimmed.find_last_not_of(" \t") + 1);

                    if (f_trimmed == "x" || f_trimmed == "X") {
                        l.mpidir = 1;
                    } else if (f_trimmed == "y" || f_trimmed == "Y") {
                        l.mpidir = 2;
                    } else if (f_trimmed == "z" || f_trimmed == "Z") {
                        l.mpidir = 3;
                    } else {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid or duplicate incoherent -mpidir option", true);
                        statuse = -1;
                        continue;
                    }
                    if (!mpidirset) {
                        std::string chain_trimmed2 = l.chain;
                        chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t"));
                        chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t") + 1);
                        
                        std::string f_trimmed2 = f;
                        f_trimmed2.erase(0, f_trimmed2.find_first_not_of(" \t"));
                        f_trimmed2.erase(f_trimmed2.find_last_not_of(" \t") + 1);

                        std::string opts_trimmed = l.opcionespararesumeo;
                        opts_trimmed.erase(0, opts_trimmed.find_first_not_of(" \t"));
                        opts_trimmed.erase(opts_trimmed.find_last_not_of(" \t") + 1);

                        l.opcionespararesumeo = opts_trimmed + " " + chain_trimmed2 + " " + f_trimmed2;
                        mpidirset = true;
                    }
                } else if (chain_trimmed == "-pause") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    
                    std::istringstream iss(f);
                    iss >> pausetime;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid pause time", true);
                    }
                    if (pausetime <= 0) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid pause time: zero or negative value", true);
                        statuse = -1;
                    }
                    
                    l.pausar = true;
#ifdef CompileWithMPI
                    MPI_Barrier(SUBCOMM_MPI, l.ierr);
#endif
                    Report_m::get_secnds(&l.time_out2);
                    l.time_begin = l.time_out2.segundos;
                    dubuf = "Paused for (secs) " + std::to_string(pausetime);
                    Report_m::print11(l.layoutnumber, dubuf);
                    
                    while (l.pausar) {
#ifdef CompileWithMPI
                        // Placeholder for pause logic if needed, currently empty in snippet
#endif
                    }
                }
                
                i = i + 1;
            }
        }
    }

    void insertalogtmp(entrada_t& l) {
        // Stub
    }

    void print_help(entrada_t& l) {
        // Stub
    }

    void print_basic_help(entrada_t& l) {
        // Stub
    }

    void print_credits(entrada_t& l) {
        // Stub
    }

    void removeintraspaces(std::string& s) {
        // Stub
    }

    void buscaswitchficheroinput(entrada_t& l) {
        // Stub
    }

    void default_flags(entrada_t& l) {
        // Stub
    }
}

MPI_Barrier(SUBCOMM_MPI, l->ierr);
#endif
                  get_secnds(l->time_out2);
                  l->time_end = l->time_out2.segundos;
                  if (l->time_end - l->time_begin > pausetime) {
                     l->pausar = false;
                  }
               }
#ifdef CompileWithMPI
               MPI_Barrier(SUBCOMM_MPI, l->ierr);
               l->l_aux = l->pausar;
               MPI_AllReduce(&l->l_aux, &l->pausar, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l->ierr);
#endif
            case '-NF2FFDecim':
               l->NF2FFDecim = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain)) + " " + trim(adjustl(f));
               break;
            case '-noNF2FF':
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               if (trim(adjustl(f)) == "back" || trim(adjustl(f)) == "BACK") {
                  l->facesNF2FF.TR = false;
               } else if (trim(adjustl(f)) == "front" || trim(adjustl(f)) == "FRONT") {
                  l->facesNF2FF.FR = false;
               } else if (trim(adjustl(f)) == "left" || trim(adjustl(f)) == "LEFT") {
                  l->facesNF2FF.IZ = false;
               } else if (trim(adjustl(f)) == "right" || trim(adjustl(f)) == "RIGHT") {
                  l->facesNF2FF.DE = false;
               } else if (trim(adjustl(f)) == "down" || trim(adjustl(f)) == "DOWN") {
                  l->facesNF2FF.AB = false;
               } else if (trim(adjustl(f)) == "up" || trim(adjustl(f)) == "UP") {
                  l->facesNF2FF.AR = false;
               } else {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid -noNF2FF option", true);
                  statuse = -1;
               }
               continue;
               // COMO LA RCS SE CALCULA SOLO AL FINAL NO OBLIGO A RESUMEAR CON IGUAL -NONFF2FF PARA PODER CALCULAR CON Y SIN ESTA OPCION resumeando
               // l->opcionespararesumeo = trim (adjustl(l->opcionespararesumeo)) + " " + trim (adjustl(l->chain)) + " " + trim (adjustl(f));
            case '-force':
               l->forcing = true;
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               try {
                  std::istringstream iss(f);
                  iss >> l->forced;
                  if (iss.fail()) throw std::runtime_error("Invalid cut");
               } catch (...) {
                  goto label_412;
               }
               goto label_312;
            label_412:
               stoponerror(l->layoutnumber, l->num_procs, "Invalid cut", true);
               statuse = -1;
            label_312:
               continue;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain)) + " " + trim(adjustl(f));
               break;
            case '-singlefile':
               l->singlefilewrite = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-ignoresamplingerrors':
               l->ignoresamplingerrors = true;
               break;
            case '-prioritizeTHINWIRE':
               l->prioritizeTHINWIRE = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               l->ignoreerrors = true;
               break;
            case '-prioritizeCOMPOoverPEC':
               l->prioritizeCOMPOoverPEC = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               l->ignoreerrors = true;
               break;
            case '-noshared':
               l->updateshared = false;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-prioritizeISOTROPICBODYoverall':
               l->prioritizeISOTROPICBODYoverall = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-wirecrank':
               l->wirecrank = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-clip':
               l->CLIPREGION = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
               // !!!!#endif
               // !

            case '-verbose':
               l->verbose = true;
               break;
            case '-ignoreerrors':
               l->ignoreerrors = true;
               break;
            case '-r':
               l->resume = true;
               l->forcesteps = true;
#ifdef CompileWithOldSaving
            case '-old':
               l->resume_fromold = true;
               break;
#endif
            case '-cpumax':
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               try {
                  std::istringstream iss(f);
                  iss >> l->maxCPUtime;
                  if (iss.fail()) throw std::runtime_error("Invalid CPU maximum time");
               } catch (...) {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid CPU maximum time", true);
               }
               if (l->maxCPUtime <= 0) {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid CPU maximum time", true);
                  statuse = -1;
               }
               break;
            case '-s':
               l->freshstart = true;
               break;
            case '-flush':
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               try {
                  std::istringstream iss(f);
                  iss >> l->flushminutesFields;
                  if (iss.fail()) throw std::runtime_error("Invalid flushing interval");
               } catch (...) {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid flushing interval", true);
               }
               if (l->flushminutesFields <= 0) {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid flushing interval", true);
                  statuse = -1;
               }
               break;
            case '-flushdata':
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               try {
                  std::istringstream iss(f);
                  iss >> l->flushminutesData;
                  if (iss.fail()) throw std::runtime_error("Invalid flushing interval");
               } catch (...) {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid flushing interval", true);
               }
            label_401:
               if (l->flushminutesData <= 0) {
                  stoponerror(l->layoutnumber, l->num_procs, "Invalid flushing interval", true);
                  statuse = -1;
               }
               break;
            case '-run':
               l->run = true;
               break;
            case '-map':
               l->createmap = true;
               break;
            case '-dontwritevtk':
               l->dontwritevtk = true;
               break;
            case '-vtkindex':
               l->vtkindex = true;
               break;
            case '-mapvtk':
               l->createmapvtk = true;
               break;
            case '-dmma':
               l->run_with_dmma = true;
               l->run_with_abrezanjas = false;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-takeintcripte':
               l->takeintcripte = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
#ifdef CompileWithNIBC
            case '-skindepthpre':
               l->skindepthpre = true;
               // l->opcionespararesumeo = trim (adjustl(l->opcionespararesumeo)) + " " + trim (adjustl(l->chain));
               break;
            case '-mibc':
            case '-skindepth':
               l->mibc = true;
               l->sgbc = false;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-conformalskin':
               l->conformalskin = true;
               l->mibc = true;
               l->sgbc = false;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-ade':
               l->ade = true;
               l->mibc = true;
               l->sgbc = false;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-NOcompomur':
               l->NOcompomur = true;
               l->mibc = true;
               l->sgbc = false;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
#endif
            case '-mur2':
               l->MurAfterPML = true;
               // l->mur_second=true;
               l->mur_first = true;
               // arreglar cuando resuelva el bug en mur segundo orden
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
               break;
            case '-mur1':
               l->MurAfterPML = true;
               l->mur_first = true;
               break;

l.mur_second = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case '-pmlalpha':
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.alphamaxpar;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
                  statuse = -1;
                  goto label_8621;
               }
               goto label_8621;
            label_7621:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
               statuse = -1;
            label_8621:
               if (l.alphamaxpar < 0.0) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
                  statuse = -1;
               }
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.alphaOrden;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML order factor", true);
                  statuse = -1;
                  goto label_8121;
               }
               goto label_8121;
            label_7121:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML order factor", true);
               statuse = -1;
            label_8121:
               if (l.alphaOrden < 0.0) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
                  statuse = -1;
               }
               break;
            case '-pmlkappa':
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.kappamaxpar;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML kappa factor", true);
                  statuse = -1;
                  goto label_8622;
               }
               goto label_8622;
            label_7622:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML kappa factor", true);
               statuse = -1;
            label_8622:
               if (l.kappamaxpar < 1.0) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML kappa factor", true);
                  statuse = -1;
               }
               break;
            case '-pmlcorr':
               l.MEDIOEXTRA.exists = true;
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.MEDIOEXTRA.sigma;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid pmlcorr sigma factor", true);
                  statuse = -1;
                  goto label_8672;
               }
               goto label_8672;
            label_7672:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid pmlcorr sigma factor", true);
               statuse = -1;
            label_8672:
               if (l.MEDIOEXTRA.sigma < 0.0) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid pmlcorr sigma factor", true);
                  statuse = -1;
               }
               l.MEDIOEXTRA.sigmam = -1.0; // voids it. later overriden
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.MEDIOEXTRA.pml_size;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid pmlcorr depth factor", true);
                  statuse = -1;
                  goto label_8662;
               }
               goto label_8662;
            label_7662:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid pmlcorr depth factor", true);
               statuse = -1;
            label_8662:
               if (l.MEDIOEXTRA.pml_size < 0) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid pmlcorr depth factor", true);
                  statuse = -1;
               }
               break;
            case '-attc':
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.attfactorc;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
                  statuse = -1;
                  goto label_866;
               }
               goto label_866;
            label_766:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
               statuse = -1;
            label_866:
               if ((l.attfactorc <= -1.0) || (l.attfactorc > 1.0)) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
                  statuse = -1;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case '-sgbcdepth':
               l.mibc = false;
               l.sgbc = true;
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.sgbcdepth;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc depth ", true);
                  statuse = -1;
                  goto label_8466;
               }
               goto label_8466;
            label_7466:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc depth ", true);
               statuse = -1;
            label_8466:
               if (l.sgbcdepth < -1) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc depth", true);
                  statuse = -1;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case '-sgbcfreq':
               l.sgbc = true;
               l.mibc = false;
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.sgbcfreq;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc freq ", true);
                  statuse = -1;
                  goto label_84616;
               }
               goto label_84616;
            label_74616:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc freq ", true);
               statuse = -1;
            label_84616:
               if (l.sgbcfreq < 0.) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc freq", true);
                  statuse = -1;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case '-sgbcresol':
               l.mibc = false;
               l.sgbc = true;
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                  std::istringstream iss(f);
                  iss >> l.sgbcresol;
                  if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc decay ", true);
                  statuse = -1;
                  goto label_84626;
               }
               goto label_84626;
            label_74626:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc decay ", true);
               statuse = -1;
            label_84626:
               if (l.sgbcresol < 0.0) {
                  stoponerror(l.layoutnumber, l.num_procs, "Invalid sgbc decay", true);
                  statuse = -1;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case '-sgbcyee':
               l.sgbc = true;
               l.mibc = false;
               l.sgbccrank = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case '-sgbccrank': // es el default. Lo mantengo por compatibilidad con lanzamientos previos
               l.sgbccrank = true;
               l.mibc = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case '-nosgbc': // opcion generica que aglutina varios switches que estan den default (l%sgbcresol, l%sgbccrank, l%sgbcfreq)
               l.sgbc = false;
               l.mibc = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;

case ("-sgbc"): // opcion generica que aglutina varios switches que estan den default (l.sgbcresol, l.sgbccrank, l.sgbcfreq)
               l.sgbc = true;
               l.mibc = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-sgbcDispersive"): // opcion generica que aglutina varios switches que estan den default (l.sgbcresol, l.sgbccrank, l.sgbcfreq)
               l.sgbc = true;
               l.mibc = false;
               l.sgbcDispersive = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-saveall"):
               l.saveall = true;
               break;
            case ("-attw"):
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                   std::istringstream iss(f);
                   iss >> l.attfactorw;
                   if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
                   statuse = -1;
                   goto label_668;
               }
               goto label_832;
            label_732:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
               statuse = -1;
               goto label_668;
            label_832:
               if ((l.attfactorw <= -1.0) || (l.attfactorw > 1.0)) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
                   statuse = -1;
                   goto label_668;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case ("-maxwireradius"):
               l.boundwireradius = true;
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                   std::istringstream iss(f);
                   iss >> l.maxwireradius;
                   if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
                   statuse = -1;
                   goto label_668;
               }
               goto label_837;
            label_737:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid dissipation factor", true);
               statuse = -1;
               goto label_668;
            label_837:
               if ((l.maxwireradius <= 0.0)) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid maximumwireradius", true);
                   statuse = -1;
                   goto label_668;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case ("-mindistwires"):
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                   std::istringstream iss(f);
                   iss >> l.mindistwires;
                   if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid minimum distance between wires", true);
                   statuse = -1;
                   goto label_668;
               }
               goto label_1832;
            label_1732:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid minimum distance between wires", true);
               statuse = -1;
               goto label_668;
            label_1832:
               if (l.mindistwires <= 0.0) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid minimum distance between wires", true);
                   statuse = -1;
                   goto label_668;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case ("-makeholes"):
               l.makeholes = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-connectendings"):
               l.connectendings = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-nostrictOLD"):
               l.strictOLD = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-nomtlnberenger"):
               l.mtlnberenger = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-stableradholland"):
               l.stableradholland = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            //  CASE ('-mtln')
            //      buff='-mtln option deprecated and ignored. Check -nomtlnberenger or -l%stableradholland'
            //      call WarnErrReport(Trim(buff),.false.)
            case ("-intrawiresimplify"):
               l.strictOLD = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-notaparrabos"):
               l.TAPARRABOS = false;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
            case ("-fieldtotl"):
               l.fieldtotl = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;
          //case ('-experimentalVideal')
          //    l.experimentalVideal=.true.
          //    l.opcionespararesumeo = trim (adjustl(l.opcionespararesumeo)) // ' ' // trim (adjustl(l.chain))

            case ("-forceresampled"): //a menos que se pida explicitamente, no se resamplea 120123
               l.forceresampled = true;
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain));
               break;

            case ("-wirethickness"):
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               // Converts the characters to real
               try {
                   std::istringstream iss(f);
                   iss >> l.wirethickness;
                   if (iss.fail()) throw std::runtime_error("Read error");
               } catch (...) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid l.wirethickness ", true);
                   statuse = -1;
                   goto label_668;
               }
               goto label_8416;
            label_7416:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid l.wirethickness ", true);
               statuse = -1;
               goto label_668;
            label_8416:
               if (l.sgbcdepth < -1) {
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid l.wirethickness", true);
                   statuse = -1;
                   goto label_668;
               }
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               break;
            case ("-wiresflavor"):
               i = i + 1;
               getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
               l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(l.chain)) + " " + trim(adjustl(f));
               try {
                   std::istringstream iss(f);
                   iss >> std::noskipws;
                   char c;
                   iss >> c;
                   // Simulating READ (f, '(a)') which reads a string
                   // In C++, we just read the whole string f again or use the variable f directly if it was formatted
                   // The Fortran READ (f, '(a)') reads the whole line/string into l.wiresflavor
                   // Assuming f contains the string to be read into l.wiresflavor
                   l.wiresflavor = f; 
               } catch (...) {
                   goto label_3621;
               }
               if (trim(adjustl(l.wiresflavor.substr(0, 1))) == "g") l.wiresflavor = "slanted";
               std::string flv_trimmed = trim(adjustl(l.wiresflavor));
               if (flv_trimmed == "holland" || flv_trimmed == "old") {
                   l.wiresflavor = "holland";
               } else if (flv_trimmed == "berenger" || flv_trimmed == "new") {
                   l.wiresflavor = "berenger";
               } else if (flv_trimmed == "slanted" || flv_trimmed == "experimental") {
                   l.wiresflavor = "slanted";
               } else if (flv_trimmed == "transition") {
                   l.wiresflavor = "transition";
               } else if (flv_trimmed == "semistructured") {
                   l.wiresflavor = "semistructured";
                   //
                   i = i + 1;
                   getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                   l.opcionespararesumeo = trim(adjustl(l.opcionespararesumeo)) + " " + trim(adjustl(f));
                   // Converts the characters to real
                   try {
                       std::istringstream iss(f);
                       iss >> l.precision;
                       if (iss.fail()) throw std::runtime_error("Read error");
                   } catch (...) {
                       goto label_2561;
                   }
                   goto label_2562;
                label_2561:
                   stoponerror(l.layoutnumber, l.num_procs, "Invalid l.precision for semistructured", true);
                   statuse = -1;
                   goto label_668;
                label_2562:
                   if (l.precision < 0) {
                       stoponerror(l.layoutnumber, l.num_procs, "Invalid l.precision for semistructured", true);
                       statuse = -1;
                       goto label_668;
                   }
                   //
               }
               goto label_4621;
            label_3621:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid wires flavor", true);
               statuse = -1;
               goto label_668;
            label_4621:
               if (((trim(adjustl(l.wiresflavor)) != "holland") &&
                    (trim(adjustl(l.wiresflavor)) != "transition") &&
                    (trim(adjustl(l.wiresflavor)) != "berenger") &&
                    (trim(adjustl(l.wiresflavor)) != "slanted") &&
                    (trim(adjustl(l.wiresflavor)) != "semistructured")) ||
                   !((trim(adjustl(l.wiresflavor)) == "holland") ^

((l.wiresflavor == "transition") ^ (l.wiresflavor == "berenger") ^ (l.wiresflavor == "slanted") ^ (l.wiresflavor == "semistructured"))) {
                stoponerror(l.layoutnumber, l.num_procs, "Invalid wires flavor->" + l.wiresflavor, true);
                statuse = -1;
            }
#ifndef CompileWithThickWires
            if (l.wiresflavor == "holland" || l.wiresflavor == "transition") {
                if (l.wirethickness != 1) {
                    stoponerror(l.layoutnumber, l.num_procs, "Holland wire flavor not available in this compilation", true);
                    statuse = -1;
                }
            }
#endif
#ifndef CompileWithThickWires
            if (l.wiresflavor == "holland") {
                if (l.wirethickness != 1) {
                    stoponerror(l.layoutnumber, l.num_procs, "Holland wire flavor thickness>1 requires recompiling", true);
                    statuse = -1;
                }
            }
#endif
            if (l.wiresflavor == "berenger" || l.wiresflavor == "slanted" || l.wiresflavor == "experimental" || l.wiresflavor == "transition") {
                if (l.wirethickness != 1) {
                    stoponerror(l.layoutnumber, l.num_procs, "Thickness>1 unsupported for this wireflavor", true);
                    statuse = -1;
                }
            }
#ifndef CompileWithBerengerWires
            if (l.wiresflavor == "berenger") {
                stoponerror(l.layoutnumber, l.num_procs, "Berenger wire flavor not available in this compilation", true);
                statuse = -1;
            }
#endif
#ifndef CompileWithSlantedWires
            if (l.wiresflavor == "slanted" || l.wiresflavor == "experimental") {
                stoponerror(l.layoutnumber, l.num_procs, "Experimental wire flavor not available in this compilation", true);
                statuse = -1;
            }
#endif
            case -isolategroupgroups:
                l.isolategroupgroups = true;
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain;
                break;
            case -groundwires:
                l.groundwires = true;
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain;
                break;
            case -noSlantedcrecepelo:
                l.noSlantedcrecepelo = true;
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain;
                break;
            case -inductance:
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                try {
                    std::istringstream iss(f);
                    iss >> l.inductance_model;
                } catch (...) {
                    goto label_361;
                }
                goto label_461;
            label_361:
                stoponerror(l.layoutnumber, l.num_procs, "Invalid inductance model", true);
                statuse = -1;
                goto label_668;
            label_461:
                if (l.inductance_model != "ledfelt" && l.inductance_model != "berenger" && l.inductance_model != "boutayeb") {
                    stoponerror(l.layoutnumber, l.num_procs, "Invalid inductance model", true);
                    statuse = -1;
                    goto label_668;
                }
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain + " " + f;
                break;
            case -inductanceorder:
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                try {
                    std::istringstream iss(f);
                    iss >> l.inductance_order;
                } catch (...) {
                    goto label_179;
                }
                goto label_180;
            label_179:
                stoponerror(l.layoutnumber, l.num_procs, "Invalid inductance order", true);
                statuse = -1;
                goto label_668;
            label_180:
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain + " " + f;
                break;
            case -prefix:
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                l.prefix = "_" + f;
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain + " " + f;
                break;
            case -cfl:
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                try {
                    std::istringstream iss(f);
                    iss >> l.cfltemp;
                } catch (...) {
                    goto label_3762;
                }
                goto label_3862;
            label_3762:
                stoponerror(l.layoutnumber, l.num_procs, "Invalid Courant Number", true);
                statuse = -1;
                goto label_668;
            label_3862:
                if (l.cfltemp <= 0.0) {
                    print11(l.layoutnumber, "------> Ignoring negative or null l%cfl Courant Number");
                    l.forcecfl = false;
                } else {
                    l.cfl = l.cfltemp;
                    l.forcecfl = true;
                    l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain + " " + f;
                }
                break;
            case -noconformalmapvtk:
                l.noconformalmapvtk = true;
                break;
            case -niapapostprocess:
                l.niapapostprocess = true;
                break;
#ifdef CompileWithPrescale
            case -pscale:
                l.permitscaling = true;
                l.saveall = true;
                i = i + 1;
                buff = "";
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                try {
                    std::istringstream iss(f);
                    iss >> buff;
                } catch (...) {
                    goto label_33762;
                }
                l.EpsMuTimeScale_input_parameters.electric = false;
                l.EpsMuTimeScale_input_parameters.electric = false;
                if (buff == "ee") {
                    l.EpsMuTimeScale_input_parameters.electric = true;
                } else if (buff == "hh") {
                    l.EpsMuTimeScale_input_parameters.magnetic = true;
                } else if (buff == "eh" || buff == "he") {
                    l.EpsMuTimeScale_input_parameters.electric = true;
                    l.EpsMuTimeScale_input_parameters.magnetic = true;
                } else {
                    goto label_33862;
                }
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                l.opcionespararesumeo = l.opcionespararesumeo + " " + l.chain + " " + f;
                try {
                    std::istringstream iss(f);
                    iss >> l.EpsMuTimeScale_input_parameters.tini;
                } catch (...) {
                    goto label_33762;
                }
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                l.opcionespararesumeo = l.opcionespararesumeo + " " + f;
                try {
                    std::istringstream iss(f);
                    iss >> l.EpsMuTimeScale_input_parameters.tend;
                } catch (...) {
                    goto label_33762;
                }
                i = i + 1;
                getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                l.opcionespararesumeo = l.opcionespararesumeo + " " + f;
                try {
                    std::istringstream iss(f);
                    iss >> l.EpsMuTimeScale_input_parameters.alpha_max;
                } catch (...) {
                    goto label_33762;
                }
                goto label_33862;
            label_33762:
                stoponerror(l.layoutnumber, l.num_procs, "Invalid pscale parameters", true);
                statuse = -1;
                goto label_668;
            label_33862:
                if (l.EpsMuTimeScale_input_parameters.checkError() != 0) {
                    stoponerror(l.layoutnumber, l.num_procs, "Invalid -pscale parameters: some parameters have to be greater than 0.0: -pscale t0(>=0) tend slope(>0)", true);
                    statuse = -1;
                    goto label_668;
                }
                break;
#endif

true); statuse = -1; // goto 668
               } else {
                  l->EpsMuTimeScale_input_parameters->are_there = true;
               }
#endif
            case ('-n'):
               l->forcesteps = true;
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               // Converts the characters to integer
               try {
                   std::istringstream iss(f);
                   iss >> l->finaltimestep;
                   if (iss.fail()) throw std::runtime_error("Invalid time step");
               } catch (...) {
                   goto label_602;
               }
               goto label_702;
label_602:
               stoponerror(l->layoutnumber, l->num_procs, "Invalid time step", true); statuse = -1; // goto 668
label_702:
               if (l->finaltimestep < -2) {
                   stoponerror(l->layoutnumber, l->num_procs, "Invalid time step", true); statuse = -1; // goto 668
               }
               // !!!!
            case ('-factorradius'):
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               // Converts the characters to integer
               try {
                   std::istringstream iss(f);
                   iss >> l->factorradius;
                   if (iss.fail()) throw std::runtime_error("Invalid l->factorradius");
               } catch (...) {
                   goto label_6032;
               }
               goto label_7032;
label_6032:
               stoponerror(l->layoutnumber, l->num_procs, "Invalid l->factorradius", true); statuse = -1; // goto 668
label_7032:
               continue;
            case ('-factordelta'):
               i = i + 1;
               getcommandargument(l->chaininput, i, f, l->length, statuse, binaryPath);
               // Converts the characters to integer
               try {
                   std::istringstream iss(f);
                   iss >> l->factordelta;
                   if (iss.fail()) throw std::runtime_error("Invalid l->factordelta");
               } catch (...) {
                   goto label_6072;
               }
               goto label_7072;
label_6072:
               stoponerror(l->layoutnumber, l->num_procs, "Invalid l->factordelta", true); statuse = -1; // goto 668
label_7072:
               continue;
               // !!!!!!!!!
            case ('-stoch'):
               l->stochastic = true;
               l->chosenyesornostochastic = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
#ifndef CompileWithMPI
               stoponerror(l->layoutnumber, l->num_procs, "l->stochastic simulation unsupported without MPI compilation", true); statuse = -1; // goto 668
#endif
            case ('-nostoch'):
               l->stochastic = false;
               l->chosenyesornostochastic = true;
               l->opcionespararesumeo = trim(adjustl(l->opcionespararesumeo)) + " " + trim(adjustl(l->chain));
            case ('-forcecreateh5bin'):
               l->createh5bin = true;
            case (''): // 100615 para evitar el crlf del .sh
               continue;
            default:
               stoponerror(l->layoutnumber, l->num_procs, "Wrong switch " + trim(adjustl(l->chain)), true); statuse = -1; // goto 668
               break;
            }
            i = i + 1;
         }
      }

      if (l->connectendings && l->strictOLD) {
         stoponerror(l->layoutnumber, l->num_procs, "l->strictOLD option not compatible with -l->connectendings", true); statuse = -1; // goto 668
      }
      if (l->TAPARRABOS && (!l->strictOLD)) {
         stoponerror(l->layoutnumber, l->num_procs, "-nostrictOLD option requires -notaparrabos ", true); statuse = -1; // goto 668
      }
      if (l->isolategroupgroups && l->strictOLD) {
         stoponerror(l->layoutnumber, l->num_procs, "-intrawiresimplify option not compatible with -l->isolategroupgroups", true); statuse = -1; // goto 668
      }

      if ((l->sgbc && l->mibc)) {
         stoponerror(l->layoutnumber, l->num_procs, "Use only one of -sgbc -l->mibc", true); statuse = -1; // goto 668
      }
      if (l->freshstart && l->resume) {
         stoponerror(l->layoutnumber, l->num_procs, "Fresh Start option -s not compatible with restarting -r", true); statuse = -1; // goto 668
      }
      if (l->freshstart && l->resume_fromold) {
         stoponerror(l->layoutnumber, l->num_procs, "Fresh Start option -s not compatible with -old", true); statuse = -1; // goto 668
      }
      if ((!l->resume) && (!l->run) && l->resume_fromold) {
         stoponerror(l->layoutnumber, l->num_procs, "l->resume option -r must be used if issuing -old", true); statuse = -1; // goto 668
      }
      if ((l->flushminutesFields != 0) && (l->deleteintermediates)) {
         stoponerror(l->layoutnumber, l->num_procs, "-delete is not compatible with -flush", true); statuse = -1; // goto 668
      }
      if (l->run_with_abrezanjas && l->run_with_dmma) {
         stoponerror(l->layoutnumber, l->num_procs, "-abrezanjas is not compatible with -dmma", true); statuse = -1; // goto 668
      }
      if (l->stochastic && (trim(adjustl(l->wiresflavor)) != "holland")) {
         stoponerror(l->layoutnumber, l->num_procs, "Old wires flavor is the only supported with l->stochastic", true); statuse = -1; // goto 668
      }
      if (l->stochastic && l->wirecrank) {
         stoponerror(l->layoutnumber, l->num_procs, "Wires Crank Nicolson is unsupported with l->stochastic", true); statuse = -1; // goto 668
      }
      // !!!si esta soportado 170719
      // !! if (l->permitscaling.and.l->resume) then
      // !!   call stoponerror (l->layoutnumber, l->num_procs, 'Resuming with Permittivity scaling unsupported',.true.); statuse=-1; !goto 668
      // !!end if
      if (l->permitscaling && (l->kappamaxpar > 1.000001_rkind)) {
         // !!!061118 no lo permito porque cpml toca los idxe, idye, idze en funcion del kappa y permittivity scaling conflicta
         stoponerror(l->layoutnumber, l->num_procs, "Unsupported CPML kappa factor since 061118 because conflicts with Idxe...in permittivity scaling", true);
      }
      if (l->stochastic) {
#ifndef CompileWithStochastic
         StopOnError(l->layoutnumber, l->num_procs, "l->stochastic without compilation support. Recompile");
#endif
#ifdef CompileWithStochastic
#ifndef CompileWithMPI
         StopOnError(l->layoutnumber, l->num_procs, "l->stochastic unsupported without MPI compilation. Recompile");
#endif
#endif
         continue;
      }

      // !!!

      //
      l->prefixopci1 = trim(adjustl(l->opcionespararesumeo));
      l->prefixopci = " ";
      for (i = 0; i < static_cast<int>(trim(adjustl(l->prefixopci1)).length()); ++i) {
         l->prefixopci[i] = l->prefixopci1[i];
         int j = static_cast<int>(static_cast<unsigned char>(l->prefixopci1[i]));
         if (j <= 47) l->prefixopci[i] = '_';
         if (j >= 123) l->prefixopci[i] = '_';
         if ((j >= 58) && (j <= 64)) l->prefixopci[i] = '_';
         if ((j >= 91) && (j <= 96)) l->prefixopci[i] = '_';
         if (j == 46) l->prefixopci[i] = 'p';
      }

      for (i = 0; i < static_cast<int>(trim(adjustl(l->prefixopci)).length()); ++i) {
         while (i + 1 < static_cast<int>(trim(adjustl(l->prefixopci)).length()) && l->prefixopci[i] == '_' && l->prefixopci[i+1] == '_') {
             l->prefixopci.erase(i, 1);
         }
      }
      if (l->prefix.length() > 0 && l->prefix[0] == '_') {
         // !!!acortado 120219  l->nEntradaRoot = trim (adjustl(l->fichin)) // trim (adjustl(prefix))// trim (adjustl(l->prefixopci))
         l->nEntradaRoot = trim(adjustl(l->fichin)) + "_" + trim(adjustl(l->prefixopci));
      } else {
         l->nEntradaRoot = trim(adjustl(l->fichin));
      }
      // !!!l->stochastic
#ifdef CompileWithStochastic
      if (l->stochastic) {
         if (l->layoutnumber <= l->num_procs / 2 - 1) { // aun no se ha dividido el l->num_procs
             l->nEntradaRoot = trim(adjustl(l->nEntradaRoot));
         } else {
             l->nEntradaRoot = trim(adjustl("devia_" + trim(adjustl(l->nEntradaRoot))));
         }
      }
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &l->ierr);
#endif
#endif
      // !!!fin l->stochastic
      // !!!   sgg->nEntradaRoot=trim (adjustl(l->nEntradaRoot))
      //
      chari = std::to_string(l->layoutnumber + 1);
      // Pad chari to 5 chars if necessary, though std::to_string doesn't pad.
      // Assuming Fortran (i5) behavior which pads with spaces.
      if (chari.length() < 5) {
          chari.insert(0, 5 - chari.length(), ' ');
      }
      l->nresumeable2 = trim(adjustl(l->nEntradaRoot)) + "_" + trim(adjustl(chari)) + ".fields";
      //

      l->geomfile = trim(adjustl(l->nEntradaRoot)) + "_" + trim(adjustl(chari));
      // warning file management
      if (statuse != -1) {

CLOSEWARNINGFILE(l.layoutnumber, l.num_procs, l.fatalerror, false, false);
            if ((!l.fatalerror) && (l.layoutnumber == 0)) {
                std::string filename = trim(adjustl(l.fichin)) + "_tmpWarnings.txt_Warnings.txt";
                std::remove(filename.c_str());
            }
            INITWARNINGFILE(l.layoutnumber, l.num_procs, l.nEntradaRoot, l.verbose, l.ignoreerrors);
        }

        if (l.resume_fromold) {
            std::string filename = trim(adjustl(l.nresumeable2)) + ".old";
            resume3 = std::filesystem::exists(filename);
        } else {
            std::string filename = trim(adjustl(l.nresumeable2));
            resume3 = std::filesystem::exists(filename);
        }

        if (l.resume) {
            if (!resume3) {
                stoponerror(l.layoutnumber, l.num_procs, "l%resume fields not present", true);
                statuse = -1;
            }
            dubuf = "RESUMING simulation " + trim(adjustl(l.nEntradaRoot)) + " until n= " + std::to_string(l.finaltimestep);
            print11(l.layoutnumber, dubuf);
        } else {
            if (resume3 && (!l.freshstart) && (!l.run)) {
                stoponerror(l.layoutnumber, l.num_procs, "Restarting file exists. Either specify -r to l%resume, -s to do a fresh START, or -run to run in whatever the case", true);
                statuse = -1;
            } else if (resume3 && l.run) {
                l.resume = true;
            } else {
                {
                    std::ofstream ofs35(trim(adjustl(l.nresumeable2)));
                    ofs35 << "!END" << std::endl;
                }
                {
                    std::ofstream ofs35(trim(adjustl(l.nresumeable2)) + ".old");
                    ofs35 << "!END" << std::endl;
                }
            }
        }

        if (((l.wiresflavor == "slanted") || (l.wiresflavor == "semistructured")) && (l.mpidir != 3)) {
            // arreglado l%mpidir slanted 2019
        }
        if (l.input_conformal_flag && (l.mpidir != 3)) {
            // arreglado l%mpidir conformal 2019
        }
        if (l.run_with_abrezanjas && (l.mpidir != 3)) {
            // arreglado l%mpidir conformal 2019
        }
        if (l.run_with_abrezanjas && l.flag_conf_sgg) {
            // se hace en otro sitio
        }

        if (((l.forcesteps) && (!l.freshstart)) && (statuse != -1)) {
            if (l.resume_fromold) {
                std::string filename = trim(adjustl(l.nresumeable2)) + ".old";
                l.resume = std::filesystem::exists(filename);
            } else {
                std::string filename = trim(adjustl(l.nresumeable2));
                l.resume = std::filesystem::exists(filename);
            }
            if (l.resume) {
                if ((l.layoutnumber == 0) || ((l.layoutnumber == l.num_procs / 2) && l.stochastic)) {
                    if (l.file11isopen) {
                        // CLOSE (11)
                        l.file11isopen = false;
                    }
                    std::string reportFile = trim(adjustl(l.nEntradaRoot)) + "_Report.txt";
                    std::ofstream ofs11(reportFile, std::ios::app);
                    l.file11isopen = true;
                    if (l.resume_fromold) {
                        print11(l.layoutnumber, "Resuming from .fields.old files");
                    } else {
                        print11(l.layoutnumber, "Resuming from .fields files");
                    }
                }
            } else {
                l.freshstart = true;
                l.resume_fromold = false;
            }
            if (((l.layoutnumber == 0) || ((l.layoutnumber == l.num_procs / 2) && l.stochastic)) && l.file11isopen) {
                // close (11)
                l.file11isopen = false;
            }
        }

        if (l.run) {
#ifdef keeppause
            std::string runningFile = "running";
            hayinput = std::filesystem::exists(runningFile);
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, l.ierr);
#endif
            if (hayinput) {
                std::ifstream ifs9(runningFile);
                std::getline(ifs9, chain4);
                chain4 = trim(adjustl(chain4));
                ifs9.close();
#ifdef CompileWithMPI
                MPI_Barrier(SUBCOMM_MPI, l.ierr);
#endif
                removeintraspaces(l.opcionespararesumeo);
                removeintraspaces(chain4);
                if (trim(adjustl(l.opcionespararesumeo)) == trim(adjustl(chain4))) {
                    existiarunningigual = true;
                }
            }
#endif
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, l.ierr);
#endif
            if (l.layoutnumber == 0) {
                std::ofstream ofs38("running");
                ofs38 << trim(adjustl(l.opcionespararesumeo)) << std::endl;
            }
        }
#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, l.ierr);
#endif

        if (((l.layoutnumber == 0) || ((l.layoutnumber == l.num_procs / 2) && l.stochastic)) && (statuse != -1)) {
            std::cout << "Opening _Report.txt file" << std::endl;
            if (l.resume) {
                if (l.file11isopen) {
                    // CLOSE (11)
                    l.file11isopen = false;
                }
                std::string reportFile = trim(adjustl(l.nEntradaRoot)) + "_Report.txt";
                std::ifstream ifs11(reportFile);
                l.file11isopen = true;
                donde = 0;
                while (donde == 0) {
                    std::getline(ifs11, l.chdummy);
                    donde = l.chdummy.find("mpirun -n");
                }
                ifs11.close();
                l.file11isopen = false;
                l.opcionesoriginales = l.chdummy;

                removeintraspaces(l.opcionespararesumeo);
                removeintraspaces(l.opcionesoriginales);
                if (trim(adjustl(l.opcionesoriginales)) != trim(adjustl(l.opcionespararesumeo))) {
                    stoponerror(l.layoutnumber, l.num_procs, "Different resumed/original switches: " + trim(adjustl(l.opcionespararesumeo)) + " <> " + trim(adjustl(l.opcionesoriginales)), true);
                    statuse = -1;
                }
                std::ifstream ifs11_2(reportFile);
                l.file11isopen = true;
                donde = 0;
                while (donde == 0) {
                    std::getline(ifs11_2, l.chdummy);
                    donde = l.chdummy.find("!SLICES");
                }
            }
        }

}
            l.slicesoriginales = trim(adjustl(l.chdummy));
            close(11);
            l.file11isopen = false;
            open(11, (trim(adjustl(l.nEntradaRoot)) + "_Report.txt").c_str(), ios::out | ios::app);
            l.file11isopen = true;
            if (l.layoutnumber == 0) insertalogtmp(l);
         } else {
            close(11);
            l.file11isopen = false;
            open(11, (trim(adjustl(l.nEntradaRoot)) + "_Report.txt").c_str(), ios::out);
            l.file11isopen = true;
            if (l.layoutnumber == 0) insertalogtmp(l);
         }
         
         get_secnds(l.time_out2);
         print_credits(l);
#ifdef CompileWithReal8
         dubuf = "Compiled with Double precision (real*8)";
         print11(l.layoutnumber, dubuf);
#endif
#ifdef CompileWithReal4
         dubuf = "Compiled with Single precision (real*4)";
         print11(l.layoutnumber, dubuf);
#endif
#ifdef CompileWithReal16
         dubuf = "Compiled with Quadruple precision (real*16)";
         print11(l.layoutnumber, dubuf);
#endif
         dubuf = SEPARADOR + SEPARADOR + SEPARADOR;
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = "Launched on              " + l.time_out2.fecha.substr(6, 2) + "/" + l.time_out2.fecha.substr(4, 2) + "/" + l.time_out2.fecha.substr(0, 4) + " " + l.time_out2.hora.substr(0, 2) + ":" + l.time_out2.hora.substr(2, 2);
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = SEPARADOR + SEPARADOR + SEPARADOR;
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = "Launched with total options ";
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = trim(adjustl(l.opcionestotales));
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = "If later resuming use compulsory options ";
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = trim(adjustl(l.opcionespararesumeo));
         // !!!call print11 (l%layoutnumber, dubuf,.true.)
         write(11, trim(adjustl(dubuf))); // a capon para que el l%stochastic pueda resumear
         dubuf = SEPARADOR + SEPARADOR + SEPARADOR;
         print11(l.layoutnumber, dubuf);
      }
      
      //
      //
      //in seconds
      l.flushsecondsFields = l.flushminutesFields * 60;
      //in seconds
      l.flushsecondsData = l.flushminutesData * 60;

      if ((!l.existeNFDE) && (!l.existeh5)) {
         stoponerror(l.layoutnumber, l.num_procs, "Some input file missing .h5/.nfde/.conf", true);
         statuse = -1;
         // goto 668
      }
      //
      //
#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, l.ierr);
#endif
      if (existiarunningigual) { // lo pongo aqui pq si no no se escribe en el report
         stoponerror(l.layoutnumber, l.num_procs, "Running flag file with same options than requested exist. ", true);
         statuse = -1;
      }

668:
      ;

      input_conformal_flag = l.input_conformal_flag; // es un flag global!!!!ojooo 051223 !devolverlo correctamente
      return; // el unico return que he dejado !240817

   } // end subroutine interpreta

   void insertalogtmp(entrada_t& l) { // para 100920
      char dubuf[BUFSIZE];
      int MYUNIT11;
      OffPrint(); // no reimprimas, esto ya estaba por pantalla
      open_new_unit(MYUNIT11, "SEMBA_FDTD_temp.log");
      while (true) {
         if (!read_line(MYUNIT11, dubuf, 1024)) break;
         dubuf[0] = '&';
         memmove(dubuf + 1, dubuf, strlen(dubuf));
         print11(l.layoutnumber, dubuf);
      }
      close_unit_delete(MYUNIT11);
      OnPrint();
      return;
   }

   void print_basic_help(entrada_t& l) {
      print_credits(l);
      print11(l.layoutnumber, "___________________________________________________________________________");
      print11(l.layoutnumber, "Basic usage: ");
      print11(l.layoutnumber, "&   For help use          -h ");
      print11(l.layoutnumber, "&   For launching use                     ");
      print11(l.layoutnumber, "&                         -i inputfile (native)");
      print11(l.layoutnumber, "___________________________________________________________________________");
      return;
   }

   void print_credits(entrada_t& l) {
      char dubuf[BUFSIZE];

      if (l.creditosyaprinteados) return;
      l.creditosyaprinteados = true;
      print11(l.layoutnumber, "=========================");
      print11(l.layoutnumber, program_name);
      print11(l.layoutnumber, "=========================");

      dubuf[0] = '\0';
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      print11(l.layoutnumber, dubuf);
      print11(l.layoutnumber, ("Compilation date: " + std::string(compilation_date)).c_str());
      print11(l.layoutnumber, ("Compiler Id: " + std::string(compiler_id)).c_str());
      print11(l.layoutnumber, ("git commit: " + std::string(git_commit)).c_str());
      print11(l.layoutnumber, ("cmake build type: " + std::string(cmake_build_type)).c_str());
      if (std::string(cmake_build_type) == "Debug") {
         print11(l.layoutnumber, ("cmake compilation flags: " + std::string(compilation_flags_debug)).c_str());
      } else if (std::string(cmake_build_type) == "Release") {
         print11(l.layoutnumber, ("cmake compilation flags: " + std::string(compilation_flags_release)).c_str());
      } else {
         print11(l.layoutnumber, ("cmake compilation flags: " + std::string(compilation_flags)).c_str());
      }
      dubuf[0] = '\0';
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      print11(l.layoutnumber, dubuf);
      dubuf[0] = '\0';
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      print11(l.layoutnumber, dubuf);
      print11(l.layoutnumber, "All rights reserved by the University of Granada (Spain)");
      print11(l.layoutnumber, "       Contact person: Luis D. Angulo <lmdiazangulo@ugr.es>");
      print11(l.layoutnumber, " ");
      // *******************************************************************************

      dubuf[0] = '\0';
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      print11(l.layoutnumber, dubuf);
#ifdef CompileWithMPI
      print11(l.layoutnumber, "Compiled WITH MPI support");
#endif
#ifdef CompileWithHDF
      print11(l.layoutnumber, "Compiled WITH .h5 HDF support");
#endif
#ifdef CompileWithMTLN
      print11(l.layoutnumber, "Compiled WITH MTLN support");
#endif
#ifdef CompileWithSMBJSON
      print11(l.layoutnumber, "Compiled WITH SMBJSON support");
#endif
      dubuf[0] = '\0';
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      strcat(dubuf, SEPARADOR);
      print11(l.layoutnumber, dubuf);
      get_secnds(l.time_out2);
      dubuf[0] = '\0';
      sprintf(dubuf, "Launched on              %s/%s/%s %s:%s", 
              l.time_out2.fecha.substr(6, 2).c_str(),
              l.time_out2.fecha.substr(4, 2).c_str(),
              l.time_out2.fecha.substr(0, 4).c_str(),
              l.time_out2.hora.substr(0, 2).c_str(),
              l.time_out2.hora.substr(2, 2).c_str());
      print11(l.layoutnumber, dubuf);
      if (l.layoutnumber == 0) {
         std::cout << "Highest integer " << std::numeric_limits<int>::max() << std::endl;
      }
      return;
   }

   void print_help(entrada_t& l) {
      char buff[BUFSIZE];
      print11(l.layoutnumber, "___________________________________________________________________________");
      print11(l.layoutnumber, "Command line arguments: ");
      print11(l.layoutnumber, "___________________________________________________________________________");
      print11(l.layoutnumber, "-i geometryfile        : Simulates the Native format input file            ");
   }

print11(l.layoutnumber, "-r                     : Restarts a previous execution until a given step. ");
        print11(l.layoutnumber, "&                        Needs -n                                          ");
        print11(l.layoutnumber, "-run                   : Uses a semaphore running file and automatically   ");
        print11(l.layoutnumber, "&                        relaunches simulation if ended or aborted (cluter)");
#ifdef CompileWithOldSaving
        print11(l.layoutnumber, "-old                   : Jointly with -r restarts from .fields.old files   ");
        print11(l.layoutnumber, "&                        instead (for safety .fields.old fields are saved  ");
        print11(l.layoutnumber, "&                        too if -flush is issued)                          ");
#endif
        print11(l.layoutnumber, "-cfl number            : Courant number (suggested<=0.8)  overriding input ");
        print11(l.layoutnumber, "-n numberoftimesteps   : Run the simulation until a specified step         ");
        print11(l.layoutnumber, "&                        either restarting if the necessary files are      ");
        print11(l.layoutnumber, "&                        present, or starting a fresh new one otherwise    ");
        print11(l.layoutnumber, "&                        Special cases: n=-1 -> Run only .h5/.nfde preproc.");
        print11(l.layoutnumber, "&                        Special cases: n=-2 -> Run only .h5 preprocessing ");
        print11(l.layoutnumber, "-s                     : Forces a fresh new simulation, erasing the        ");
        print11(l.layoutnumber, "&                        restarting files if they are present              ");
        print11(l.layoutnumber, "&                        Jointly with -n, it enforces a fresh restart      ");
        print11(l.layoutnumber, "&                        (erases .fields files from previous simulations)  ");
        print11(l.layoutnumber, "___________________________________________________________________________");
        print11(l.layoutnumber, "-pause seconds         : Wait seconds to start simulation                  ");
        print11(l.layoutnumber, "-prefix string         : Adds a string to the output filenames             ");
        print11(l.layoutnumber, "-saveall               : Saves all the observation time steps              ");
        print11(l.layoutnumber, "&                        (default saves only the specified windows of time)");
        print11(l.layoutnumber, "-singlefile            : Compacts E, H, J probes in single files to        ");
        print11(l.layoutnumber, "&                        overcome a large number of file openings          ");
        // !!#ifdef CompileWithMPI
        // !!          print11 (l.layoutnumber, "-maxmessages number    : Buffer of messages for MPI Warnings file. Just    ");
        // !!          print11 (l.layoutnumber, "&                        increase if requested at runtime                  ");
        // !!#endif

#ifdef CompileWithNIBC
        print11(l.layoutnumber, "-skindepthpre          : Pre-processor for sgbc metals including skin depth.");
        print11(l.layoutnumber, "-mibc                  : Uses pure l.mibc to deal with composites.  ");
        print11(l.layoutnumber, "-ade                   : Uses l.ade-l.mibc to deal with composites. ");
        print11(l.layoutnumber, "&                        Alternative to -l.mibc.");
        print11(l.layoutnumber, "-conformalskin         : Uses a conformal l.mibc to deal with skin-depth");
        print11(l.layoutnumber, "&                        Do not use this switch if the problem also involves ");
        print11(l.layoutnumber, "&                        traditional composites, since these do not hold the right ");
        print11(l.layoutnumber, "&                        thickness parameter. Only use it if the problem only ");
        print11(l.layoutnumber, "&                        contains metals for which both the conductivity and ");
        print11(l.layoutnumber, "&                        thickness are CORRECTLY specified in the .nfde file. ");
        print11(l.layoutnumber, "-NOcompomur            : Uses OLD (possibly unstable) upwinding scheme to deal with composites, ");
        print11(l.layoutnumber, "&                        instead of the NEW default, which uses a causal time-domain extrapolation ");
        print11(l.layoutnumber, "&                        of magnetic fields at the surface, by using the one-way ");
        print11(l.layoutnumber, "&                        advection equation (similar to 1D Mur ABCs) for its ");
        print11(l.layoutnumber, "&                        superior stability of the default new Mur formulation");
        print11(l.layoutnumber, "-attc   dissipation    : Positive factor (under 1) for stable composites,   ");
        print11(l.layoutnumber, "&                        permits to solve some instabilities in the simulation of l.mibc materials.");
        print11(l.layoutnumber, "&                        It just adds a 1 cell lossy magnetic coating to the l.mibc composite.");
        print11(l.layoutnumber, "&                        The dissipation factor is used to find the magnetic conductivity ");
        print11(l.layoutnumber, "&                        from the coefficient updating the current magnetic ");
        print11(l.layoutnumber, "&                        field from the previous one.  ");
        std::ostringstream buff_stream;
        buff_stream << std::setprecision(2) << std::scientific << "&                        Default= " << l.attfactorc;
        std::string buff = buff_stream.str();
        print11(l.layoutnumber, buff);
#endif
        print11(l.layoutnumber, "-prioritizeCOMPOoverPEC: Uses Composites instead of PEC in conflicts.       ");
        print11(l.layoutnumber, "-prioritizeISOTROPICBODYoverall: Uses ISOTROPIC BODY FOR conflicts (JUST FOR SIVA).       ");
        print11(l.layoutnumber, "-sgbc               : Enables the defaults sgbc model for composites. Default sgbc:");
        print11(l.layoutnumber, "-nosgbc             : Disables the defaults sgbc model for composites. Default sgbc:");
        print11(l.layoutnumber, "&                        -sgbfreq 3e9 -sgbresol 1 -sgbcrank      ");
        print11(l.layoutnumber, "-sgbcfreq           : Maximum frequency to consider the skin-depth       ");
        print11(l.layoutnumber, "-sgbcresol          : Number of cells per skin-depth a the Maximum frequency");
        print11(l.layoutnumber, "-sgbcyee            : Uses pure Yee ETD sgbc instead of Crank-Nicolson");
        print11(l.layoutnumber, "-sgbccrank          : Uses sgbc Crank-Nicolson (default)        ");
        print11(l.layoutnumber, "-sgbcdepth number   : Overrides automatic calculation of number of cells ");
        print11(l.layoutnumber, "&                        within sgbc                              ");

        print11(l.layoutnumber, "-pmlalpha factor order : CPML Alpha factor (>=0, <1 sug.) & polyn. grading.");
        print11(l.layoutnumber, "&                        alpha=factor * maximum_PML_sigma , order=polynom. ");
        std::ostringstream buff_stream2;
        buff_stream2 << std::setprecision(2) << std::scientific << "&                        Default= " << l.alphamaxpar << " " << l.alphaOrden;
        std::string buff2 = buff_stream2.str();
        print11(l.layoutnumber, buff2);
        std::ostringstream buff_stream3;
        buff_stream3 << std::setprecision(2) << std::scientific << "-pmlkappa number       : CPML Kappa (>=1). Default= " << l.kappamaxpar;
        std::string buff3 = buff_stream3.str();
        print11(l.layoutnumber, buff3);
        print11(l.layoutnumber, "-pmlcorr factor depth  : Factor for CPML enhanced stability (default none).");
        print11(l.layoutnumber, "&                        sigma=factor * maximum_PML_sigma, depth= # layers ");
        print11(l.layoutnumber, "-mur1                  : Supplement PMLs with 1st order Mur ABCs           ");
        print11(l.layoutnumber, "-mur2                  : Supplement PMLs with 2nd order Mur ABCs           ");
        print11(l.layoutnumber, "-wiresflavor {holland.or.old} : model for the wires    ");
#ifdef CompileWithBerengerWires

print11(l.layoutnumber, "-wiresflavor {berenger} : model for the wires    ");
#endif
#ifdef CompileWithSlantedWires
        print11(l.layoutnumber, "-wiresflavor {new/Slanted.or.experimental.or.slanted/transition/semistructured l.precision} : model for the wires    ");
#endif
        print11(l.layoutnumber, "&                        (default " + trim(adjustl(l.wiresflavor)) + ")   ");
        print11(l.layoutnumber, "-notaparrabos          : Do not remove extra double tails at the end of the wires ");
        print11(l.layoutnumber, "&                        only available for the native format.             ");
        print11(l.layoutnumber, "-intrawiresimplify     : Disable strict interpretation of .NFDE topology.  ");
        print11(l.layoutnumber, "&                        Collapse internal parallel wires and create       ");
        print11(l.layoutnumber, "&                        intra-wire junctions.                             ");
        print11(l.layoutnumber, "-nomtlnberenger        : Disables MTLN improvements for Berenger l.wiresflavor");
        print11(l.layoutnumber, "-stableradholland             : Automatic correction of radii for Holland l.wiresflavor");
        print11(l.layoutnumber, "&                        Use only in case of instabilities.  (experimental)");
        print11(l.layoutnumber, "-groundwires           : Ground wires touching/embedded/crossing PEC/Lossy.");
        print11(l.layoutnumber, "&                        Use with CAUTION. Revise *Warnings.txt file!      ");
        print11(l.layoutnumber, "-noSlantedcrecepelo : Ground open nodes. Experimental. Do not use.");
        print11(l.layoutnumber, "-connectendings        : Joins ohmicly endings nodes of adjacent segments  ");
        print11(l.layoutnumber, "&                        from multiwires (segments do no collapse).        ");
        print11(l.layoutnumber, "&                        regardless of whether they are actually connected ");
        print11(l.layoutnumber, "&                        through the LeftEnd/RightEnd numbering ");
        print11(l.layoutnumber, "&                        Automatic with -a                                 ");
        print11(l.layoutnumber, "&                        Use with CAUTION. Revise *Warnings.txt file!      ");
        print11(l.layoutnumber, "-isolategroupgroups    : Detach ohmicly endings nodes of adjacent segments ");
        print11(l.layoutnumber, "&                        from multiwires if they are in different          ");
        print11(l.layoutnumber, "-makeholes             : Create a void 2-cell area around wire segments    ");
        print11(l.layoutnumber, "&                        Use with CAUTION. Revise *Warnings.txt (experim.) ");
        print11(l.layoutnumber, "-mindistwires dist     : Specify the min distance between wires in a       ");
        print11(l.layoutnumber, "&                        multiwire in new and experimental wires flavors   ");
        write(buff, "(a,e10.2e3)")("&                        Default= ", l.mindistwires);
        print11(l.layoutnumber, buff);
        print11(l.layoutnumber, "-inductance {ledfelt/berenger/boutayeb} : model for the self-inductance    ");
        print11(l.layoutnumber, "&                        (default " + trim(adjustl(l.inductance_model)) + ")   ");
        print11(l.layoutnumber, "-inductanceorder order : order for the self-inductance calculation for     ");
        print11(l.layoutnumber, "&                        slanted wires in experimental l.wiresflavor         ");
        write(buff, "(a,i8)")("&                        Default= ", l.inductance_order);
        print11(l.layoutnumber, buff);
        print11(l.layoutnumber, "-attw   dissipation    : Positive factor (under 1) for stability in wires, ");
        write(buff, "(a,e10.2e3)")("&                        Default= ", l.attfactorw);
        print11(l.layoutnumber, buff);
        print11(l.layoutnumber, "-maxwireradius number  : Bounds globally the wire radius                   ");
        print11(l.layoutnumber, "-clip                  : Permits to clip a bigger problem truncating wires.");
        print11(l.layoutnumber, "-wirecrank             : Uses Crank-Nicolson for wires (development)       ");
        print11(l.layoutnumber, "-noNF2FF string        : Supress a NF2FF plane for calculation             ");
        print11(l.layoutnumber, "&                        String can be: up, down, left, right, back , front");
        print11(l.layoutnumber, "-NF2FFDecim            : Uses decimation in NF2FF calculation (faster).    ");
        print11(l.layoutnumber, "&                        WARNING: High-freq aliasing may occur             ");
        print11(l.layoutnumber, "-vtkindex              : Output index instead of real point in 3D slices.  ");
        print11(l.layoutnumber, "-ignoreerrors          : Run even if errors reported in *Warnings.txt file.");
        print11(l.layoutnumber, "___________________________________________________________________________");
        print11(l.layoutnumber, "-cpumax minutes        : CPU runtime (useful for limited CPU queuing       ");
        print11(l.layoutnumber, "-noshared              : Do not waste time with shared fields              ");
        print11(l.layoutnumber, "-flush minutes         : Minutes between data flush of restarting fields   ");
        print11(l.layoutnumber, "&                        (default 0=No flush)                              ");
        print11(l.layoutnumber, "-flushdata minutes     : Minutes between flushing observation data         ");
        print11(l.layoutnumber, "&                        (default is every 5 minutes)                      ");
        print11(l.layoutnumber, "-map                   : Creates map ASCII files of the geometry           ");
        print11(l.layoutnumber, "&                        with wires and PEC                ");
        print11(l.layoutnumber, "&                        (in conjunction with -n 0 only creates the maps)  ");
        print11(l.layoutnumber, "-mapvtk                : Creates .VTK map of the PEC/wires/Surface geometry");
        print11(l.layoutnumber, "-dmma                  : Thin-gaps treated in DMMA manner  ");
#ifdef CompileWithMPI
        print11(l.layoutnumber, "-mpidir {x,y,z}        : Rotate model to force MPI along z be the largest  ");
        print11(l.layoutnumber, "-force    cutplane     : Force a MPI layout to begin at cutplane (debug!)  ");
#endif
        print11(l.layoutnumber, "___________________________________________________________________________");
        print11(l.layoutnumber, "Control through signaling files during the simulation: (after erased)      ");
        print11(l.layoutnumber, "&  stop         : (void) Forces a graceful end (it Cannot be resumed)      ");
        print11(l.layoutnumber, "&                 No restarting data is flushed, only observation data     ");
        print11(l.layoutnumber, "&  stopflushing : (void) Forces a graceful end (it can be resumed)         ");
        print11(l.layoutnumber, "&  flush        : (void) Forces a flush of resuming fields and observation ");
        print11(l.layoutnumber, "&                 data in 1 minute time approx.                            ");
        print11(l.layoutnumber, "&  flushdata    : (void) Forces a flush only of the observation data in    ");
        print11(l.layoutnumber, "&                 1 minute time approx.                                    ");
        print11(l.layoutnumber, "&                 Both restarting and observation data are flushed         ");
        print11(l.layoutnumber, "&  stop_only         : Forces a graceful end (cannot be resumed) only of a ");
        print11(l.layoutnumber, "&                      given problem name (without the .nfde extension)    ");
        print11(l.layoutnumber, "&                      No restarting data is flushed, only observation data");

print11(l.layoutnumber, "&  stopflushing_only : Forces a graceful end (it can be resumed) only of a ");
      print11(l.layoutnumber, "&                      give problem name (without the .nfde extension)     ");
      print11(l.layoutnumber, "&                      Both restarting and observation data is flushed     ");
      print11(l.layoutnumber, "&  flush_only   : Forces flush of resuming fields and observation data only");
      print11(l.layoutnumber, "&                 of a given problem name (without the .nfde extension)    ");
      print11(l.layoutnumber, "&                 in 1 minute time approx.                                 ");
      print11(l.layoutnumber, "&  flushdata_only : Forces a flush only of the observation data only of a  ");
      print11(l.layoutnumber, "&                   given problem name (without the .nfde extension)       ");
      print11(l.layoutnumber, "&                   in 1 minute time approx.                               ");
      print11(l.layoutnumber, "&                   Both restarting and observation data are flushed       ");
      print11(l.layoutnumber, "&  pause        : (void) While this field exist no simulation is started   ");
      print11(l.layoutnumber, "&  unpack       : (void) Unpacks on-the-fly .bin probes files created      ");
      print11(l.layoutnumber, "&                 with the -singlefile packaging option                    ");
      print11(l.layoutnumber, "&  postprocess  : (void) Do frequency domain and transfer function         ");
      print11(l.layoutnumber, "&                 postprocess on-the-fly                                   ");
      print11(l.layoutnumber, "&  flushxdmf    : (void) Flush .xdmf animation probes on the fly           ");
      print11(l.layoutnumber, "&  flushvtk     : (void) Flush .vtk  animation probes on the fly           ");
      print11(l.layoutnumber, "&  snap         : Creates a .h5 and .xdmf snapshot per MPI layout if the   ");
      print11(l.layoutnumber, "&                 field value is over the first number found in this file  ");
      print11(l.layoutnumber, "&                 in space steps by the 2nd integer number                 ");
      print11(l.layoutnumber, "&                 in time steps by the 3rd integer number (1-minute lapse) ");
      print11(l.layoutnumber, "&  relaunch     : Relaunches the simulation upon termination with the      ");
      print11(l.layoutnumber, "&                 switches read from this file. Used jointly with a        ");
      print11(l.layoutnumber, "&                 stop file permits to launch simulations on-demand        ");
      print11(l.layoutnumber, "___________________________________________________________________________");
      //
      sprintf(buff, "Max CPU time is %14d seconds (can be overriden by -cpumax)", topCPUtime);
      print11(l.layoutnumber, buff);
#ifdef CompileWithOpenMP
      print11(l.layoutnumber, "SUPPORTED:   MultiCPU parallel simulation (OpenMP)");
#endif
#ifdef CompileWithMPI
      print11(l.layoutnumber, "SUPPORTED:   MultiCPU/Multinode parallel simulation (MPI)");
#endif
      print11(l.layoutnumber, "SUPPORTED:   Near-to-Far field probes");
      print11(l.layoutnumber, "SUPPORTED:   Lossy anistropic materials, both electric and magnetic");
      print11(l.layoutnumber, "SUPPORTED:   Thin Slots ");
      print11(l.layoutnumber, "SUPPORTED:   Electric and Magnetic Dispersive materials ");
      print11(l.layoutnumber, "SUPPORTED:   Isotropic Multilayer Skin-depth Materials (sgbc)");
#ifdef CompileWithNIBC
      print11(l.layoutnumber, "SUPPORTED:   Isotropic Multilayer Skin-depth Materials (l.mibc)");
#endif
      print11(l.layoutnumber, "SUPPORTED:   Loaded and grounded thin-wires with juntions");
      print11(l.layoutnumber, "SUPPORTED:   Nodal hard/soft electric and magnetic sources");
#ifdef CompileWithHDF
      print11(l.layoutnumber, "SUPPORTED:   .xdmf+.h5 probes ");
#endif
#ifdef CompileWithOldSaving
      print11(l.layoutnumber, "SUPPORTED:   .fields.old files created (fail-safe)");
#endif
#ifdef CompileWithStochastic
      print11(l.layoutnumber, "SUPPORTED:   l.stochastic analysis");
#endif
#ifdef CompileWithPrescale
      print11(l.layoutnumber, "SUPPORTED:   Permittivity scaling accelerations");
#endif
      print11(l.layoutnumber, "SUPPORTED:   Holland Wires");
#ifdef CompileWithBerengerWires
      print11(l.layoutnumber, "SUPPORTED:   Multi-Wires");
#endif
#ifdef CompileWithSlantedWires
      print11(l.layoutnumber, "SUPPORTED:   Slanted Wires");
#endif
      //
#ifdef CompileWithReal4
      print11(l.layoutnumber, "Single precission simulations (reals are 4-byte)");
#endif
#ifdef CompileWithReal8
      print11(l.layoutnumber, "Double precission simulations (reals are 8-byte)");
#endif
#ifdef CompileWithInt4
      print11(l.layoutnumber, "Media matrices are 4 bytes");
#endif
#ifdef CompileWithInt2
      print11(l.layoutnumber, "Media matrices are 2 bytes");
#endif
#ifdef CompileWithInt1
      print11(l.layoutnumber, "Media matrices are 1 byte");
#endif
#ifdef CompileWithMPI
      MPI_FINALIZE(l.ierr);
#endif
      return;
   }

   void removeintraspaces(std::string& a) {
      int i, longi;
      bool correc;
      correc = true;
      while (correc) {
         correc = false;
         // trim left
         size_t start = a.find_first_not_of(' ');
         if (start == std::string::npos) {
            a = "";
            return;
         }
         a = a.substr(start);
         // trim right
         size_t end = a.find_last_not_of(' ');
         if (end != std::string::npos) {
            a = a.substr(0, end + 1);
         } else {
            a = "";
         }
         
         longi = static_cast<int>(a.length());
         for (i = 0; i < longi - 1; ++i) {
            if (a[i] == ' ' && a[i + 1] == ' ') {
               a = a.substr(0, i) + a.substr(i + 1);
               correc = true;
               break;
            }
         }
      }
      return;
   }

   //
   //
   //

   void buscaswitchficheroinput(entrada_t& l) {
      //
      char dato[BUFSIZE], buff[BUFSIZE], f[BUFSIZE], binaryPath[BUFSIZE];
      int i, n, statuse, NUM_NFDES, TEMP_NUMNFDES, p;
      char NFDEEXTENSION[6], CONFEXTENSION[6], CMSHEXTENSION[6];

      //
      strcpy(NFDEEXTENSION, ".nfde");
      strcpy(CONFEXTENSION, ".conf");
      strcpy(CMSHEXTENSION, ".cmsh");
      statuse = 0;
      //
      getBinaryPath(binaryPath);
      n = commandargumentcount(l.chain2, binaryPath);
      if (n == 0) {
         print_basic_help(l);
         stoponerror(l.layoutnumber, l.num_procs, "Error: NO arguments neither command line nor in launch file. Correct and remove pause...", true);
         statuse = -1;
         goto label_667;
      }

      if (n > 0) {
         num_nfdes = 0;
         i = 2;
         while (i <= n) {
            getcommandargument(l.chain2, i, l.chain, l.length, statuse, binaryPath);
            if (statuse != 0) {
               stoponerror(l.layoutnumber, l.num_procs, "Reading input", true);
               goto label_667;
            }
            //
            std::string chain_trimmed = l.chain;
            // trim and adjustl equivalent
            size_t start = chain_trimmed.find_first_not_of(' ');
            if (start != std::string::npos) {
               chain_trimmed = chain_trimmed.substr(start);
               size_t end = chain_trimmed.find_last_not_of(' ');
               if (end != std::string::npos) {
                  chain_trimmed = chain_trimmed.substr(0, end + 1);
               }
            } else {
               chain_trimmed = "";
            }

            if (chain_trimmed == "-mpidir") {
               i = i + 1;
               getcommandargument(l.chain2, i, f, l.length, statuse, binaryPath);
               std::string f_trimmed = f;
               size_t start_f = f_trimmed.find_first_not_of(' ');
               if (start_f != std::string::npos) {
                  f_trimmed = f_trimmed.substr(start_f);
                  size_t end_f = f_trimmed.find_last_not_of(' ');
                  if (end_f != std::string::npos) {
                     f_trimmed = f_trimmed.substr(0, end_f + 1);
                  }
               } else {
                  f_trimmed = "";
               }
               
               if (f_trimmed == "x" || f_trimmed == "X") {
                  l.mpidir = 1;
               } else if (f_trimmed == "y" || f_trimmed == "Y") {
                  l.mpidir = 2;
               } else if (f_trimmed == "z" || f_trimmed == "Z") {
                  l.mpidir = 3;
               } else {
                  goto label_1762;
               }
               goto label_2762;
label_1762:
               stoponerror(l.layoutnumber, l.num_procs, "Invalid -l.mpidir option", true);
               statuse = -1;
               goto label_667;
label_2762:
               continue;
            } else if (chain_trimmed == "-h") {
               print_credits(l);
               print_help(l);
               print_credits(l);
               exit(0); // STOP equivalent in C++
            } else if (chain_trimmed == "-i") {
               num_nfdes = num_nfdes + 1;
            }
            i = i + 1;
         }
      }
      
label_667:;
   }

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <filesystem>

// Assuming necessary types and functions are available in the global scope or a namespace
// Based on the Fortran code, we assume:
// - l is a reference to a struct/class 'entrada_t'
// - getcommandargument, stoponerror, INITWARNINGFILE are functions
// - NFDEEXTENSION is a global string constant
// - RKIND is a type alias for double (or similar)
// - SUBCOMM_MPI, l.ierr are MPI related (if CompileWithMPI is defined)

// Helper to trim whitespace
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

// Placeholder for getcommandargument
// Fortran: call getcommandargument(l%chain2, i, l%chain, l%length, statuse, binaryPath)
// Assuming l.chain2 is a vector of strings (arguments), l.chain is output string, l.length is output length
void getcommandargument(const std::vector<std::string>& args, int i, std::string& out, int& out_len, int& statuse, const std::string& binaryPath);

// Placeholder for stoponerror
void stoponerror(int layoutnumber, int num_procs, const std::string& message, bool fatal);

// Placeholder for INITWARNINGFILE
void INITWARNINGFILE(int layoutnumber, int num_procs, const std::string& filename, bool verbose, bool ignoreerrors);

// Global constant
extern const std::string NFDEEXTENSION;

// Forward declaration of entrada_t
struct entrada_t;

void buscaswitchficheroinput_impl(entrada_t& l, int n, int num_nfdes) {
    int statuse = 0;
    std::string binaryPath = ""; // Assuming binaryPath is passed or global
    std::string f = "";
    std::string dato = "";
    std::string buff = "";
    int p = 0;
    int i = 0;
    int temp_numnfdes = 0;

    if (num_nfdes > 1) {
        temp_numnfdes = 0;
        i = 2; // se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
        while (i <= n) {
            getcommandargument(l.chain2, i, l.chain, l.length, statuse, binaryPath);
            if (statuse != 0) {
                stoponerror(l.layoutnumber, l.num_procs, "Reading input", true);
                goto label_667;
            }
            
            std::string trimmed_chain = trim(l.chain);
            
            if (trimmed_chain == "-i") {
                temp_numnfdes = temp_numnfdes + 1;
                i = i + 1;
                getcommandargument(l.chain2, i, f, l.length, statuse, binaryPath);
                p = trim(f).length();
                
                if ((p - 4) >= 1) {
                    if (f[(p - 4)] == NFDEEXTENSION[0]) {
                        // Extract extension
                        std::string ext = f.substr(p - 4, 4);
                        NFDEEXTENSION = ext; // Assuming NFDEEXTENSION can be modified or is a reference
                        l.extension = NFDEEXTENSION;
                        l.fichin = f.substr(0, p - 5);
                    } else {
                        l.fichin = f.substr(0, p);
                    }
                } else if (p >= 1) {
                    l.fichin = f.substr(0, p);
                } else {
                    stoponerror(l.layoutnumber, l.num_procs, "There is not a .nfde file for input", true);
                    statuse = -1;
                    goto label_667;
                }
                
                // Inquire file existence
                std::string full_filename = trim(l.fichin) + NFDEEXTENSION;
                l.existeNFDE = std::filesystem::exists(full_filename);
                
                if (!l.existeNFDE) {
                    buff = "The input file was not found " + trim(l.fichin) + NFDEEXTENSION;
                    stoponerror(l.layoutnumber, l.num_procs, buff, true);
                    statuse = -1;
                    goto label_667;
                }

                if (temp_numnfdes == 1) { // solo el primero
                    if (l.layoutnumber == 0) {
                        // Open file for writing merged content
                        std::string multi_filename = "multi_" + trim(l.fichin) + NFDEEXTENSION;
                        std::ofstream file194(multi_filename, std::ios::out);
                        if (!file194.is_open()) {
                             // Handle error appropriately
                        }
                        // Store file handle or path for later use. 
                        // Since we can't easily map Fortran unit numbers in C++, we'll manage files by name/path
                        // This part is tricky to translate directly without more context on how file units are managed.
                        // We will assume a helper class or global map manages these files.
                        // For simplicity, let's assume we open them and keep them open or re-open.
                        // Given the complexity, we might need to refactor this part significantly.
                        // However, sticking to the prompt, we translate logic.
                        
                        // Note: Direct translation of Fortran unit numbers to C++ file streams is complex.
                        // We will assume a simplified approach where we manage files by name.
                        // Let's assume l has members to store file streams or paths for units 194 and 196.
                        // Since we don't have the struct definition, we'll use local variables and assume they are managed.
                        
                        // To make this compile, we need to assume some structure.
                        // Let's assume l has:
                        // std::ofstream* file194;
                        // std::ofstream* file196;
                        // But since we can't modify the struct, we'll use a different approach.
                        // We will open files and write to them. If we need to keep them open, we need a way to track them.
                        // Given the constraints, I will assume a global or class-level file manager exists.
                        // For the sake of this exercise, I will open and close files as needed, or assume they are managed elsewhere.
                        // Actually, looking at the code, it reads from 196 and writes to 194.
                        // It opens 196, reads line by line, writes to 194.
                        // Then closes 196.
                        // Then if temp_numnfdes == num_nfdes, it writes '!END' to 194 and closes 194.
                        
                        // We need to keep 194 open across iterations.
                        // Let's assume l has a member std::ofstream* merged_file;
                        // And we open it here if it's the first time.
                        
                        // Since I cannot change the struct, I will assume the existence of a helper function or global state.
                        // Let's assume there is a function open_unit(int unit, const std::string& filename)
                        // and close_unit(int unit).
                        
                        // For now, I will write the logic assuming we can open files.
                        // We will use a map to store file streams.
                        static std::map<int, std::ofstream> file_streams;
                        
                        if (file_streams.find(194) == file_streams.end()) {
                            file_streams[194].open(multi_filename, std::ios::out);
                        }
                    }
                }
                
                if (l.layoutnumber == 0) {
                    std::string input_filename = trim(l.fichin) + NFDEEXTENSION;
                    std::ifstream file196(input_filename, std::ios::in);
                    if (!file196.is_open()) {
                        stoponerror(l.layoutnumber, l.num_procs, "Cannot open input file", true);
                        statuse = -1;
                        goto label_667;
                    }
                    
                    std::string line;
                    while (std::getline(file196, line)) {
                        dato = trim(line);
                        if (dato != "!END") {
                            if (l.layoutnumber == 0) {
                                // Write to 194
                                if (file_streams.find(194) != file_streams.end()) {
                                    file_streams[194] << dato << "\n";
                                }
                            }
                        } else {
                            dato = "***** End merging file: " + trim(l.fichin) + NFDEEXTENSION + " ********";
                            if (l.layoutnumber == 0) {
                                if (file_streams.find(194) != file_streams.end()) {
                                    file_streams[194] << dato << "\n";
                                }
                            }
                        }
                    }
                    file196.close();
                }
                
                if (temp_numnfdes == num_nfdes) { // solo el primero
                    if (l.layoutnumber == 0) {
                        if (file_streams.find(194) != file_streams.end()) {
                            file_streams[194] << "!END" << "\n";
                            file_streams[194].close();
                            file_streams.erase(194);
                        }
                    }
                }
            }
            i = i + 1;
        }
    }
    
#ifdef CompileWithMPI
    // MPI_Barrier call
    // Assuming MPI is available
    MPI_Barrier(SUBCOMM_MPI, &l.ierr);
#endif

    temp_numnfdes = 0;
    if (n > 0) {
        i = 2; // se empieza en 2 porque el primer argumento es siempre el nombre del ejecutable
        while (i <= n) {
            getcommandargument(l.chain2, i, l.chain, l.length, statuse, binaryPath);
            if (statuse != 0) {
                stoponerror(l.layoutnumber, l.num_procs, "Reading input", true);
                goto label_667;
            }
            
            std::string trimmed_chain = trim(l.chain);
            
            if (trimmed_chain == "-i") {
                temp_numnfdes = temp_numnfdes + 1;
                i = i + 1;
                if (temp_numnfdes == 1) {
                    getcommandargument(l.chain2, i, f, l.length, statuse, binaryPath);
                    p = trim(f).length();
                    
                    if ((p - 4) >= 1) {
                        if (f[(p - 4)] == NFDEEXTENSION[0]) {
                            std::string ext = f.substr(p - 4, 4);
                            NFDEEXTENSION = ext;
                            l.extension = NFDEEXTENSION;
                            l.fichin = f.substr(0, p - 5);
                        } else {
                            l.fichin = f.substr(0, p);
                        }
                    } else if (p >= 1) {
                        l.fichin = f.substr(0, p);
                    } else {
                        stoponerror(l.layoutnumber, l.num_procs, "There is not a .nfde file for input", true);
                        statuse = -1;
                        goto label_667;
                    }
                    
                    // Get current working directory
                    std::string cwd = std::filesystem::current_path().string();
                    std::cout << cwd << std::endl;
                    
                    // Inquire file existence
                    std::string full_filename = trim(l.fichin) + NFDEEXTENSION;
                    l.existeNFDE = std::filesystem::exists(full_filename);
                    
                    if (!l.existeNFDE) {
                        buff = "The input file was not found " + trim(l.fichin) + NFDEEXTENSION;
                        stoponerror(l.layoutnumber, l.num_procs, buff, true);
                        statuse = -1;
                        goto label_667;
                    }
                } else if (temp_numnfdes == 2) {
                    l.fichin = "multi_" + trim(l.fichin);
                } else {
                    l.fichin = l.fichin;
                }
            }
            i = i + 1;
        }
    }
    
    // If no input is present we stop
    if (trim(l.fichin).length() <= 0) {
        stoponerror(l.layoutnumber, l.num_procs, "ERROR! -> No input file was specified. Use -i ****.fdtd.json", true);
        statuse = -1;
        goto label_667;
    }

    l.fileFDE = trim(l.fichin) + NFDEEXTENSION;
    l.fileH5 = trim(l.fichin) + ".h5";
    INITWARNINGFILE(l.layoutnumber, l.num_procs, trim(l.fichin) + "_tmpWarnings.txt", l.verbose, l.ignoreerrors);

label_667:
    return;
}

void default_flags_impl(entrada_t& l) {
    l.noconformalmapvtk = false;
    l.forced = -1;
    l.sgbcdepth = -1;
    l.statuse = 0;
    l.time_begin = 0;
    
    l.precision = 0; // redondeo del semiestructurado
    l.stochastic = false;
    l.chosenyesornostochastic = false; // es un flag informativo que debe inicializarse a .false. a pesar de qu el sentido comun diga lo contrario
    l.simu_devia = false;

#ifdef CompileWithHDF
    l.createh5bin = false;
#else
    l.createh5bin = true;
#endif
    l.createh5filefromsinglebin = false;
    l.permitscaling = false;
    l.niapapostprocess = false;
    l.prioritizeCOMPOoverPEC = false; // pec has default more priority than compo (para siva hay que cambiarlo)
    l.prioritizeTHINWIRE = false; // solo para visualizacion y experimentacion 231024
    l.prioritizeISOTROPICBODYoverall = false; // PARA EL SIVA SE CAMBIA POR LINEA DE COMANDO
    l.mpidir = 3; // DEFAULT do NOT ROTATE GEOMETRY !JUST TO TAKE PROFIT OF MPI
    l.maxwireradius = -1.0; // Assuming RKIND is double
    l.boundwireradius = false;
    l.wirecrank = false;
    l.ignoreerrors = false;
    l.ignoresamplingerrors = false;
    l.vtkindex = false; // SOLO AFECTA A LOS VTK (SACA INDICES EN VEZ DE POSICION FISICA)
    l.CLIPREGION = false;
    l.NF2FFDecim = false;
    l.facesNF2FF.tr = true;
    l.facesNF2FF.fr = true;
    l.facesNF2FF.iz = true;
    l.facesNF2FF.de = true;
    l.facesNF2FF.ab = true;
    l.facesNF2FF.ar = true;
    // defaults
    l.read_command_line = true;
    l.hay_slanted_wires = false;
    l.forcing = false;
    l.resume_fromold = false;
    l.singlefilewrite = false;
    l.updateshared = true;

    l.finaltimestep = 0;
    l.cfltemp = 1.0; // dummy
    l.cfl = 1.0; // default courant number !no tocarlo 310715 solo afecta si se usa -l%cfl
}

namespace interpreta_switches_m {

void default_flags(InterpretaSwitches& l) {
    l.forcecfl = false;
    // PML default
    // cpml stretching maximum parameters !!l%alphamaxpar=StaticFrequency*2*pi*Eps0
    l.alphamaxpar = 0.0; // 0.24 !expresion 7.78 taflove 3 edic)
    l.alphaOrden = 1.0;
    l.kappamaxpar = 1.0; // 15.0 !061118 mantener a 1 por conflictos cpml and permittivity scaling
    // and final layer electric sigma
    l.MEDIOEXTRA.exists = false;
    l.MEDIOEXTRA.index = -7; // void
    l.MEDIOEXTRA.pml_size = -1; // void
    l.MEDIOEXTRA.sigma = -1e20; // void
    //
    l.MurAfterPML = false;
    l.mur_second = false;
    l.mur_first = false;
    l.mur_exist = false;
    // !!!!!!!!!!!!
    l.takeintcripte = false; // a peticion de OLD, redondear los nodos de cripte a la baja
    l.attfactorc = 1.0; // default dissipation factor for composites
    l.attfactorw = 1.0; // default dissipation factor for wires

    l.mibc = false;
    l.ade = false; // auxiliary differential equation for composites
    l.conformalskin = false;
    l.sgbc = true; // default is false unless required
    l.sgbcDispersive = false; // default is false unless required
    l.skindepthpre = false;
    l.sgbcdepth = -1; // se calcula automaticamente a menos que se use el switch
    l.sgbcfreq = 1e9; // default es cazar el skin depth hasta 1e9
    l.sgbcresol = 1.0; // numero de celdas por skin depth (caida a exp(-1))
    l.sgbccrank = true; // default es l%sgbccrank

    l.fatalerror = false;
    l.fatalerrornfde2sgg = false;

    l.dontwritevtk = false;

    l.NOcompomur = false; // DEFAULT mi formulacion

    // default no join the wires which are adjacent (ORIGINAL election)
    // do not connect endings unless specified in ORIGINAL
    l.makeholes = false;
    l.connectendings = false;
    l.isolategroupgroups = false;
    l.strictOLD = true; // default is strict ORIGINAL overriden manually
    l.TAPARRABOS = true; // default since 101116 !cortar los multizigzag rabitos
    l.mtlnberenger = true; // solo actua si se invoca con l%wiresflavor berenger esto va a ser siempre true a menos que tambien se invoque con -nomtlnberenger (solo para debugeo y que coincida con Holland) 020719
    l.stableradholland = false; // solo actua si se invoca con l%wiresflavor holland
    l.fieldtotl = false;
    l.experimentalVideal = false;
    l.thereare_stoch = false;
    l.forceresampled = false;
    l.factorradius = 1.0e+30; // para evitar division por cero 120123
    l.factordelta = 1.0e+30; // para evitar division por cero 120123
    // default
    l.groundwires = false;
    l.noSlantedcrecepelo = false; // 131219 experimental niapa ojoooo
    l.inductance_model = "boutayeb";
    l.inductance_order = 8;
    l.wiresflavor = "holland";
    l.wirethickness = 1;
    l.mindistwires = 0.5;
    //
    l.MurAfterPML = false;
    //
    l.createmap = false;
    l.createmapvtk = false;
    l.verbose = false;
    l.saveall = false;
    l.forcesteps = false;
    l.resume = false;
    l.freshstart = false;
    l.run = false; // si hay .fields restartea y si no comienza
    l.deleteintermediates = false;
    //
    l.existeNFDE = false;
    l.existeh5 = false;
    //
    // default is NO flush fields
    l.flushminutesFields = 0;
    // default is to flush data when the buffer is filled up
    // si se pone cada tantos minutos y se guardan las sondas en trancos!puede haber errores de redondeo porque el buffer se limpia tras cada flusheo
    l.flushminutesData = topCPUtime;
    //
    // maximum runtime
    l.maxCPUtime = topCPUtime;
    l.input_conformal_flag = false;
    l.file11isopen = false;
    l.relaunching = false;
    l.forcestop = false;
    l.input_conformal_flag = false;
    l.run_with_dmma = true;
    l.run_with_abrezanjas = false;
    return;
}

} // namespace interpreta_switches_m

