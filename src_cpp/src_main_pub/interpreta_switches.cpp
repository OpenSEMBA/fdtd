```cpp
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <sstream>

// Forward declarations for external modules/types used in the Fortran code
// These would typically be defined in other headers included by this module

// Placeholder for external types and functions referenced in the Fortran code
// These need to be implemented or included from other translation units

namespace FDETYPES_m {
    // Placeholder for type definitions
    using RKIND = float; // Assuming RKIND is single precision float based on context
    using RKIND_wires = float;
    const int BUFSIZE = 256;
    const int BUFSIZE_LONG = 1024;
}

namespace Getargs_m {
    // Placeholder for argument parsing functions
    std::string getBinaryPath();
    int commandargumentcount(const std::string& chaininput, const std::string& binaryPath);
    void getcommandargument(const std::string& chaininput, int i, std::string& result, int& length, int& statuse, const std::string& binaryPath);
}

namespace EpsMuTimeScale_m {
    struct EpsMuTimeScale_input_parameters_t {
        bool electric = false;
        bool magnetic = false;
        double tini = 0.0;
        double tend = 0.0;
        double alpha_max = 0.0;
        
        int checkError() const {
            if (tini < 0.0 || tend <= 0.0 || alpha_max <= 0.0) return -1;
            return 0;
        }
    };
}

namespace Report_m {
    void print11(int layoutnumber, const std::string& msg);
    void stoponerror(int layoutnumber, int num_procs, const std::string& msg, bool fatal);
    void CLOSEWARNINGFILE(int layoutnumber, int num_procs, bool& fatalerror, bool dummy1, bool dummy2);
    void INITWARNINGFILE(int layoutnumber, int num_procs, const std::string& nEntradaRoot, bool verbose, bool ignoreerrors);
    void OffPrint();
    void OnPrint();
}

namespace version_m {
    extern std::string program_name;
    extern std::string compilation_date;
    extern std::string compiler_id;
    extern std::string git_commit;
    extern std::string cmake_build_type;
    extern std::string compilation_flags_debug;
    extern std::string compilation_flags_release;
    extern std::string compilation_flags;
    extern std::string SEPARADOR;
    extern int topCPUtime;
}

// External global variables referenced in the code
extern bool input_conformal_flag;
extern int SUBCOMM_MPI;
extern int MPI_LOGICAL;
extern int MPI_LOR;
extern int MPI_DOUBLE_PRECISION;
extern int MPI_REAL;

// External types referenced in the code
struct nf2ff_T {
    bool TR = true;
    bool FR = true;
    bool IZ = true;
    bool DE = true;
    bool AB = true;
    bool AR = true;
};

struct MedioExtra_t {
    bool exists = false;
    int index = -7;
    double pml_size = -1.0;
    double sigma = -1e20;
};

struct tiempo_t {
    std::string fecha;
    std::string hora;
    double segundos = 0.0;
};

void get_secnds(tiempo_t& time_out2);

namespace interpreta_switches_m {

    struct entrada_t {
        // Logical fields
        bool forcing = false;
        bool singlefilewrite = false;
        bool ignoresamplingerrors = false;
        bool ignoreerrors = false;
        bool updateshared = true;
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
        bool run_with_dmma = true;
        bool run_with_abrezanjas = false;
        bool input_conformal_flag = false;
        bool pausar = false;
        bool relaunching = false;
        bool forcestop = false;
        bool l_aux = false;
        bool flag_conf_sgg = false;
        bool takeintcripte = false;
        bool skindepthpre = false;
        bool sgbc = true;
        bool conformalskin = false;
        bool ade = false;
        bool mibc = false;
        bool NOcompomur = false;
        bool MurAfterPML = false;
        bool sgbccrank = true;
        bool sgbcDispersive = false;
        bool saveall = false;
        bool boundwireradius = false;
        bool hay_slanted_wires = false;
        bool makeholes = false;
        bool mur_first = false;
        bool mur_second = false;
        bool mur_exist = false;
        bool connectendings = false;
        bool strictOLD = true;
        bool mtlnberenger = true;
        bool stableradholland = false;
        bool TAPARRABOS = true;
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
        bool read_command_line = true;

        // Integer fields
        int wirethickness = 1;
        int inductance_order = 8;
        int finaltimestep = 0;
        int ierr = 0;
        int layoutnumber = 0;
        int num_procs = 1;
        int length = 0;
        int mpidir = 3;
        int flushminutesFields = 0;
        int flushminutesData = 0;
        int flushsecondsFields = 0;
        int flushsecondsData = 0;
        int forced = -1;
        int maxCPUtime = 0;
        int sgbcdepth = -1;
        int precision = 0;
        int statuse = 0;

        // Real fields (RKIND)
        double maxwireradius = -1.0;
        double mindistwires = 0.5;
        double attfactorc = 1.0;
        double attfactorw = 1.0;
        double cfltemp = 1.0;
        double cfl = 1.0;
        double sgbcfreq = 1e9;
        double sgbcresol = 1.0;
        double alphamaxpar = 0.0;
        double kappamaxpar = 1.0;
        double alphaOrden = 1.0;

        // Real fields (RKIND=8)
        double time_begin = 0.0;
        double time_end = 0.0;

        // Real fields (RKIND_wires)
        double factorradius = 1.0e+30;
        double factordelta = 1.0e+30;

        // Complex types
        nf2ff_T facesNF2FF;
        MedioExtra_t MEDIOEXTRA;
        EpsMuTimeScale_m::EpsMuTimeScale_input_parameters_t EpsMuTimeScale_input_parameters;
        tiempo_t time_out2;

        // Character strings
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

    void interpreta(entrada_t& l, int& statuse);
    void insertalogtmp(entrada_t& l);
    void print_help(entrada_t& l);
    void print_basic_help(entrada_t& l);
    void print_credits(entrada_t& l);
    void removeintraspaces(std::string& a);
    void buscaswitchficheroinput(entrada_t& l);
    void default_flags(entrada_t& l);

} // namespace interpreta_switches_m

// Implementation of functions

namespace interpreta_switches_m {

    void interpreta(entrada_t& l, int& statuse) {
        std::string chari, f, dubuf, buff, binaryPath;
        bool existiarunningigual = false;
        bool mpidirset = false;
        bool resume3 = false;
        int i, j, donde, n, newmpidir;
        double pausetime;
        int iostatus = 0;

        l.input_conformal_flag = input_conformal_flag;

        mpidirset = false;
        existiarunningigual = false;
        statuse = 0;

        binaryPath = Getargs_m::getBinaryPath();
        n = Getargs_m::commandargumentcount(l.chaininput, binaryPath);
        if (n == 0) {
            print_basic_help(l);
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
            // Simulate trim(adjustl())
            std::string trimmed_opciones = l.opcionestotales;
            trimmed_opciones.erase(0, trimmed_opciones.find_first_not_of(" \t\n\r\f\v"));
            trimmed_opciones.erase(trimmed_opciones.find_last_not_of(" \t\n\r\f\v") + 1);
            
            std::string trimmed_chain = l.chain;
            trimmed_chain.erase(0, trimmed_chain.find_first_not_of(" \t\n\r\f\v"));
            trimmed_chain.erase(trimmed_chain.find_last_not_of(" \t\n\r\f\v") + 1);
            
            l.opcionestotales = trimmed_opciones + " " + trimmed_chain;
        }
        Report_m::print11(l.layoutnumber, "Switches " + l.opcionestotales);

        if (n > 0) {
            i = 2; // Start from 2 because first argument is executable name
            while (i <= n) {
                Getargs_m::getcommandargument(l.chaininput, i, l.chain, l.length, statuse, binaryPath);
                if (statuse != 0) {
                    Report_m::stoponerror(l.layoutnumber, l.num_procs, "Reading input", true);
                    statuse = -1;
                    return;
                }
                
                std::string chain_trimmed = l.chain;
                chain_trimmed.erase(0, chain_trimmed.find_first_not_of(" \t\n\r\f\v"));
                chain_trimmed.erase(chain_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);

                if (chain_trimmed == "-i") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    continue;
                } else if (chain_trimmed == "-a") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    continue;
                } else if (chain_trimmed == "-mpidir") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::string f_trimmed = f;
                    f_trimmed.erase(0, f_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    f_trimmed.erase(f_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    if (f_trimmed == "x" || f_trimmed == "X") {
                        l.mpidir = 1;
                    } else if (f_trimmed == "y" || f_trimmed == "Y") {
                        l.mpidir = 2;
                    } else if (f_trimmed == "z" || f_trimmed == "Z") {
                        l.mpidir = 3;
                    } else {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid or duplicate incoherent -l%mpidir option", true);
                        statuse = -1;
                        continue;
                    }
                    if (!mpidirset) {
                        std::string opciones_trimmed = l.opcionespararesumeo;
                        opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                        opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                        
                        std::string chain_trimmed2 = l.chain;
                        chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                        chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                        
                        std::string f_trimmed2 = f;
                        f_trimmed2.erase(0, f_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                        f_trimmed2.erase(f_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                        
                        l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2 + " " + f_trimmed2;
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
                    // MPI Barrier placeholder
                    // call MPI_Barrier(SUBCOMM_MPI, l%ierr)
                    
                    get_secnds(l.time_out2);
                    l.time_begin = l.time_out2.segundos;
                    dubuf = "Paused for (secs) " + std::to_string(pausetime);
                    Report_m::print11(l.layoutnumber, dubuf);
                    
                    while (l.pausar) {
                        // MPI Barrier placeholder
                        // call MPI_Barrier(SUBCOMM_MPI, l%ierr)
                        
                        get_secnds(l.time_out2);
                        l.time_end = l.time_out2.segundos;
                        if (l.time_end - l.time_begin > pausetime) {
                            l.pausar = false;
                        }
                    }
                    // MPI Barrier placeholder
                    // call MPI_Barrier(SUBCOMM_MPI, l%ierr)
                    // l%l_aux = l%pausar
                    // call MPI_AllReduce(l%l_aux, l%pausar, 1_4, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, l%ierr)
                } else if (chain_trimmed == "-NF2FFDecim") {
                    l.NF2FFDecim = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string f_trimmed2 = f;
                    f_trimmed2.erase(0, f_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    f_trimmed2.erase(f_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2 + " " + f_trimmed2;
                } else if (chain_trimmed == "-noNF2FF") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::string f_trimmed = f;
                    f_trimmed.erase(0, f_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    f_trimmed.erase(f_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    if (f_trimmed == "back" || f_trimmed == "BACK") {
                        l.facesNF2FF.TR = false;
                    } else if (f_trimmed == "front" || f_trimmed == "FRONT") {
                        l.facesNF2FF.FR = false;
                    } else if (f_trimmed == "left" || f_trimmed == "LEFT") {
                        l.facesNF2FF.IZ = false;
                    } else if (f_trimmed == "right" || f_trimmed == "RIGHT") {
                        l.facesNF2FF.DE = false;
                    } else if (f_trimmed == "down" || f_trimmed == "DOWN") {
                        l.facesNF2FF.AB = false;
                    } else if (f_trimmed == "up" || f_trimmed == "UP") {
                        l.facesNF2FF.AR = false;
                    } else {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid -noNF2FF option", true);
                        statuse = -1;
                    }
                    continue;
                } else if (chain_trimmed == "-force") {
                    l.forcing = true;
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss(f);
                    iss >> l.forced;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid cut", true);
                        statuse = -1;
                        goto label_312;
                    }
                    label_312:;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string f_trimmed2 = f;
                    f_trimmed2.erase(0, f_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    f_trimmed2.erase(f_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2 + " " + f_trimmed2;
                } else if (chain_trimmed == "-singlefile") {
                    l.singlefilewrite = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-ignoresamplingerrors") {
                    l.ignoresamplingerrors = true;
                } else if (chain_trimmed == "-prioritizeTHINWIRE") {
                    l.prioritizeTHINWIRE = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                    l.ignoreerrors = true;
                } else if (chain_trimmed == "-prioritizeCOMPOoverPEC") {
                    l.prioritizeCOMPOoverPEC = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                    l.ignoreerrors = true;
                } else if (chain_trimmed == "-noshared") {
                    l.updateshared = false;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-prioritizeISOTROPICBODYoverall") {
                    l.prioritizeISOTROPICBODYoverall = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-wirecrank") {
                    l.wirecrank = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-clip") {
                    l.CLIPREGION = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-verbose") {
                    l.verbose = true;
                } else if (chain_trimmed == "-ignoreerrors") {
                    l.ignoreerrors = true;
                } else if (chain_trimmed == "-r") {
                    l.resume = true;
                    l.forcesteps = true;
                } else if (chain_trimmed == "-cpumax") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss(f);
                    iss >> l.maxCPUtime;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid CPU maximum time", true);
                    }
                    if (l.maxCPUtime <= 0) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid CPU maximum time", true);
                        statuse = -1;
                    }
                } else if (chain_trimmed == "-s") {
                    l.freshstart = true;
                } else if (chain_trimmed == "-flush") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss(f);
                    iss >> l.flushminutesFields;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid flushing interval", true);
                    }
                    if (l.flushminutesFields <= 0) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid flushing interval", true);
                        statuse = -1;
                    }
                } else if (chain_trimmed == "-flushdata") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss(f);
                    iss >> l.flushminutesData;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid flushing interval", true);
                    }
                    if (l.flushminutesData <= 0) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid flushing interval", true);
                        statuse = -1;
                    }
                } else if (chain_trimmed == "-run") {
                    l.run = true;
                } else if (chain_trimmed == "-map") {
                    l.createmap = true;
                } else if (chain_trimmed == "-dontwritevtk") {
                    l.dontwritevtk = true;
                } else if (chain_trimmed == "-vtkindex") {
                    l.vtkindex = true;
                } else if (chain_trimmed == "-mapvtk") {
                    l.createmapvtk = true;
                } else if (chain_trimmed == "-dmma") {
                    l.run_with_dmma = true;
                    l.run_with_abrezanjas = false;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-takeintcripte") {
                    l.takeintcripte = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-mur2") {
                    l.MurAfterPML = true;
                    l.mur_first = true;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-mur1") {
                    l.MurAfterPML = true;
                    l.mur_first = true;
                    l.mur_second = false;
                    std::string opciones_trimmed = l.opcionespararesumeo;
                    opciones_trimmed.erase(0, opciones_trimmed.find_first_not_of(" \t\n\r\f\v"));
                    opciones_trimmed.erase(opciones_trimmed.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    std::string chain_trimmed2 = l.chain;
                    chain_trimmed2.erase(0, chain_trimmed2.find_first_not_of(" \t\n\r\f\v"));
                    chain_trimmed2.erase(chain_trimmed2.find_last_not_of(" \t\n\r\f\v") + 1);
                    
                    l.opcionespararesumeo = opciones_trimmed + " " + chain_trimmed2;
                } else if (chain_trimmed == "-pmlalpha") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss(f);
                    iss >> l.alphamaxpar;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
                        statuse = -1;
                        goto label_8621;
                    }
                    label_8621:;
                    if (l.alphamaxpar < 0.0) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
                        statuse = -1;
                    }
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss2(f);
                    iss2 >> l.alphaOrden;
                    if (iss2.fail()) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML order factor", true);
                        statuse = -1;
                        goto label_8121;
                    }
                    label_8121:;
                    if (l.alphaOrden < 0.0) {
                        Report_m::stoponerror(l.layoutnumber, l.num_procs, "Invalid CPML alpha factor", true);
                        statuse = -1;
                    }
                } else if (chain_trimmed == "-pmlkappa") {
                    i = i + 1;
                    Getargs_m::getcommandargument(l.chaininput, i, f, l.length, statuse, binaryPath);
                    std::istringstream iss(f);
                    iss >> l.kappamaxpar;
                    if (iss.fail()) {
                        Report_m::stoponerror(l.layout