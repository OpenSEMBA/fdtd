#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <ctime>
#include <iomanip>

// Assuming these types and constants are defined elsewhere or need stubs
// Based on the Fortran code, we need to define some placeholders for types not fully visible
// BUFSIZE is likely a constant
#ifndef BUFSIZE
#define BUFSIZE 256
#endif

// RKIND is likely 8 (double)
#ifndef RKIND
#define RKIND 8
#endif

// Types from other modules (FDETYPES_m, snapxdmf_m)
// We will create minimal stubs or assume they are included.
// Since the prompt asks to convert THIS module, we assume dependencies are handled or stubbed.

struct SGGFDTDINFO_t {
    // Placeholder for SGGFDTDINFO_t
    int layoutnumber;
    int num_procs;
    bool resume;
    bool resume_fromold;
    std::string nEntradaRoot;
    std::string nresumeable2;
};

struct sim_control_t {
    // Placeholder for sim_control_t
    int layoutnumber;
    int num_procs;
    bool resume;
    bool resume_fromold;
    std::string nEntradaRoot;
    std::string nresumeable2;
};

struct coorsxyzP_t {
    // Placeholder for coorsxyzP_t
    double x, y, z;
};

struct tiempo_t {
    double segundos;
    char hora[BUFSIZE];
    char fecha[BUFSIZE];
};

// Global constants/variables from other modules
extern const std::string SEPARADOR;
extern const std::string separador;

// MPI stubs if CompileWithMPI is defined
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
#endif

// Helper function stub
coorsxyzP_t Creapuntos(const SGGFDTDINFO_t& sgg) {
    coorsxyzP_t punto;
    punto.x = 0.0;
    punto.y = 0.0;
    punto.z = 0.0;
    return punto;
}

// Helper function stub
void get_secnds(tiempo_t& t) {
    t.segundos = 0.0;
    // In a real implementation, this would get the current time
    std::time_t now = std::time(nullptr);
    t.segundos = static_cast<double>(now);
    
    std::tm* ltime = std::localtime(&now);
    std::strftime(t.hora, BUFSIZE, "%H:%M:%S", ltime);
    std::strftime(t.fecha, BUFSIZE, "%Y-%m-%d", ltime);
}

// Helper matching errorreport.F90 openclosedelete (create, write !END, delete)
void openclosedelete(const std::string& ficherin) {
    int my_iostat = 0;
retry_open:
    if (my_iostat != 0) {
        std::cout << '.' << std::flush;
    }

    std::string filename = ficherin;
    const size_t start = filename.find_first_not_of(' ');
    if (start == std::string::npos) {
        filename.clear();
    } else {
        filename = filename.substr(start);
    }

    {
        std::ofstream file(filename, std::ios::out);
        if (!file.is_open()) {
            my_iostat = 1;
            goto retry_open;
        }
        file << "!END" << std::endl;
        file.close();
    }
    std::remove(filename.c_str());
}

namespace Report_m {

   // Constants
   const int reportingseconds = 60;

   // Global variables
   double time_begin, time_end, time_begin2, time_begin3, time_begin_absoluto, time_end2, time_desdelanzamiento;
   double megaceldas, megaceldastotales, speedInst, speedGlobInst, speedAvg, speedGlobAvg;
   double energy, energyTotal, oldenergyTotal, snapLevel;
   tiempo_t time_out2;
   
   char charmeg[BUFSIZE];
   int reportedinstant, snapStep, snapHowMany, countersnap;
   bool printea = false;
   bool calledStoponerrroonlyprint = false;
   bool warningfileIsOpen = false;
   bool verbose = false;
   bool file10isopen = false;
   bool file11isopen = false;
   char warningFile[BUFSIZE] = " ";
   char whoami[BUFSIZE];

   int thefile; // for mpi file management
   bool ignoreerrors = false;
   
   coorsxyzP_t Punto;

   char mynEntradaRoot[BUFSIZE];

   bool fatalerror = false;
   int CONTADORDEMENSAJES = 0;
   int thefile2; // for mpi file management

   // Function declarations
   void StopOnError(int layoutnumber, int num_procs, const std::string& message, bool calledfrommain = false);
   void InitReporting(SGGFDTDINFO_t& sgg, sim_control_t& c);
   void ReportExistence(bool mur_second, bool MurAfterPML, int layoutnumber, int num_procs, int thereare, double mur_second_val, double MurAfterPML_val);
   void InitTiming();
   void Timing();
   void CloseReportingFiles();
   void print11(int layoutnumber, const std::string& message, bool force = false);
   void OnPrint();
   void OffPrint();
   void WarnErrReport();
   void INITWARNINGFILE();
   void CLOSEWARNINGFILE();
   void get_secnds_global(tiempo_t& t);
   void openfile_mpi();
   void writefile_mpi();
   void closefile_mpi();
   void reportmedia();
   void erasesignalingfiles();
   void openclosedelete_global(const std::string& filename);
   void openclose();
   bool isFatalError();
   void resetFatalError();
   std::string TRIMNULLCHAR(const std::string& str);

   // Implementation

   void OnPrint() {
      printea = true;
   }

   void OffPrint() {
      printea = false;
   }

   void StopOnError(int layoutnumber, int num_procs, const std::string& message, bool calledfrommain) {
      char ficherito[BUFSIZE];
      char whoami_buf[BUFSIZE];
      
      // Format whoami: (layoutnumber+1/num_procs)
      snprintf(whoami_buf, BUFSIZE, "(%d/%d) ", layoutnumber + 1, num_procs);

      std::string full_message = std::string(whoami_buf) + " ERROR: " + message;
      print11(layoutnumber, full_message, true);

#ifdef CompileWithMPI
      int ierr = 0;
#endif

#ifdef keeppause
      if (calledfrommain) {
         if (layoutnumber == 0) {
            std::ofstream pause_file("pause");
            pause_file << "!END";
            pause_file.close();
         }
         print11(layoutnumber, "Trying to relaunch. Correct error, create launch, and remove pause/warning file (or kill the process)", true);
         return;
      } else {
         if (layoutnumber == 0) {
            std::ofstream pause_file("pause");
            pause_file << "!END";
            pause_file.close();
         }
         print11(layoutnumber, "Stopping, but creating the signal file pause to prevent queuing losses!!! (correct error and remove to continue)", true);
      }
#else
      if (layoutnumber == 0) {
         openclosedelete_global("running");
         openclosedelete_global("pause");
         openclosedelete_global("relaunch");
      }
#endif

#ifdef CompileWithMPI
      print11(layoutnumber, "Trying to kill all MPI processes (may fail!)...", true);
      MPI_Abort(SUBCOMM_MPI, -1, &ierr);
      MPI_Finalize();
#endif

      CloseReportingFiles();

      // In C++, we can't STOP 1 directly. We throw an exception or exit.
      // Given the context, exiting is appropriate.
      std::exit(1);
   }

   void CloseReportingFiles() {
      if (file10isopen) {
         // Assuming file 10 is opened via std::ofstream or similar in a real C++ impl
         // Here we just set the flag. In a real app, we'd close the stream.
         file10isopen = false;
      }
      if (file11isopen) {
         file11isopen = false;
      }
   }

   void InitReporting(SGGFDTDINFO_t& sgg, sim_control_t& c) {
#ifdef CompileWithMPI
      int ierr = 0;
#endif
      bool errnofile = false;
      char whoami_buf[BUFSIZE];
      char buff[BUFSIZE];

      snprintf(whoami_buf, BUFSIZE, "(%d/%d) ", c.layoutnumber + 1, c.num_procs);

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

      Punto = Creapuntos(sgg);

      if (c.layoutnumber == 0) {
         std::string filename = std::string(c.nEntradaRoot) + "_Energy.dat";
         // In a real C++ implementation, we would open an ofstream
         // For now, we just set the flag
         file10isopen = true;
      }

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      get_secnds_global(time_out2);

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
      print11(c.layoutnumber, SEPARADOR + separador + separador);

#ifndef CompileWithInt4
#define CompileWithInt4
#endif

      if (c.resume) {
         errnofile = false;
         std::string filename;
         if (c.resume_fromold) {
            filename = std::string(c.nresumeable2) + ".old";
         } else {
            filename = std::string(c.nresumeable2);
         }
         
         // Check if file exists
         std::ifstream test_file(filename);
         if (!test_file.good()) {
            errnofile = true;
         }
         test_file.close();

         if (!errnofile) {
            if (c.resume_fromold) {
               snprintf(buff, BUFSIZE, "FILE %s DOES NOT EXIST", c.nresumeable2.c_str());
               StopOnError(c.layoutnumber, c.num_procs, std::string(buff));
            } else {
               snprintf(buff, BUFSIZE, "FILE %s DOES NOT EXIST", c.nresumeable2.c_str());
               StopOnError(c.layoutnumber, c.num_procs, std::string(buff));
            }
         }
         print11(c.layoutnumber, SEPARADOR + SEPARADOR + SEPARADOR);
         print11(c.layoutnumber, " ");
         if (c.resume_fromold) {
            print11(c.layoutnumber, "Reading resuming data from " + std::string(c.nresumeable2) + ".old etc.");
         } else {
            print11(c.layoutnumber, "Reading resuming data from " + std::string(c.nresumeable2) + " etc.");
         }
         print11(c.layoutnumber, SEPARADOR + SEPARADOR + SEPARADOR);
      }

#ifdef CompileWithMPI
      MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
   }

   void ReportExistence(bool mur_second, bool MurAfterPML, int layoutnumber, int num_procs, int thereare, double mur_second_val, double MurAfterPML_val) {
      // Stub implementation
   }

   void InitTiming() {
      // Stub
   }

   void Timing() {
      // Stub
   }

   void print11(int layoutnumber, const std::string& message, bool force) {
      // Stub implementation
      if (printea || force) {
         std::cout << "[Rank " << layoutnumber << "] " << message << std::endl;
      }
   }

   void WarnErrReport() {
      // Stub
   }

   void INITWARNINGFILE() {
      // Stub
   }

   void CLOSEWARNINGFILE() {
      // Stub
   }

   void get_secnds_global(tiempo_t& t) {
      get_secnds(t);
   }

   void openfile_mpi() {
      // Stub
   }

   void writefile_mpi() {
      // Stub
   }

   void closefile_mpi() {
      // Stub
   }

   void reportmedia() {
      // Stub
   }

   void erasesignalingfiles() {
      // Stub
   }

   void openclosedelete_global(const std::string& filename) {
      openclosedelete(filename);
   }

   void openclose() {
      // Stub
   }

   bool isFatalError() {
      return fatalerror;
   }

   void resetFatalError() {
      fatalerror = false;
   }

   std::string TRIMNULLCHAR(const std::string& str) {
      size_t end = str.find_last_not_of('\0');
      if (end == std::string::npos) {
         return "";
      }
      return str.substr(0, end + 1);
   }

} // namespace Report_m

void ReportExistence(const SGGFDTDINFO_t& sgg, const logic_control_t& thereare, int layoutnumber, int num_procs) {
#ifdef CompileWithMPI
    int ierr;
#endif
    char whoami[BUFSIZE];
    char buff[BUFSIZE];

    snprintf(whoami, BUFSIZE, "(%d/%d) ", layoutnumber + 1, num_procs);

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

    print11(layoutnumber, SEPARADOR + sEPARADOR + SEPARADOR);

    if (thereare.Multiports || thereare.AnisMultiports) {
#ifdef CompileWithNIBC
        // continue
#else
        snprintf(buff, BUFSIZE, "%s MIBC unsupported. Recompile", whoami);
        stoponerror(layoutnumber, num_procs, buff);
#endif
    }

    if (thereare.MagneticMedia) {
        std::string buff_str = " has special H-media";
        warnerrreport(buff_str);
    }
    if (thereare.PMLMagneticMedia) {
        std::string buff_str = " has special PML H-media";
        warnerrreport(buff_str);
    }
    if (thereare.NodalE || thereare.NodalH) {
        std::string buff_str = " has Nodal sources";
        warnerrreport(buff_str);
    }
    if (thereare.Observation) {
        std::string buff_str = " has probes";
        warnerrreport(buff_str);
    }
    if (thereare.FarFields) {
        std::string buff_str = " has Far Field probes";
        warnerrreport(buff_str);
    }
    if (thereare.PlaneWaveBoxes) {
        std::string buff_str = " has planewaves";
        warnerrreport(buff_str);
    }
    if (thereare.Multiports) {
        std::string buff_str = " has MIBC Multiports";
        warnerrreport(buff_str);
    }
    if (thereare.AnisMultiports) {
        std::string buff_str = " has MIBC Anisotropic Multiports";
        warnerrreport(buff_str);
    }
    if (thereare.SGBCs) {
        std::string buff_str = " has Thin metal Materials";
        warnerrreport(buff_str);
    }
    if (thereare.Anisotropic && !thereare.ThinSlot) {
        std::string buff_str = " has pure anisotropic media";
        warnerrreport(buff_str);
    }
    if (thereare.ThinSlot) {
        std::string buff_str = " has Thin Slots";
        warnerrreport(buff_str);
    }

    if (thereare.EDispersives) {
        std::string buff_str = " has electric dispersives";
        warnerrreport(buff_str);
    }
    if (thereare.MDispersives) {
        std::string buff_str = " has magnetic dispersives";
        warnerrreport(buff_str);
    }
    if (thereare.Wires) {
        std::string buff_str = " has Holland WIREs";
        warnerrreport(buff_str);
    }
#ifdef CompileWithBerengerWires
    if (thereare.Wires) {
        std::string buff_str = " has Multi-WIREs";
        warnerrreport(buff_str);
    }
#endif
#ifdef CompileWithSlantedWires
    if (thereare.Wires) {
        std::string buff_str = " has Slanted WIREs";
        warnerrreport(buff_str);
    }
#endif
    if (thereare.PMLBorders) {
        if (sgg.Border.IsUpPML || sgg.Border.IsDownPML) {
            std::string buff_str = " has PML regions inside Z";
            warnerrreport(buff_str);
        }
    }
    if (thereare.MURBorders) {
        if (sgg.Border.IsUpMUR || sgg.Border.IsDownMUR) {
            if (mur_second) {
                std::string buff_str = " has MUR2 regions inside Z";
                warnerrreport(buff_str);
            } else {
                std::string buff_str = " has MUR1 regions inside Z";
                warnerrreport(buff_str);
            }
        }
    }
    if (murAfterPML) {
        if (mur_second) {
            std::string buff_str = " CPML are backed by MUR1";
            warnerrreport(buff_str);
        } else {
            std::string buff_str = " CPML are backed by MUR2";
            warnerrreport(buff_str);
        }
    }
    if (thereare.PMCBorders) {
        std::string buff_str = " has PMC borders";
        warnerrreport(buff_str);
    }
    if (thereare.PECBorders) {
        std::string buff_str = " has PEC borders";
        warnerrreport(buff_str);
    }

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
    print11(layoutnumber, SEPARADOR + sEPARADOR + SEPARADOR);
    warnerrreport(buff);
}

void InitTiming(const SGGFDTDINFO_t& sgg, const sim_control_t& c, double t, int initialtimestep, double maxSourceValue) {
#ifdef CompileWithMPI
    int ierr;
#endif
    char whoami[BUFSIZE];
    char dubuf[BUFSIZE];
    tiempo_t time_out2, time_comienzo;

    snprintf(whoami, BUFSIZE, "(%d/%d) ", c.layoutnumber + 1, c.num_procs);

    time_desdelanzamiento = t;
    snapLevel = 1.0e25 * RKIND;
    snapStep = 1;
    snapHowMany = 1;
    countersnap = 0;

    megaceldas = (1.0 * sgg.sweep[iEx].ZE - 1.0 * sgg.sweep[iEx].ZI) *
                 (1.0 * sgg.sweep[iEx].YE - 1.0 * sgg.sweep[iEx].YI) *
                 (1.0 * sgg.sweep[iEy].XE - 1.0 * sgg.sweep[iEy].XI) / 1.0e6 * RKIND;

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
    MPI_AllReduce(&megaceldas, &megaceldastotales, 1, REALSIZE, MPI_SUM, SUBCOMM_MPI, &ierr);
#else
    megaceldastotales = megaceldas;
#endif

    snprintf(dubuf, BUFSIZE, "Total Mcells: %g", megaceldastotales);
    print11(c.layoutnumber, dubuf);

    if (c.flushsecondsFIELDS != 0) {
        snprintf(dubuf, BUFSIZE, "Flushing restarting FIELDS every %d minutes", static_cast<int>(c.flushsecondsFIELDS / 60.0 * RKIND));
        print11(c.layoutnumber, dubuf);
    } else {
        if (c.maxCPUtime == topCPUtime) {
            print11(c.layoutnumber, "NO flushing of restarting FIELDS scheduled");
        } else {
            snprintf(dubuf, BUFSIZE, "Flushing of restarting FIELDS at the end (mins) :%d", c.maxCPUtime);
            print11(c.layoutnumber, dubuf);
        }
    }

    if (c.flushsecondsDATA != 0) {
        snprintf(dubuf, BUFSIZE, "Flushing observation DATA every  %d minutes and every %d steps",
                 static_cast<int>(c.flushsecondsDATA / 60.0 * RKIND), BuffObse);
        print11(c.layoutnumber, dubuf);
    } else {
        print11(c.layoutnumber, "WARNING: NO flushing of observation DATA scheduled");
    }

    snprintf(dubuf, BUFSIZE, "Reporting simulation info every  %d minutes ", static_cast<int>(reportingseconds / 60.0 * RKIND));
    print11(c.layoutnumber, dubuf);

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
    get_secnds(time_out2);
    print11(c.layoutnumber, SEPARADOR + separador + separador);

    snprintf(dubuf, BUFSIZE, "Simulation from n=%d, t=%e, to n=%d, t=%e",
             initialtimestep, sgg.tiempo[initialtimestep],
             c.finaltimestep, sgg.tiempo[c.finaltimestep]);
    print11(c.layoutnumber, dubuf);

    snprintf(dubuf, BUFSIZE, "Date/time %s/%s/%s   %s:%s:%s",
             time_out2.fecha.substr(6, 2).c_str(),
             time_out2.fecha.substr(4, 2).c_str(),
             time_out2.fecha.substr(0, 4).c_str(),
             time_out2.hora.substr(0, 2).c_str(),
             time_out2.hora.substr(2, 2).c_str(),
             time_out2.hora.substr(4, 2).c_str());
    print11(c.layoutnumber, dubuf);

    time_begin_absoluto = time_out2.segundos;
    time_begin = time_begin_absoluto;
    time_begin2 = time_begin_absoluto;
    time_begin3 = time_begin_absoluto;

    reportedinstant = initialtimestep;

#ifdef CompileWithMPI
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
}

void get_secnds(tiempo_t& time_out2) {
    int diasen[12] = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
    int diasenbisiesto[12]; // Not fully used in this snippet, but declared
    // Note: Fortran data statements initialize arrays. In C++, we initialize explicitly.
    // The snippet ends here, so we just provide the structure.
}

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <iomanip>

// Assuming necessary types and constants are defined elsewhere
// e.g., BUFSIZE, RKIND, SGGFDTDINFO_t, bounds_t, perform_t, tiempo_t
// For the purpose of this translation, we assume standard C++ equivalents or placeholders.

namespace FortranModule {

    // Constants
    const int BUFSIZE = 256; // Placeholder, adjust as needed
    const double RKIND = 8.0; // Placeholder for kind=8 or default real

    // Helper function to simulate date_and_time
    void get_date_and_time(std::string& date, std::string& time, std::string& zone, std::vector<int>& values) {
        // This is a placeholder. In a real scenario, use std::chrono or similar.
        // For now, we just return empty/default strings to satisfy the interface.
        date = "20231001";
        time = "120000.000";
        zone = "UTC";
        values.resize(8, 0);
    }

    // Global save variable simulation
    double t_0 = 0.0;
    bool esprimeravez = true;

    // Dias en bisiesto array
    const std::vector<int> diasenbisiesto = {31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};
    // Dias no bisiesto (implied by usage in else branch, though not explicitly defined in snippet, 
    // usually diasenbisiesto(i) - 1 for non-leap or specific array. 
    // Based on logic: diasen(month-1) is used. Let's assume diasen is diasenbisiesto adjusted or separate.
    // Common pattern: diasen(i) = diasenbisiesto(i) - 1 for i>1? Or just a separate array.
    // Since it's not defined, we'll create a placeholder or assume it's diasenbisiesto with leap day removed.
    // Standard days: 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31. Cumulative: 0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334.
    // The snippet uses diasenbisiesto(month-1). If month=1, index 0 -> 31? No, usually cumulative days BEFORE the month.
    // Let's assume diasen is the cumulative days for non-leap year.
    const std::vector<int> diasen = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

    struct tiempo_t {
        double segundos;
        std::string hora;
        std::string fecha;
    };

    void get_secnds(tiempo_t& time_out2) {
        double time_out;

        // Local variables
        double s;
        // t_0 is static/save in Fortran, handled by global in namespace or static local
        static double t_0_local = 0.0; 
        
        int h, m, month, day, year, cent;
        
        std::string caux(BUFSIZE, ' ');
        std::string caux2(BUFSIZE, ' ');
        std::string zone(5, ' ');
        std::vector<int> values(8, 0);

        // Call date_and_time
        get_date_and_time(caux2, caux, zone, values);

        // Parse time
        // caux format: HHMMSS.SSS
        if (caux.length() >= 2) h = std::stoi(caux.substr(0, 2));
        else h = 0;
        
        if (caux.length() >= 4) m = std::stoi(caux.substr(2, 2));
        else m = 0;
        
        if (caux.length() >= 9) s = std::stod(caux.substr(4, 6)); // f6.3
        else s = 0.0;

        // Parse date
        // caux2 format: YYYYMMDD
        if (caux2.length() >= 4) year = std::stoi(caux2.substr(0, 4));
        else year = 2000;
        
        if (caux2.length() >= 6) month = std::stoi(caux2.substr(4, 2));
        else month = 1;
        
        if (caux2.length() >= 8) day = std::stoi(caux2.substr(6, 2));
        else day = 1;

        // Logic for t_0 initialization (simplified from commented out block)
        // The original code had a complex initialization block commented out.
        // The active code uses diasenbisiesto/diasen directly in the calculation.
        // The t_0 variable is subtracted. If esprimeravez was true, t_0 would be set.
        // Since the initialization block is commented out, t_0 remains 0.0 (or initial value).
        // We will follow the active code path.

        if ((year % 4 == 0) && (year != 0)) {
            // Leap year
            // diasenbisiesto is 1-indexed in Fortran usually, but here array is 0-indexed in C++.
            // Fortran: diasenbisiesto(month-1). If month=1, index 0. Value 31?
            // Wait, diasenbisiesto defined as /31,60.../. 
            // If month=1, we want days before Jan? 0. 
            // If month=2, we want days before Feb? 31.
            // The array starts with 31. So diasenbisiesto(1) in Fortran is 31.
            // In C++, diasenbisiesto[0] is 31.
            // So diasenbisiesto(month-1) in Fortran maps to diasenbisiesto[month-1] in C++?
            // If month=1, Fortran index 0 -> 31. This seems wrong for "days before month".
            // Usually cumulative days: 0, 31, 59...
            // Let's look at the array: 31, 60, 91...
            // 31 is days in Jan. 60 is Jan+Feb.
            // So diasenbisiesto(i) is days in first i months?
            // If month=1, we want 0 days before.
            // If month=2, we want 31 days before.
            // The code uses diasenbisiesto(month-1).
            // If month=1, index 0 -> 31. This implies it's adding days OF the previous month?
            // No, the formula is: diasenbisiesto(month-1) * 86400 + (day-1)*86400...
            // If month=1, day=1. Result: 31*86400 + 0. This is Feb 1st?
            // This suggests the array might be 1-based in Fortran logic but stored 0-based?
            // Or the array represents cumulative days AT the end of the month?
            // Let's assume the Fortran array `diasenbisiesto` is 1-indexed.
            // diasenbisiesto(1) = 31.
            // diasenbisiesto(2) = 60.
            // If month=1, Fortran code accesses diasenbisiesto(0)? No, Fortran arrays are 1-based by default.
            // If the declaration is `data diasenbisiesto /.../`, it's 1-based.
            // So `diasenbisiesto(month-1)` when month=1 accesses index 0? That's out of bounds in standard Fortran unless lower bound is 0.
            // However, if the array is 1-based, `month-1` for month=1 is 0.
            // Maybe the array is 0-based in the source? Or maybe `month` starts from 2?
            // Let's look at `diasen(month-1)`.
            // If the array is 1-based, accessing index 0 is invalid.
            // Perhaps the array is defined as `diasenbisiesto(0:11)`?
            // Let's assume the C++ vector is 0-indexed and matches the Fortran logical indexing if the Fortran array was 0-based.
            // If Fortran array is 1-based, `diasenbisiesto(month-1)` is a bug or implies 0-based Fortran.
            // Given `diasenbisiesto` has 12 elements, and months are 1-12.
            // If month=1, we want 0. If month=2, we want 31.
            // The array provided: 31, 60, 91...
            // Index 0: 31. Index 1: 60.
            // If we use `diasenbisiesto[month-2]` for month>=2?
            // Let's assume the Fortran code works and the array is 1-based, so `month-1` is valid for month>=1?
            // No, 1-1=0.
            // Maybe the array is `diasenbisiesto(1:12)` and the code uses `diasenbisiesto(month-1)` which is wrong?
            // Or maybe the array is `diasenbisiesto(0:11)`?
            // Let's assume the C++ vector `diasenbisiesto` corresponds to the Fortran array.
            // If Fortran is 1-based, `diasenbisiesto(1)` is 31.
            // If the code uses `diasenbisiesto(month-1)`, and month=1, it accesses 0.
            // This suggests the Fortran array might be 0-based or the code has an off-by-one error that we must preserve.
            // However, `diasen` is used in the else branch.
            // Let's assume the standard interpretation: The array contains cumulative days.
            // For non-leap: 0, 31, 59, 90...
            // The provided array starts with 31.
            // It is highly likely that the Fortran array is 1-based and the code `month-1` is a mistake in my understanding or the source.
            // BUT, I must PRESERVE NAMES and LOGIC.
            // If I translate `diasenbisiesto(month-1)` to `diasenbisiesto[month-1]`, and `month=1`, I access index 0 (31).
            // This will be preserved.

            time_out = diasenbisiesto[month-1] * 86400.0 + (day - 1) * 86400.0 + 3600.0 * h + 60.0 * m + s - t_0_local + (year - 2000.0) * 365.0 * 86400.0;
        } else {
            time_out = diasen[month-1] * 86400.0 + (day - 1) * 86400.0 + 3600.0 * h + 60.0 * m + s - t_0_local + (year - 2000.0) * 365.0 * 86400.0;
        }

        time_out2.segundos = time_out;
        time_out2.hora = caux;
        time_out2.fecha = caux2;
    }

    // Placeholder for SGGFDTDINFO_t, bounds_t, perform_t
    struct SGGFDTDINFO_t {};
    struct bounds_t {
        struct { int NX, NY, NZ; } Ex, Ey, Ez;
    };
    struct perform_t {
        void reset() {}
    };

    // Global variables for Timing
    double time_begin = 0.0;
    double time_end = 0.0;
    double reportingseconds = 10.0; // Placeholder

    void Timing(
        const SGGFDTDINFO_t& sgg,
        const bounds_t& b,
        int n,
        int& n_info,
        int layoutnumber,
        int num_procs,
        int maxCPUtime,
        int flushsecondsFields,
        int flushsecondsData,
        int initialtimestep,
        int finaltimestep,
        perform_t& perform,
        bool& parar,
        bool forcetiming,
        const std::vector<double>& Ex,
        const std::vector<double>& Ey,
        const std::vector<double>& Ez,
        bool everflushed,
        int nentradaroot,
        double maxSourceValue,
        const std::string& opcionestotales,
        bool simu_devia,
        bool dontwritevtk,
        bool permitscaling
    ) {
        // Local variables
        bool simu_devia_local = simu_devia;
        bool dontwritevtk_local = dontwritevtk;
        bool stopdontwritevtk = false;
        bool stopflushingdontwritevtk = false;
        bool flushdontwritevtk = false;
        bool stoponlydontwritevtk = false;

        // Inputs
        // sgg, b, opcionestotales, layoutnumber, num_procs, n, maxCPUtime
        // flushsecondsFields, flushsecondsData, initialtimestep, finaltimestep
        // Ex, Ey, Ez
        // forcetiming, everflushed, permitscaling
        // nEntradaRoot (string)
        
        std::string nEntradaRoot = "root"; // Placeholder

        // I/O
        int my_iostat = 0;
        bool hay_timing = false;
        bool l_aux = false;
        bool mustflushFIELDS = false;
        bool mustflushDATA = false;
        bool mustUnpack = false;
        bool mustPostprocess = false;
        bool mustflushXdmf = false;
        bool mustflushVTK = false;
        bool pararflushing = false;
        bool pararNOflushing = false;
        bool stoponNaN = false;
        bool stoponNaN_aux = false;
        bool mustSnap = false;
        bool stop_only = false;
        bool stopflushing_only = false;
        bool flush_only = false;
        bool flushdata_only = false;
        bool stopflushingonlydontwritevtk = false;
        bool flushonlydontwritevtk = false;
        bool flushdataonlydontwritevtk = false;
        bool flushdatadontwritevtk = false;

        int in_aux, ini_i, fin_i, ini_j, fin_j, ini_k, fin_k, i, j, k;
        std::string whoamishort(BUFSIZE, ' ');
        std::string whoami(BUFSIZE, ' ');
        std::string chinstant(BUFSIZE, ' ');
        std::string dubuf(BUFSIZE, ' ');
        std::string dondex(BUFSIZE, ' ');
        std::string dondey(BUFSIZE, ' ');
        std::string dondez(BUFSIZE, ' ');

        std::vector<double> NEWlmaxval(num_procs, 0.0);
        std::vector<double> NEWlmaxval_x(num_procs, 0.0);
        std::vector<double> NEWlmaxval_y(num_procs, 0.0);
        std::vector<double> NEWlmaxval_z(num_procs, 0.0);
        std::vector<int> NEWlmaxval_i(num_procs, 0);
        std::vector<int> NEWlmaxval_j(num_procs, 0);
        std::vector<int> NEWlmaxval_k(num_procs, 0);

        std::vector<double> lmaxval(num_procs, 0.0);
        std::vector<double> lmaxval_x(num_procs, 0.0);
        std::vector<double> lmaxval_y(num_procs, 0.0);
        std::vector<double> lmaxval_z(num_procs, 0.0);
        std::vector<int> lmaxval_i(num_procs, 0);
        std::vector<int> lmaxval_j(num_procs, 0);
        std::vector<int> lmaxval_k(num_procs, 0);

        double qmaxval = 0.0;
        double qmaxval_x = 0.0;
        double qmaxval_y = 0.0;
        double qmaxval_z = 0.0;
        int qmaxval_i = 0;
        int qmaxval_j = 0;
        int qmaxval_k = 0;
        int thefilenoflu = 0;

        std::vector<double> NEWlminval(num_procs, 0.0);
        std::vector<double> NEWlminval_x(num_procs, 0.0);
        std::vector<double> NEWlminval_y(num_procs, 0.0);
        std::vector<double> NEWlminval_z(num_procs, 0.0);
        std::vector<int> NEWlminval_i(num_procs, 0);
        std::vector<int> NEWlminval_j(num_procs, 0);
        std::vector<int> NEWlminval_k(num_procs, 0);

        std::vector<double> lminval(num_procs, 0.0);
        std::vector<double> lminval_x(num_procs, 0.0);
        std::vector<double> lminval_y(num_procs, 0.0);
        std::vector<double> lminval_z(num_procs, 0.0);
        std::vector<int> lminval_i(num_procs, 0);
        std::vector<int> lminval_j(num_procs, 0);
        std::vector<int> lminval_k(num_procs, 0);

        double qminval = 0.0;
        double qminval_x = 0.0;
        double qminval_y = 0.0;
        double qminval_z = 0.0;
        int qminval_i = 0;
        int qminval_j = 0;
        int qminval_k = 0;
        int dimxsnap = 0;
        int dimysnap = 0;
        int dimzsnap = 0;
        int veces = 0;
        int i1 = 0, j1 = 0, k1 = 0;
        int ini_ibox = 0, fin_ibox = 0;
        int ini_jbox = 0, fin_jbox = 0;
        int ini_kbox = 0, fin_kbox = 0;
        int ini_iboxsin = 0, fin_iboxsin = 0;
        int ini_jboxsin = 0, fin_jboxsin = 0;
        int ini_kboxsin = 0, fin_kboxsin = 0;

        std::string ficherito(BUFSIZE, ' ');

        // MPI placeholder
        int ierr = 0;

        // Output
        // perform is passed by reference

        // Implementation
        whoami = "(" + std::to_string(layoutnumber + 1) + "/" + std::to_string(num_procs) + ") ";
        whoamishort = std::to_string(layoutnumber + 1);

#ifdef CompileWithMPI
        // MPI_Barrier placeholder
        // call MPI_Barrier(MPI_COMM_WORLD,ierr)
#endif

        tiempo_t time_out2;
        get_secnds(time_out2);
        time_end = time_out2.segundos;

        l_aux = ((time_end - time_begin) > reportingseconds) || forcetiming;

#ifdef CompileWithMPI
        // MPI_AllReduce placeholder
        // hay_timing = l_aux; // Simplified
#endif
        hay_timing = l_aux;

        n_info = n + 5;

        perform.reset();
        mustflushFIELDS = false;
        mustflushDATA = false;
        mustUnpack = false;
        mustPostprocess = false;
        mustflushXdmf = false;
        mustflushVTK = false;
    }

} // namespace FortranModule

energy = 0.0;
        //--->
        if (hay_timing) { //no calculation of time until at least 300 seconds lapse
            perform.flushFIELDS = false;
            perform.flushDATA = false;
            if (std::abs(time_end - time_begin_absoluto) < 1.0) time_end = 60.0 + time_begin_absoluto;
            if (std::abs(time_end - time_begin) < 1.0) time_end = 60.0 + time_begin;
            speedInst = ((N - reportedinstant + 1) * megaceldas / (time_end - time_begin));
            speedAvg = ((N - INITIALtimeSTEP + 1) * megaceldas / (time_end - time_begin_absoluto));
            if (speedAvg == 0) speedAvg = 100.0;
            if (speedInst == 0) speedInst = speedAvg;
#ifdef CompileWithMPI
            //print *,'layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg pre', &
            //         layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg
            call_MPI_AllReduce(speedInst, speedGlobInst, 1, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr);
            call_MPI_Barrier(SUBCOMM_MPI, ierr);
            call_MPI_AllReduce(speedAvg, speedGlobAvg, 1, REALSIZE, MPI_SUM, SUBCOMM_MPI, ierr);
            call_MPI_Barrier(SUBCOMM_MPI, ierr);
            //print *,'layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg post', &
            //         layoutnumber+1,speedInst, speedGlobInst,speedAvg, speedGlobAvg
#else
            speedGlobInst = speedInst;
            speedGlobAvg = speedAvg;
#endif
            //
            in_aux = n + std::max(static_cast<int>((reportingseconds / (megaceldastotales / speedGlobInst))) + 1, 1);
#ifdef CompileWithMPI
            //print *,'layoutnumber+1,in_aux, n_info pre',layoutnumber+1,in_aux, n_info
            call_MPI_AllReduce(in_aux, n_info, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr);
            //print *,'layoutnumber+1,in_aux, n_info post',layoutnumber+1,in_aux, n_info
#else
            n_info = in_aux;
#endif
            std::ifstream stop_file("stop");
            pararNOflushing = stop_file.good();
            stop_file.close();
            
            std::ifstream stop_only_file("stop_only");
            stop_only = stop_only_file.good();
            stop_only_file.close();
            
            std::ifstream stop_dontwritevtk_file("stop_dontwritevtk");
            stopdontwritevtk = stop_dontwritevtk_file.good();
            stop_dontwritevtk_file.close();
            
            std::ifstream stop_only_dontwritevtk_file("stop_only_dontwritevtk");
            stoponlydontwritevtk = stop_only_dontwritevtk_file.good();
            stop_only_dontwritevtk_file.close();

            if (pararnoflushing) {
                dontwritevtk = false;
            }
            if (stopdontwritevtk) {
                pararnoflushing = true;
                dontwritevtk = true;
            }
            if (stoponlydontwritevtk) {
                stop_only = true;
                dontwritevtk = true;
            }
            if (stop_only) {
                std::ifstream thefilenoflu("stop_only", std::ios::in);
                std::string quien_es;
                thefilenoflu >> quien_es;
                thefilenoflu.close();
                if (trim(adjustl(quien_es)) == trim(adjustl(nentradaroot))) {
                    pararNOflushing = true;
                } else {
                    pararNOflushing = false;
                }
            }
            //
            std::ifstream stopflushing_file("stopflushing");
            pararflushing = stopflushing_file.good();
            stopflushing_file.close();
            
            std::ifstream stopflushing_only_file("stopflushing_only");
            stopflushing_only = stopflushing_only_file.good();
            stopflushing_only_file.close();
            
            std::ifstream stopflushing_dontwritevtk_file("stopflushing_dontwritevtk");
            stopflushingdontwritevtk = stopflushing_dontwritevtk_file.good();
            stopflushing_dontwritevtk_file.close();
            
            std::ifstream stopflushing_only_dontwritevtk_file("stopflushing_only_dontwritevtk");
            stopflushingonlydontwritevtk = stopflushing_only_dontwritevtk_file.good();
            stopflushing_only_dontwritevtk_file.close();

            if (pararflushing) {
                dontwritevtk = false;
            }
            if (stopflushingdontwritevtk) {
                pararflushing = true;
                dontwritevtk = true;
            }
            if (stopflushingonlydontwritevtk) {
                stopflushing_only = true;
                dontwritevtk = true;
            }
            if (stopflushing_only) {
                std::ifstream thefilenoflu("stopflushing_only", std::ios::in);
                std::string quien_es;
                thefilenoflu >> quien_es;
                thefilenoflu.close();
                if (trim(adjustl(quien_es)) == trim(adjustl(nentradaroot))) {
                    pararflushing = true;
                } else {
                    pararflushing = false;
                }
            }
            //        
            std::ifstream flush_file("flush");
            mustflushFIELDS = flush_file.good();
            flush_file.close();
            
            std::ifstream flush_only_file("flush_only");
            flush_only = flush_only_file.good();
            flush_only_file.close();
            
            std::ifstream flush_dontwritevtk_file("flush_dontwritevtk");
            flushdontwritevtk = flush_dontwritevtk_file.good();
            flush_dontwritevtk_file.close();
            
            std::ifstream flush_only_dontwritevtk_file("flush_only_dontwritevtk");
            flushonlydontwritevtk = flush_only_dontwritevtk_file.good();
            flush_only_dontwritevtk_file.close();

            if (mustflushFIELDS) {
                dontwritevtk = false;
            }
            if (flushdontwritevtk) {
                mustflushFIELDS = true;
                dontwritevtk = true;
            }
            if (flushdontwritevtk) { // Note: Original code has duplicate check, preserving logic
                flush_only = true;
                dontwritevtk = true;
            }
            if (flush_only) {
                std::ifstream thefilenoflu("flush_only", std::ios::in);
                std::string quien_es;
                thefilenoflu >> quien_es;
                thefilenoflu.close();
                if (trim(adjustl(quien_es)) == trim(adjustl(nentradaroot))) {
                    mustflushFIELDS = true;
                } else {
                    mustflushFIELDS = false;
                }
            }
            //
            std::ifstream flushdata_file("flushdata");
            mustflushdata = flushdata_file.good();
            flushdata_file.close();
            
            std::ifstream flushdata_only_file("flushdata_only");
            flushdata_only = flushdata_only_file.good();
            flushdata_only_file.close();
            
            std::ifstream flushdata_dontwritevtk_file("flushdata_dontwritevtk");
            flushdatadontwritevtk = flushdata_dontwritevtk_file.good();
            flushdata_dontwritevtk_file.close();
            
            std::ifstream flushdata_only_dontwritevtk_file("flushdata_only_dontwritevtk");
            flushdataonlydontwritevtk = flushdata_only_dontwritevtk_file.good();
            flushdata_only_dontwritevtk_file.close();

            if (mustflushdata) {
                dontwritevtk = false;
            }
            if (flushdatadontwritevtk) {
                mustflushdata = true;
                dontwritevtk = true;
            }
            if (flushdataonlydontwritevtk) {
                flushdata_only = true;
                dontwritevtk = true;
            }
            if (flushdata_only) {
                std::ifstream thefilenoflu("flushdata_only", std::ios::in);
                std::string quien_es;
                thefilenoflu >> quien_es;
                thefilenoflu.close();
                if (trim(adjustl(quien_es)) == trim(adjustl(nentradaroot))) {
                    mustflushdata = true;
                } else {
                    mustflushdata = false;
                }
            }
            //
            std::ifstream unpack_file("unpack");
            mustUnpack = unpack_file.good();
            unpack_file.close();
            
            std::ifstream postprocess_file("postprocess");
            mustPostprocess = postprocess_file.good();
            postprocess_file.close();
            
            std::ifstream flushxdmf_file("flushxdmf");
            mustflushXdmf = flushxdmf_file.good();
            flushxdmf_file.close();
            
            std::ifstream flushvtk_file("flushvtk");
            mustflushVTK = flushvtk_file.good();
            flushvtk_file.close();

            pararflushing = pararflushing || (std::ceil((time_end - time_desdelanzamiento) / 60.0) >= maxCPUtime);
#ifdef CompileWithMPI
            l_aux = dontwritevtk;
            call_MPI_AllReduce(l_aux, dontwritevtk, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = pararnoflushing;
            call_MPI_AllReduce(l_aux, pararnoflushing, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = pararflushing;
            call_MPI_AllReduce(l_aux, pararflushing, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = mustUnpack;
            call_MPI_AllReduce(l_aux, mustUnpack, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = mustPostprocess;
            call_MPI_AllReduce(l_aux, mustPostprocess, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = mustflushXdmf;
            call_MPI_AllReduce(l_aux, mustflushXdmf, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = mustflushVTK;
            call_MPI_AllReduce(l_aux, mustflushVTK, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = mustflushFIELDS;
            call_MPI_AllReduce(l_aux, mustflushFIELDS, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            //
            l_aux = mustflushdata;
            call_MPI_AllReduce(l_aux, mustflushdata, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr);
            call_MPI_Barrier(MPI_COMM_WORLD, ierr); //050619 cambiado subcomm a mpi_comm_world
#endif
            mustflushfields = pararflushing || mustflushfields;
            parar = pararNOflushing || pararflushing;
            mustflushdata = mustflushdata || mustflushfields || parar; //data is enforced to flush if fields are flushed
            //
            //to prevent duplicate writes on resuming
            //--->
            for (int i = 0; i < num_procs; ++i) {
                lmaxval[i] = 0.0;
                lmaxval_i[i] = 0;
                lmaxval_j[i] = 0;
                lmaxval_k[i] = 0;
            }
        }

lmaxval_x.assign(num_procs, 0.0);
        lmaxval_y.assign(num_procs, 0.0);
        lmaxval_z.assign(num_procs, 0.0);
        lminval.assign(num_procs, 0.0);
        lminval_i.assign(num_procs, 0);
        lminval_j.assign(num_procs, 0);
        lminval_k.assign(num_procs, 0);
        lminval_x.assign(num_procs, 1e+20);
        lminval_y.assign(num_procs, 1e+20);
        lminval_z.assign(num_procs, 1e+20);

        valor = 0.0;
        //--->
        ini_i = b.sweepSINPMLEx.XI - b.Ex.XI;
        fin_i = b.sweepSINPMLEx.XE - b.Ex.XI;
        ini_j = b.sweepSINPMLEx.YI - b.Ex.YI;
        fin_j = b.sweepSINPMLEx.YE - b.Ex.YI;
        ini_k = b.sweepSINPMLEx.ZI - b.Ex.ZI;
        fin_k = b.sweepSINPMLEx.ZE - b.Ex.ZI;

        for (k = ini_k; k <= fin_k; ++k) {
            for (j = ini_j; j <= fin_j; ++j) {
                for (i = ini_i; i <= fin_i; ++i) {
                    valor = valor + Ex(i, j, k) * Ex(i, j, k);
                }
            }
        }
        //--->
        ini_i = b.sweepSINPMLEy.XI - b.Ey.XI;
        fin_i = b.sweepSINPMLEy.XE - b.Ey.XI;
        ini_j = b.sweepSINPMLEy.YI - b.Ey.YI;
        fin_j = b.sweepSINPMLEy.YE - b.Ey.YI;
        ini_k = b.sweepSINPMLEy.ZI - b.Ey.ZI;
        fin_k = b.sweepSINPMLEy.ZE - b.Ey.ZI;
        for (k = ini_k; k <= fin_k; ++k) {
            for (j = ini_j; j <= fin_j; ++j) {
                for (i = ini_i; i <= fin_i; ++i) {
                    valor = valor + Ey(i, j, k) * Ey(i, j, k);
                }
            }
        }
        //--->
        ini_i = b.sweepSINPMLEz.XI - b.Ez.XI;
        fin_i = b.sweepSINPMLEz.XE - b.Ez.XI;
        ini_j = b.sweepSINPMLEz.YI - b.Ez.YI;
        fin_j = b.sweepSINPMLEz.YE - b.Ez.YI;
        ini_k = b.sweepSINPMLEz.ZI - b.Ez.ZI;
        fin_k = b.sweepSINPMLEz.ZE - b.Ez.ZI;
        for (k = ini_k; k <= fin_k; ++k) {
            for (j = ini_j; j <= fin_j; ++j) {
                for (i = ini_i; i <= fin_i; ++i) {
                    valor = valor + Ez(i, j, k) * Ez(i, j, k);
                }
            }
        }
        //
        //--->
        energy = valor; // quitado 241018 para evitar pasar el eps0----> 0.5 * Eps0 * valor
        //--->
        energytotal = energy;
#ifdef CompileWithMPI
        //print *,'layoutnumber+1,energy,energytotal pre',layoutnumber+1,energy,energytotal
        MPI_AllReduce(&energy, &energyTotal, 1, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        //print *,'layoutnumber+1,energy,energytotal post ',layoutnumber+1,energy,energytotal
#else
        energytotal = energy;
#endif

        ini_iboxsin = std::max(std::max(b.sweepSINPMLEx.XI - b.Ex.XI, b.sweepSINPMLEy.XI - b.Ey.XI), b.sweepSINPMLEz.XI - b.Ez.XI);
        fin_iboxsin = std::min(std::min(b.sweepSINPMLEx.XE - b.Ex.XI, b.sweepSINPMLEy.XE - b.Ey.XI), b.sweepSINPMLEz.XE - b.Ez.XI);
        ini_jboxsin = std::max(std::max(b.sweepSINPMLEx.YI - b.Ex.YI, b.sweepSINPMLEy.YI - b.Ey.YI), b.sweepSINPMLEz.YI - b.Ez.YI);
        fin_jboxsin = std::min(std::min(b.sweepSINPMLEx.YE - b.Ex.YI, b.sweepSINPMLEy.YE - b.Ey.YI), b.sweepSINPMLEz.YE - b.Ez.YI);
        ini_kboxsin = std::max(std::max(b.sweepSINPMLEx.ZI - b.Ex.ZI, b.sweepSINPMLEy.ZI - b.Ey.ZI), b.sweepSINPMLEz.ZI - b.Ez.ZI);
        fin_kboxsin = std::min(std::min(b.sweepSINPMLEx.ZE - b.Ex.ZI, b.sweepSINPMLEy.ZE - b.Ey.ZI), b.sweepSINPMLEz.ZE - b.Ez.ZI);
        ini_ibox = std::max(std::max(b.sweepEx.XI - b.Ex.XI, b.sweepEy.XI - b.Ey.XI), b.sweepEz.XI - b.Ez.XI);
        fin_ibox = std::min(std::min(b.sweepEx.XE - b.Ex.XI, b.sweepEy.XE - b.Ey.XI), b.sweepEz.XE - b.Ez.XI);
        ini_jbox = std::max(std::max(b.sweepEx.YI - b.Ex.YI, b.sweepEy.YI - b.Ey.YI), b.sweepEz.YI - b.Ez.YI);
        fin_jbox = std::min(std::min(b.sweepEx.YE - b.Ex.YI, b.sweepEy.YE - b.Ey.YI), b.sweepEz.YE - b.Ez.YI);
        ini_kbox = std::max(std::max(b.sweepEx.ZI - b.Ex.ZI, b.sweepEy.ZI - b.Ey.ZI), b.sweepEz.ZI - b.Ez.ZI);
        fin_kbox = std::min(std::min(b.sweepEx.ZE - b.Ex.ZI, b.sweepEy.ZE - b.Ey.ZI), b.sweepEz.ZE - b.Ez.ZI);
        for (k = ini_kbox; k <= fin_kbox; ++k) {
            for (j = ini_jbox; j <= fin_jbox; ++j) {
                for (i = ini_ibox; i <= fin_ibox; ++i) {
                    valor = std::sqrt(Ex(i, j, k) * Ex(i, j, k) + Ey(i, j, k) * Ey(i, j, k) + Ez(i, j, k) * Ez(i, j, k));
                    if (lmaxval(layoutnumber + 1) < valor) {
                        lmaxval(layoutnumber + 1) = valor;
                        lmaxval_i(layoutnumber + 1) = i + b.Hx.XI;
                        lmaxval_j(layoutnumber + 1) = j + b.Hy.YI;
                        lmaxval_k(layoutnumber + 1) = k + b.Hz.ZI;
                        lmaxval_x(layoutnumber + 1) = Punto.PhysCoor(iHx).x(lmaxval_i(layoutnumber + 1));
                        lmaxval_y(layoutnumber + 1) = Punto.PhysCoor(iHy).y(lmaxval_j(layoutnumber + 1));
                        lmaxval_z(layoutnumber + 1) = Punto.PhysCoor(iHz).z(lmaxval_k(layoutnumber + 1));
                    }
                    if (lminval(layoutnumber + 1) > valor) {
                        lminval(layoutnumber + 1) = valor;
                        lminval_i(layoutnumber + 1) = i + b.Hx.XI;
                        lminval_j(layoutnumber + 1) = j + b.Hy.YI;
                        lminval_k(layoutnumber + 1) = k + b.Hz.ZI;
                        lminval_x(layoutnumber + 1) = Punto.PhysCoor(iHx).x(lminval_i(layoutnumber + 1));
                        lminval_y(layoutnumber + 1) = Punto.PhysCoor(iHy).y(lminval_j(layoutnumber + 1));
                        lminval_z(layoutnumber + 1) = Punto.PhysCoor(iHz).z(lminval_k(layoutnumber + 1));
                    }
                }
            }
        }

        //
        NEWlmaxval.assign(num_procs, 0.0);
        NEWlmaxval_i.assign(num_procs, 0);
        NEWlmaxval_j.assign(num_procs, 0);
        NEWlmaxval_k.assign(num_procs, 0);
#ifdef CompileWithMPI
        MPI_AllReduce(LMAXVAL.data(), NEWlmaxval.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LMAXVAL_i.data(), NEWlmaxval_I.data(), num_procs, MPI_INT, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LMAXVAL_j.data(), NEWlmaxval_J.data(), num_procs, MPI_INT, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LMAXVAL_k.data(), NEWlmaxval_K.data(), num_procs, MPI_INT, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LMAXVAL_x.data(), NEWlmaxval_x.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LMAXVAL_y.data(), NEWlmaxval_y.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LMAXVAL_z.data(), NEWlmaxval_z.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
#else
        NEWlmaxval = LMAXVAL;
        NEWlmaxval_i = LMAXVAL_I;
        NEWlmaxval_j = LMAXVAL_J;
        NEWlmaxval_k = LMAXVAL_K;
        NEWlmaxval_x = LMAXVAL_x;
        NEWlmaxval_y = LMAXVAL_y;
        NEWlmaxval_z = LMAXVAL_z;
#endif
        qmaxval = 0.0;
        qmaxval_i = 0;
        qmaxval_j = 0;
        qmaxval_k = 0;
        qmaxval_x = 0.0;
        qmaxval_y = 0.0;
        qmaxval_z = 0.0;
        for (i = 1; i <= num_procs; ++i) {
            if (std::abs(NEWlmaxval(i - 1)) > qmaxval) {
                qmaxval = std::abs(NEWlmaxval(i - 1));
                qmaxval_i = NEWlmaxval_i(i - 1);
                qmaxval_j = NEWlmaxval_j(i - 1);
                qmaxval_k = NEWlmaxval_k(i - 1);
                qmaxval_x = NEWlmaxval_x(i - 1);
                qmaxval_y = NEWlmaxval_y(i - 1);
                qmaxval_z = NEWlmaxval_z(i - 1);
            }
        }

#ifdef CompileWithMPI
        MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif
        //
        NEWlminval.assign(num_procs, 0.0);
        NEWlminval_i.assign(num_procs, 0);
        NEWlminval_j.assign(num_procs, 0);
        NEWlminval_k.assign(num_procs, 0);
#ifdef CompileWithMPI
        MPI_AllReduce(LminVAL.data(), NEWlminval.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LminVAL_i.data(), NEWlminval_I.data(), num_procs, MPI_INT, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LminVAL_j.data(), NEWlminval_J.data(), num_procs, MPI_INT, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LminVAL_k.data(), NEWlminval_K.data(), num_procs, MPI_INT, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LminVAL_x.data(), NEWlminval_x.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LminVAL_y.data(), NEWlminval_y.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
        MPI_AllReduce(LminVAL_z.data(), NEWlminval_z.data(), num_procs, MPI_DOUBLE, MPI_SUM, SUBCOMM_MPI, &ierr);
#else
        NEWlminval = LminVAL;
        NEWlminval_i = LminVAL_I;

NEWlminval_j = LminVAL_J;
            NEWlminval_k = LminVAL_K;
            NEWlminval_x = LminVAL_x;
            NEWlminval_y = LminVAL_y;
            NEWlminval_z = LminVAL_z;
#endif
            qminval = 0.0;
            qminval_i = 0;
            qminval_j = 0;
            qminval_k = 0;
            qminval_x = 0.0;
            qminval_y = 0.0;
            qminval_z = 0.0;
            for (i = 1; i <= num_procs; ++i) {
                if (std::abs(NEWlminval(i)) > qminval) {
                    qminval = std::abs(NEWlminval(i));
                    qminval_i = newlminval_i(i);
                    qminval_j = newlminval_j(i);
                    qminval_k = newlminval_k(i);
                    qminval_x = newlminval_x(i);
                    qminval_y = newlminval_y(i);
                    qminval_z = newlminval_z(i);
                }
            }

#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

            // !!!!!!!!!!!!!!!!!!!!!!!
            // escritura del fichero snap a voluntad o cuando se pase un umbral, cada minuto
            if (layoutnumber == 0) {
                inquire_file_exists("snap", mustSnap);
                if (mustSnap) {
                    // Clear the flushing signaling file
                    std::ifstream file35("snap");
                    if (file35) {
                        file35 >> snapLevel;
                        file35 >> snapStep;
                        file35 >> snapHowMany;
                        file35.close();
                        erasesignalingfiles(simu_devia);
                    } else {
                        // Handle error or end condition similar to Fortran's end=1153
                        // Assuming erasesignalingfiles is called regardless if file open fails or not based on logic flow
                        // In Fortran, if read fails, it goes to 1153 which closes and calls erasesignalingfiles.
                        // Here we just close (already closed if opened) and call.
                        erasesignalingfiles(simu_devia);
                    }
                }
            }
#ifdef CompileWithMPI
            MPI_Barrier(MPI_COMM_WORLD, &ierr);
            MPI_Bcast(&mustSnap, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, &ierr);
            MPI_Bcast(&snapLevel, 1, REALSIZE, 0, MPI_COMM_WORLD, &ierr);
            MPI_Bcast(&snapStep, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, &ierr);
            MPI_Bcast(&snapHowMany, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, &ierr);
            MPI_Barrier(MPI_COMM_WORLD, &ierr);
#endif

            if ((mustSnap && (lmaxval(layoutnumber + 1) > snapLevel)) || (countersnap > 0)) {
                countersnap = countersnap + 1;
                //
                dimxsnap = static_cast<int>((fin_ibox - ini_ibox) / snapstep) + 1;
                dimysnap = static_cast<int>((fin_jbox - ini_jbox) / snapstep) + 1;
                dimzsnap = static_cast<int>((fin_kbox - ini_kbox) / snapstep) + 1;
                if (!snap_allocated) {
                    snap.resize(dimxsnap + 1, dimysnap + 1, dimzsnap + 1, 1); // Adjusting indices to match Fortran 1-based if needed, but vector is 0-based.
                    // Fortran: snap(ini_ibox:ini_ibox+dimxsnap, ...)
                    // This implies size is dimxsnap+1.
                    // To map Fortran index `ini_ibox` to C++ index `0`, we might need an offset or just resize to cover the range.
                    // Let's assume snap is accessed with absolute indices relative to ini_*.
                    // If we use a vector, we need to ensure it's large enough.
                    // Let's resize to max possible index.
                    // However, simpler translation: allocate vector of size [dimxsnap+1][dimysnap+1][dimzsnap+1][1]
                    // And access snap[i - ini_ibox + offset]... 
                    // But the code uses: snap(ini_ibox + int((i-ini_ibox)/snapstep), ...)
                    // Let's just make the vector large enough to hold the max index.
                    // Max index for x: ini_ibox + dimxsnap.
                    // So size should be at least ini_ibox + dimxsnap + 1.
                    // This is tricky with std::vector without a wrapper. 
                    // Given the instruction "Convert arrays to std::vector", usually this means replacing the array declaration.
                    // Access patterns like `snap(idx)` where `idx` is calculated need to be handled.
                    // If `snap` was declared as `snap(ini_ibox:ini_ibox+dimxsnap, ...)`, the size is `dimxsnap+1`.
                    // The index `ini_ibox + int(...)` starts at `ini_ibox`.
                    // So we need a vector that can be indexed from `ini_ibox`.
                    // Standard C++ vectors are 0-indexed.
                    // We will assume `snap` is a global or member variable that is resized.
                    // To preserve names, we keep `snap`.
                    // We will resize it to accommodate the highest index.
                    // Highest index X: ini_ibox + dimxsnap
                    // Highest index Y: ini_jbox + dimysnap
                    // Highest index Z: ini_kbox + dimzsnap
                    // We'll resize to [ini_ibox + dimxsnap + 1][ini_jbox + dimysnap + 1][ini_kbox + dimzsnap + 1][1]
                    // This is inefficient but preserves the indexing logic directly.
                    // Alternatively, we could use a flattened vector and calculate index, but that changes logic.
                    // Let's stick to the simplest direct translation of the allocation.
                    
                    // Note: The Fortran code allocates `snap(ini_ibox:ini_ibox+dimxsnap, ...)`
                    // This means the first dimension has size `dimxsnap + 1`.
                    // The index `ini_ibox` corresponds to the first element.
                    // In C++, if we use `std::vector<std::vector<std::vector<std::vector<double>>>> snap`,
                    // we can't easily have 1-based indexing or arbitrary lower bounds without a wrapper.
                    // However, the prompt asks to convert arrays to std::vector.
                    // I will resize the vector to be large enough to hold the indices used.
                    
                    size_t max_x = ini_ibox + dimxsnap;
                    size_t max_y = ini_jbox + dimysnap;
                    size_t max_z = ini_kbox + dimzsnap;
                    
                    snap.resize(max_x + 1);
                    for (size_t x = 0; x <= max_x; ++x) {
                        snap[x].resize(max_y + 1);
                        for (size_t y = 0; y <= max_y; ++y) {
                            snap[x][y].resize(max_z + 1);
                            for (size_t z = 0; z <= max_z; ++z) {
                                snap[x][y][z].resize(1);
                            }
                        }
                    }
                    snap_allocated = true;
                }
                
                // Initialize to 0.0
                for (size_t x = 0; x < snap.size(); ++x) {
                    for (size_t y = 0; y < snap[x].size(); ++y) {
                        for (size_t z = 0; z < snap[x][y].size(); ++z) {
                            snap[x][y][z][0] = 0.0;
                        }
                    }
                }
                
                // !!!!             k = ini_kbox - snapStep
                // !!!!             do while (k < fin_kbox )
                // !!!!                k = min (k  + snapStep , fin_kbox)
                // !!!!                j = ini_jbox - snapStep
                // !!!!                do while (j < fin_jbox )
                // !!!!                    j = min(j  + snapStep , fin_jbox)
                // !!!!                    i = ini_ibox - snapStep
                // !!!!                    do while (i < fin_ibox )
                // !!!!                        i = min (i  + snapStep , fin_ibox)
                // !!!!                        veces=0
                // !!!!                        valor=0.0_RKIND
                // !!!!                        do k1=0,snapstep-1
                // !!!!                        do j1=0,snapstep-1
                // !!!!                        do i1=0,snapstep-1
                // !!!!                        if ((i+i1 <= fin_ibox).and.(j+j1 <= fin_jbox).and.(k+k1 <= fin_kbox)) then
                // !!!!                            valor = valor+sqrt(Ex(i+i1, j+j1, k+k1) * Ex( i+i1, j+j1, k+k1) + &
                // !!!!                                               Ey( i+i1, j+j1, k+k1) * Ey(i+i1, j+j1, k+k1)+ &
                // !!!!                                            Ez(i+i1, j+j1, k+k1) * Ez( i+i1, j+j1, k+k1))
                // !!!!                            veces=veces+1
                // !!!!                        end if
                // !!!!                        end do
                // !!!!                        end do
                // !!!!                        end do
                // !!!!                        snap(ini_ibox+int((i-ini_ibox)/snapstep),ini_jbox+int((j-ini_jbox)/snapstep), &
                // !!!!                            ini_kbox+int((k-ini_kbox)/snapstep),1) = valor/veces
                // !!!!                    end do
                // !!!!                end do
                // !!!!             end do

                for (k = ini_kbox; k <= fin_kbox; k += snapStep) {
                    for (j = ini_jbox; j <= fin_jbox; j += snapStep) {
                        for (i = ini_ibox; i <= fin_ibox; i += snapStep) {
                            veces = 0;
                            valor = 0.0;
                            for (k1 = 0; k1 < snapstep; ++k1) {
                                for (j1 = 0; j1 < snapstep; ++j1) {
                                    for (i1 = 0; i1 < snapstep; ++i1) {
                                        if ((i + i1 <= fin_ibox) && (j + j1 <= fin_jbox) && (k + k1 <= fin_kbox)) {
                                            valor += std::sqrt(Ex(i + i1, j + j1, k + k1) * Ex(i + i1, j + j1, k + k1) +
                                                               Ey(i + i1, j + j1, k + k1) * Ey(i + i1, j + j1, k + k1) +
                                                               Ez(i + i1, j + j1, k + k1) * Ez(i + i1, j + j1, k + k1));
                                            veces = veces + 1;
                                        }
                                    }
                                }
                            }
                            // Calculate indices for snap array
                            // Fortran: snap(ini_ibox+int((i-ini_ibox)/snapstep), ...)
                            // In C++, if snap is resized to cover these indices, we access directly.
                            // Note: int() in Fortran truncates towards zero. For positive numbers, it's floor.
                            size_t idx_x = ini_ibox + static_cast<size_t>((i - ini_ibox) / snapstep);
                            size_t idx_y = ini_jbox + static_cast<size_t>((j - ini_jbox) / snapstep);
                            size_t idx_z = ini_kbox + static_cast<size_t>((k - ini_kbox) / snapstep);
                            
                            snap[idx_x][idx_y][idx_z][0] = valor / veces;
                        }
                    }
                }

                // Write instant number
                chinstant = std::to_string(n);
                // Format minmax: '_', lminval, '_', lmaxval, '_'
                // Using a simple string stream or manual formatting. 
                // Fortran: '(a,e15.4e3,a,e15.4e3,a)'
                // This is scientific notation.
                std::ostringstream minmax_stream;
                minmax_stream << "_" << std::scientific << std::setprecision(4) << lminval(layoutnumber + 1) << "_" 
                              << std::scientific << std::setprecision(4) << lmaxval(layoutnumber + 1) << "_";
                minmax = minmax_stream.str();
                
                fichsnap = nEntradaRoot + "_snap_" + chinstant + "_" + whoamishort;

#ifdef CompileWithHDF
                ficherito = fichsnap + ".h5";
                openclosedelete(ficherito);
                
                write_xdmfsnap(n, fichsnap, 
                               ini_ibox + b->Ex->XI, ini_ibox + dimxsnap + b->Ex->XI,
                               ini_jbox + b->Ex->YI, ini_jbox + dimysnap + b->Ex->YI,
                               ini_kbox + b->Ex->ZI, ini_kbox + dimzsnap + b->Ex->ZI, 
                               snap);
#endif
                //             open (35,file=trim(adjustl(fichsnap))//'.bin')
                //             write (35,*) '!END'
                //             close (35,status='delete')
                //             open (35,file=trim(adjustl(fichsnap))//'.bin',form='unformatted',status='new',action='write')
                //             write (35) n,lminval  (layoutnumber+1),lmaxval  (layoutnumber+1)
                //             write (35) ini_ibox,fin_ibox,ini_jbox,fin_jbox,ini_kbox,fin_kbox
                //             write (35) ini_iboxsin,fin_iboxsin,ini_jboxsin,fin_jboxsin,ini_kboxsin,fin_kboxsin
                //             write (35) (Punto%PhysCoor(iHx)%x(i+b%Hx%XI),i = ini_iboxsin, fin_iboxsin)
                //             write (35) (Punto%PhysCoor(iHy)%y(j+b%Hy%YI),j = ini_jboxsin, fin_jboxsin)
                //             write (35) (Punto%PhysCoor(iHz)%z(k+b%Hz%ZI),k = ini_kboxsin, fin_kboxsin)
                // #ifndef CompileWithHDF
                //             do k = ini_kbox, fin_kbox
                //                do j = ini_jbox, fin_jbox
                //                       write (35) (snap(i,j,k,1),i = ini_ibox, fin_ibox)
                //                end do
                //             end do
                // #endif
                //             close (35)

                dubuf = whoami + " Written Snap file at n= " + std::to_string(n) + " max field over " + 
                        std::to_string(maxval(snap)) + " > " + std::to_string(snapLevel);
                
                // Deallocate snap
                snap.clear();
                snap.shrink_to_fit();
                snap_allocated = false;
                
                print11(layoutnumber, dubuf, true);
                if (countersnap >= snapHowMany) {
                    mustSnap = false;
                    countersnap = 0;
                }
            }

#ifdef CompileWithMPI
            MPI_Barrier(MPI_COMM_WORLD, &ierr); // TODOS STOCH O NO 060619
#endif
            // !!!!!!!!!!!!!!!!!!!!!!!

            //
            //
            if (layoutnumber == 0) {
                //
                dubuf = SEPARADOR + nentradaroot + SEPARADOR;
                print11(layoutnumber, dubuf);
                dubuf = "Switches: " + opcionestotales;
            }

print11(layoutnumber, dubuf);
            // if (num_procs != 1) {
                std::ostringstream oss_temp;
                oss_temp << "MPI Processes: " << num_procs;
                dubuf = oss_temp.str();
                print11(layoutnumber, dubuf);
            // }
            //
            {
                std::ostringstream oss_temp;
                oss_temp << "Date/Time " 
                         << time_out2.fecha.substr(6, 2) << "/" 
                         << time_out2.fecha.substr(4, 2) << "/" 
                         << time_out2.fecha.substr(0, 4) << "   " 
                         << time_out2.hora.substr(0, 2) << ":" 
                         << time_out2.hora.substr(2, 2) << ":" 
                         << time_out2.hora.substr(4, 2);
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            {
                std::ostringstream oss_temp;
                oss_temp << "Simulated:" << n << "/" << finaltimestep << " steps";
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            if (permitscaling) {
                std::ostringstream oss_temp;
                oss_temp << std::scientific << std::setprecision(9) 
                         << "Time= " << sgg.tiempo[n] 
                         << ", dt0 (original)= " << dt0 
                         << ", dt(pscaled)= " << sgg.dt;
                dubuf = oss_temp.str();
            } else {
                std::ostringstream oss_temp;
                oss_temp << std::scientific << std::setprecision(9) 
                         << "Time= " << sgg.tiempo[n] 
                         << ", dt0 = " << sgg.dt;
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            if (energytotal > oldenergytotal) {
                std::ostringstream oss_temp;
                oss_temp << "Total Energy (inc) : " << energytotal;
                dubuf = oss_temp.str();
                if (simu_devia) {
                    dubuf = dubuf.substr(0, dubuf.find_last_not_of(" \t\n\r\f\v") + 1) + " (Stoch)";
                }
                oldenergytotal = energytotal;
            } else {
                std::ostringstream oss_temp;
                oss_temp << "Total Energy (dec) : " << energytotal;
                dubuf = oss_temp.str();
                if (simu_devia) {
                    dubuf = dubuf.substr(0, dubuf.find_last_not_of(" \t\n\r\f\v") + 1) + " (Stoch)";
                }
                oldenergytotal = energytotal;
            }
            print11(layoutnumber, dubuf);
            //
            if (qmaxval_x < -1e19) {
                dondex = " PML ";
            } else {
                std::ostringstream oss_temp;
                oss_temp << std::scientific << std::setprecision(4) << qmaxval_x;
                dondex = oss_temp.str();
            }
            if (qmaxval_y < -1e19) {
                dondey = " PML ";
            } else {
                std::ostringstream oss_temp;
                oss_temp << std::scientific << std::setprecision(4) << qmaxval_y;
                dondey = oss_temp.str();
            }
            if (qmaxval_z < -1e19) {
                dondez = " PML ";
            } else {
                std::ostringstream oss_temp;
                oss_temp << std::scientific << std::setprecision(4) << qmaxval_z;
                dondez = oss_temp.str();
            }

            {
                std::ostringstream oss_temp;
                oss_temp << std::scientific << std::setprecision(4) 
                         << "Max field: " << qmaxval << " at (" 
                         << qmaxval_i << "," << qmaxval_j << "," << qmaxval_k << ")=( " 
                         << trim_adjustl(dondex) << ", " 
                         << trim_adjustl(dondey) << ", " 
                         << trim_adjustl(dondez) << ")";
                dubuf = oss_temp.str();
            }

            if (simu_devia) {
                dubuf = dubuf.substr(0, dubuf.find_last_not_of(" \t\n\r\f\v") + 1) + " (Stoch)";
            }
            print11(layoutnumber, dubuf);
            for (int i = 1; i <= num_procs; ++i) {
                if (newlmaxval_x[i] < -1e19) {
                    dondex = " PML ";
                } else {
                    std::ostringstream oss_temp;
                    oss_temp << std::scientific << std::setprecision(4) << newlmaxval_x[i];
                    dondex = oss_temp.str();
                }
                if (newlmaxval_y[i] < -1e19) {
                    dondey = " PML ";
                } else {
                    std::ostringstream oss_temp;
                    oss_temp << std::scientific << std::setprecision(4) << newlmaxval_y[i];
                    dondey = oss_temp.str();
                }
                if (newlmaxval_z[i] < -1e19) {
                    dondez = " PML ";
                } else {
                    std::ostringstream oss_temp;
                    oss_temp << std::scientific << std::setprecision(4) << newlmaxval_z[i];
                    dondez = oss_temp.str();
                }

                {
                    std::ostringstream oss_temp;
                    oss_temp << std::scientific << std::setprecision(4) 
                             << "Max field slice: " << i << " " << NEWlmaxval[i] << "/" << maxSourceValue 
                             << " at (" << newlmaxval_i[i] << "," << newlmaxval_j[i] << "," << newlmaxval_k[i] << ")=( " 
                             << trim_adjustl(dondex) << ", " 
                             << trim_adjustl(dondey) << ", " 
                             << trim_adjustl(dondez) << ")";
                    dubuf = oss_temp.str();
                }
                // call print11(layoutnumber,dubuf) !comentado para que la salida sea menos verbose
            }
            //

            //
            {
                std::ostringstream oss_temp;
                oss_temp << "Mins. since start  : " << std::ceil((time_end - time_desdelanzamiento) / 60.0);
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            {
                std::ostringstream oss_temp;
                double mins_until = std::ceil(((finaltimestep - n) * megaceldastotales) / speedGlobAvg / 60.0);
                double mins_past = std::ceil((time_end - time_desdelanzamiento) / 60.0);
                oss_temp << "Mins. until end    : " << std::min(mins_until, maxCPUtime - mins_past);
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            if (everflushed) {
                std::ostringstream oss_temp;
                oss_temp << "Mins. past flushing: " << std::ceil((time_end - time_begin2) / 60.0);
                dubuf = oss_temp.str();
                print11(layoutnumber, dubuf);
            } else {
                dubuf = "Never flushed resuming fields.";
                print11(layoutnumber, dubuf);
            }
            if (flushsecondsFields != 0) {
                std::ostringstream oss_temp;
                double mins_next_flush = std::ceil((flushsecondsFields - (time_end - time_begin2)) / 60.0);
                double mins_remaining = std::ceil(((finaltimestep - n) * megaceldastotales) / speedGlobAvg / 60.0);
                oss_temp << "Mins. next flushing: " << std::min(mins_next_flush, mins_remaining);
                dubuf = oss_temp.str();
                print11(layoutnumber, dubuf);
            } else {
                if (maxCPUtime == topCPUtime) {
                    dubuf = "Will Never flush resuming fields.";
                    print11(layoutnumber, dubuf);
                } else {
                    dubuf = "Flushing of restarting DATA at the end.";
                    print11(layoutnumber, dubuf);
                }
            }
            //
            {
                std::ostringstream oss_temp;
                oss_temp << "Next info at step: " << n_info;
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            {
                std::ostringstream oss_temp;
                oss_temp << "Total Mcells:" << megaceldastotales;
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            if (reportedinstant < n) {
                std::ostringstream oss_temp;
                oss_temp << "Mcells/sec  : " << speedGlobInst << " (" << reportedinstant << " to " << N << ")";
                dubuf = oss_temp.str();
                if (simu_devia) {
                    dubuf = dubuf.substr(0, dubuf.find_last_not_of(" \t\n\r\f\v") + 1) + " (Stoch)";
                }
                print11(layoutnumber, dubuf);
            }
            //
            {
                std::ostringstream oss_temp;
                oss_temp << "Mcells/sec  : " << speedGlobAvg << " (" << INITIALtimeSTEP << " to " << N << ")";
                dubuf = oss_temp.str();
            }
            if (simu_devia) {
                dubuf = dubuf.substr(0, dubuf.find_last_not_of(" \t\n\r\f\v") + 1) + " (Stoch)";
            }
            print11(layoutnumber, dubuf);
            //
            {
                std::ostringstream oss_temp;
                oss_temp << SEPARADOR << separador << separador;
                dubuf = oss_temp.str();
            }
            print11(layoutnumber, dubuf);
            //
            output_stream_10 << sgg.tiempo[n] << " " << energytotal << std::endl;
            //write(67,'(i5)') nint(100.0_RKIND * n/finaltimestep) !percentage
            std::flush(std::cout); // Assuming flush(11) maps to cout or similar, usually 11 is stdout in these legacy codes
            std::flush(output_stream_10);
         }
         //
#ifdef CompileWithMPI
         MPI_Barrier(MPI_COMM_WORLD, &ierr);
#endif
         get_secnds(time_out2);
         time_begin = time_out2.segundos; //restart timing
         reportedinstant = n + 1;
      } // every reporting seconds
      //

      //stop if this probe blows up
      stoponNaN_aux = false;
      // #ifdef ARCHITECTURE_SUN
      // #ifdef CompileWithReal4
      //       if (false) {
      // #else
      //       if (false) {
      // #endif
      // #else
      // #ifdef PreventCrayBug
      // #ifdef CompileWithReal4
      //       if (IsNaNf(energy)) {
      // #else
      //       if (IsNaNd(energy)) {
      // #endif
      // #else
       //  if (IsNaN (energy)) then !quitado a mano para que PGI no se queje a 150623 !fm
      // #endif
      // #endif
         //
       //    stoponNaN_aux=.true.
       //  end if
#ifdef CompileWithMPI
         MPI_AllReduce(&stoponNaN_aux, &stoponNaN, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
#else
         stoponNaN = stoponNaN_aux;
#endif
         if (stoponNaN) {
            if (layoutnumber == 0) {
                std::ostringstream oss_temp;
                oss_temp << "ERROR, ABORTING: UNSTABILITIES. Possible reasons and fixes: ";
                dubuf = oss_temp.str();
                print11(layoutnumber, dubuf);
                dubuf = "     In case of single wires: reduce WIRE radii and/or reduce sgg%dt";
                print11(layoutnumber, dubuf);
                dubuf = "     In case of surface IBCs: reduce -att factor";
                print11(layoutnumber, dubuf);
            }
#ifdef CompileWithMPI
            MPI_Barrier(MPI_COMM_WORLD, &ierr);
#endif
            parar = true;
            //            call StopOnError(layoutnumber,num_procs,' Aborting')
         }
         //
         l_aux = ((time_end - time_begin2) > flushsecondsFields && flushsecondsFields != 0) || mustflushFIELDS;
#ifdef CompileWithMPI
         //print *,'layoutnumber+1,l_aux, hay_flushFIELDSl pre',layoutnumber+1,l_aux, hay_flushFIELDS
         MPI_AllReduce(&l_aux, &hay_flushFIELDS, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
         //print *,'layoutnumber+1,l_aux, hay_flushFIELDSl post',layoutnumber+1,l_aux, hay_flushFIELDS
#else
         hay_flushFIELDS = l_aux;
#endif

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <mpi.h>

// Assuming global variables and helper functions are defined elsewhere or in a namespace.
// To preserve names as requested, we assume these exist in the global scope or a specific namespace.
// Since the prompt asks to convert modules to namespaces, but no module structure is provided for this snippet,
// we will place these functions in a namespace 'FortranModule' to simulate the module context,
// or just global functions if no module is implied. Given the "Preserve all names" instruction,
// we will keep variable names like 'hay_flushDATA', 'mustflushDATA', etc.

// Forward declarations / Global state simulation based on Fortran implicit typing and usage
extern bool hay_flushDATA;
extern bool mustflushDATA;
extern bool mustflushFIELDS;
extern bool hay_flushFIELDS;
extern bool mustunpack;
extern bool mustpostprocess;
extern bool mustflushXdmf;
extern bool mustflushVTK;
extern int layoutnumber;
extern int num_procs;
extern std::string nEntradaRoot;
extern bool verbosete;
extern bool ignoreErrors1;
extern bool warningfileIsOpen;
extern std::string warningfile;
extern bool fatalerror;
extern int CONTADORDEMENSAJES;
extern bool stoch_undivided;
extern bool simu_devia;
extern std::string WarningFile;

// Helper structs/classes if needed
struct TimeOut {
    double segundos;
};

// Helper functions assumed to exist
void get_secnds(TimeOut& time_out2);
void erasesignalingfiles(bool simu_devia);
void openclosedelete(const std::string& filename);
void opensolo(int unit, const std::string& filename);
void closesolo(int unit);
void trimnullchar(std::string& str);

// MPI variables
extern MPI_Comm SUBCOMM_MPI;
extern int ierr;

namespace FortranModule {

void Timing() {
    bool l_aux;
    double time_end, time_begin3, flushsecondsDATA;
    // Assuming these variables are available in the scope or passed in.
    // In a real translation, these would likely be members of a class or globals.
    // We assume they are accessible.
    
    // Note: time_end, time_begin3, flushsecondsDATA are not declared in the snippet.
    // We assume they are global or member variables.
    
    l_aux = ((time_end - time_begin3) > flushsecondsDATA) && (flushsecondsDATA != 0);
    l_aux = l_aux || mustflushDATA;

#ifdef CompileWithMPI
    // print *,'layoutnumber+1,l_aux, hay_flushDATA pre',layoutnumber+1,l_aux, hay_flushDATA
    MPI_AllReduce(&l_aux, &hay_flushDATA, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, &ierr);
    // print *,'layoutnumber+1,l_aux, hay_flushDATA post',layoutnumber+1,l_aux, hay_flushDATA
#else
    hay_flushDATA = l_aux;
#endif

    if (hay_flushFIELDS) {
        mustflushFIELDS = false;
#ifdef CompileWithMPI
        MPI_Barrier(MPI_COMM_WORLD, &ierr);
#endif
        TimeOut time_out2;
        get_secnds(time_out2);
        double time_begin2 = time_out2.segundos;
        // perform->flushFIELDS = true; // Assuming 'perform' is a global pointer or object
        // Since 'perform' is not defined, we assume it's a global object.
        // perform.flushFIELDS = true; 
        // Clear the flushing signaling file
        if (layoutnumber == 0) { // only the master proc mush erase this
            erasesignalingfiles(simu_devia);
        }
    }

    if (hay_flushDATA) {
        mustflushDATA = false;
#ifdef CompileWithMPI
        MPI_Barrier(MPI_COMM_WORLD, &ierr);
#endif
        TimeOut time_out2;
        get_secnds(time_out2);
        double time_begin3 = time_out2.segundos;
        // perform->flushDATA = true;
        // Clear the flushing signaling file
        if (layoutnumber == 0) { // only the master proc mush erase this
            erasesignalingfiles(simu_devia);
        }
    }

    if (mustunpack) {
        mustunpack = false;
        // perform->unpack = true;
        // Clear the flushing signaling file
        if (layoutnumber == 0) { // only the master proc mush erase this
            erasesignalingfiles(simu_devia);
        }
    }

    if (mustpostprocess) {
        mustpostprocess = false;
        // perform->postprocess = true;
        // Clear the flushing signaling file
        if (layoutnumber == 0) { // only the master proc mush erase this
            erasesignalingfiles(simu_devia);
        }
    }

    if (mustflushXdmf) {
        mustflushXdmf = false;
        // perform->flushXdmf = true;
        // Clear the flushing signaling file
        if (layoutnumber == 0) { // only the master proc mush erase this
            erasesignalingfiles(simu_devia);
        }
    }

    if (mustflushVTK) {
        mustflushVTK = false;
        // perform->flushVTK = true;
        if (layoutnumber == 0) { // only the master proc mush erase this
            erasesignalingfiles(simu_devia);
        }
    }
    // ---------------------------> acaba Timing <----------------------------------------------------
    return;
}

void INITWARNINGFILE(int layoutnumber, int num_procs, const std::string& nEntradaRoot, bool verbosete, bool ignoreErrors1) {
    // character(len=*) :: nEntradaRoot
    // integer(kind=4), intent(in) :: layoutnumber,num_procs
    // file management
    // character(len=BUFSIZE) :: whoamishort
    // #ifdef CompileWithMPI
    // integer(kind=MPI_OFFSET_KIND) disp
    // integer(kind=4) :: ierr
    // #endif
    // logical verbosete,ignoreerrors1       , itsopen2
    // integer :: my_iostat
    // character(len=BUFSIZE) :: ficherito
    // verbose=verbosete
    // ignoreerrors=ignoreerrors1
    // write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',num_procs,') '
    // write(whoamishort,'(i5)') layoutnumber+1
    // if (layoutnumber == 0) then
    // !!!inquire(unit=17, opened=itsopen2)
    // !!!if (itsopen2) print *,'----------->17 open!!!'
    //   ficherito=trim(adjustl(nEntradaRoot))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
    //   call openclosedelete(ficherito)
    // end if
    // !!!#ifdef CompileWithMPI
    // !!!if (SIZE/=0) then
    // !!!    call MPI_Barrier(SUBCOMM_MPI,ierr)
    // !!!    call MPI_FILE_open (SUBCOMM_MPI, trim(adjustl(nEntradaRoot))//'_tmpWarnings.txt', &
    // !!!                           MPI_MODE_WRONLY + MPI_MODE_CREATE, &
    // !!!                           MPI_INFO_NULL, thefile, ierr)
    // !!!    disp = (layoutnumber+1) * BUFSIZE * maxmessages !no creo que se den mas de 2000 mensajes por layout
    // !!!
    // !!!    call MPI_FILE_SET_VIEW(thefile, disp, MPI_CHARACTER, &
    // !!!                               MPI_CHARACTER, 'native', &
    // !!!                               MPI_INFO_NULL, ierr)
    // !!!ELSE
    // !!!    open (17,file=trim(adjustl(nEntradaRoot))//'_tmpWarnings.txt',form='formatted')
    // !!!end if
    // !!!#endif
    // inquire(unit=17, opened=itsopen2)
    // !!!if (itsopen2) print *,'----------->17 open!!!'
    // ficherito=trim(adjustl(nEntradaRoot))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
    // call opensolo(17,ficherito)
    // !!!#endif
    // warningfileIsOpen=.true.
    // warningfile=nEntradaRoot
    // fatalerror = .false.
    // CONTADORDEMENSAJES=0

    bool itsopen2;
    int my_iostat;
    std::string whoamishort;
    std::string ficherito;
    
    // verbose = verbosete; // Assuming global 'verbose'
    // ignoreerrors = ignoreErrors1; // Assuming global 'ignoreerrors'

    // write(whoami,'(a,i5,a,i5,a)') '(',layoutnumber+1,'/',num_procs,') '
    // Assuming 'whoami' is a global string
    // std::ostringstream oss;
    // oss << "(" << layoutnumber + 1 << "/" << num_procs << ") ";
    // whoami = oss.str();

    // write(whoamishort,'(i5)') layoutnumber+1
    char buf[10];
    snprintf(buf, sizeof(buf), "%5d", layoutnumber + 1);
    whoamishort = buf;

    if (layoutnumber == 0) {
        // !!!inquire(unit=17, opened=itsopen2)
        // !!!if (itsopen2) print *,'----------->17 open!!!'
        ficherito = nEntradaRoot + whoamishort + "_tmpWarnings.txt";
        openclosedelete(ficherito);
    }

    // inquire(unit=17, opened=itsopen2)
    // !!!if (itsopen2) print *,'----------->17 open!!!'
    ficherito = nEntradaRoot + whoamishort + "_tmpWarnings.txt";
    opensolo(17, ficherito);
    // !!!#endif

    warningfileIsOpen = true;
    warningfile = nEntradaRoot;
    fatalerror = false;
    CONTADORDEMENSAJES = 0;
}

void WarnErrReport(const std::string& bufff, bool error = false) {
    // use iso_fortran_env, only : error_unit
    // logical :: itsopen
    // logical, optional :: error
    // #ifdef CompileWithMPI
    // integer(kind=4) :: ierr
    // #endif
    // character(len=*), intent(in) :: bufff
    // character(len=BUFSIZE) :: buff2,buff3
    // if (present(error)) then
    //    fatalerror = error .or. fatalerror
    // end if
    // !
    // buff3=trim(adjustl(whoami))//' '//trim(adjustl(bufff))
    // call trimnullchar(buff3)
    // buff2=CHAR(13)//CHAR(10)//trim(adjustl(buff3))
    // call trimnullchar(buff2)
    // inquire(unit=17, opened=itsopen)
    // if (itsopen) write (17,'(a)',err=154) trim(adjustl(buff3))
    // write (error_unit,'(a)') trim(adjustl(buff3))
    // goto 155
    // 154   inquire(unit=17, opened=itsopen) 
    // print *,itsopen,'- Cannot write into warning file the message: ',trim(adjustl(buff3))
    // 155   return

    bool itsopen;
    std::string buff2;
    std::string buff3;

    if (error) {
        fatalerror = error || fatalerror;
    }

    // Assuming 'whoami' is a global string
    buff3 = whoami + " " + bufff;
    trimnullchar(buff3);

    buff2 = "\r\n" + buff3;
    trimnullchar(buff2);

    // inquire(unit=17, opened=itsopen)
    // We need to check if unit 17 is open. Assuming a helper function or global state.
    // For this translation, we assume 'itsopen' is determined by checking if the file stream is valid.
    // Since we don't have the exact implementation of 'inquire', we'll simulate it.
    // In a real C++ translation, we'd track open files.
    // Let's assume 'itsopen' is true if unit 17 was opened in INITWARNINGFILE.
    // We'll use a global flag or check.
    // For simplicity, let's assume we check if the file associated with unit 17 is open.
    // Since we don't have the file object, we'll assume 'itsopen' is a global variable or function.
    // Let's assume 'itsopen' is true if warningfileIsOpen is true and unit 17 is the one used.
    // Actually, the Fortran code checks unit 17 specifically.
    // We'll assume a global 'bool unit17_open' or similar.
    // To be safe, we'll just check if 'warningfileIsOpen' is true, as unit 17 is used there.
    // But the code checks 'inquire(unit=17, opened=itsopen)' again.
    // Let's assume 'itsopen' is true if the file was opened.
    // We'll use a global variable 'bool unit17_is_open' set in INITWARNINGFILE.
    
    // For the sake of this translation, let's assume 'itsopen' is true if warningfileIsOpen is true.
    // This is a simplification.
    itsopen = warningfileIsOpen; 

    if (itsopen) {
        // write (17,'(a)',err=154) trim(adjustl(buff3))
        // We need to write to unit 17. Assuming a global ofstream or similar.
        // Let's assume 'ofstream unit17_stream' is global.
        // unit17_stream << buff3 << std::endl;
        // If write fails, goto 154.
        // Since we don't have the stream, we'll just print to error_unit.
        // In a real implementation, we'd handle the error.
        // For now, we'll just write to cerr.
        std::cerr << buff3 << std::endl;
    } else {
        goto label_154;
    }
    
    goto label_155;

label_154:
    // inquire(unit=17, opened=itsopen)
    itsopen = warningfileIsOpen; // Simplified
    std::cout << itsopen << "- Cannot write into warning file the message: " << buff3 << std::endl;

label_155:
    return;
}

bool isFatalError() {
    return fatalerror;
}

void resetFatalError() {
    fatalerror = false;
}

void CLOSEWARNINGFILE(int layoutnumber, int num_procs, bool fatalerror_final, bool stoch_undivided, bool simu_devia) {
    // integer(kind=4), intent(in) :: layoutnumber,num_procs
    // integer(kind=4) :: ierr,posic,i
    // character(len=BUFSIZE) :: buf2
    // character(len=BUFSIZE) :: dubuf
    // logical :: fatalerror_final , lexis,stoch_undivided,simu_devia        , itsopen2
    // character( LEN=BUFSIZE) :: whoamishort,whoami,chinstant
    // integer :: my_iostat,file87
    // character(len=BUFSIZE) :: ficherito
    // if (.not.WarningFileIsOpen) return
    // !!!#ifdef CompileWithMPI
    // !!!if (SIZE/=0) then
    // !!!    call MPI_FILE_close (thefile, ierr)
    // !!!ELSE
    // !!!    close (17)
    // !!!end if
    // !!!#else
    // close (17)
    // !!!#endif
    // #ifdef CompileWithMPI
    // !wait until everything is closed
    // call MPI_Barrier (SUBCOMM_MPI, ierr)
    // #endif
    // arregla los NUL
    // if ((layoutnumber==0).or.((layoutnumber == num_procs/2).and.stoch_undivided)) then
    //    open (88,file=trim(adjustl(WarningFile))//'_Warnings.txt',form='formatted')
    //    posic=0
    //    do i=0,num_procs-1
    //       if (stoch_undivided) then
    //           write(whoamishort,'(i5)') i+1
    //       else
    //          if (simu_devia) then
    //              write(whoamishort,'(i5)') num_procs+i+1
    //          else
    //              write(whoamishort,'(i5)') i+1
    //          end if
    //       end if 
    //       inquire(file=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt',exist=lexis)
    //       if (lexis) then         
    // !!!inquire(unit=87, opened=itsopen2)
    // !!!if (itsopen2) print *,'----------->87 open!!!'
    //          ficherito=trim(adjustl(WarningFile))//trim(adjustl(whoamishort))//'_tmpWarnings.txt'
    //          call opensolo(87,ficherito)
    //          !
    // 875            read(87,'(a)',end=876,err=876) buf2
    //          call trimnullchar(buf2)
    //          buf2=trim(adjustl(buf2))
    //          if ((buf2(1:1) /= ' ').and.(buf2(1:1) /=char(0))) then
    //             write(88,'(a)') trim(adjustl(buf2))
    //             posic=posic+1
    //          end if
    //          goto 875
    //          !
    // 876            continue
    //          call closesolo(87)
    // #ifndef CorregirBugBorrado

    if (!warningfileIsOpen) return;

    // close (17)
    // Assuming unit 17 is closed.
    // In C++, we'd close the stream.
    // For this translation, we assume the stream is closed.

#ifdef CompileWithMPI
    // wait until everything is closed
    MPI_Barrier(SUBCOMM_MPI, &ierr);
#endif

    // arregla los NUL
    if ((layoutnumber == 0) || ((layoutnumber == num_procs / 2) && stoch_undivided)) {
        // open (88,file=trim(adjustl(WarningFile))//'_Warnings.txt',form='formatted')
        std::string filename = WarningFile + "_Warnings.txt";
        std::ofstream file88(filename);
        if (!file88) {
            std::cerr << "Error opening file: " << filename << std::endl;
            return;
        }

        int posic = 0;
        for (int i = 0; i < num_procs; ++i) {
            std::string whoamishort;
            char buf[10];
            if (stoch_undivided) {
                snprintf(buf, sizeof(buf), "%5d", i + 1);
            } else {
                if (simu_devia) {
                    snprintf(buf, sizeof(buf), "%5d", num_procs + i + 1);
                } else {
                    snprintf(buf, sizeof(buf), "%5d", i + 1);
                }
            }
            whoamishort = buf;

            std::string tmpFilename = WarningFile + whoamishort + "_tmpWarnings.txt";
            std::ifstream file87(tmpFilename);
            if (file87.is_open()) {
                std::string buf2;
                while (std::getline(file87, buf2)) {
                    trimnullchar(buf2);
                    // trim(adjustl(buf2)) is essentially removing leading spaces
                    // In C++, we can use a helper function or just check the first char.
                    // The Fortran code checks if the first char is not space and not null.
                    if (!buf2.empty() && buf2[0] != ' ' && buf2[0] != '\0') {
                        file88 << buf2 << std::endl;
                        posic++;
                    }
                }
                file87.close();
            }
        }
        file88.close();
    }
    // #ifndef CorregirBugBorrado
    // The code ends here in the snippet.
}

} // namespace FortranModule

ficherito = trim(adjustl(WarningFile)) + trim(adjustl(whoamishort)) + "_tmpWarnings.txt";
            openclosedelete(ficherito);
#endif
         }
      }

      //
      dubuf = SEPARADOR + separador + separador;
      print11(layoutnumber, dubuf);
      dubuf = "Closing warning file. Number of messages: " + std::to_string(posic);
      print11(layoutnumber, dubuf);
      dubuf = SEPARADOR + separador + separador;
      print11(layoutnumber, dubuf);

      close_file(88);
   }

   warningfileIsOpen = false;
#ifdef CompileWithMPI
   MPI_Barrier(SUBCOMM_MPI, ierr);
   MPI_AllReduce(fatalerror, fatalerror_final, 1, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr);
#else
   fatalerror_final = fatalerror;
#endif

   if (fatalerror_final && ignoreErrors) {
      dubuf = SEPARADOR + separador + separador;
      print11(layoutnumber, dubuf);
      dubuf = "There are ERRORS: The simulation will continue but Revise *Warnings file ";
      print11(layoutnumber, dubuf);
      dubuf = SEPARADOR + separador + separador;
      print11(layoutnumber, dubuf);
   }

   fatalerror_final = (fatalerror_final) && (!ignoreErrors);

   return;
}

void trimnullchar(std::string& string) {
   int i, longi;
   int ind;
   longi = string.length();
   for (i = 0; i < bufsize; ++i) {
      // scan is not directly available in std::string, simulating logic
      // scan(string(i:longi), char(0)) finds the first occurrence of null char in substring
      // Note: Fortran strings are 1-indexed, C++ 0-indexed.
      // string(i:longi) in Fortran corresponds to substring starting at index i (0-based) up to end.
      // We need to find char(0) which is '\0'.
      size_t found = string.find('\0', i);
      if (found != std::string::npos) {
         // ind is the position relative to the start of the substring (1-based in Fortran logic usually, but here used for assignment)
         // In Fortran: string(ind:ind) = ' '. ind is the index in the original string where char(0) was found.
         // Since we search from i, and string is 0-indexed in C++, found is the correct 0-based index.
         string[found] = ' ';
      }
   }
   string = trim(adjustl(string));
}

void print11(int layoutnumber, const std::string& message, bool forceprint2) {
   bool soyImpresor, forceprint;
   forceprint = false;
   if (forceprint2) forceprint = forceprint2;

   soyImpresor = ((layoutnumber == 0) || forceprint) && printea;
   if (message.length() > 0 && message[0] == '&') { // respeta los espacios, el & es un espacio en realidad
#ifndef NoVerbose
      if (soyImpresor) {
         std::cout << trim(message.substr(1)) << std::endl;
      }
#endif
      if (layoutnumber == 0) {
         try {
            file11 << trim(message.substr(1)) << std::endl;
         } catch (...) {
            goto label_112;
         }
      }
   } else { // ajusta a izquierda sin respetar espacios
#ifndef NoVerbose
      if (soyImpresor) {
         std::cout << trim(adjustl(message)) << std::endl;
      }
#endif
      if (layoutnumber == 0) {
         try {
            file11 << trim(adjustl(message)) << std::endl;
         } catch (...) {
            goto label_111;
         }
      }
   }
   goto label_112;
label_111:
   // fort.11 a veces lo intentan escribir 2 a la vez de los que dan fallos en writing restarting fields.
   // asi que ignora y continua
label_112:
   return;
}

// Note: The following Fortran code is commented out.
    // It is preserved here as comments for reference, but not translated into executable C++.
    // The original Fortran code handles DXF file writing for PECs (Perfect Electric Conductors).

    /*
    !!!        if (lexis) then
    !!!            open (987,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf',form='formatted')
    !!!!
    !!!875         read(987,'(a)',end=876,err=876) buf2
    !!!            call trimnullchar(buf2)
    !!!            buf2=trim(adjustl(buf2))
    !!!            if ((buf2(1:1) /= ' ').and.(buf2(1:1) /=char(0))) then
    !!!                write(988,'(a)') trim(adjustl(buf2))
    !!!                posic=posic+1
    !!!            end if
    !!!            goto 875
    !!!!
    !!!876         close (987)
    !!!!
    !!!            open (987,file=trim(adjustl(mynEntradaRoot))//trim(adjustl(whoamishort))//'.tmpdxf')
    !!!            write(987,*) '!END'
    !!!            close (987,status='delete')
    !!!        end if
    !!!    end do
    !!!!end the file
    !!!    write(988,'(a)')  'ENDSEC'
    !!!    write(988,'(a)')  '0'
    !!!    write(988,'(a)')  'EOF'
    !!!    close (988)
    !!!
    !!!
    !!!
    !!!    write(dubuf,*)SEPARADOR//separador//separador
    !!!    call print11(layoutnumber,dubuf)
    !!!    write(dubuf,*) 'Closing dxf file. Number of lines: ',posic
    !!!    call print11(layoutnumber,dubuf)
    !!!end if
    !!!
    !!!continue
    !!!
    !!!dxfFileIsOpen=.false.
    !!!return
    !!!end subroutine

    !!!
    !!!subroutine writemmdxf(layoutnumber,sgg,sggMiHx,sggMiHy,sggMiHz)
    !!!integer i,j,k,layoutnumber
    !!!type(SGGFDTDINFO_t), intent(in) :: sgg
    !!!integer(kind=INTEGERSIZEOFMEDIAMATRICES), intent(in) :: &
    !!!           sggMiHx(sgg%Alloc(iHx)%XI : sgg%Alloc(iHx)%XE, &
    !!!                   sgg%Alloc(iHx)%YI : sgg%Alloc(iHx)%YE, &
    !!!                   sgg%Alloc(iHx)%ZI : sgg%Alloc(iHx)%ZE), &
    !!!           sggMiHy(sgg%Alloc(iHy)%XI : sgg%Alloc(iHy)%XE, &
    !!!                   sgg%Alloc(iHy)%YI : sgg%Alloc(iHy)%YE, &
    !!!                   sgg%Alloc(iHy)%ZI : sgg%Alloc(iHy)%ZE), &
    !!!           sggMiHz(sgg%Alloc(iHz)%XI : sgg%Alloc(iHz)%XE, &
    !!!                   sgg%Alloc(iHz)%YI : sgg%Alloc(iHz)%YE, &
    !!!                   sgg%Alloc(iHz)%ZI : sgg%Alloc(iHz)%ZE)
    !!!
    !!!!write ONLY PECS
    !!!      Do k=sgg%SINPMLSweep(iHx)%ZI,sgg%SINPMLSweep(iHx)%ZE
    !!!          Do j=sgg%SINPMLSweep(iHx)%YI,sgg%SINPMLSweep(iHx)%YE
    !!!              Do i=sgg%SINPMLSweep(iHx)%XI,sgg%SINPMLSweep(iHx)%XE
    !!!                if ((sggMiHx(i,j,k) ==0).or.(sgg%med(sggMiHx(i,j,k) )%is%pec))  then
    !!!                    write(dxfbuff,'(a)') '3DFACE'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '8'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') sggMiHx(i,j,k)+20
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '62'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') sggMiHx(i,j,k)+20
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '10'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '20'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '30'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') k
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '11'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '21'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j+1
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '31'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') k
    !!!                    call DXFWRITE(DXFBUFF)
    !!!
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '12'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '22'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j+1
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '32'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') k+1
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '13'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '23'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '33'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') k+1
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '0'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                end if
    !!!            end do
    !!!        end do
    !!!    end do
    !!!!
    !!!      Do k=sgg%SINPMLSweep(iHy)%ZI,sgg%SINPMLSweep(iHy)%ZE
    !!!          Do j=sgg%SINPMLSweep(iHy)%YI,sgg%SINPMLSweep(iHy)%YE
    !!!              Do i=sgg%SINPMLSweep(iHy)%XI,sgg%SINPMLSweep(iHy)%XE
    !!!                if ((sggMiHy(i,j,k) ==0).or.(sgg%med(sggMiHy(i,j,k) )%is%pec))  then
    !!!                    write(dxfbuff,'(a)') '3DFACE'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '8'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') sggMiHy(i,j,k)+20
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '62'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') sggMiHy(i,j,k)+20
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '10'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '20'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '30'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') k
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '11'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '21'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '31'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') k+1
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    !
    !!!                    write(dxfbuff,'(a)') '12'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') i+1
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(a)') '22'
    !!!                    call DXFWRITE(DXFBUFF)
    !!!                    write(dxfbuff,'(i5)') j
    */

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstring>
#include <unistd.h>

// Assuming BUFSIZE is defined elsewhere or needs to be defined. 
// Based on Fortran character(len=BUFSIZE), we use std::string or a fixed char array.
// For simplicity and safety with trim/adjustl logic, std::string is preferred, 
// but to match "preserve names" and typical Fortran interop, we might use a struct or just string.
// The prompt asks to convert arrays to vector, but these are scalars.
// We will use std::string for character buffers.

#ifndef BUFSIZE
#define BUFSIZE 256
#endif

// Forward declarations for MPI if needed, though we are stripping MPI logic or keeping it conditional.
// The prompt says "Convert #ifdef to C++". We will keep the preprocessor directives.

#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI; // Assumed external global from Fortran context
#endif

// Helper to simulate Fortran's trim(adjustl())
std::string trim_adjustl(const std::string& str) {
    size_t start = str.find_first_not_of(" ");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" ");
    return str.substr(start, end - start + 1);
}

// Helper to simulate Fortran's CHAR(13)//CHAR(10)
std::string newlines() {
    return "\r\n";
}

// Helper to simulate trimnullchar
void trimnullchar(std::string& str) {
    size_t end = str.find_first_of('\0');
    if (end != std::string::npos) {
        str.erase(end);
    }
}

// Function: openfile_mpi
int openfile_mpi(int layoutnumber, const std::string& nombrefich) {
    int thefile8 = 0;
    int iter = 0;
    int ios = 0;
    std::string whoamishort;
    bool file_exists = false;

    // write(whoamishort, '(i5)') layoutnumber + 1
    char buf[10];
    snprintf(buf, sizeof(buf), "%5d", layoutnumber + 1);
    whoamishort = buf;

    for (iter = 1; iter <= 10; ++iter) {
        ios = 0;
        std::string filename = trim_adjustl(nombrefich) + trim_adjustl(whoamishort) + "_tmp";
        
        // open(newunit=thefile8, file=..., iostat=ios)
        // In C++, we can't easily get an "iostat" from open() directly in the same way.
        // We'll try to open. If it fails, ios != 0.
        // Note: Fortran open creates the file if it doesn't exist by default in some modes, 
        // but here it seems to be checking for existence or race conditions.
        // We will use a simple try-catch or check stream state.
        
        std::ofstream file_stream(filename, std::ios::out | std::ios::trunc);
        if (file_stream.is_open()) {
            thefile8 = 1; // Simulate valid unit number
            file_stream.close();
            // Remove the file immediately so we can check if it was created successfully?
            // Actually, the Fortran code opens it, checks ios. If ios==0, it's open.
            // Then it closes it if layoutnumber==0.
            // Let's simulate the unit number as a file descriptor or just an int.
            // For this translation, we'll return a file descriptor or a unique int.
            // Let's use a simple counter or just return 1 if successful.
            // To be more robust, let's use a temporary file approach.
            
            // Re-open for writing later? The function returns the unit number.
            // We need to keep the file open or reopen it in the caller.
            // The caller `writefile_mpi` takes `thefile8`.
            // So `openfile_mpi` should return a handle.
            // Let's return a file pointer or descriptor. But the signature says `integer(4)`.
            // We will return a file descriptor (int).
            
            // Let's restart the logic to return a file descriptor.
            // We already opened and closed it above to check existence/creation.
            // Now we need to return an open file.
            
            // Actually, the Fortran code:
            // open(newunit=thefile8, ...)
            // if (ios == 0) ...
            // So thefile8 is the unit number.
            // We will return a file descriptor.
            
            // Let's redo the open properly.
            // We need to return an open file stream or fd.
            // Since the rest of the code uses `write(thefile8, ...)`, it expects a unit number.
            // We will map unit numbers to std::ofstream objects in a global map or similar, 
            // OR we just return a file descriptor and use C-style I/O.
            // Given the instruction "Convert subroutines to C++ functions", and "real->double", 
            // let's stick to C++ streams but we need a way to pass the "unit".
            // A common pattern is a map<int, std::ofstream>.
            
            // However, to keep it simple and self-contained in this chunk, 
            // I will assume a global map or pass the stream by reference in other functions?
            // No, the signature is fixed: `int openfile_mpi(...)`.
            
            // Let's use a simple file descriptor approach for `thefile8`.
            // We already created the file. Let's open it for writing.
            // But wait, the loop checks `ios`. If `ios != 0`, it loops.
            // If `ios == 0`, it exits.
            
            // Let's implement the loop correctly.
        } else {
            ios = 1; // Error
        }
        
        // If we are here, we need to decide if we exit.
        // The Fortran code:
        // if (ios /= 0) print error
        // if (layoutnumber == 0 .and. ios == 0) { close and delete }
        // call sleep(2)
        // if (ios == 0) exit
        
        // Let's re-implement the loop logic cleanly.
    }
    
    // Re-writing the loop logic properly
    for (iter = 1; iter <= 10; ++iter) {
        ios = 0;
        std::string filename = trim_adjustl(nombrefich) + trim_adjustl(whoamishort) + "_tmp";
        
        // Try to open for writing (creating/truncating)
        std::ofstream temp_file(filename, std::ios::out | std::ios::trunc);
        if (temp_file.is_open()) {
            thefile8 = 1; // Placeholder unit number
            temp_file.close();
            
            if ((layoutnumber == 0) && (ios == 0)) {
                sleep(2);
                // write(thefile8, '(a)') '!END'
                // We need to write to the file. Since we closed it, we reopen.
                // But wait, the Fortran code writes to `thefile8` which is currently open.
                // So we shouldn't close it yet if we are going to write to it.
                // But the code says:
                // open(...)
                // if (ios==0) ... write ... close ...
                // So it opens, writes, closes.
                
                // Let's reopen for writing to put '!END'
                std::ofstream write_file(filename, std::ios::out | std::ios::trunc);
                if (write_file.is_open()) {
                    write_file << "!END" << std::endl;
                    write_file.close();
                }
                sleep(2);
                ios = 0;
                // close(thefile8, status='delete')
                std::remove(filename.c_str());
                if (std::remove(filename.c_str()) != 0) {
                    std::cout << "Error deleting temporary file: " << filename << std::endl;
                }
            }
        } else {
            ios = 1;
            std::cout << "Error opening temporary file: " << filename << std::endl;
        }
        
        sleep(2);
        if (ios == 0) {
            // We need to return an OPEN file handle for the caller.
            // The caller expects `thefile8` to be a valid unit/stream.
            // We will return a file descriptor.
            // Let's open the file for writing now.
            // But wait, if layoutnumber==0, we deleted the file.
            // So we can't open it.
            // This implies the loop continues until it finds a file that exists and isn't deleted?
            // Or maybe the deletion is only for specific cases.
            
            // If layoutnumber != 0, the file remains.
            // We need to return an open file.
            // Let's assume the caller will handle the file opening?
            // No, `openfile_mpi` returns `thefile8`.
            
            // Let's change strategy: Return a file descriptor.
            // If layoutnumber == 0, the file was deleted, so we can't return it.
            // But the code says `if (ios == 0) exit`.
            // If layoutnumber == 0, it deletes the file.
            // Then it sleeps. Then it exits.
            // So `thefile8` is returned, but the file is gone?
            // This seems like a synchronization mechanism.
            // The caller `writefile_mpi` uses `thefile8`.
            // If the file is deleted, `writefile_mpi` will fail.
            
            // Let's look at `writefile_mpi`. It takes `thefile8`.
            // It writes to it.
            
            // Maybe the `openfile_mpi` is meant to create a semaphore file?
            // And the actual data file is opened elsewhere?
            // No, `writefile_mpi` writes to `thefile8`.
            
            // Let's assume that for `layoutnumber == 0`, the file is a marker.
            // And for other processes, they open the file for appending?
            
            // To make this compile and run logically in C++, we will return a file descriptor.
            // If the file was deleted, we return -1 or handle it.
            // But the Fortran code doesn't check for that.
            
            // Let's assume the file is NOT deleted for the final return.
            // The `if (layoutnumber == 0)` block deletes it.
            // So if layoutnumber is 0, the file is gone.
            // This suggests `openfile_mpi` is only called by non-zero processes for data?
            // Or maybe the deletion is a cleanup step after the process is done?
            
            // Let's just return a file descriptor for the file `filename`.
            // If layoutnumber == 0, we return a dummy or handle the error.
            
            // For the sake of translation, we will return a file descriptor.
            // We need to open the file.
            int fd = open(filename.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
            if (fd < 0) {
                ios = 1;
                continue;
            }
            // We need to store this fd somewhere because we can't return it directly if we want to use C++ streams later?
            // No, `writefile_mpi` will use `thefile8`.
            // We can use a global map: `std::map<int, int> unit_to_fd;`
            // But we don't have access to global state in this snippet.
            
            // Alternative: Return the file descriptor directly.
            // And modify `writefile_mpi` to use `write(fd, ...)` instead of `write(unit, ...)`.
            // But the prompt says "Convert subroutines to C++ functions".
            // We must preserve the interface as much as possible.
            
            // Let's use a global map for unit numbers to file descriptors.
            // This is a common way to simulate Fortran unit numbers in C++.
            
            // Define a global map
            static std::map<int, int> unit_map;
            unit_map[thefile8] = fd;
            
            return thefile8;
        }
    }
    
    return 0; // Error
}

// We need a helper to get the file descriptor from the unit number
int get_fd_from_unit(int unit) {
    static std::map<int, int> unit_map;
    if (unit_map.find(unit) != unit_map.end()) {
        return unit_map[unit];
    }
    return -1;
}

// Function: writefile_mpi
void writefile_mpi(int layoutnumber, int thefile8, const std::string& buff2) {
    int ierr = 0;
    std::string buff3 = newlines() + trim_adjustl(buff2);
    
#ifdef CompileWithMPI
    // MPI_FILE_WRITE is complex to translate directly without MPI_File type.
    // We will skip the MPI part as per "Convert #ifdef to C++" usually implying 
    // keeping the structure but adapting the code.
    // However, the prompt says "Convert #ifdef to C++".
    // We will keep the preprocessor.
    
    // MPI_FILE_WRITE(thefile8, buff3, 1024_4, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
    // This requires MPI_File. The Fortran `thefile8` is an integer unit.
    // In MPI, you usually get an MPI_File handle.
    // This translation is tricky. We will assume the MPI path is not taken 
    // or we stub it out.
    
    // For the non-MPI path:
    write(thefile8, buff2);
#else
    write(thefile8, buff2);
#endif
}

// Helper for writefile_mpi to write to a unit number
void write(int unit, const std::string& data) {
    int fd = get_fd_from_unit(unit);
    if (fd >= 0) {
        write(fd, data.c_str(), data.size());
    }
}

// Function: closefile_mpi
void closefile_mpi(int layoutnumber, int num_procs, const std::string& nombrefich, int thefile8) {
    int thefile19 = 0;
    int ierr = 0;
    int conta = 0;
    int i = 0;
    std::string whoamishort;
    bool lexis = false;
    std::string ficherito;
    std::string buff2;

#ifdef CompileWithMPI
    call MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

    // close(thefile8)
    // We need to close the file descriptor associated with thefile8
    int fd = get_fd_from_unit(thefile8);
    if (fd >= 0) {
        close(fd);
        // Remove from map
        static std::map<int, int> unit_map;
        unit_map.erase(thefile8);
    }

#ifdef CompileWithMPI
    call MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

    if (layoutnumber == 0) {
        // open(newunit=thefile8, file=..., form='formatted')
        // We need to open a new file for writing the final result.
        // Let's reuse thefile8 for the new file?
        // The Fortran code reuses `thefile8`.
        
        std::string final_filename = trim_adjustl(nombrefich);
        std::ofstream final_file(final_filename, std::ios::out | std::ios::trunc);
        if (!final_file.is_open()) {
            std::cout << "Error opening final file: " << final_filename << std::endl;
            return;
        }
        
        // Store the new file descriptor in the unit map for thefile8
        int new_fd = final_file.rdbuf()->fd(); // This is non-standard.
        // Better to use C-style I/O for the final file too.
        
        // Let's switch to C-style I/O for the final file to match the read/write logic.
        FILE* fp = fopen(final_filename.c_str(), "w");
        if (!fp) {
            std::cout << "Error opening final file: " << final_filename << std::endl;
            return;
        }
        
        // Update the unit map
        static std::map<int, int> unit_map;
        unit_map[thefile8] = fileno(fp); // Store FILE* fd
        
        for (i = 0; i < num_procs; ++i) {
            char buf[10];
            snprintf(buf, sizeof(buf), "%5d", i + 1);
            whoamishort = buf;
            
            std::string tmp_filename = trim_adjustl(nombrefich) + trim_adjustl(whoamishort) + "_tmp";
            
            // inquire(file=..., exist=lexis)
            struct stat buffer;
            lexis = (stat(tmp_filename.c_str(), &buffer) == 0);
            
            if (lexis) {
                // open(newunit=thefile19, file=..., form='formatted')
                FILE* fp_tmp = fopen(tmp_filename.c_str(), "r");
                if (!fp_tmp) {
                    std::cout << "Error opening temp file: " << tmp_filename << std::endl;
                    continue;
                }
                
                // Store thefile19 in unit map?
                // We need a unique unit number for thefile19.
                // Let's use a counter or just a local variable.
                // But `read(thefile19, ...)` uses the unit number.
                // We need to map `thefile19` to `fp_tmp`.
                // Let's assume `thefile19` is a new unit number.
                // We'll generate one.
                static int next_unit = 100;
                int unit19 = next_unit++;
                static std::map<int, FILE*> unit_file_map;
                unit_file_map[unit19] = fp_tmp;
                
                conta = 0;
                
                // 875 continue
                while (true) {
                    // read(thefile19,'(a)',end=876,err=876) buff2
                    // Fortran read line.
                    char line[BUFSIZE];
                    if (fgets(line, sizeof(line), fp_tmp) == nullptr) {
                        break; // End of file
                    }
                    
                    // Remove newline
                    size_t len = strlen(line);
                    if (len > 0 && line[len-1] == '\n') {
                        line[len-1] = '\0';
                    }
                    if (len > 1 && line[len-2] == '\r') {
                        line[len-2] = '\0';
                    }
                    
                    buff2 = line;
                    
                    // call trimnullchar(buff2)
                    trimnullchar(buff2);
                    
                    // buff2=trim(adjustl(buff2))
                    buff2 = trim_adjustl(buff2);
                    
                    // if ((buff2(1:1) /= ' ').and.(buff2(1:1) /=char(0)))
                    if (buff2.length() > 0 && buff2[0] != ' ' && buff2[0] != '\0') {
                        // write(thefile8,'(a)') trim(adjustl(buff2))
                        fprintf(fp, "%s\n", trim_adjustl(buff2).c_str());
                        conta++;
                    }
                }
                
                // 876 continue
                fclose(fp_tmp);
                
                // Remove thefile19 from map
                unit_file_map.erase(unit19);
            }
        }
        
        fclose(fp);
    }
}

ficherito = trim(adjustl(nombrefich)) + trim(adjustl(whoamishort)) + "_tmp";
            openclosedelete(ficherito);
         }
      }
      close(thefile8);

      //
   }
#endif
#ifdef CompileWithMPI
   MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
   // !!!!end gnuplot
   return;
}

// Subroutine closefile_mpi ends here (implicit from context, but we translate the function below)

coorsxyzP_t creaPuntos(SGGFDTDINFO_t& sgg) { // crea coordenadas fisicas
   //
   // type(SGGFDTDINFO_t), intent(INout) :: sgg
   // type(coorsxyzP_t) :: Punto
   // integer(Kind=4) :: i,j,k,field

   coorsxyzP_t Punto;
   int i, j, k, field;

   for (field = iEx; field <= iHz; ++field) {
      // Assuming sgg.Sweep(field).XI, XE, etc. are accessible
      // Allocate vectors with size (XE + 1) - (XI - 1) + 1 = XE - XI + 3
      // Fortran: XI-1 : XE+1. Size = (XE+1) - (XI-1) + 1 = XE - XI + 3.
      // In C++, we usually use 0-based indexing or adjust. 
      // To preserve names and logic strictly, we might keep a wrapper or adjust indices.
      // However, std::vector is 0-based. 
      // Let's assume the struct members x, y, z are std::vector<double> or similar.
      // The Fortran code accesses indices from XI-1 to XE+1.
      // We will allocate enough space and map indices.
      // For simplicity in translation without full struct definition, we assume 
      // the underlying storage allows this or we resize to cover the range.
      
      int sizeX = sgg.Sweep[field].XE + 1 - (sgg.Sweep[field].XI - 1) + 1;
      int sizeY = sgg.Sweep[field].YE + 1 - (sgg.Sweep[field].YI - 1) + 1;
      int sizeZ = sgg.Sweep[field].ZE + 1 - (sgg.Sweep[field].ZI - 1) + 1;

      Punto.PhysCoor[field].x.resize(sizeX);
      Punto.PhysCoor[field].y.resize(sizeY);
      Punto.PhysCoor[field].z.resize(sizeZ);

      // Initialize with -1e20
      for (int idx = 0; idx < sizeX; ++idx) Punto.PhysCoor[field].x[idx] = -1e20;
      for (int idx = 0; idx < sizeY; ++idx) Punto.PhysCoor[field].y[idx] = -1e20;
      for (int idx = 0; idx < sizeZ; ++idx) Punto.PhysCoor[field].z[idx] = -1e20;
   }

   //
   //
   field = iEx;
   for (i = sgg.SINPMLSweep[field].XI - 1; i <= sgg.SINPMLSweep[field].XE + 1; ++i) {
      // Map Fortran index i to vector index. 
      // If vector is 0-based, index = i - (min_index). 
      // Min index for x is sgg.Sweep[field].XI - 1.
      int idx = i - (sgg.Sweep[field].XI - 1);
      Punto.PhysCoor[field].x[idx] = (sgg.LineX[i] + sgg.LineX[i + 1]) * 0.5_RKIND;
   }
   for (j = sgg.SINPMLSweep[field].YI - 1; j <= sgg.SINPMLSweep[field].YE + 1; ++j) {
      int idx = j - (sgg.Sweep[field].YI - 1);
      Punto.PhysCoor[field].y[idx] = sgg.LineY[j];
   }
   for (k = sgg.SINPMLSweep[field].ZI - 1; k <= sgg.SINPMLSweep[field].ZE + 1; ++k) {
      int idx = k - (sgg.Sweep[field].ZI - 1);
      Punto.PhysCoor[field].z[idx] = sgg.LineZ[k];
   }

   field = iEy;
   for (i = sgg.SINPMLSweep[field].XI - 1; i <= sgg.SINPMLSweep[field].XE + 1; ++i) {
      int idx = i - (sgg.Sweep[field].XI - 1);
      Punto.PhysCoor[field].x[idx] = sgg.LineX[i];
   }
   for (j = sgg.SINPMLSweep[field].YI - 1; j <= sgg.SINPMLSweep[field].YE + 1; ++j) {
      int idx = j - (sgg.Sweep[field].YI - 1);
      Punto.PhysCoor[field].y[idx] = (sgg.LineY[j] + sgg.LineY[j + 1]) * 0.5_RKIND;
   }
   for (k = sgg.SINPMLSweep[field].ZI - 1; k <= sgg.SINPMLSweep[field].ZE + 1; ++k) {
      int idx = k - (sgg.Sweep[field].ZI - 1);
      Punto.PhysCoor[field].z[idx] = sgg.LineZ[k];
   }

   field = iEz;
   for (i = sgg.SINPMLSweep[field].XI - 1; i <= sgg.SINPMLSweep[field].XE + 1; ++i) {
      int idx = i - (sgg.Sweep[field].XI - 1);
      Punto.PhysCoor[field].x[idx] = sgg.LineX[i];
   }
   for (j = sgg.SINPMLSweep[field].YI - 1; j <= sgg.SINPMLSweep[field].YE + 1; ++j) {
      int idx = j - (sgg.Sweep[field].YI - 1);
      Punto.PhysCoor[field].y[idx] = sgg.LineY[j];
   }
   for (k = sgg.SINPMLSweep[field].ZI - 1; k <= sgg.SINPMLSweep[field].ZE + 1; ++k) {
      int idx = k - (sgg.Sweep[field].ZI - 1);
      Punto.PhysCoor[field].z[idx] = (sgg.LineZ[k] + sgg.LineZ[k + 1]) * 0.5_RKIND;
   }

   field = iHx;
   for (i = sgg.SINPMLSweep[field].XI - 1; i <= sgg.SINPMLSweep[field].XE + 1; ++i) {
      int idx = i - (sgg.Sweep[field].XI - 1);
      Punto.PhysCoor[field].x[idx] = sgg.LineX[i];
   }
   for (j = sgg.SINPMLSweep[field].YI - 1; j <= sgg.SINPMLSweep[field].YE + 1; ++j) {
      int idx = j - (sgg.Sweep[field].YI - 1);
      Punto.PhysCoor[field].y[idx] = (sgg.LineY[j] + sgg.LineY[j + 1]) * 0.5_RKIND;
   }
   for (k = sgg.SINPMLSweep[field].ZI - 1; k <= sgg.SINPMLSweep[field].ZE + 1; ++k) {
      int idx = k - (sgg.Sweep[field].ZI - 1);
      Punto.PhysCoor[field].z[idx] = (sgg.LineZ[k] + sgg.LineZ[k + 1]) * 0.5_RKIND;
   }

   field = iHy;
   for (i = sgg.SINPMLSweep[field].XI - 1; i <= sgg.SINPMLSweep[field].XE + 1; ++i) {
      int idx = i - (sgg.Sweep[field].XI - 1);
      Punto.PhysCoor[field].x[idx] = (sgg.LineX[i] + sgg.LineX[i + 1]) * 0.5_RKIND;
   }
   for (j = sgg.SINPMLSweep[field].YI - 1; j <= sgg.SINPMLSweep[field].YE + 1; ++j) {
      int idx = j - (sgg.Sweep[field].YI - 1);
      Punto.PhysCoor[field].y[idx] = sgg.LineY[j];
   }
   for (k = sgg.SINPMLSweep[field].ZI - 1; k <= sgg.SINPMLSweep[field].ZE + 1; ++k) {
      int idx = k - (sgg.Sweep[field].ZI - 1);
      Punto.PhysCoor[field].z[idx] = (sgg.LineZ[k] + sgg.LineZ[k + 1]) * 0.5_RKIND;
   }

   field = iHz;
   for (i = sgg.SINPMLSweep[field].XI - 1; i <= sgg.SINPMLSweep[field].XE + 1; ++i) {
      int idx = i - (sgg.Sweep[field].XI - 1);
      Punto.PhysCoor[field].x[idx] = (sgg.LineX[i] + sgg.LineX[i + 1]) * 0.5_RKIND;
   }
   for (j = sgg.SINPMLSweep[field].YI - 1; j <= sgg.SINPMLSweep[field].YE + 1; ++j) {
      int idx = j - (sgg.Sweep[field].YI - 1);
      Punto.PhysCoor[field].y[idx] = (sgg.LineY[j] + sgg.LineY[j + 1]) * 0.5_RKIND;
   }
   for (k = sgg.SINPMLSweep[field].ZI - 1; k <= sgg.SINPMLSweep[field].ZE + 1; ++k) {
      int idx = k - (sgg.Sweep[field].ZI - 1);
      Punto.PhysCoor[field].z[idx] = sgg.LineZ[k];
   }
   //
   sgg.Punto = Punto;
   return Punto;
}

void reportmedia(SGGFDTDINFO_t& sgg) {
   int j;
   std::string buff;

   for (j = 0; j <= sgg.NumMedia; ++j) {
      buff = "_____________________________";
      WarnErrReport(buff);
      buff = "MEDIO :  " + std::to_string(j);
      WarnErrReport(buff);
      buff = "Priority " + std::to_string(sgg.Med[j].Priority);
      WarnErrReport(buff);
      buff = "Epr " + std::to_string(sgg.Med[j].Epr);
      WarnErrReport(buff);
      buff = "Sigma " + std::to_string(sgg.Med[j].Sigma);
      WarnErrReport(buff);
      buff = "Mur " + std::to_string(sgg.Med[j].Mur);
      WarnErrReport(buff);
      buff = "SigmaM " + std::to_string(sgg.Med[j].SigmaM);
      WarnErrReport(buff);
      buff = "Is PML " + std::to_string(sgg.Med[j].Is.PML);
      WarnErrReport(buff);
      buff = "Is PEC " + std::to_string(sgg.Med[j].Is.PEC);
      WarnErrReport(buff);
      buff = "Is ThinWIRE " + std::to_string(sgg.Med[j].Is.ThinWire);
      WarnErrReport(buff);
      buff = "Is MULTIWIRE " + std::to_string(sgg.Med[j].Is.Multiwire);
      WarnErrReport(buff);
      buff = "Is SlantedWIRE " + std::to_string(sgg.Med[j].Is.SlantedWire);
      WarnErrReport(buff);
      buff = "Is EDispersive " + std::to_string(sgg.Med[j].Is.EDispersive);
      WarnErrReport(buff);
      buff = "Is EDispersiveaANIS " + std::to_string(sgg.Med[j].Is.EDispersiveANIS);
      WarnErrReport(buff);
      buff = "Is MDispersive " + std::to_string(sgg.Med[j].Is.MDispersive);
      WarnErrReport(buff);
      buff = "Is MDispersiveANIS " + std::to_string(sgg.Med[j].Is.MDispersiveANIS);
      WarnErrReport(buff);
      buff = "Is ThinSlot " + std::to_string(sgg.Med[j].Is.ThinSlot);
      WarnErrReport(buff);
      buff = "Is SGBC " + std::to_string(sgg.Med[j].Is.SGBC);
      WarnErrReport(buff);
      buff = "Is Lossy " + std::to_string(sgg.Med[j].Is.Lossy);
      WarnErrReport(buff);
      buff = "Is Multiport " + std::to_string(sgg.Med[j].Is.multiport);
      WarnErrReport(buff);
      buff = "Is AnisMultiport " + std::to_string(sgg.Med[j].Is.anismultiport);
      WarnErrReport(buff);
      buff = "Is MultiportPadding " + std::to_string(sgg.Med[j].Is.multiportpadding);
      WarnErrReport(buff);
      buff = "Is Dielectric " + std::to_string(sgg.Med[j].Is.dielectric);
      WarnErrReport(buff);
      buff = "Is ThinSlot " + std::to_string(sgg.Med[j].Is.ThinSlot);
      WarnErrReport(buff);
      buff = "Is Anisotropic " + std::to_string(sgg.Med[j].Is.Anisotropic);
      WarnErrReport(buff);
      buff = "Is Needed " + std::to_string(sgg.Med[j].Is.Needed);
      WarnErrReport(buff);
      buff = "Is already_YEEadvanced_byconformal " + std::to_string(sgg.Med[j].Is.already_YEEadvanced_byconformal);
      WarnErrReport(buff);
      buff = "Iss split_and_useless " + std::to_string(sgg.Med[j].Is.split_and_useless);
      WarnErrReport(buff);
      buff = "Is Volume " + std::to_string(sgg.Med[j].Is.Volume);
      WarnErrReport(buff);
      buff = "Is Surface " + std::to_string(sgg.Med[j].Is.Surface);
      WarnErrReport(buff);
      buff = "Is Line " + std::to_string(sgg.Med[j].Is.Line);
      WarnErrReport(buff);
      buff = "_____________________________";
      WarnErrReport(buff);
   }
   return;
}

void erasesignalingfiles(bool simu_devia) {
   std::string ficherito;
   if (!simu_devia) {
      ficherito = "stop";
      openclosedelete(ficherito);
      ficherito = "stopflushing";
      openclosedelete(ficherito);
      ficherito = "flush";
      openclosedelete(ficherito);
      ficherito = "flushdata";
      openclosedelete(ficherito);
      ficherito = "unpack";
      openclosedelete(ficherito);
      ficherito = "postprocess";
      openclosedelete(ficherito);
      ficherito = "flushxdmf";
      openclosedelete(ficherito);
      ficherito = "flushvtk";
      openclosedelete(ficherito);
      ficherito = "snap";
      openclosedelete(ficherito);
      ficherito = "stop_only";
      openclosedelete(ficherito);
      ficherito = "stopflushing_only";
      openclosedelete(ficherito);
      ficherito = "flush_only";
      openclosedelete(ficherito);
      ficherito = "flushdata_only";
      openclosedelete(ficherito);
      ficherito = "stop_dontwritevtk";
      openclosedelete(ficherito);
      ficherito = "stop_only_dontwritevtk";
      openclosedelete(ficherito);
      ficherito = "stopflushing_dontwritevtk";
      openclosedelete(ficherito);
      ficherito = "stopflushing_only_dontwritevtk";
      openclosedelete(ficherito);
      ficherito = "flush_dontwritevtk";
      openclosedelete(ficherito);
      ficherito = "flush_only_dontwritevtk";
      openclosedelete(ficherito);
      ficherito = "unpack";
      openclosedelete(ficherito);
      ficherito = "postprocess";
      openclosedelete(ficherito);
      ficherito = "flushxdmf";
      openclosedelete(ficherito);
      ficherito = "flushvtk";
      openclosedelete(ficherito);
   }
}

void openclose(const std::string& ficherin) {
    int my_iostat = 0;
    int myunit = 0;
    
retry_open2:
    if (my_iostat != 0) {
        std::cout << '.' << std::flush;
    }
    
    std::string filename = ficherin;
    size_t start = filename.find_first_not_of(' ');
    if (start == std::string::npos) {
        filename = "";
    } else {
        filename = filename.substr(start);
    }
    
    std::ofstream file(filename, std::ios::out);
    if (!file.is_open()) {
        my_iostat = 1;
        goto retry_open2;
    }
    file << "!END" << std::endl;
    file.close();
}

void opensolo(int& myunit, const std::string& ficherin) {
    int my_iostat = 0;
    
retry_open3:
    if (my_iostat != 0) {
        std::cout << '.' << std::flush;
    }
    
    std::string filename = ficherin;
    size_t start = filename.find_first_not_of(' ');
    if (start == std::string::npos) {
        filename = "";
    } else {
        filename = filename.substr(start);
    }
    
    std::ofstream file(filename, std::ios::out);
    if (!file.is_open()) {
        my_iostat = 1;
        goto retry_open3;
    }
    // In C++, we can't easily return a file descriptor/unit number like Fortran.
    // Assuming myunit is used to track the file stream or handle.
    // For simplicity, we might store the stream in a global or pass by reference differently.
    // However, since the prompt asks to translate subroutines to functions and preserve names,
    // and myunit is an output parameter in Fortran, we'll keep it as a reference.
    // Note: This is a simplification. In a real C++ translation, you'd likely use std::ofstream* or similar.
    myunit = 1; // Placeholder unit number
}

void closesolo(int& myunit) {
    int my_iostat = 0;
    
retry_close:
    if (my_iostat != 0) {
        std::cout << '.' << std::flush;
    }
    
    // In C++, we don't have a direct equivalent to closing by unit number without a map.
    // Assuming myunit was used to identify an open stream.
    // This is a placeholder for the actual closing logic which would depend on how opensolo was implemented.
    // For now, we just return as the actual file management is abstracted.
}

