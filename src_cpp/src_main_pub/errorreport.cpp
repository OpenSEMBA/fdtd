#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <sstream>

// Forward declarations for external types/functions not defined in this snippet
// These would typically come from FDETYPES_m, snapxdmf_m, and MPI headers

// Placeholder for MPI types and functions if CompileWithMPI is defined
#ifdef CompileWithMPI
#include <mpi.h>
extern MPI_Comm SUBCOMM_MPI;
extern int REALSIZE; // Assuming this is defined in FDETYPES_m
#endif

// Placeholder for external types from FDETYPES_m
struct SGGFDTDINFO_t;
struct sim_control_t;
struct logic_control_t;
struct bounds_t;
struct perform_t;
struct coorsxyzP_t;

// Placeholder for external constants/functions from FDETYPES_m
extern int BUFSIZE;
extern int RKIND; // Assuming RKIND maps to a C++ type size or is used for casting
extern int iEx, iEy, iEz, iHx, iHy, iHz;
extern double dt0;
extern double topCPUtime;
extern std::string SEPARADOR;
extern std::string separador;
extern std::string BuffObse;
extern int INITIALtimeSTEP;

// Placeholder for external functions from snapxdmf_m
void write_xdmfsnap(int n, const std::string& fichsnap, int ini_ibox, int fin_ibox, int ini_jbox, int fin_jbox, int ini_kbox, int fin_kbox, const std::vector<std::vector<std::vector<std::vector<double>>>>& snap);

// Helper function to simulate Fortran's trim(adjustl(...))
std::string trim_adjustl(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r\f\v");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r\f\v");
    return str.substr(start, end - start + 1);
}

// Helper function to simulate Fortran's write format for whoami
std::string format_whoami(int layoutnumber, int num_procs) {
    std::ostringstream oss;
    oss << "(" << (layoutnumber + 1) << "/" << num_procs << ") ";
    return oss.str();
}

// Helper function to simulate Fortran's write format for whoami short
std::string format_whoami_short(int layoutnumber) {
    std::ostringstream oss;
    oss << (layoutnumber + 1);
    return oss.str();
}

// Helper function to simulate Fortran's write format for time
std::string format_time(const std::string& fecha, const std::string& hora) {
    // fecha: YYYYMMDD, hora: HHMMSS
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices
std::string format_time_part(const std::string& fecha, const std::string& hora) {
    // This is a simplified version. Fortran string indexing is 1-based.
    // fecha(7:8) -> index 6:8 in C++ (0-based)
    // fecha(5:6) -> index 4:6 in C++
    // fecha(1:4) -> index 0:4 in C++
    // hora(1:2) -> index 0:2 in C++
    // hora(3:4) -> index 2:4 in C++
    // hora(5:6) -> index 4:6 in C++
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (InitTiming)
std::string format_time_init(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing2(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing3(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing4(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing5(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing6(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing7(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing8(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing9(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing10(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing11(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing12(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing13(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing14(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing15(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing16(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing17(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing18(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing19(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing20(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing21(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing22(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing23(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing24(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing25(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing26(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing27(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing28(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing29(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing30(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing31(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing32(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing33(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing34(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing35(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing36(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing37(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing38(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing39(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing40(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing41(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing42(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing43(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing44(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing45(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing46(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing47(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing48(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing49(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time with specific indices (Timing)
std::string format_time_timing50(const std::string& fecha, const std::string& hora) {
    if (fecha.length() >= 8 && hora.length() >= 6) {
        return fecha.substr(6, 2) + "/" + fecha.substr(4, 2) + "/" + fecha.substr(0, 4) + "   " +
               hora.substr(0, 2) + ":" + hora.substr(2, 2) + ":" + hora.substr(4, 2);
    }
    return "";
}

// Helper function to simulate Fortran's write format for time