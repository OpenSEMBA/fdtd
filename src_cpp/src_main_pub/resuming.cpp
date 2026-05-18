#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <algorithm>

// Forward declarations for external modules/types used in the Fortran code
// These would typically be defined in their respective headers included via the 'use' statements

// Placeholder for Report_m
void print11(int unit, const std::string& msg);
const std::string SEPARADOR = "========================================";

// Placeholder for FDETYPES_m
using RKIND = double;
using RKIND_tiempo = double;
constexpr int RKIND_size = sizeof(double);

// Placeholder for types defined in other modules
struct XYZlimit_t {
    int XI, XE, YI, YE, ZI, ZE;
};

// Constants for indices, assumed to be defined in FDETYPES_m or similar
extern const int iEx, iEy, iEz, iHx, iHy, iHz;

struct SGGFDTDINFO_t {
    std::vector<double> tiempo;
    double dt;
};

struct bounds_t {
    struct { int NX, NY, NZ; } Ex, Ey, Ez, Hx, Hy, Hz;
};

struct logic_control_t {
    bool PMLBorders;
    bool PMLbodies;
    bool MURBorders;
    bool wires;
    bool Lumpeds;
    bool SGBCs;
    bool Multiports;
    bool EDispersives;
    bool MDispersives;
    bool PlaneWaveBoxes;
    bool FarFields;
};

// Placeholder for StoreFields functions from other modules
void StoreFieldsCPMLBorders();
void StorefieldsPMLbodies();
void StoreFieldsMURBorders();
void StoreFieldsWires();
void StoreFieldsWires_Berenger();
void StoreFieldsWires_Slanted();
void StoreFieldsLumpeds(bool stochastic);
void StoreFieldsSGBCs(bool stochastic);
void StoreFieldsMultiports();
void StoreFieldsEDispersives();
void StoreFieldsMDispersives();
void StorePlaneWaves(const SGGFDTDINFO_t& sgg);
void StoreFarFields(const bounds_t& b);

// Placeholder for MPI functions
#ifdef CompileWithMPI
void MPI_Barrier(int comm, int& ierr);
extern int SUBCOMM_MPI;
void newFlushWiresMPI(int layoutnumber, int num_procs);
void FlushWiresMPI_Berenger(int layoutnumber, int num_procs);
#endif

// Placeholder for Stochastic functions
#ifdef CompileWithStochastic
void syncstoch_mpi_wires(bool simu_devia, int layoutnumber, int num_procs);
void syncstoch_mpi_lumped(bool simu_devia, int layoutnumber, int num_procs);
void syncstoch_mpi_SGBCs(bool simu_devia, int layoutnumber, int num_procs);
#endif

// Constants
constexpr int BLOCK_SIZE = 1024;
constexpr int BUFSIZE = 256;

namespace resuming_m {

    // Global variables from the module
    double zvac = 0.0;
    double cluz = 0.0;
    double eps0 = 0.0;
    double mu0 = 0.0;

    void ReadFields(
        const std::vector<XYZlimit_t>& sggalloc,
        int& lastexecutedtimestep,
        double& lastexecutedtime,
        double& ultimodt,
        double eps00,
        double mu00,
        std::vector<std::vector<std::vector<double>>>& Ex,
        std::vector<std::vector<std::vector<double>>>& Ey,
        std::vector<std::vector<std::vector<double>>>& Ez,
        std::vector<std::vector<std::vector<double>>>& Hx,
        std::vector<std::vector<std::vector<double>>>& Hy,
        std::vector<std::vector<std::vector<double>>>& Hz
    ) {
        eps0 = eps00;
        mu0 = mu00;
        zvac = std::sqrt(mu0 / eps0);
        cluz = 1.0 / std::sqrt(mu0 * eps0);

        std::ifstream file(14, std::ios::in); // Assuming unit 14 maps to a specific file stream in context
        // Note: In a real translation, unit 14 would be managed by a file handling abstraction.
        // Here we assume a global or context-specific file stream for unit 14.
        // For strict translation, we might need a helper class for Fortran units.
        // Simplified: using a dummy file stream for demonstration of logic.
        // In practice, 'open(14, ...)' would be handled by a wrapper.
        
        // Since we cannot easily map Fortran unit numbers to C++ streams without context,
        // we will assume a helper function get_stream(14) exists or use a global map.
        // For this translation, we will use a placeholder file stream 'file14'.
        // The actual implementation would require a file manager.
        
        // Let's assume 'file14' is a global ifstream* or similar managed elsewhere.
        // To make this compile, we'll use a local dummy file for structure, 
        // but in reality, this needs the actual file I/O context.
        
        // Re-reading the prompt: "Convert Fortran I/O ... to C++ equivalents".
        // We will assume a helper `get_file_stream(int unit)` returns an `std::istream&`.
        
        auto& stream14 = get_file_stream(14); 

        stream14 >> lastexecutedtimestep >> lastexecutedtime >> ultimodt >> eps0 >> mu0;

        // Helper lambda to read a 1D block of data
        auto read_block = [&](std::vector<double>& data, int start, int end) {
            for (int i = start; i <= end; ++i) {
                stream14 >> data[i];
            }
        };

        // Process Ex
        for (int k = sggalloc[iEx].ZI; k <= sggalloc[iEx].ZE; ++k) {
            for (int j = sggalloc[iEx].YI; j <= sggalloc[iEx].YE; ++j) {
                int n_block = ((sggalloc[iEx].XE - sggalloc[iEx].XI + 1) / BLOCK_SIZE);
                int ini = sggalloc[iEx].XI;
                for (int i_block = 0; i_block < n_block; ++i_block) {
                    int fin = ini - 1 + BLOCK_SIZE;
                    read_block(Ex[i][j][k], ini, fin); // Note: Ex is 0-indexed in C++, but Fortran indices are used for bounds.
                    // The vector Ex passed in is likely sized to match the Fortran array dimensions.
                    // We assume Ex is accessed as Ex[i][j][k] where i,j,k are the Fortran indices.
                    // To preserve 1-based indexing logic if necessary, the vector should be sized appropriately.
                    // However, std::vector is 0-based. We will assume the caller passes vectors sized to accommodate
                    // the max index used, or we adjust indices.
                    // Given the complexity, we assume Ex is accessed directly with Fortran indices.
                    ini += BLOCK_SIZE;
                }
                read_block(Ex[i][j][k], ini, sggalloc[iEx].XE);
            }
        }

        // Process Ey
        for (int k = sggalloc[iEy].ZI; k <= sggalloc[iEy].ZE; ++k) {
            for (int j = sggalloc[iEy].YI; j <= sggalloc[iEy].YE; ++j) {
                int n_block = ((sggalloc[iEy].XE - sggalloc[iEy].XI + 1) / BLOCK_SIZE);
                int ini = sggalloc[iEy].XI;
                for (int i_block = 0; i_block < n_block; ++i_block) {
                    int fin = ini - 1 + BLOCK_SIZE;
                    read_block(Ey[i][j][k], ini, fin);
                    ini += BLOCK_SIZE;
                }
                read_block(Ey[i][j][k], ini, sggalloc[iEy].XE);
            }
        }

        // Process Ez
        for (int k = sggalloc[iEz].ZI; k <= sggalloc[iEz].ZE; ++k) {
            for (int j = sggalloc[iEz].YI; j <= sggalloc[iEz].YE; ++j) {
                int n_block = ((sggalloc[iEz].XE - sggalloc[iEz].XI + 1) / BLOCK_SIZE);
                int ini = sggalloc[iEz].XI;
                for (int i_block = 0; i_block < n_block; ++i_block) {
                    int fin = ini - 1 + BLOCK_SIZE;
                    read_block(Ez[i][j][k], ini, fin);
                    ini += BLOCK_SIZE;
                }
                read_block(Ez[i][j][k], ini, sggalloc[iEz].XE);
            }
        }

        // Process Hx
        for (int k = sggalloc[iHx].ZI; k <= sggalloc[iHx].ZE; ++k) {
            for (int j = sggalloc[iHx].YI; j <= sggalloc[iHx].YE; ++j) {
                int n_block = ((sggalloc[iHx].XE - sggalloc[iHx].XI + 1) / BLOCK_SIZE);
                int ini = sggalloc[iHx].XI;
                for (int i_block = 0; i_block < n_block; ++i_block) {
                    int fin = ini - 1 + BLOCK_SIZE;
                    read_block(Hx[i][j][k], ini, fin);
                    ini += BLOCK_SIZE;
                }
                read_block(Hx[i][j][k], ini, sggalloc[iHx].XE);
            }
        }

        // Process Hy
        for (int k = sggalloc[iHy].ZI; k <= sggalloc[iHy].ZE; ++k) {
            for (int j = sggalloc[iHy].YI; j <= sggalloc[iHy].YE; ++j) {
                int n_block = ((sggalloc[iHy].XE - sggalloc[iHy].XI + 1) / BLOCK_SIZE);
                int ini = sggalloc[iHy].XI;
                for (int i_block = 0; i_block < n_block; ++i_block) {
                    int fin = ini - 1 + BLOCK_SIZE;
                    read_block(Hy[i][j][k], ini, fin);
                    ini += BLOCK_SIZE;
                }
                read_block(Hy[i][j][k], ini, sggalloc[iHy].XE);
            }
        }

        // Process Hz
        for (int k = sggalloc[iHz].ZI; k <= sggalloc[iHz].ZE; ++k) {
            for (int j = sggalloc[iHz].YI; j <= sggalloc[iHz].YE; ++j) {
                int n_block = ((sggalloc[iHz].XE - sggalloc[iHz].XI + 1) / BLOCK_SIZE);
                int ini = sggalloc[iHz].XI;
                for (int i_block = 0; i_block < n_block; ++i_block) {
                    int fin = ini - 1 + BLOCK_SIZE;
                    read_block(Hz[i][j][k], ini, fin);
                    ini += BLOCK_SIZE;
                }
                read_block(Hz[i][j][k], ini, sggalloc[iHz].XE);
            }
        }
    }

    void flush_and_save_resume(
        const SGGFDTDINFO_t& sgg,
        const bounds_t& b,
        int layoutnumber,
        int num_procs,
        const std::string& nentradaroot,
        const std::string& nresumeable2,
        const logic_control_t& thereare,
        int fin,
        double eps00,
        double mu00,
        bool& everflushed,
        const std::vector<std::vector<std::vector<double>>>& Ex,
        const std::vector<std::vector<std::vector<double>>>& Ey,
        const std::vector<std::vector<std::vector<double>>>& Ez,
        const std::vector<std::vector<std::vector<double>>>& Hx,
        const std::vector<std::vector<std::vector<double>>>& Hy,
        const std::vector<std::vector<std::vector<double>>>& Hz,
        const std::string& wiresflavor,
        bool simu_devia,
        bool stochastic
    ) {
        eps0 = eps00;
        mu0 = mu00;
        zvac = std::sqrt(mu0 / eps0);
        cluz = 1.0 / std::sqrt(mu0 * eps0);

        char whoami[BUFSIZE];
        snprintf(whoami, BUFSIZE, "(%d/%d) ", layoutnumber + 1, num_procs);
        
        everflushed = true;

        if (layoutnumber == 0) {
            flush(11);
            flush(10);
        }

#ifdef CompileWithOldSaving
        bool existe = false;
        // inquire equivalent
        std::ifstream test_file(nresumeable2 + ".old");
        existe = test_file.good();
        test_file.close();

        if (existe) {
            int my_iostat = 0;
            // Label 8766
            // Error handling loop for open
            std::ofstream file_old(nresumeable2 + ".old", std::ios::out);
            if (!file_old) {
                // Error handling
                my_iostat = 1;
                // goto 8766 logic would be a loop or goto, simplified here
            } else {
                file_old << "!END" << std::endl;
                file_old.close();
                // close with status delete
                std::remove((nresumeable2 + ".old").c_str());
                // rename
                std::rename(nresumeable2.c_str(), (nresumeable2 + ".old").c_str());
            }
        }
#endif

#ifdef CompileWithMPI
        int ierr = 0;
        MPI_Barrier(SUBCOMM_MPI, ierr);
#endif

        int my_iostat = 0;
        // Label 8776
        std::ofstream file1(nresumeable2, std::ios::out);
        if (!file1) {
            my_iostat = 1;
        } else {
            file1 << "!END" << std::endl;
            file1.close();
            std::remove(nresumeable2.c_str());
        }

        // Label 8777
        std::ofstream file2(nresumeable2, std::ios::out | std::ios::binary);
        if (!file2) {
            my_iostat = 1;
        } else {
            StoreFields(sgg, fin, eps0, mu0, b, Ex, Ey, Ez, Hx, Hy, Hz);

            if (thereare.PMLBorders) StoreFieldsCPMLBorders();
            if (thereare.PMLbodies) StorefieldsPMLbodies();
            if (thereare.MURBorders) StoreFieldsMURBorders();

#ifdef CompileWithMPI
            if (num_procs > 1) {
                if (wiresflavor == "holland" || wiresflavor == "transition") {
                    if (num_procs > 1 && thereare.wires) {
                        newFlushWiresMPI(layoutnumber, num_procs);
                    }
#ifdef CompileWithStochastic
                    if (stochastic) {
                        syncstoch_mpi_wires(simu_devia, layoutnumber, num_procs);
                    }
#endif
                }
#ifdef CompileWithBerengerWires
                if (wiresflavor == "berenger") {
                    FlushWiresMPI_Berenger(layoutnumber, num_procs);
                }
#endif
            }
#endif

            if (thereare.wires) {
                if (wiresflavor == "holland" || wiresflavor == "transition") {
                    StoreFieldsWires();
                }
#ifdef CompileWithBerengerWires
                if (wiresflavor == "berenger") {
                    StoreFieldsWires_Berenger();
                }
#endif
#ifdef CompileWithSlantedWires
                if (wiresflavor == "slanted" || wiresflavor == "semistructured") {
                    StoreFieldsWires_Slanted();
                }
#endif
            }

#ifdef CompileWithMPI
#ifdef CompileWithStochastic
            if (stochastic) {
                syncstoch_mpi_lumped(simu_devia, layoutnumber, num_procs);
            }
#endif
#endif
            if (thereare.Lumpeds) StoreFieldsLumpeds(stochastic);

#ifdef CompileWithMPI
#ifdef CompileWithStochastic
            if (stochastic) {
                syncstoch_mpi_SGBCs(simu_devia, layoutnumber, num_procs);
            }
#endif
#endif
            if (thereare.SGBCs) StoreFieldsSGBCs(stochastic);

#ifdef CompileWithNIBC
            if (thereare.Multiports) StoreFieldsMultiports();
#endif
            if (thereare.EDispersives) StoreFieldsEDispersives();
            if (thereare.MDispersives) StoreFieldsMDispersives();
            if (thereare.PlaneWaveBoxes) StorePlaneWaves(sgg);
            if (thereare.FarFields) StoreFarFields(b);

#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
            file2.close();
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
#ifdef CompileWithMPI
            MPI_Barrier(SUBCOMM_MPI, ierr);
#endif
            goto label_635;
        }

label_634:
        print11(0, SEPARADOR + SEPARADOR + SEPARADOR);
        print11(0, "RESUMING FLUSHSAVEANDRESUME: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
        print11(0, SEPARADOR + SEPARADOR + SEPARADOR);

label_635:
        return;
    }

    void StoreFields(
        const SGGFDTDINFO_t& sgg,
        int finaltimestep,
        double eps0,
        double mu0,
        const bounds_t& b,
        const std::vector<std::vector<std::vector<double>>>& Ex,
        const std::vector<std::vector<std::vector<double>>>& Ey,
        const std::vector<std::vector<std::vector<double>>>& Ez,
        const std::vector<std::vector<std::vector<double>>>& Hx,
        const std::vector<std::vector<std::vector<double>>>& Hy,
        const std::vector<std::vector<std::vector<double>>>& Hz
    ) {
        auto& stream14 = get_file_stream(14);
        stream14 << finaltimestep << " " << sgg.tiempo[finaltimestep] << " " << sgg.dt << " " << eps0 << " " << mu0 << std::endl;

        auto write_block = [&](const std::vector<std::vector<std::vector<double>>>& data, int NX, int NY, int NZ) {
            for (int k = 0; k < NZ; ++k) {
                for (int j = 0; j < NY; ++j) {
                    int n_block = (NX / BLOCK_SIZE);
                    int ini = 0;
                    for (int i_block = 0; i_block < n_block; ++i_block) {
                        int fin = ini - 1 + BLOCK_SIZE;
                        for (int i = ini; i <= fin; ++i) {
                            stream14 << data[i][j][k] << " ";
                        }
                        stream14 << std::endl;
                        ini += BLOCK_SIZE;
                    }
                    for (int i = ini; i < NX; ++i) {
                        stream14 << data[i][j][k] << " ";
                    }
                    stream14 << std::endl;
                }
            }
        };

        write_block(Ex, b.Ex.NX, b.Ex.NY, b.Ex.NZ);
        write_block(Ey, b.Ey.NX, b.Ey.NY, b.Ey.NZ);
        write_block(Ez, b.Ez.NX, b.Ez.NY, b.Ez.NZ);
        write_block(Hx, b.Hx.NX, b.Hx.NY, b.Hx.NZ);
        write_block(Hy, b.Hy.NX, b.Hy.NY, b.Hy.NZ);
        write_block(Hz, b.Hz.NX, b.Hz.NY, b.Hz.NZ);

        goto label_635;

label_634:
        print11(0, SEPARADOR + SEPARADOR + SEPARADOR);
        print11(0, "RESUMING STOREFIELDS: ERROR WRITING RESTARTING FIELDS. IGNORING AND CONTINUING");
        print11(0, SEPARADOR + SEPARADOR + SEPARADOR);

label_635:
        return;
    }

} // namespace resuming_m