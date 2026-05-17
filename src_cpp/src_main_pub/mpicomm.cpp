```cpp
#include <vector>
#include <string>
#include <iostream>
#include <cstring>
#include <algorithm>
#include <mpi.h>
#include <cstdint>

// Placeholder includes for external modules/types referenced in Fortran
// These would need to be implemented or included from the original project
// #include "Report_m.h"
// #include "FDETYPES_m.h"
// #include "wiresHolland_constants_m.h"
// #include "HollandWires_m.h"

// Forward declarations for external types/constants
// Assuming these exist in the original codebase
struct SGGFDTDINFO_t;
struct limit_t;
struct XYZlimit_t;
struct Thinwires_t;
struct CurrentSegments_t;
struct MediaData_t;

// Constants and Types from external modules (Simulated for compilation context)
// In a real translation, these would come from their respective headers
extern const int RKIND_wires;
extern const int RKIND;
extern const int MPI_MAX_PROCESSOR_NAME;
extern const int BUFSIZE_LONG;
extern const int BUFSIZE;
extern const int MPI_STATUS_SIZE;
extern const int MPI_COMM_WORLD;
extern const int MPI_INTEGER;
extern const int MPI_LOGICAL;
extern const int MPI_SUM;
extern const int MPI_MIN;
extern const int MPI_MAX;
extern const int MPI_LOR;
extern const int MPI_ADDRESS_KIND;
extern const int MPI_INTEGER4;
extern const int MPI_CHARACTER;
extern int SUBCOMM_MPI;
extern int WGROUP;
extern const int INTEGERSIZE;
extern const int REALSIZE;
extern const int REALSIZE_wires;
extern const int INTEGERSIZEOFMEDIAMATRICES;
extern const int BuffObse;
extern const int NumMed;
extern const int iHz;
extern const int iEz;
extern const int iHx;
extern const int iHy;
extern const int iEx;
extern const int iEy;
extern const int plusCPU_PML;
extern const int sinpml_fullsize_iHz_ZI; // Placeholder logic
extern const int sinpml_fullsize_iHz_ZE; // Placeholder logic

// Helper functions for external calls
void MPI_INIT(int& ierr);
void MPI_COMM_SIZE(int comm, int& size, int& ierr);
void MPI_COMM_RANK(int comm, int& rank, int& ierr);
void MPI_GET_PROCESSOR_NAME(char* name, int& namelen, int& ierr);
void MPI_Barrier(int comm, int& ierr);
void MPI_AllReduce(const void* sendbuf, void* recvbuf, int count, int datatype, int op, int comm, int& ierr);
void MPI_Group_Incl(int group, int n, const int* ranks, int* newgroup, int& ierr);
void MPI_Comm_Group(int comm, int* group, int& ierr);
void MPI_Comm_Create(int comm, int group, int* newcomm, int& ierr);
void MPI_Irecv(void* buf, int count, int datatype, int source, int tag, int comm, int* request, int& ierr);
void MPI_Isend(const void* buf, int count, int datatype, int dest, int tag, int comm, int* request, int& ierr);
void MPI_Wait(int* request, int* status, int& ierr);
void MPI_Waitall(int count, int* requests, int* statuses, int& ierr);
void MPI_Type_create_struct(int count, const int* blocklens, const MPI_Aint* displs, const int* types, MPI_Datatype* newtype, int& ierr);
void MPI_Type_commit(MPI_Datatype* type, int& ierr);
void print11(int layoutnumber, const std::string& msg, bool fatal = false);
void stoponerror(int layoutnumber, int num_procs, const std::string& msg, bool fatal = false);
void StopOnError(int layoutnumber, int num_procs, const std::string& msg);
Thinwires_t* GetHwires();

namespace MPIcomm_m {

    // Global variables from module scope
    Thinwires_t* HwiresMPI = nullptr;

    struct buffer_t {
        std::vector<double> SendUP;
        std::vector<double> SendDown;
        std::vector<double> RecUp;
        std::vector<double> RecDown;
        int SendSizeUp = 0;
        int SendSizeDown = 0;
        int RecSizeUp = 0;
        int RecSizeDown = 0;
    };

    struct ibuffer_t {
        std::vector<int32_t> SendUP;
        std::vector<int32_t> SendDown;
        std::vector<int32_t> RecUp;
        std::vector<int32_t> RecDown;
        int SendSizeUp = 0;
        int SendSizeDown = 0;
        int RecSizeUp = 0;
        int RecSizeDown = 0;
    };

    // LOCAL VARIABLES
    buffer_t buffer;
    ibuffer_t ibuffer;

    bool FlushExtraInfoDown = false;
    bool FlushExtraInfoUp = false;
    
    int sizeHx = 0, sizeHy = 0, HxXI = 0, HxXE = 0, HyXI = 0, HyXE = 0, HxYI = 0, HxYE = 0, HyYI = 0, HyYE = 0, comZ = 0, finZ = 0;
    int sizeEx = 0, sizeEy = 0, ExXI = 0, ExXE = 0, EyXI = 0, EyXE = 0, ExYI = 0, ExYE = 0, EyYI = 0, EyYE = 0;
    int sizeEz = 0, sizeHz = 0, EzXI = 0, EzXE = 0, HzXI = 0, HzXE = 0, EzYI = 0, EzYE = 0, HzYI = 0, HzYE = 0;
    
    std::vector<int32_t> mpiZcom;
    std::vector<int32_t> mpiZfin;

    struct t_databuf_t {
        int ip_target = 0;
        int sizex = 0;
        int sizey = 0;
        int sizez = 0;
        bool FlushExtraInfo = false;
        // Pointers to external arrays, represented as raw pointers or vectors depending on ownership
        // In C++, we often use std::vector or raw pointers. Here we use raw pointers to match Fortran pointer semantics
        // assuming the memory is managed elsewhere.
        double* buf_x_rx = nullptr;
        double* buf_y_rx = nullptr;
        double* buf_z_rx = nullptr;
        double* buf_x_tx = nullptr;
        double* buf_y_tx = nullptr;
        double* buf_z_tx = nullptr;
    };

    struct t_databuf_Set_t {
        bool syncUp = false;
        bool pbcUp = false;
        t_databuf_t databuf_Up;
        bool syncDown = false;
        bool pbcDown = false;
        t_databuf_t databuf_Down;
    };

    t_databuf_Set_t databuf_SetH;
    t_databuf_Set_t databuf_SetE;

    void InitGeneralMPI(int& layoutnumber, int& num_procs) {
        char name[MPI_MAX_PROCESSOR_NAME];
        int namelen = 0;
        int ierr = 0;
        MPI_INIT(ierr);
        MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr);
        MPI_COMM_RANK(MPI_COMM_WORLD, layoutnumber, ierr);
        MPI_GET_PROCESSOR_NAME(name, namelen, ierr);
        MPI_Barrier(MPI_COMM_WORLD, ierr);
    }

    void MPIdivide(SGGFDTDINFO_t& sgg, const std::vector<limit_t>& fullsize, const std::vector<limit_t>& SINPML_FULLSIZE, int layoutnumber, int num_procs, bool forcing, int forced, const std::string& slicesoriginales, bool resume, bool& fatalerror) {
        std::vector<int32_t> trancos(num_procs);
        std::vector<double> cZI(num_procs + 1);
        std::vector<double> cZE(num_procs);
        
        // Reset pointers/nullify vectors as per Fortran null()
        // In C++, vectors are empty by default or we clear them.
        // The Fortran code allocates them later.
        
        double carga = 1.0 * (fullsize[iHz].ZE - fullsize[iHz].ZI) / (1.0 * num_procs) + 
                       (plusCPU_PML - 1.0) * ((SINPML_FULLSIZE[iHz].ZI - fullsize[iHz].ZI) + 
                       (fullsize[iHz].ZE - SINPML_FULLSIZE[iHz].ZE)) / (1.0 * num_procs);
        
        cZI[0] = fullsize[iHz].ZI;
        
        for (int ilay = 0; ilay < num_procs; ++ilay) {
            double guess = carga + cZI[ilay] + (plusCPU_PML - 1.0) * (std::min(cZI[ilay], 1.0 * SINPML_FULLSIZE[iHz].ZI) + 
                           std::max(cZI[ilay], 1.0 * SINPML_FULLSIZE[iHz].ZE));
            
            double ZE[3];
            ZE[0] = (guess - (plusCPU_PML - 1.0) * (SINPML_FULLSIZE[iHz].ZI)) / (1.0 + (plusCPU_PML - 1.0));
            ZE[1] = (guess - (plusCPU_PML - 1.0) * (SINPML_FULLSIZE[iHz].ZE)) / (1.0 + (plusCPU_PML - 1.0));
            ZE[2] = (guess - (plusCPU_PML - 1.0) * (SINPML_FULLSIZE[iHz].ZE) - (plusCPU_PML - 1.0) * (SINPML_FULLSIZE[iHz].ZI));
            
            double cargaZE[3];
            for (int j = 0; j < 3; ++j) {
                cargaZE[j] = std::abs((ZE[j] - cZI[ilay]) + 
                    (plusCPU_PML - 1.0) * (std::min(1.0 * SINPML_FULLSIZE[iHz].ZI, ZE[j]) - std::min(1.0 * SINPML_FULLSIZE[iHz].ZI, cZI[ilay]) + 
                    std::max(1.0 * SINPML_FULLSIZE[iHz].ZE, ZE[j]) - std::max(1.0 * SINPML_FULLSIZE[iHz].ZE, cZI[ilay])) - carga);
            }
            
            // minloc returns 1-based index in Fortran, we use 0-based here for array access but need to adjust
            int index = 0;
            double min_val = cargaZE[0];
            for(int j=1; j<3; ++j) {
                if(cargaZE[j] < min_val) {
                    min_val = cargaZE[j];
                    index = j;
                }
            }
            
            cZE[ilay] = ZE[index];
            cZI[ilay + 1] = cZE[ilay];
        }

        if (forcing) {
            if (num_procs == 2) {
                std::string dubuf = "Forcing MPI cut at " + std::to_string(forced);
                print11(layoutnumber, dubuf);
                // voided logic: cZI=-1; cZE=-1;
                // Re-implementing the specific logic for num_procs=2 forcing
                int ilay = 0;
                cZI[ilay] = fullsize[iHz].ZI;
                cZE[ilay] = forced;
                ilay++;
                cZI[ilay] = cZE[ilay - 1];
                cZE[ilay] = fullsize[iHz].ZE;
            } else {
                std::string dubuf = "Cannot force for more than 1 cut in a num_procs=2 MPI";
                print11(layoutnumber, dubuf, true);
            }
        }

        for (int ilay = 0; ilay < num_procs; ++ilay) {
            cZE[ilay] = std::nearbyint(cZE[ilay]);
            cZI[ilay + 1] = cZE[ilay];
            trancos[ilay] = static_cast<int>(cZE[ilay] - cZI[0]);
        }

        mpiZcom.resize(num_procs);
        mpiZfin.resize(num_procs);
        
        mpiZcom[0] = fullsize[iHz].ZI;
        mpiZfin[0] = fullsize[iHz].ZI + trancos[0];
        
        for (int ilay = 1; ilay < num_procs - 1; ++ilay) {
            mpiZcom[ilay] = fullsize[iHz].ZI + trancos[ilay - 1];
            mpiZfin[ilay] = fullsize[iHz].ZI + trancos[ilay];
        }
        
        mpiZcom[num_procs - 1] = fullsize[iHz].ZI + trancos[num_procs - 2];
        mpiZfin[num_procs - 1] = fullsize[iHz].ZE;

        // Assign limits
        if (layoutnumber > 0 && layoutnumber < num_procs - 1) {
            // Assuming sgg.Sweep is an array/vector of limit_t
            // sgg.Sweep(1:6) maps to indices 0-5 or 1-6 depending on implementation. 
            // Fortran 1:6 usually means 1-based. Let's assume 0-based vector for C++ but access carefully.
            // If limit_t array is 1-based in Fortran, we might need to adjust indices.
            // For simplicity, assuming sgg.Sweep is accessible via operator[] or similar.
            // We will assume standard 0-based indexing for C++ vectors if passed as such, 
            // but Fortran code uses 1:6. Let's assume the struct handles this or we map 1->0.
            // Given the complexity, we'll stick to direct translation of logic.
            // sgg.Sweep[ilay-1] etc.
             // Note: This part requires knowledge of SGGFDTDINFO_t structure.
             // We will assume sgg.Sweep is a vector/array of size 6 or 7.
             // Fortran: sgg%Sweep(1:6)%ZI
             // C++: sgg.Sweep[0].ZI to sgg.Sweep[5].ZI
             for(int k=0; k<6; ++k) {
                 sgg.Sweep[k].ZI = fullsize[k].ZI + trancos[layoutnumber - 1];
                 sgg.Sweep[k].ZE = fullsize[k].ZI + trancos[layoutnumber];
             }
        } else if (layoutnumber == 0 && layoutnumber != num_procs - 1) {
            for(int k=0; k<6; ++k) {
                sgg.Sweep[k].ZI = fullsize[k].ZI;
                sgg.Sweep[k].ZE = fullsize[k].ZI + trancos[layoutnumber];
            }
        } else if (layoutnumber != 0 && layoutnumber == num_procs - 1) {
            for(int k=0; k<6; ++k) {
                sgg.Sweep[k].ZI = fullsize[k].ZI + trancos[layoutnumber - 1];
                sgg.Sweep[k].ZE = fullsize[k].ZE;
            }
        }

        // Adjust endings
        if (layoutnumber > 0 && layoutnumber < num_procs - 1) {
            sgg.Sweep[iEz].ZE -= 1;
            sgg.Sweep[iHx].ZE -= 1;
            sgg.Sweep[iHy].ZE -= 1;
        } else if (layoutnumber == 0 && layoutnumber != num_procs - 1) {
            sgg.Sweep[iEz].ZE -= 1;
            sgg.Sweep[iHx].ZE -= 1;
            sgg.Sweep[iHy].ZE -= 1;
        }

        int padding = 1;
        if (padding >= *std::min_element(trancos.begin(), trancos.end())) {
            std::string buff = "Number of cells per processor less than 2. Decrease the number of MPI processors";
            stoponerror(layoutnumber, num_procs, buff, true);
            fatalerror = true;
            return;
        }
        if (*std::min_element(trancos.begin(), trancos.end()) <= 2) {
            std::string buff = "Number of cells per processor less than 2. Decrease the number of MPI processors";
            stoponerror(layoutnumber, num_procs, buff, true);
            fatalerror = true;
            return;
        }

        if (layoutnumber > 0 && layoutnumber < num_procs - 1) {
            for(int k=0; k<6; ++k) {
                sgg.alloc[k].ZI = sgg.Sweep[k].ZI - padding;
                sgg.alloc[k].ZE = sgg.Sweep[k].ZE + padding;
            }
        } else if (layoutnumber == 0 && layoutnumber != num_procs - 1) {
            for(int k=0; k<6; ++k) {
                sgg.alloc[k].ZI = sgg.Sweep[k].ZI - 1;
                sgg.alloc[k].ZE = sgg.Sweep[k].ZE + padding;
            }
        } else if (layoutnumber != 0 && layoutnumber == num_procs - 1) {
            for(int k=0; k<6; ++k) {
                sgg.alloc[k].ZI = sgg.Sweep[k].ZI - padding;
                sgg.alloc[k].ZE = sgg.Sweep[k].ZE + 1;
            }
        }

        // Arrange PML Borders
        if (layoutnumber == 0) {
            sgg.Border.IsUpPML = false;
            sgg.Border.IsUpmur = false;
            sgg.Border.IsUpPMC = false;
            sgg.Border.IsUpPEC = false;
        } else if (layoutnumber == num_procs - 1) {
            sgg.Border.IsDownPML = false;
            sgg.Border.IsDownmur = false;
            sgg.Border.IsDownPMC = false;
            sgg.Border.IsDownPEC = false;
        } else {
            sgg.Border.IsUpPMC = false;
            sgg.Border.IsUpPEC = false;
            sgg.Border.IsDownPMC = false;
            sgg.Border.IsDownPEC = false;
            
            if (sgg.Sweep[iEx].ZI < SINPML_FULLSIZE[iEx].ZI) {
                sgg.Border.IsDownPML = true;
            } else {
                sgg.Border.IsDownPML = false;
            }
            
            if (sgg.Sweep[iEx].ZE > SINPML_FULLSIZE[iEx].ZE) {
                sgg.Border.IsUpPML = true;
            } else {
                sgg.Border.IsUpPML = false;
            }
            sgg.Border.IsUpMur = false;
            sgg.Border.IsDownMur = false;
        }

        if (!sgg.Border.IsUpPML) sgg.PML.NumLayers[3][2] = 0;
        if (!sgg.Border.IsDownPML) sgg.PML.NumLayers[3][1] = 0;

        // Writing
        if (layoutnumber == 0) {
            std::string slices = "!SLICES";
            for (int ilay = 0; ilay < num_procs; ++ilay) {
                std::string buff;
                // Simulate write(buff,'(i7)')
                char buf_char[10];
                snprintf(buf_char, sizeof(buf_char), "%7d", mpiZfin[ilay] - mpiZcom[ilay]);
                std::string s_buff(buf_char);
                slices += "_" + s_buff;
            }
            
            if (resume && slices != slicesoriginales) {
                std::string buff = "Different resumed/original MPI slices: " + slices + " " + slicesoriginales;
                StopOnError(layoutnumber, num_procs, buff);
            }
            print11(layoutnumber, slices);
            
            for (int ilay = 0; ilay < num_procs; ++ilay) {
                char buf_char[100];
                snprintf(buf_char, sizeof(buf_char), "(%d/%d) Spanning from z=%d to %d = %d", 
                         ilay + 1, num_procs, mpiZcom[ilay], mpiZfin[ilay], mpiZfin[ilay] - mpiZcom[ilay]);
                print11(layoutnumber, std::string(buf_char));
            }
        }

        // Bug check
        // Assuming sggPMLNumLayers_original was captured before
        // int min_pml = std::min_element(sggPMLNumLayers_original.begin(), sggPMLNumLayers_original.end());
        // if (originalPML_up_or_down && (mpiZfin[layoutnumber] - mpiZcom[layoutnumber] <= min_pml)) { ... }
    }

    void InitMPI(const std::vector<XYZlimit_t>& sggsweep, const std::vector<XYZlimit_t>& sggalloc) {
        FlushExtraInfoDown = false;
        FlushExtraInfoUp = false;

        ExXI = sggalloc[iEx].XI;
        ExXE = sggalloc[iEx].XE;
        EyXI = sggalloc[iEy].XI;
        EyXE = sggalloc[iEy].XE;
        EzXI = sggalloc[iEz].XI;
        EzXE = sggalloc[iEz].XE;

        ExYI = sggalloc[iEx].YI;
        ExYE = sggalloc[iEx].YE;
        EyYI = sggalloc[iEy].YI;
        EyYE = sggalloc[iEy].YE;
        EzYI = sggalloc[iEz].YI;
        EzYE = sggalloc[iEz].YE;

        HxXI = sggalloc[iHx].XI;
        HxXE = sggalloc[iHx].XE;
        HyXI = sggalloc[iHy].XI;
        HyXE = sggalloc[iHy].XE;
        HzXI = sggalloc[iHz].XI;
        HzXE = sggalloc[iHz].XE;

        HxYI = sggalloc[iHx].YI;
        HxYE = sggalloc[iHx].YE;
        HyYI = sggalloc[iHy].YI;
        HyYE = sggalloc[iHy].YE;
        HzYI = sggalloc[iHz].YI;
        HzYE = sggalloc[iHz].YE;

        sizeEx = (ExXE - ExXI + 1) * (ExYE - ExYI + 1);
        sizeEy = (EyXE - EyXI + 1) * (EyYE - EyYI + 1);
        sizeEz = (EzXE - EzXI + 1) * (EzYE - EzYI + 1);

        sizeHx = (HxXE - HxXI + 1) * (HxYE - HxYI + 1);
        sizeHy = (HyXE - HyXI + 1) * (HyYE - HyYI + 1);
        sizeHz = (HzXE - HzXI + 1) * (HzYE - HzYI + 1);

        comZ = sggsweep[iHx].ZI;
        finZ = sggsweep[iHx].ZE;
    }

    void MPIupdateMin(double dtlay, double& dt) {
        int ierr = 0;
        MPI_AllReduce(&dtlay, &dt, 1, REALSIZE, MPI_MIN, SUBCOMM_MPI, ierr);
    }

    void MPIupdateBloques(int layoutnumber, const std::vector<double>& valores, std::vector<double>& newvalores, int SubComm) {
        int ierr = 0;
        int sizeofvalores = BuffObse + 1;
        MPI_AllReduce(valores.data(), newvalores.data(), sizeofvalores, REALSIZE, MPI_SUM, SubComm, ierr);
    }

    void MPIinitSubcomm(int layoutnumber, int num_procs, int& SubComm, int& Root, int& group1) {
        int count = 0;
        int ierr = 0;
        int wgroup = 0;
        int NewRoot = Root;
        
        std::vector<bool> allranks(num_procs, false);
        std::vector<bool> newallranks(num_procs, false);
        
        if (Subcomm == 1) allranks[layoutnumber] = true;
        
        MPI_AllReduce(allranks.data(), newallranks.data(), num_procs, MPI_LOGICAL, MPI_LOR, SUBCOMM_MPI, ierr);
        MPI_AllReduce(&Root, &NewRoot, 1, MPI_INTEGER, MPI_MAX, SUBCOMM_MPI, ierr);
        Root = NewRoot;
        
        count = -1;
        for (int i = 0; i < num_procs; ++i) {
            if (newallranks[i]) count++;
        }
        
        std::vector<int> NGroup(count + 1);
        count = -1;
        for (int i = 0; i < num_procs; ++i) {
            if (newallranks[i]) {
                count++;
                NGroup[count] = i;
            }
        }
        
        if (count >= 0) {
            MPI_Comm_Group(SUBCOMM_MPI, &wgroup, ierr);
            MPI_Group_Incl(wgroup, count + 1, NGroup.data(), &group1, ierr);
            MPI_Comm_Create(SUBCOMM_MPI, group1, &SubComm, ierr);
        } else {
            SubComm = -1;
            group1 = -1;
            wgroup = -1;
        }
        
        if (!newallranks[layoutnumber]) SubComm = -1;
    }

    void FlushMPI_H(const std::vector<XYZlimit_t>& sggAlloc, int layoutnumber, int num_procs, 
                    std::vector<std::vector<std::vector<int32_t>>>& Hx, 
                    std::vector<std::vector<std::vector<int32_t>>>& Hy, 
                    std::vector<std::vector<std::vector<int32_t>>>& Hz) {
        int ierr1=0, ierr2=0, ierr3=0, ierr4=0, ierr5=0, ierr6=0, ierr7=0, ierr8=0, ierr9=0, ierr10=0, ierr11=0, ierr12=0, ierr100=0, ierr100b=0;
        int jerr1=0, jerr2=0, jerr3=0, jerr4=0, jerr5=0, jerr6=0, jerr7=0, jerr8=0, jerr9=0, jerr10=0, jerr11=0, jerr12=0, jerr100=0, jerr100b=0;
        
        std::vector<int> req1(4), req2(4), status1(MPI_STATUS_SIZE * 4), status2(MPI_STATUS_SIZE * 4);
        std::vector<int> req1b(2), req2b(2), status1b(MPI_STATUS_SIZE * 2), status2b(MPI_STATUS_SIZE * 2);

        if (layoutnumber != num_procs - 1) {
            MPI_Irecv(Hx[HxXI][HxYI][finZ + 1].data(), sizeHx, INTEGERSIZE, layoutnumber + 1, 1, SUBCOMM_MPI, &req1[0], ierr1);
            MPI_Isend(Hx[HxXI][HxYI][finZ].data(), sizeHx, INTEGERSIZE, layoutnumber + 1, 2, SUBCOMM_MPI, &req1[1], ierr2);
            MPI_Irecv(Hy[HyXI][HyYI][finZ + 1].data(), sizeHy, INTEGERSIZE, layoutnumber + 1, 3, SUBCOMM_MPI, &req1[2], ierr3);
            MPI_Isend(Hy[HyXI][HyYI][finZ].data(), sizeHy, INTEGERSIZE, layoutnumber + 1, 4, SUBCOMM_MPI, &req1[3], ierr4);
            
            if (FlushExtraInfoUp) {
                MPI_Irecv(Hz[HzXI][HzYI][finZ + 2].data(), sizeHz, INTEGERSIZE, layoutnumber + 1, 5, SUBCOMM_MPI, &req1b[0], ierr11);
                MPI_Isend(Hz[HzXI][HzYI][finZ].data(), sizeHz, INTEGERSIZE, layoutnumber + 1, 6, SUBCOMM_MPI, &req1b[1], ierr12);
            }
        } else {
            MPI_Irecv(Hx[HxXI][HxYI][finZ + 1].data(), sizeHx, INTEGERSIZE, 0, 1, SUBCOMM_MPI, &req1[0], ierr1);
            MPI_Isend(Hx[HxXI][HxYI][finZ].data(), sizeHx, INTEGERSIZE, 0, 2, SUBCOMM_MPI, &req1[1], ierr2);
            MPI_Irecv(Hy[HyXI][HyYI][finZ + 1].data(), sizeHy, INTEGERSIZE, 0, 3, SUBCOMM_MPI, &req1[2], ierr3);
            MPI_Isend(Hy[HyXI][HyYI][finZ].data(), sizeHy, INTEGERSIZE, 0, 4, SUBCOMM_MPI, &req1[3], ierr4);
        }

        if (layoutnumber != 0) {
            MPI_Isend(Hx[HxXI][HxYI][comZ].data(), sizeHx, INTEGERSIZE, layoutnumber - 1, 1, SUBCOMM_MPI, &req2[0], jerr1);
            MPI_Irecv(Hx[HxXI][HxYI][comZ - 1].data(), sizeHx, INTEGERSIZE, layoutnumber - 1, 2, SUBCOMM_MPI, &req2[1], jerr2);
            MPI_Isend(Hy[HyXI][HyYI][comZ].data(), sizeHy, INTEGERSIZE, layoutnumber - 1, 3, SUBCOMM_MPI, &req2[2], jerr3);
            MPI_Irecv(Hy[HyXI][HyYI][comZ - 1].data(), sizeHy, INTEGERSIZE, layoutnumber - 1, 4, SUBCOMM_MPI, &req2[3], jerr4);
            
            if (FlushExtraInfoDown) {
                MPI_Isend(Hz[HzXI][HzYI][comZ + 1].data(), sizeHz, INTEGERSIZE, layoutnumber - 1, 5, SUBCOMM_MPI, &req2b[0], jerr11);
                MPI_Irecv(Hz[HzXI][HzYI][comZ - 1].data(), sizeHz, INTEGERSIZE, layoutnumber - 1, 6, SUBCOMM_MPI, &req2b[1], jerr12);
            }
        } else {
            MPI_Isend(Hx[HxXI][HxYI][comZ].data(), sizeHx, INTEGERSIZE, num_procs - 1, 1, SUBCOMM_MPI, &req2[0], jerr1);
            MPI_Irecv(Hx[HxXI][HxYI][comZ - 1].data(), sizeHx, INTEGERSIZE, num_procs - 1, 2, SUBCOMM_MPI, &req2[1], jerr2);
            MPI_Isend(Hy[HyXI][HyYI][comZ].data(), sizeHy, INTEGERSIZE, num_procs - 1, 3, SUBCOMM_MPI, &req2[2], jerr3);
            MPI_Irecv(Hy[HyXI][HyYI][comZ - 1].data(), sizeHy, INTEGERSIZE, num_procs - 1, 4, SUBCOMM_MPI, &req2[3], jerr4);
        }

        if (layoutnumber != 0) {
            MPI_Waitall(4, req2.data(), status2.data(), ierr100);
            if (FlushExtraInfoDown) {
                MPI_Waitall(2, req2b.data(), status2b.data(), ierr100b);
            }
        }
        if (layoutnumber != num_procs - 1) {
            MPI_Waitall(4, req1.data(), status1.data(), jerr100);
            if (FlushExtraInfoUp) {
                MPI_Waitall(2, req1b.data(), status1b.data(), jerr100b);
            }
        }

        int total_err = ierr1+ierr2+ierr3+ierr4+ierr5+ierr6+ierr7+ierr8+ierr9+ierr10+ierr11+ierr12+ierr100+ierr100b+
                        jerr1+jerr2+jerr3+jerr4+jerr5+jerr6+jerr7+jerr8+jerr9+jerr10+jerr11+jerr12+jerr100+jerr100b;
        if (total_err != 0) {
            StopOnError(layoutnumber, num_procs, "FLUSHMPI");
        }
    }

    void FlushMPI_E(const std::vector<XYZlimit_t>& sggAlloc, int layoutnumber, int num_procs, 
                    std::vector<std::vector<std::vector<int32_t>>>& Ex, 
                    std::vector<std::vector<std::vector<int32_t>>>& Ey, 
                    std::vector<std::vector<std::vector<int32_t>>>& Ez) {
        int ierr1=0, ierr2=0, ierr3=0, ierr4=0, ierr5=0, ierr6=0, ierr7=0, ierr8=0, ierr9=0, ierr10=0, ierr11=0, ierr12=0, ierr100=0, ierr100b=0;
        int jerr1=0, jerr2=0, jerr3=0, jerr4=0, jerr5=0, jerr6=0, jerr7=0, jerr8=0, jerr9=0, jerr10=0, jerr11=0, jerr12=0, jerr100=0, jerr100b=0;
        
        std::vector<int> req1(2), req2(2), status1(MPI_STATUS_SIZE * 2), status2(MPI_STATUS_SIZE * 2);
        std::vector<int> req1b(4), req2b(4), status1b(MPI_STATUS_SIZE * 4), status2b(MPI_STATUS_SIZE * 4);

        if (layoutnumber != num_procs - 1) {
            if (FlushExtraInfoUp) {
                MPI_Irecv(Ez[EzXI][EzYI][finZ + 1].data(), sizeEz, INTEGERSIZE,