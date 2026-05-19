#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <complex>

// Forward declarations and includes for external modules/types
// #include "FDETYPES_m.h"
// #include "Observa_m.h"
// #include "Report_m.h"

// Assuming these types/constants are defined in included headers
// extern int IKINDMTAG;
// extern int RKIND;
// extern int RKIND_tiempo;
// extern int REALSIZE;
// extern int MPI_INTEGER;
// extern int MPI_SUM;
// extern int iCur, iCurX, iCurY, iCurZ, mapvtk;
constexpr int iCur = 53, iCurX = 54, iCurY = 55, iCurZ = 56, mapvtk = 57;
// extern int BUFSIZE;

// Placeholder for external functions/types if not fully defined here
// In a real translation, these would be in their respective .h files
struct SGGFDTDINFO_t {
    int NumberRequest;
    struct {
        bool Volumic;
        int nP;
        struct {
            int What[10]; // Assuming size, needs actual definition
            int XI, XE, YI, YE, ZI, ZE;
        } P[10]; // Assuming size
        bool done;
        bool flushed;
        bool Begun;
        bool TimeDomain;
        bool FreqDomain;
    } observation[100]; // Assuming size
    struct {
        int XI, XE, YI, YE, ZI, ZE;
    } Alloc[10]; // Assuming size
};

struct output_t {
    struct {
        struct {
            std::string path;
            int UNIT;
            int unit; // Note: Fortran uses lowercase often for variables, keeping name
            int columnas;
            int MPISubComm;
            int ZIorig, ZEorig;
            int MPIRoot;
        } item[10]; // Assuming size
        int TimesWritten;
    } item[100]; // Assuming size
};

struct Serialized_t {
    std::vector<int> eI;
    std::vector<int> eJ;
    std::vector<int> eK;
    std::vector<int> currentType;
    std::vector<int> sggMtag;
    std::vector<double> Valor;
    std::vector<double> Valor_x;
    std::vector<double> Valor_y;
    std::vector<double> Valor_z;
    std::vector<double> ValorE;
    std::vector<double> Valor_Ex;
    std::vector<double> Valor_Ey;
    std::vector<double> Valor_Ez;
    std::vector<double> ValorH;
    std::vector<double> Valor_Hx;
    std::vector<double> Valor_Hy;
    std::vector<double> Valor_Hz;
    std::vector<std::complex<double>> ValorComplex_x;
    std::vector<std::complex<double>> ValorComplex_y;
    std::vector<std::complex<double>> ValorComplex_z;
    std::vector<std::complex<double>> ValorComplex_Ex;
    std::vector<std::complex<double>> ValorComplex_Ey;
    std::vector<std::complex<double>> ValorComplex_Ez;
    std::vector<std::complex<double>> ValorComplex_Hx;
    std::vector<std::complex<double>> ValorComplex_Hy;
    std::vector<std::complex<double>> ValorComplex_Hz;

    void allocate_for_time_domain(int n) {
        eI.resize(n); eJ.resize(n); eK.resize(n);
        currentType.resize(n); sggMtag.resize(n);
        Valor.resize(n); Valor_x.resize(n); Valor_y.resize(n); Valor_z.resize(n);
        ValorE.resize(n); Valor_Ex.resize(n); Valor_Ey.resize(n); Valor_Ez.resize(n);
        ValorH.resize(n); Valor_Hx.resize(n); Valor_Hy.resize(n); Valor_Hz.resize(n);
    }

    void allocate_for_frequency_domain(int n) {
        eI.resize(n); eJ.resize(n); eK.resize(n);
        currentType.resize(n); sggMtag.resize(n);
        ValorComplex_x.resize(n); ValorComplex_y.resize(n); ValorComplex_z.resize(n);
        ValorComplex_Ex.resize(n); ValorComplex_Ey.resize(n); ValorComplex_Ez.resize(n);
        ValorComplex_Hx.resize(n); ValorComplex_Hy.resize(n); ValorComplex_Hz.resize(n);
    }

    void allocate_current_value(int n) {
        // This method seems to be a placeholder or specific to the original logic
        // In C++, vectors are resized in the allocate_for_* methods above
    }
};

// External functions placeholders
output_t* GetOutput();
void stoponerror(int layoutnumber, int num_procs, const std::string& msg);
void print11(int layoutnumber, const std::string& buff);
void creaUnstructData(Serialized_t& Serialized, int numberOfSerialized, SGGFDTDINFO_t& sgg, 
                      std::vector<double>& Nodes, int& NumNodes, 
                      std::vector<int>& Elems, int& NumEdges, int& NumQuads, bool vtkindex);

// MPI placeholders
extern "C" {
    void MPI_Barrier(int comm, int& ierr);
    void MPI_AllReduce(const void* sendbuf, void* recvbuf, int count, int datatype, int op, int comm, int& ierr);
}

namespace VTK_m {

    void createVTK(int layoutnumber, int num_procs, SGGFDTDINFO_t& sgg, 
                   const std::vector<std::vector<int>>& sggMtag, 
                   bool& somethingdone, int mpidir, bool dontwritevtk) {
        
        // Note: sggMtag in Fortran is a 3D array. In C++, we pass it as a flattened vector or a custom struct.
        // For simplicity, assuming a helper to access it or converting the interface.
        // The original signature: integer(kind=IKINDMTAG) :: sggMtag (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, ...)
        // We will assume a helper function or struct to access this data if needed, 
        // but for this translation, we'll keep the variable name and assume access logic is handled or passed differently.
        // To strictly preserve names, we keep 'sggMtag' but its type needs adaptation.
        // Let's assume a 1D vector representing the 3D array for simplicity in this snippet, 
        // or better, pass a reference to a data structure. 
        // Given the complexity, we'll stick to the variable name and assume the caller handles the mapping.
        
        bool yacreado = false;
        std::string filename;
        std::string fichero;
        std::string fichero_input;
        std::string char_i_sub_time;
        int k;
        std::vector<std::string> suffFile = {"_current.vtk", "_efield.vtk ", "_hfield.vtk "};
        std::vector<std::string> suffTag = {"cu", "ef", "hf"};
        
        int ierr;
        int posicionMPI;
        int conta;
        int ecurrentType;
        int eei, eej, eek, esggMtag;
        std::vector<int> sizeofvalores;
        std::vector<int> NewsizeOfValores;
        
        double time;
        double rdum;
        
        output_t* output = GetOutput();
        int iroot;
        
        Serialized_t NewSerialized;
        Serialized_t Serialized;
        std::vector<int> PosiMPI;
        std::vector<int> NewPosiMPI;
        int indi;
        int numberOfSerialized;
        std::vector<double> att;
        double att_rkind;
        double att_rkind_tiempo;
        
        int ii, i1, finalstep;
        bool lexis;
        bool freqdomain;
        std::string dubuf;
        int minXabs, maxXabs, minYabs, maxYabs, minZabs, maxZabs;
        std::string pathroot;
        std::string chari, charj, chark, chari2, charj2, chark2;
        std::string extpoint;
        std::string buff;
        std::string charc;
        std::string tag;
        std::string whoami, whoamishort;
        int numNodes, numEdges, numQuads, iroot2, iroot1, i_sub_time, total_sub_times;
        const int time_phases_param = 35;
        std::vector<double> Nodes;
        std::vector<int> Elems;
        int coldummy;
        std::vector<int> volumicCurrentFlags = {iCur, iCurX, iCurY, iCurZ, mapvtk};
        
        yacreado = false;
        numNodes = 0;
        numEdges = 0;
        numQuads = 0;

        // Format whoamishort: '(i5)'
        char buffer[10];
        snprintf(buffer, sizeof(buffer), "%5d", layoutnumber + 1);
        whoamishort = std::string(buffer);
        
        // Format whoami: '(a,i5,a,i5,a)'
        snprintf(buffer, sizeof(buffer), "(%d/%d) ", layoutnumber + 1, num_procs);
        whoami = std::string(buffer);

        output = GetOutput(); // Re-assign if GetOutput returns a pointer to a singleton or similar

        somethingdone = false;

        for (ii = 0; ii < sgg.NumberRequest; ++ii) { // Fortran 1-based, assuming ii starts at 1 in original, but loop is 1 to NumberRequest
            // Adjust for 0-based index if sgg.observation is 0-based in C++ struct
            // Original: do ii = 1, sgg%NumberRequest
            // Accessing sgg.observation[ii-1] if 0-based, or sgg.observation[ii] if 1-based padding.
            // Let's assume the struct members are accessed with 1-based indexing logic mapped to 0-based vector/array.
            // For strict translation, we'll use ii+1 if the struct is 1-based, but C++ structs are 0-based.
            // We will assume sgg.observation is a vector of size NumberRequest+1 for 1-based access.
            
            if (ii >= sgg.NumberRequest) break; // Safety check

            if (sgg.observation[ii].Volumic && sgg.observation[ii].nP == 1) {
                bool found = false;
                for (int flag : volumicCurrentFlags) {
                    if (sgg.observation[ii].P[0].What[0] == flag) { // Assuming What is an array, accessing first element or checking all
                        // Original: any(sgg%observation(ii)%P(1)%What == volumicCurrentFlags)
                        // This checks if any element in What matches any in volumicCurrentFlags
                        // Simplified: check if the first element of What is in the flags list
                        // A proper 'any' check would iterate through What and volumicCurrentFlags
                        // For now, assuming a simple match
                        found = true;
                        break;
                    }
                }
                // More accurate 'any' check:
                found = false;
                for (int w : sgg.observation[ii].P[0].What) {
                    for (int v : volumicCurrentFlags) {
                        if (w == v) {
                            found = true;
                            break;
                        }
                    }
                    if (found) break;
                }

                if (found) {
                    if (sgg.observation[ii].done) {
                        if (sgg.observation[ii].flushed) {
                            continue;
                        } else {
                            sgg.observation[ii].flushed = true;
                            continue;
                        }
                    } else {
                        if (sgg.observation[ii].Begun) {
                            continue;
                        } else {
                            continue;
                        }
                    }
                } else {
                    continue;
                }
            } else {
                continue;
            }

            // sondas Volumic traducelas a VTK
            if (sgg.observation[ii].Volumic) {
                if (sgg.observation[ii].nP == 1) {
                    bool found2 = false;
                    for (int w : sgg.observation[ii].P[0].What) {
                        for (int v : volumicCurrentFlags) {
                            if (w == v) {
                                found2 = true;
                                break;
                            }
                        }
                        if (found2) break;
                    }
                    if (found2) {
                        // inquire(FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                        // C++ doesn't have direct inquire, use file existence check
                        std::ifstream test_file(output->item[ii].item[0].path);
                        lexis = test_file.good();
                        test_file.close();

                        if (lexis && output->item[ii].TimesWritten != 0) {
                            minXabs = sgg.observation[ii].P[0].XI;
                            maxXabs = sgg.observation[ii].P[0].XE;
                            minYabs = sgg.observation[ii].P[0].YI;
                            maxYabs = sgg.observation[ii].P[0].YE;
                            
                            // #ifdef CompileWithMPI
                            minZabs = output->item[ii].item[0].ZIorig;
                            maxZabs = output->item[ii].item[0].ZEorig;
                            // #else
                            // minZabs = sgg.observation[ii].P[0].zI;
                            // maxZabs = sgg.observation[ii].P[0].zE;
                            // #endif

                            snprintf(buffer, sizeof(buffer), "%7d", minXabs);
                            chari = std::string(buffer);
                            snprintf(buffer, sizeof(buffer), "%7d", minYabs);
                            charj = std::string(buffer);
                            snprintf(buffer, sizeof(buffer), "%7d", minZabs);
                            chark = std::string(buffer);
                            snprintf(buffer, sizeof(buffer), "%7d", maxXabs);
                            chari2 = std::string(buffer);
                            snprintf(buffer, sizeof(buffer), "%7d", maxYabs);
                            charj2 = std::string(buffer);
                            snprintf(buffer, sizeof(buffer), "%7d", maxZabs);
                            chark2 = std::string(buffer);

                            // mpidir logic
                            if (mpidir == 3) {
                                extpoint = chari + "_" + charj + "_" + chark + "__" + chari2 + "_" + charj2 + "_" + chark2;
                            } else if (mpidir == 2) {
                                extpoint = charj + "_" + chark + "_" + chari + "__" + charj2 + "_" + chark2 + "_" + chari2;
                            } else if (mpidir == 1) {
                                extpoint = chark + "_" + chari + "_" + charj + "__" + chark2 + "_" + chari2 + "_" + charj2;
                            } else {
                                stoponerror(layoutnumber, num_procs, "Buggy error in mpidir. ");
                            }

                            // iroot=index(output(ii)%item(1)%path,'__',.true.)
                            size_t iroot_pos = output->item[ii].item[0].path.rfind("__");
                            if (iroot_pos != std::string::npos) {
                                pathroot = output->item[ii].item[0].path.substr(0, iroot_pos);
                            } else {
                                pathroot = output->item[ii].item[0].path;
                            }
                            
                            // iroot = index (pathroot, '_',.true.)
                            iroot_pos = pathroot.rfind('_');
                            if (iroot_pos != std::string::npos) {
                                pathroot = pathroot.substr(0, iroot_pos);
                            }
                            
                            iroot_pos = pathroot.rfind('_');
                            if (iroot_pos != std::string::npos) {
                                pathroot = pathroot.substr(0, iroot_pos);
                            }
                            
                            iroot_pos = pathroot.rfind('_');
                            if (iroot_pos != std::string::npos) {
                                pathroot = pathroot.substr(0, iroot_pos);
                            }
                            
                            pathroot += "_" + extpoint;
                            filename = pathroot;

#ifdef CompileWithMPI
                            if (num_procs > 1) {
                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    continue;
                                }
                            }
#endif
                            finalstep = output->item[ii].TimesWritten;
                            att.resize(finalstep + 1); // 1-based indexing in Fortran, so size is finalstep+1, using 1..finalstep

                            numberOfSerialized = 0;
                            sizeofvalores.resize(num_procs);
                            std::fill(sizeofvalores.begin(), sizeofvalores.end(), 0);
                            sizeofvalores[layoutnumber] = output->item[ii].item[0].columnas;

#ifdef CompileWithMPI
                            if (num_procs > 1) {
                                NewsizeOfValores.resize(num_procs);
                                std::fill(NewsizeOfValores.begin(), NewsizeOfValores.end(), 0);
                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(sizeofvalores.data(), NewsizeOfValores.data(), num_procs, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                sizeofvalores = NewsizeOfValores;
                            }
#endif

                            for (i1 = 0; i1 < num_procs; ++i1) {
                                numberOfSerialized += sizeofvalores[i1];
                            }

                            PosiMPI.resize(numberOfSerialized + 1); // 1-based

                            if (sgg.observation[ii].TimeDomain) {
                                Serialized.allocate_for_time_domain(numberOfSerialized);
                                freqdomain = false;
                            } else if (sgg.observation[ii].FreqDomain) {
                                Serialized.allocate_for_frequency_domain(numberOfSerialized);
                                freqdomain = true;
                            }
                            Serialized.allocate_current_value(numberOfSerialized);
                            std::fill(PosiMPI.begin(), PosiMPI.end(), 0);

                            posicionMPI = 0;
#ifdef CompileWithMPI
                            if (num_procs > 1) {
                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    for (i1 = 0; i1 < layoutnumber; ++i1) {
                                        posicionMPI += sizeofvalores[i1];
                                    }
                                }
                                NewPosiMPI.resize(numberOfSerialized + 1);
                                if (sgg.observation[ii].TimeDomain) {
                                    NewSerialized.allocate_for_time_domain(numberOfSerialized);
                                } else if (sgg.observation[ii].FreqDomain) {
                                    NewSerialized.allocate_for_frequency_domain(numberOfSerialized);
                                }
                                NewSerialized.allocate_current_value(numberOfSerialized);
                                std::fill(NewPosiMPI.begin(), NewPosiMPI.end(), 0);
                            }
#endif

                            // open(output(ii)%item(1)%UNIT, FILE=trim(adjustl(output(ii)%item(1)%path)), FORM='unformatted')
                            // C++ fstream doesn't support unformatted binary read/write easily in the same way.
                            // Assuming standard binary read.
                            std::ifstream file(output->item[ii].item[0].path, std::ios::binary);
                            if (!file.is_open()) {
                                // Handle error
                            }
                            
                            file.read(reinterpret_cast<char*>(&coldummy), sizeof(int));
                            if (coldummy != output->item[ii].item[0].columnas) {
                                snprintf(buffer, sizeof(buffer), "ERROR: Buggy error creating .vtk%9d%9d", coldummy, output->item[ii].item[0].columnas);
                                print11(0, std::string(buffer));
                            }

                            for (conta = 1; conta <= output->item[ii].item[0].columnas; ++conta) {
                                file.read(reinterpret_cast<char*>(&eei), sizeof(int));
                                file.read(reinterpret_cast<char*>(&eej), sizeof(int));
                                file.read(reinterpret_cast<char*>(&eek), sizeof(int));
                                file.read(reinterpret_cast<char*>(&ecurrentType), sizeof(int));
                                file.read(reinterpret_cast<char*>(&esggMtag), sizeof(int));
                                
                                PosiMPI[posicionMPI + conta] = posicionMPI + conta;
                                Serialized.eI[posicionMPI + conta - 1] = eei; // Adjust for 0-based vector
                                Serialized.eJ[posicionMPI + conta - 1] = eej;
                                Serialized.eK[posicionMPI + conta - 1] = eek;
                                Serialized.currentType[posicionMPI + conta - 1] = ecurrentType;
                                Serialized.sggMtag[posicionMPI + conta - 1] = esggMtag;
                            }

                            if (sgg.observation[ii].FreqDomain) {
                                file.read(reinterpret_cast<char*>(&rdum), sizeof(double));
                            }

#ifdef CompileWithMPI
                            if (num_procs > 1) {
                                std::fill(NewPosiMPI.begin(), NewPosiMPI.end(), -1);
                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(PosiMPI.data(), NewPosiMPI.data(), numberOfSerialized, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                PosiMPI = NewPosiMPI;

                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(Serialized.eI.data(), NewSerialized.eI.data(), numberOfSerialized, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                Serialized.eI = NewSerialized.eI;

                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(Serialized.eJ.data(), NewSerialized.eJ.data(), numberOfSerialized, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                Serialized.eJ = NewSerialized.eJ;

                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(Serialized.eK.data(), NewSerialized.eK.data(), numberOfSerialized, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                Serialized.eK = NewSerialized.eK;

                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(Serialized.currentType.data(), NewSerialized.currentType.data(), numberOfSerialized, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                Serialized.currentType = NewSerialized.currentType;

                                if (output->item[ii].item[0].MPISubComm != -1) {
                                    int err;
                                    MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                    MPI_AllReduce(Serialized.sggMtag.data(), NewSerialized.sggMtag.data(), numberOfSerialized, MPI_INTEGER, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                }
                                Serialized.sggMtag = NewSerialized.sggMtag;
                            }
#endif

#ifdef CompileWithMPI
                            if (layoutnumber == output->item[ii].item[0].MPIRoot) {
#else
                            if (layoutnumber == 0) {
#endif
                                creaUnstructData(Serialized, numberOfSerialized, sgg, Nodes, numNodes, Elems, numEdges, numQuads, false); // vtkindex is input, but not used in call in snippet
                            }

                            for (indi = 1; indi <= finalstep; ++indi) {
                                if (sgg.observation[ii].TimeDomain) {
                                    std::fill(Serialized.Valor.begin(), Serialized.Valor.end(), 0.0);
                                    std::fill(Serialized.Valor_x.begin(), Serialized.Valor_x.end(), 0.0);
                                    std::fill(Serialized.Valor_y.begin(), Serialized.Valor_y.end(), 0.0);
                                    std::fill(Serialized.Valor_z.begin(), Serialized.Valor_z.end(), 0.0);
                                    std::fill(Serialized.ValorE.begin(), Serialized.ValorE.end(), 0.0);
                                    std::fill(Serialized.Valor_Ex.begin(), Serialized.Valor_Ex.end(), 0.0);
                                    std::fill(Serialized.Valor_Ey.begin(), Serialized.Valor_Ey.end(), 0.0);
                                    std::fill(Serialized.Valor_Ez.begin(), Serialized.Valor_Ez.end(), 0.0);
                                    std::fill(Serialized.ValorH.begin(), Serialized.ValorH.end(), 0.0);
                                    std::fill(Serialized.Valor_Hx.begin(), Serialized.Valor_Hx.end(), 0.0);
                                    std::fill(Serialized.Valor_Hy.begin(), Serialized.Valor_Hy.end(), 0.0);
                                    std::fill(Serialized.Valor_Hz.begin(), Serialized.Valor_Hz.end(), 0.0);

                                    file.read(reinterpret_cast<char*>(&att_rkind_tiempo), sizeof(double));
                                    att[indi] = att_rkind_tiempo;

                                    if (output->item[ii].item[0].columnas != 0) {
                                        for (conta = 1; conta <= output->item[ii].item[0].columnas; ++conta) {
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_x[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_y[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_z[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.ValorE[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_Ex[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_Ey[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_Ez[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.ValorH[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_Hx[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_Hy[posicionMPI + conta - 1]), sizeof(double));
                                            file.read(reinterpret_cast<char*>(&Serialized.Valor_Hz[posicionMPI + conta - 1]), sizeof(double));
                                        }
                                    }
                                } else if (sgg.observation[ii].FreqDomain) {
                                    std::fill(Serialized.ValorComplex_x.begin(), Serialized.ValorComplex_x.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_y.begin(), Serialized.ValorComplex_y.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_z.begin(), Serialized.ValorComplex_z.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_Ex.begin(), Serialized.ValorComplex_Ex.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_Ey.begin(), Serialized.ValorComplex_Ey.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_Ez.begin(), Serialized.ValorComplex_Ez.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_Hx.begin(), Serialized.ValorComplex_Hx.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_Hy.begin(), Serialized.ValorComplex_Hy.end(), std::complex<double>(0,0));
                                    std::fill(Serialized.ValorComplex_Hz.begin(), Serialized.ValorComplex_Hz.end(), std::complex<double>(0,0));

                                    file.read(reinterpret_cast<char*>(&att_rkind), sizeof(double));
                                    att[indi] = att_rkind;

                                    if (output->item[ii].item[0].columnas != 0) {
                                        for (conta = 1; conta <= output->item[ii].item[0].columnas; ++conta) {
                                            file.read(reinterpret_cast<char*>(&Serialized.ValorComplex_x[posicionMPI + conta - 1]), sizeof(std::complex<double>));
                                            file.read(reinterpret_cast<char*>(&Serialized.ValorComplex_y[posicionMPI + conta - 1]), sizeof(std::complex<double>));
                                            file.read(reinterpret_cast<char*>(&Serialized.ValorComplex_z[posicionMPI + conta - 1]), sizeof(std::complex<double>));
                                        }
                                    }
                                }

#ifdef CompileWithMPI
                                int err;
                                MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                if (num_procs > 1) {
                                    if (output->item[ii].item[0].MPISubComm != -1) {
                                        if (sgg.observation[ii].TimeDomain) {
                                            MPI_AllReduce(Serialized.Valor.data(), NewSerialized.Valor.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor = NewSerialized.Valor;

                                            MPI_AllReduce(Serialized.Valor_x.data(), NewSerialized.Valor_x.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_x = NewSerialized.Valor_x;

                                            MPI_AllReduce(Serialized.Valor_y.data(), NewSerialized.Valor_y.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_y = NewSerialized.Valor_y;

                                            MPI_AllReduce(Serialized.Valor_z.data(), NewSerialized.Valor_z.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_z = NewSerialized.Valor_z;

                                            MPI_AllReduce(Serialized.ValorE.data(), NewSerialized.ValorE.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.ValorE = NewSerialized.ValorE;

                                            MPI_AllReduce(Serialized.Valor_Ex.data(), NewSerialized.Valor_Ex.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_Ex = NewSerialized.Valor_Ex;

                                            MPI_AllReduce(Serialized.Valor_Ey.data(), NewSerialized.Valor_Ey.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_Ey = NewSerialized.Valor_Ey;

                                            MPI_AllReduce(Serialized.Valor_Ez.data(), NewSerialized.Valor_Ez.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_Ez = NewSerialized.Valor_Ez;

                                            MPI_AllReduce(Serialized.ValorH.data(), NewSerialized.ValorH.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.ValorH = NewSerialized.ValorH;

                                            MPI_AllReduce(Serialized.Valor_Hx.data(), NewSerialized.Valor_Hx.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_Hx = NewSerialized.Valor_Hx;

                                            MPI_AllReduce(Serialized.Valor_Hy.data(), NewSerialized.Valor_Hy.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_Hy = NewSerialized.Valor_Hy;

                                            MPI_AllReduce(Serialized.Valor_Hz.data(), NewSerialized.Valor_Hz.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            MPI_Barrier(output->item[ii].item[0].MPISubComm, err);
                                            Serialized.Valor_Hz = NewSerialized.Valor_Hz;
                                        } else if (sgg.observation[ii].FreqDomain) {
                                            std::fill(Serialized.Valor_x.begin(), Serialized.Valor_x.end(), 0.0);
                                            std::fill(NewSerialized.Valor_x.begin(), NewSerialized.Valor_x.end(), 0.0);
                                            for (conta = 1; conta <= output->item[ii].item[0].columnas; ++conta) {
                                                Serialized.Valor_x[posicionMPI + conta - 1] = std::real(Serialized.ValorComplex_x[posicionMPI + conta - 1]);
                                            }
                                            MPI_AllReduce(Serialized.Valor_x.data(), NewSerialized.Valor_x.data(), numberOfSerialized, REALSIZE, MPI_SUM, output->item[ii].item[0].MPISubComm, err);
                                            // The rest of the frequency domain MPI reduction is cut off in the source
                                        }
                                    }
                                }
#endif
                            }
                        }
                    }
                }
            }
        }
    }

    void createVTKOnTheFly() {
        // Stub for the second public subroutine
    }

} // namespace VTK_m
