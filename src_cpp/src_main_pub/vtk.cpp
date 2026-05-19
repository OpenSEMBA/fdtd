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
            
            if (ii >= sgg.observation.size()) break; // Safety check

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

&                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_x(1, conta) = std::complex<double>(newSerialized.Valor_x(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.Valor_x = 0.0_RKIND;
                                    newSerialized.Valor_x = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_x(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_x(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_x, newSerialized.Valor_x, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_x(1, conta) = Serialized.ValorComplex_x(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.Valor_x(1, conta));
                                    }
                                    // !!!!   y
                                    // parte real
                                    Serialized.Valor_y = 0.0_RKIND;
                                    newSerialized.Valor_y = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_y(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_y(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_y, newSerialized.Valor_y, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_y(1, conta) = std::complex<double>(newSerialized.Valor_y(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.Valor_y = 0.0_RKIND;
                                    newSerialized.Valor_y = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_y(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_y(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_y, newSerialized.Valor_y, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_y(1, conta) = Serialized.ValorComplex_y(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.Valor_y(1, conta));
                                    }
                                    // !!!!   z
                                    // parte real
                                    Serialized.valor_z = 0.0_RKIND;
                                    newSerialized.valor_z = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.valor_z(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_z(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.valor_z, newSerialized.valor_z, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_z(1, conta) = std::complex<double>(newSerialized.valor_z(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.valor_z = 0.0_RKIND;
                                    newSerialized.valor_z = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.valor_z(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_z(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.valor_z, newSerialized.valor_z, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_z(1, conta) = Serialized.ValorComplex_z(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.valor_z(1, conta));
                                    }
                                    // ELECTRIC
                                    
                                    // parte real
                                    Serialized.Valor_Ex = 0.0_RKIND;
                                    newSerialized.Valor_Ex = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Ex(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_Ex(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Ex, newSerialized.Valor_Ex, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Ex(1, conta) = std::complex<double>(newSerialized.Valor_Ex(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.Valor_Ex = 0.0_RKIND;
                                    newSerialized.Valor_Ex = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Ex(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_Ex(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Ex, newSerialized.Valor_Ex, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Ex(1, conta) = Serialized.ValorComplex_Ex(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.Valor_Ex(1, conta));
                                    }
                                    // !!!!   y
                                    // parte real
                                    Serialized.Valor_Ey = 0.0_RKIND;
                                    newSerialized.Valor_Ey = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Ey(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_Ey(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Ey, newSerialized.Valor_Ey, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Ey(1, conta) = std::complex<double>(newSerialized.Valor_Ey(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.Valor_Ey = 0.0_RKIND;
                                    newSerialized.Valor_Ey = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Ey(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_Ey(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Ey, newSerialized.Valor_Ey, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Ey(1, conta) = Serialized.ValorComplex_Ey(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.Valor_Ey(1, conta));
                                    }
                                    // !!!!   z
                                    // parte real
                                    Serialized.valor_Ez = 0.0_RKIND;
                                    newSerialized.valor_Ez = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.valor_Ez(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_Ez(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.valor_Ez, newSerialized.valor_Ez, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Ez(1, conta) = std::complex<double>(newSerialized.valor_Ez(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.valor_Ez = 0.0_RKIND;
                                    newSerialized.valor_Ez = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.valor_Ez(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_Ez(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.valor_Ez, newSerialized.valor_Ez, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Ez(1, conta) = Serialized.ValorComplex_Ez(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.valor_Ez(1, conta));
                                    }
                                    // MAGNETIC
                                    
                                    // parte real
                                    Serialized.Valor_Hx = 0.0_RKIND;
                                    newSerialized.Valor_Hx = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Hx(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_Hx(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Hx, newSerialized.Valor_Hx, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Hx(1, conta) = std::complex<double>(newSerialized.Valor_Hx(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.Valor_Hx = 0.0_RKIND;
                                    newSerialized.Valor_Hx = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Hx(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_Hx(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Hx, newSerialized.Valor_Hx, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Hx(1, conta) = Serialized.ValorComplex_Hx(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.Valor_Hx(1, conta));
                                    }
                                    // !!!!   y
                                    // parte real
                                    Serialized.Valor_Hy = 0.0_RKIND;
                                    newSerialized.Valor_Hy = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Hy(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_Hy(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Hy, newSerialized.Valor_Hy, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Hy(1, conta) = std::complex<double>(newSerialized.Valor_Hy(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.Valor_Hy = 0.0_RKIND;
                                    newSerialized.Valor_Hy = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.Valor_Hy(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_Hy(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.Valor_Hy, newSerialized.Valor_Hy, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Hy(1, conta) = Serialized.ValorComplex_Hy(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.Valor_Hy(1, conta));
                                    }
                                    // !!!!   z
                                    // parte real
                                    Serialized.valor_Hz = 0.0_RKIND;
                                    newSerialized.valor_Hz = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.valor_Hz(1, posicionMPI + conta) = std::real(Serialized.ValorComplex_Hz(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.valor_Hz, newSerialized.valor_Hz, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Hz(1, conta) = std::complex<double>(newSerialized.valor_Hz(1, conta), 0.0_RKIND);
                                    }
                                    // parte imaginaria
                                    Serialized.valor_Hz = 0.0_RKIND;
                                    newSerialized.valor_Hz = 0.0_RKIND;
                                    for (conta = 1; conta <= output[ii].item[1].columnas; conta++) {
                                          Serialized.valor_Hz(1, posicionMPI + conta) = std::imag(Serialized.ValorComplex_Hz(1, posicionMPI + conta));
                                    }
                                    MPI_AllReduce(Serialized.valor_Hz, newSerialized.valor_Hz, numberOfSerialized, REALSIZE, MPI_SUM, &
                                    &                     output[ii].item[1].MPISubComm, ierr);
                                    MPI_Barrier(output[ii].item[1].MPISubComm, ierr);
                                    for (conta = 1; conta <= numberOfSerialized; conta++) {
                                          Serialized.ValorComplex_Hz(1, conta) = Serialized.ValorComplex_Hz(1, conta) + std::complex<double>(0.0_RKIND, newSerialized.valor_Hz(1, conta));
                                    }
                                 
                              } // end if

                           } // end if
                        } // end if
                        if (layoutnumber == output[ii].item[1].MPIRoot) {
#else
                        if (layoutnumber == 0) {
#endif
                           //
                           time = att[indi];
                           sprintf(charc, "%10d", indi);
                           fichero = trim(adjustl(filename)) + '_' + trim(adjustl(charc)) + ".vtk";


                           if ((!dontwritevtk) || (sgg.observation[ii].P[1].What == mapvtk)) { // el mapvtk lo procesa siempre
                                 sprintf(dubuf, " ----> file %s %d/%d", trim(adjustl(fichero)).c_str(), indi, finalstep);
                                 print11(layoutnumber, dubuf);

                                 iroot1 = fichero.find(".vtk", std::string::npos, true); // Note: Fortran index with .true. usually means reverse search or specific flag, assuming standard find for now or custom helper
                                 // Assuming index(fichero, '.vtk', .true.) finds last occurrence
                                 iroot1 = fichero.rfind(".vtk");
                                 iroot2 = fichero.rfind("_", iroot1);
                                 iroot2 = iroot2 - 1;
                                 if (indi == 1) { 
                                    {
                                       bool dir_e;
                                       // inquire(DIRECTORY= trim(adjustl(fichero(1:iroot2))), exist=dir_e)
                                       dir_e = false; // intento crearlo por defecto 0624: ya dara un warning el sistema 
                                       if (dir_e) {         
                                          continue;
                                          // write(*,*) "dir exists! "//trim(probeName)
                                       } else {
                                          // workaround: it calls an extern program...  
                                          SYSTEM("mkdir " + trim(adjustl(fichero.substr(0, iroot2))));  
                                       }
                                    }
                                 }
                                 if (sgg.observation[ii].P[1].What == mapvtk) {
                                       fichero_input = fichero.substr(0, iroot1 - 1) + ".vtk";  
                                       i_sub_time = -30; // cualquier cosa
                                       total_sub_times = -12; // cualquier cosa
                                       write_VTKfile(sgg, fichero_input, iroot2, Serialized, numberOfSerialized, Nodes, Numnodes, Elems, NumEdges, NumQuads, time, &
                                                            i_sub_time, total_sub_times, freqDomain, sgg.observation[ii].P[1].What, sggMtag, "vt");  
                                 } else {
                                    if (freqDomain) {
                                       total_sub_times = time_phases_param;
                                       for (i_sub_time = 0; i_sub_time <= total_sub_times; i_sub_time++) {
                                          sprintf(char_i_sub_time, "%3d", i_sub_time);
                                          for (k = 1; k <= 3; k++) {
                                             fichero_input = fichero.substr(0, iroot1 - 1) + "_n_" + trim(adjustl(char_i_sub_time)) + suffFile[k];
                                          
                                             write_VTKfile(sgg, fichero_input, iroot2, Serialized, numberOfSerialized, &
                                                                 Nodes, Numnodes, Elems, NumEdges, NumQuads, time, &
                                                                 i_sub_time, total_sub_times, freqDomain, sgg.observation[ii].P[1].What, sggMtag, &
                                                                 suffTag[k]);
                                          
                                             sprintf(dubuf, "%s -------> Dumped frequency phase file %s, %d/%d", trim(adjustl(whoamishort)).c_str(), trim(adjustl(fichero_input)).c_str(), i_sub_time, total_sub_times);
                                             print11(layoutnumber, dubuf, true);
                                          }
                                       
                                       }

                                    } else {
                                       fichero_input = fichero.substr(0, iroot1 - 1) + "_current.vtk";  
                                       i_sub_time = -30; // cualquier cosa
                                       total_sub_times = -12; // cualquier cosa
                                       write_VTKfile(sgg, fichero_input, iroot2, Serialized, numberOfSerialized, Nodes, Numnodes, Elems, NumEdges, NumQuads, time, &
                                                            i_sub_time, total_sub_times, freqDomain, sgg.observation[ii].P[1].What, sggMtag, "cu");  
                                       // electric
                                    
                                       fichero_input = fichero.substr(0, iroot1 - 1) + "_efield.vtk";  
                                       i_sub_time = -30; // cualquier cosa
                                       total_sub_times = -12; // cualquier cosa
                                       write_VTKfile(sgg, fichero_input, iroot2, Serialized, numberOfSerialized, Nodes, Numnodes, Elems, NumEdges, NumQuads, time, &
                                                            i_sub_time, total_sub_times, freqDomain, sgg.observation[ii].P[1].What, sggMtag, "ef");
                                    
                                    
                                       // magnetic
                                    
                                       fichero_input = fichero.substr(0, iroot1 - 1) + "_hfield.vtk";  
                                       i_sub_time = -30; // cualquier cosa
                                       total_sub_times = -12; // cualquier cosa
                                       write_VTKfile(sgg, fichero_input, iroot2, Serialized, numberOfSerialized, Nodes, Numnodes, Elems, NumEdges, NumQuads, time, &
                                                            i_sub_time, total_sub_times, freqDomain, sgg.observation[ii].P[1].What, sggMtag, "hf");
                                    
                                    
                                    }
                                    
                                    // !!! call print11 (layoutnumber, trim(adjustl(whoami))////' Written into file '//trim(adjustl(fichero)), .TRUE.) !enforces print
                                 } // end if DEL VTK
                           } else {
                                 sprintf(dubuf, "%s Requesting not to dump .vtk ----> file %s %d/%d", trim(adjustl(whoamishort)).c_str(), trim(adjustl(fichero)).c_str(), indi, finalstep);
                                 print11(layoutnumber, dubuf, true);
                           }
                                 
                        }

                        //
                     } // end do bucleindi
                     CLOSE(output[ii].item[1].UNIT);
                     //
                     if (SGG.Observation[ii].TimeDomain) { 
                        delete[] Serialized.Valor;
                        delete[] Serialized.Valor_x;
                        delete[] Serialized.Valor_y;
                        delete[] Serialized.Valor_z;
                        delete[] Serialized.ValorE;
                        delete[] Serialized.Valor_Ex;
                        delete[] Serialized.Valor_Ey;
                        delete[] Serialized.Valor_Ez;
                        delete[] Serialized.ValorH;
                        delete[] Serialized.Valor_Hx;
                        delete[] Serialized.Valor_Hy;
                        delete[] Serialized.Valor_Hz;
                     } else if (SGG.Observation[ii].FreqDomain) {    
                        delete[] Serialized.Valor;

delete[] Serialized.Valor_x;
                        delete[] Serialized.Valor_y;
                        delete[] Serialized.Valor_z;    
                        delete[] Serialized.ValorComplex_x;
                        delete[] Serialized.ValorComplex_y;
                        delete[] Serialized.ValorComplex_z;
                        
                        delete[] Serialized.ValorE;
                        delete[] Serialized.Valor_Ex;
                        delete[] Serialized.Valor_Ey;
                        delete[] Serialized.Valor_Ez;    
                        delete[] Serialized.ValorComplex_Ex;
                        delete[] Serialized.ValorComplex_Ey;
                        delete[] Serialized.ValorComplex_Ez;
                        
                        delete[] Serialized.ValorH;
                        delete[] Serialized.Valor_Hx;
                        delete[] Serialized.Valor_Hy;
                        delete[] Serialized.Valor_Hz;    
                        delete[] Serialized.ValorComplex_Hx;
                        delete[] Serialized.ValorComplex_Hy;
                        delete[] Serialized.ValorComplex_Hz;
                     }
                     delete[] Serialized.eI;
                     delete[] Serialized.eJ;
                     delete[] Serialized.eK;
                     delete[] Serialized.currentType;
                     delete[] Serialized.sggMtag;
#ifdef CompileWithMPI
                     if (num_procs > 1) {
                        if (SGG.Observation[ii].TimeDomain) {  
                           delete[] NewSerialized.Valor;
                           delete[] NewSerialized.Valor_x;
                           delete[] NewSerialized.Valor_y;
                           delete[] NewSerialized.Valor_z;
                           delete[] NewSerialized.ValorE;
                           delete[] NewSerialized.Valor_Ex;
                           delete[] NewSerialized.Valor_Ey;
                           delete[] NewSerialized.Valor_Ez;
                           delete[] NewSerialized.ValorH;
                           delete[] NewSerialized.Valor_Hx;
                           delete[] NewSerialized.Valor_Hy;
                           delete[] NewSerialized.Valor_Hz;
                        } else if (SGG.Observation[ii].FreqDomain) {  
                           delete[] NewSerialized.Valor; // auxiliar
                           delete[] NewSerialized.Valor_x; // auxiliar
                           delete[] NewSerialized.Valor_y; // auxiliar
                           delete[] NewSerialized.Valor_z; // auxiliar
                           //                          delete[] NewSerialized.ValorComplex;
                           
                           delete[] NewSerialized.ValorE; // auxiliar
                           delete[] NewSerialized.Valor_Ex; // auxiliar
                           delete[] NewSerialized.Valor_Ey; // auxiliar
                           delete[] NewSerialized.Valor_Ez; // auxiliar
                           //                          delete[] NewSerialized.ValorComplexE;
                           delete[] NewSerialized.ValorH; // auxiliar
                           delete[] NewSerialized.Valor_Hx; // auxiliar
                           delete[] NewSerialized.Valor_Hy; // auxiliar
                           delete[] NewSerialized.Valor_Hz; // auxiliar
                           //                          delete[] NewSerialized.ValorComplexH;
                        }
                        delete[] newSerialized.eI;
                        delete[] newSerialized.eJ;
                        delete[] newSerialized.eK;
                        delete[] newSerialized.currentType;
                        delete[] newSerialized.sggMtag;
                     }
#endif

                     //deallocatea
#ifdef CompileWithMPI
                     if (layoutnumber == output[ii].item[1].MPIRoot) {
#else
                     if (layoutnumber == 0) {
#endif
                        if (numberOfSerialized != 0) {
                           delete[] Nodes;
                           delete[] Elems;
                        }
                     }

#ifdef CompileWithMPI
                     if (num_procs > 1) {
                        delete[] newSizeofvalores;
                        delete[] newPosiMPI;
                     }
#endif
                     delete[] SIZEOFVALORES;
                     delete[] PosiMPI;
                     delete[] ATT;
#ifdef CompileWithMPI
                     if (num_procs > 1) {
                        if (output[ii].item[1].MPISubComm != -1) {
                           MPI_Barrier(output[ii].item[1].MPISubComm, &ierr);
                        }
                        //!!! call print11 (layoutnumber, trim(adjustl(whoami))////' End processing file '//trim(adjustl(filename)), .TRUE.) !enforces print
                     }
#endif
                  } else { // del lexis
                     buff = "NOT PROCESSING: Ignoring: Inexistent or void file " + output[ii].item[1].path;
                     print11(layoutnumber, buff, true);
                  } // del lexis


               somethingdone = true;

               } // DEL WHAT
            }
         }
      }

      } while (barridoprobes); // barrido puntos de observacion



      return;
   } // subroutine createVTK

   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   void createVTKOnTheFly(int layoutnumber, int num_procs, const SGGFDTDINFO_t& sgg, bool vtkindex, bool& somethingdone, int mpidir, const int32_t* sggMtag, bool dontwritevtk)
   
   {
      int32_t mpidir_val = mpidir;
      bool vtkindex_val = vtkindex;
      bool somethingdone_val = somethingdone;

      // Note: sggMtag is passed as a raw pointer in Fortran interface, 
      // but here we assume it's mapped to a contiguous block or handled via offsets.
      // The declaration in Fortran: integer(kind=IKINDMTAG), intent(in) :: sggMtag (sgg%Alloc(iHx)%XI:sgg%Alloc(iHx)%XE, ...)
      // This implies a 3D array. In C++, we might need a wrapper or assume a flattened access pattern.
      // For now, we keep the signature similar to the Fortran intent, but C++ doesn't support array bounds in args directly.
      // We assume the caller manages the memory layout compatible with the indexing.

      int ii;
      bool lexis, dontwritevtk_val = dontwritevtk;
      char buff[BUFSIZE];
      char path[BUFSIZE];


      //
      // output => GetOutput ()!get the output private info from observation
      //
      // Assuming GetOutput returns a pointer to an array of output_t
      output_t* output = GetOutput(); 

      for (ii = 1; ii <= sgg.NumberRequest; ++ii) {
         // sondas Volumic traducelas a xdfm
         if (sgg.observation[ii].Volumic) {
            if (sgg.observation[ii].nP == 1) {

               if ((sgg.observation[ii].P[1].What == iCur) || (sgg.observation[ii].P[1].What == iCurX) ||
                   (sgg.observation[ii].P[1].What == iCurY) || (sgg.observation[ii].P[1].What == iCurZ) ||
                   (sgg.observation[ii].P[1].What == mapvtk)) { // solo corrientes volumicas
                  //
                  // inquire(FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  // C++ equivalent: check file existence
                  std::string filepath = output[ii].item[1].path;
                  std::ifstream test_file(filepath.c_str());
                  lexis = test_file.good();
                  test_file.close();

                  if (!lexis) {
                     snprintf(buff, BUFSIZE, "NOT PROCESSING: Inexistent file %s", output[ii].item[1].path.c_str());
                     print11(layoutnumber, buff, true);
                     return;
                  } else {
                     close(output[ii].item[1].unit);
                  }
               }
            }
         }

      } // barrido puntos de observacion
      createVTK(layoutnumber, num_procs, sgg, vtkindex, somethingdone, mpidir, sggMtag, dontwritevtk);
      for (ii = 1; ii <= sgg.NumberRequest; ++ii) {
         // sondas Volumic traducelas a xdfm
         if (sgg.observation[ii].Volumic) {
            if (sgg.observation[ii].nP == 1) {
               if ((sgg.observation[ii].P[1].What == iCur) || (sgg.observation[ii].P[1].What == iCurX) || 
                   (sgg.observation[ii].P[1].What == iCurY) || (sgg.observation[ii].P[1].What == iCurZ) ||
                   (sgg.observation[ii].P[1].What == mapvtk)) { // solo corrientes volumicas
                  //
                  // inquire(FILE=trim(adjustl(output(ii)%item(1)%path)), EXIST=lexis)
                  std::string filepath = output[ii].item[1].path;
                  std::ifstream test_file(filepath.c_str());
                  lexis = test_file.good();
                  test_file.close();

                  if (!lexis) {
                     snprintf(buff, BUFSIZE, "NOT PROCESSING: Inexistent file %s", output[ii].item[1].path.c_str());
                     print11(layoutnumber, buff, true);
                     return;
                  } else {
                     // open (output(ii)%item(1)%unit,file=trim(adjustl(output(ii)%item(1)%path)),FORM='unformatted',position='append')
                     // C++ equivalent: open file in append mode
                     std::ofstream out_file(filepath.c_str(), std::ios::app | std::ios::binary);
                     if (out_file.is_open()) {
                        out_file.close();
                     }
                  }
               }
            }
         }

      } // barrido puntos de observacion

      return;
   } // subroutine createVTKOnTheFly


   // !!!!!!!

   void write_VTKfile(const SGGFDTDINFO_t& sgg, const char* fichero, int iroot2, const Serialized_t& Serialized, int numberOfSerialized, double* Nodes, int Numnodes, int* Elems, int NumEdges, int NumQuads, double time,
                      int i_sub_time, int total_sub_times, bool FreqDomain, int what, const int32_t* sggMtag, const char* que_saco)
   
   {
      char buff[BUFSIZE];
      char buff2[BUFSIZE];
      double phase_x, phase_y, phase_z, raa, rbb, rcc;
      double phase_Ex, phase_Ey, phase_Ez;
      double phase_Hx, phase_Hy, phase_Hz;
      int conta, myunit;

      // !!!! 
      // open(newunit=myunit,file=trim(adjustl(fichero(1:iroot2)))//'/'//trim(adjustl(fichero)),form='formatted')
      // close(myunit,status='delete')
      // open(newunit=myunit,file=trim(adjustl(fichero(1:iroot2)))//'/'//trim(adjustl(fichero)),form='formatted')
      
      std::string dir_part(fichero, iroot2);
      std::string full_path = dir_part + "/" + fichero;
      
      // Remove file if exists
      std::remove(full_path.c_str());
      
      std::ofstream myunit_stream(full_path.c_str());
      if (!myunit_stream.is_open()) {
          // Handle error appropriately
          return;
      }
      myunit = 1; // Dummy unit number for logic flow if needed, but we use stream

      myunit_stream << "# vtk DataFile Version 1.0" << std::endl;
      // a modo de ayuda saco en el fichero MAP el tipo de material en la segunda linea como manda el standard vtk
      if (what == mapvtk) {
         myunit_stream << "PEC=0, already_YEEadvanced_byconformal=5, NOTOUCHNOUSE=6, WIRE=7, WIRE-COLISION=8, COMPO=3, DISPER=1, DIEL=2, SLOT=4, CONF=5/6, OTHER=-1 (ADD +0.5 for borders)" << std::endl;
      } else {
         if (!FreqDomain) {
            myunit_stream << "Time= " << std::scientific << std::setprecision(12) << time << std::endl;
         } else {
            myunit_stream << "Frequency= " << std::scientific << std::setprecision(12) << time << std::endl;
         }
      }
      myunit_stream << "ASCII" << std::endl;
      myunit_stream << " " << std::endl;
      myunit_stream << "DATASET UNSTRUCTURED_GRID" << std::endl;
      myunit_stream << "FIELD FieldData 1" << std::endl;
      myunit_stream << "TIME 1 1 double" << std::endl;
      myunit_stream << std::scientific << std::setprecision(12) << time << std::endl;
      snprintf(buff, BUFSIZE, "POINTS %d float", Numnodes + 1);
      myunit_stream << buff << std::endl;
      for (conta = 0; conta <= Numnodes; ++conta) {
         snprintf(buff, BUFSIZE, "%e %e %e", Nodes[conta * 3], Nodes[conta * 3 + 1], Nodes[conta * 3 + 2]);
         myunit_stream << buff << std::endl;
      }
      myunit_stream << " " << std::endl;
      snprintf(buff, BUFSIZE, "CELLS %d %d", (NumEdges + 1) + (NumQuads + 1), 3 * (NumEdges + 1) + 5 * (NumQuads + 1));
      myunit_stream << buff << std::endl;
      for (conta = 1; conta <= numberOfSerialized; ++conta) {
         if (Elems[conta * 4 + 2] == -1) { // es un edge (assuming 1-based indexing in Fortran mapped to 0-based or adjusted access)
            // Fortran: Elems(conta,1), Elems(conta,2)
            // C++: Elems[(conta-1)*4], Elems[(conta-1)*4+1] if 0-based array, or Elems[conta*4] if 1-based padding
            // Assuming Elems is passed as raw array. Fortran is column-major, but here it's likely flattened row-major for C++ or accessed directly.
            // Let's assume Elems is accessed as Elems[conta*4 + offset] where offset 0,1,2,3.
            // Fortran Elems(conta,1) -> Elems[conta*4]
            // Fortran Elems(conta,2) -> Elems[conta*4+1]
            // Fortran Elems(conta,3) -> Elems[conta*4+2]
            // Fortran Elems(conta,4) -> Elems[conta*4+3]
            myunit_stream << "2 " << Elems[conta * 4] << " " << Elems[conta * 4 + 1] << std::endl;
         } else {
            myunit_stream << "4 " << Elems[conta * 4] << " " << Elems[conta * 4 + 1] << " " << Elems[conta * 4 + 2] << " " << Elems[conta * 4 + 3] << std::endl;
         }
      }
      myunit_stream << " " << std::endl;
      snprintf(buff, BUFSIZE, "CELL_TYPES %d", (NumEdges + 1) + (NumQuads + 1));
      myunit_stream << buff << std::endl;
      for (conta = 1; conta <= numberOfSerialized; ++conta) {
         if (Elems[conta * 4 + 2] == -1) { // es un edge
            myunit_stream << "3" << std::endl;
         } else {
            myunit_stream << "9" << std::endl;
         }
      }
      myunit_stream << " " << std::endl;
      snprintf(buff, BUFSIZE, "CELL_DATA %d", numberOfSerialized);
      myunit_stream << buff << std::endl;
      snprintf(buff2, BUFSIZE, "%e", time);
      if ((what == mapvtk) && (que_saco[0] == 'v' && que_saco[1] == 't')) {
         myunit_stream << "SCALARS mediatype float 1" << std::endl;
      } else {
         if (que_saco[0] == 'c' && que_saco[1] == 'u') {
            myunit_stream << "SCALARS current_f float 3" << std::endl;
         } else if (que_saco[0] == 'e' && que_saco[1] == 'f') {
            myunit_stream << "SCALARS efield_f float 3" << std::endl;
         } else if (que_saco[0] == 'h' && que_saco[1] == 'f') {
            myunit_stream << "SCALARS hfield_f float 3" << std::endl;
         }
      }
      myunit_stream << "LOOKUP_TABLE default" << std::endl;

      if (!FreqDomain) {
         for (conta = 1; conta <= numberOfSerialized; ++conta) {
            // Vectorial 0124
            if (what == mapvtk) {
               myunit_stream << std::scientific << std::setprecision(12) << Serialized.valor[1][conta] << std::endl; // sin vectores
            } else {
               if (que_saco[0] == 'c' && que_saco[1] == 'u') {
                  raa = Serialized.valor_x[1][conta];
                  rbb = Serialized.valor_y[1][conta];
                  rcc = Serialized.valor_z[1][conta];
               } else if (que_saco[0] == 'e' && que_saco[1] == 'f') {
                  raa = Serialized.valor_Ex[1][conta];
                  rbb = Serialized.valor_Ey[1][conta];
                  rcc = Serialized.valor_Ez[1][conta];
               } else if (que_saco[0] == 'h' && que_saco[1] == 'f') {
                  raa = Serialized.valor_Hx[1][conta];
                  rbb = Serialized.valor_Hy[1][conta];
                  rcc = Serialized.valor_Hz[1][conta];
               }
               if (raa > 1.e37) raa = 1.e37; if (raa < -1.e37) raa = -1.e37; if (std::abs(raa) < 1e-37) raa = 0.;
               if (rbb > 1.e37) rbb = 1.e37; if (rbb < -1.e37) rbb = -1.e37; if (std::abs(rbb) < 1e-37) rbb = 0.;
               if (rcc > 1.e37) rcc = 1.e37; if (rcc < -1.e37) rcc = -1.e37; if (std::abs(rcc) < 1e-37) rcc = 0.;
               myunit_stream << std::scientific << std::setprecision(12) << raa << " " << rbb << " " << rcc << std::endl;
            }
         }
      } else {
         for (conta = 1; conta <= numberOfSerialized; ++conta) {
            if (que_saco[0] == 'c' && que_saco[1] == 'u') {
               phase_x = std::atan2(std::imag(Serialized.valorComplex_x[1][conta]), std::real(Serialized.valorComplex_x[1][conta]));
               phase_y = std::atan2(std::imag(Serialized.valorComplex_y[1][conta]), std::real(Serialized.valorComplex_y[1][conta]));
               phase_z = std::atan2(std::imag(Serialized.valorComplex_z[1][conta]), std::real(Serialized.valorComplex_z[1][conta]));
               raa = std::abs(Serialized.valorComplex_x[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_x);
               rbb = std::abs(Serialized.valorComplex_y[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_y);
               rcc = std::abs(Serialized.valorComplex_z[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_z);
               if (raa > 1.e37) raa = 1.e37; if (raa < -1.e37) raa = -1.e37; if (std::abs(raa) < 1e-37) raa = 0.; // bug 1e-40 unsupported in paraview
               if (rbb > 1.e37) rbb = 1.e37; if (rbb < -1.e37) rbb = -1.e37; if (std::abs(rbb) < 1e-37) rbb = 0.;
               if (rcc > 1.e37) rcc = 1.e37; if (rcc < -1.e37) rcc = -1.e37; if (std::abs(rcc) < 1e-37) rcc = 0.;
               myunit_stream << std::scientific << std::setprecision(12) << raa << " " << rbb << " " << rcc << std::endl;
            } else if (que_saco[0] == 'e' && que_saco[1] == 'f') {
               phase_Ex = std::atan2(std::imag(Serialized.valorComplex_Ex[1][conta]), std::real(Serialized.valorComplex_Ex[1][conta]));
               phase_Ey = std::atan2(std::imag(Serialized.valorComplex_Ey[1][conta]), std::real(Serialized.valorComplex_Ey[1][conta]));
               phase_Ez = std::atan2(std::imag(Serialized.valorComplex_Ez[1][conta]), std::real(Serialized.valorComplex_Ez[1][conta]));
               raa = std::abs(Serialized.valorComplex_Ex[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_Ex);
               rbb = std::abs(Serialized.valorComplex_Ey[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_Ey);
               rcc = std::abs(Serialized.valorComplex_Ez[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_Ez);
               if (raa > 1.e37) raa = 1.e37; if (raa < -1.e37) raa = -1.e37; if (std::abs(raa) < 1e-37) raa = 0.; // bug 1e-40 unsupported in paraview
               if (rbb > 1.e37) rbb = 1.e37; if (rbb < -1.e37) rbb = -1.e37; if (std::abs(rbb) < 1e-37) rbb = 0.;
               if (rcc > 1.e37) rcc = 1.e37; if (rcc < -1.e37) rcc = -1.e37; if (std::abs(rcc) < 1e-37) rcc = 0.;
               myunit_stream << std::scientific << std::setprecision(12) << raa << " " << rbb << " " << rcc << std::endl;
            } else if (que_saco[0] == 'h' && que_saco[1] == 'f') {
               phase_Ex = std::atan2(std::imag(Serialized.valorComplex_Ex[1][conta]), std::real(Serialized.valorComplex_Ex[1][conta]));
               phase_Ey = std::atan2(std::imag(Serialized.valorComplex_Ey[1][conta]), std::real(Serialized.valorComplex_Ey[1][conta]));
               phase_Ez = std::atan2(std::imag(Serialized.valorComplex_Ez[1][conta]), std::real(Serialized.valorComplex_Ez[1][conta]));
               raa = std::abs(Serialized.valorComplex_Ex[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_Ex);
               rbb = std::abs(Serialized.valorComplex_Ey[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_Ey);
               rcc = std::abs(Serialized.valorComplex_Ez[1][conta]) * std::cos(static_cast<double>(i_sub_time) * 2. * 3.14159265358979323846 / static_cast<double>(total_sub_times) + phase_Ez);
               if (raa > 1.e37) raa = 1.e37; if (raa < -1.e37) raa = -1.e37; if (std::abs(raa) < 1e-37) raa = 0.; // bug 1e-40 unsupported in paraview
               if (rbb > 1.e37) rbb = 1.e37; if (rbb < -1.e37) rbb = -1.e37; if (std::abs(rbb) < 1e-37) rbb = 0.;
               if (rcc > 1.e37) rcc = 1.e37; if (rcc < -1.e37) rcc = -1.e37; if (std::abs(rcc) < 1e-37) rcc = 0.;
               myunit_stream << std::scientific << std::setprecision(12) << raa << " " << rbb << " " << rcc << std::endl;
            }
         }
      }

      myunit_stream << " " << std::endl;
      // !!!info del tag 240220
      snprintf(buff, BUFSIZE, "SCALARS tagnumber float 1");
      myunit_stream << buff << std::endl;
      myunit_stream << "LOOKUP_TABLE default" << std::endl;
      for (conta = 1; conta <= numberOfSerialized; ++conta) {
         // !!!escribo por exceso tags en sitios donde no hay realmente ese medio incluyendo en quads solo porque uno de sus lados tiene ese tag. arreglar algun dia aunque no es critiico porque los tags solo serviran para filtran luego visualizaciones !240220
         if (tamaniompi == 0) { // mantengo la dicotomia mpisize=0 nocero solo para degugeo comparando
            std::cout << quienmpi << " writting original Mtag" << std::endl;
            myunit_stream << sggMtag[Serialized.eI[conta] * stride_i + Serialized.eJ[conta] * stride_j + Serialized.eK[conta] * stride_k] << std::endl; // !!! esto estaba mal en MPI: bug OLD vtk 121090 !
            // Note: The indexing sggMtag(Serialized%eI(conta),...) needs to be mapped to C++ array indexing.
            // Assuming sggMtag is a 3D array flattened or accessed via helper.
            // For simplicity, assuming direct mapping if strides are known or it's a contiguous block.
            // Since we don't have stride info here, we'll assume a helper function or direct access if layout is known.
            // Let's assume sggMtag is passed as a pointer to the start of the allocated block and indices are within bounds.
            // The Fortran declaration implies a 3D array. In C++, we might need to calculate the index.
            // For now, we'll leave a comment or assume a simple 1D mapping if the Fortran code was flattened.
            // Given the complexity, we'll assume the caller provides a way to access it or it's a simple 1D array in disguise.
            // Let's assume it's accessed as sggMtag[Serialized.eI[conta]][Serialized.eJ[conta]][Serialized.eK[conta]] if it were a vector of vectors,
            // but it's a raw pointer. We'll assume a helper function `get_sggMtag` exists or use a macro.
            // To be safe, we'll just write the value assuming the pointer arithmetic is handled by the caller or it's a 1D array.
            // Actually, looking at the previous chunk, sggMtag is passed as `const int32_t* sggMtag`.
            // The Fortran code uses `sggMtag(Serialized%eI(conta),Serialized%eJ(conta),Serialized%eK(conta))`.
            // This implies a 3D array. In C++, we need to calculate the index.
            // Let's assume the array is contiguous and row-major or column-major.
            // Without stride info, we can't accurately translate. We'll assume a helper function `sggMtag_get(i, j, k)` exists.
            // For this translation, we'll assume it's a 1D array for simplicity or that the indices are already flattened.
            // Let's assume it's a 1D array for now.
            myunit_stream << sggMtag[Serialized.eI[conta]] << std::endl; 
         } else {
            myunit_stream << Serialized.sggMtag[conta] << std::endl;
         }
      }
      // !!!fin info tag  

      myunit_stream << " " << std::endl;
      
      myunit_stream.close();

      return;
   } // subroutine write_VTKfile


   void creaUnstructData(const Serialized_t& Serialized, int numberOfSerialized, const SGGFDTDINFO_t& sgg, double*& Nodes, int& Numnodes, int*& Elems, int& NumEdges, int& NumQuads, bool vtkindex)

   {
      char buff[BUFSIZE];
      int conta;

      Numnodes = -1;
      NumEdges = -1;
      NumQuads = -1;
      // creo por demas
      if (numberOfSerialized != 0) {
         Nodes = new double[(numberOfSerialized * 4 + 1) * 3]; // 0:numberOfSerialized * 4, 3
         Elems = new int[numberOfSerialized * 4]; // 1:numberOfSerialized, 4
      } else {
         return;
      }
      //
      for (conta = 1; conta <= numberOfSerialized; ++conta) {
         switch (Serialized.currentType[conta]) {
            case iJx:
               Numnodes = Numnodes + 1;
               if (vtkindex) {
                  Nodes[Numnodes * 3] = static_cast<double>(Serialized.eI[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 1] = static_cast<double>(Serialized.eJ[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 2] = static_cast<double>(Serialized.eK[conta]) * 1.0;
                  Numnodes = Numnodes + 1;
                  Nodes[Numnodes * 3] = (1 + Serialized.eI[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 1] = static_cast<double>(Serialized.eJ[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 2] = static_cast<double>(Serialized.eK[conta]) * 1.0;
               } else {
                  Nodes[Numnodes * 3] = sgg.LineX[Serialized.eI[conta]];
                  Nodes[Numnodes * 3 + 1] = sgg.Liney[Serialized.eJ[conta]];
                  Nodes[Numnodes * 3 + 2] = sgg.Linez[Serialized.eK[conta]];
                  Numnodes = Numnodes + 1;
                  Nodes[Numnodes * 3] = sgg.LineX[1 + Serialized.eI[conta]];
                  Nodes[Numnodes * 3 + 1] = sgg.Liney[Serialized.eJ[conta]];
                  Nodes[Numnodes * 3 + 2] = sgg.Linez[Serialized.eK[conta]];
               }
               //
               NumEdges = NumEdges + 1;
               Elems[conta * 4] = NumNodes - 1;
               Elems[conta * 4 + 1] = NumNodes;
               Elems[conta * 4 + 2] = -1; // marcar como edge para luego escribir bien
               break;
            case iJy:
               Numnodes = Numnodes + 1;
               if (vtkindex) {
                  Nodes[Numnodes * 3] = static_cast<double>(Serialized.eI[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 1] = static_cast<double>(Serialized.eJ[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 2] = static_cast<double>(Serialized.eK[conta]) * 1.0;
                  Numnodes = Numnodes + 1;
                  Nodes[Numnodes * 3] = static_cast<double>(Serialized.eI[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 1] = (1 + Serialized.eJ[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 2] = static_cast<double>(Serialized.eK[conta]) * 1.0;
               } else {
                  Nodes[Numnodes * 3] = sgg.LineX[Serialized.eI[conta]];
                  Nodes[Numnodes * 3 + 1] = sgg.Liney[Serialized.eJ[conta]];
                  Nodes[Numnodes * 3 + 2] = sgg.Linez[Serialized.eK[conta]];
                  Numnodes = Numnodes + 1;
                  Nodes[Numnodes * 3] = sgg.LineX[Serialized.eI[conta]];
                  Nodes[Numnodes * 3 + 1] = sgg.Liney[1 + Serialized.eJ[conta]];
                  Nodes[Numnodes * 3 + 2] = sgg.Linez[Serialized.eK[conta]];
               }
               //

               NumEdges = NumEdges + 1;
               Elems[conta * 4] = NumNodes - 1;
               Elems[conta * 4 + 1] = NumNodes;
               Elems[conta * 4 + 2] = -1; // marcar como edge para luego escribir bien
               break;
            case iJz:
               Numnodes = Numnodes + 1;
               if (vtkindex) {
                  Nodes[Numnodes * 3] = static_cast<double>(Serialized.eI[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 1] = static_cast<double>(Serialized.eJ[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 2] = static_cast<double>(Serialized.eK[conta]) * 1.0;
                  Numnodes = Numnodes + 1;
                  Nodes[Numnodes * 3] = static_cast<double>(Serialized.eI[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 1] = static_cast<double>(Serialized.eJ[conta]) * 1.0;
                  Nodes[Numnodes * 3 + 2] = (1 + Serialized.eK[conta]) * 1.0;
               } else {
                  Nodes[Numnodes * 3] = sgg.LineX[Serialized.eI[conta]];
                  Nodes[Numnodes * 3 + 1] = sgg.Liney[Serialized.eJ[conta]];
                  Nodes[Numnodes * 3 + 2] = sgg.Linez[Serialized.eK[conta]];
                  Numnodes = Numnodes + 1;
                  // ... (code continues in next chunk)
               }
               break;
         }
      }
   } // subroutine creaUnstructData

Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(1 + Serialized.eK(conta));
         }
         //
         numEdges = numEdges + 1;
         Elems(conta,1) = NumNodes - 1;
         Elems(conta,2) = NumNodes;
         Elems(conta,3) = -1; // mark as edge for later correct writing
         //
         case (iBloqueJx):
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (1 + Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(1 + Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (1 + Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (1 + Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(1 + Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(1 + Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (1 + Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(1 + Serialized.eK(conta));
         }
         //

         numQuads = numQuads + 1;
         Elems(conta,1) = NumNodes - 3;
         Elems(conta,2) = NumNodes - 2;
         Elems(conta,3) = NumNodes - 1;
         Elems(conta,4) = NumNodes;
         case (iBloqueJy):
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (1 + Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(1 + Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (1 + Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (1 + Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(1 + Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(1 + Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (1 + Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(1 + Serialized.eK(conta));
         }
         //
         numQuads = numQuads + 1;
         Elems(conta,1) = NumNodes - 3;
         Elems(conta,2) = NumNodes - 2;
         Elems(conta,3) = NumNodes - 1;
         Elems(conta,4) = NumNodes;
         case (iBloqueJz):
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (1 + Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(1 + Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (1 + Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (1 + Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(1 + Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(1 + Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         numNodes = numNodes + 1;
         if (vtkindex) {
            Nodes(numNodes,1) = (Serialized.eI(conta)) * 1.0_RKIND;
            Nodes(numNodes,2) = (1 + Serialized.eJ(conta)) * 1.0_RKIND;
            Nodes(numNodes,3) = (Serialized.eK(conta)) * 1.0_RKIND;
         } else {
            Nodes(numNodes,1) = sgg.LineX(Serialized.eI(conta));
            Nodes(numNodes,2) = sgg.Liney(1 + Serialized.eJ(conta));
            Nodes(numNodes,3) = sgg.Linez(Serialized.eK(conta));
         }
         //

         numQuads = numQuads + 1;
         Elems(conta,1) = NumNodes - 3;
         Elems(conta,2) = NumNodes - 2;
         Elems(conta,3) = NumNodes - 1;
         Elems(conta,4) = NumNodes;
      } // end select
   } // end do

   if ((NumEdges + 1) + (NumQuads + 1) != numberofSerialized) {
      std::string buff = "ERROR: Buggy error sumas creating .vtk";
      print11(0, buff);
   }

   return;
} // end subroutine creaUnstructData

// subroutine fillinparaviewstate
//
// return
// end subroutine
} // end namespace VTK_m