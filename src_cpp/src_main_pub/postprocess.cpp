#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <sstream>

// Assuming these types and functions are defined in other modules
// FDETYPES_m
using RKIND = double;
using RKIND_tiempo = double;
using CKIND = std::complex<double>;
using integer = int;

// Report_m
void print11(int layoutnumber, const std::string& msg) {
    std::cout << "[Process " << layoutnumber << "] " << msg << std::endl;
}

// Observa_m
// Forward declarations for types and functions used
struct SGGFDTDINFO_t {
    int NumberRequest;
    struct {
        bool Volumic;
        bool FreqDomain;
        bool Transfer;
        bool TimeDomain;
        RKIND InitialFreq;
        RKIND FinalFreq;
        RKIND FreqStep;
        RKIND TimeStep;
        std::string FileNormalize;
        int dt; // Assuming dt is accessible here or via a method
    } Observation[100]; // Assuming max size or dynamic, simplified for translation
    RKIND dt;
};

struct output_item_t {
    std::string path;
    int unit;
    int columnas;
    int MPIRoot;
    int Trancos;
};

struct output_t {
    std::vector<output_item_t> item;
};

output_t GetOutput(); // Placeholder for actual implementation

// Helper function to simulate Fortran's index function for strings
int index(const std::string& str, const std::string& sub) {
    size_t pos = str.find(sub);
    if (pos != std::string::npos) {
        return static_cast<int>(pos + 1); // Fortran is 1-based
    }
    return 0;
}

// Placeholder for dtft subroutine
void dtft(std::vector<std::complex<double>>& fqValues, 
          const std::vector<double>& fqPos, 
          int fqLength, 
          const std::vector<double>& samplingTime, 
          const std::vector<double>& signal, 
          int timesteps) {
    // Simplified DTFT calculation placeholder
    fqValues.assign(fqLength, std::complex<double>(0.0, 0.0));
    for (int i = 0; i < fqLength; ++i) {
        for (int n = 0; n < timesteps; ++n) {
            double angle = -2.0 * M_PI * fqPos[i] * samplingTime[n];
            fqValues[i] += signal[n] * std::exp(std::complex<double>(0, angle));
        }
    }
}

// Placeholder for conviertecabecera
void conviertecabecera(const std::string& cabecera, std::string& cabeceraNew, int numComp, RKIND rinstant) {
    cabeceraNew = "Time";
    for(int i=0; i<numComp-1; ++i) {
        cabeceraNew += " Value" + std::to_string(i+1);
    }
    cabeceraNew += "\n";
}

namespace PostProcessing_m {

    void PostProcess(int layoutnumber, int num_procs, const SGGFDTDINFO_t& sgg, const std::string& nEntradaRoot, 
                     RKIND rinstant, bool& somethingdone, bool niapapostprocess, bool forceresampled) {

        output_t output = GetOutput(); 
        
        somethingdone = false;
        
        std::string whoamishort;
        std::ostringstream oss_short;
        oss_short << std::setw(5) << std::setfill(' ') << (layoutnumber + 1);
        whoamishort = oss_short.str();

        std::string whoami;
        std::ostringstream oss_whoami;
        oss_whoami << "(" << (layoutnumber + 1) << "/" << num_procs << ") ";
        whoami = oss_whoami.str();

        std::string filename = nEntradaRoot + "_Outputrequests_" + whoamishort + ".txt";
        std::ifstream infile(filename);
        if (!infile.is_open()) {
            print11(layoutnumber, "Could not open output requests file");
            return;
        }

        std::string cabecera;
        while (std::getline(infile, cabecera)) {
            // Trim whitespace
            size_t start = cabecera.find_first_not_of(" \t\r\n");
            size_t end = cabecera.find_last_not_of(" \t\r\n");
            if (start != std::string::npos) {
                cabecera = cabecera.substr(start, end - start + 1);
            } else {
                cabecera = "";
            }
            
            if (cabecera == "!END") {
                break;
            }
        }
        // In Fortran, backspace moves back one record. In stream/text files, this is tricky.
        // We assume the loop stopped at !END, so we need to re-read or handle state.
        // For simplicity in C++ stream processing, we might need to reopen or manage state.
        // Here we assume the file pointer is conceptually at the start of the data block.
        // Re-opening to simulate backspace behavior for line-based files is safer.
        infile.close();
        infile.open(filename);
        if (!infile.is_open()) return;

        // Skip header until !END
        while (std::getline(infile, cabecera)) {
             size_t start = cabecera.find_first_not_of(" \t\r\n");
             size_t end = cabecera.find_last_not_of(" \t\r\n");
             if (start != std::string::npos) {
                 cabecera = cabecera.substr(start, end - start + 1);
             } else {
                 cabecera = "";
             }
             if (cabecera == "!END") {
                 break;
             }
        }

        // Now read the actual data blocks
        // Note: The original Fortran code structure implies a specific format for the .txt file
        // which is not fully detailed here. We will iterate based on sgg.NumberRequest.
        
        for (int ii = 0; ii < sgg.NumberRequest; ++ii) {
            // Assuming sgg.Observation is an array/vector. 
            // In C++, we need to access it properly. Let's assume it's a vector or array.
            // Since the struct definition above was simplified, we access directly.
            // Note: Fortran arrays are 1-based, C++ 0-based. Adjusting indices if necessary.
            // The loop `do i=1,sgg%Observation(ii)%nP` suggests Np items per observation.
            // We need access to Np. Let's assume it's part of the structure.
            
            // Placeholder for Np access. In real code, this would be sgg.Observation[ii].nP
            int Np = 1; // Placeholder
            for (int i = 0; i < Np; ++i) {
                int field = 0; // Placeholder for sgg.Observation[ii].P[i].What
                // if (field != nothing) { ... }
                
                // Simplified logic flow based on the provided snippet
                // The snippet cuts off, so we translate what is visible.
                
                bool escribirBloque = true;
                bool escribir = true;
                
                if (escribir) {
                    if (sgg.Observation[ii].FreqDomain) {
                        int pasadas = 1;
                        int pp = 1;
                        std::string path = output.item[i].path;
                        
                        bool existe = false;
                        // Check file existence
                        std::ifstream test_file(path);
                        existe = test_file.good();
                        test_file.close();

                        if (!existe) {
                            std::string buff = "Not processing: Inexistent file " + path;
                            print11(layoutnumber, buff);
                        } else {
                            int numComp = output.item[i].columnas;
                            bool neverprecounted = true;
                            
                            std::vector<double> valores;
                            std::vector<double> tiempo;
                            std::vector<double> signal;
                            std::vector<double> samplingtime;
                            
                            if (neverprecounted) {
                                neverprecounted = false;
                                std::ifstream datafile(path);
                                if (!datafile.is_open()) {
                                    print11(layoutnumber, "Could not open data file: " + path);
                                    continue;
                                }
                                
                                std::string header_line;
                                std::getline(datafile, header_line);
                                
                                int timesteps = 0;
                                double dummy;
                                while (datafile >> dummy >> dummy) {
                                    timesteps++;
                                }
                                datafile.close();
                                timesteps--; // Adjust for last step issue mentioned in comment
                                
                                if (timesteps <= 0) {
                                    print11(layoutnumber, "No timesteps found in file: " + path);
                                    continue;
                                }

                                valores.resize(timesteps * numComp);
                                tiempo.resize(timesteps);
                                signal.resize(timesteps);
                                samplingtime.resize(timesteps);
                                
                                // Read data
                                std::ifstream datafile2(path);
                                std::getline(datafile2, header_line);
                                for (int ns = 0; ns < timesteps; ++ns) {
                                    datafile2 >> tiempo[ns];
                                    for (int compo = 0; compo < numComp; ++compo) {
                                        datafile2 >> valores[ns * numComp + compo];
                                    }
                                }
                                datafile2.close();
                            }
                            
                            if (niapapostprocess) {
                                std::cout << "Correcting in FreqDomain postprocess " << timesteps << " " << path << std::endl;
                                for (int ns = 0; ns < timesteps; ++ns) {
                                    tiempo[ns] = static_cast<RKIND_tiempo>(ns * sgg.dt);
                                }
                            }
                            
                            if (forceresampled) {
                                std::string path_resampled = path.substr(0, index(path, ".dat") - 1) + "_resampled_time.dat";
                                int columna = 1; // 1-based index in Fortran, mapped to 0-based logic if needed, but here used as column index
                                // Fortran column index 1 corresponds to C++ index 0
                                int col_idx = columna - 1;
                                
                                std::ofstream resampled_file(path_resampled);
                                if (!resampled_file.is_open()) {
                                    print11(layoutnumber, "Could not open resampled file: " + path_resampled);
                                    continue;
                                }
                                
                                RKIND t_pedido = tiempo[0];
                                resampled_file << std::setprecision(15) << t_pedido << " " << valores[col_idx] << std::endl;
                                
                                for (int iii = 1; iii < timesteps; ++iii) {
                                    while (t_pedido <= tiempo[iii]) {
                                        t_pedido += sgg.Observation[ii].TimeStep;
                                        bool found = false;
                                        for (int jjj = iii - 1; jjj < timesteps - 1; ++jjj) {
                                            if (t_pedido >= tiempo[jjj] && t_pedido < tiempo[jjj + 1]) {
                                                RKIND value_interp = (valores[(jjj + 1) * numComp + col_idx] - valores[jjj * numComp + col_idx]) / 
                                                                     (tiempo[jjj + 1] - tiempo[jjj]) * (t_pedido - tiempo[jjj]) + valores[jjj * numComp + col_idx];
                                                resampled_file << std::setprecision(15) << t_pedido << " " << value_interp << std::endl;
                                                found = true;
                                                break;
                                            }
                                        }
                                        if (!found) break;
                                    }
                                }
                                resampled_file.close();
                            }
                            
                            RKIND fmin = std::min(sgg.Observation[ii].FinalFreq, sgg.Observation[ii].InitialFreq);
                            RKIND fmax = std::max(sgg.Observation[ii].FinalFreq, sgg.Observation[ii].InitialFreq);
                            
                            RKIND fstep;
                            if (sgg.Observation[ii].FreqStep == 0.0 || sgg.Observation[ii].FreqStep > fmax - fmin) {
                                fstep = fmax - fmin;
                            } else {
                                fstep = sgg.Observation[ii].FreqStep;
                            }
                            
                            int fqLength = static_cast<int>((fmax - fmin) / fstep) + 2;
                            
                            std::vector<double> fqPos(fqLength);
                            std::vector<std::complex<double>> fqValues(fqLength);
                            std::vector<std::complex<double>> valoresDF(fqLength * numComp);
                            
                            int pozi = index(path, "_log_");
                            if (pozi != 0) {
                                fmin = std::max(1.0, std::log10(fmin));
                                fmax = std::log10(fmax);
                                fstep = (fmax - fmin) / (fqLength - 2.0);
                            }
                            
                            for (int i1 = 0; i1 < fqLength; ++i1) {
                                fqPos[i1] = fmin + (i1) * fstep;
                            }
                            
                            if (pozi != 0) {
                                for (int i1 = 0; i1 < fqLength; ++i1) {
                                    fqPos[i1] = std::pow(10.0, fqPos[i1]);
                                }
                            }
                            
                            signal.assign(timesteps, 0.0);
                            samplingtime.assign(timesteps, 0.0);
                            samplingtime = tiempo;
                            
                            for (int i1 = 0; i1 < numComp; ++i1) {
                                for (int ns = 0; ns < timesteps; ++ns) {
                                    signal[ns] = valores[ns * numComp + i1];
                                }
                                dtft(fqValues, fqPos, fqLength, samplingtime, signal, timesteps);
                                for (int k = 0; k < fqLength; ++k) {
                                    valoresDF[k * numComp + i1] = fqValues[k];
                                }
                            }
                            
                            std::string path2 = path.substr(0, index(path, ".dat") - 1) + "_df.dat";
                            std::ofstream df_file(path2);
                            if (!df_file.is_open()) {
                                print11(layoutnumber, "Could not open DF file: " + path2);
                            } else {
                                std::string cabeceraNew;
                                conviertecabecera(header_line, cabeceraNew, numComp + 1, rinstant);
                                df_file << cabeceraNew;
                                
                                for (int i1 = 0; i1 < fqLength; ++i1) {
                                    df_file << std::setprecision(15) << fqPos[i1];
                                    for (int j1 = 0; j1 < numComp; ++j1) {
                                        RKIND abs_val = std::abs(valoresDF[i1 * numComp + j1]);
                                        RKIND phase = std::atan2(std::imag(valoresDF[i1 * numComp + j1]), std::real(valoresDF[i1 * numComp + j1]));
                                        df_file << " " << std::setprecision(15) << abs_val << " " << std::setprecision(15) << phase;
                                    }
                                    df_file << std::endl;
                                }
                                df_file.close();
                            }
                            
                            // Update output list file
                            std::ofstream list_file(filename, std::ios::app);
                            if (list_file.is_open()) {
                                list_file << path2 << std::endl;
                                list_file.close();
                            }
                            
                            // Cleanup
                            // In C++, vectors are automatically cleaned up when going out of scope or reassigned
                        }
                    }
                }
                somethingdone = true;
            }
        }
    }

    void postprocessonthefly() {
        // Placeholder as the full implementation wasn't provided in the snippet
    }

} // namespace PostProcessing_m

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

// Assuming necessary types and functions are defined in headers included previously
// e.g., SGGFDTDINFO_t, output_t, GetOutput(), print11(), BUFSIZE, RKIND, RKIND_tiempo, fmt, etc.
// Also assuming constants like iBloqueJx, iBloqueJy, iBloqueMx, iBloqueMy, CompileWithMPI

// Helper function to simulate Fortran's trim(adjustl(str))
inline std::string trim_adjustl(const std::string& str) {
    size_t start = str.find_first_not_of(' ');
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(' ');
    return str.substr(start, end - start + 1);
}

// Helper function to simulate Fortran's write to stream with a format string
// Note: The original code uses a specific 'fmt' variable which is likely a character string format.
// For simplicity in translation, we assume a standard stream output or a helper function.
// Since 'fmt' is used in write statements like `write(unit, fmt) val1, val2`, 
// we will assume a generic print function or direct stream insertion if the format is simple.
// However, to be safe and preserve behavior, we'll use a placeholder or standard cout/cerr logic 
// if 'fmt' is not explicitly passed. Given the context, it's likely a specific format string.
// Let's assume a helper function `write_formatted` exists or use standard stream formatting.
// For this translation, I will use standard stream insertion for simplicity, 
// assuming the 'fmt' handles basic numeric output. If 'fmt' is complex, a parser would be needed.
// Given the constraints, I'll use a simple stream write.

void write_formatted(std::ostream& os, const std::string& fmt, double val1, double val2) {
    // Simplified: Just write the values. In a real scenario, 'fmt' would be parsed.
    // Fortran default format for reals is often G12.5 or similar.
    os << std::scientific << std::setprecision(6) << val1 << " " << val2 << "\n";
}

// Function: almostequal
bool almostequal(double a, double b) {
    double ratio;
    const double tolerancia = 0.01;

    if (b == 0.0) {
        // Avoid division by zero, though Fortran might handle it differently or assume non-zero
        return false; 
    }
    
    ratio = a / b;
    if (std::abs(ratio) < 1.0) {
        ratio = 1.0 / ratio;
    }
    
    if ((ratio > 0.0) && (ratio < 1.0 + tolerancia)) {
        return true;
    } else {
        return false;
    }
}

// Subroutine: conviertecabecera
void conviertecabecera(int columnas, std::string& c, std::string& cNew, double rinstant) {
    const int BUFSIZE = 256; // Assumed size
    std::string chninstant(BUFSIZE, ' ');
    std::string c2(columnas, ' '); // Vector of strings, initialized to spaces
    int i, j, k, longi, ii;
    
    // Simulate Fortran write to chninstant with fmt
    // Assuming fmt outputs rinstant into chninstant
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << rinstant;
    chninstant = oss.str();
    // Trim chninstant to remove trailing spaces if any, though Fortran char is fixed length
    // Fortran write to char usually pads or truncates. Let's assume it fits.
    
    // Initialize c2 with spaces
    for (int idx = 0; idx < columnas; ++idx) {
        c2[idx] = " ";
    }

    longi = static_cast<int>(trim_adjustl(c).length()) + 1;
    k = 1;
    i = 1;
    
    // Fortran do while loop
    while (i <= longi) {
        if (c[i-1] != ' ') { // Fortran is 1-indexed, C++ is 0-indexed
            // interno do j=i,longi
            for (j = i; j <= longi; ++j) {
                if (c[j-1] == ' ') {
                    if (k == 1) {
                        c2[k-1] = " f_at_" + trim_adjustl(chninstant);
                    } else {
                        // c(i:j-1) corresponds to substring from i to j-1
                        // In C++, substring is (start, length)
                        // start index is i-1, length is (j-1) - (i-1) = j - i
                        std::string sub = c.substr(i-1, j - i);
                        c2[k-1] = sub + "_Mod   " + sub + "_Phase   ";
                    }
                    k++;
                    if (k > columnas) {
                        goto end_buscam;
                    }
                    // do ii=j,longi
                    for (ii = j; ii <= longi; ++ii) {
                        if ((c[ii-1] != ' ') || (ii == longi)) {
                            i = ii;
                            goto end_interno;
                        }
                    }
                }
            }
            end_interno:;
        } else {
            i++;
        }
    }
    
    end_buscam:;

    cNew = " ";
    for (i = 1; i <= columnas; ++i) {
        cNew = trim_adjustl(cNew) + "   " + trim_adjustl(c2[i-1]);
    }
}

// Subroutine: postprocessonthefly
void postprocessonthefly(int layoutnumber, int num_procs, const SGGFDTDINFO_t& sgg, const std::string& nEntradaRoot, double rinstant, bool& somethingdone, bool niapapostprocess, bool forceresampled) {
    std::vector<std::unique_ptr<output_t>> output; // Assuming GetOutput returns a vector or similar
    // Note: The original code uses `Output => GetOutput()`, implying a pointer assignment.
    // In C++, we might get a reference or a pointer. Let's assume GetOutput returns a pointer or reference.
    // For translation safety, let's assume GetOutput() returns a pointer to a vector of output_t or similar structure.
    // However, without the definition of GetOutput, we'll assume it returns a pointer to the first element of a vector,
    // or we manage the vector ourselves.
    // Let's assume `output` is a vector of pointers to output_t objects.
    
    // Re-reading: `type(output_t), pointer, dimension( : ) :: output`
    // `Output => GetOutput()`
    // This suggests GetOutput() returns a pointer to an array of output_t.
    // We will assume GetOutput() returns `std::vector<output_t>*` or similar.
    // To keep it simple and safe, let's assume GetOutput() returns a pointer to the first element of a dynamically allocated array,
    // or we use a vector. Given the complexity, I'll assume GetOutput() returns a vector of output_t by value or reference.
    // Let's assume `output` is a `std::vector<output_t>` obtained from GetOutput.
    
    // Actually, looking at the usage `output(ii)%item(i)%unit`, it's an array of structs.
    // Let's assume GetOutput() returns a `std::vector<output_t>*`.
    
    std::vector<output_t>* output_ptr = GetOutput();
    if (!output_ptr) {
        return; // Handle error
    }
    std::vector<output_t>& output = *output_ptr;

    std::string path;
    bool existe, escribir, escribirBloque;
    int ii, i, field;
    std::string buff;

    for (ii = 1; ii <= sgg.NumberRequest; ++ii) {
        for (i = 1; i <= sgg.Observation[ii-1].nP; ++i) {
            field = sgg.Observation[ii-1].P[i-1].What;
            
#ifdef CompileWithMPI
            escribirBloque = ((field == iBloqueJx) || (field == iBloqueJy) || 
                              (field == iBloqueMx) || (field == iBloqueMy)) && 
                             (layoutnumber == output[ii-1].item[i-1].MPIRoot);
#else
            escribirBloque = true;
#endif
            
            escribir = (((field != iBloqueJx) && (field != iBloqueJy) && 
                         (field != iBloqueMx) && (field != iBloqueMy)) || 
                         escribirBloque) && 
                        (sgg.Observation[ii-1].FreqDomain || sgg.Observation[ii-1].Transfer);
            
            if (escribir) {
                if (sgg.Observation[ii-1].FreqDomain) {
                    path = trim_adjustl(output[ii-1].item[i-1].path);
                    // inquire(file=..., exist=existe)
                    std::ifstream test_file(path);
                    existe = test_file.good();
                    test_file.close();
                    
                    if (!existe) {
                        buff = "Not processing: Inexistent file " + path;
                        print11(layoutnumber, buff);
                    } else {
                        // close (output(ii)%item(i)%unit)
                        // Assuming output[ii-1].item[i-1].unit is a file stream
                        if (output[ii-1].item[i-1].unit.is_open()) {
                            output[ii-1].item[i-1].unit.close();
                        }
                    }
                }
            }
        }
    }
    
    // Call postprocess
    postprocess(layoutnumber, num_procs, sgg, nEntradaRoot, rinstant, somethingdone, niapapostprocess, forceresampled);
    
    for (ii = 1; ii <= sgg.NumberRequest; ++ii) {
        for (i = 1; i <= sgg.Observation[ii-1].nP; ++i) {
            field = sgg.Observation[ii-1].P[i-1].What;
            
#ifdef CompileWithMPI
            escribirBloque = ((field == iBloqueJx) || (field == iBloqueJy) || 
                              (field == iBloqueMx) || (field == iBloqueMy)) && 
                             (layoutnumber == output[ii-1].item[i-1].MPIRoot);
#else
            escribirBloque = true;
#endif
            
            escribir = (((field != iBloqueJx) && (field != iBloqueJy) && 
                         (field != iBloqueMx) && (field != iBloqueMy)) || 
                         escribirBloque) && 
                        (sgg.Observation[ii-1].FreqDomain || sgg.Observation[ii-1].Transfer);
            
            if (escribir) {
                if (sgg.Observation[ii-1].FreqDomain) {
                    path = trim_adjustl(output[ii-1].item[i-1].path);
                    // inquire(file=..., exist=existe)
                    std::ifstream test_file(path);
                    existe = test_file.good();
                    test_file.close();
                    
                    if (!existe) {
                        buff = "ERROR: Inexistent file. Creating new " + path;
                        print11(layoutnumber, buff);
                    } else {
                        // open (output(ii)%item(i)%unit, file=..., form='formatted', position='append')
                        if (output[ii-1].item[i-1].unit.is_open()) {
                            output[ii-1].item[i-1].unit.close();
                        }
                        output[ii-1].item[i-1].unit.open(path, std::ios::out | std::ios::app | std::ios::trunc); // 'append' usually means ios::app
                        // Note: Fortran 'append' keeps existing content and adds to end.
                        // C++ ios::app does this.
                        // 'formatted' is default for text files.
                    }
                }
            }
        }
    }
}