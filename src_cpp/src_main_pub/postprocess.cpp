#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <iomanip>

// Forward declarations for external modules/types
// These would typically be in their respective headers
namespace FDETYPES_m {
    using RKIND = double;
    using RKIND_tiempo = double;
    using CKIND = std::complex<double>;
    const int BUFSIZE = 256;
}

namespace Report_m {
    void print11(int layoutnumber, const std::string& buff);
    void conviertecabecera(const std::string& cabecera, std::string& cabeceraNew, int numComp, double rinstant);
}

namespace Observa_m {
    // Placeholder for SGGFDTDINFO_t and output_t structures
    // These need to be defined based on Observa_m.f90 content
    
    enum FieldTypes {
        nothing = 0,
        iBloqueJx = 1,
        iBloqueJy = 2,
        iBloqueMx = 3,
        iBloqueMy = 4
    };

    struct Item_t {
        std::string path;
        int unit;
        int columnas;
        int MPIRoot;
        int Trancos;
    };

    struct Observation_t {
        int nP;
        bool Volumic;
        bool FreqDomain;
        bool Transfer;
        double InitialFreq;
        double FinalFreq;
        double FreqStep;
        double TimeStep;
        std::vector<int> P; // Indices or pointers to P structures
        std::string FileNormalize;
    };

    struct P_t {
        int What;
    };

    struct SGGFDTDINFO_t {
        int NumberRequest;
        double dt;
        std::vector<Observation_t> Observation;
    };

    struct output_t {
        std::vector<Item_t> item;
        int Trancos;
    };

    output_t GetOutput();
    
    // External function declaration for DTFT
    void dtft(std::vector<std::complex<double>>& fqValues, 
              std::vector<double>& fqPos, 
              int fqLength, 
              const std::vector<double>& samplingTime, 
              const std::vector<double>& signal, 
              int timesteps);
}

namespace PostProcessing_m {

    using namespace FDETYPES_m;
    using namespace Report_m;
    using namespace Observa_m;

    // Helper to mimic Fortran's trim(adjustl(...))
    std::string trim_adjustl(const std::string& str) {
        size_t first = str.find_first_not_of(' ');
        if (std::string::npos == first) {
            return str;
        }
        size_t last = str.find_last_not_of(' ');
        return str.substr(first, (last - first + 1));
    }

    // Format string for output, mimicking Fortran default formatting
    // Note: Fortran default formatting can vary. Using a generic scientific or fixed format.
    // The original code uses 'fmt' which is likely defined elsewhere or implicit.
    // Assuming a standard double precision format.
    const std::string fmt = "%.15e"; 

    void PostProcess(int layoutnumber, int num_procs, const SGGFDTDINFO_t& sgg, 
                     const std::string& nEntradaRoot, double rinstant, 
                     bool& somethingdone, bool niapapostprocess, bool forceresampled) {

        output_t* output = new output_t(); // Mimicking pointer assignment from GetOutput()
        *output = GetOutput();

        std::string whoami, whoamishort;
        std::string cabecera, cabeceraNew, path, path2, path3, path_resampled;
        int numComp;
        
        // Allocatable arrays
        std::vector<double> valores;
        std::vector<std::complex<double>> valoresDF, valoresDF2;
        std::vector<double> tiempo, samplingtime;
        
        int pozi;
        bool existe, neverprecounted, escribir, escribirBloque;
        int fqLength, ii, i, i1, j1, field, ns, timesteps, compo, iii, pp, pasadas;
        double dummy;
        
        std::vector<double> fqPos, signal;
        double fmin, fmax, fstep, value_interp;
        std::vector<std::complex<double>> fqValues;
        
        std::string buff, dubuf;
        int columna, jjj;
        int my_iostat = 0; // Placeholder for iostat

        if (niapapostprocess) {
            std::cout << "Copiar a mano los .dat en tiempo para que se postrocesen bien..." << std::endl;
            // pause is not directly available in C++, usually handled by debugger or exit
            std::cout << "Continuing..." << std::endl;
        }

        // Format whoamishort: (i5)
        std::ostringstream oss_short;
        oss_short << std::setw(5) << std::setfill(' ') << (layoutnumber + 1);
        whoamishort = oss_short.str();

        // Format whoami: (a,i5,a,i5,a)
        std::ostringstream oss_whoami;
        oss_whoami << "(" << std::setw(5) << std::setfill(' ') << (layoutnumber + 1) 
                   << "/" << std::setw(5) << std::setfill(' ') << num_procs << ") ";
        whoami = oss_whoami.str();

        // Open output request list file
        std::string output_requests_file = trim_adjustl(nEntradaRoot) + "_Outputrequests_" + trim_adjustl(whoamishort) + ".txt";
        std::ifstream file119(output_requests_file);
        if (!file119.is_open()) {
            std::cerr << "Error opening file: " << output_requests_file << std::endl;
            return;
        }

        cabecera = " ";
        while (trim_adjustl(cabecera) != "!END") {
            std::getline(file119, cabecera);
        }
        // backspace(119) in Fortran moves the file position back one record.
        // In text files, this is tricky. We'll assume the next read will get the last line read.
        // However, since we just read "!END", we need to re-read the previous line.
        // A simpler approach in C++ for this specific logic:
        // The loop reads until "!END". The line before "!END" is the first record of interest.
        // We need to rewind or re-read.
        file119.clear();
        file119.seekg(0, std::ios::beg);
        std::string line;
        std::vector<std::string> records;
        while (std::getline(file119, line)) {
            if (trim_adjustl(line) == "!END") break;
            records.push_back(line);
        }
        file119.close();

        // Process each request
        for (ii = 0; ii < sgg.NumberRequest; ++ii) {
            for (i = 0; i < sgg.Observation[ii].nP; ++i) {
                // Accessing P(i). What. Assuming P is a vector of P_t or similar.
                // The Fortran code uses sgg%Observation(ii)%P(i)%What
                // We need to ensure P is accessible.
                if (i >= sgg.Observation[ii].P.size()) continue; 
                
                field = sgg.Observation[ii].P[i].What;
                
                if (field != nothing) {
                    if (!sgg.Observation[ii].Volumic) {
                        
                        // MPI Logic
#ifdef CompileWithMPI
                        bool isBlockField = (field == iBloqueJx) || (field == iBloqueJy) || 
                                            (field == iBloqueMx) || (field == iBloqueMy);
                        escribirBloque = isBlockField && (layoutnumber == output->item[i].MPIRoot);
#else
                        escribirBloque = true;
#endif

                        bool isNonBlockField = (field != iBloqueJx) && (field != iBloqueJy) && 
                                               (field != iBloqueMx) && (field != iBloqueMy);
                        
                        escribir = ((isNonBlockField || escribirBloque) && 
                                    (sgg.Observation[ii].FreqDomain || sgg.Observation[ii].Transfer));

                        if (escribir) {
                            if (sgg.Observation[ii].FreqDomain) {
                                pasadas = 1;
                                pp = 1;
                                if (pp == 1) {
                                    path = trim_adjustl(output->item[i].path);
                                }

                                // Check if file exists
                                std::ifstream test_file(path);
                                existe = test_file.good();
                                test_file.close();

                                if (!existe) {
                                    buff = "Not processing: Inexistent file " + trim_adjustl(path);
                                    print11(layoutnumber, buff);
                                } else {
                                    numComp = output->item[i].columnas;
                                    neverprecounted = true;

                                    if (neverprecounted) {
                                        neverprecounted = false;
                                        std::ifstream infile(path);
                                        if (!infile.is_open()) {
                                            std::cerr << "Error opening file for precounting: " << path << std::endl;
                                            continue;
                                        }
                                        
                                        std::string header_line;
                                        std::getline(infile, header_line); // Read header
                                        
                                        timesteps = 0;
                                        double d1, d2;
                                        while (infile >> d1 >> d2) {
                                            timesteps++;
                                        }
                                        infile.close();
                                        
                                        // Fortran logic: timesteps = timesteps - 1
                                        timesteps--; 
                                        
                                        if (valores.size() > 0) {
                                            valores.clear();
                                            tiempo.clear();
                                            signal.clear();
                                        }
                                        
                                        valores.resize(timesteps, 0.0);
                                        tiempo.resize(timesteps, 0.0);
                                        signal.resize(timesteps, 0.0);
                                        samplingtime.resize(timesteps, 0.0);
                                    }

                                    // Read data
                                    std::ifstream infile(path);
                                    if (!infile.is_open()) {
                                        std::cerr << "Error opening file for reading: " << path << std::endl;
                                        continue;
                                    }
                                    
                                    std::string header_line;
                                    std::getline(infile, header_line); // Read header
                                    
                                    for (ns = 0; ns < timesteps; ++ns) {
                                        infile >> tiempo[ns];
                                        for (compo = 0; compo < numComp; ++compo) {
                                            infile >> valores[ns * numComp + compo];
                                        }
                                    }
                                    infile.close();

                                    if (niapapostprocess) {
                                        std::cout << "Correcting in FreqDomain postprocess " << timesteps << " " << trim_adjustl(path) << std::endl;
                                        for (ns = 0; ns < timesteps; ++ns) {
                                            tiempo[ns] = static_cast<double>(ns + 1) * sgg.dt; // Fortran ns starts at 1, C++ at 0. 
                                                                                            // Fortran: tiempo(ns)=real(ns*sgg%dt). If ns=1..timesteps.
                                                                                            // Here ns is 0..timesteps-1. So index is ns+1.
                                        }
                                    }

                                    if (forceresampled) {
                                        // Resample logic
                                        size_t dot_pos = path.find(".dat");
                                        std::string base_path = path.substr(0, dot_pos);
                                        path_resampled = base_path + "_resampled_time.dat";
                                        
                                        columna = 1; // Fortran 1-based index, used as column index
                                        std::ofstream outfile_resampled(path_resampled);
                                        if (!outfile_resampled.is_open()) {
                                            std::cerr << "Error opening resampled file: " << path_resampled << std::endl;
                                        } else {
                                            double t_pedido = tiempo[0];
                                            outfile_resampled << std::fixed << std::setprecision(15) << t_pedido << " " << valores[0] << std::endl;
                                            
                                            for (iii = 1; iii < timesteps; ++iii) {
                                                while (t_pedido <= tiempo[iii]) {
                                                    t_pedido += sgg.Observation[ii].TimeStep;
                                                    
                                                    bool found = false;
                                                    for (jjj = iii - 1; jjj < timesteps - 1; ++jjj) {
                                                        if (t_pedido >= tiempo[jjj] && t_pedido < tiempo[jjj + 1]) {
                                                            value_interp = (valores[(jjj + 1) * numComp + (columna - 1)] - 
                                                                            valores[jjj * numComp + (columna - 1)]) / 
                                                                           (tiempo[jjj + 1] - tiempo[jjj]) * 
                                                                           (t_pedido - tiempo[jjj]) + 
                                                                           valores[jjj * numComp + (columna - 1)];
                                                            outfile_resampled << std::fixed << std::setprecision(15) << t_pedido << " " << value_interp << std::endl;
                                                            found = true;
                                                            break;
                                                        }
                                                    }
                                                    if (!found) break; // Should not happen if logic is correct
                                                }
                                            }
                                            outfile_resampled.close();
                                        }
                                    }

                                    // DFT Calculation
                                    fmin = std::min(sgg.Observation[ii].FinalFreq, sgg.Observation[ii].InitialFreq);
                                    fmax = std::max(sgg.Observation[ii].FinalFreq, sgg.Observation[ii].InitialFreq);

                                    if (sgg.Observation[ii].FreqStep == 0.0 || sgg.Observation[ii].FreqStep > (fmax - fmin)) {
                                        fstep = fmax - fmin;
                                    } else {
                                        fstep = sgg.Observation[ii].FreqStep;
                                    }

                                    fqLength = static_cast<int>((fmax - fmin) / fstep) + 2;
                                    
                                    fqPos.resize(fqLength, 0.0);
                                    fqValues.resize(fqLength, std::complex<double>(0.0, 0.0));
                                    valoresDF.resize(fqLength * numComp, 0.0); // Flattened or 2D? Fortran is (1:fqLength, 1:numComp)
                                    // Using vector of vectors for 2D logic or flattened. Let's use flattened for simplicity in C++ vector
                                    // But Fortran access is valoresDF(i1, j1).
                                    // Let's use std::vector<std::vector<double>> for clarity if numComp is small, or flattened.
                                    // Given the complexity, let's stick to flattened: index = (i1-1)*numComp + (j1-1)
                                    
                                    // Re-allocating specifically for complex values if needed, but Fortran uses valoresDF as real?
                                    // Fortran: complex(kind=CKIND), allocatable, dimension(:,:) :: valoresDF
                                    // So valoresDF is complex.
                                    std::vector<std::complex<double>> valoresDF_complex(fqLength * numComp, std::complex<double>(0.0, 0.0));

                                    pozi = path.find("_log_");
                                    if (pozi != std::string::npos) {
                                        fmin = std::max(1.0, std::log10(fmin));
                                        fmax = std::log10(fmax);
                                        fstep = (fmax - fmin) / (fqLength - 2.0);
                                    }

                                    for (i1 = 0; i1 < fqLength; ++i1) {
                                        fqPos[i1] = fmin + i1 * fstep;
                                    }

                                    if (pozi != std::string::npos) {
                                        for (i1 = 0; i1 < fqLength; ++i1) {
                                            fqPos[i1] = std::pow(10.0, fqPos[i1]);
                                        }
                                    }

                                    signal.assign(timesteps, 0.0);
                                    samplingtime.assign(timesteps, 0.0);
                                    samplingtime = tiempo;

                                    for (i1 = 0; i1 < numComp; ++i1) {
                                        for (ns = 0; ns < timesteps; ++ns) {
                                            signal[ns] = valores[ns * numComp + i1];
                                        }
                                        
                                        std::vector<std::complex<double>> fqValues_local(fqLength, std::complex<double>(0.0, 0.0));
                                        dtft(fqValues_local, fqPos, fqLength, samplingtime, signal, timesteps);
                                        
                                        for (int k = 0; k < fqLength; ++k) {
                                            valoresDF_complex[k * numComp + i1] = fqValues_local[k];
                                        }
                                    }

                                    path2 = path.substr(0, path.find(".dat")) + "_df.dat";
                                    std::ofstream outfile_df(path2);
                                    if (!outfile_df.is_open()) {
                                        std::cerr << "Error opening DF file: " << path2 << std::endl;
                                    } else {
                                        // Update output list file
                                        std::ofstream file119_update(output_requests_file, std::ios::app);
                                        if (file119_update.is_open()) {
                                            file119_update << trim_adjustl(path2) << std::endl;
                                            file119_update.close();
                                        }

                                        conviertecabecera(cabecera, cabeceraNew, numComp + 1, rinstant);
                                        outfile_df << trim_adjustl(cabeceraNew) << std::endl;
                                        
                                        for (i1 = 0; i1 < fqLength; ++i1) {
                                            outfile_df << std::fixed << std::setprecision(15) << fqPos[i1];
                                            for (j1 = 0; j1 < numComp; ++j1) {
                                                double abs_val = std::abs(valoresDF_complex[i1 * numComp + j1]);
                                                double phase = std::atan2(std::imag(valoresDF_complex[i1 * numComp + j1]), 
                                                                          std::real(valoresDF_complex[i1 * numComp + j1]));
                                                outfile_df << " " << abs_val << " " << phase;
                                            }
                                            outfile_df << std::endl;
                                        }
                                        outfile_df.close();
                                    }

                                    valores.clear();
                                    tiempo.clear();
                                    signal.clear();
                                    samplingtime.clear();

                                    if (sgg.Observation[ii].Transfer) {
                                        path3 = trim_adjustl(sgg.Observation[ii].FileNormalize);
                                        std::ifstream test_file3(path3);
                                        existe = test_file3.good();
                                        test_file3.close();

                                        if (!existe) {
                                            buff = "Not processing: Inexistent file " + trim_adjustl(path3);
                                            print11(layoutnumber, buff);
                                        } else {
                                            neverprecounted = true;
                                            if (neverprecounted) {
                                                neverprecounted = false;
                                                std::ifstream infile3(path3);
                                                if (!infile3.is_open()) {
                                                    std::cerr << "Error opening normalize file for precounting: " << path3 << std::endl;
                                                } else {
                                                    timesteps = 0;
                                                    double d1, d2;
                                                    while (infile3 >> d1 >> d2) {
                                                        timesteps++;
                                                    }
                                                    infile3.close();
                                                    timesteps--;
                                                    
                                                    if (valores.size() > 0) {
                                                        valores.clear();
                                                        tiempo.clear();
                                                        signal.clear();
                                                    }
                                                    valores.resize(timesteps * 1, 0.0);
                                                    tiempo.resize(timesteps, 0.0);
                                                    signal.resize(timesteps, 0.0);
                                                    samplingtime.resize(timesteps, 0.0);
                                                }
                                            }

                                            std::ifstream infile3(path3);
                                            if (!infile3.is_open()) {
                                                std::cerr << "Error opening normalize file: " << path3 << std::endl;
                                            } else {
                                                for (ns = 0; ns < timesteps; ++ns) {
                                                    infile3 >> tiempo[ns] >> valores[ns];
                                                }
                                                infile3.close();
                                            }

                                            if (niapapostprocess) {
                                                std::cout << "Correcting in Transfer postprocess " << timesteps << " " << trim_adjustl(path) << std::endl;
                                                for (ns = 0; ns < timesteps; ++ns) {
                                                    tiempo[ns] = static_cast<double>(ns + 1) * sgg.dt;
                                                }
                                            }

                                            if (forceresampled) {
                                                std::string base_path = path.substr(0, path.find(".dat"));
                                                std::string resampled_name = base_path + "_resampled_time.dat";
                                                std::string path2_norm = trim_adjustl(sgg.Observation[ii].FileNormalize) + "_at_" + resampled_name;
                                                
                                                columna = 1;
                                                std::ofstream outfile_norm(path2_norm);
                                                if (!outfile_norm.is_open()) {
                                                    std::cerr << "Error opening normalization resampled file: " << path2_norm << std::endl;
                                                } else {
                                                    double t_pedido = tiempo[0];
                                                    outfile_norm << std::fixed << std::setprecision(15) << t_pedido << " " << valores[0] << std::endl;
                                                    
                                                    for (iii = 1; iii < timesteps; ++iii) {
                                                        while (t_pedido <= tiempo[iii]) {
                                                            t_pedido += sgg.Observation[ii].TimeStep;
                                                            bool found = false;
                                                            for (jjj = iii - 1; jjj < timesteps - 1; ++jjj) {
                                                                if (t_pedido >= tiempo[jjj] && t_pedido < tiempo[jjj + 1]) {
                                                                    value_interp = (valores[jjj + 1] - valores[jjj]) / 
                                                                                   (tiempo[jjj + 1] - tiempo[jjj]) * 
                                                                                   (t_pedido - tiempo[jjj]) + 
                                                                                   valores[jjj];
                                                                    outfile_norm << std::fixed << std::setprecision(15) << t_pedido << " " << value_interp << std::endl;
                                                                    found = true;
                                                                    break;
                                                                }
                                                            }
                                                            if (!found) break;
                                                        }
                                                    }
                                                    outfile_norm.close();
                                                }
                                            }

                                            std::vector<std::complex<double>> valoresDF2(fqLength * numComp, std::complex<double>(0.0, 0.0));
                                            
                                            signal.assign(timesteps, 0.0);
                                            samplingtime.assign(timesteps, 0.0);
                                            for (ns = 0; ns < timesteps; ++ns) {
                                                signal[ns] = valores[ns];
                                                samplingtime[ns] = tiempo[ns];
                                            }
                                            
                                            std::vector<std::complex<double>> fqValues_local(fqLength, std::complex<double>(0.0, 0.0));
                                            dtft(fqValues_local, fqPos, fqLength, samplingtime, signal, timesteps);
                                            
                                            for (iii = 0; iii < fqLength; ++iii) {
                                                for (j1 = 0; j1 < numComp; ++j1) {
                                                    valoresDF2[iii * numComp + j1] = valoresDF_complex[iii * numComp + j1] / fqValues_local[iii];
                                                }
                                            }

                                            path2 = path.substr(0, path.find(".dat")) + "_normalization_df.dat";
                                            std::ofstream outfile_norm_df(path2);
                                            if (!outfile_norm_df.is_open()) {
                                                std::cerr << "Error opening normalization DF file: " << path2 << std::endl;
                                            } else {
                                                outfile_norm_df << "Freq Mod_Normalization_data Phase_Normalization_data" << std::endl;
                                                for (i1 = 0; i1 < fqLength; ++i1) {
                                                    outfile_norm_df << std::fixed << std::setprecision(15) << fqPos[i1] << " " 
                                                                    << std::abs(fqValues_local[i1]) << " " 
                                                                    << std::atan2(std::imag(fqValues_local[i1]), std::real(fqValues_local[i1])) 
                                                                    << std::endl;
                                                }
                                                outfile_norm_df.close();
                                            }

                                            path2 = path.substr(0, path.find(".dat")) + "_tr.dat";
                                            std::ofstream outfile_tr(path2);
                                            if (!outfile_tr.is_open()) {
                                                std::cerr << "Error opening TR file: " << path2 << std::endl;
                                            } else {
                                                std::ofstream file119_update(output_requests_file, std::ios::app);
                                                if (file119_update.is_open()) {
                                                    file119_update << trim_adjustl(path2) << std::endl;
                                                    file119_update.close();
                                                }

                                                conviertecabecera(cabecera, cabeceraNew, numComp + 1, rinstant);
                                                pozi = path.find("_log_");
                                                if (pozi != std::string::npos) {
                                                    cabeceraNew += "_(In_dB_(20_log10)_and_radians)";
                                                }
                                                outfile_tr << trim_adjustl(cabeceraNew) << std::endl;

                                                for (i1 = 0; i1 < fqLength; ++i1) {
                                                    outfile_tr << std::fixed << std::setprecision(15) << fqPos[i1];
                                                    for (j1 = 0; j1 < numComp; ++j1) {
                                                        if (std::abs(valoresDF2[i1 * numComp + j1]) < 1e-30) {
                                                            valoresDF2[i1 * numComp + j1] = 1e-30;
                                                        }
                                                        double db = 20.0 * std::log10(std::abs(valoresDF2[i1 * numComp + j1]));
                                                        double phase = std::atan2(std::imag(valoresDF2[i1 * numComp + j1]), 
                                                                                  std::real(valoresDF2[i1 * numComp + j1]));
                                                        outfile_tr << " " << db << " " << phase;
                                                    }
                                                    outfile_tr << std::endl;
                                                }
                                                outfile_tr.close();
                                            }
                                            
                                            valoresDF2.clear();
                                            valores.clear();
                                            tiempo.clear();
                                            signal.clear();
                                            samplingtime.clear();
                                        }
                                    }
                                    
                                    fqPos.clear();
                                    fqValues.clear();
                                    valoresDF_complex.clear();
                                }
                                
                                somethingdone = true;
                            }
                        }
                    }
                }
                
                // TimeDomain / FreqDomain / Transfer check for resampling
                if ((sgg.Observation[ii].TimeDomain) || (sgg.Observation[ii].FreqDomain) || (sgg.Observation[ii].Transfer)) {
                    path = trim_adjustl(output->item[i].path);
                    
                    std::ifstream test_file(path);
                    existe = test_file.good();
                    test_file.close();

                    if (!existe) {
                        buff = "Not processing: Inexistent file " + trim_adjustl(path);
                        print11(layoutnumber, buff);
                    } else {
                        numComp = output->item[i].columnas;
                        neverprecounted = true;

                        if (neverprecounted) {
                            neverprecounted = false;
                            std::ifstream infile(path);
                            if (!infile.is_open()) {
                                std::cerr << "Error opening file for precounting: " << path << std::endl;
                            } else {
                                std::string header_line;
                                std::getline(infile, header_line);
                                
                                timesteps = 0;
                                double d1, d2;
                                while (infile >> d1 >> d2) {
                                    timesteps++;
                                }
                                infile.close();
                                timesteps--;
                                
                                if (valores.size() > 0) {
                                    valores.clear();
                                    tiempo.clear();
                                    signal.clear();
                                }
                                valores.resize(timesteps * numComp, 0.0);
                                tiempo.resize(timesteps, 0.0);
                                signal.resize(timesteps, 0.0);
                                samplingtime.resize(timesteps, 0.0);
                            }
                        }

                        std::ifstream infile(path);
                        if (!infile.is_open()) {
                            std::cerr << "Error opening file for reading: " << path << std::endl;
                        } else {
                            std::string header_line;
                            std::getline(infile, header_line);
                            
                            for (ns = 0; ns < timesteps; ++ns) {
                                infile >> tiempo[ns];
                                for (compo = 0; compo < numComp; ++compo) {
                                    infile >> valores[ns * numComp + compo];
                                }
                            }
                            infile.close();
                        }

                        if (niapapostprocess) {
                            std::cout << "Correcting in FreqDomain postprocess " << timesteps << " " << trim_adjustl(path) << std::endl;
                            for (ns = 0; ns < timesteps; ++ns) {
                                tiempo[ns] = static_cast<double>(ns + 1) * sgg.dt;
                            }
                        }

                        if (forceresampled) {
                            size_t dot_pos = path.find(".dat");
                            std::string base_path = path.substr(0, dot_pos);
                            path_resampled = base_path + "_resampled_time.dat";
                            
                            columna = 1;
                            std::ofstream outfile_resampled(path_resampled);
                            if (!outfile_resampled.is_open()) {
                                std::cerr << "Error opening resampled file: " << path_resampled << std::endl;
                            } else {
                                double t_pedido = tiempo[0];
                                outfile_resampled << std::fixed << std::setprecision(15) << t_pedido << " " << valores[0] << std::endl;
                                
                                for (iii = 1; iii < timesteps; ++iii) {
                                    while (t_pedido <= tiempo[iii]) {
                                        t_pedido += sgg.Observation[ii].TimeStep;
                                        
                                        bool found = false;
                                        for (jjj = iii - 1; jjj < timesteps - 1; ++jjj) {
                                            if (t_pedido >= tiempo[jjj] && t_pedido < tiempo[jjj + 1]) {
                                                value_interp = (valores[(jjj + 1) * numComp + (columna - 1)] - 
                                                                valores[jjj * numComp + (columna - 1)]) / 
                                                               (tiempo[jjj + 1] - tiempo[jjj]) * 
                                                               (t_pedido - tiempo[jjj]) + 
                                                               valores[jjj * numComp + (columna - 1)];
                                                outfile_resampled << std::fixed << std::setprecision(15) << t_pedido << " " << value_interp << std::endl;
                                                found = true;
                                                break;
                                            }
                                        }
                                        if (!found) break;
                                    }
                                }
                                outfile_resampled.close();
                            }
                        }
                    }
                }
            }
        }
        
        delete output;
    }
    
    void postprocessonthefly() {
        // Stub implementation as per original code structure
    }

}

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstring>
#include <cstdint>

// Assuming necessary headers and types are defined in previous chunks or included here
// e.g., SGGFDTDINFO_t, output_t, GetOutput, print11, BUFSIZE, RKIND, RKIND_tiempo, fmt, etc.
// Forward declarations or includes would go here based on the full project structure.

// Helper function to simulate trim and adjustl
inline std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, (last - first + 1));
}

// Helper function to simulate adjustl (left align, pad with spaces on right if needed, though usually just trim leading)
// In Fortran, ADJUSTL removes leading spaces.
inline std::string adjustl(const std::string& str) {
    return trim(str); // Simplified: usually just removing leading whitespace is enough for file paths/names
}

// External function declaration assumed from previous context
// std::vector<output_t>& GetOutput();
// void print11(int layoutnumber, const std::string& buff);

// Global or extern variable for format string, assumed defined elsewhere
extern const char* fmt; 

// Enum constants assumed defined elsewhere
// enum { iBloqueJx, iBloqueJy, iBloqueMx, iBloqueMy };

bool almostequal(double a, double b) {
    double ratio;
    const double tolerancia = 0.01;

    if (b == 0.0) {
        // Avoid division by zero, though logic might imply b != 0
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

void conviertecabecera(std::string& c, std::string& cNew, int columnas, double rinstant) {
    char chninstant[BUFSIZE];
    std::string c2[BUFSIZE]; // Using array of strings, max size BUFSIZE for columns
    int i, j, k, longi, ii;
    
    // Simulate write to chninstant with fmt
    // Assuming fmt is a C-style format string compatible with printf-like behavior
    // Since 'fmt' is extern, we assume it's something like "%.6f"
    snprintf(chninstant, BUFSIZE, fmt, rinstant);
    
    // Initialize c2 with spaces
    for (int idx = 0; idx < columnas; ++idx) {
        c2[idx] = " ";
    }
    
    longi = static_cast<int>(c.length()) + 1;
    
    k = 1;
    i = 1;
    
    while (i <= longi) {
        if (c[i-1] != ' ') { // Fortran is 1-based, C++ is 0-based
            int j_start = i;
            int j_end = longi;
            bool found_space = false;
            
            for (j = j_start; j <= j_end; ++j) {
                if (c[j-1] == ' ') {
                    found_space = true;
                    if (k == 1) {
                        c2[k-1] = " f_at_" + trim(std::string(chninstant));
                    } else {
                        // Extract substring from i to j-1
                        std::string sub = c.substr(i-1, j-i);
                        c2[k-1] = sub + "_Mod   " + sub + "_Phase   ";
                    }
                    k++;
                    if (k > columnas) {
                        goto end_buscam;
                    }
                    
                    // Inner loop to find next non-space or end
                    for (ii = j; ii <= longi; ++ii) {
                        if ((c[ii-1] != ' ') || (ii == longi)) {
                            i = ii;
                            goto end_interno;
                        }
                    }
                    // If loop finishes without break/exit, i should be updated? 
                    // Fortran logic: do ii=j,longi ... if exit interno. 
                    // If it doesn't exit, the loop ends. But the structure implies it always exits or i is updated.
                    // Let's assume it finds a char or end.
                    i = longi + 1; // Force exit outer loop if not found? Or just continue?
                    // Actually, if no space found until end, the inner do j loop finishes.
                    // Then we are still in the if (c(i:i) /= ' ') block.
                    // We need to ensure i advances.
                    break; 
                }
            }
            
            if (!found_space) {
                // Reached end of string without finding a space delimiter
                if (k == 1) {
                     c2[k-1] = " f_at_" + trim(std::string(chninstant));
                } else {
                     std::string sub = c.substr(i-1, longi-i+1);
                     c2[k-1] = sub + "_Mod   " + sub + "_Phase   ";
                }
                k++;
                i = longi + 1;
            }
            
        } else {
            i++;
        }
        end_interno:;
    }
    
    end_buscam:;
    
    cNew = " ";
    for (i = 1; i <= columnas; ++i) {
        cNew = trim(cNew) + "   " + trim(c2[i-1]);
    }
}

void postprocessonthefly(int layoutnumber, int num_procs, const SGGFDTDINFO_t& sgg, const std::string& nEntradaRoot, double rinstant, bool& somethingdone, bool& niapapostprocess, bool forceresampled) {
    std::vector<output_t>& output = GetOutput();
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
                    path = trim(adjustl(output[ii-1].item[i-1].path));
                    // Inquire file existence
                    std::ifstream test_file(path);
                    existe = test_file.good();
                    test_file.close();
                    
                    if (!existe) {
                        buff = "Not processing: Inexistent file " + path;
                        print11(layoutnumber, buff);
                    } else {
                        // Close the unit if it was open
                        // Assuming output[ii-1].item[i-1].unit is an ofstream or similar
                        // We need to check if it's open before closing
                        // Since we don't have the exact type of 'unit', we assume it has a close method or is an int file descriptor
                        // Based on previous chunk, it seems to be an ofstream or similar stream object wrapped or an int.
                        // If it's an int (file descriptor), we'd use close(). If ofstream, .close().
                        // Let's assume it's a stream object with a close() method or a file descriptor.
                        // Given "open (output(ii)%item(i)%unit,file=...)", it looks like a Fortran unit number or a stream.
                        // In the previous chunk translation, we likely mapped this to an ofstream or similar.
                        // If it's an ofstream:
                        if (output[ii-1].item[i-1].unit.is_open()) {
                            output[ii-1].item[i-1].unit.close();
                        }
                    }
                }
            }
        }
    }
    
    // Call the main postprocess subroutine
    // Assuming postprocess is defined in the same namespace/module
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
                    path = trim(adjustl(output[ii-1].item[i-1].path));
                    // Inquire file existence
                    std::ifstream test_file(path);
                    existe = test_file.good();
                    test_file.close();
                    
                    if (!existe) {
                        buff = "ERROR: Inexistent file. Creating new " + path;
                        print11(layoutnumber, buff);
                    } else {
                        // Open file in append mode
                        // Assuming output[ii-1].item[i-1].unit is an ofstream
                        output[ii-1].item[i-1].unit.open(path, std::ios::out | std::ios::app);
                    }
                }
            }
        }
    }
}