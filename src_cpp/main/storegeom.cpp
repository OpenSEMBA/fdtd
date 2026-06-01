#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

// Assuming these types and constants are defined in FDETYPES_m
// Since the prompt asks to translate the provided code, we assume the 
// dependent types (media_matrices_t, SGGFDTDINFO_t) and constants 
// (INTEGERSIZEOFMEDIAMATRICES, iEx, iEy, etc.) are available in the 
// global namespace or included headers. For the sake of a standalone 
// translation of the logic, we will define placeholders or assume 
// they are included. However, to strictly follow "Preserve ALL names", 
// we will use the names as given.

// Placeholder definitions for external dependencies to make the code 
// syntactically valid C++ if compiled in isolation. In a real project, 
// these would come from FDETYPES_m and other modules.

#ifndef INTEGERSIZEOFMEDIAMATRICES
#define INTEGERSIZEOFMEDIAMATRICES int
#endif

// Enumerations for field components
enum FieldComponent {
    iEx = 1,
    iEy = 2,
    iEz = 3,
    iHx = 4,
    iHy = 5,
    iHz = 6
};

// Placeholder for Is_t structure
struct Is_t {
    bool PML;
    bool PEC;
    bool ThinWire;
    bool SlantedWire;
    bool EDispersive;
    bool MDispersive;
    bool ThinSlot;
    bool SGBC;
    bool Lossy;
    bool multiport;
    bool anismultiport;
    bool multiportpadding;
    bool dielectric;
    bool Anisotropic;
    bool Needed;
    bool already_YEEadvanced_byconformal;
    bool split_and_useless;
    bool Volume;
    bool Surface;
    bool Line;
};

// Placeholder for Media_t structure
struct Media_t {
    int Priority;
    double Epr;
    double Sigma;
    double Mur;
    Is_t Is;
    double SigmaM;
};

// Placeholder for Sweep_t structure
struct Sweep_t {
    int XI;
    int XE;
    int YI;
    int YE;
    int ZI;
    int ZE;
};

// Placeholder for Alloc_t structure
struct Alloc_t {
    int XI;
    int XE;
};

// Placeholder for SGGFDTDINFO_t
struct SGGFDTDINFO_t {
    int NumMedia;
    std::vector<Media_t> Med;
    std::vector<Sweep_t> sweep;
    std::vector<Sweep_t> SINPMLsweep;
    std::vector<Alloc_t> Alloc;
};

// Placeholder for media_matrices_t
struct media_matrices_t {
    // Assuming 3D arrays for media matrices. 
    // Dimensions are not specified in the snippet, so we use a generic 
    // vector of vectors of vectors. The indices i, j, k are used.
    std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>> sggMiEx;
    std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>> sggMiEy;
    std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>> sggMiEz;
    std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>> sggMiHx;
    std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>> sggMiHy;
    std::vector<std::vector<std::vector<INTEGERSIZEOFMEDIAMATRICES>>> sggMiHz;
};

namespace storeData_m {

    // Helper function to translate media indexes into characters
    char chartranslate(INTEGERSIZEOFMEDIAMATRICES entero) {
        if (entero == 1) {
            return '_';
        } else if (entero == 0) {
            return '0';
        } else if (entero == -1) {
            return '#';
        } else {
            // Fortran char(48+Abs(entero)) converts ASCII value. 
            // 48 is '0'. So if entero is 2, it returns '2'.
            // Note: If entero > 9, this might produce non-digit characters.
            // We assume entero is small enough or behaves as expected in Fortran context.
            return static_cast<char>(48 + std::abs(entero));
        }
    }

    void store_geomData(const media_matrices_t& media, const SGGFDTDINFO_t& sgg, const std::string& fileFDE) {
        // Open files
        // Fortran: open(20, FILE=trim(adjustl(fileFDE))//'_MapEx.txt')
        // C++: We use fstream. Note: Fortran unit numbers 20-25 are used.
        // We will map them to file streams.
        
        std::string baseName = fileFDE;
        // trim and adjustl equivalent in C++:
        // Assuming fileFDE is already clean or we just append.
        
        std::ofstream fileEx(baseName + "_MapEx.txt");
        std::ofstream fileEy(baseName + "_MapEy.txt");
        std::ofstream fileEz(baseName + "_MapEz.txt");
        std::ofstream fileHx(baseName + "_MapHx.txt");
        std::ofstream fileHy(baseName + "_MapHy.txt");
        std::ofstream fileHz(baseName + "_MapHz.txt");

        if (!fileEx.is_open() || !fileEy.is_open() || !fileEz.is_open() ||
            !fileHx.is_open() || !fileHy.is_open() || !fileHz.is_open()) {
            std::cerr << "Error opening files for writing." << std::endl;
            return;
        }

        // Map unit numbers to streams for easy access in the loop
        // 20: Ex, 21: Ey, 22: Ez, 23: Hx, 24: Hy, 25: Hz
        std::ofstream* files[6] = {&fileEx, &fileEy, &fileEz, &fileHx, &fileHy, &fileHz};

        for (int campo = 1; campo <= 6; ++campo) {
            int i = 19 + campo;
            int q = 19 + campo;
            
            // Write header
            // Fortran: write(q,*) '____ 1-Sustrato, -n PML_______'
            // q is the unit number. In C++, we use the corresponding file stream.
            // Since q is 20-25, we use files[campo-1]
            if (q >= 20 && q <= 25) {
                (*files[q - 20]) << "____ 1-Sustrato, -n PML_______" << std::endl;
            }

            // Loop over media
            for (int j = 0; j <= sgg.NumMedia; ++j) {
                INTEGERSIZEOFMEDIAMATRICES INTJ = static_cast<INTEGERSIZEOFMEDIAMATRICES>(j);
                
                if (q >= 20 && q <= 25) {
                    (*files[q - 20]) << "_____________________________" << std::endl;
                    (*files[q - 20]) << "MEDIO :  " << INTJ << std::endl;
                    (*files[q - 20]) << "Priority " << sgg.Med[j].Priority << std::endl;
                    (*files[q - 20]) << "Epr " << sgg.Med[j].Epr << std::endl;
                    (*files[q - 20]) << "Sigma " << sgg.Med[j].Sigma << std::endl;
                    (*files[q - 20]) << "Mur " << sgg.Med[j].Mur << std::endl;
                    (*files[q - 20]) << "Is PML " << sgg.Med[j].Is.PML << std::endl;
                    (*files[q - 20]) << "Is PEC " << sgg.Med[j].Is.PEC << std::endl;
                    (*files[q - 20]) << "SigmaM " << sgg.Med[j].SigmaM << std::endl;
                    (*files[q - 20]) << "Is ThinWIRE " << sgg.Med[j].Is.ThinWire << std::endl;
                    (*files[q - 20]) << "Is SlantedWIRE " << sgg.Med[j].Is.SlantedWire << std::endl;
                    (*files[q - 20]) << "Is EDispersive " << sgg.Med[j].Is.EDispersive << std::endl;
                    (*files[q - 20]) << "Is MDispersive " << sgg.Med[j].Is.MDispersive << std::endl;
                    (*files[q - 20]) << "Is ThinSlot " << sgg.Med[j].Is.ThinSlot << std::endl;
                    (*files[q - 20]) << "Is SGBC " << sgg.Med[j].Is.SGBC << std::endl;
                    (*files[q - 20]) << "Is Lossy " << sgg.Med[j].Is.Lossy << std::endl;
                    (*files[q - 20]) << "Is Multiport " << sgg.Med[j].Is.multiport << std::endl;
                    (*files[q - 20]) << "Is AnisMultiport " << sgg.Med[j].Is.anismultiport << std::endl;
                    (*files[q - 20]) << "Is MultiportPadding " << sgg.Med[j].Is.multiportpadding << std::endl;
                    (*files[q - 20]) << "Is Dielectric " << sgg.Med[j].Is.dielectric << std::endl;
                    (*files[q - 20]) << "Is ThinSlot " << sgg.Med[j].Is.ThinSlot << std::endl;
                    (*files[q - 20]) << "Is Anisotropic " << sgg.Med[j].Is.Anisotropic << std::endl;
                    (*files[q - 20]) << "Is Needed " << sgg.Med[j].Is.Needed << std::endl;
                    (*files[q - 20]) << "Is already_YEEadvanced_byconformal " << sgg.Med[j].Is.already_YEEadvanced_byconformal << std::endl;
                    (*files[q - 20]) << "Is split_and_useless " << sgg.Med[j].Is.split_and_useless << std::endl;
                    (*files[q - 20]) << "Is Volume " << sgg.Med[j].Is.Volume << std::endl;
                    (*files[q - 20]) << "Is Surface " << sgg.Med[j].Is.Surface << std::endl;
                    (*files[q - 20]) << "Is Line " << sgg.Med[j].Is.Line << std::endl;
                }
            }

            // Write PML info
            if (i >= 20 && i <= 25) {
                (*files[i - 20]) << campo << " con PML IINIC, IFIN " << sgg.sweep[campo].XI << " " << sgg.sweep[campo].XE << std::endl;
                (*files[i - 20]) << campo << " con PML JINIC, JFIN " << sgg.sweep[campo].YI << " " << sgg.sweep[campo].YE << std::endl;
                (*files[i - 20]) << campo << " con PML KINIC, KFIN " << sgg.sweep[campo].ZI << " " << sgg.sweep[campo].ZE << std::endl;
                (*files[i - 20]) << campo << " sin PML IINIC, IFIN " << sgg.SINPMLsweep[campo].XI << " " << sgg.SINPMLsweep[campo].XE << std::endl;
                (*files[i - 20]) << campo << " sin PML JINIC, JFIN " << sgg.SINPMLsweep[campo].YI << " " << sgg.SINPMLsweep[campo].YE << std::endl;
                (*files[i - 20]) << campo << " sin PML KINIC, KFIN " << sgg.SINPMLsweep[campo].ZI << " " << sgg.SINPMLsweep[campo].ZE << std::endl;
            }

            // Write grid maps
            for (int k = sgg.sweep[campo].ZI; k <= sgg.sweep[campo].ZE; ++k) {
                i = 19 + campo;
                if (i >= 20 && i <= 25) {
                    (*files[i - 20]) << "_______________________________________________________________________" << std::endl;
                    (*files[i - 20]) << "!!!!!!** k=" << k << std::endl;
                    
                    // Header line with indices
                    // Fortran: write(19+campo, '(A,400a)') 'I=  |', ('0123456789', i=sgg%Alloc(campo)%XI, sgg%Alloc(campo)%XE+10, 10)
                    // This prints digits 0-9 repeatedly.
                    std::string header = "I=  |";
                    int startIdx = sgg.Alloc[campo].XI;
                    int endIdx = sgg.Alloc[campo].XE + 10;
                    for (int idx = startIdx; idx < endIdx; ++idx) {
                        header += std::to_string(idx % 10);
                    }
                    (*files[i - 20]) << header << std::endl;
                    
                    (*files[i - 20]) << "J______________________________________________________________________" << std::endl;
                }

                // Loop J from YE down to YI
                for (int j = sgg.sweep[campo].YE; j >= sgg.sweep[campo].YI; --j) {
                    i = 19 + campo;
                    if (i >= 20 && i <= 25) {
                        INTEGERSIZEOFMEDIAMATRICES val = 0;
                        switch (campo) {
                            case iEx:
                                val = media.sggMiEx[i][j][k];
                                break;
                            case iEy:
                                val = media.sggMiEy[i][j][k];
                                break;
                            case iEz:
                                val = media.sggMiEz[i][j][k];
                                break;
                            case iHx:
                                val = media.sggMiHx[i][j][k];
                                break;
                            case iHy:
                                val = media.sggMiHy[i][j][k];
                                break;
                            case iHz:
                                val = media.sggMiHz[i][j][k];
                                break;
                            default:
                                val = 0;
                                break;
                        }
                        
                        std::string charVal = chartranslate(val);
                        // Fortran: write(19+campo, '(I3,A,4000a)') j, ' |', (chartranslate(...))
                        // This prints J (3 digits), " |", and then a string of characters.
                        // The string length is determined by the loop range XI to XE.
                        std::string rowStr = "";
                        int xi = sgg.sweep[campo].XI;
                        int xe = sgg.sweep[campo].XE;
                        for (int idx = xi; idx <= xe; ++idx) {
                            INTEGERSIZEOFMEDIAMATRICES v = 0;
                            switch (campo) {
                                case iEx: v = media.sggMiEx[idx][j][k]; break;
                                case iEy: v = media.sggMiEy[idx][j][k]; break;
                                case iEz: v = media.sggMiEz[idx][j][k]; break;
                                case iHx: v = media.sggMiHx[idx][j][k]; break;
                                case iHy: v = media.sggMiHy[idx][j][k]; break;
                                case iHz: v = media.sggMiHz[idx][j][k]; break;
                                default: v = 0; break;
                            }
                            rowStr += chartranslate(v);
                        }
                        
                        (*files[i - 20]) << std::setw(3) << j << " |" << rowStr << std::endl;
                    }
                }
            }
        }

        // Close files
        for (int i = 20; i <= 25; ++i) {
            if (i >= 20 && i <= 25) {
                files[i - 20]->close();
            }
        }
    }

} // namespace storeData_m