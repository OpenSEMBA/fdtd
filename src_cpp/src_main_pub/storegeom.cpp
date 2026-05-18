#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdint>

// Forward declarations for types defined in FDETYPES_m
// Assuming FDETYPES_m defines these types and constants.
// In a real translation, these would be included from a header.

// Placeholder for INTEGERSIZEOFMEDIAMATRICES kind
using INTEGERSIZEOFMEDIAMEDIAMATRICES = int32_t;

// Placeholder for iEx, iEy, etc. constants
enum FieldComponent {
    iEx = 1,
    iEy = 2,
    iEz = 3,
    iHx = 4,
    iHy = 5,
    iHz = 6
};

// Placeholder for media_matrices_t
struct media_matrices_t {
    // Assuming these are 3D arrays. Fortran is column-major, C++ is row-major.
    // We will use std::vector with a flattened index or 3D vector.
    // For simplicity and preserving logic, we'll use a 3D vector structure.
    // Note: The original code accesses media%sggMiEx(i, j, k).
    // We assume these are allocated appropriately.
    std::vector<std::vector<std::vector<std::string>>> sggMiEx;
    std::vector<std::vector<std::vector<std::string>>> sggMiEy;
    std::vector<std::vector<std::vector<std::string>>> sggMiEz;
    std::vector<std::vector<std::vector<std::string>>> sggMiHx;
    std::vector<std::vector<std::vector<std::string>>> sggMiHy;
    std::vector<std::vector<std::vector<std::string>>> sggMiHz;
};

// Placeholder for SGGFDTDINFO_t
struct SGGFDTDINFO_t {
    int NumMedia;
    struct {
        int Priority;
        double Epr;
        double Sigma;
        double Mur;
        struct {
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
        } Is;
        double SigmaM;
    }* Med; // Array of media info, 1-based in Fortran, so we might pad or adjust indexing

    struct {
        int XI, XE;
        int YI, YE;
        int ZI, ZE;
    }* sweep; // Array of sweep info, 1-based

    struct {
        int XI, XE;
        int YI, YE;
        int ZI, ZE;
    }* SINPMLsweep; // Array of sinpml sweep info, 1-based

    struct {
        int XI, XE;
    }* Alloc; // Array of alloc info, 1-based
};

namespace storeData_m {

    void store_geomData(const media_matrices_t& media, const SGGFDTDINFO_t& sgg, const std::string& fileFDE) {
        std::ofstream files[6];
        std::string suffixes[] = {"_MapEx.txt", "_MapEy.txt", "_MapEz.txt", "_MapHx.txt", "_MapHy.txt", "_MapHz.txt"};
        
        for (int i = 0; i < 6; ++i) {
            files[i].open(fileFDE + suffixes[i]);
            if (!files[i].is_open()) {
                std::cerr << "Error opening file: " << fileFDE << suffixes[i] << std::endl;
                return;
            }
        }

        for (int campo = 1; campo <= 6; ++campo) {
            int i = 19 + campo;
            int q = 19 + campo;
            
            // Map 1-based campo to 0-based index for files array
            int file_idx = campo - 1;
            
            files[file_idx] << "____ 1-Sustrato, -n PML_______" << std::endl;
            
            for (int j = 0; j <= sgg.NumMedia; ++j) {
                // Fortran uses 1-based indexing for Med. 
                // Assuming sgg.Med is allocated with size NumMedia+1 or accessed 1-based.
                // In C++, if Med is a raw array or vector, we need to handle indexing.
                // Let's assume Med is a pointer to an array where index 1 is valid.
                // We'll access sgg.Med[j] directly if it's 1-based in memory layout or adjusted.
                // For this translation, we assume the struct members are accessible as in Fortran.
                
                int INTJ = j;
                files[file_idx] << "_____________________________" << std::endl;
                files[file_idx] << "MEDIO :  " << chartranslate(INTJ) << std::endl;
                files[file_idx] << "Priority " << sgg.Med[j].Priority << std::endl;
                files[file_idx] << "Epr " << sgg.Med[j].Epr << std::endl;
                files[file_idx] << "Sigma " << sgg.Med[j].Sigma << std::endl;
                files[file_idx] << "Mur " << sgg.Med[j].Mur << std::endl;
                files[file_idx] << "Is PML " << sgg.Med[j].Is.PML << std::endl;
                files[file_idx] << "Is PEC " << sgg.Med[j].Is.PEC << std::endl;
                files[file_idx] << "SigmaM " << sgg.Med[j].SigmaM << std::endl;
                files[file_idx] << "Is ThinWIRE " << sgg.Med[j].Is.ThinWire << std::endl;
                files[file_idx] << "Is SlantedWIRE " << sgg.Med[j].Is.SlantedWire << std::endl;
                files[file_idx] << "Is EDispersive " << sgg.Med[j].Is.EDispersive << std::endl;
                files[file_idx] << "Is MDispersive " << sgg.Med[j].Is.MDispersive << std::endl;
                files[file_idx] << "Is ThinSlot " << sgg.Med[j].Is.ThinSlot << std::endl;
                files[file_idx] << "Is SGBC " << sgg.Med[j].Is.SGBC << std::endl;
                files[file_idx] << "Is Lossy " << sgg.Med[j].Is.Lossy << std::endl;
                files[file_idx] << "Is Multiport " << sgg.Med[j].Is.multiport << std::endl;
                files[file_idx] << "Is AnisMultiport " << sgg.Med[j].Is.anismultiport << std::endl;
                files[file_idx] << "Is MultiportPadding " << sgg.Med[j].Is.multiportpadding << std::endl;
                files[file_idx] << "Is Dielectric " << sgg.Med[j].Is.dielectric << std::endl;
                files[file_idx] << "Is ThinSlot " << sgg.Med[j].Is.ThinSlot << std::endl;
                files[file_idx] << "Is Anisotropic " << sgg.Med[j].Is.Anisotropic << std::endl;
                files[file_idx] << "Is Needed " << sgg.Med[j].Is.Needed << std::endl;
                files[file_idx] << "Is already_YEEadvanced_byconformal " << sgg.Med[j].Is.already_YEEadvanced_byconformal << std::endl;
                files[file_idx] << "Is split_and_useless " << sgg.Med[j].Is.split_and_useless << std::endl;
                files[file_idx] << "Is Volume " << sgg.Med[j].Is.Volume << std::endl;
                files[file_idx] << "Is Surface " << sgg.Med[j].Is.Surface << std::endl;
                files[file_idx] << "Is Line " << sgg.Med[j].Is.Line << std::endl;
            }
            
            files[file_idx] << campo << " con PML IINIC, IFIN " << sgg.sweep[campo].XI << " " << sgg.sweep[campo].XE << std::endl;
            files[file_idx] << campo << " con PML JINIC, JFIN " << sgg.sweep[campo].YI << " " << sgg.sweep[campo].YE << std::endl;
            files[file_idx] << campo << " con PML KINIC, KFIN " << sgg.sweep[campo].ZI << " " << sgg.sweep[campo].ZE << std::endl;
            files[file_idx] << campo << " sin PML IINIC, IFIN " << sgg.SINPMLsweep[campo].XI << " " << sgg.SINPMLsweep[campo].XE << std::endl;
            files[file_idx] << campo << " sin PML JINIC, JFIN " << sgg.SINPMLsweep[campo].YI << " " << sgg.SINPMLsweep[campo].YE << std::endl;
            files[file_idx] << campo << " sin PML KINIC, KFIN " << sgg.SINPMLsweep[campo].ZI << " " << sgg.SINPMLsweep[campo].ZE << std::endl;
            
            for (int k = sgg.sweep[campo].ZI; k <= sgg.sweep[campo].ZE; ++k) {
                files[file_idx] << "_______________________________________________________________________" << std::endl;
                files[file_idx] << "!!!!!!** k=" << k << std::endl;
                
                // Header line: I= | 0123456789...
                std::string header = "I=  |";
                for (int ii = sgg.Alloc[campo].XI; ii <= sgg.Alloc[campo].XE + 10; ii += 10) {
                    header += "0123456789";
                }
                files[file_idx] << header << std::endl;
                files[file_idx] << "J______________________________________________________________________" << std::endl;
                
                for (int j = sgg.sweep[campo].YE; j >= sgg.sweep[campo].YI; --j) {
                    std::string val;
                    switch (campo) {
                        case iEx:
                            val = chartranslate(media.sggMiEx[campo][j][k]); // Assuming 1-based indexing for media arrays too
                            break;
                        case iEy:
                            val = chartranslate(media.sggMiEy[campo][j][k]);
                            break;
                        case iEz:
                            val = chartranslate(media.sggMiEz[campo][j][k]);
                            break;
                        case iHx:
                            val = chartranslate(media.sggMiHx[campo][j][k]);
                            break;
                        case iHy:
                            val = chartranslate(media.sggMiHy[campo][j][k]);
                            break;
                        case iHz:
                            val = chartranslate(media.sggMiHz[campo][j][k]);
                            break;
                        default:
                            val = "?";
                            break;
                    }
                    files[file_idx] << j << " |" << val << std::endl;
                }
            }
        }
        
        for (int i = 0; i < 6; ++i) {
            if (files[i].is_open()) {
                files[i].close();
            }
        }
    }

    std::string chartranslate(int entero) {
        if (entero == 1) {
            return "_";
        } else if (entero == 0) {
            return "0";
        } else if (entero == -1) {
            return "#";
        } else {
            // Fortran: char(48+Abs(entero))
            // ASCII 48 is '0'. So this converts integer to its digit character representation.
            // Note: This assumes entero is a single digit. If not, it might produce non-digit chars.
            // Fortran's CHAR function returns a character with the given ASCII code.
            char c = static_cast<char>(48 + std::abs(entero));
            return std::string(1, c);
        }
    }

} // namespace storeData_m