#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <memory>
#include <algorithm>

// Assuming necessary type definitions from included modules are available
// or defined here for compilation context.
// RKIND is typically double
using RKIND = double;
const RKIND RKIND_VAL = 1.0;

// Placeholder types to match Fortran structures
struct media_matrices_t {
    // ... members ...
};

struct MedioExtra_t {
    // ... members ...
};

struct limit_t {
    // ... members ...
};

struct SGGFDTDINFO_t {
    bool thereAreMagneticMedia;
    bool thereArePMLMagneticMedia;
    double dt;
    std::vector<double> DX;
    std::vector<double> DY;
    std::vector<double> DZ;
    
    struct {
        int Conta;
        int MaxConta;
        std::vector<void*> elem; // Placeholder for complex type
    } EShared;
};

struct Parseador_t {
    struct {
        std::vector<void*> volumes; // Placeholder
    } conformalRegs;
};

struct taglist_t {
    // ... members ...
};

struct tagtype_t {
    // ... members ...
};

struct FreqDepenMaterial_t {
    // ... members ...
};

struct ConformalMedia_t {
    double time_step_scale_factor;
};

struct side_tris_map_t {
    // ... members ...
};

struct XYZlimit_t {
    // ... members ...
};

struct xyzlimit_scaled_t {
    // ... members ...
};

// Global variables from the module
namespace Preprocess_m {

    RKIND cluz = 0.0;
    RKIND zvac = 0.0;
    RKIND eps0 = 0.0;
    RKIND mu0 = 0.0;

    // Forward declarations for functions called within this module
    void cuentatags(Parseador_t& this_, tagtype_t& tagtype, int layoutnumber, const std::string& fichin);
    void prepro_skindepth(Parseador_t& this_, const std::string& fichin);
    void print11(int layoutnumber, const std::string& msg);
    std::vector<ConformalMedia_t> buildConformalMedia(void* conformalRegs); // Placeholder signature
    std::vector<side_tris_map_t> buildSideMaps(void* conformalRegs); // Placeholder signature

    void read_geomData(media_matrices_t& media, taglist_t& tag_numbers, const std::string& fichin, 
                       int layoutnumber, int num_procs, 
                       std::vector<limit_t>& SINPML_fullsize, std::vector<limit_t>& fullsize, 
                       Parseador_t& this_, 
                       std::vector<bool>& groundwires, double attfactor, 
                       bool& mibc, bool& SGBC, bool& SGBCDispersive, 
                       MedioExtra_t& MEDIOEXTRA, double maxSourceValue, 
                       double skindepthpre, bool& createmapvtk, 
                       bool input_conformal_flag, bool& CLIPREGION, 
                       double boundwireradius, double maxwireradius, 
                       bool updateshared, bool run_with_dmma, 
                       double eps00, double mu00, bool simu_devia, 
                       bool hay_slanted_wires, bool verbose, 
                       bool ignoresamplingerrors, tagtype_t& tagtype, 
                       const std::string& wiresflavor) {
        
        // Local variables
        int tama, tama2, tama3, tama4, tama5, tama6, i, j, k, tipotemp, tamaSonda;
        int tamaoldSONDA, tamaBloquePrb, tamaScrPrb, pozi, tama2bis, numeroasignaciones, ci;
        std::string probenumber;
        RKIND ex, ey, ez, px, py, pz, amplitud, minSpaceStep;
        std::string tag;
        XYZlimit_t punto, BoundingBox, conf_bounding_box;
        xyzlimit_scaled_t punto_s;
        int orientacion, orientacionL, orientacionR, direccion, contamedia, oldcontamedia;
        int maxcontamedia, mincontamedia, inicontamedia, i1, j1, field, k1, pecmedio, ii, medio1, medio2;
        int sondas, CONTACURR, CONTAVOLT, I_, J_;
        bool isathinwire, VALIDO, existia, medioespecial, nodo_cazado;
        bool errnofile, errnofile1, errnofile2, errnofile3, errnofile4;
        RKIND tiempo1, tiempo2, field1, field2, rdummy;
        int nsurfs, numus, OrigIndex, numminus;
        RKIND delta, del, sig_max;
        std::vector<int> contapuntos;
        int conta1, conta2, MEDIO, imenos1, jmenos1, kmenos1, o, p, puntoxi, puntoyi, puntozi;
        int bboxwirXI, dummy_bboxwirXI, bboxwirYI, dummy_bboxwirYI, bboxwirzI, dummy_bboxwirzI;
        int bboxwirXE, dummy_bboxwirXE, bboxwirYE, dummy_bboxwirYE, bboxwirZE, dummy_bboxwirZE, IERR;
        long long memo;
        std::string MultiportFile, MultiportFile2, buff;
        std::vector<RKIND> dummy_px, dummy_py, dummy_pz, dummy_ex, dummy_ey, dummy_ez, dummy_INCERT;
        std::string whoami, whoamishort, ext, extpoint, chari, charj, chark, chari2, charj2, chark2;
        bool paraerrhilo, islossy, DENTRO;
        RKIND width;
        std::vector<RKIND> dir(3);
        RKIND epr1, mur1;
        bool oriX, oriY, oriZ, oriX2, oriY2, oriZ2, oriX3, oriY3, oriZ3, iguales;
        bool oriX4, oriY4, oriZ4;
        std::vector<std::vector<RKIND>> EprSlot(3, std::vector<RKIND>(3));
        std::vector<std::vector<RKIND>> MurSlot(3, std::vector<RKIND>(3));
        int indicemedio, i11, j11;
        FreqDepenMaterial_t* fdgeom = nullptr;
        int numertag;
        int Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE;
        int Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE;
        int Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE;
        int Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE;
        int Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE;
        int Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE;

        std::vector<ConformalMedia_t> conformal_media;
        std::vector<RKIND> edge_ratios;
        std::vector<side_tris_map_t> side_to_triangles_maps;

        // Initialize globals
        eps0 = eps00;
        mu0 = mu00;
        cluz = 1.0_RKIND_VAL / std::sqrt(eps0 * mu0);
        zvac = std::sqrt(mu0 / eps0);

        call cuentatags(this_, tagtype, layoutnumber, fichin);

        delta = -1.0; // Initialize to avoid warnings

        sgg.thereAreMagneticMedia = true;
        sgg.thereArePMLMagneticMedia = true;

        if (skindepthpre) {
            if (layoutnumber == 0) {
                print11(layoutnumber, "Preprocessing SGBC materials to include skin-depth effects....");
                prepro_skindepth(this_, fichin);
                print11(layoutnumber, "Finished preprocessing for skin-depth.");
            }
            // MPI_Barrier placeholder
            return;
        }

        whoami = "(" + std::to_string(layoutnumber + 1) + "/" + std::to_string(num_procs) + ") ";
        whoamishort = std::to_string(layoutnumber + 1);

        // Create space for etangential shared info
        sgg.EShared.Conta = 0;
        sgg.EShared.MaxConta = 10;
        sgg.EShared.elem.resize(sgg.EShared.MaxConta);

        {
            int m;
            RKIND min_scale_factor = 1.0;
            RKIND dt;
            if (!this_.conformalRegs.volumes.empty()) { // Placeholder check for associated
                conformal_media = buildConformalMedia(this_.conformalRegs.volumes); // Placeholder cast
                side_to_triangles_maps = buildSideMaps(this_.conformalRegs.volumes); // Placeholder cast
                
                for (m = 0; m < conformal_media.size(); ++m) {
                    if (conformal_media[m].time_step_scale_factor < min_scale_factor) {
                        min_scale_factor = conformal_media[m].time_step_scale_factor;
                    }
                }
                
                // Calculate dt
                RKIND inv_dx_sq = 0.0, inv_dy_sq = 0.0, inv_dz_sq = 0.0;
                if (!sgg.DX.empty()) inv_dx_sq = (1.0_RKIND_VAL / sgg.DX[0]) * (1.0_RKIND_VAL / sgg.DX[0]);
                if (!sgg.DY.empty()) inv_dy_sq = (1.0_RKIND_VAL / sgg.DY[0]) * (1.0_RKIND_VAL / sgg.DY[0]);
                if (!sgg.DZ.empty()) inv_dz_sq = (1.0_RKIND_VAL / sgg.DZ[0]) * (1.0_RKIND_VAL / sgg.DZ[0]);
                
                dt = (1.0_RKIND_VAL / (cluz * std::sqrt(inv_dx_sq + inv_dy_sq + inv_dz_sq)));
                
                if (sgg.dt > dt * min_scale_factor) {
                    std::cout << "-- Conformal geometry requires a time step change" << std::endl;
                    std::cout << "Previous time step: " << sgg.dt << std::endl;
                    sgg.dt = dt * min_scale_factor;
                    std::cout << "New time step: " << sgg.dt << std::endl;
                }
            } else {
                conformal_media.resize(0);
            }
        }

        // Count media, calculate sizes, reserve space
        // ... (Rest of the logic would follow here)
    }

    void read_limits_nogeom() {
        // Stub
    }

    void AssigLossyOrPECtoNodes() {
        // Stub
    }

    void searchtag() {
        // Stub
    }

    void checkDielectricTagForDuplicate() {
        // Stub
    }

    void checkAnimatedTagForDuplicate() {
        // Stub
    }

    void checkLossyTagForDuplicate() {
        // Stub
    }
}

// el medio 0 se reserva para PEC
      // regiones PEC
      // el medio 1 se reserva para sustrato  y saltamos
      contamedia = 1;
      if ((this->pmcregs.nvols) + (this->pmcregs.nsurfs) + (this->pmcregs.nLINS) != 0) {
         // los PMC empiezan en 2
         contamedia = 2;
         // fin regions PMC
      }
      // materialList
      // NonMetalREgions   and frequencydependent media
      // Anisotropic
      //
      contamedia = contamedia + (this->DielRegs.nvols) + (this->DielRegs.nsurfs) + (this->DielRegs.nLINS) + this->FRQDEPMATS.nvols +
      this->FRQDEPMATS.nsurfs + this->FRQDEPMATS.nLINS + (this->ANIMATS.nvols + this->ANIMATS.nsurfs + this->ANIMATS.nLINS);
      // Multiports
      // worst case 6 orientations per surface plus the the lossy padding
      contamedia = contamedia + this->LossyThinSurfs.length * 7;
      // wires
      // nueva formulacion que almacena also the lenghts
      contamedia = contamedia + this->twires.n_tw;
      contamedia = contamedia + this->swires.n_sw;
      // echo por demas, habria que precontar pero es complicado porque depende del procesamiento
      // thin Slots

      if (run_with_dmma) {
         for (j = 0; j < this->tSlots.n_tg; ++j) {
            contamedia = contamedia + this->tSlots.Tg[j].N_tgc;
         }
      }

      // end thin Slots
      // PARA LA CAPA EXTRA 2013
      if (medioextra.exists) {
         CONTAMEDIA = CONTAMEDIA + 1;
         MEDIOEXTRA.index = CONTAMEDIA;
      }
      // para modulos que necesien senialar con already_YEEadvanced_byconformal y split_and_useless (eg. conformal)
      // se crea siempre por defecto
      contamedia = contamedia + 2; // para acomodar los no_use no_use_notouch
      // !!!!!!!!!!!!
      contamedia = contamedia + 1; // para acomodar los nodal sources como caso especial de linea vacia

      // contamedia = contamedia + this->conformalRegs.nEdges + this->conformalRegs.nFaces
      
      edge_ratios = getDifferentEdgeRatios(conformal_media);
      face_ratios = getDifferentFaceRatios(conformal_media);
      contamedia = contamedia + static_cast<int>(edge_ratios.size()) + static_cast<int>(face_ratios.size());
      if (std::find(edge_ratios.begin(), edge_ratios.end(), 0.0) != edge_ratios.end()) contamedia = contamedia - 1;
      if (std::find(face_ratios.begin(), face_ratios.end(), 0.0) != face_ratios.end()) contamedia = contamedia - 1;

#ifdef CompileWithMTLN
      contamedia = contamedia + this->mtln.n_unsh;
#endif

      sgg.NumMedia = contamedia;
      sgg.AllocMed = contamedia;
      // reserva espacio
     sgg.Med.resize(sgg.NumMedia + 1);
      // comienzo barrido resto :  medios y observaciones
      BoundingBox.XI = sgg.Alloc[iHx].XI;
      BoundingBox.XE = sgg.Alloc[iHx].XE;
      BoundingBox.YI = sgg.Alloc[iHy].YI;
      BoundingBox.YE = sgg.Alloc[iHy].YE;
      BoundingBox.ZI = sgg.Alloc[iHz].ZI;
      BoundingBox.ZE = sgg.Alloc[iHz].ZE;
      //
      Alloc_iEx_XI = sgg.Alloc[iEx].XI;
      Alloc_iEx_XE = sgg.Alloc[iEx].XE;
      Alloc_iEx_YI = sgg.Alloc[iEx].YI;
      Alloc_iEx_YE = sgg.Alloc[iEx].YE;
      Alloc_iEx_ZI = sgg.Alloc[iEx].ZI;
      Alloc_iEx_ZE = sgg.Alloc[iEx].ZE;
      Alloc_iEy_XI = sgg.Alloc[iEy].XI;
      Alloc_iEy_XE = sgg.Alloc[iEy].XE;
      Alloc_iEy_YI = sgg.Alloc[iEy].YI;
      Alloc_iEy_YE = sgg.Alloc[iEy].YE;
      Alloc_iEy_ZI = sgg.Alloc[iEy].ZI;
      Alloc_iEy_ZE = sgg.Alloc[iEy].ZE;
      Alloc_iEz_XI = sgg.Alloc[iEz].XI;
      Alloc_iEz_XE = sgg.Alloc[iEz].XE;
      Alloc_iEz_YI = sgg.Alloc[iEz].YI;
      Alloc_iEz_YE = sgg.Alloc[iEz].YE;
      Alloc_iEz_ZI = sgg.Alloc[iEz].ZI;
      Alloc_iEz_ZE = sgg.Alloc[iEz].ZE;
      Alloc_iHx_XI = sgg.Alloc[iHx].XI;
      Alloc_iHx_XE = sgg.Alloc[iHx].XE;
      Alloc_iHx_YI = sgg.Alloc[iHx].YI;
      Alloc_iHx_YE = sgg.Alloc[iHx].YE;
      Alloc_iHx_ZI = sgg.Alloc[iHx].ZI;
      Alloc_iHx_ZE = sgg.Alloc[iHx].ZE;
      Alloc_iHy_XI = sgg.Alloc[iHy].XI;
      Alloc_iHy_XE = sgg.Alloc[iHy].XE;
      Alloc_iHy_YI = sgg.Alloc[iHy].YI;
      Alloc_iHy_YE = sgg.Alloc[iHy].YE;
      Alloc_iHy_ZI = sgg.Alloc[iHy].ZI;
      Alloc_iHy_ZE = sgg.Alloc[iHy].ZE;
      Alloc_iHz_XI = sgg.Alloc[iHz].XI;
      Alloc_iHz_XE = sgg.Alloc[iHz].XE;
      Alloc_iHz_YI = sgg.Alloc[iHz].YI;
      Alloc_iHz_YE = sgg.Alloc[iHz].YE;
      Alloc_iHz_ZI = sgg.Alloc[iHz].ZI;
      Alloc_iHz_ZE = sgg.Alloc[iHz].ZE;
      //
      //
      field = 1;
      //
      numertag = 0;
     media.sggMtag.resize(Alloc_iHx_XE - Alloc_iHx_XI + 1, std::vector<std::vector<int>>(Alloc_iHy_YE - Alloc_iHy_YI + 1, std::vector<int>(Alloc_iHz_ZE - Alloc_iHz_ZI + 1)));
     media.sggMiNo.resize(Alloc_iHx_XE - Alloc_iHx_XI + 1, std::vector<std::vector<int>>(Alloc_iHy_YE - Alloc_iHy_YI + 1, std::vector<int>(Alloc_iHz_ZE - Alloc_iHz_ZI + 1)));

     tag_numbers.edge.x.resize(Alloc_iEx_XE - Alloc_iEx_XI + 1, std::vector<std::vector<int>>(Alloc_iEx_YE - Alloc_iEx_YI + 1, std::vector<int>(Alloc_iEx_ZE - Alloc_iEx_ZI + 1)));
     tag_numbers.edge.y.resize(Alloc_iEy_XE - Alloc_iEy_XI + 1, std::vector<std::vector<int>>(Alloc_iEy_YE - Alloc_iEy_YI + 1, std::vector<int>(Alloc_iEy_ZE - Alloc_iEy_ZI + 1)));
     tag_numbers.edge.z.resize(Alloc_iEz_XE - Alloc_iEz_XI + 1, std::vector<std::vector<int>>(Alloc_iEz_YE - Alloc_iEz_YI + 1, std::vector<int>(Alloc_iEz_ZE - Alloc_iEz_ZI + 1)));
     tag_numbers.face.x.resize(Alloc_iHx_XE - Alloc_iHx_XI + 1, std::vector<std::vector<int>>(Alloc_iHx_YE - Alloc_iHx_YI + 1, std::vector<int>(Alloc_iHx_ZE - Alloc_iHx_ZI + 1)));
     tag_numbers.face.y.resize(Alloc_iHy_XE - Alloc_iHy_XI + 1, std::vector<std::vector<int>>(Alloc_iHy_YE - Alloc_iHy_YI + 1, std::vector<int>(Alloc_iHy_ZE - Alloc_iHy_ZI + 1)));
     tag_numbers.face.z.resize(Alloc_iHz_XE - Alloc_iHz_XI + 1, std::vector<std::vector<int>>(Alloc_iHz_YE - Alloc_iHz_YI + 1, std::vector<int>(Alloc_iHz_ZE - Alloc_iHz_ZI + 1)));

      !!!nodos materiales: se precisan para el conformal !sgg310715
     media.sggMiEx.resize(Alloc_iEx_XE - Alloc_iEx_XI + 1, std::vector<std::vector<int>>(Alloc_iEx_YE - Alloc_iEx_YI + 1, std::vector<int>(Alloc_iEx_ZE - Alloc_iEx_ZI + 1)));
     media.sggMiEy.resize(Alloc_iEy_XE - Alloc_iEy_XI + 1, std::vector<std::vector<int>>(Alloc_iEy_YE - Alloc_iEy_YI + 1, std::vector<int>(Alloc_iEy_ZE - Alloc_iEy_ZI + 1)));
     media.sggMiEz.resize(Alloc_iEz_XE - Alloc_iEz_XI + 1, std::vector<std::vector<int>>(Alloc_iEz_YE - Alloc_iEz_YI + 1, std::vector<int>(Alloc_iEz_ZE - Alloc_iEz_ZI + 1)));
     media.sggMiHx.resize(Alloc_iHx_XE - Alloc_iHx_XI + 1, std::vector<std::vector<int>>(Alloc_iHx_YE - Alloc_iHx_YI + 1, std::vector<int>(Alloc_iHx_ZE - Alloc_iHx_ZI + 1)));
     media.sggMiHy.resize(Alloc_iHy_XE - Alloc_iHy_XI + 1, std::vector<std::vector<int>>(Alloc_iHy_YE - Alloc_iHy_YI + 1, std::vector<int>(Alloc_iHy_ZE - Alloc_iHy_ZI + 1)));
     media.sggMiHz.resize(Alloc_iHz_XE - Alloc_iHz_XI + 1, std::vector<std::vector<int>>(Alloc_iHz_YE - Alloc_iHz_YI + 1, std::vector<int>(Alloc_iHz_ZE - Alloc_iHz_ZI + 1)));


      // el tag esta voided porque luego el numero va con el del tag
      for (int i = 0; i < media.sggMtag.size(); ++i) {
         for (int j = 0; j < media.sggMtag[i].size(); ++j) {
            for (int k = 0; k < media.sggMtag[i][j].size(); ++k) {
               media.sggMtag[i][j][k] = 0; // LO VOIDEO A 0 EN VEZ DE A -1 PORQUE EL TAG 0 NO VA A EXISTIR NUNCA 141020
            }
         }
      }
      for (int i = 0; i < tag_numbers.edge.x.size(); ++i) {
         for (int j = 0; j < tag_numbers.edge.x[i].size(); ++j) {
            for (int k = 0; k < tag_numbers.edge.x[i][j].size(); ++k) {
               tag_numbers.edge.x[i][j][k] = 0;
            }
         }
      }
      for (int i = 0; i < tag_numbers.edge.y.size(); ++i) {
         for (int j = 0; j < tag_numbers.edge.y[i].size(); ++j) {
            for (int k = 0; k < tag_numbers.edge.y[i][j].size(); ++k) {
               tag_numbers.edge.y[i][j][k] = 0;
            }
         }
      }
      for (int i = 0; i < tag_numbers.edge.z.size(); ++i) {
         for (int j = 0; j < tag_numbers.edge.z[i].size(); ++j) {
            for (int k = 0; k < tag_numbers.edge.z[i][j].size(); ++k) {
               tag_numbers.edge.z[i][j][k] = 0;
            }
         }
      }
      for (int i = 0; i < tag_numbers.face.x.size(); ++i) {
         for (int j = 0; j < tag_numbers.face.x[i].size(); ++j) {
            for (int k = 0; k < tag_numbers.face.x[i][j].size(); ++k) {
               tag_numbers.face.x[i][j][k] = 0;
            }
         }
      }
      for (int i = 0; i < tag_numbers.face.y.size(); ++i) {
         for (int j = 0; j < tag_numbers.face.y[i].size(); ++j) {
            for (int k = 0; k < tag_numbers.face.y[i][j].size(); ++k) {
               tag_numbers.face.y[i][j][k] = 0;
            }
         }
      }
      for (int i = 0; i < tag_numbers.face.z.size(); ++i) {
         for (int j = 0; j < tag_numbers.face.z[i].size(); ++j) {
            for (int k = 0; k < tag_numbers.face.z[i][j].size(); ++k) {
               tag_numbers.face.z[i][j][k] = 0;
            }
         }
      }
      // todo sustrato por defecto
      for (int i = 0; i < media.sggMiNo.size(); ++i) {
         for (int j = 0; j < media.sggMiNo[i].size(); ++j) {
            for (int k = 0; k < media.sggMiNo[i][j].size(); ++k) {
               media.sggMiNo[i][j][k] = 1;
            }
         }
      }
      for (int i = 0; i < media.sggMiEx.size(); ++i) {
         for (int j = 0; j < media.sggMiEx[i].size(); ++j) {
            for (int k = 0; k < media.sggMiEx[i][j].size(); ++k) {
               media.sggMiEx[i][j][k] = 1;
            }
         }
      }
      for (int i = 0; i < media.sggMiEy.size(); ++i) {
         for (int j = 0; j < media.sggMiEy[i].size(); ++j) {
            for (int k = 0; k < media.sggMiEy[i][j].size(); ++k) {
               media.sggMiEy[i][j][k] = 1;
            }
         }
      }
      for (int i = 0; i < media.sggMiEz.size(); ++i) {
         for (int j = 0; j < media.sggMiEz[i].size(); ++j) {
            for (int k = 0; k < media.sggMiEz[i][j].size(); ++k) {
               media.sggMiEz[i][j][k] = 1;
            }
         }
      }
      for (int i = 0; i < media.sggMiHx.size(); ++i) {
         for (int j = 0; j < media.sggMiHx[i].size(); ++j) {
            for (int k = 0; k < media.sggMiHx[i][j].size(); ++k) {
               media.sggMiHx[i][j][k] = 1;
            }
         }
      }
      for (int i = 0; i < media.sggMiHy.size(); ++i) {
         for (int j = 0; j < media.sggMiHy[i].size(); ++j) {
            for (int k = 0; k < media.sggMiHy[i][j].size(); ++k) {
               media.sggMiHy[i][j][k] = 1;
            }
         }
      }
      for (int i = 0; i < media.sggMiHz.size(); ++i) {
         for (int j = 0; j < media.sggMiHz[i].size(); ++j) {
            for (int k = 0; k < media.sggMiHz[i][j].size(); ++k) {
               media.sggMiHz[i][j][k] = 1;
            }
         }
      }
      

      // planeWaves
      //
      tama = (this->plnSrc.nc);
!!!      write(buff,*) 'More than 1 Huygens box unsupported'
!!!      if (tama > 1) call STOPONERROR(layoutnumber,num_procs,buff)
      // LO PONGO A MANO ojo
      amplitud = 1.0;
      sgg.NumPlaneWaves = tama;
     sgg.PlaneWave.resize(sgg.NumPlaneWaves + 1);
      for (i = 0; i < sgg.NumPlaneWaves; ++i) {
         punto.XI = std::min(this->plnSrc.collection[i].coor1[0], this->plnSrc.collection[i].coor2[0]);
         punto.XE = std::max(this->plnSrc.collection[i].coor1[0], this->plnSrc.collection[i].coor2[0]);
         punto.YI = std::min(this->plnSrc.collection[i].coor1[1], this->plnSrc.collection[i].coor2[1]);
         punto.YE = std::max(this->plnSrc.collection[i].coor1[1], this->plnSrc.collection[i].coor2[1]);
         punto.ZI = std::min(this->plnSrc.collection[i].coor1[2], this->plnSrc.collection[i].coor2[2]);
         punto.ZE = std::max(this->plnSrc.collection[i].coor1[2], this->plnSrc.collection[i].coor2[2]);
         // just for the sake of peace of my mind
         // readjust Huygens surface CLEARLY in/out in case of coincidente
         if ((punto.XI == SINPML_fullsize[iHx].XI)) {
            punto.XI = SINPML_fullsize[iHx].XI - 5;
         }
         if ((punto.XE == SINPML_fullsize[iHx].XE)) {
            punto.XE = SINPML_fullsize[iHx].XE + 5;
         }
         if ((punto.YI == SINPML_fullsize[iHy].YI)) {
            punto.YI = SINPML_fullsize[iHy].YI - 5;
         }
         if ((punto.YE == SINPML_fullsize[iHy].YE)) {
            punto.YE = SINPML_fullsize[iHy].YE + 5;
         }
         if ((punto.ZI == SINPML_fullsize[iHz].ZI)) {
            punto.ZI = SINPML_fullsize[iHz].ZI - 5;
         }

}
         if (punto.ZE == SINPML_fullsize[iHz].ZE) {
            punto.ZE = SINPML_fullsize[iHz].ZE + 5;
         }
         //
         sgg.PlaneWave[i].isRC = this->plnSrc.collection[i].isRC;
         sgg.PlaneWave[i].numModes = this->plnSrc.collection[i].numModes;
         sgg.PlaneWave[i].incertmax = this->plnSrc.collection[i].incertmax;
         if (sgg.PlaneWave[i].isRC) {
            sgg.PlaneWave[i].px.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].py.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].pz.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].ex.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].ey.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].ez.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].INCERT.resize(sgg.PlaneWave[i].numModes);
            dummy_px.resize(sgg.PlaneWave[i].numModes);
            dummy_py.resize(sgg.PlaneWave[i].numModes);
            dummy_pz.resize(sgg.PlaneWave[i].numModes);
            dummy_ex.resize(sgg.PlaneWave[i].numModes);
            dummy_ey.resize(sgg.PlaneWave[i].numModes);
            dummy_ez.resize(sgg.PlaneWave[i].numModes);
            dummy_INCERT.resize(sgg.PlaneWave[i].numModes);
            std::fill(sgg.PlaneWave[i].px.begin(), sgg.PlaneWave[i].px.end(), 0.0);
            std::fill(sgg.PlaneWave[i].py.begin(), sgg.PlaneWave[i].py.end(), 0.0);
            std::fill(sgg.PlaneWave[i].pz.begin(), sgg.PlaneWave[i].pz.end(), 0.0);
            std::fill(sgg.PlaneWave[i].ex.begin(), sgg.PlaneWave[i].ex.end(), 0.0);
            std::fill(sgg.PlaneWave[i].ey.begin(), sgg.PlaneWave[i].ey.end(), 0.0);
            std::fill(sgg.PlaneWave[i].ez.begin(), sgg.PlaneWave[i].ez.end(), 0.0);
            std::fill(sgg.PlaneWave[i].INCERT.begin(), sgg.PlaneWave[i].INCERT.end(), 0.0);
            if (layoutnumber == 0) populatePlaneWaveRC(sgg.PlaneWave[i], simu_devia); // only the master populates
#ifdef CompileWithMPI
            MPI_Barrier(MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].px.data(), dummy_px.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].py.data(), dummy_py.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].pz.data(), dummy_pz.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].ex.data(), dummy_ex.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].ey.data(), dummy_ey.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].ez.data(), dummy_ez.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_AllReduce(sgg.PlaneWave[i].INCERT.data(), dummy_INCERT.data(), sgg.PlaneWave[i].numModes, REALSIZE, MPI_SUM, MPI_COMM_WORLD, &ierr);
            MPI_Barrier(MPI_COMM_WORLD, &ierr);
            sgg.PlaneWave[i].px = dummy_px;
            sgg.PlaneWave[i].py = dummy_py;
            sgg.PlaneWave[i].pz = dummy_pz;
            sgg.PlaneWave[i].ex = dummy_ex;
            sgg.PlaneWave[i].ey = dummy_ey;
            sgg.PlaneWave[i].ez = dummy_ez;
            sgg.PlaneWave[i].INCERT = dummy_INCERT;
#endif
            dummy_px.clear();
            dummy_py.clear();
            dummy_pz.clear();
            dummy_ex.clear();
            dummy_ey.clear();
            dummy_ez.clear();
            dummy_INCERT.clear();
         } else {
            sgg.PlaneWave[i].numModes = 1;
            sgg.PlaneWave[i].px.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].py.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].pz.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].ex.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].ey.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].ez.resize(sgg.PlaneWave[i].numModes);
            sgg.PlaneWave[i].INCERT.resize(sgg.PlaneWave[i].numModes);
            //
            ez = amplitud * std::cos(this->plnSrc.collection[i].alpha);
            ey = amplitud * std::sin(this->plnSrc.collection[i].alpha) * std::sin(this->plnSrc.collection[i].beta);
            ex = amplitud * std::sin(this->plnSrc.collection[i].alpha) * std::cos(this->plnSrc.collection[i].beta);
            pz = std::cos(this->plnSrc.collection[i].theta);
            py = std::sin(this->plnSrc.collection[i].theta) * std::sin(this->plnSrc.collection[i].phi);
            px = std::sin(this->plnSrc.collection[i].theta) * std::cos(this->plnSrc.collection[i].phi);
            //ojo con estos redondeos.
            //!!!if (std::abs(ex/amplitud) < 1e-4) ex = 0.0;
            //!!!if (std::abs(ey/amplitud) < 1e-4) ey = 0.0;
            //!!!if (std::abs(ez/amplitud) < 1e-4) ez = 0.0;
            //!!!if (std::abs(px) < 1e-4) px = 0.0;
            //!!!if (std::abs(py) < 1e-4) py = 0.0;
            //!!!if (std::abs(pz) < 1e-4) pz = 0.0;
            if (std::abs(px * ex + py * ey + pz * ez) >= 1e-4) {
               std::ostringstream buff_stream;
               buff_stream << "NO TEM PLANEWAVE " << ex << " " << ey << " " << ez << " " << px << " " << py << " " << pz << " " << (px * ex + py * ey + pz * ez) << " " << this->plnSrc.collection[i].alpha << " "
                  << this->plnSrc.collection[i].beta << " " << this->plnSrc.collection[i].theta << " " << this->plnSrc.collection[i].phi;
               std::string buff = buff_stream.str();
               STOPONERROR(layoutnumber, num_procs, buff);
            }
            //
            sgg.PlaneWave[i].px[0] = px;
            sgg.PlaneWave[i].py[0] = py;
            sgg.PlaneWave[i].pz[0] = pz;
            sgg.PlaneWave[i].ex[0] = ex;
            sgg.PlaneWave[i].ey[0] = ey;
            sgg.PlaneWave[i].ez[0] = ez;
            sgg.PlaneWave[i].INCERT[0] = 0.0;
         }
         sgg.PlaneWave[i].fichero.name = trim(adjustl(this->plnSrc.collection[i].nombre_fichero));
         sgg.PlaneWave[i].esqx1 = std::min(punto.XI, punto.XE);
         sgg.PlaneWave[i].esqy1 = std::min(punto.YI, punto.YE);
         sgg.PlaneWave[i].esqz1 = std::min(punto.ZI, punto.ZE);
         sgg.PlaneWave[i].esqx2 = std::max(punto.XI, punto.XE);
         sgg.PlaneWave[i].esqy2 = std::max(punto.YI, punto.YE);
         sgg.PlaneWave[i].esqz2 = std::max(punto.ZI, punto.ZE);
      }
      //Media parsing
      //Default
      //background
      sgg.Med.Priority = prior_BV;
      sgg.Med.Epr = 1.0;
      sgg.Med.Sigma = 0.0;
      sgg.Med.Sigmareasignado = false; // solo afecta a un chequeo de errores en lumped 120123
      sgg.Med.Mur = 1.0;
      sgg.Med.SigmaM = 0.0;
      sgg.Med.Is.Interfase = false;
      sgg.Med.Is.PMLbody = false;
      sgg.Med.Is.Needed = true;
      sgg.Med.Is.Anisotropic = false;
      sgg.Med.Is.Dielectric = false;
      sgg.Med.Is.EDispersive = false;
      sgg.Med.Is.EDispersiveAnis = false;
      sgg.Med.Is.MDispersive = false;
      sgg.Med.Is.MDispersiveAnis = false;
      sgg.Med.Is.Lumped = false;
      sgg.Med.Is.SGBC = false;
      sgg.Med.Is.SGBCDispersive = false;
      sgg.Med.Is.Lossy = false;
      sgg.Med.Is.multiport = false;
      sgg.Med.Is.multiportpadding = false;
      sgg.Med.Is.AnisMultiport = false;
      sgg.Med.Is.ThinWire = false;
      sgg.Med.Is.Multiwire = false;
      sgg.Med.Is.SlantedWire = false;
      sgg.Med.Is.ThinSlot = false;
      sgg.Med.Is.PEC = false;
      sgg.Med.Is.ConformalPEC = false;
      sgg.Med.Is.PMC = false;
      sgg.Med.Is.PML = false;
      sgg.Med.Is.Volume = false;
      sgg.Med.Is.Surface = false;
      sgg.Med.Is.Line = false;
      sgg.Med.Is.already_YEEadvanced_byconformal = false;
      sgg.Med.Is.split_and_useless = false;
      //ojo tocar tambien en el readjust de healing si se crean nuevos flags
      //
      //medio PEC y PML es intrascendente si es surface o volume
      //son los de prioridad mas alta y siempre contienen a sus campos tangenciales electricos
      //Background    only differences from default are needed
      sgg.Med[1].Priority = prior_BV;
      sgg.Med[1].Epr = this->mats.mats[1].eps / Eps0;
      sgg.Med[1].Sigma = this->mats.mats[1].Sigma;
      sgg.Med[1].Mur = this->mats.mats[1].mu / Mu0;
      sgg.Med[1].SigmaM = this->mats.mats[1].SigmaM;
      sgg.Med[1].Is.Dielectric = false; // considero el vacio como NO dielectrico '251114
      sgg.Med[1].Is.Volume = false; // considero el vacio como no volumic false '251114
      //
      sgg.Med[0].Is.PEC = true;
      sgg.Med[0].Is.Needed = true;
      sgg.Med[0].Priority = prior_PEC;

sgg.Med[0].Epr = this->mats.mats[1].eps / Eps0;
sgg.Med[0].Sigma = 1.0e29;
sgg.Med[0].Mur = this->mats.mats[1].mu / Mu0;
sgg.Med[0].SigmaM = 0.0;

// CAPA EXTRA
// Background    only differences from default are needed

if (medioextra.exists) {
    // estimate in terms of percentage of the maximum PML conductivity the conductivity of the extra medium
    // This info is available from read_limits_nogeom
    // the calculus is taken from borderscpml.F90
    sig_max = 0.0;
    for (int o = 1; o <= 3; ++o) {
        for (int p = 1; p <= 2; ++p) {
            double del;
            if ((o == 1) && (p == 1)) del = sgg.dx(SINPML_fullsize[iHx].XI);
            if ((o == 1) && (p == 2)) del = sgg.dx(SINPML_fullsize[iHx].XE - 1);
            if ((o == 2) && (p == 1)) del = sgg.dy(SINPML_fullsize[iHy].YI);
            if ((o == 2) && (p == 2)) del = sgg.dy(SINPML_fullsize[iHy].YE - 1);
            if ((o == 3) && (p == 1)) del = sgg.dz(SINPML_fullsize[iHz].ZI);
            if ((o == 3) && (p == 2)) del = sgg.dz(SINPML_fullsize[iHz].ZE - 1);

            if (sgg.PML.NumLayers[o][p] != 0) {
                if ((sgg.PML.NumLayers[o][p] == 10) || (sgg.PML.NumLayers[o][p] == 5)) {
                    sig_max = std::max(sig_max, 0.8 * (sgg.PML.orden[o][p] + 1) / (zvac * del));
                } else {
                    if (sgg.PML.CoeffReflPML[o][p] == 1.0) {
                        // realmente en el borderscpml
                        // sig_max(sig_max,-((log( 0.99999d0                 )*(sgg%PML%orden(o,p)+1))/ &
                        //    (2.0_RKIND *sqrt(Mu0/eps0)*sgg%PML%NumLayers(o,p)*del)))
                        // trampa para que entonces tome la conductividad autentica que se especifique y poder anular las PML y solo dejar capa fisica !!?!?
                        sig_max = 1.0;
                    } else {
                        sig_max = std::max(sig_max, -((std::log(sgg.PML.CoeffReflPML[o][p]) * (sgg.PML.orden[o][p] + 1)) /
                            (2.0 * std::sqrt(Mu0 / eps0) * sgg.PML.NumLayers[o][p] * del)));
                    }
                }
            }
        }
    }
    MEDIOEXTRA.sigma = MEDIOEXTRA.sigma * sig_max; // la especificacion se da en terminos de tanto por uno en la linea de comandos

    sgg.Med[MEDIOEXTRA.index].Epr = this->mats.mats[1].eps / Eps0; // luego se machaca este valor
    sgg.Med[MEDIOEXTRA.index].Sigma = MEDIOEXTRA.sigma; // luego se machaca este valor
    sgg.Med[MEDIOEXTRA.index].Mur = this->mats.mats[1].mu / Mu0; // luego se machaca este valor
    sgg.Med[MEDIOEXTRA.index].SigmaM = 0.0; // solo lo creo para las tangenciales electricas
    sgg.Med[MEDIOEXTRA.index].Priority = prior_PEC;
    sgg.Med[MEDIOEXTRA.index].Is.Dielectric = true;
    sgg.Med[MEDIOEXTRA.index].Is.Volume = true;
    sgg.Med[MEDIOEXTRA.index].Is.PML = true;
}

// barre los medios
// Primero todos los pec
// PECRegions
// volumenes
// el medio 0 se reserva para PEC
// regiones PEC

if ((this->pecregs.nvols) + (this->pecregs.nsurfs) + (this->pecregs.nLINS) != 0) {
    pecmedio = 0;
    tama = (this->pecregs.nvols);
    // BODYes
    for (int i = 1; i <= tama; ++i) {
        punto.XI = this->pecregs.vols[i].XI;
        punto.XE = this->pecregs.vols[i].XE;
        punto.YI = this->pecregs.vols[i].YI;
        punto.YE = this->pecregs.vols[i].YE;
        punto.ZI = this->pecregs.vols[i].ZI;
        punto.ZE = this->pecregs.vols[i].ZE;
        numertag = searchtag(tagtype, this->pecregs.vols[i].tag);
        CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                       media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                       Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                       Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                       Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                       Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                       Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, pecmedio);
    }
    // SURFs
    tama = (this->pecregs.nsurfs);
    for (int i = 1; i <= tama; ++i) {
        punto.XI = this->pecregs.surfs[i].XI;
        punto.XE = this->pecregs.surfs[i].XE;
        punto.YI = this->pecregs.surfs[i].YI;
        punto.YE = this->pecregs.surfs[i].YE;
        punto.ZI = this->pecregs.surfs[i].ZI;
        punto.ZE = this->pecregs.surfs[i].ZE;
        orientacion = this->pecregs.surfs[i].or;
        numertag = searchtag(tagtype, this->pecregs.surfs[i].tag);
        CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                        media.sggMiHx, media.sggMiHy, media.sggMiHz,
                        Alloc_iEx_XI, Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                        Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE,
                        Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, Alloc_iEz_ZE,
                        Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE,
                        Alloc_iHy_XI, Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE,
                        Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE,
                        sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, pecmedio);
    }
    // LINs
    tama = (this->pecregs.nLINS);
    for (int i = 1; i <= tama; ++i) {
        punto.XI = this->pecregs.lins[i].XI;
        punto.XE = this->pecregs.lins[i].XE;
        punto.YI = this->pecregs.lins[i].YI;
        punto.YE = this->pecregs.lins[i].YE;
        punto.ZI = this->pecregs.lins[i].ZI;
        punto.ZE = this->pecregs.lins[i].ZE;
        orientacion = this->pecregs.lins[i].or;
        isathinwire = false;
        numertag = searchtag(tagtype, this->pecregs.lins[i].tag);
        CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                     media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                     Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                     Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                     Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                     Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                     Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                     pecmedio, isathinwire, verbose, numeroasignaciones);
    }
    // regiones PEC
}

// el medio 1 se reserva para sustrato  y saltamos
contamedia = 1;

// para el conformal !debe ser tipicamente contamedia =1+1=2 pq el 0 es pec y el 1 es vacio. Ojo cambiado de sitio el PMC porque podia hacer que fuesen 3 y 4. 130220!!! y puede haber error pq por ahi se comprueba el 2 y el 3
contamedia = contamedia + 1;
sgg.Med[contamedia].Is.already_YEEadvanced_byconformal = true;
// debe ser contamedia =2+1=3
contamedia = contamedia + 1;
sgg.Med[contamedia].Is.split_and_useless = true;

// cambiado aqui 130220

// materialList
// regiones PMC
if ((this->pmcregs.nvols) + (this->pmcregs.nsurfs) + (this->pmcregs.nLINS) != 0) {
    // los PMC de existir tienen todos indice 2
}

contamedia = contamedia + 1; // !!!!contamedia = 2 !!!ufff. cambiado a 130220 por posible bug con conformal si algun dia habia regiones PMC
            sgg.Med[contamedia].Epr = sgg.Med[1].Epr;
            sgg.Med[contamedia].Mur = sgg.Med[1].Mur;
            sgg.Med[contamedia].Sigma = 0.0;
            sgg.Med[contamedia].SigmaM = 1.0e29;
            sgg.Med[contamedia].Priority = prior_PMC;
            sgg.Med[contamedia].Is.PMC = true;
            //BODYes
            tama = this->pmcregs.nvols;
            for (i = 1; i <= tama; ++i) {
                punto.XI = this->pmcregs.vols[i].XI;
                punto.XE = this->pmcregs.vols[i].XE;
                punto.YI = this->pmcregs.vols[i].YI;
                punto.YE = this->pmcregs.vols[i].YE;
                punto.ZI = this->pmcregs.vols[i].ZI;
                punto.ZE = this->pmcregs.vols[i].ZE;
                //
                //
                numertag = searchtag(tagtype, this->pmcregs.vols[i].tag);
                CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
            }
            //SURFs
            tama = this->pmcregs.nsurfs;
            for (i = 1; i <= tama; ++i) {
                punto.XI = this->pmcregs.surfs[i].XI;
                punto.XE = this->pmcregs.surfs[i].XE;
                punto.YI = this->pmcregs.surfs[i].YI;
                punto.YE = this->pmcregs.surfs[i].YE;
                punto.ZI = this->pmcregs.surfs[i].ZI;
                punto.ZE = this->pmcregs.surfs[i].ZE;
                orientacion = this->pmcregs.surfs[i].or;
                numertag = searchtag(tagtype, this->pmcregs.surfs[i].tag);
                CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                    contamedia);
            }
            //LINs
            tama = this->pmcregs.nLINS;
            for (i = 1; i <= tama; ++i) {
                punto.XI = this->pmcregs.lins[i].XI;
                punto.XE = this->pmcregs.lins[i].XE;
                punto.YI = this->pmcregs.lins[i].YI;
                punto.YE = this->pmcregs.lins[i].YE;
                punto.ZI = this->pmcregs.lins[i].ZI;
                punto.ZE = this->pmcregs.lins[i].ZE;
                orientacion = this->pmcregs.lins[i].or;
                isathinwire = false;
                numertag = searchtag(tagtype, this->pmcregs.lins[i].tag);
                CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                    contamedia, isathinwire, verbose, numeroasignaciones);
                //
            }
            //fin regions PMC
        }
        !!!!fin cambiado 130220

        //NonMetalREgions
        //BODYes
        tama = this->DielRegs.nvols;
        for (i = 1; i <= tama; ++i) {
            contamedia = contamedia + 1;
            sgg.Med[contamedia].Is.Dielectric = true;
            sgg.Med[contamedia].Priority = prior_IB;
            sgg.Med[contamedia].Epr = this->DielRegs.vols[i].eps / Eps0;
            sgg.Med[contamedia].Sigma = this->DielRegs.vols[i].Sigma;
            sgg.Med[contamedia].Mur = this->DielRegs.vols[i].mu / Mu0;
            sgg.Med[contamedia].SigmaM = this->DielRegs.vols[i].SigmaM;
            !!!!pmlbody
            if (this->DielRegs.vols[i].PMLbody) {
                sgg.Med[contamedia].Priority = prior_pmlbody; // machaca con una prioridad superior a la de thin wires y backgroud !prueba HOLD 251019 coax
                sgg.Med[contamedia].Is.PMLbody = true;
                sgg.Med[contamedia].PMLbody.resize(1);
                sgg.Med[contamedia].PMLbody[0].orient = this->DielRegs.vols[i].orient;
            }
            !!!!!
            tama2 = this->DielRegs.vols[i].n_c2P;
            for (j = 1; j <= tama2; ++j) {
                if ((j == 1) && (this->DielRegs.vols[i].PMLbody)) sgg.Med[contamedia].PMLbody[0].orient = this->DielRegs.vols[i].c2P[j].OR; // ES IGUAL PARA TODOS
                punto.XI = this->DielRegs.vols[i].c2P[j].XI;
                punto.XE = this->DielRegs.vols[i].c2P[j].XE;
                punto.YI = this->DielRegs.vols[i].c2P[j].YI;
                punto.YE = this->DielRegs.vols[i].c2P[j].YE;
                punto.ZI = this->DielRegs.vols[i].c2P[j].ZI;
                punto.ZE = this->DielRegs.vols[i].c2P[j].ZE;
                numertag = searchtag(tagtype, this->DielRegs.vols[i].c2P[j].tag);
                CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
            }
            tama3 = this->DielRegs.vols[i].n_c1P;
            for (j = 1; j <= tama3; ++j) {
                if ((j == 1) && (this->DielRegs.vols[i].PMLbody)) sgg.Med[contamedia].PMLbody[0].orient = this->DielRegs.vols[i].c1P[j].OR; // ES IGUAL PARA TODOS
                punto.XI = this->DielRegs.vols[i].c1P[j].XI;
                punto.XE = this->DielRegs.vols[i].c1P[j].XI;
                punto.YI = this->DielRegs.vols[i].c1P[j].YI;
                punto.YE = this->DielRegs.vols[i].c1P[j].YI;
                punto.ZI = this->DielRegs.vols[i].c1P[j].ZI;
                punto.ZE = this->DielRegs.vols[i].c1P[j].ZI;
                //
                numertag = searchtag(tagtype, this->DielRegs.vols[i].c1P[j].tag);
                CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,

Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
         }
      }
      //SURFs
      tama = (this->DielRegs.nsurfs);
      for (i = 0; i < tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med[contamedia].Is.Dielectric = true;
         sgg.Med[contamedia].Priority = prior_IS;
         sgg.Med[contamedia].Epr = this->DielRegs.surfs[i].eps / Eps0;
         sgg.Med[contamedia].Sigma = this->DielRegs.surfs[i].Sigma;
         sgg.Med[contamedia].Mur = this->DielRegs.surfs[i].mu / Mu0;
         sgg.Med[contamedia].SigmaM = this->DielRegs.surfs[i].SigmaM;
         tama2 = (this->DielRegs.surfs[i].n_c2P);
         for (j = 0; j < tama2; ++j) {
            punto.XI = this->DielRegs.surfs[i].c2P[j].XI;
            punto.XE = this->DielRegs.surfs[i].c2P[j].XE;
            punto.YI = this->DielRegs.surfs[i].c2P[j].YI;
            punto.YE = this->DielRegs.surfs[i].c2P[j].YE;
            punto.ZI = this->DielRegs.surfs[i].c2P[j].ZI;
            punto.ZE = this->DielRegs.surfs[i].c2P[j].ZE;
            orientacion = this->DielRegs.surfs[i].c2P[j].or;
            numertag = searchtag(tagtype, this->DielRegs.surfs[i].c2P[j].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
         tama3 = (this->DielRegs.surfs[i].n_c1P);

         for (j = 0; j < tama3; ++j) {
            punto.XI = this->DielRegs.surfs[i].c1P[j].XI;
            punto.XE = this->DielRegs.surfs[i].c1P[j].XI;
            punto.YI = this->DielRegs.surfs[i].c1P[j].YI;
            punto.YE = this->DielRegs.surfs[i].c1P[j].YI;
            punto.ZI = this->DielRegs.surfs[i].c1P[j].ZI;
            punto.ZE = this->DielRegs.surfs[i].c1P[j].ZI;
            orientacion = this->DielRegs.surfs[i].c1P[j].or;
            numertag = searchtag(tagtype, this->DielRegs.surfs[i].c1P[j].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
      }
      //LINs
      tama = (this->DielRegs.nLINS);
      for (i = 0; i < tama; ++i) {
         numeroasignaciones = 0; //solo lo usa lumped para echarselo al primer y el resto ponerlo a PEC
         contamedia = contamedia + 1;
         sgg.Med[contamedia].Is.Dielectric = true;
         sgg.Med[contamedia].Priority = prior_IL;
         sgg.Med[contamedia].Epr = this->DielRegs.lins[i].eps / Eps0;
         sgg.Med[contamedia].Sigma = this->DielRegs.lins[i].Sigma;
         sgg.Med[contamedia].Mur = this->DielRegs.lins[i].mu / Mu0;
         sgg.Med[contamedia].SigmaM = this->DielRegs.lins[i].SigmaM;
         //lumped
         if (this->DielRegs.lins[i].resistor) {
            sgg.Med[contamedia].Is.Lumped = true;
            sgg.Med[contamedia].Is.lossy = true; //importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med[contamedia].lumped.resize(1);
            sgg.Med[contamedia].lumped[0].resistor = true;
            sgg.Med[contamedia].lumped[0].inductor = false;
            sgg.Med[contamedia].lumped[0].capacitor = false;
            sgg.Med[contamedia].lumped[0].diodo = false;
            sgg.Med[contamedia].lumped[0].R = this->DielRegs.lins[i].R;
            sgg.Med[contamedia].lumped[0].L = 0.0;
            sgg.Med[contamedia].lumped[0].C = 0.0;
            sgg.Med[contamedia].lumped[0].R_devia = this->DielRegs.lins[i].R_devia;
            sgg.Med[contamedia].lumped[0].L_devia = 0.0;
            sgg.Med[contamedia].lumped[0].C_devia = 0.0;
            sgg.Med[contamedia].lumped[0].Rtime_on = this->DielRegs.lins[i].Rtime_on;
            sgg.Med[contamedia].lumped[0].Rtime_off = this->DielRegs.lins[i].Rtime_off;
            sgg.Med[contamedia].lumped[0].DiodB = 0.0;
            sgg.Med[contamedia].lumped[0].DiodIsat = 0.0;
            sgg.Med[contamedia].lumped[0].orient = this->DielRegs.lins[i].DiodOri;
         } else if (this->DielRegs.lins[i].inductor) {
            sgg.Med[contamedia].Is.Lumped = true;
            sgg.Med[contamedia].Is.lossy = true; //importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med[contamedia].lumped.resize(1);
            sgg.Med[contamedia].lumped[0].resistor = false;
            sgg.Med[contamedia].lumped[0].inductor = true;
            sgg.Med[contamedia].lumped[0].capacitor = false;
            sgg.Med[contamedia].lumped[0].diodo = false;
            sgg.Med[contamedia].lumped[0].R = this->DielRegs.lins[i].R;
            sgg.Med[contamedia].lumped[0].L = this->DielRegs.lins[i].L;
            sgg.Med[contamedia].lumped[0].C = 0.0;
            sgg.Med[contamedia].lumped[0].R_devia = this->DielRegs.lins[i].R_devia;
            sgg.Med[contamedia].lumped[0].L_devia = this->DielRegs.lins[i].L_devia;
            sgg.Med[contamedia].lumped[0].C_devia = 0.0;
            sgg.Med[contamedia].lumped[0].Rtime_on = 0.0; //irrelevant
            sgg.Med[contamedia].lumped[0].Rtime_off = 0.0; //irrelevant
            sgg.Med[contamedia].lumped[0].DiodB = 0.0;
            sgg.Med[contamedia].lumped[0].DiodIsat = 0.0;
            sgg.Med[contamedia].lumped[0].orient = this->DielRegs.lins[i].DiodOri;
         } else if (this->DielRegs.lins[i].capacitor) {
            sgg.Med[contamedia].Is.Lumped = true;
            sgg.Med[contamedia].Is.lossy = true; //importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med[contamedia].lumped.resize(1);
            sgg.Med[contamedia].lumped[0].resistor = false;
            sgg.Med[contamedia].lumped[0].inductor = false;
            sgg.Med[contamedia].lumped[0].capacitor = true;
            sgg.Med[contamedia].lumped[0].diodo = false;
            sgg.Med[contamedia].lumped[0].R = this->DielRegs.lins[i].R;
            sgg.Med[contamedia].lumped[0].L = 0.0;
            sgg.Med[contamedia].lumped[0].C = this->DielRegs.lins[i].C;
            sgg.Med[contamedia].lumped[0].R_devia = this->DielRegs.lins[i].R_devia;

sgg.Med[contamedia].lumped[0].L_devia = 0.0;
            sgg.Med[contamedia].lumped[0].C_devia = this->DielRegs.lins[i].C_devia;
            sgg.Med[contamedia].lumped[0].Rtime_on = 0.0; //irrelevant
            sgg.Med[contamedia].lumped[0].Rtime_off = 0.0; //irrelevant
            sgg.Med[contamedia].lumped[0].DiodB = 0.0;
            sgg.Med[contamedia].lumped[0].DiodIsat = 0.0;
            sgg.Med[contamedia].lumped[0].orient = this->DielRegs.lins[i].DiodOri;
         } else if (this->DielRegs.lins[i].diodo) {
            //!!!27/08/15 diodos aun no soportados
            std::string buff = "Lumped Diodes currently unsupported. .";
            STOPONERROR(layoutnumber, num_procs, buff);
            //!!!
            sgg.Med[contamedia].Is.Lumped = true;
            sgg.Med[contamedia].Is.lossy = true; //importante que si es lumped esto se ponga a lossy para que thin-wires haga bien el bonding !bug agb 120123 test_GGGbugresis_wire_stoch_foragasconbug
            sgg.Med[contamedia].lumped.resize(1);
            sgg.Med[contamedia].lumped[0].resistor = false;
            sgg.Med[contamedia].lumped[0].inductor = false;
            sgg.Med[contamedia].lumped[0].capacitor = false;
            sgg.Med[contamedia].lumped[0].diodo = true;
            sgg.Med[contamedia].lumped[0].R = this->DielRegs.lins[i].R;
            sgg.Med[contamedia].lumped[0].Rtime_on = 0.0; //irrelevant
            sgg.Med[contamedia].lumped[0].Rtime_off = 0.0; //irrelevant
            sgg.Med[contamedia].lumped[0].L = 0.0;
            sgg.Med[contamedia].lumped[0].C = 0.0;
            sgg.Med[contamedia].lumped[0].DiodB = this->DielRegs.lins[i].DiodB;
            sgg.Med[contamedia].lumped[0].DiodIsat = this->DielRegs.lins[i].DiodIsat;
            sgg.Med[contamedia].lumped[0].orient = this->DielRegs.lins[i].DiodOri;
         } else {
            sgg.Med[contamedia].Is.Lumped = false;
            if (!this->DielRegs.lins[i].plain) {
               std::string buff = "Buggy error 1 in preprocess lumped. .";
               STOPONERROR(layoutnumber, num_procs, buff);
            }
         }
         //!!!fin lumped
         tama2 = (this->DielRegs.lins[i].n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->DielRegs.lins[i].c2P[j].XI;
            punto.XE = this->DielRegs.lins[i].c2P[j].XE;
            punto.YI = this->DielRegs.lins[i].c2P[j].YI;
            punto.YE = this->DielRegs.lins[i].c2P[j].YE;
            punto.ZI = this->DielRegs.lins[i].c2P[j].ZI;
            punto.ZE = this->DielRegs.lins[i].c2P[j].ZE;
            orientacion = this->DielRegs.lins[i].c2P[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->DielRegs.lins[i].c2P[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, &punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
         tama3 = (this->DielRegs.lins[i].n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->DielRegs.lins[i].c1P[j].XI;
            punto.XE = this->DielRegs.lins[i].c1P[j].XI;
            punto.YI = this->DielRegs.lins[i].c1P[j].YI;
            punto.YE = this->DielRegs.lins[i].c1P[j].YI;
            punto.ZI = this->DielRegs.lins[i].c1P[j].ZI;
            punto.ZE = this->DielRegs.lins[i].c1P[j].ZI;
            orientacion = this->DielRegs.lins[i].c1P[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->DielRegs.lins[i].c1P[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, &punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
         }
      }

      //Anisotropic materials
      //materialList
      //BODYes
      tama = (this->ANIMATS.nvols);
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med[contamedia].Anisotropic.resize(1);
         sgg.Med[contamedia].Is.Anisotropic = true;
         sgg.Med[contamedia].Priority = prior_AB;
         sgg.Med[contamedia].Anisotropic[0].Epr = this->ANIMATS.vols[i].eps / Eps0;
         sgg.Med[contamedia].Anisotropic[0].Sigma = this->ANIMATS.vols[i].Sigma;
         sgg.Med[contamedia].Anisotropic[0].Mur = this->ANIMATS.vols[i].mu / Mu0;
         sgg.Med[contamedia].Anisotropic[0].SigmaM = this->ANIMATS.vols[i].SigmaM;
         tama2 = (this->ANIMATS.vols[i].n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->ANIMATS.vols[i].c2P[j].XI;
            punto.XE = this->ANIMATS.vols[i].c2P[j].XE;
            punto.YI = this->ANIMATS.vols[i].c2P[j].YI;
            punto.YE = this->ANIMATS.vols[i].c2P[j].YE;
            punto.ZI = this->ANIMATS.vols[i].c2P[j].ZI;
            punto.ZE = this->ANIMATS.vols[i].c2P[j].ZE;
            numertag = searchtag(tagtype, this->ANIMATS.vols[i].c2P[j].tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, &punto, contamedia);
         }
         tama3 = (this->ANIMATS.vols[i].n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->ANIMATS.vols[i].c1P[j].XI;
            punto.XE = this->ANIMATS.vols[i].c1P[j].XI;
            punto.YI = this->ANIMATS.vols[i].c1P[j].YI;
            punto.YE = this->ANIMATS.vols[i].c1P[j].YI;
            punto.ZI = this->ANIMATS.vols[i].c1P[j].ZI;
            punto.ZE = this->ANIMATS.vols[i].c1P[j].ZI;
            //
            numertag = searchtag(tagtype, this->ANIMATS.vols[i].c1P[j].tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, &punto, contamedia);
         }
      }

& Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
         }
      }
      //SURFs
      tama = (this->ANIMATS.nsurfs);
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med[contamedia].Anisotropic.resize(1);
         sgg.Med[contamedia].Is.Anisotropic = true;
         sgg.Med[contamedia].Priority = prior_IS;
         sgg.Med[contamedia].Anisotropic[0].Epr = this->ANIMATS.surfs[i].eps / Eps0;
         sgg.Med[contamedia].Anisotropic[0].Sigma = this->ANIMATS.surfs[i].Sigma;
         sgg.Med[contamedia].Anisotropic[0].Mur = this->ANIMATS.surfs[i].mu / Mu0;
         sgg.Med[contamedia].Anisotropic[0].SigmaM = this->ANIMATS.surfs[i].SigmaM;
         tama2 = (this->ANIMATS.surfs[i].n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->ANIMATS.surfs[i].c2P[j].XI;
            punto.XE = this->ANIMATS.surfs[i].c2P[j].XE;
            punto.YI = this->ANIMATS.surfs[i].c2P[j].YI;
            punto.YE = this->ANIMATS.surfs[i].c2P[j].YE;
            punto.ZI = this->ANIMATS.surfs[i].c2P[j].ZI;
            punto.ZE = this->ANIMATS.surfs[i].c2P[j].ZE;
            orientacion = this->ANIMATS.surfs[i].c2P[j].or;
            numertag = searchtag(tagtype, this->ANIMATS.surfs[i].c2P[j].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz, &
            & media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, &
            & contamedia);
         }
         tama3 = (this->ANIMATS.surfs[i].n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->ANIMATS.surfs[i].c1P[j].XI;
            punto.XE = this->ANIMATS.surfs[i].c1P[j].XI;
            punto.YI = this->ANIMATS.surfs[i].c1P[j].YI;
            punto.YE = this->ANIMATS.surfs[i].c1P[j].YI;
            punto.ZI = this->ANIMATS.surfs[i].c1P[j].ZI;
            punto.ZE = this->ANIMATS.surfs[i].c1P[j].ZI;
            orientacion = this->ANIMATS.surfs[i].c1P[j].or;
            numertag = searchtag(tagtype, this->ANIMATS.surfs[i].c1P[j].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz, &
            & media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, &
            & contamedia);
         }
      }
      //LINs
      tama = (this->ANIMATS.nLINS);
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         sgg.Med[contamedia].Anisotropic.resize(1);
         sgg.Med[contamedia].Is.Anisotropic = true;
         sgg.Med[contamedia].Priority = prior_IL;
         sgg.Med[contamedia].Anisotropic[0].Epr = this->ANIMATS.lins[i].eps / Eps0;
         sgg.Med[contamedia].Anisotropic[0].Sigma = this->ANIMATS.lins[i].Sigma;
         sgg.Med[contamedia].Anisotropic[0].Mur = this->ANIMATS.lins[i].mu / Mu0;
         sgg.Med[contamedia].Anisotropic[0].SigmaM = this->ANIMATS.lins[i].SigmaM;
         tama2 = (this->ANIMATS.lins[i].n_c2P);
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->ANIMATS.lins[i].c2P[j].XI;
            punto.XE = this->ANIMATS.lins[i].c2P[j].XE;
            punto.YI = this->ANIMATS.lins[i].c2P[j].YI;
            punto.YE = this->ANIMATS.lins[i].c2P[j].YE;
            punto.ZI = this->ANIMATS.lins[i].c2P[j].ZI;
            punto.ZE = this->ANIMATS.lins[i].c2P[j].ZE;
            orientacion = this->ANIMATS.lins[i].c2P[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->ANIMATS.lins[i].c2P[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz, &
            & media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire, verbose, numeroasignaciones);
         }
         tama3 = (this->ANIMATS.lins[i].n_c1P);
         for (j = 1; j <= tama3; ++j) {
            punto.XI = this->ANIMATS.lins[i].c1P[j].XI;
            punto.XE = this->ANIMATS.lins[i].c1P[j].XI;
            punto.YI = this->ANIMATS.lins[i].c1P[j].YI;
            punto.YE = this->ANIMATS.lins[i].c1P[j].YI;
            punto.ZI = this->ANIMATS.lins[i].c1P[j].ZI;
            punto.ZE = this->ANIMATS.lins[i].c1P[j].ZI;
            orientacion = this->ANIMATS.lins[i].c1P[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->ANIMATS.lins[i].c1P[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz, &
            & media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI, &
            & Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI, &
            & Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI, &
            & Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI, &
            & Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI, &
            & Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion, &
            & contamedia, isathinwire, verbose, numeroasignaciones);
         }
      }
      //frequency dependent materials
      //bodies
      //
      tama = this->FRQDEPMATS.nvols;
      for (i = 1; i <= tama; ++i) {
         contamedia = contamedia + 1;
         fdgeom.reset();
         fdgeom = this->FRQDEPMATS.vols[i];
         asignadisper(fdgeom);
         //geometry
         //!!!!!!!!

         tama2 = this->FRQDEPMATS.vols[i].n_C;
         for (j = 1; j <= tama2; ++j) {
            punto.XI = this->FRQDEPMATS.vols[i].C[j].XI;
            punto.XE = this->FRQDEPMATS.vols[i].C[j].XE;
            punto.YI = this->FRQDEPMATS.vols[i].C[j].YI;
            punto.YE = this->FRQDEPMATS.vols[i].C[j].YE;
            punto.ZI = this->FRQDEPMATS.vols[i].C[j].ZI;
            punto.ZE = this->FRQDEPMATS.vols[i].C[j].ZE;

numertag = searchtag(tagtype, this->FRQDEPMATS.vols[i].C[j].tag);
            CreateVolumeMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                           media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                           Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                           Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                           Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                           Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                           Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, contamedia);
        }
    }

    // SURFs
    tama = this->FRQDEPMATS.nsurfs;
    for (i = 1; i <= tama; ++i) {
        contamedia = contamedia + 1;
        fdgeom = nullptr;
        fdgeom = &this->FRQDEPMATS.surfs[i];
        asignadisper(fdgeom);

        tama2 = this->FRQDEPMATS.surfs[i].n_C;
        for (j = 1; j <= tama2; ++j) {
            punto.XI = this->FRQDEPMATS.surfs[i].C[j].XI;
            punto.XE = this->FRQDEPMATS.surfs[i].C[j].XE;
            punto.YI = this->FRQDEPMATS.surfs[i].C[j].YI;
            punto.YE = this->FRQDEPMATS.surfs[i].C[j].YE;
            punto.ZI = this->FRQDEPMATS.surfs[i].C[j].ZI;
            punto.ZE = this->FRQDEPMATS.surfs[i].C[j].ZE;
            orientacion = this->FRQDEPMATS.surfs[i].C[j].or;
            numertag = searchtag(tagtype, this->FRQDEPMATS.surfs[i].C[j].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
        }
    }

    // LINs
    tama = this->FRQDEPMATS.nLINS;
    for (i = 1; i <= tama; ++i) {
        contamedia = contamedia + 1;
        fdgeom = nullptr;
        fdgeom = &this->FRQDEPMATS.lins[i];
        asignadisper(fdgeom);

        tama2 = this->FRQDEPMATS.lins[i].n_C;
        for (j = 1; j <= tama2; ++j) {
            punto.XI = this->FRQDEPMATS.lins[i].C[j].XI;
            punto.XE = this->FRQDEPMATS.lins[i].C[j].XE;
            punto.YI = this->FRQDEPMATS.lins[i].C[j].YI;
            punto.YE = this->FRQDEPMATS.lins[i].C[j].YE;
            punto.ZI = this->FRQDEPMATS.lins[i].C[j].ZI;
            punto.ZE = this->FRQDEPMATS.lins[i].C[j].ZE;
            orientacion = this->FRQDEPMATS.lins[i].C[j].or;
            isathinwire = false;
            numertag = searchtag(tagtype, this->FRQDEPMATS.lins[i].C[j].tag);
            CreateLineMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                         media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                         Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                         Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                         Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                         Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                         Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                         contamedia, isathinwire, verbose, numeroasignaciones);
        }
    }

    //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // ISOTROPIC Multiports
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    inicontamedia = contamedia + 1;
    maxcontamedia = contamedia;
    tama = this->LossyThinSurfs.length;
    for (j = 1; j <= tama; ++j) {
        // carbon Multiports a guevo
        if (this->LossyThinSurfs.cs[j].numcapas == 0) {
            this->LossyThinSurfs.cs[j].numcapas = 1;
            this->LossyThinSurfs.cs[j].SigmaM.resize(1);
            this->LossyThinSurfs.cs[j].Sigma.resize(1);
            this->LossyThinSurfs.cs[j].EPS.resize(1);
            this->LossyThinSurfs.cs[j].MU.resize(1);
            this->LossyThinSurfs.cs[j].thk.resize(1);

            // _for_devia 090519
            this->LossyThinSurfs.cs[j].SigmaM_devia.resize(1);
            this->LossyThinSurfs.cs[j].Sigma_devia.resize(1);
            this->LossyThinSurfs.cs[j].EPS_devia.resize(1);
            this->LossyThinSurfs.cs[j].MU_devia.resize(1);
            this->LossyThinSurfs.cs[j].thk_devia.resize(1);
            // !!!

            this->LossyThinSurfs.cs[j].SigmaM[0] = 0.0; // TRUCO PARA QUE CUANDO NO TENGA CAPAS (LECTURA DESDE FICHERO DE POLOS /RESIDUOS) NO PETE
            this->LossyThinSurfs.cs[j].Sigma[0] = 0.0;
            this->LossyThinSurfs.cs[j].EPS[0] = EPS0;
            this->LossyThinSurfs.cs[j].MU[0] = MU0;
            this->LossyThinSurfs.cs[j].thk[0] = -1.0;

            // _for_devia 090519
            this->LossyThinSurfs.cs[j].SigmaM_devia[0] = 0.0;
            this->LossyThinSurfs.cs[j].Sigma_devia[0] = 0.0;
            this->LossyThinSurfs.cs[j].EPS_devia[0] = 0.0;
            this->LossyThinSurfs.cs[j].MU_devia[0] = 0.0;
            this->LossyThinSurfs.cs[j].thk_devia[0] = 0.0;
            // !!!
            // !!!comentado el 120219 pq no se lleva bien con Semba !no entiendo ahora el comentario de malonyedispersive!!!120219
            // !!!     write(buff, '(a)')    'pre1_Error:  SGBC materials must have at least one layyer even in dummy for malonyedispersive'
            // !!!     call WarnErrReport (buff,.true.)
        }

        if (std::abs(this->LossyThinSurfs.cs[j].SigmaM[0]) <= 1.0e-2) { // !!!ojoooo a 210319 manda guevos que tengamos que estar con el flag de la conductidad magnetica para llamar a SGBC todavia en 2015!!!
            this->LossyThinSurfs.cs[j].SigmaM[0] = 0.0;
            if (!mibc) {
                // if (this->LossyThinSurfs.cs[j].numcapas >1) then
                //    write(buff, '(a)')    'pre1_Warning:  SGBC materials are just averaged for multilayered structures.'// &
                //    ' Use preferably -mibc instead.'
                //    call WarnErrReport (buff)
                // end if
                SGBC = true; // si la conductividad es 0.0 (o casi) utiliza directamente SGBC
                mibc = false;
            }
        }
        if (this->LossyThinSurfs.cs[j].SigmaM[0] >= 0.0) {
            // SURFs (siempre son surfs)
            //
            tama2 = this->LossyThinSurfs.cs[j].nc;
            mincontamedia = maxcontamedia + 1;
            MultiportFile = trim(adjustl(this->LossyThinSurfs.cs[j].files)) + "_z11.txt";
            for (i = 1; i <= tama2; ++i) {
                orientacion = this->LossyThinSurfs.cs[j].C[i].or;
                punto.XI = this->LossyThinSurfs.cs[j].C[i].XI;
                punto.XE = this->LossyThinSurfs.cs[j].C[i].XE;
                punto.YI = this->LossyThinSurfs.cs[j].C[i].YI;
                punto.YE = this->LossyThinSurfs.cs[j].C[i].YE;
                punto.ZI = this->LossyThinSurfs.cs[j].C[i].ZI;
                punto.ZE = this->LossyThinSurfs.cs[j].C[i].ZE;
                existia = false;
                for (k = inicontamedia; k <= maxcontamedia; ++k) {

if (trim(adjustl(sgg.Med[k].multiport[1].multiportFileZ11)) == trim(adjustl(MultiportFile))) {
                     if (sgg.Med[k].multiport[1].Multiportdir == orientacion) {
                        contamedia = k;
                        existia = true;
                        goto doexis_exit;
                     }
                  }
               }
doexis_exit:;

               if (!existia) {
                  maxcontamedia = maxcontamedia + 1;
                  contamedia = maxcontamedia;
                  sgg.Med.resize(contamedia + 1);
                  sgg.Med[contamedia].multiport.resize(2);
                  //
                  if ((this->LossyThinSurfs.cs[j].numcapas > 1) && SGBCDispersive) {
                     write(buff, *) 'ERROR in SGBCs Number of layers >1 still unsupported for SGBCDispersive. ';
                     StopOnError(0, 0, buff);
                  }
                  //
                  sgg.Med[contamedia].Multiport.resize(2);
                  sgg.Med[contamedia].Multiport[1].epr.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].mur.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].sigma.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].sigmam.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].width.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  // _for_devia 090519
                  sgg.Med[contamedia].Multiport[1].epr_devia.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].mur_devia.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].sigma_devia.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].sigmaM_devia.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  sgg.Med[contamedia].Multiport[1].width_devia.resize(this->LossyThinSurfs.cs[j].numcapas + 1);
                  // !!!
                  puntoXI = Max(punto.XI, Min(BoundingBox.XI, BoundingBox.XE));
                  puntoYI = Max(punto.YI, Min(BoundingBox.YI, BoundingBox.YE));
                  puntoZI = Max(punto.ZI, Min(BoundingBox.ZI, BoundingBox.ZE));

                  // !!!!!!!!estaba antes  maaaal. bug 140815verano
                  if (!((puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE))) {
                     puntoXI = sgg.allocDxI;
                     write(buff, *) 'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)', puntoXI, puntoYI, puntoZI, sgg.allocDxI, sgg.allocDyI, sgg.allocDzI;
                     WarnErrReport(buff, true);
                  }
                  if (!((puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE))) {
                     puntoYI = sgg.allocDyI;
                     write(buff, *) 'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)', puntoXI, puntoYI, puntoZI, sgg.allocDxI, sgg.allocDyI, sgg.allocDzI;
                     WarnErrReport(buff, true);
                  }
                  if (!((puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE))) {
                     puntoZI = sgg.allocDzI;
                     write(buff, *) 'ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)', puntoXI, puntoYI, puntoZI, sgg.allocDxI, sgg.allocDyI, sgg.allocDzI;
                     WarnErrReport(buff, true);
                  }
                  dentro = (puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE) &&
                     (puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE) &&
                     (puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE);
                  delta = -1.0;
                  if (DENTRO) {
                     switch (abs(this->LossyThinSurfs.cs[j].C[i].or)) {
                      case iEx:
                        delta = (sgg.DX[puntoXI] + sgg.DX[puntoXI - 1]) / 2.0;
                        break;
                      case iEy:
                        delta = (sgg.DY[puntoYI] + sgg.Dy[puntoYI - 1]) / 2.0;
                        break;
                      case iEz:
                        delta = (sgg.DZ[puntoZI] + sgg.Dz[puntoZI - 1]) / 2.0;
                        break;
                      default:
                        write(buff, '(a)') 'Buggy error 1 in preprocess composites. .';
                        STOPONERROR(layoutnumber, num_procs, buff);
                     }
                  } else {
                     write(buff, '(a)') 'Buggy error 2 in preprocess composites. .';
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }
                  sgg.Med[contamedia].Multiport[1].numcapas = this->LossyThinSurfs.cs[j].numcapas;
                  // el especificado
                  sgg.Med[contamedia].Multiport[1].Multiportdir = this->LossyThinSurfs.cs[j].C[i].or;
                  for (I_ = 1; I_ <= sgg.Med[contamedia].Multiport[1].numcapas; I_++) {
                     if (sgg.Med[contamedia].Multiport[1].Multiportdir > 0) {
                        j_ = I_;
                     } else {
                        j_ = sgg.Med[contamedia].Multiport[1].numcapas - I_ + 1; // dale la vuelta (medios no simetricos) !0121
                     }
                     sgg.Med[contamedia].Multiport[1].epr[j_] = this->LossyThinSurfs.cs[j].eps[I_] / Eps0;
                     sgg.Med[contamedia].Multiport[1].mur[j_] = this->LossyThinSurfs.cs[j].mu[I_] / mu0;
                     sgg.Med[contamedia].Multiport[1].sigma[j_] = this->LossyThinSurfs.cs[j].Sigma[I_];
                     sgg.Med[contamedia].Multiport[1].sigmam[j_] = abs(this->LossyThinSurfs.cs[j].Sigmam[I_]);
                     sgg.Med[contamedia].Multiport[1].width[j_] = this->LossyThinSurfs.cs[j].thk[I_];

                     // _for_devia 090519
                     sgg.Med[contamedia].Multiport[1].epr_devia[j_] = this->LossyThinSurfs.cs[j].eps_devia[I_] / Eps0;
                     sgg.Med[contamedia].Multiport[1].mur_devia[j_] = this->LossyThinSurfs.cs[j].MU_devia[I_] / mu0;
                     sgg.Med[contamedia].Multiport[1].sigma_devia[j_] = this->LossyThinSurfs.cs[j].Sigma_devia[I_];
                     sgg.Med[contamedia].Multiport[1].sigmaM_devia[j_] = abs(this->LossyThinSurfs.cs[j].SigmaM_devia[I_]);
                     sgg.Med[contamedia].Multiport[1].width_devia[j_] = this->LossyThinSurfs.cs[j].thk_devia[I_];
                  }

                  rdummy = maxval(abs(this->LossyThinSurfs.cs[j].MU_devia)) + maxval(abs(this->LossyThinSurfs.cs[j].SigmaM_devia));
                  if (rdummy > 1.0e-15) {
                     write(buff, '(a)') 'Non null deviations found in sigmam or mu in composites. Still unsupported.';
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }

                  // !!!

                  // !!old pre 17/07/15
                  // !!!                           sgg.Med[contamedia].Multiport[1].transversalSpaceDelta = delta;
                  sgg.Med[contamedia].Priority = prior_CS;
                  sgg.Med[contamedia].Epr = this->LossyThinSurfs.cs[j].eps[1] / Eps0;
                  sgg.Med[contamedia].Sigma = this->LossyThinSurfs.cs[j].Sigma[1];
                  sgg.Med[contamedia].Mur = this->LossyThinSurfs.cs[j].mu[1] / Mu0;
                  sgg.Med[contamedia].SigmaM = abs(this->LossyThinSurfs.cs[j].SigmaM[1]); // may be negative
                  if (mibc) {
                     sgg.Med[contamedia].Is.multiport = true;
                     sgg.Med[contamedia].Is.Lossy = true;
                  } else if (SGBC) {
                     sgg.Med[contamedia].Is.SGBC = true;
                     sgg.Med[contamedia].Is.Lossy = true;
                  }

if (SGBCDispersive) sgg.Med(contamedia).Is.SGBCDispersive = true;
                  else {
                     write(buff, '(a)') "Some -mibc -sgbc switch should be used for Composites.";
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }
                  sgg.Med(contamedia).Is.Dielectric = false;

                  sgg.Med(contamedia).multiport(1).multiportFileZ11 = trim(adjustl(this.LossyThinSurfs.cs(j).files)) + "_z11.txt";
                  sgg.Med(contamedia).multiport(1).multiportFileZ22 = trim(adjustl(this.LossyThinSurfs.cs(j).files)) + "_z22.txt";
                  sgg.Med(contamedia).multiport(1).multiportFileZ12 = trim(adjustl(this.LossyThinSurfs.cs(j).files)) + "_z12.txt";
                  sgg.Med(contamedia).multiport(1).multiportFileZ21 = trim(adjustl(this.LossyThinSurfs.cs(j).files)) + "_z12.txt";

                  if (mibc) {
                     sgg.Med(contamedia).Is.SGBC = false;
                     sgg.Med(contamedia).Is.SGBCDispersive = false;
                     sgg.Med(contamedia).Is.Lossy = true;

                     errnofile1 = false;
                     errnofile2 = false;
                     errnofile3 = false;
                     errnofile4 = false;
                     pozi = index(sgg.Med(contamedia).Multiport(1).multiportFileZ11, "_z11.txt");
                     multiportfile2 = trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ11.substr(0, pozi - 1)));
                     inquire(file=trim(adjustl(multiportfile2)), EXIST=errnofile1);
                     inquire(file=trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ11)), EXIST=errnofile2);
                     inquire(file=trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ12)), EXIST=errnofile3);
                     inquire(file=trim(adjustl(sgg.Med(contamedia).Multiport(1).multiportFileZ22)), EXIST=errnofile4);

                     if (!(errnofile1 || (errnofile2 && errnofile3 && errnofile4))) {
                        buff = "Neither New nor Old style mibc FILE " + trim(adjustl(multiportfile2)) + " EXISTS.";
                        WarnErrReport(buff, true);
                     }
                  }

               }
               numertag = searchtag(tagtype, this.LossyThinSurfs.cs(j).C(i).tag);
               CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                               media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                               Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                               Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                               Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                               Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                               Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                               contamedia);
            }
         }
      }
      contamedia = maxcontamedia;

      inicontamedia = contamedia + 1;
      maxcontamedia = contamedia;
      tama = this.LossyThinSurfs.length;
      for (j = 1; j <= tama; ++j) {
         if (this.LossyThinSurfs.cs(j).SigmaM(1) < 0.0_RKIND) {
            tama2 = this.LossyThinSurfs.cs(j).nc;
            mincontamedia = maxcontamedia + 1;
            MultiportFile = trim(adjustl(this.LossyThinSurfs.cs(j).files)) + "_z11.txt";
            for (i = 1; i <= tama2; ++i) {
               orientacion = this.LossyThinSurfs.cs(j).C(i).or;
               punto.XI = this.LossyThinSurfs.cs(j).C(i).XI;
               punto.XE = this.LossyThinSurfs.cs(j).C(i).XE;
               punto.YI = this.LossyThinSurfs.cs(j).C(i).YI;
               punto.YE = this.LossyThinSurfs.cs(j).C(i).YE;
               punto.ZI = this.LossyThinSurfs.cs(j).C(i).ZI;
               punto.ZE = this.LossyThinSurfs.cs(j).C(i).ZE;
               existia = false;
               for (k = inicontamedia; k <= maxcontamedia; ++k) {
                  if (trim(adjustl(sgg.Med(k).AnisMultiport(1).multiportFileZ11)) == trim(adjustl(MultiportFile))) {
                     if (sgg.Med(k).AnisMultiport(1).Multiportdir == orientacion) {
                        contamedia = k;
                        existia = true;
                        break;
                     }
                  }
               }
               if (!existia) {
                  maxcontamedia = maxcontamedia + 1;
                  contamedia = maxcontamedia;
                  sgg.Med(contamedia).AnisMultiport.resize(1);
                  sgg.Med(contamedia).AnisMultiport[0] = Multiport(); // Assuming default constructor or appropriate init

                  if (this.LossyThinSurfs.cs(j).numcapas > 1) {
                     write(buff, '(a)') "pre1_ERROR:  Anisotropic multiport materials unsupported for multilayered structures.";
                     WarnErrReport(buff, true);
                  }
                  puntoXI = Max(punto.XI, Min(BoundingBox.XI, BoundingBox.XE));
                  puntoYI = Max(punto.YI, Min(BoundingBox.YI, BoundingBox.YE));
                  puntoZI = Max(punto.ZI, Min(BoundingBox.ZI, BoundingBox.ZE));

                  if (!((puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE))) puntoXI = sgg.allocDxI;
                  if (!((puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE))) puntoYI = sgg.allocDyI;
                  if (!((puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE))) puntoZI = sgg.allocDzI;
                  write(buff, '(a)') "ERROR: precompo 2: Readjusting composite init point. Only ignore if parts of the geometry fall out of the the domain deliberately (only if manual clipping)";
                  WarnErrReport(buff, true);

                  dentro = (puntoXI >= sgg.allocDxI) && (puntoXI <= sgg.allocDxE) &&
                           (puntoYI >= sgg.allocDyI) && (puntoYI <= sgg.allocDyE) &&
                           (puntoZI >= sgg.allocDzI) && (puntoZI <= sgg.allocDzE);
                  delta = -1.0_RKIND;
                  if (dentro) {
                     switch (abs(this.LossyThinSurfs.cs(j).C(i).or)) {
                        case iEx:
                           delta = (sgg.DX(puntoXI) + sgg.DX(puntoXI - 1)) / 2.0_RKIND;
                           break;
                        case iEy:
                           delta = (sgg.DY(puntoYI) + sgg.Dy(puntoYI - 1)) / 2.0_RKIND;
                           break;
                        case iEz:
                           delta = (sgg.DZ(puntoZI) + sgg.Dz(puntoZI - 1)) / 2.0_RKIND;
                           break;
                        default:
                           write(buff, '(a)') "Buggy error 1 in preprocess composites. .";
                           STOPONERROR(layoutnumber, num_procs, buff);
                     }
                  } else {
                     write(buff, '(a)') "Buggy error 2 in preprocess composites. .";
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }
                  sgg.Med(contamedia).AnisMultiport(1).Multiportdir = this.LossyThinSurfs.cs(j).C(i).or;
                  sgg.Med(contamedia).AnisMultiport(1).epr = this.LossyThinSurfs.cs(j).eps / Eps0;

sgg.Med(contamedia).AnisMultiport[1].mur = this->LossyThinSurfs.cs[j].mu / mu0;
                  sgg.Med(contamedia).AnisMultiport[1].sigma = this->LossyThinSurfs.cs[j].Sigma;
                  sgg.Med(contamedia).AnisMultiport[1].sigmam = std::abs(this->LossyThinSurfs.cs[j].Sigmam);
                  sgg.Med(contamedia).AnisMultiport[1].width = this->LossyThinSurfs.cs[j].thk;
                  sgg.Med(contamedia).Priority = prior_CS;
                  sgg.Med(contamedia).Epr = this->LossyThinSurfs.cs[j].eps[0] / Eps0;
                  sgg.Med(contamedia).Sigma = this->LossyThinSurfs.cs[j].Sigma[0];
                  sgg.Med(contamedia).Mur = this->LossyThinSurfs.cs[j].mu[0] / Mu0;
                  sgg.Med(contamedia).SigmaM = std::abs(this->LossyThinSurfs.cs[j].SigmaM[0]); //may be negative

                  if (mibc) {
                     sgg.Med(contamedia).Is.Anismultiport = true;
                     sgg.Med(contamedia).Is.Lossy = true;
                  } else {
                     std::string buff = "Some -mibc -sgbc switch should be used for Anisotropic Composites.";
                     STOPONERROR(layoutnumber, num_procs, buff);
                  }

                  sgg.Med(contamedia).Is.Dielectric = false;

                  sgg.Med(contamedia).AnisMultiport[1].multiportFileZ11 = trim(adjustl(this->LossyThinSurfs.cs[j].files)) + "_z11.txt";
                  sgg.Med(contamedia).AnisMultiport[1].multiportFileZ22 = trim(adjustl(this->LossyThinSurfs.cs[j].files)) + "_z22.txt";
                  sgg.Med(contamedia).AnisMultiport[1].multiportFileZ12 = trim(adjustl(this->LossyThinSurfs.cs[j].files)) + "_z12.txt";
                  sgg.Med(contamedia).AnisMultiport[1].multiportFileZ21 = trim(adjustl(this->LossyThinSurfs.cs[j].files)) + "_z12.txt";
               }
            }

            numertag = searchtag(tagtype, this->LossyThinSurfs.cs[j].C[i].tag);
            CreateSurfaceMM(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy, media.sggMiEz,
                            media.sggMiHx, media.sggMiHy, media.sggMiHz, Alloc_iEx_XI,
                            Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                            Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                            Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                            Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                            Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared, BoundingBox, punto, orientacion,
                            contamedia);
         }
      }
   }
   //Multiports
}
contamedia = maxcontamedia;
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//end ANISOTROPIC multiports
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

//
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Multiports lossy padding
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if (std::abs(attfactor - 1.0) > 1.0e-12) {
   tama = this->LossyThinSurfs.length;
   for (j = 1; j <= tama; ++j) {
      tama2 = this->LossyThinSurfs.cs[j].nc;
      contamedia = contamedia + 1;
      for (i = 1; i <= tama2; ++i) {
         orientacion = this->LossyThinSurfs.cs[j].C[i].or;
         //ES UN free-space multiportpadding CON LA PRIORIDAD DE UN MULTIPORT con conductividad magnetica que luego se desanulara
         sgg.Med(contamedia).Priority = prior_CS;
         sgg.Med(contamedia).Is.multiport = false;
         sgg.Med(contamedia).Is.ANISmultiport = false;
         sgg.Med(contamedia).Is.MultiportPadding = true;
         sgg.Med(contamedia).Is.Lossy = true;
         sgg.Med(contamedia).Is.Dielectric = false;
         sgg.Med(contamedia).Epr = 1.0; //this->LossyThinSurfs.cs[j].eps / Eps0
         sgg.Med(contamedia).Sigma = 0.0; //abs(this->LossyThinSurfs.cs[j].Sigma) !may be negative
         sgg.Med(contamedia).Mur = 1.0; //this->LossyThinSurfs.cs[j].mu / Mu0
         sgg.Med(contamedia).Is.Dielectric = true;
         //provisionalmente (luego se retocara con el sigmam correcto)
         sgg.Med(contamedia).SigmaM = 0.0; //abs(this->LossyThinSurfs.cs[j].Sigma) !may be negative

         punto.XI = this->LossyThinSurfs.cs[j].C[i].XI;
         punto.XE = this->LossyThinSurfs.cs[j].C[i].XE;
         punto.YI = this->LossyThinSurfs.cs[j].C[i].YI;
         punto.YE = this->LossyThinSurfs.cs[j].C[i].YE;
         punto.ZI = this->LossyThinSurfs.cs[j].C[i].ZI;
         punto.ZE = this->LossyThinSurfs.cs[j].C[i].ZE;
         //!!

         switch (std::abs(orientacion)) {
         case iEx:
            punto.XI = this->LossyThinSurfs.cs[j].C[i].XI;
            punto.XE = this->LossyThinSurfs.cs[j].C[i].XI;
            numertag = searchtag(tagtype, this->LossyThinSurfs.cs[j].C[i].tag);
            CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                                  media.sggMiEz,
                                  media.sggMiHx, media.sggMiHy, media.sggMiHz,
                                  Alloc_iEx_XI,
                                  Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE,
                                  Alloc_iEy_XI, Alloc_iEy_XE, Alloc_iEy_YI,
                                  Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE,
                                  Alloc_iEz_YI, Alloc_iEz_YE, Alloc_iEz_ZI,
                                  Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE,
                                  Alloc_iHx_ZI, Alloc_iHx_ZE, Alloc_iHy_XI,
                                  Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE,
                                  Alloc_iHz_XI, Alloc_iHz_XE, Alloc_iHz_YI,
                                  Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia,
                                  sgg.EShared, BoundingBox, punto, orientacion,
                                  contamedia);
            punto.XI = this->LossyThinSurfs.cs[j].C[i].XI - 1;
            punto.XE = this->LossyThinSurfs.cs[j].C[i].XI - 1;
            numertag = searchtag(tagtype, this->LossyThinSurfs.cs[j].C[i].tag);
            CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                                  media.sggMiEz,
                                  media.sggMiHx, media.sggMiHy, media.sggMiHz,
                                  Alloc_iEx_XI,
                                  Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,
                                  Alloc_iEy_XE, Alloc_iEy_YI,
                                  Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                                  Alloc_iEz_YE, Alloc_iEz_ZI,
                                  Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                                  Alloc_iHx_ZE, Alloc_iHy_XI,
                                  Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                                  Alloc_iHz_XE, Alloc_iHz_YI,
                                  Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                                  BoundingBox, punto, orientacion,
                                  contamedia);
            break;
         case iEy:
            punto.YI = this->LossyThinSurfs.cs[j].C[i].YI;
            punto.YE = this->LossyThinSurfs.cs[j].C[i].YI;
            numertag = searchtag(tagtype, this->LossyThinSurfs.cs[j].C[i].tag);

CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                    media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz,
                    Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,
                    Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                    Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                    Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                    Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                    BoundingBox, punto, orientacion,
                    contamedia);
                punto.YI = this->LossyThinSurfs.cs(j).C(i).YI - 1;
                punto.YE = this->LossyThinSurfs.cs(j).C(i).YI - 1;
                numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                    media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz,
                    Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,
                    Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                    Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                    Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                    Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                    BoundingBox, punto, orientacion,
                    contamedia);
                break;
            case iEz:
                punto.ZI = this->LossyThinSurfs.cs(j).C(i).ZI;
                punto.ZE = this->LossyThinSurfs.cs(j).C(i).ZI;
                numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                    media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz,
                    Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,
                    Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                    Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                    Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                    Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                    BoundingBox, punto, orientacion,
                    contamedia);
                punto.ZI = this->LossyThinSurfs.cs(j).C(i).ZI - 1;
                punto.ZE = this->LossyThinSurfs.cs(j).C(i).ZI - 1;
                numertag = searchtag(tagtype, this->LossyThinSurfs.cs(j).C(i).tag);
                CreateMagneticSurface(layoutnumber, media.sggMtag, tag_numbers, numertag, media.sggMiEx, media.sggMiEy,
                    media.sggMiEz,
                    media.sggMiHx, media.sggMiHy, media.sggMiHz,
                    Alloc_iEx_XI,
                    Alloc_iEx_XE, Alloc_iEx_YI, Alloc_iEx_YE, Alloc_iEx_ZI, Alloc_iEx_ZE, Alloc_iEy_XI,
                    Alloc_iEy_XE, Alloc_iEy_YI,
                    Alloc_iEy_YE, Alloc_iEy_ZI, Alloc_iEy_ZE, Alloc_iEz_XI, Alloc_iEz_XE, Alloc_iEz_YI,
                    Alloc_iEz_YE, Alloc_iEz_ZI,
                    Alloc_iEz_ZE, Alloc_iHx_XI, Alloc_iHx_XE, Alloc_iHx_YI, Alloc_iHx_YE, Alloc_iHx_ZI,
                    Alloc_iHx_ZE, Alloc_iHy_XI,
                    Alloc_iHy_XE, Alloc_iHy_YI, Alloc_iHy_YE, Alloc_iHy_ZI, Alloc_iHy_ZE, Alloc_iHz_XI,
                    Alloc_iHz_XE, Alloc_iHz_YI,
                    Alloc_iHz_YE, Alloc_iHz_ZI, Alloc_iHz_ZE, sgg.Med, sgg.NumMedia, sgg.EShared,
                    BoundingBox, punto, orientacion,
                    contamedia);
                break;
            }
        }
    }
}
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// wires
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
oldcontamedia = contamedIA;
tama = this->twires.n_tw;
for (j = 1; j <= tama; ++j) {
    contamedia = contamedia + 1;
    sgg.Med.resize(contamedia);
    sgg.Med[contamedia - 1].wire.resize(1);
    sgg.Med[contamedia - 1].Priority = prior_TW;

    //
    // background
    //
    sgg.Med[contamedia - 1].Epr = sgg.Med[0].Epr;
    sgg.Med[contamedia - 1].Sigma = sgg.Med[0].Sigma;
    sgg.Med[contamedia - 1].Mur = sgg.Med[0].Mur;
    sgg.Med[contamedia - 1].SigmaM = sgg.Med[0].SigmaM;
    sgg.Med[contamedia - 1].Is.ThinWire = true;
    sgg.Med[contamedia - 1].Is.Dielectric = false;
    sgg.Med[contamedia - 1].wire[0].radius = this->twires.TW[j - 1].RAD;
    sgg.Med[contamedia - 1].wire[0].radius_devia = this->twires.TW[j - 1].RAD_devia;
    if (boundwireradius) {
        if (sgg.Med[contamedia - 1].wire[0].radius > maxwireradius) sgg.Med[contamedia - 1].wire[0].radius = maxwireradius;
    }
    sgg.Med[contamedia - 1].wire[0].R = this->twires.TW[j - 1].RES;
    sgg.Med[contamedia - 1].wire[0].l = this->twires.TW[j - 1].IND;
    sgg.Med[contamedia - 1].wire[0].C = this->twires.TW[j - 1].CAP;
    sgg.Med[contamedia - 1].wire[0].P_R = this->twires.TW[j - 1].P_RES;
    sgg.Med[contamedia - 1].wire[0].P_l = this->twires.TW[j - 1].P_IND;

    sgg.Med[contamedia - 1].wire[0].P_C = this->twires.TW[j - 1].P_CAP;
    if (this->twires.TW[j - 1].disp) {
        sgg.Med[contamedia - 1].wire[0].disp.resize(1);
        asignawiredisper(sgg.Med[contamedia - 1].wire[0].disp[0],
            this->twires.TW[j - 1].dispfile);
    }
    sgg.Med[contamedia - 1].wire[0].LeftEnd = this->twires.TW[j - 1].LeftEnd;
    sgg.Med[contamedia - 1].wire[0].RightEnd = this->twires.TW[j - 1].RightEnd;
    sgg.Med[contamedia - 1].wire[0].VsourceExists = false;
    sgg.Med[contamedia - 1].wire[0].IsourceExists = false;
    //
    sgg.Med[contamedia - 1].wire[0].HasAbsorbing_RightEnd = false;
    sgg.Med[contamedia - 1].wire[0].HasAbsorbing_LeftEnd = false;
    sgg.Med[contamedia - 1].wire[0].HasParallel_RightEnd = false;
    sgg.Med[contamedia - 1].wire[0].HasParallel_LeftEnd = false;
    sgg.Med[contamedia - 1].wire[0].HasSeries_RightEnd = false;
    sgg.Med[contamedia - 1].wire[0].HasSeries_LeftEnd = false;
    sgg.Med[contamedia - 1].wire[0].Parallel_R_RightEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Parallel_R_LeftEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Series_R_RightEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Series_R_LeftEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Parallel_L_RightEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Parallel_L_LeftEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Series_L_RightEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Series_L_LeftEnd = 0.0;
    sgg.Med[contamedia - 1].wire[0].Parallel_C_RightEnd = 0.0;
}

sgg.Med(contamedia).wire[1].Parallel_C_LeftEnd = 0.0;
            sgg.Med(contamedia).wire[1].Series_C_RightEnd = 2.0e7; //en corto 14/2/14
            sgg.Med(contamedia).wire[1].Series_C_LeftEnd = 2.0e7; //en corto 14/2/14
//stoch
            sgg.Med(contamedia).wire[1].R_devia = this.twires.TW[j].RES_devia;
            sgg.Med(contamedia).wire[1].l_devia = this.twires.TW[j].IND_devia;
            sgg.Med(contamedia).wire[1].C_devia = this.twires.TW[j].CAP_devia;

            sgg.Med(contamedia).wire[1].Parallel_R_RightEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Parallel_R_LeftEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Parallel_L_RightEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Parallel_L_LeftEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Parallel_C_RightEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Parallel_C_LeftEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Series_R_RightEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Series_R_LeftEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Series_L_RightEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Series_L_LeftEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Series_C_RightEnd_devia = 0.0;
            sgg.Med(contamedia).wire[1].Series_C_LeftEnd_devia = 0.0;
//fin stoch
            //
            if (this.twires.TW[j].TL == MATERIAL_absorbing) {
                sgg.Med(contamedia).wire[1].HasAbsorbing_LeftEnd = true;
            } else if (this.twires.TW[j].TL == Parallel_CONS) {
                sgg.Med(contamedia).wire[1].HasParallel_LeftEnd = true;
                sgg.Med(contamedia).wire[1].Parallel_R_LeftEnd = this.twires.TW[j].R_LeftEnd;
                sgg.Med(contamedia).wire[1].Parallel_L_LeftEnd = this.twires.TW[j].L_LeftEnd;
                sgg.Med(contamedia).wire[1].Parallel_C_LeftEnd = this.twires.TW[j].C_LeftEnd;

                sgg.Med(contamedia).wire[1].Parallel_R_LeftEnd_devia = this.twires.TW[j].R_LeftEnd_devia;
                sgg.Med(contamedia).wire[1].Parallel_L_LeftEnd_devia = this.twires.TW[j].L_LeftEnd_devia;
                sgg.Med(contamedia).wire[1].Parallel_C_LeftEnd_devia = this.twires.TW[j].C_LeftEnd_devia;

            } else if (this.twires.TW[j].TL == SERIES_CONS) {
                sgg.Med(contamedia).wire[1].HasSeries_LeftEnd = true;
                sgg.Med(contamedia).wire[1].Series_R_LeftEnd = this.twires.TW[j].R_LeftEnd;
                sgg.Med(contamedia).wire[1].Series_L_LeftEnd = this.twires.TW[j].L_LeftEnd;
                sgg.Med(contamedia).wire[1].Series_C_LeftEnd = this.twires.TW[j].C_LeftEnd;

                sgg.Med(contamedia).wire[1].Series_R_LeftEnd_devia = this.twires.TW[j].R_LeftEnd_devia;
                sgg.Med(contamedia).wire[1].Series_L_LeftEnd_devia = this.twires.TW[j].L_LeftEnd_devia;
                sgg.Med(contamedia).wire[1].Series_C_LeftEnd_devia = this.twires.TW[j].C_LeftEnd_devia;
            } else if (this.twires.TW[j].TL == DISPERSIVE_CONS) {
                sgg.Med(contamedia).wire[1].disp_LeftEnd.resize(1);
                asignawiredisper(sgg.Med(contamedia).wire[1].disp_LeftEnd[0],
                    this.twires.TW[j].dispfile_LeftEnd);
            }
            //

            if (this.twires.TW[j].TR == MATERIAL_absorbing) {
                sgg.Med(contamedia).wire[1].HasAbsorbing_RightEnd = true;
            } else if (this.twires.TW[j].TR == Parallel_CONS) {
                sgg.Med(contamedia).wire[1].HasParallel_RightEnd = true;
                sgg.Med(contamedia).wire[1].Parallel_R_RightEnd = this.twires.TW[j].R_RightEnd;
                sgg.Med(contamedia).wire[1].Parallel_L_RightEnd = this.twires.TW[j].L_RightEnd;
                sgg.Med(contamedia).wire[1].Parallel_C_RightEnd = this.twires.TW[j].C_RightEnd;

                sgg.Med(contamedia).wire[1].Parallel_R_RightEnd_devia = this.twires.TW[j].R_RightEnd_devia;
                sgg.Med(contamedia).wire[1].Parallel_L_RightEnd_devia = this.twires.TW[j].L_RightEnd_devia;
                sgg.Med(contamedia).wire[1].Parallel_C_RightEnd_devia = this.twires.TW[j].C_RightEnd_devia;
            } else if (this.twires.TW[j].TR == SERIES_CONS) {
                sgg.Med(contamedia).wire[1].HasSeries_RightEnd = true;
                sgg.Med(contamedia).wire[1].Series_R_RightEnd = this.twires.TW[j].R_RightEnd;
                sgg.Med(contamedia).wire[1].Series_L_RightEnd = this.twires.TW[j].L_RightEnd;
                sgg.Med(contamedia).wire[1].Series_C_RightEnd = this.twires.TW[j].C_RightEnd;

                sgg.Med(contamedia).wire[1].Series_R_RightEnd_devia = this.twires.TW[j].R_RightEnd_devia;
                sgg.Med(contamedia).wire[1].Series_L_RightEnd_devia = this.twires.TW[j].L_RightEnd_devia;
                sgg.Med(contamedia).wire[1].Series_C_RightEnd_devia = this.twires.TW[j].C_RightEnd_devia;

            } else if (this.twires.TW[j].TR == DISPERSIVE_CONS) {
                sgg.Med(contamedia).wire[1].disp_RightEnd.resize(1);
                asignawiredisper(sgg.Med(contamedia).wire[1].disp_RightEnd[0],
                    this.twires.TW[j].dispfile_RightEnd);
            }

//stoch
            rdummy = std::abs(sgg.Med(contamedia).wire[1].radius_devia)
                + std::abs(sgg.Med(contamedia).wire[1].l_devia)
                + std::abs(sgg.Med(contamedia).wire[1].C_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Parallel_L_RightEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Parallel_L_LeftEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Parallel_C_RightEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Parallel_C_LeftEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Series_L_RightEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Series_L_LeftEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Series_C_RightEnd_devia)
                + std::abs(sgg.Med(contamedia).wire[1].Series_C_LeftEnd_devia);
            if (rdummy > 1.0e-15) {
                sprintf(buff, "%s", "Non null deviations found in L, C or radius in wires stoch. Still unsupported.");
                STOPONERROR(layoutnumber, num_procs, buff);
            }
//fin stoch
            //
            //esto se soportaba desde versiones antiguas (hilos de un solo segmento. Por error se descomento en la R2417 cuando se trabajo en lo del strictnfde tras vuelta de madrid
            //vuelvo a comentarlo porque si que tenemos la capacidad de hilos de un solo segmento
            //
            //        tama2 = this.twires.TW[j].N_TWC;
            //        if (tama2 == 1) {
            //           stoponerror(layoutnumber, num_procs, "A WIRE must have at least two segments");
            //        }
            //
            //esto no es ya necesario porque lo calculo yo luego en el wires
            //!record the LeftEnd and RightEnd coordinates (first and last points)
            //!
            //        sgg.Med(contamedia).wire[1].LextremoI = this.twires.TW[j].TWC[1].i;
            //        sgg.Med(contamedia).wire[1].LextremoJ = this.twires.TW[j].TWC[1].j;
            //        sgg.Med(contamedia).wire[1].LextremoK = this.twires.TW[j].TWC[1].k;
            //        orientacionL = this.twires.TW[j].TWC[1].D;
            //        sgg.Med(contamedia).wire[1].RextremoI = this.twires.TW[j].TWC[tama2].i;
            //        sgg.Med(contamedia).wire[1].RextremoJ = this.twires.TW[j].TWC[tama2].j;
            //        sgg.Med(contamedia).wire[1].RextremoK = this.twires.TW[j].TWC[tama2].k;
            //        orientacionR = this.twires.TW[j].TWC[tama2].D;
            //        //
            //!correct each ending
            //        numminus = 0;
            //        for (i = 2; i < tama2 - 1; i++) { //bug OLD 12/09/13 Model_unidos.nfde segmentos finales duplicados internamente
            //            punto.XI = this.twires.TW[j].TWC[i].i;
            //            punto.YI = this.twires.TW[j].TWC[i].j;
            //            punto.ZI = this.twires.TW[j].TWC[i].k;
            //            orientacion = this.twires.TW[j].TWC[i].D;
            //            switch (orientacion) {
            //            case iEx:
            //                if (punto.YI == sgg.Med(contamedia).wire[1].LextremoJ && ...) {

// !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%LextremoK)) then
            // !                    if ((orientacion /= orientacionL).and.(punto%XI   == sgg%Med(contamedia)%wire(1)%LextremoI)) numminus=numminus +1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
            // !                    if                                    (punto%XI+1 == sgg%Med(contamedia)%wire(1)%LextremoI) numminus =numminus  +1
            // !                end if
            // !            case (iEy)
            // !                if (                                      (punto%XI   == sgg%Med(contamedia)%wire(1)%LextremoI).and.  &
            // !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%LextremoK)) then
            // !                    if ((orientacion /= orientacionL).and.(punto%YI   == sgg%Med(contamedia)%wire(1)%LextremoJ)) numminus=numminus +1
            // !                    if                                    (punto%YI+1 == sgg%Med(contamedia)%wire(1)%LextremoJ) numminus =numminus  +1
            // !                end if
            // !            case (iEz)
            // !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%LextremoJ).and.  &
            // !                                                          (punto%XI   == sgg%Med(contamedia)%wire(1)%LextremoI)) then
            // !                    if ((orientacion /= orientacionL).and.(punto%ZI   == sgg%Med(contamedia)%wire(1)%LextremoK)) numminus=numminus +1
            // !                    if                                    (punto%ZI+1 == sgg%Med(contamedia)%wire(1)%LextremoK) numminus =numminus  +1
            // !                end if
            // !            end select
            // !
            // !        end do
            // !        if (numminus >= 1) then !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
            // !              select case (this%twires%TW(j)%TWC(1)%D)
            // !              case (iEx) !si son iguales a 2 es cerrado
            // !                  sgg%Med(contamedia)%wire(1)%LextremoI = sgg%Med(contamedia)%wire(1)%LextremoI + 1
            // !              case (iEy)
            // !                  sgg%Med(contamedia)%wire(1)%LextremoJ = sgg%Med(contamedia)%wire(1)%LextremoJ + 1
            // !              case (iEz)
            // !                  sgg%Med(contamedia)%wire(1)%LextremoK = sgg%Med(contamedia)%wire(1)%LextremoK + 1
            // !              end select
            // !        end if
            // !        !
            // !!correct each ending
            // !        numminus=0
            // !        do i = 2, tama2-1 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
            // !            punto%XI = this%twires%TW(j)%TWC(i)%i
            // !            punto%YI = this%twires%TW(j)%TWC(i)%j
            // !            punto%ZI = this%twires%TW(j)%TWC(i)%k
            // !            orientacion = this%twires%TW(j)%TWC(i)%D
            // !            select case (orientacion)
            // !            case (iEx)
            // !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%RextremoJ).and.  &
            // !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%RextremoK)) then
            // !                    if ((orientacion /= orientacionR).and.(punto%XI   == sgg%Med(contamedia)%wire(1)%RextremoI)) numminus=numminus +1
            // !                    if                                    (punto%XI+1 == sgg%Med(contamedia)%wire(1)%RextremoI) numminus =numminus  +1
            // !                end if
            // !            case (iEy)
            // !                if (                                      (punto%XI   == sgg%Med(contamedia)%wire(1)%RextremoI).and.  &
            // !                                                          (punto%ZI   == sgg%Med(contamedia)%wire(1)%RextremoK)) then
            // !                    if ((orientacion /= orientacionR).and.(punto%YI   == sgg%Med(contamedia)%wire(1)%RextremoJ)) numminus=numminus +1
            // !                    if                                    (punto%YI+1 == sgg%Med(contamedia)%wire(1)%RextremoJ) numminus =numminus  +1
            // !                end if
            // !            case (iEz)
            // !                if (                                      (punto%YI   == sgg%Med(contamedia)%wire(1)%RextremoJ).and.  &
            // !                                                          (punto%XI   == sgg%Med(contamedia)%wire(1)%RextremoI)) then
            // !                    if ((orientacion /= orientacionR).and.(punto%ZI   == sgg%Med(contamedia)%wire(1)%RextremoK)) numminus=numminus +1
            // !                    if                                    (punto%ZI+1 == sgg%Med(contamedia)%wire(1)%RextremoK) numminus =numminus  +1
            // !                end if
            // !            end select
            // !        end do
            // !        if ((numminus >= 1).or.(tama2 == 1)) then  !bug ca295 !bug OLD 12/09/13  Model_unidos.nfde segmentos finales duplicados internamente
            // !              select case (this%twires%TW(j)%TWC(tama2)%D)
            // !              case (iEx)
            // !                  sgg%Med(contamedia)%wire(1)%RextremoI = sgg%Med(contamedia)%wire(1)%RextremoI  + 1
            // !!si son iguales a 2 es cerrado
            // !              case (iEy)
            // !                  sgg%Med(contamedia)%wire(1)%RextremoJ = sgg%Med(contamedia)%wire(1)%RextremoJ  + 1
            // !              case (iEz)
            // !                  sgg%Med(contamedia)%wire(1)%RextremoK = sgg%Med(contamedia)%wire(1)%RextremoK  + 1
            // !              end select
            // !        end if
         } // end do !del tama

         // preanalisis de hilos embeddeds en materiales  antes de asignarlos
         tama = this->twires->n_tw;
         paraerrhilo = false;
         for (int j1 = 1; j1 <= tama; ++j1) {
            tama2 = this->twires->TW[j1]->N_TWC;
            for (int i1 = 1; i1 <= tama2; ++i1) {
               i = this->twires->TW[j1]->TWC[i1]->i;
               j = this->twires->TW[j1]->TWC[i1]->j;
               k = this->twires->TW[j1]->TWC[i1]->k;
               orientacion = this->twires->TW[j1]->TWC[i1]->D;
               OrigIndex = this->twires->TW[j1]->TWC[i1]->nd;
               if ((i >= BoundingBox->XI) && (i < BoundingBox->XE) &&
                   (j >= BoundingBox->YI) && (j < BoundingBox->YE) &&
                   (k >= BoundingBox->ZI) && (k < BoundingBox->ZE)) {
                  if (i > BoundingBox->XI) {
                     imenos1 = i - 1;
                  } else {
                     imenos1 = i;
                  }
                  if (j > BoundingBox->YI) {
                     jmenos1 = j - 1;
                  } else {
                     jmenos1 = j;
                  }
                  if (k > BoundingBox->ZI) {
                     kmenos1 = k - 1;
                  } else {
                     kmenos1 = k;
                  }
                  switch (orientacion) {
                  case iEx:
                     if ((media->sggMiEx[i][j][k] == 0) || (sgg->med(media->sggMiEx[i][j][k])->is->pec)) {
                        paraerrhilo = true;
                        sprintf(buff, "pre1_WARNING:   x-WIRE at %7d %5d %5d %5d embedded within PEC", OrigIndex, i, j, k);
                        if (verbose) WarnErrReport(buff);
                     } else if (media->sggMiEx[i][j][k] != 1) {
                        islossy = (sgg->Med(media->sggMiEx[i][j][k])->Sigma != 0.0_RKIND);
                        if (islossy) {
                           paraerrhilo = true;
                           sprintf(buff, "pre1_WARNING: x-WIRE at %7d %5d %5d %5d embedded within LOSSY medium %5d", OrigIndex, i, j, k, media->sggMiEx[i][j][k]);
                        } else {
                           sprintf(buff, "pre1_WARNING: x-WIRE at %7d %5d %5d %5d embedded within medium %5d", OrigIndex, i, j, k, media->sggMiEx[i][j][k]);
                        }
                        if (verbose) WarnErrReport(buff);
                     }
                     if ((((media->sggMiEy[i][j][k] == 0) || (sgg->med(media->sggMiEy[i][j][k])->is->pec)) ||

} else if (
                        ((media.sggMiEz(i, j, k) == 0) || (sgg.med(media.sggMiEz(i, j, k)).is.pec)) ||
                        ((media.sggMiEy(i, jmenos1, k) == 0) || (sgg.med(media.sggMiEy(i, jmenos1, k)).is.pec)) ||
                        ((media.sggMiEz(i, j, kmenos1) == 0) || (sgg.med(media.sggMiEz(i, j, kmenos1)).is.pec))
                    ) &&
                    ((media.sggMiEx(i, j, k) != 0) && (!sgg.med(media.sggMiEx(i, j, k)).is.pec))
                    ) {
                        if ((i1 != 1) && (i1 != tama2)) { // solo en LeftEnd y RightEnd pueden tocar
                            paraerrhilo = true;
                            std::ostringstream buff_stream;
                            buff_stream << "pre1_WARNING:   intermediate node of x-WIRE at  " << OrigIndex << " " << i << " " << j << " " << k << " touching PEC";
                            std::string buff = buff_stream.str();
                            if (verbose) WarnErrReport(buff);
                        } else {
                            // continue
                            // write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC'
                            // if (verbose) call WarnErrReport (buff)
                        }
                    } else if (
                        (
                            ((media.sggMiEy(i+1, j, k) == 0) || (sgg.med(media.sggMiEy(i+1, j, k)).is.pec)) ||
                            ((media.sggMiEz(i+1, j, k) == 0) || (sgg.med(media.sggMiEz(i+1, j, k)).is.pec)) ||
                            ((media.sggMiEy(i+1, jmenos1, k) == 0) || (sgg.med(media.sggMiEy(i+1, jmenos1, k)).is.pec)) ||
                            ((media.sggMiEz(i+1, j, kmenos1) == 0) || (sgg.med(media.sggMiEz(i+1, j, kmenos1)).is.pec))
                        ) &&
                        ((media.sggMiEx(i, j, k) != 0) && (!sgg.med(media.sggMiEx(i, j, k)).is.pec))
                    ) {
                        if ((i1 != 1) && (i1 != tama2)) { // solo en LeftEnd y RightEnd pueden tocar
                            paraerrhilo = true;
                            std::ostringstream buff_stream;
                            buff_stream << "pre1_WARNING:   intermediate node of x-WIRE at  " << OrigIndex << " " << i+1 << " " << j << " " << k << " touching PEC";
                            std::string buff = buff_stream.str();
                            if (verbose) WarnErrReport(buff);
                        } else {
                            // continue
                            // write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i+1, j, k,' touching PEC'
                            // if (verbose) call WarnErrReport (buff)
                        }
                    } else if ((media.sggMiEy(i, j, k) != 1) && (media.sggMiEx(i, j, k) == 1)) {
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING: x-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " touching medium " << media.sggMiEy(i, j, k);
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    } else if ((media.sggMiEz(i, j, k) != 1) && (media.sggMiEx(i, j, k) == 1)) {
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING: x-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " touching medium " << media.sggMiEz(i, j, k);
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    } else if ((media.sggMiEy(i+1, j, k) != 1) && (media.sggMiEx(i, j, k) == 1)) {
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING: x-WIRE at " << OrigIndex << " " << i+1 << " " << j << " " << k << " touching medium " << media.sggMiEy(i+1, j, k);
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    } else if ((media.sggMiEz(i+1, j, k) != 1) && (media.sggMiEx(i, j, k) == 1)) {
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING: x-WIRE at " << OrigIndex << " " << i+1 << " " << j << " " << k << " touching medium " << media.sggMiEz(i+1, j, k);
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    }
                    break;

                case iEy:
                    if ((media.sggMiEy(i, j, k) == 0) || (sgg.med(media.sggMiEy(i, j, k)).is.pec)) {
                        paraerrhilo = true;
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING:   y-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " embedded within PEC";
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    } else if (media.sggMiEy(i, j, k) != 1) {
                        islossy = (sgg.Med(media.sggMiEy(i, j, k)).Sigma != 0.0_RKIND);
                        if (islossy) {
                            paraerrhilo = true;
                            std::ostringstream buff_stream;
                            buff_stream << "pre1_WARNING: Y-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " embedded within LOSSY medium " << media.sggMiEY(i, j, k);
                            std::string buff = buff_stream.str();
                            if (verbose) WarnErrReport(buff);
                        } else {
                            std::ostringstream buff_stream;
                            buff_stream << "pre1_WARNING: y-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " embedded within medium " << media.sggMiEy(i, j, k);
                            std::string buff = buff_stream.str();
                            if (verbose) WarnErrReport(buff);
                        }
                    }
                    if (
                        (
                            ((media.sggMiEx(i, j, k) == 0) || (sgg.med(media.sggMiEx(i, j, k)).is.pec)) ||
                            ((media.sggMiEz(i, j, k) == 0) || (sgg.med(media.sggMiEz(i, j, k)).is.pec)) ||
                            ((media.sggMiEx(imenos1, j, k) == 0) || (sgg.med(media.sggMiEx(imenos1, j, k)).is.pec)) ||
                            ((media.sggMiEz(i, j, kmenos1) == 0) || (sgg.med(media.sggMiEz(i, j, kmenos1)).is.pec))
                        ) &&
                        ((media.sggMiEy(i, j, k) != 0) && (!sgg.med(media.sggMiEy(i, j, k)).is.pec))
                    ) {
                        if ((i1 != 1) && (i1 != tama2)) { // solo en LeftEnd y RightEnd pueden tocar
                            paraerrhilo = true;
                            std::ostringstream buff_stream;
                            buff_stream << "pre1_WARNING:   intermediate node of y-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " touching PEC";
                            std::string buff = buff_stream.str();
                            if (verbose) WarnErrReport(buff);
                        } else {
                            // continue
                            // write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j, k,' touching PEC'
                            // if (verbose) call WarnErrReport (buff)
                        }
                    } else if (
                        (
                            ((media.sggMiEx(i, j+1, k) == 0) || (sgg.med(media.sggMiEx(i, j+1, k)).is.pec)) ||
                            ((media.sggMiEz(i, j+1, k) == 0) || (sgg.med(media.sggMiEz(i, j+1, k)).is.pec)) ||
                            ((media.sggMiEx(imenos1, j+1, k) == 0) || (sgg.med(media.sggMiEx(imenos1, j+1, k)).is.pec)) ||
                            ((media.sggMiEz(i, j+1, kmenos1) == 0) || (sgg.med(media.sggMiEz(i, j+1, kmenos1)).is.pec))
                        ) &&
                        ((media.sggMiEy(i, j, k) != 0) && (!sgg.med(media.sggMiEy(i, j, k)).is.pec))
                    ) {
                        if ((i1 != 1) && (i1 != tama2)) { // solo en LeftEnd y RightEnd pueden tocar
                            paraerrhilo = true;
                            std::ostringstream buff_stream;
                            buff_stream << "pre1_WARNING:   intermediate node of y-WIRE at " << OrigIndex << " " << i << " " << j+1 << " " << k << " touching PEC";
                            std::string buff = buff_stream.str();
                            if (verbose) WarnErrReport(buff);
                        } else {
                            // continue
                            // write(buff, '(a,i7,3i5,a)')    'A node of terminal x-WIRE at ',OrigIndex, i, j+1, k,' touching PEC'
                            // if (verbose) call WarnErrReport (buff)
                        }
                    } else if ((media.sggMiEx(i, j, k) != 1) && (media.sggMiEy(i, j, k) == 1)) {
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING: y-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " touching medium " << media.sggMiEx(i, j, k);
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    } else if ((media.sggMiEz(i, j, k) != 1) && (media.sggMiEy(i, j, k) == 1)) {
                        std::ostringstream buff_stream;
                        buff_stream << "pre1_WARNING: y-WIRE at " << OrigIndex << " " << i << " " << j << " " << k << " touching medium " << media.sggMiEz(i, j, k);
                        std::string buff = buff_stream.str();
                        if (verbose) WarnErrReport(buff);
                    }
                    break;

