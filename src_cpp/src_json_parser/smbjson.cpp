#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <algorithm>
#include <cstring>

// Forward declarations and includes for external modules/types
// These would correspond to the Fortran use statements
// #include "NFDETypes_m.h"
// #include "NFDETypes_extension_m.h"
// #include "smbjson_labels_m.h"
// #include "mesh_m.h"
// #include "parser_tools_m.h"
// #include "idchildtable_m.h"
// #include "Report_m.h"
// #include "json_module.h"
// #include "json_kinds.h"
// #include "conformal_types_m.h"

// Assuming these types and functions exist in the included headers
// struct NFDEGeneral_t;
// struct MatrizMedios_t;
// struct Desplazamiento_t;
// struct Frontera_t;
// struct FronteraPML_t;
// struct Materials_t;
// struct PECRegions_t;
// struct ConformalPECRegions_t;
// struct ConformalPECElements_t;
// struct conformal_region_t;
// struct cell_interval_t;
// struct interval_t;
// struct coords_t;
// struct triangle_t;

// External functions assumed to be available
// void WarnErrReport(const std::string& msg, bool fatal);
// void initializeProblemDescription(Parseador_t& res);
// void addCoordinates(Mesh_t& mesh);
// void addElements(Mesh_t& mesh);
// void appendRegion(std::vector<coords_t>& resCoords, int& resNCoords, int& resNCoordsMax, const std::vector<coords_t>& cs);
// void appendRegion(std::vector<ConformalPECElements_t>& regions, const conformal_region_t& region, const std::string& tagName);
// std::vector<interval_t> copyIntervals(const std::vector<cell_interval_t>& intervals);
// void fillDielectricsOfCellType(std::vector<...>& vols, int cellType);
// void fillDielectricsOfCellType(std::vector<...>& surfs, int cellType);

// Constants from labels module (placeholders)
// #define J_MATERIALS "materials"
// #define J_MESH "mesh"
// #define J_ELEMENTS "elements"
// #define J_COORDINATES "coordinates"
// #define J_ID "id"
// #define J_COORDINATE_POS "position"
// #define J_TYPE "type"
// #define J_ELEM_TYPE_NODE "node"
// #define J_ELEM_TYPE_POLYLINE "polyline"
// #define J_ELEM_TYPE_CELL "cell"
// #define J_CONF_VOLUME_TRIANGLES "triangles"
// #define J_CELL_INTERVALS "intervals"
// #define J_SUBTYPE "subtype"
// #define J_CONF_SUBTYPE_VOLUME "volume"
// #define J_CONF_SUBTYPE_SURFACE "surface"
// #define J_GENERAL "general"
// #define J_GEN_ADDITIONAL_ARGUMENTS "additionalArguments"
// #define J_GEN_TIME_STEP "timeStep"
// #define J_GEN_NUMBER_OF_STEPS "numberOfSteps"
// #define J_GEN_MTLN_PROBLEM "mtlnProblem"
// #define J_GRID "grid"
// #define J_GRID_NUMBER_OF_CELLS "numberOfCells"
// #define J_GRID_STEPS "steps"
// #define J_GRID_ORIGIN "origin"
// #define J_BOUNDARY "boundary"
// #define J_BND_ALL "all"
// #define J_BND_XL "xl"
// #define J_BND_XU "xu"
// #define J_BND_YL "yl"
// #define J_BND_YU "yu"
// #define J_BND_ZL "zl"
// #define J_BND_ZU "zu"
// #define J_BND_PML_LAYERS "layers"
// #define J_BND_PML_ORDER "order"
// #define J_BND_PML_REFLECTION "reflection"
// #define J_BND_TYPE_PEC "pec"
// #define J_BND_TYPE_PMC "pmc"
// #define J_BND_TYPE_PERIODIC "periodic"
// #define J_BND_TYPE_MUR "mur"
// #define J_BND_TYPE_PML "pml"
// #define J_BACKGROUND "background"
// #define J_BKG_ABS_PERMITTIVITY "permittivity"
// #define J_BKG_ABS_PERMEABILITY "permeability"
// #define J_MAT_TYPE_PEC "pec"
// #define J_MAT_TYPE_PMC "pmc"
// #define J_CONF_SUBTYPE_VOLUME "volume"
// #define J_CONF_SUBTYPE_SURFACE "surface"

// Enumerations/Constants for boundary types and cell types
// enum { F_PML, F_PEC, F_PMC, F_PER, F_MUR, F_XL, F_XU, F_YL, F_YU, F_ZL, F_ZU };
// enum { CELL_TYPE_LINEL, CELL_TYPE_SURFEL, CELL_TYPE_VOXEL };
// enum { REGION_TYPE_VOLUME, REGION_TYPE_SURFACE };
// enum { NP_T1_PLAIN, NP_T2_TIME };

// RKIND and CK are kind parameters
using RKIND = double;
using CK = char;

namespace smbjson_m {

    // Placeholder for BUFSIZE if not defined elsewhere
    #ifndef BUFSIZE
    const int BUFSIZE = 1024;
    #endif

    // Placeholder for NP_T1_PLAIN, NP_T2_TIME if not defined elsewhere
    #ifndef NP_T1_PLAIN
    const int NP_T1_PLAIN = 0;
    #endif
    #ifndef NP_T2_TIME
    const int NP_T2_TIME = 1;
    #endif

    struct parser_t {
    private:
        std::string filename;
        std::shared_ptr<json_file> jsonfile; // Assuming json_file is a class with pointer semantics
        std::shared_ptr<json_core> core;
        std::shared_ptr<json_value> root;
        mesh_t mesh;
        IdChildTable_t matTable, elementTable;

    public:
        bool isInitialized = false;

        parser_t() : jsonfile(std::make_shared<json_file>()), core(std::make_shared<json_core>()), root(std::make_shared<json_value>()) {}

        void readProblemDescription(Parseador_t& res);
        Mesh_t readMesh();

#ifdef CompileWithMTLN
        void readMTLN(MTLN_t& res); // Placeholder return type
#endif

    private:
        std::string readAdditionalArguments();
        NFDEGeneral_t readGeneral();
        Desplazamiento_t readGrid();
        MatrizMedios_t readMediaMatrix();
        void readBackgroundMaterial(Materials_t& mats);
        PECRegions_t readPECRegions();
        PECRegions_t readPMCRegions();
        DielectricRegions_t readDielectricRegions();
        void readLossyThinSurfaces(LossyThinSurfaces_t& res); // Placeholder
        Frontera_t readBoundary();
        std::vector<Planewave_t> readPlanewaves(); // Placeholder
        std::vector<NodalSource_t> readNodalSources(); // Placeholder
        std::vector<Probe_t> readProbes(); // Placeholder
        std::vector<Probe_t> readMoreProbes(); // Placeholder
        std::vector<BlockProbe_t> readBlockProbes(); // Placeholder
        std::vector<VolumicProbe_t> readVolumicProbes(); // Placeholder
        std::vector<ThinWire_t> readThinWires(); // Placeholder
        std::vector<ThinSlot_t> readThinSlots(); // Placeholder
        ConformalPECRegions_t readConformalRegions();

        // Helper methods
        bool getLogicalAt(const std::shared_ptr<json_value>& val, const std::string& key, bool default_val = false);
        int getIntAt(const std::shared_ptr<json_value>& val, const std::string& key);
        std::vector<int> getIntsAt(const std::shared_ptr<json_value>& val, const std::string& key);
        double getRealAt(const std::shared_ptr<json_value>& val, const std::string& key, double default_val = 0.0);
        std::vector<double> getRealsAt(const std::shared_ptr<json_value>& val, const std::string& key);
        // Matrix handling omitted for brevity, assuming similar pattern
        std::string getStrAt(const std::shared_ptr<json_value>& val, const std::string& key, const std::string& default_val = "");
        bool existsAt(const std::shared_ptr<json_value>& val, const std::string& key);
        int dimensionAt(const std::shared_ptr<json_value>& val, const std::string& key);
        std::shared_ptr<json_value> getDomain(const std::string& path);
        void buildPECPMCRegions(); // Logic moved to buildPECPMCRegions function below
        std::vector<materialAssociation_t> getMaterialAssociations(const std::vector<std::string>& matTypes, const std::vector<std::string>& elemTypes);
        void parseMaterialAssociation(const materialAssociation_t& ma);
        void matAssToCoords(std::vector<coords_t>& cs, const materialAssociation_t& ma, int cellType);
        std::string buildTagName(int materialId, int elementId);
        std::vector<std::shared_ptr<json_value>> jsonValueFilterByKeyValue(const std::vector<std::shared_ptr<json_value>>& values, const std::string& key, const std::string& value);
        std::vector<std::shared_ptr<json_value>> jsonValueFilterByKeyValues(const std::vector<std::shared_ptr<json_value>>& values, const std::string& key, const std::vector<std::string>& values_list);
        std::vector<int> getSingleVolumeInElementsIds(int elementId);
    };

    struct thinwiretermination_t {
        int terminationType;
        double r, l, c;
    };

    struct generator_description_t {
        std::string srctype, srcfile;
        double multiplier;
    };

    struct materialAssociation_t {
        std::string name;
        int materialId;
        std::vector<int> elementIds;
        std::string matAssType;
        // Cable specific fields.
        int initialTerminalId = -1;
        int endTerminalId = -1;
        int initialConnectorId = -1;
        int endConnectorId = -1;
        int containedWithinElementId = -1;
        // Optional total resistance override.
        std::vector<double> totalResistance;
        bool hasTotalResistance = false;
    };

    struct domain_t {
        double tstart = 0.0, tstop = 0.0, tstep = 0.0;
        double fstart = 0.0, fstop = 0.0;
        double fstep = 0.0;
        std::string filename;
        int type1 = NP_T1_PLAIN, type2 = NP_T2_TIME;
        bool isLogarithmicFrequencySpacing = false;
    };

    // Implementation of parser_t methods

    void parser_t::readProblemDescription(Parseador_t& res) {
        mesh = readMesh();
        // Assuming IdChildTable_t constructor takes core, root, and key
        matTable = IdChildTable_t(core, root, J_MATERIALS);
        elementTable = IdChildTable_t(core, root, J_MESH + "." + J_ELEMENTS);
        
        initializeProblemDescription(res);
        
        res.switches = readAdditionalArguments();
        
        // Basics
        res.general = readGeneral();
        res.matriz = readMediaMatrix();
        res.despl = readGrid();
        res.front = readBoundary();
        
        // Materials
        readBackgroundMaterial(res.mats);
        res.pecRegs = readPECRegions();
        res.pmcRegs = readPMCRegions();
        res.dielRegs = readDielectricRegions();
        readLossyThinSurfaces(res.lossyThinSurfs);
        
        // Sources
        res.plnSrc = readPlanewaves();
        res.nodSrc = readNodalSources();
        
        // Probes
        res.oldSonda = readProbes();
        res.sonda = readMoreProbes();
        res.BloquePrb = readBlockProbes();
        res.VolPrb = readVolumicProbes();
        
        // Conformal elements
        res.conformalRegs = readConformalRegions();

        // Thin elements
#ifdef CompileWithMTLN 
        readMTLN(res.mtln);
#else
        // Placeholder for readThinWires logic which modifies res.tWires and res.sonda
        // This part is complex and depends on internal state, simplified here
        std::vector<ThinWire_t> tWires = readThinWires();
        // res.tWires = tWires; 
        // res.sonda update logic omitted
#endif
        res.tSlots = readThinSlots();
    }

    Mesh_t parser_t::readMesh() {
        Mesh_t res;
        addCoordinates(res);
        addElements(res);
        return res;
    }

    void parser_t::addCoordinates(Mesh_t& mesh) {
        bool found = false;
        auto jcs = getDomain(J_MESH + "." + J_COORDINATES);
        if (jcs && existsAt(jcs, "")) { // Simplified check
             // Assuming core->count works on the value
             int numberOfCoordinates = core->count(jcs);
             mesh.allocateCoordinates(50 * numberOfCoordinates);
             for (int i = 1; i <= numberOfCoordinates; ++i) {
                 auto jc = core->get_child(jcs, i);
                 int id = getIntAt(jc, J_ID);
                 auto pos = getRealsAt(jc, J_COORDINATE_POS);
                 coordinate_t c;
                 c.position = pos;
                 mesh.addCoordinate(id, c);
             }
        }
    }

    void parser_t::addElements(Mesh_t& mesh) {
        bool found = false;
        auto jes = getDomain(J_MESH + "." + J_ELEMENTS);
        int numberOfElements = core->count(jes);
        mesh.allocateElements(50 * numberOfElements);
        
        if (jes) {
            for (int i = 1; i <= numberOfElements; ++i) {
                auto je = core->get_child(jes, i);
                int id = getIntAt(je, J_ID);
                std::string elementType = getStrAt(je, J_TYPE);
                
                if (elementType == J_ELEM_TYPE_NODE) {
                    auto coordIds = getIntsAt(je, J_COORDINATE_IDS);
                    node_t node;
                    node.coordIds = coordIds;
                    mesh.addElement(id, node);
                } else if (elementType == J_ELEM_TYPE_POLYLINE) {
                    auto coordIds = getIntsAt(je, J_COORDINATE_IDS);
                    polyline_t polyline;
                    polyline.coordIds = coordIds;
                    mesh.addElement(id, polyline);
                } else if (elementType == J_ELEM_TYPE_CELL) {
                    bool isConformal = false;
                    auto triangles = core->get(je, J_CONF_VOLUME_TRIANGLES, isConformal);
                    if (!isConformal) {
                        cell_region_t cR;
                        cR.intervals = readCellIntervals(je, J_CELL_INTERVALS);
                        mesh.addCellRegion(id, cR);
                    } else {
                        conformal_region_t cV;
                        cV.triangles = readTriangles(je, J_CONF_VOLUME_TRIANGLES);
                        for (size_t k = 0; k < cV.triangles.size(); ++k) {
                            for (int j = 0; j < 3; ++j) {
                                auto c = mesh.getCoordinate(cV.triangles[k].vertices[j].id);
                                cV.triangles[k].vertices[j].position[0] = c.position[0];
                                cV.triangles[k].vertices[j].position[1] = c.position[1];
                                cV.triangles[k].vertices[j].position[2] = c.position[2];
                            }
                        }
                        cV.intervals = readCellIntervals(je, J_CELL_INTERVALS);
                        std::string subtype = getStrAt(je, J_SUBTYPE);

                        if (subtype == J_CONF_SUBTYPE_VOLUME) cV.type = REGION_TYPE_VOLUME;
                        if (subtype == J_CONF_SUBTYPE_SURFACE) cV.type = REGION_TYPE_SURFACE;

                        mesh.addConformalRegion(id, cV);
                    }
                } else {
                    WarnErrReport("Invalid element type", true);
                }
            }
        }
    }

    std::vector<cell_interval_t> parser_t::readCellIntervals(const std::shared_ptr<json_value>& place, const std::string& path) {
        bool containsInterval = false;
        auto intervalsPlace = core->get(place, path, containsInterval);
        if (!containsInterval) {
            return {};
        }
        int nIntervals = core->count(intervalsPlace);
        std::vector<cell_interval_t> res(nIntervals);
        for (int i = 0; i < nIntervals; ++i) {
            auto interval = core->get_child(intervalsPlace, i + 1); // 1-based indexing in Fortran
            auto cellIni = getRealsAt(interval, "(1)");
            auto cellEnd = getRealsAt(interval, "(2)");
            res[i].ini.cell[0] = cellIni[0];
            res[i].ini.cell[1] = cellIni[1];
            res[i].ini.cell[2] = cellIni[2];
            res[i].end.cell[0] = cellEnd[0];
            res[i].end.cell[1] = cellEnd[1];
            res[i].end.cell[2] = cellEnd[2];
        }
        return res;
    }

    std::vector<triangle_t> parser_t::readTriangles(const std::shared_ptr<json_value>& place, const std::string& path) {
        bool containsTriangles = false;
        auto triangles = core->get(place, path, containsTriangles);
        if (!containsTriangles) {
            return {};
        }
        int nTriangles = core->count(triangles);
        std::vector<triangle_t> res(nTriangles);
        for (int i = 0; i < nTriangles; ++i) {
            auto triangle_ptr = core->get_child(triangles, i + 1);
            auto triangle = getRealsAt(triangle_ptr, ""); // Assuming direct array access
            for (int j = 0; j < 3; ++j) {
                res[i].vertices[j].id = static_cast<int>(triangle[j]);
            }
        }
        return res;
    }

    std::string parser_t::readAdditionalArguments() {
        return getStrAt(root, J_GENERAL + "." + J_GEN_ADDITIONAL_ARGUMENTS, "");
    }

    NFDEGeneral_t parser_t::readGeneral() {
        NFDEGeneral_t res;
        res.dt = getRealAt(root, J_GENERAL + "." + J_GEN_TIME_STEP, 0.0);
        if (res.dt < 0) WarnErrReport("timStep cannot be negative", true);
        res.nmax = getRealAt(root, J_GENERAL + "." + J_GEN_NUMBER_OF_STEPS);
        if (res.nmax <= 0) WarnErrReport("numberOfSteps has to be positive", true);
        res.mtlnProblem = getLogicalAt(root, J_GENERAL + "." + J_GEN_MTLN_PROBLEM, false);
        return res;
    }

    MatrizMedios_t parser_t::readMediaMatrix() {
        MatrizMedios_t res;
        std::string P = J_MESH + "." + J_GRID + "." + J_GRID_NUMBER_OF_CELLS;
        res.totalX = getIntAt(root, P + "(1)") + 1;
        res.totalY = getIntAt(root, P + "(2)") + 1;
        res.totalZ = getIntAt(root, P + "(3)") + 1;
        return res;
    }

    Desplazamiento_t parser_t::readGrid() {
        Desplazamiento_t res;
        std::string P = J_MESH + "." + J_GRID;
        int nX = getIntAt(root, P + "." + J_GRID_NUMBER_OF_CELLS + "(1)");
        int nY = getIntAt(root, P + "." + J_GRID_NUMBER_OF_CELLS + "(2)");
        int nZ = getIntAt(root, P + "." + J_GRID_NUMBER_OF_CELLS + "(3)");

        res.nX = nX;
        res.nY = nY;
        res.nZ = nZ;

        // Helper lambda for assignDes logic
        auto assignDes = [&](const std::string& path, std::vector<double>& dest, int& n) {
            bool found = false;
            auto vec = getRealsAt(root, path, found);
            if (!found) {
                WarnErrReport("Error reading grid: steps not found.", true);
            }
            if (vec.size() != 1 && vec.size() != static_cast<size_t>(n)) {
                WarnErrReport("Error reading grid: steps must be arrays of size 1 (for regular grids) or size equal to the number of cells.", true);
            }

            if (vec.size() == 1) {
                n = 1;
                dest.resize(1);
                dest[0] = vec[0];
            } else {
                dest.resize(n);
                dest = vec;
            }
        };

        assignDes(P + "." + J_GRID_STEPS + ".x", res.desX, res.nX);
        assignDes(P + "." + J_GRID_STEPS + ".y", res.desY, res.nY);
        assignDes(P + "." + J_GRID_STEPS + ".z", res.desZ, res.nZ);

        res.originx = getRealAt(root, P + "." + J_GRID_ORIGIN + "(1)", 0.0);
        res.originy = getRealAt(root, P + "." + J_GRID_ORIGIN + "(2)", 0.0);
        res.originz = getRealAt(root, P + "." + J_GRID_ORIGIN + "(3)", 0.0);

        res.mx1 = 0;
        res.my1 = 0;
        res.mz1 = 0;
        res.mx2 = nX;
        res.my2 = nY;
        res.mz2 = nZ;

        return res;
    }

    Frontera_t parser_t::readBoundary() {
        Frontera_t res;
        std::string bdrType;
        bool found = false;
        auto bdrs = core->get(root, J_BOUNDARY, found);
        if (!found) {
            WarnErrReport("Error reading boundary: " + std::string(J_BOUNDARY) + " not found.", true);
        }

        {
            bdrType = getStrAt(bdrs, J_BND_ALL + "." + J_TYPE, found);
            if (found) {
                // Assuming labelToBoundaryType returns an int and res.tipoFrontera is an array/vector
                // This part is highly dependent on specific struct definitions
                // Placeholder logic
                // res.tipoFrontera.fill(labelToBoundaryType(bdrType));
                // if (all(res.tipoFrontera == F_PML)) {
                //     res.propiedadesPML = readPMLProperties(J_BOUNDARY "." J_BND_ALL);
                // }
                return res; // Simplified return
            }
        }
        
        {
            std::vector<std::string> placeLabels = {J_BND_XL, J_BND_XU, J_BND_YL, J_BND_YU, J_BND_ZL, J_BND_ZU};
            for (int i = 0; i < 6; ++i) {
                bdrType = getStrAt(bdrs, placeLabels[i] + "." + J_TYPE, found);
                if (!found) {
                    WarnErrReport("Error reading boundary: " + placeLabels[i] + " or " + J_BND_ALL + " not found.", true);
                }
                int j = labelToBoundaryPlace(placeLabels[i]);
                // res.tipoFrontera[j] = labelToBoundaryType(bdrType);
                // if (res.tipoFrontera[j] == F_PML) {
                //     res.propiedadesPML[j] = readPMLProperties(J_BOUNDARY "." placeLabels[i]);
                // }
            }
        }
        return res;
    }

    void parser_t::readBackgroundMaterial(Materials_t& mats) {
        bool found = false;
        double val = getRealAt(root, J_BACKGROUND + "." + J_BKG_ABS_PERMITTIVITY, found);
        if (found) mats.mats[0].eps = val; // Assuming 1-based indexing in Fortran maps to 0-based in C++ or mats is 1-based

        found = false;
        val = getRealAt(root, J_BACKGROUND + "." + J_BKG_ABS_PERMEABILITY, found);
        if (found) mats.mats[0].mu = val;
    }

    PECRegions_t parser_t::readPECRegions() {
        return buildPECPMCRegions(J_MAT_TYPE_PEC);
    }

    PECRegions_t parser_t::readPMCRegions() {
        return buildPECPMCRegions(J_MAT_TYPE_PMC);
    }

    PECRegions_t parser_t::buildPECPMCRegions(const std::string& matType) {
        PECRegions_t res;
        std::vector<materialAssociation_t> mAs = getMaterialAssociations(
            {matType}, 
            {"-" + J_CONF_SUBTYPE_SURFACE, J_ELEM_TYPE_CELL + "    ", "-" + J_CONF_SUBTYPE_VOLUME + " "}
        );
        
        if (mAs.empty()) {
            std::vector<coords_t> emptyCoords;
            // appendRegion(res.lins, res.nLins, res.nLins_max, emptyCoords);
            // appendRegion(res.surfs, res.nSurfs, res.nSurfs_max, emptyCoords);
            // appendRegion(res.vols, res.nVols, res.nVols_max, emptyCoords);
            return res;
        }
        
        for (size_t i = 0; i < mAs.size(); ++i) {
            std::vector<coords_t> cs;
            matAssToCoords(cs, mAs[i], CELL_TYPE_LINEL);
            // appendRegion(res.lins, res.nLins, res.nLins_max, cs);
            cs.clear();
            matAssToCoords(cs, mAs[i], CELL_TYPE_SURFEL);
            // appendRegion(res.surfs, res.nSurfs, res.nSurfs_max, cs);
            cs.clear();
            matAssToCoords(cs, mAs[i], CELL_TYPE_VOXEL);
            // appendRegion(res.vols, res.nVols, res.nVols_max, cs);
        }
        return res;
    }

    ConformalPECRegions_t parser_t::readConformalRegions() {
        ConformalPECRegions_t res;
        std::vector<materialAssociation_t> mAs = getMaterialAssociations(
            {J_MAT_TYPE_PEC}, 
            {J_CONF_SUBTYPE_VOLUME, J_CONF_SUBTYPE_SURFACE}
        );

        for (size_t i = 0; i < mAs.size(); ++i) {
            for (size_t j = 0; j < mAs[i].elementIds.size(); ++j) {
                bool found = false;
                conformal_region_t cR = mesh.getConformalRegion(mAs[i].elementIds[j], found);
                if (found) {
                    std::string tagName = buildTagName(mAs[i].materialId, mAs[i].elementIds[j]);
                    if (cR.type == REGION_TYPE_VOLUME) {
                        // appendRegion(res.volumes, cR, tagName);
                    }
                    if (cR.type == REGION_TYPE_SURFACE) {
                        // appendRegion(res.surfaces, cR, tagName);
                    }
                }
            }
        }
        return res;
    }

    DielectricRegions_t parser_t::readDielectricRegions() {
        DielectricRegions_t res;
        fillDielectricsOfCellType(res.vols, CELL_TYPE_VOXEL);
        fillDielectricsOfCellType(res.surfs, CELL_TYPE_SURFEL);
        return res;
    }

    // Helper method implementations
    bool parser_t::getLogicalAt(const std::shared_ptr<json_value>& val, const std::string& key, bool default_val) {
        // Placeholder implementation
        return default_val;
    }

    int parser_t::getIntAt(const std::shared_ptr<json_value>& val, const std::string& key) {
        // Placeholder implementation
        return 0;
    }

    std::vector<int> parser_t::getIntsAt(const std::shared_ptr<json_value>& val, const std::string& key) {
        // Placeholder implementation
        return {};
    }

    double parser_t::getRealAt(const std::shared_ptr<json_value>& val, const std::string& key, double default_val) {
        // Placeholder implementation
        return default_val;
    }

    std::vector<double> parser_t::getRealsAt(const std::shared_ptr<json_value>& val, const std::string& key) {
        // Placeholder implementation
        return {};
    }

    std::string parser_t::getStrAt(const std::shared_ptr<json_value>& val, const std::string& key, const std::string& default_val) {
        // Placeholder implementation
        return default_val;
    }

    bool parser_t::existsAt(const std::shared_ptr<json_value>& val, const std::string& key) {
        // Placeholder implementation
        return false;
    }

    int parser_t::dimensionAt(const std::shared_ptr<json_value>& val, const std::string& key) {
        // Placeholder implementation
        return 0;
    }

    std::shared_ptr<json_value> parser_t::getDomain(const std::string& path) {
        // Placeholder implementation
        return root;
    }

    std::vector<materialAssociation_t> parser_t::getMaterialAssociations(const std::vector<std::string>& matTypes, const std::vector<std::string>& elemTypes) {
        // Placeholder implementation
        return {};
    }

    void parser_t::matAssToCoords(std::vector<coords_t>& cs, const materialAssociation_t& ma, int cellType) {
        // Placeholder implementation
    }

    std::string parser_t::buildTagName(int materialId, int elementId) {
        return "tag_" + std::to_string(materialId) + "_" + std::to_string(elementId);
    }

    // Helper functions for boundary
    int labelToBoundaryPlace(const std::string& str) {
        if (str == J_BND_XL) return F_XL;
        if (str == J_BND_XU) return F_XU;
        if (str == J_BND_YL) return F_YL;
        if (str == J_BND_YU) return F_YU;
        if (str == J_BND_ZL) return F_ZL;
        if (str == J_BND_ZU) return F_ZU;
        return 0;
    }

    int labelToBoundaryType(const std::string& str) {
        if (str == J_BND_TYPE_PEC) return F_PEC;
        if (str == J_BND_TYPE_PMC) return F_PMC;
        if (str == J_BND_TYPE_PERIODIC) return F_PER;
        if (str == J_BND_TYPE_MUR) return F_MUR;
        if (str == J_BND_TYPE_PML) return F_PML;
        return 0;
    }

} // namespace smbjson_m

fillDielectricsOfCellType(res.lins, CELL_TYPE_LINEL);
      
      res.nVols = static_cast<int>(res.vols.size());
      res.nSurfs = static_cast<int>(res.Surfs.size());
      res.nLins = static_cast<int>(res.Lins.size());

      res.nVols_max = res.nVols;
      res.nSurfs_max = res.nSurfs;
      res.nLins_max = res.nLins;
   }

   void fillDielectricsOfCellType(std::vector<dielectric_t>& res, int cellType) {
      auto mAs = this->getMaterialAssociations({J_MAT_TYPE_ISOTROPIC, J_MAT_TYPE_LUMPED});
      if (mAs.empty()) {
         res.clear();
         return;
      }

      // Precounts
      int nDielectrics = 0;
      for (const auto& mA : mAs) {           
         if (containsCellRegionsWithType(mA, cellType)) {
            nDielectrics++;
         } 
      }

      // Fills
      res.resize(nDielectrics);
      
      if (nDielectrics == 0) return;

      int j = 0;
      mAs = this->getMaterialAssociations({J_MAT_TYPE_ISOTROPIC});
      for (const auto& mA : mAs) {       
         if (!containsCellRegionsWithType(mA, cellType)) continue;
         j++;
         res[j-1] = readDielectric(mA, cellType);
      }

      mAs = this->getMaterialAssociations({J_MAT_TYPE_LUMPED});
      for (const auto& mA : mAs) {
         if (!containsCellRegionsWithType(mA, cellType)) continue;
         j++;
         res[j-1] = readLumped(mA, cellType);
      }
   }

   Dielectric_t readDielectric(const materialAssociation_t& mA, int cellType) {
      Dielectric_t res;
      cell_region_t cR;
      std::vector<coords_t> coords;
      json_value_ptr_t matPtr;
      int e, j;

      res.c1P.clear();
      res.n_c1p = 0;
      this->matAssToCoords(res.c2p, mA, cellType);
      res.n_c2p = static_cast<int>(res.c2p.size());
      
      matPtr = this->matTable.getId(mA.materialId);
      // Fills rest of dielectric data.
      res.sigma  = this->getRealAt(matPtr.p, J_MAT_ELECTRIC_CONDUCTIVITY, 0.0_RKIND);
      res.sigmam = this->getRealAt(matPtr.p, J_MAT_MAGNETIC_CONDUCTIVITY, 0.0_RKIND);
      res.eps    = this->getRealAt(matPtr.p, J_MAT_REL_PERMITTIVITY, 1.0_RKIND) * EPSILON_VACUUM;
      res.mu     = this->getRealAt(matPtr.p, J_MAT_REL_PERMEABILITY, 1.0_RKIND) * MU_VACUUM;

      return res;
   }

   Dielectric_t readLumped(const materialAssociation_t& mA, int cellType) {
      Dielectric_t res;
      cell_region_t cR;
      std::vector<coords_t> coords;
      json_value_ptr_t matPtr;
      int e, j;
      std::string model;
      bool found;
      const std::string errorMsgInit = "ERROR reading lumped material: ";
      char errorMsg[BUFSIZE];

      res.c1P.clear();
      res.n_c1p = 0;
      this->matAssToCoords(res.c2p, mA, cellType);
      res.n_c2p = static_cast<int>(res.c2p.size());
      
      matPtr = this->matTable.getId(mA.materialId);
      
      // Get the model type
      model = this->getStrAt(matPtr.p, J_MAT_LUMPED_MODEL, found);
      if (!found) {
         snprintf(errorMsg, BUFSIZE, "%s%d model not found.", errorMsgInit.c_str(), mA.materialId);
         WarnErrReport(errorMsg, true);
      }
      
      // Not really needed for resistor, inductor, or capacitor. 
      // But avoids error in lumped initialization.
      res.orient = 1;
      res.DiodOri = 1;

      res.eps = EPSILON_VACUUM;
      res.mu = MU_VACUUM;

      // Handle resistor model
      if (model == J_MAT_LUMPED_MODEL_RESISTOR) {
         res.resistor = true;
         res.R = this->getRealAt(matPtr.p, J_MAT_LUMPED_RESISTANCE, found);
         if (!found) {
            snprintf(errorMsg, BUFSIZE, "%s%d resistance not found.", errorMsgInit.c_str(), mA.materialId);
            WarnErrReport(errorMsg, true);
         }
         res.Rtime_on = this->getRealAt(matPtr.p, J_MAT_LUMPED_STARTING_TIME, 0.0_RKIND);
         res.Rtime_off = this->getRealAt(matPtr.p, J_MAT_LUMPED_END_TIME, 1.0_RKIND);
         if (res.Rtime_on < 0.0 || res.Rtime_off < 0.0) { 
            snprintf(errorMsg, BUFSIZE, "%s%d starting or end time is negative", errorMsgInit.c_str(), mA.materialId);
            WarnErrReport("", true);
         }
      } else if (model == J_MAT_LUMPED_MODEL_INDUCTOR) {
         res.inductor = true;
         res.L = this->getRealAt(matPtr.p, J_MAT_LUMPED_INDUCTANCE, found);
         if (!found) {
            snprintf(errorMsg, BUFSIZE, "%s%d inductance not found.", errorMsgInit.c_str(), mA.materialId);
            WarnErrReport(errorMsg, true);
         }
         res.R = this->getRealAt(matPtr.p, J_MAT_LUMPED_RESISTANCE, 0.0_RKIND);
      } else if (model == J_MAT_LUMPED_MODEL_CAPACITOR) {
         res.capacitor = true;
         res.C = this->getRealAt(matPtr.p, J_MAT_LUMPED_CAPACITANCE, found);
         if (!found) {
            snprintf(errorMsg, BUFSIZE, "%s%d capacitance not found.", errorMsgInit.c_str(), mA.materialId);
            WarnErrReport(errorMsg, true);
         }
         res.R = this->getRealAt(matPtr.p, J_MAT_LUMPED_RESISTANCE, found);
         if (!found) {
            snprintf(errorMsg, BUFSIZE, "%s%d resistance not found.", errorMsgInit.c_str(), mA.materialId);
            WarnErrReport(errorMsg, true);
         }
      } else {
         snprintf(errorMsg, BUFSIZE, "%s%d invalid model.", errorMsgInit.c_str(), mA.materialId);
         WarnErrReport(errorMsg, true);
      }

      return res;
   }

   bool containsCellRegionsWithType(const materialAssociation_t& mA, int cellType) {
      for (int e = 0; e < static_cast<int>(mA.elementIds.size()); e++) {
         cR = this->mesh.getCellRegion(mA.elementIds[e]);
         if (!cellRegionToCoords(cR, cellType).empty()) {
            return true;
         }
      }

      return false;
   }

   void matAssToCoords(materialAssociation_t& mA, int cellType, std::vector<coords_t>& res) {
      std::string tagName;
      std::vector<coords_t> newCoords;
      cell_region_t cR;
      int nCs;
      int e, jIni, jEnd;
      
      // Precount
      nCs = 0;
      for (e = 0; e < static_cast<int>(mA.elementIds.size()); e++) {
         cR = this->mesh.getCellRegion(mA.elementIds[e]);
         newCoords = cellRegionToCoords(cR, cellType);
         nCs += static_cast<int>(newCoords.size());
      }

      // Fills coords
      jIni = 0;
      res.resize(nCs);
      for (e = 0; e < static_cast<int>(mA.elementIds.size()); e++) {
         cR = this->mesh.getCellRegion(mA.elementIds[e]);
         tagName = this->buildTagName(mA.materialId, mA.elementIds[e]);
         newCoords = cellRegionToCoords(cR, cellType, tagName);
         if (newCoords.empty()) continue;
         jEnd = jIni + static_cast<int>(newCoords.size()) - 1;
         for (int k = 0; k < static_cast<int>(newCoords.size()); k++) {
            res[jIni + k] = newCoords[k];
         }
         jIni = jEnd + 1; 
      }
   }

   LossyThinSurfaces_t readLossyThinSurfaces() {
      LossyThinSurfaces_t res;
      std::vector<materialAssociation_t> mAs;
      json_value_ptr_t mat;
      int nLossySurfaces;
      bool found;
      int i, j, k;
      std::vector<coords_t>* cs;

      mAs = this->getMaterialAssociations({J_MAT_TYPE_MULTILAYERED_SURFACE});
      
      // Precounts
      nLossySurfaces = 0;
      for (i = 0; i < static_cast<int>(mAs.size()); i++) {
         this->matAssToCoords(mAs[i], CELL_TYPE_SURFEL, cs);
         if (!cs->empty()) nLossySurfaces++;
      }

      // Fills
      if (nLossySurfaces == 0) {
         res = emptyLossyThinSurfaces();
         return res;
      }

      res.cs.resize(nLossySurfaces);
      res.length = nLossySurfaces;
      res.length_max = nLossySurfaces;

      k = 0;
      for (i = 0; i < static_cast<int>(mAs.size()); i++) {
         this->matAssToCoords(mAs[i], CELL_TYPE_SURFEL, cs);
         if (cs->empty()) continue;
         res.cs[k] = readLossyThinSurface(mAs[i]);
         k++;
      }

      for (i = 0; i < nLossySurfaces; i++) {
         if (res.nC_max < static_cast<int>(res.cs[i].c.size())) {
            res.nC_max = static_cast<int>(res.cs[i].c.size());
         }
      }
      
      return res;
   }

   LossyThinSurface_t readLossyThinSurface(const materialAssociation_t& mA) {
      LossyThinSurface_t res;
      bool found, hasAbsPermittivity, hasAbsPermeability;
      const std::string errorMsgInit = "ERROR reading lossy thin surface: ";
      int i;
      json_value_ptr_t mat;
      json_value* layer;
      json_value* layers;
      
      this->matAssToCoords(mA, CELL_TYPE_SURFEL, res.c);
      res.nc = static_cast<int>(res.c.size());

      mat = this->matTable.getId(mA.materialId);
      res.files = trim(adjustl(this->getStrAt(mat.p, J_NAME, " ")));
      this->core.get(mat.p, J_MAT_MULTILAYERED_SURF_LAYERS, layers);

      res.numcapas = this->core.count(layers);
      res.sigma.resize(res.numcapas);
      res.eps.resize(res.numcapas);
      res.mu.resize(res.numcapas);
      res.sigmam.resize(res.numcapas);
      res.thk.resize(res.numcapas);
      res.sigma_devia.resize(res.numcapas);
      res.eps_devia.resize(res.numcapas);
      res.mu_devia.resize(res.numcapas);
      res.sigmam_devia.resize(res.numcapas);
      res.thk_devia.resize(res.numcapas);
      for (i = 0; i < res.numcapas; i++) {
         this->core.get_child(layers, i + 1, layer);
         res.sigma[i]  = this->getRealAt(layer, J_MAT_ELECTRIC_CONDUCTIVITY, 0.0_RKIND);
         res.sigmam[i] = this->getRealAt(layer, J_MAT_MAGNETIC_CONDUCTIVITY, 0.0_RKIND);
         res.eps[i]    = this->getRealAt(layer, J_MAT_ABS_PERMITTIVITY, hasAbsPermittivity);
         if (!hasAbsPermittivity) {
            res.eps[i] = this->getRealAt(layer, J_MAT_REL_PERMITTIVITY, 1.0_RKIND) * EPSILON_VACUUM;
         }
         res.mu[i]     = this->getRealAt(layer, J_MAT_ABS_PERMEABILITY, hasAbsPermeability);
         if (!hasAbsPermeability) {
            res.mu[i] = this->getRealAt(layer, J_MAT_REL_PERMEABILITY, 1.0_RKIND) * MU_VACUUM;
         }
         res.thk[i]    = this->getRealAt(layer, J_MAT_MULTILAYERED_SURF_THICKNESS, found);
         if (!found) {
            WarnErrReport(errorMsgInit + J_MAT_MULTILAYERED_SURF_THICKNESS + " in layer not found.", true);
         }
         res.sigma_devia[i] = 0.0_RKIND;
         res.eps_devia[i] = 0.0_RKIND;
         res.mu_devia[i] = 0.0_RKIND;
         res.sigmam_devia[i] = 0.0_RKIND;
         res.thk_devia[i] = 0.0_RKIND;
      }
      return res;
   }

   LossyThinSurfaces_t emptyLossyThinSurfaces() {
      LossyThinSurfaces_t res;
      res.cs.clear();
      res.length = 0;
      res.length_max = 0;
      res.nC_max = 0;
      return res;
   }

   PlaneWaves_t readPlanewaves() {
      PlaneWaves_t res;
      json_value* sources;
      std::vector<json_value_ptr_t> pws;
      int i;
      bool found;

      this->core.get(this->root, J_SOURCES, sources, found);
      
      if (!found) {
         res.collection.clear();
         res.nc = static_cast<int>(res.collection.size());
         res.nc_max = static_cast<int>(res.collection.size());
         return res;
      }

      pws = this->jsonValueFilterByKeyValue(sources, J_TYPE, J_SRC_TYPE_PW);

      res.collection.resize(pws.size());
      for (i = 0; i < static_cast<int>(pws.size()); i++) {
         res.collection[i] = readPlanewave(pws[i].p);
      }
      res.nc = static_cast<int>(res.collection.size());
      res.nc_max = static_cast<int>(res.collection.size());

      return res;
   }

   PlaneWave_t readPlanewave(json_value* pw) {
      PlaneWave_t res;
      std::string label;
      bool found;

      res.nombre_fichero = trim(adjustl(this->getStrAt(pw, J_SRC_MAGNITUDE_FILE)));

      res.atributo = "LOCKED";

      res.theta = this->getRealAt(pw, J_SRC_PW_DIRECTION + "." + J_SRC_PW_THETA);
      res.phi = this->getRealAt(pw, J_SRC_PW_DIRECTION + "." + J_SRC_PW_PHI);
      res.alpha = this->getRealAt(pw, J_SRC_PW_POLARIZATION + "." + J_SRC_PW_THETA);
      res.beta = this->getRealAt(pw, J_SRC_PW_POLARIZATION + "." + J_SRC_PW_PHI);

      {
         std::vector<cell_interval_t> cellIntervals;
         std::vector<coords_t> nfdeCoords;
         cellIntervals = this->getSingleVolumeInElementsIds(pw);
         if (cellIntervals.empty()) return res;
         nfdeCoords = cellIntervalsToCoords(cellIntervals);
         res.coor1 = {nfdeCoords[0].Xi, nfdeCoords[0].Yi, nfdeCoords[0].Zi};
         res.coor2 = {nfdeCoords[0].Xe, nfdeCoords[0].Ye, nfdeCoords[0].Ze};
      }

      res.isRC = false;
      res.nummodes = 1;
      res.incertmax = 0.0_RKIND;
      return res;
   }

   NodSource_t readNodalSources() {
      NodSource_t res;
      json_value* sources;
      std::vector<json_value_ptr_t> nodSrcs;
      bool found;
      int i;

      this->core.get(this->root, J_SOURCES, sources, found);
      if (!found) {
         res.NodalSource.clear();
         return res;
      }

      nodSrcs = this->jsonValueFilterByKeyValues(sources, J_TYPE, {J_SRC_TYPE_NS});
      if (nodSrcs.empty()) {
         res.NodalSource.clear();
         return res;
      }

      res.NodalSource.resize(nodSrcs.size());
      res.n_nodSrc = static_cast<int>(nodSrcs.size());
      res.n_nodSrc_max = static_cast<int>(nodSrcs.size());
      for (i = 0; i < static_cast<int>(nodSrcs.size()); i++) {
         res.NodalSource[i] = readField(nodSrcs[i].p);
      }
      for (i = 0; i < static_cast<int>(nodSrcs.size()); i++) {
         res.n_C1P_max = std::max(res.n_C1P_max, res.NodalSource[i].n_C1P);
         res.n_C2P_max = std::max(res.n_C2P_max, res.NodalSource[i].n_C2P);
      }
      if (res.n_nodSrc > 0) {
         if (res.n_C1P_max == 0) res.n_C1P_max = 1;
         if (res.n_C2P_max == 0) res.n_C2P_max = 1;
      }

      return res;
   }

   Curr_Field_Src_t readField(json_value* jns) {
      Curr_Field_Src_t res;
      json_value* entry;
      std::vector<int> elementIds;
      char nodalSourceName[BUFSIZE];
      std::vector<coords_scaled_t>* allCoords;
      int j, cnt_c1p, cnt_c2p;

      if (this->getStrAt(jns, J_FIELD, J_FIELD_CURRENT) == J_FIELD_CURRENT) {
         res.isElec = true;
      } else {
         WarnErrReport("Error reading current field source. Field label not recognized.", true);
      }
      
      if (this->getStrAt(jns, J_SRC_NS_HARDNESS, J_SRC_NS_HARDNESS_SOFT) == J_SRC_NS_HARDNESS_SOFT) {
         res.isHard = false;
      } else if (this->getStrAt(jns, J_SRC_NS_HARDNESS, J_SRC_NS_HARDNESS_SOFT) == J_SRC_NS_HARDNESS_HARD) {
         res.isHard = true;
      } else {
         WarnErrReport("Error reading current field source. Hardness label not recognized.", true);
      }
      
      res.isInitialValue = false;

      res.nombre = trim(adjustl(this->getStrAt(jns, J_SRC_MAGNITUDE_FILE)));

      std::string nodalSourceNameStr = this->getStrAt(jns, J_NAME, " ");
      strncpy(nodalSourceName, nodalSourceNameStr.c_str(), BUFSIZE - 1);
      nodalSourceName[BUFSIZE - 1] = '\0';

      elementIds = this->getIntsAt(jns, J_ELEMENTIDS);
      cellRegionsToScaledCoords(allCoords, this->mesh.getCellRegions(elementIds));

      cnt_c1p = 0;
      cnt_c2p = 0;
      for (j = 0; j < static_cast<int>(allCoords->size()); j++) {
         if (allCoords->at(j).Xi == allCoords->at(j).Xe &&
             allCoords->at(j).Yi == allCoords->at(j).Ye &&
             allCoords->at(j).Zi == allCoords->at(j).Ze) {
            cnt_c1p++;
         } else {
            cnt_c2p++;
         }
      }

      if (cnt_c1p > 0) res.c1P.resize(cnt_c1p);
      if (cnt_c2p > 0) res.c2P.resize(cnt_c2p);
      cnt_c1p = 0;
      cnt_c2p = 0;
      for (j = 0; j < static_cast<int>(allCoords->size()); j++) {
         if (allCoords->at(j).Xi == allCoords->at(j).Xe &&
             allCoords->at(j).Yi == allCoords->at(j).Ye &&
             allCoords->at(j).Zi == allCoords->at(j).Ze) {
            cnt_c1p++;
            res.c1P[cnt_c1p - 1] = allCoords->at(j);
         } else {
            cnt_c2p++;
            res.c2P[cnt_c2p - 1] = allCoords->at(j);
         }
      }
      res.n_C1P = cnt_c1p;
      res.n_C2P = cnt_c2p;

      delete allCoords;

      return res;
   }

   Sondas_t readProbes() {
      Sondas_t res;
      json_value* allProbes;
      std::vector<json_value_ptr_t> ps;
      // The only oldProbe present in the format is the far field.
      const char* validTypes[] = {J_PR_TYPE_FARFIELD};
      int i;
      bool found;

      this->core.get(this->root, J_PROBES, allProbes, found);
      if (!found) {
         res.probes.clear();
         res.n_probes = static_cast<int>(res.probes.size());
         res.n_probes_max = static_cast<int>(res.probes.size());
         return res;
      }

      ps = this->jsonValueFilterByKeyValues(allProbes, J_TYPE, validTypes, 1);

      res.n_probes = static_cast<int>(ps.size());
      res.n_probes_max = static_cast<int>(ps.size());
      res.probes.resize(ps.size());
      for (i = 0; i < static_cast<int>(ps.size()); i++) {
         res.probes[i] = readFarFieldProbe(ps[i].p);
      }

      return res;
   }

   abstractSonda_t readFarFieldProbe(json_value* p) {
      abstractSonda_t res;
      Sonda_t* ff;
      std::string outputName;
      bool transferFunctionFound;
      domain_t domain;

      res.n_FarField = 1;
      res.n_FarField_max = 1;
      res.FarField.resize(1);
      ff = &res.FarField[0].probe;

      ff->grname = " ";
      outputName = this->getStrAt(p, J_NAME);
      ff->outputrequest = trim(adjustl(outputName));

      // Far fields only accept frequency domains.
      domain = this->getDomain(p, J_PR_DOMAIN);
      if (domain.type2 != NP_T2_FREQ) {
         WarnErrReport("Only frequency domain is accepted for far field probes.", true);
      }
      ff->tstart = 0.0_RKIND;
      ff->tstop = 0.0_RKIND;
      ff->tstep = 0.0_RKIND;
      ff->fstart = domain.fstart;
      ff->fstop = domain.fstop;
      ff->fstep = domain.fstep;

      {
         bool sourcesFound;
         json_value* sources, * src;
         std::string fn;

         fn = this->getStrAt(p, J_PR_DOMAIN + J_PR_DOMAIN_MAGNITUDE_FILE, transferFunctionFound);
         if (!transferFunctionFound) {
            this->core.get(this->root, J_SOURCES, sources, sourcesFound);
            if (sourcesFound) {
               if (this->core.count(sources) == 1) {
                  this->core.get_child(sources, 1, src);
                  fn = this->getStrAt(src, J_SRC_MAGNITUDE_FILE, transferFunctionFound);
               }
            }
         }

         if (transferFunctionFound) {
            ff->FileNormalize = trim(adjustl(fn));
         } else {
            ff->FileNormalize = " ";
         }
      }

      if (domain.isLogarithmicFrequencySpacing) {
         appendLogSufix(ff->outputrequest);
      }

      {
         std::vector<coords_t> nfdeCoords;
         nfdeCoords = cellIntervalsToCoords(this->getSingleVolumeInElementsIds(p));
         ff->n_cord = 2;
         ff->n_cord_max = 2;
         ff->i.resize(2);
         ff->j.resize(2);
         ff->k.resize(2);
         ff->node.clear();
         ff->i[0] = nfdeCoords[0].Xi;
         ff->i[1] = nfdeCoords[0].Xe;
         ff->j[0] = nfdeCoords[0].Yi;
         ff->j[1] = nfdeCoords[0].Ye;
         ff->k[0] = nfdeCoords[0].Zi;
         ff->k[1] = nfdeCoords[0].Ze;
      }

      {
         readDirection(p, J_PR_FAR_FIELD_PHI, ff->phistart, ff->phistop, ff->phistep);
         readDirection(p, J_PR_FAR_FIELD_THETA, ff->thetastart, ff->thetastop, ff->thetastep);
      }
      return res;
   }

   void readDirection(json_value* p, const std::string& label, real_kind& initial, real_kind& final, real_kind& step) {
      json_value* dir;
      bool found;

      this->core.get(p, label, dir, found);
      if (!found) {
         WarnErrReport("Error reading far field probe. Direction label not found.", true);
      }
      initial = this->getRealAt(dir, J_PR_FAR_FIELD_DIR_INITIAL);
      final   = this->getRealAt(dir, J_PR_FAR_FIELD_DIR_FINAL);
      step    = this->getRealAt(dir, J_PR_FAR_FIELD_DIR_STEP);
   }

   MasSondas_t readMoreProbes() {
      MasSondas_t res;
      json_value* allProbes;
      std::vector<json_value_ptr_t> ps;

      int i;
      const char* validTypes[] = {J_PR_TYPE_POINT, J_PR_TYPE_LINE};
      bool found;
      std::string probeLbl;
      int filtered_size, n;

      this->core.get(this->root, J_PROBES, allProbes, found);
      if (!found) {
         res.collection.clear();
         res.length = static_cast<int>(res.collection.size());
         res.length_max = static_cast<int>(res.collection.size());
         res.len_cor_max = 0;
         return res;
      }

      ps = this->jsonValueFilterByKeyValues(allProbes, J_TYPE, validTypes, 2);
      
      filtered_size = 0;
      for (i = 0; i < static_cast<int>(ps.size()); i++) {
         if (isMoreProbe(ps[i].p)) { 
            filtered_size++;
         }
      }

      n = 0;
      res.collection.resize(filtered_size);
      for (i = 0; i < static_cast<int>(ps.size()); i++) {
         if (isMoreProbe(ps[i].p)) { 
            probeLbl = this->getStrAt(ps[i].p, J_TYPE, J_FIELD_ELECTRIC);
            if (probeLbl == J_PR_TYPE_POINT) { 
               res.collection[n] = readPointProbe(ps[i].p);
               n++;
            } else if (probeLbl == J_PR_TYPE_LINE) {
               res.collection[n] = readLineProbe(ps[i].p);
               n++;
            }
         }
      }

      res.length = static_cast<int>(res.collection.size());
      res.length_max = static_cast<int>(res.collection.size());
      for (i = 0; i < static_cast<int>(res.collection.size()); i++) {
          if (static_cast<int>(res.collection[i].cordinates.size()) > res.len_cor_max) {
              res.len_cor_max = static_cast<int>(res.collection[i].cordinates.size());
          }
      }
      return res;
   }

   bool isMoreProbe(json_value* p) {
      return isPointProbe(p) || isLineProbe(p);
   }

   bool isLineProbe(json_value* p) {
      return this->getStrAt(p, J_TYPE) == J_PR_TYPE_LINE;
   }

   bool isPointProbe(json_value* p) {
      std::string typeLabel, fieldLabel;
      bool found;

      typeLabel = this->getStrAt(p, J_TYPE, found);
      if (!found) {
         WarnErrReport("Point probe type label not found.", true);
      }
      if (typeLabel != J_PR_TYPE_POINT) {
         return false;
      }

      fieldLabel = this->getStrAt(p, J_FIELD, J_FIELD_ELECTRIC);
      if (fieldLabel == J_FIELD_ELECTRIC ||
            fieldLabel == J_FIELD_MAGNETIC) {
         return true;
      } else {
         return false;
      }
   }

   MasSonda_t readLineProbe(json_value* p) {
      MasSonda_t res;
      int i;
      std::string outputName;
      std::vector<linel_t> linels;
      polyline_t polyline;

      std::vector<int> elemIds;
      bool elementIdsFound, nameFound;

      outputName = this->getStrAt(p, J_NAME, nameFound);
      // Note: The Fortran code cuts off here. 
      // Assuming the rest of the function logic would follow similar patterns.
      // Since the input is truncated, we provide the translation up to the cut-off point.
      // In a real scenario, the rest of the function would be translated.
      
      // Placeholder for truncated code
      return res; 
   }

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Forward declarations and includes for types used in this chunk
// Assuming these are defined in previous chunks or headers
// #include "parser_t.h"
// #include "MasSonda_t.h"
// #include "BloqueProbes_t.h"
// #include "VolProbes_t.h"
// #include "ThinSlots_t.h"
// #include "ThinWires_t.h"
// #include "json_value.h"
// #include "mesh_t.h"
// #include "material_t.h"
// #include "constants.h" // For J_*, NP_*, iEx, etc.

// Helper function declarations (assumed to exist in the translation unit or headers)
void WarnErrReport(const std::string& msg, bool fatal);
std::string trim(const std::string& str);
std::string adjustl(const std::string& str);
int strToFieldType(const std::string& fieldLabel, const std::string& dirLabel);
std::vector<std::string> buildDirLabels(const json_value* dirLabelsPtr, const parser_t& core);
pixel_t getPixelFromElementId(const mesh_t* mesh, int elemId);
std::vector<coords_t> cellIntervalsToCoords(const std::vector<interval_t>& intervals);
int buildVolProbeType(const std::string& fieldType, const std::string& component);
std::vector<generator_description_t> readGeneratorOnThinWire(const std::vector<linel_t>& linels, const std::vector<int>& elementIds);
int getOrAssignNodeIndex(int coordId, int& nGlobal, int& nNodes, std::vector<int>& nodeCoordIds, std::vector<int>& nodeNodeIdx);

// Constants assumed to be defined elsewhere
// extern const int NP_COR_LINE;
// extern const int NP_COR_EX;
// ... etc

namespace parser_namespace { // Or whatever namespace the parser_t belongs to

    // ... (Previous code would be here)

    // Continuation of parser_t methods

    // Note: The Fortran code ends abruptly inside getOrAssignNodeIndex. 
    // I will complete the function and close the class/namespace structure appropriately 
    // based on the context of a class method.

    int getOrAssignNodeIndex(int coordId, int& nGlobal, int& nNodes, std::vector<int>& nodeCoordIds, std::vector<int>& nodeNodeIdx) {
        for (int k = 0; k < nNodes; ++k) {
            if (nodeCoordIds[k] == coordId) {
                return nodeNodeIdx[k];
            }
        }
        nGlobal = nGlobal + 1;
        nNodes = nNodes + 1;
        if (nNodes > static_cast<int>(nodeCoordIds.size())) {
            // This case should ideally not happen if pre-allocated correctly, 
            // but Fortran allocatable arrays grow. In C++, we might need to resize if logic differs.
            // Based on Fortran: allocate(nodeCoordIds(2 * res%n_tw))
            // If nNodes exceeds size, we have an issue. Assuming pre-allocation is sufficient.
            // If not, we'd need to push_back or resize.
            // For strict translation of logic:
            // nodeCoordIds.resize(nNodes + 1); 
            // nodeNodeIdx.resize(nNodes + 1);
            // But let's assume the pre-allocation in readThinWires was correct.
            throw std::runtime_error("Node index array overflow in getOrAssignNodeIndex");
        }
        nodeCoordIds[nNodes - 1] = coordId; // 0-based index for vector
        nodeNodeIdx[nNodes - 1] = nGlobal;
        return nGlobal;
    }

    // The Fortran code snippet ended inside the function. 
    // The rest of readThinWire and readThinWires would follow if they were in this chunk.
    // Since the chunk ends at `nodeNodeIdx(nNodes) = nGlobal`, I have completed that function.
    // The `readThinWire` function and `readThinWires` subroutine are not closed in the provided snippet,
    // but the snippet ends inside a helper function. 
    // I will output the completed helper function. 
    // If this is a continuation, the previous chunk likely defined the class and earlier methods.
    
    // Note: The provided Fortran code ends with `nodeNodeIdx(nNodes) = nGlobal` inside `getOrAssignNodeIndex`.
    // It does not show the closing braces for `readThinWire`, `readThinWires`, or the class.
    // I will provide the C++ translation for the code present in the chunk.

} // namespace parser_namespace

nodeIdx = nGlobal;
   }

   std::vector<generator_description_t> parser_t::readGeneratorOnThinWire(const std::vector<linel_t>& linels, const std::vector<int>& plineElemIds) {
      std::vector<generator_description_t> res;
      json_value* sources = nullptr;
      std::vector<json_value_ptr_t> genSrcs;
      bool found = false;
      int i;

      res.resize(linels.size());
      for (i = 0; i < linels.size(); i++) {
         res[i].srcfile = "None";
         res[i].srctype = "None";
         res[i].multiplier = 0.0_RKIND;
      }

      this->core->get(this->root, J_SOURCES, sources, found);
      if (!found) {
         return res;
      }

      genSrcs = this->jsonValueFilterByKeyValues(sources, J_TYPE, {J_SRC_TYPE_GEN});
      if (genSrcs.size() == 0) {
         return res;
      }

      {
         std::vector<int> sourceElemIds;
         int position;
         node_t srcCoord;
         polyline_t polylineCoords;
         for (i = 0; i < genSrcs.size(); i++) {
            sourceElemIds = this->getIntsAt(genSrcs[i].p, J_ELEMENTIDS);
            srcCoord = this->mesh->getNode(sourceElemIds[0]);
            polylineCoords = this->mesh->getPolyline(plineElemIds[0]);
            if (!std::any_of(polylineCoords.coordIds.begin(), polylineCoords.coordIds.end(), [&](int id){ return id == srcCoord.coordIds[0]; })) {
               continue; // generator is not in this polyline
            }

            position = findSourcePositionInLinels(sourceElemIds, linels);

            if (!this->existsAt(genSrcs[i].p, J_SRC_MAGNITUDE_FILE)) {
               WarnErrReport("magnitudeFile of source missing", true);
            }

            switch (this->getStrAt(genSrcs[i].p, J_FIELD)) {
               case J_FIELD_VOLTAGE:
                  res[position].srctype = F_SOURCE_VOLTAGE;
                  break;
               case J_FIELD_CURRENT:
                  res[position].srctype = F_SOURCE_CURRENT;
                  break;
               default:
                  WarnErrReport("Field block of source of type generator must be current or voltage", true);
            }
            res[position].srcfile = this->getStrAt(genSrcs[i].p, J_SRC_MAGNITUDE_FILE);
            res[position].multiplier = 1.0_RKIND * orientFieldFromGenerator(linels, position);

         }
      }
      return res;
   }

   int parser_t::orientFieldFromGenerator(const std::vector<linel_t>& linels, int position) {
      int res;
      if (position == 1) { 
         res = std::copysign(1.0, linels[position-1].orientation);
      } else if (position == linels.size()) { 
         res = -std::copysign(1.0, linels[position-1].orientation);
      } else {
         res = std::copysign(1.0, linels[position-1].orientation);
      }
      return res;
   }

   int parser_t::findSourcePositionInLinels(const std::vector<int>& srcElemIds, const std::vector<linel_t>& linels) {
      pixel_t pixel;
      int res;
      int i;
      
      pixel = this->mesh->nodeToPixel(this->mesh->getNode(srcElemIds[0]));
      for (i = 0; i < linels.size(); i++) {
         if (linels[i].tag == pixel.tag) {
            res = i + 1;
            return res;
         }
      }

      WarnErrReport("Source could not be found in linels.", true);
      return 0; // Should not reach here due to WarnErrReport
   }

   thinwiretermination_t parser_t::readWireTermination(json_value* terminal) {
      thinwiretermination_t res;
      json_value* tms = nullptr;
      json_value* tm = nullptr;
      std::string label;
      bool found = false;

      this->core->get(terminal, J_MAT_TERM_TERMINATIONS, tms, found);

      if (!found) {
         WarnErrReport("Error reading wire terminal. terminations not found.", true);
      }
      if (this->core->count(tms) != 1) {
         WarnErrReport("Only terminals with a single termination are allowed for a wire.", true);
      }

      this->core->get_child(tms, 1, tm);

      label = this->getStrAt(tm, J_TYPE, found);
      res.terminationType = strToTerminationType(label);
      if (!found) {
         WarnErrReport("Error reading wire terminal. termination must specify a type.", true);
      }

      if (label == J_MAT_TERM_TYPE_OPEN) {
         res.r = 0.0_RKIND;
         res.l = 0.0_RKIND;
         res.c = 0.0_RKIND;
      } else if (label == J_MAT_TERM_TYPE_SHORT) {
         res.r = 0.0_RKIND;
         res.l = 0.0_RKIND;
         res.c = 0.0_RKIND;
      } else if (label == J_MAT_TERM_TYPE_SERIES) {
         res.r = this->getRealAt(tm, J_MAT_TERM_RESISTANCE, 0.0_RKIND);
         res.l = this->getRealAt(tm, J_MAT_TERM_INDUCTANCE, 0.0_RKIND);
         res.c = this->getRealAt(tm, J_MAT_TERM_CAPACITANCE, 1e22_RKIND);
      } else {
         WarnErrReport("Error reading wire terminal. Holland wires only support open, short, and series terminations", true);
      }

      return res;
   }

   int parser_t::strToTerminationType(const std::string& label) {
      int res;
      if (label == J_MAT_TERM_TYPE_OPEN) {
         res = MATERIAL_CONS;
      } else if (label == J_MAT_TERM_TYPE_SERIES) {
         res = SERIES_CONS;
      } else if (label == J_MAT_TERM_TYPE_SHORT) {
         res = MATERIAL_CONS;
      } else {
         res = MATERIAL_CONS; // Default or handle error
      }
      return res;
   }

   bool parser_t::isThinWire(materialAssociation_t mA) {
      json_value_ptr_t mat;
      polyline_t pl;
      bool found = false;
      bool isThinWireRes = false;

      isThinWireRes = false;

      if (mA.elementIds.size() != 1) {
         WarnErrReport("Thin wires must be defined by a single element id.", true);
      }

      pl = this->mesh->getPolyline(mA.elementIds[0]);
      if (!this->mesh->arePolylineSegmentsStructured(pl)) {
         WarnErrReport("Thin wires must be defined by a structured polyline.", true);
      }
   
      isThinWireRes = true;
      return isThinWireRes;
   }

   MasSonda_t parser_t::readWireProbe(json_value* p) {
      MasSonda_t res;
      std::string outputName;
      std::string fieldLabel;
      node_t node;
      coordinate_t probe_coord;
      std::vector<int> elemIds;
      bool nameFound = false;
      bool elementIdsFound = false;

      outputName = this->getStrAt(p, J_NAME, nameFound);
      if (!nameFound) {
         WarnErrReport("Wire probes must define a name.", true);
      }
      res.outputrequest = outputName; // trim(adjustl(outputName)) equivalent in C++ string handling usually just copy if no leading/trailing spaces expected or use std::algorithm

      setDomainOfWireProbe(res, this->getDomain(p, J_PR_DOMAIN));

      elemIds = this->getIntsAt(p, J_ELEMENTIDS, elementIdsFound);
      if (!elementIdsFound) {
         WarnErrReport("Element ids entry not found for wire probe.", true);
      }
      if (elemIds.size() != 1) {
         WarnErrReport("Wire probe must contain a single element id.", true);
      }

      node = this->mesh->getNode(elemIds[0]);
      probe_coord = this->mesh->getCoordinate(node.coordIds[0]);
      fieldLabel = this->getStrAt(p, J_FIELD, "voltage"); // default=J_FIELD_VOLTAGE

      res.cordinates.resize(1);
      res.cordinates[0].tag = outputName;
      res.cordinates[0].Xi = getSegmentNdWhichMatchesCoord(node.coordIds[0], probe_coord);
      res.cordinates[0].Yi = 0;
      res.cordinates[0].Zi = 0;
      if (fieldLabel == "current") { // J_FIELD_CURRENT
         res.cordinates[0].Or = NP_COR_WIRECURRENT;
      } else if (fieldLabel == "voltage") { // J_FIELD_VOLTAGE
         res.cordinates[0].Or = NP_COR_DDP;
      } else {
         WarnErrReport("Invalid field label for wire probe.", true);
      }

      res.len_cor = 1;
      return res;
   }

   void parser_t::setDomainOfWireProbe(MasSonda_t& res, const domain_t& domain) {
      res.tstart = domain.tstart;
      res.tstep  = domain.tstep;
      res.tstop  = domain.tstop;
      res.fstart = domain.fstart;
      res.fstep  = domain.fstep;
      res.fstop  = domain.fstop;
      if (domain.filename != "") {
         res.filename = domain.filename;
      } else {
         res.filename = " ";
      }
      res.type1  = domain.type1;
      res.type2  = domain.type2;

      if (domain.isLogarithmicFrequencySpacing) {
         appendLogSufix(res.outputrequest);
      }
   }

   int parser_t::getSegmentNdWhichMatchesCoord(int coordId, const coordinate_t& probe_coord) {
      int nd_index;
      std::vector<linel_t> linels;
      polyline_t polyline;
      std::vector<coordinate_t> linelCoords;
      std::vector<double> distance_to_linel_cell;
      int k, j, i, or_val, n_linels, local_idx;
      std::vector<int> m(1);

      nd_index = 0;
      for (k = 0; k < mAs.size(); k++) {
         polyline = this->mesh->getPolyline(mAs[k].elementIds[0]);
         for (j = 0; j < polyline.coordIds.size(); j++) {
            if (polyline.coordIds[j] == coordId) {
               linels = this->mesh->polylineToLinels(polyline);
               n_linels = linels.size();
               linelCoords.resize(n_linels + 1);
               for (i = 0; i < n_linels; i++) {
                  linelCoords[i].position[0] = linels[i].cell[0];
                  linelCoords[i].position[1] = linels[i].cell[1];
                  linelCoords[i].position[2] = linels[i].cell[2];
                  if (linels[i].orientation < 0) {
                     or_val = std::abs(linels[i].orientation);
                     linelCoords[i].position[or_val-1] = linelCoords[i].position[or_val-1] + 1; // 1-based to 0-based adjustment if position is 1-based array
                  }
               }
               or_val = linels[n_linels-1].orientation;
               linelCoords[n_linels].position = linelCoords[n_linels-1].position;
               linelCoords[n_linels].position[std::abs(or_val)-1] = linelCoords[n_linels].position[std::abs(or_val)-1] + (or_val > 0 ? 1 : -1);
               distance_to_linel_cell.resize(n_linels + 1);
               for (i = 0; i < n_linels + 1; i++) {
                  distance_to_linel_cell[i] = std::sqrt(
                     std::pow(linelCoords[i].position[0] - probe_coord.position[0], 2) +
                     std::pow(linelCoords[i].position[1] - probe_coord.position[1], 2) +
                     std::pow(linelCoords[i].position[2] - probe_coord.position[2], 2)
                  );
               }
               // minloc returns 1-based index in Fortran
               int min_idx = 0;
               double min_val = distance_to_linel_cell[0];
               for (i = 1; i < distance_to_linel_cell.size(); i++) {
                  if (distance_to_linel_cell[i] < min_val) {
                     min_val = distance_to_linel_cell[i];
                     min_idx = i;
                  }
               }
               m[0] = min_idx + 1; // Convert to 1-based for logic below if needed, but local_idx logic handles it
               local_idx = std::min(m[0], n_linels);
               nd_index = res.tw[k].twc[local_idx-1].nd; // Adjusting for 0-based vector access if twc is vector
               return nd_index;
            }
         }
      }
      return nd_index;
   }

   domain_t parser_t::getDomain(json_value* place, const std::string& path) {
      domain_t res;
      json_value* domain = nullptr;
      std::string fn;
      std::string domainType;
      std::string freqSpacing;
      bool found = false;
      bool transferFunctionFound = false;
      double val;

      this->core->get(place, path, domain, found);
      if (!found) {
         res.filename = " ";
         return res;
      }

      fn = this->getStrAt(domain, J_PR_DOMAIN_MAGNITUDE_FILE, transferFunctionFound, " ");
      if (transferFunctionFound) {
         res.filename = fn; // trim(adjustl(fn))
      }

      res.type1 = NP_T1_PLAIN;

      domainType = this->getStrAt(domain, J_TYPE, "time"); // default=J_PR_DOMAIN_TYPE_TIME
      res.type2 = getNPDomainType(domainType, transferFunctionFound);

      res.tstart = this->getRealAt(domain, J_PR_DOMAIN_TIME_START, 0.0_RKIND);
      res.tstop = this->getRealAt(domain, J_PR_DOMAIN_TIME_STOP, 0.0_RKIND);
      res.tstep = this->getRealAt(domain, J_PR_DOMAIN_TIME_STEP, 0.0_RKIND);
      if (res.tstart < 0.0 || res.tstop < 0.0 || res.tstep < 0.0) {
         std::string errorMsg;
         std::string p_name;
         bool nameFound = false;
         p_name = this->getStrAt(place, J_NAME, nameFound);
         errorMsg = "Probe named " + p_name + " has negative times in its domain definition";
         WarnErrReport(errorMsg, true);
      }
      res.fstart = this->getRealAt(domain, J_PR_DOMAIN_FREQ_START, 0.0_RKIND);
      res.fstop = this->getRealAt(domain, J_PR_DOMAIN_FREQ_STOP, 0.0_RKIND);

      int numberOfFrequencies = this->getIntAt(domain, J_PR_DOMAIN_FREQ_NUMBER, 0);
      if (numberOfFrequencies == 0) {
         res.fstep = 0.0_RKIND;
      } else {
         res.fstep = (res.fstop - res.fstart) / numberOfFrequencies;
      }

      freqSpacing = this->getStrAt(domain, J_PR_DOMAIN_FREQ_SPACING, "linear"); // default=J_PR_DOMAIN_FREQ_SPACING_LINEAR
      if (freqSpacing == "linear") { // J_PR_DOMAIN_FREQ_SPACING_LINEAR
         res.isLogarithmicFrequencySpacing = false;
      } else if (freqSpacing == "logarithmic") { // J_PR_DOMAIN_FREQ_SPACING_LOGARITHMIC
         res.isLogarithmicFrequencySpacing = true;
      }

      return res;
   }

   int parser_t::getNPDomainType(const std::string& typeLabel, bool hasTransferFunction) {
      int res;
      bool isTime = false;
      bool isFrequency = false;
      std::string errorMsg;

      if (typeLabel == "time") { // J_PR_DOMAIN_TYPE_TIME
         isTime = true;
         isFrequency = false;
      } else if (typeLabel == "freq") { // J_PR_DOMAIN_TYPE_FREQ
         isTime = false;
         isFrequency = true;
      } else if (typeLabel == "timefreq") { // J_PR_DOMAIN_TYPE_TIMEFREQ
         isTime = true;
         isFrequency = true;
      }

      if (isTime && !isFrequency && !hasTransferFunction) {
         res = NP_T2_TIME;
         return res;
      } else if (!isTime && isFrequency && !hasTransferFunction) {
         res = NP_T2_FREQ;
         return res;
      } else if (!isTime && !isFrequency && hasTransferFunction) {
         res = NP_T2_TRANSFER;
         return res;
      } else if (isTime && isFrequency && !hasTransferFunction) {
         res = NP_T2_TIMEFREQ;
         return res;
      } else if (isTime && !isFrequency && hasTransferFunction) {
         res = NP_T2_TIMETRANSF;
         return res;
      } else if (!isTime && isFrequency && hasTransferFunction) {
         res = NP_T2_FREQTRANSF;
         return res;
      } else if (isTime && isFrequency && hasTransferFunction) {
         res = NP_T2_TIMEFRECTRANSF;
         return res;
      }

      errorMsg = "Invalid domain type: " + typeLabel;
      WarnErrReport(errorMsg, true);
      return NP_T2_TIME; // Default
   }

   materialAssociation_t parser_t::parseMaterialAssociation(json_value* matAss) {
      materialAssociation_t res;
      bool found = false;
      bool isMultiwire = false;
      bool isWireOrMultiwire = false;
      std::string errorMsg = "ERROR reading material association: ";
      
      // Fills material association.
      res.materialId = this->getIntAt(matAss, J_MATERIAL_ID, found);
      if (!found) showLabelNotFoundError(J_MATERIAL_ID);
      
      res.elementIds = this->getIntsAt(matAss, J_ELEMENTIDS, found);
      if (!found) showLabelNotFoundError(J_ELEMENTIDS);
      
      res.name = this->getStrAt(matAss, J_NAME, found);
      if (!found) {
         res.name = "";
      }

      res.initialTerminalId  = this->getIntAt(matAss, J_MAT_ASS_CAB_INI_TERM_ID, -1);
      res.endTerminalId      = this->getIntAt(matAss, J_MAT_ASS_CAB_END_TERM_ID, -1);
      res.initialConnectorId = this->getIntAt(matAss, J_MAT_ASS_CAB_INI_CONN_ID, -1);
      res.endConnectorId     = this->getIntAt(matAss, J_MAT_ASS_CAB_END_CONN_ID, -1);
      res.containedWithinElementId = this->getIntAt(matAss, J_MAT_ASS_CAB_CONTAINED_WITHIN_ID, -1);

      res.hasTotalResistance = this->existsAt(matAss, J_MAT_ASS_TOTAL_RESISTANCE);
      if (res.hasTotalResistance) {
         if (this->dimensionAt(matAss, J_MAT_ASS_TOTAL_RESISTANCE) == 0) {
            res.totalResistance.resize(1);
            res.totalResistance[0] = this->getRealAt(matAss, J_MAT_ASS_TOTAL_RESISTANCE);
         } else {
            res.totalResistance = this->getRealsAt(matAss, J_MAT_ASS_TOTAL_RESISTANCE);
         }
      }

      // Checks validity of associations.
      if (this->matTable->checkId(res.materialId) != 0) {
         errorMsg = "ERROR reading material association: material with id " + std::to_string(res.materialId) + " not found.";
         WarnErrReport(errorMsg, true);
      }
      
      if (res.elementIds.size() == 0) {
         errorMsg = "ERROR reading material association: " + std::string(J_ELEMENTIDS) + " must not be empty.";
         WarnErrReport(errorMsg, true);
      }
      {
         int i;
         for (i = 0; i < res.elementIds.size(); i++) {
            if (this->mesh->checkElementId(res.elementIds[i]) != 0) {
               errorMsg = "ERROR reading material association: element with id " + std::to_string(res.elementIds[i]) + " not found.";
               WarnErrReport(errorMsg, true);
            }
         }
      }

      json_value_ptr_t mat = this->matTable->getId(res.materialId);
      std::string matType = this->getStrAt(mat.p, J_TYPE);
      isMultiwire = (matType == J_MAT_TYPE_SHIELDED_MULTIWIRE || matType == J_MAT_TYPE_UNSHIELDED_MULTIWIRE);
      isWireOrMultiwire = (matType == J_MAT_TYPE_WIRE || isMultiwire); 
      
      if (isWireOrMultiwire) {
         if (res.initialTerminalId == -1 || res.endTerminalId == -1) {
            errorMsg = "ERROR reading material association: wire associations must include terminals.";
            WarnErrReport(errorMsg, true);
         }
         if (!isMaterialIdOfType(res.initialTerminalId, J_MAT_TYPE_TERMINAL)) {
            errorMsg = "ERROR reading material association: material with id " + std::to_string(res.materialId) + " must be a terminal.";
            WarnErrReport(errorMsg, true);
         }
         if (!isMaterialIdOfType(res.endTerminalId, J_MAT_TYPE_TERMINAL)) {
            errorMsg = "ERROR reading material association: material with id " + std::to_string(res.materialId) + " must be a terminal.";
            WarnErrReport(errorMsg, true);
         }
         if (res.initialConnectorId != -1) {
            if (!isMaterialIdOfType(res.initialConnectorId, J_MAT_TYPE_CONNECTOR)) {
               errorMsg = "ERROR reading material association: material with id " + std::to_string(res.materialId) + " must be a connector.";
               WarnErrReport(errorMsg, true);
            }
         } 
         if (res.endConnectorId != -1) {
            if (!isMaterialIdOfType(res.endConnectorId, J_MAT_TYPE_CONNECTOR)) {
               errorMsg = "ERROR reading material association: material with id " + std::to_string(res.materialId) + " must be a connector.";
               WarnErrReport(errorMsg, true);
            }
         }
      }
      if (isMultiwire) {
         // Not defininign a containedWithinElementId is an error if the simulation is a 3D-FDTD one. 
         // For pure MTLN mode it is not an error.
         // if (res.containedWithinElementId == -1) then
         //    write(error_unit, *) errorMsgInit, "multiwire associations must include: ", J_MAT_ASS_CAB_CONTAINED_WITHIN_ID
         // end if
         // if (.not. (this%getLogicalAt(this%root, J_GENERAL//'.'//J_GEN_MTLN_PROBLEM, default = .false.)) .and. &
         //    (this%mesh%checkElementId(res%containedWithinElementId) /= 0)) then
         //    write(error_unit, *) errorMsgInit, "element with id ", res%containedWithinElementId, " not found."
         // end if
      }
      
      return res;
   }

   std::vector<materialAssociation_t> parser_t::getMaterialAssociations(const std::vector<std::string>& materialTypes, const std::vector<std::string>* elementLabels) {
      std::vector<materialAssociation_t> res;
      json_value* allMatAss = nullptr;
      json_value* mAPtr = nullptr;
      int i, j, k, e;
      int nMaterials = 0;
      bool found = false;
      this->core->get(this->root, J_MATERIAL_ASSOCIATIONS, allMatAss, found);
      if (!found) {
         res.resize(0);
         return res;
      }

      for (i = 0; i < this->core->count(allMatAss); i++) {
         this->core->get_child(allMatAss, i + 1, mAPtr); // Fortran 1-based indexing for get_child
         for (j = 0; j < materialTypes.size(); j++) {
            if (isAssociatedWithMaterial(mAPtr, materialTypes[j])) {
               if (elementLabels != nullptr) { 
                  if (isAssociatedWithElementLabel(mAPtr, *elementLabels)) nMaterials++;
               } else {
                  nMaterials++;
               }
            }
         }
      }

      res.resize(nMaterials);
      j = 0;
      for (i = 0; i < this->core->count(allMatAss); i++) {
         this->core->get_child(allMatAss, i + 1, mAPtr);
         for (k = 0; k < materialTypes.size(); k++) {
            if (isAssociatedWithMaterial(mAPtr, materialTypes[k])) {
               if (elementLabels != nullptr) { 
                  if (isAssociatedWithElementLabel(mAPtr, *elementLabels)) { 
                     res[j] = this->parseMaterialAssociation(mAPtr);
                     res[j].matAssType = materialTypes[k];
                     j++;
                  }
               } else {
                  res[j] = this->parseMaterialAssociation(mAPtr);
                  res[j].matAssType = materialTypes[k];
                  j++;
               }
            }
         }
      }

      return res;
   }

   std::string parser_t::buildTagName(int matId, int elementId) {
      std::string res;
      // Implementation missing in Fortran chunk, just returning empty or placeholder
      return "";
   }

#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cctype>

// Forward declarations and includes for types used in this chunk
// Assuming these are defined in previous chunks or headers
// #include "parser_t.h"
// #include "mtln_t.h"
// #include "cable_t.h"
// #include "materialAssociation_t.h"
// #include "connector_t.h"
// #include "terminal_network_t.h"
// #include "network_circuit_t.h"
// #include "aux_node_t.h"
// #include "coordinate_t.h"
// #include "polyline_t.h"
// #include "json_value_ptr_t.h"
// #include "json_value.h"
// #include "fhash_tbl_t.h"
// #include "cable_abstract_t.h"
// #include "shielded_multiwire_t.h"
// #include "terminal_connection_t.h"
// #include "node_t.h"
// #include "terminal_model_t.h"
// #include "generator_t.h"

// Constants and Enums assumed to be defined elsewhere
// enum class TerminalNodeSide { TERMINAL_NODE_SIDE_INI, TERMINAL_NODE_SIDE_END };
// enum class TerminationType { TERMINATION_SHORT, TERMINATION_OPEN, TERMINATION_NETWORK, ... };
// const char* J_NAME = "name";
// const char* J_TYPE = "type";
// const char* J_ID = "id";
// const char* J_MAT_TYPE_SHIELDED_MULTIWIRE = "SHIELDED_MULTIWIRE";
// const char* J_MAT_TYPE_UNSHIELDED_MULTIWIRE = "UNSHIELDED_MULTIWIRE";
// const char* J_MAT_TYPE_WIRE = "WIRE";
// const char* J_MAT_TYPE_CONNECTOR = "CONNECTOR";
// const char* J_MAT_CONN_RESISTANCES = "resistances";
// const char* J_MAT_CONN_TRANSFER_IMPEDANCES = "transfer_impedances";
// const char* J_GENERAL = "general";
// const char* J_GEN_TIME_STEP = "time_step";
// const char* J_GEN_NUMBER_OF_STEPS = "number_of_steps";
// const char* J_MATERIALS = "materials";
// const char* J_MAT_TERM_TERMINATIONS = "terminations";
// const char* J_MAT_TERM_CAPACITANCE = "capacitance";
// const char* J_MAT_TERM_RESISTANCE = "resistance";
// const char* J_MAT_TERM_INDUCTANCE = "inductance";
// const char* TAG_MATERIAL = "Material";
// const char* TAG_LAYER = "Layer";
// const int BUFSIZE = 256; // Example size
// const double RKIND = 1.0; // Example kind

// Helper functions assumed to be available
// std::string adaptName(const std::string& str);
// void WarnErrReport(const std::string& msg, bool isFatal);
// std::string intToStr(int i);
// std::string sideToStr(int side); // Assuming side is mapped to string
// std::vector<std::string> splitLineIntoWords(const std::string& line);
// std::string to_upper(const std::string& s);
// double readTransferImpedance(const json_value& z);
// std::vector<int> getCableElemIds(const json_value_ptr_t& cable);
// void addConnIdToConnectorMap(fhash_tbl_t& map, const std::vector<connector_t>& connectors);
// void addElemIdToCableMap(fhash_tbl_t& map, const std::vector<int>& elemIds, int idx);
// void addElemIdToPositionMap(fhash_tbl_t& map, const std::vector<int>& elemIds);
// std::vector<connector_t> readConnectors(); // This is actually a member function in the block below
// std::vector<shielded_multiwire_t> readWireGenerators(); // Placeholder
// std::vector<probe_t> readMultiwireProbes(); // Placeholder
// bool touchesOtherWire(int cId, int ownElemId);
// bool touchesNonVacuumMaterial(int cId, const coordinate_t& relPos);

namespace parser_namespace { // Replace with actual namespace if known

// Assuming parser_t is a class
// class parser_t {
// public:
//     // ... members ...
//     std::string getStrAt(const json_value_ptr_t& ptr, const char* key, bool& found);
//     json_value_ptr_t getId(int id);
//     std::vector<materialAssociation_t> getMaterialAssociations(const std::vector<std::string>& types);
//     double getRealAt(const json_value_ptr_t& ptr, const char* key, double defaultVal = 0.0);
//     int getIntAt(const json_value_ptr_t& ptr, const char* key);
//     bool existsAt(const json_value_ptr_t& ptr, const char* key);
//     std::vector<json_value_ptr_t> jsonValueFilterByKeyValue(const json_value& val, const char* key, const char* value);
//     void get(const json_value_ptr_t& ptr, const char* key, json_value& out, bool& found);
//     void get_child(const json_value& val, int index, json_value& out);
//     int count(const json_value& val);
//     mesh_t* mesh; // Assuming mesh is a member
//     json_value_ptr_t root;
//     json_core_t* core; // Assuming core is a member
//     material_table_t* matTable; // Assuming matTable is a member
//     // ...
// };

// Helper to simulate Fortran string concatenation and trimming
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, (last - first + 1));
}

std::string adjustl(const std::string& str) {
    return trim(str); // Simplified, Fortran adjustl only removes leading spaces
}

std::string trim_left(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    return str.substr(first);
}

std::string trim_right(const std::string& str) {
    size_t last = str.find_last_not_of(" \t\n\r");
    if (last == std::string::npos) return "";
    return str.substr(0, last + 1);
}

std::string trim_all(const std::string& str) {
    return trim(str);
}

// Function to check if a name is valid
void checkIsValidName(const std::string& str, std::string& errorMsg) {
    std::string notAllowedChars = "@";
    for (size_t i = 0; i < notAllowedChars.length(); ++i) {
        if (str.find(notAllowedChars[i]) != std::string::npos) {
            std::ostringstream oss;
            oss << "ERROR in name: " << str << " contains invalid character " << notAllowedChars[i];
            errorMsg = oss.str();
            WarnErrReport(errorMsg, true);
            return;
        }
    }
}

// Function to adapt name
std::string adaptName(const std::string& str) {
    std::string res = trim_left(str);
    for (size_t i = 0; i < res.length(); ++i) {
        if (res[i] == ' ') {
            res[i] = '_';
        }
    }
    return res;
}

// Assuming this is a member function of parser_t
// std::string parser_t::getLayerName(int matId, int elementId) {
//     std::string res;
//     std::string matName;
//     std::string layerName;
//     bool found;
//     std::string errorMsg;
//
//     {
//         json_value_ptr_t mat = this->matTable->getId(matId);
//         matName = this->getStrAt(mat, J_NAME, found);
//         if (!found) {
//             matName = std::string(TAG_MATERIAL) + std::to_string(matId);
//         }
//         matName = adaptName(matName);
//     }
//
//     {
//         json_value_ptr_t elem = this->elementTable->getId(elementId);
//         layerName = this->getStrAt(elem, J_NAME, found);
//         if (!found) {
//             layerName = std::string(TAG_LAYER) + std::to_string(elementId);
//         }
//         layerName = adaptName(layerName);
//     }
//
//     checkIsValidName(matName, errorMsg);
//     checkIsValidName(layerName, errorMsg);
//     res = matName + "@" + layerName;
//     return res;
// }

// Assuming this is a member function of parser_t
// mtln_t parser_t::readMTLN() {
//     mtln_t mtln_res;
//     fhash_tbl_t elemIdToPosition;
//     fhash_tbl_t elemIdToCable;
//     fhash_tbl_t connIdToConnector;
//     std::vector<materialAssociation_t> cables;
//
//     std::vector<std::string> types;
//     types.push_back(std::string(J_MAT_TYPE_SHIELDED_MULTIWIRE) + "  ");
//     types.push_back(std::string(J_MAT_TYPE_UNSHIELDED_MULTIWIRE) + "    ");
//     types.push_back(std::string(J_MAT_TYPE_WIRE) + "               ");
//     cables = this->getMaterialAssociations(types);
//
//     mtln_res.connectors = readConnectors();
//     addConnIdToConnectorMap(connIdToConnector, mtln_res.connectors);
//
//     if (cables.empty()) {
//         mtln_res.time_step = 0;
//         mtln_res.number_of_steps = 0;
//         mtln_res.cables = std::vector<std::unique_ptr<cable_abstract_t>>(0);
//         mtln_res.probes = std::vector<probe_t>(0);
//         mtln_res.networks = std::vector<terminal_network_t>(0);
//         return mtln_res;
//     }
//
//     {
//         std::vector<materialAssociation_t> unshielded = this->getMaterialAssociations({J_MAT_TYPE_UNSHIELDED_MULTIWIRE, J_MAT_TYPE_WIRE + "             "});
//         mtln_res.n_unsh = unshielded.size();
//         std::vector<materialAssociation_t> shielded = this->getMaterialAssociations({J_MAT_TYPE_SHIELDED_MULTIWIRE});
//         mtln_res.n_sh = shielded.size();
//     }
//
//     mtln_res.time_step = this->getRealAt(this->root, std::string(J_GENERAL) + "." + J_GEN_TIME_STEP, 0.0);
//     mtln_res.number_of_steps = this->getRealAt(this->root, std::string(J_GENERAL) + "." + J_GEN_NUMBER_OF_STEPS);
//
//     mtln_res.cables.resize(cables.size());
//     for (int i = 0; i < cables.size(); ++i) {
//         std::unique_ptr<cable_t> read_cable = readMTLNCable(cables[i]);
//         stopOnRepeteadName(*read_cable, mtln_res.cables, i);
//         mtln_res.cables[i].ptr = std::move(read_cable);
//         addElemIdToCableMap(elemIdToCable, cables[i].elementIds, i);
//         addElemIdToPositionMap(elemIdToPosition, cables[i].elementIds);
//     }
//
//     for (int i = 0; i < cables.size(); ++i) {
//         cable_abstract_t* ptr = mtln_res.cables[i].ptr.get();
//         if (auto* shielded_ptr = dynamic_cast<shielded_multiwire_t*>(ptr)) {
//             shielded_ptr->parent_cable = assignParentCable(cables[i], mtln_res.cables);
//             shielded_ptr->conductor_in_parent = assignConductorInParent(cables[i]);
//         }
//     }
//
//     mtln_res.wireGenerators = readWireGenerators();
//     mtln_res.probes = readMultiwireProbes();
//     mtln_res.networks = buildNetworks();
//
//     return mtln_res;
// }

// Helper functions for readMTLN

// std::unique_ptr<cable_t> parser_t::readMTLNCable(const materialAssociation_t& cable) {
//     // Implementation depends on specific cable types
//     // This is a placeholder
//     return nullptr;
// }

// std::unique_ptr<cable_abstract_t> parser_t::assignParentCable(const materialAssociation_t& cable, const std::vector<std::unique_ptr<cable_abstract_t>>& cables) {
//     json_value_ptr_t mat = this->matTable->getId(cable.materialId);
//     std::string type = this->getStrAt(mat, J_TYPE);
//
//     if (type == J_MAT_TYPE_SHIELDED_MULTIWIRE) {
//         int parentId = cable.containedWithinElementId;
//         if (parentId == -1) {
//             return nullptr;
//         } else {
//             return getPointerToParentCable(cables, parentId);
//         }
//     } else if (type == J_MAT_TYPE_UNSHIELDED_MULTIWIRE) {
//         return nullptr;
//     } else if (type == J_MAT_TYPE_WIRE) {
//         return nullptr;
//     } else {
//         WarnErrReport("ERROR: Material type not recognized", true);
//         return nullptr;
//     }
// }

// int parser_t::assignConductorInParent(const materialAssociation_t& cable) {
//     json_value_ptr_t mat = this->matTable->getId(cable.materialId);
//     std::string type = this->getStrAt(mat, J_TYPE);
//
//     if (type == J_MAT_TYPE_SHIELDED_MULTIWIRE) {
//         int parentId = cable.containedWithinElementId;
//         if (parentId == -1) {
//             return 0;
//         } else {
//             return getParentPositionInMultiwire(parentId);
//         }
//     } else if (type == J_MAT_TYPE_UNSHIELDED_MULTIWIRE) {
//         return 0;
//     } else if (type == J_MAT_TYPE_WIRE) {
//         return 0;
//     } else {
//         WarnErrReport("ERROR: Material type not recognized", true);
//         return 0;
//     }
// }

// void parser_t::stopOnRepeteadName(const cable_t& cable, const std::vector<std::unique_ptr<cable_abstract_t>>& cables, int n) {
//     bool unique = true;
//     for (int i = 0; i < n; ++i) {
//         if (cable.name == cables[i].ptr->name) {
//             unique = false;
//             break;
//         }
//     }
//     if (!unique) {
//         std::ostringstream oss;
//         oss << "Cable name " << cable.name << " has already been used";
//         WarnErrReport(oss.str(), true);
//     }
// }

// std::vector<connector_t> parser_t::readConnectors() {
//     json_value* mat = nullptr;
//     bool materialsFound = false;
//     this->core->get(this->root, J_MATERIALS, *mat, materialsFound);
//
//     if (!materialsFound) {
//         return std::vector<connector_t>(0);
//     }
//
//     std::vector<json_value_ptr_t> connectors = this->jsonValueFilterByKeyValue(*mat, J_TYPE, J_MAT_TYPE_CONNECTOR);
//     std::vector<connector_t> res(connectors.size());
//
//     if (!connectors.empty()) {
//         for (size_t i = 0; i < connectors.size(); ++i) {
//             res[i].id = this->getIntAt(connectors[i], J_ID);
//             if (this->existsAt(connectors[i], J_MAT_CONN_RESISTANCES)) {
//                 res[i].resistances = this->getRealsAt(connectors[i], J_MAT_CONN_RESISTANCES);
//             } else {
//                 res[i].resistances = std::vector<double>(0);
//             }
//
//             if (this->existsAt(connectors[i], J_MAT_CONN_TRANSFER_IMPEDANCES)) {
//                 json_value zs;
//                 this->core->get(connectors[i], J_MAT_CONN_TRANSFER_IMPEDANCES, zs);
//                 int n = this->core->count(zs);
//                 res[i].transfer_impedances_per_meter.resize(n);
//                 for (int j = 0; j < n; ++j) {
//                     json_value z;
//                     this->core->get_child(zs, j, z);
//                     res[i].transfer_impedances_per_meter[j] = readTransferImpedance(z);
//                 }
//             } else {
//                 res[i].transfer_impedances_per_meter = std::vector<double>(0);
//             }
//         }
//     }
//
//     return res;
// }

// int parser_t::findMaxElemId(const std::vector<json_value_ptr_t>& cables) {
//     int res = 0;
//     if (!cables.empty()) {
//         for (size_t i = 0; i < cables.size(); ++i) {
//             std::vector<int> elemIds = getCableElemIds(cables[i]);
//             if (elemIds.empty()) return 0;
//             int m = *std::max_element(elemIds.begin(), elemIds.end());
//             if (m > res) {
//                 res = m;
//             }
//         }
//     }
//     return res;
// }

// std::vector<terminal_network_t> parser_t::buildNetworks() {
//     std::vector<aux_node_t> aux_nodes;
//     std::vector<coordinate_t> networks_coordinates;
//     std::vector<materialAssociation_t> cables;
//
//     cables.push_back(this->getMaterialAssociations({J_MAT_TYPE_UNSHIELDED_MULTIWIRE}));
//     cables.push_back(this->getMaterialAssociations({J_MAT_TYPE_SHIELDED_MULTIWIRE}));
//     cables.push_back(this->getMaterialAssociations({J_MAT_TYPE_WIRE}));
//
//     for (size_t i = 0; i < cables.size(); ++i) {
//         std::vector<int> elemIds = cables[i].elementIds;
//         json_value_ptr_t cableMat = this->matTable->getId(cables[i].materialId);
//         std::string cableType = this->getStrAt(cableMat, J_TYPE);
//         bool isShieldedCable = (cableType == J_MAT_TYPE_SHIELDED_MULTIWIRE);
//
//         json_value* terminations_ini = getTerminationsOnSide(cables[i].initialTerminalId);
//         json_value* terminations_end = getTerminationsOnSide(cables[i].endTerminalId);
//
//         for (size_t j = 0; j < elemIds.size(); ++j) {
//             aux_nodes.push_back(buildNode(*terminations_ini, TERMINAL_NODE_SIDE_INI, j + 1, elemIds[j], isShieldedCable));
//             aux_nodes.push_back(buildNode(*terminations_end, TERMINAL_NODE_SIDE_END, j + 1, elemIds[j], isShieldedCable));
//             updateListOfNetworksCoordinates(networks_coordinates, elemIds[j]);
//         }
//     }
//
//     std::vector<terminal_network_t> res(networks_coordinates.size());
//     for (size_t i = 0; i < networks_coordinates.size(); ++i) {
//         res[i] = buildNetwork(networks_coordinates[i], aux_nodes, i + 1);
//     }
//
//     return res;
// }

// std::vector<coordinate_t> parser_t::buildListOfCoordinates(const std::vector<int>& elemIds) {
//     std::vector<coordinate_t> res;
//     for (int id : elemIds) {
//         updateListOfNetworksCoordinates(res, id);
//     }
//     return res;
// }

// terminal_network_t parser_t::buildNetwork(const coordinate_t& network_coordinate, const std::vector<aux_node_t>& aux_nodes, int network_index) {
//     std::vector<aux_node_t> network_nodes = filterNetworkNodesByCoordinate(aux_nodes, network_coordinate);
//     std::vector<int> node_ids = buildListOfNodeIds(network_nodes);
//     std::vector<network_circuit_t> network_circuits = buildNetworkCircuits(network_nodes, node_ids, network_index);
//
//     terminal_network_t res;
//     for (int node_id : node_ids) {
//         res.add_connection(buildConnection(node_id, network_nodes, network_circuits));
//     }
//
//     return res;
// }

// std::vector<network_circuit_t> parser_t::buildNetworkCircuits(const std::vector<aux_node_t>& nodes, const std::vector<int>& node_ids, int network_index) {
//     std::vector<aux_node_t> subckt_filtered_nodes = filterNetworkNodesByNetworkCircuit(nodes);
//     int n = 0;
//     for (int node_id : node_ids) {
//         std::vector<aux_node_t> id_filtered_nodes = filterNetworkNodesById(subckt_filtered_nodes, node_id);
//         if (!id_filtered_nodes.empty()) n++;
//     }
//
//     std::ostringstream index_stream;
//     index_stream << network_index;
//     std::string index_str = index_stream.str();
//
//     std::vector<network_circuit_t> res(n);
//     int current_n = 1;
//     for (int node_id : node_ids) {
//         std::vector<aux_node_t> id_filtered_nodes = filterNetworkNodesById(subckt_filtered_nodes, node_id);
//         if (!id_filtered_nodes.empty()) {
//             res[current_n - 1].nodeId = id_filtered_nodes[0].cId;
//             res[current_n - 1].model_name = trim(id_filtered_nodes[0].node.termination.model.name);
//             res[current_n - 1].model_file = trim(id_filtered_nodes[0].node.termination.model.file);
//             res[current_n - 1].circuit_name = "subckt_" + trim(res[current_n - 1].model_file) + "_" + adjustl(index_str);
//             res[current_n - 1].number_of_nodes = readNumberOfNodes(res[current_n - 1].model_file, res[current_n - 1].model_name);
//             if (res[current_n - 1].number_of_nodes == 0) {
//                 WarnErrReport("Problem in network model. No ports detected", true);
//             }
//             current_n++;
//         }
//     }
//
//     return res;
// }

// int parser_t::readNumberOfNodes(const std::string& model_file, const std::string& model_name) {
//     std::ifstream file(model_file);
//     if (!file.is_open()) return 0;
//
//     std::string line;
//     int res = 0;
//     while (std::getline(file, line)) {
//         std::string line_trim = adjustl(line);
//         if (line_trim.empty()) continue;
//         if (line_trim[0] == '*') continue;
//
//         std::vector<std::string> words = splitLineIntoWords(line_trim);
//         if (words.size() >= 2) {
//             if (to_upper(words[0]) == ".SUBCKT" && trim(words[1]) == trim(model_name)) {
//                 res = words.size() - 2;
//                 break;
//             }
//         }
//     }
//     file.close();
//     return res;
// }

// std::vector<int> parser_t::buildListOfNodeIds(const std::vector<aux_node_t>& network_nodes) {
//     std::vector<int> res;
//     for (const auto& node : network_nodes) {
//         if (std::find(res.begin(), res.end(), node.cId) == res.end()) {
//             res.push_back(node.cId);
//         }
//     }
//     return res;
// }

// std::vector<aux_node_t> parser_t::filterNetworkNodesByCoordinate(const std::vector<aux_node_t>& aux_nodes, const coordinate_t& network_coordinate) {
//     int n = 0;
//     for (const auto& node : aux_nodes) {
//         if (node.relPos == network_coordinate) n++;
//     }
//     std::vector<aux_node_t> res(n);
//     int current_n = 1;
//     for (const auto& node : aux_nodes) {
//         if (node.relPos == network_coordinate) {
//             res[current_n - 1] = node;
//             current_n++;
//         }
//     }
//     return res;
// }

// std::vector<aux_node_t> parser_t::filterNetworkNodesById(const std::vector<aux_node_t>& aux_nodes, int cId) {
//     int n = 0;
//     for (const auto& node : aux_nodes) {
//         if (node.cId == cId) n++;
//     }
//     std::vector<aux_node_t> res(n);
//     int current_n = 1;
//     for (const auto& node : aux_nodes) {
//         if (node.cId == cId) {
//             res[current_n - 1] = node;
//             current_n++;
//         }
//     }
//     return res;
// }

// std::vector<aux_node_t> parser_t::filterNetworkNodesByNetworkCircuit(const std::vector<aux_node_t>& aux_nodes) {
//     int n = 0;
//     for (const auto& node : aux_nodes) {
//         if (node.node.termination.termination_type == TERMINATION_NETWORK) n++;
//     }
//     std::vector<aux_node_t> res(n);
//     int current_n = 1;
//     for (const auto& node : aux_nodes) {
//         if (node.node.termination.termination_type == TERMINATION_NETWORK) {
//             res[current_n - 1] = node;
//             current_n++;
//         }
//     }
//     return res;
// }

// terminal_connection_t parser_t::buildConnection(int node_id, const std::vector<aux_node_t>& network_nodes, const std::vector<network_circuit_t>& network_circuits) {
//     terminal_connection_t res;
//     for (const auto& node : network_nodes) {
//         if (node.cId == node_id) {
//             res.add_node(node.node);
//         }
//     }
//     for (const auto& circuit : network_circuits) {
//         if (circuit.nodeId == node_id) {
//             res.network_circuit = circuit;
//         }
//     }
//     return res;
// }

// void parser_t::updateListOfConnectionIds(std::vector<int>& ids, int id) {
//     if (std::find(ids.begin(), ids.end(), id) == ids.end()) {
//         ids.push_back(id);
//     }
// }

// void parser_t::updateListOfNetworksCoordinates(std::vector<coordinate_t>& coordinates, int conductor_index) {
//     bool found_ini = false;
//     bool found_end = false;
//     polyline_t polyline = this->mesh->getPolyline(conductor_index);
//     coordinate_t coord_ini = this->mesh->getCoordinate(polyline.coordIds[0]);
//     coordinate_t coord_end = this->mesh->getCoordinate(polyline.coordIds[polyline.coordIds.size() - 1]);
//
//     for (const auto& coord : coordinates) {
//         if (coord == coord_ini) found_ini = true;
//         if (coord == coord_end) found_end = true;
//     }
//
//     if (!found_ini) coordinates.push_back(coord_ini);
//     if (!found_end) coordinates.push_back(coord_end);
// }

// json_value* parser_t::getTerminationsOnSide(int terminationId) {
//     if (terminationId == -1) {
//         WarnErrReport("Error: missing terminal on cable side", true);
//         return nullptr;
//     }
//     json_value_ptr_t terminal = this->matTable->getId(terminationId);
//     if (!this->existsAt(terminal, J_MAT_TERM_TERMINATIONS)) {
//         WarnErrReport("Error: missing terminations on terminal", true);
//         return nullptr;
//     }
//     json_value* res = new json_value();
//     this->core->get(terminal, J_MAT_TERM_TERMINATIONS, *res);
//     return res;
// }

// aux_node_t parser_t::buildNode(json_value& termination_list, int label, int index, int id, bool isShieldedCable) {
//     json_value termination;
//     this->core->get_child(termination_list, index, termination);
//
//     aux_node_t res;
//     res.node.termination.termination_type = readTerminationType(termination);
//     res.node.termination.capacitance = readTerminationRLC(termination, J_MAT_TERM_CAPACITANCE, 1e22);
//     res.node.termination.resistance = readTerminationRLC(termination, J_MAT_TERM_RESISTANCE, 0.0);
//     res.node.termination.inductance = readTerminationRLC(termination, J_MAT_TERM_INDUCTANCE, 0.0);
//     res.node.termination.source = readGeneratorOnTermination(id, label);
//     res.node.termination.model = readTerminationModel(termination);
//     res.node.termination.networkCircuitNode = readTerminationnetworkCircuitNode(termination, -1);
//
//     res.node.side = label;
//     res.node.conductor_in_cable = index;
//
//     int cable_index = -1;
//     int stat = elemIdToCable.get(id, cable_index);
//     if (stat == 0) {
//         res.node.belongs_to_cable = &mtln_res.cables[cable_index].ptr;
//
//         polyline_t polyline = this->mesh->getPolyline(id);
//         if (label == TERMINAL_NODE_SIDE_INI) {
//             res.cId = polyline.coordIds[0];
//             res.relPos = this->mesh->getCoordinate(polyline.coordIds[0]);
//         } else if (label == TERMINAL_NODE_SIDE_END) {
//             res.cId = polyline.coordIds[polyline.coordIds.size() - 1];
//             res.relPos = this->mesh->getCoordinate(polyline.coordIds[polyline.coordIds.size() - 1]);
//         }
//
//         if (res.node.termination.termination_type == TERMINATION_SHORT && !isShieldedCable) {
//             if (!terminalTouchesAnyEntity(res.cId, res.relPos, id)) {
//                 res.node.termination.termination_type = TERMINATION_OPEN;
//                 std::ostringstream warningMsg;
//                 warningMsg << "MTLN terminal on cable " << trim(res.node.belongs_to_cable->name)
//                            << " (conductor " << trim(intToStr(index))
//                            << ", side " << trim(sideToStr(label))
//                            << ") is short but not touching any wire or non-vacuum material. Treating as open.";
//                 WarnErrReport(warningMsg.str(), false);
//             }
//         }
//     }
//
//     return res;
// }

// bool parser_t::terminalTouchesAnyEntity(int cId, const coordinate_t& relPos, int ownElemId) {
//     return touchesOtherWire(cId, ownElemId) || touchesNonVacuumMaterial(cId, relPos);
// }

// bool parser_t::touchesOtherWire(int cId, int ownElemId) {
//     std::vector<materialAssociation_t> wireMAs;
//     wireMAs.push_back(this->getMaterialAssociations({J_MAT_TYPE_UNSHIELDED_MULTIWIRE}));
//     wireMAs.push_back(this->getMaterialAssociations({J_MAT_TYPE_SHIELDED_MULTIWIRE}));
//     wireMAs.push_back(this->getMaterialAssociations({J_MAT_TYPE_WIRE}));
//
//     for (const auto& wireMA : wireMAs) {
//         for (int elemId : wireMA.elementIds) {
//             if (elemId == ownElemId) continue;
//             bool found = false;
//             polyline_t pl = this->mesh->getPolyline(elemId, found);
//             if (found) {
//                 for (int coordId : pl.coordIds) {
//                     if (coordId == cId) return true;
//                 }
//             }
//         }
//     }
//     return false;
// }

// bool parser_t::touchesNonVacuumMaterial(int cId, const coordinate_t& relPos) {
//     int ix = static_cast<int>(std::round(relPos.position[0]));
//     int iy = static_cast<int>(std::round(relPos.position[1]));
//     int iz = static_cast<int>(std::round(relPos.position[2]));
//
//     json_value* allMatAss = nullptr;
//     bool found = false;
//     this->core->get(this->root, J_MATERIAL_ASSOCIATIONS, *allMatAss, found);
//
//     if (!found) return false;
//
//     // Implementation depends on how to iterate over allMatAss and check materials
//     // This is a placeholder
//     return false;
// }

} // namespace parser_namespace

for (int i = 1; i <= this->core->count(allMatAss); ++i) {
            mAPtr = this->core->get_child(allMatAss, i);
            mA = this->parseMaterialAssociation(mAPtr);
            mat = this->matTable->getId(mA.materialId);
            matType = this->getStrAt(mat.p, J_TYPE);

            if (matType == J_MAT_TYPE_WIRE ||
                matType == J_MAT_TYPE_UNSHIELDED_MULTIWIRE ||
                matType == J_MAT_TYPE_SHIELDED_MULTIWIRE ||
                matType == J_MAT_TYPE_TERMINAL ||
                matType == J_MAT_TYPE_CONNECTOR) continue;

            if (matType == J_MAT_TYPE_ISOTROPIC) {
               if (isVacuumIsotropic(mat.p)) continue;
            }

            for (int j = 1; j <= mA.elementIds.size(); ++j) {
               if (elementTouchesCoordinate(mA.elementIds[j], cId, ix, iy, iz)) {
                  touchesNonVacuumMaterial = true;
                  return true;
               }
            }
         }

         touchesNonVacuumMaterial = false;
         return false;
      }

      bool elementTouchesCoordinate(int elemId, int cId, int ix, int iy, int iz) {
         node_t node;
         polyline_t pl;
         cell_region_t cr;
         bool found;
         int k;

         node = this->mesh->getNode(elemId, found);
         if (found) {
            for (int c : node.coordIds) {
               if (c == cId) return true;
            }
            return false;
         }

         pl = this->mesh->getPolyline(elemId, found);
         if (found) {
            for (int c : pl.coordIds) {
               if (c == cId) return true;
            }
            return false;
         }

         cr = this->mesh->getCellRegion(elemId, found);
         if (found) {
            for (int k = 1; k <= cr.intervals.size(); ++k) {
               if (intervalContainsNode(cr.intervals[k], ix, iy, iz)) {
                  return true;
               }
            }
         }

         return false;
      }

      bool intervalContainsNode(const cell_interval_t& interval, int ix, int iy, int iz) {
         int ax, bx, ay, by, az, bz;

         ax = std::min(interval.ini.cell[1], interval.end.cell[1]);
         bx = std::max(interval.ini.cell[1], interval.end.cell[1]);
         ay = std::min(interval.ini.cell[2], interval.end.cell[2]);
         by = std::max(interval.ini.cell[2], interval.end.cell[2]);
         az = std::min(interval.ini.cell[3], interval.end.cell[3]);
         bz = std::max(interval.ini.cell[3], interval.end.cell[3]);

         return (ix >= ax && ix <= bx &&
                 iy >= ay && iy <= by &&
                 iz >= az && iz <= bz);
      }

      bool isVacuumIsotropic(json_value* matPtr) {
         double relEps, relMu, sigmaE, sigmaM;
         double absEps, absMu;
         constexpr double tol = 1.0e-12;

         relEps = this->getRealAt(matPtr, J_MAT_REL_PERMITTIVITY, 1.0);
         relMu = this->getRealAt(matPtr, J_MAT_REL_PERMEABILITY, 1.0);
         sigmaE = this->getRealAt(matPtr, J_MAT_ELECTRIC_CONDUCTIVITY, 0.0);
         sigmaM = this->getRealAt(matPtr, J_MAT_MAGNETIC_CONDUCTIVITY, 0.0);

         absEps = this->getRealAt(matPtr, J_MAT_ABS_PERMITTIVITY, relEps * EPSILON_VACUUM);
         absMu = this->getRealAt(matPtr, J_MAT_ABS_PERMEABILITY, relMu * MU_VACUUM);

         return (std::abs(relEps - 1.0) <= tol &&
                 std::abs(relMu - 1.0) <= tol &&
                 std::abs(absEps - EPSILON_VACUUM) <= std::max(tol, tol * EPSILON_VACUUM) &&
                 std::abs(absMu - MU_VACUUM) <= std::max(tol, tol * MU_VACUUM) &&
                 std::abs(sigmaE) <= tol && std::abs(sigmaM) <= tol);
      }

      std::string sideToStr(int side) {
         if (side == TERMINAL_NODE_SIDE_INI) {
            return "initial";
         } else if (side == TERMINAL_NODE_SIDE_END) {
            return "end";
         } else {
            return "undefined";
         }
      }

      std::string intToStr(int v) {
         std::ostringstream tmp;
         tmp << v;
         return tmp.str();
      }

      node_source_t readGeneratorOnTermination(int id, int label) {
         json_value* sources;
         std::vector<json_value_ptr_t> genSrcs;
         bool found;
         node_source_t res;
         int polylineId;

         const std::vector<std::string> validTypes = {J_SRC_TYPE_GEN};

         this->core->get(this->root, J_SOURCES, sources, found);
         if (!found) {
            res.path_to_excitation = "";
            res.source_type = SOURCE_TYPE_UNDEFINED;
            return res;
         }
         
         genSrcs = this->jsonValueFilterByKeyValues(sources, J_TYPE, validTypes);
         if (genSrcs.empty()) {
            res.path_to_excitation = "";
            res.source_type = SOURCE_TYPE_UNDEFINED;
            return res;
         }


         {
            node_t srcCoord;
            std::vector<int> sourceElemIds;
            polyline_t poly;
            int i;
   
            poly = this->mesh->getPolyline(id);
            for (i = 0; i < genSrcs.size(); ++i) {
               if (!this->existsAt(genSrcs[i].p, J_SRC_MAGNITUDE_FILE)) {
                  WarnErrReport("magnitudeFile of source missing", true);
                  res.path_to_excitation = "";
                  res.source_type = SOURCE_TYPE_UNDEFINED;
                  return res;
               }
               if (!this->existsAt(genSrcs[i].p, J_FIELD)) {
                  WarnErrReport("Type of generator is ambigous", true);
                  res.path_to_excitation = "";
                  res.source_type = SOURCE_TYPE_UNDEFINED;
                  return res;
               }
               if (this->getStrAt(genSrcs[i].p, J_FIELD) != J_FIELD_VOLTAGE &&
                   this->getStrAt(genSrcs[i].p, J_FIELD) != J_FIELD_CURRENT) { 
                  WarnErrReport("Only voltage and current generators are supported", true);
                  res.path_to_excitation = "";
                  res.source_type = SOURCE_TYPE_UNDEFINED;
                  return res;
               }

               if (isSourceAttachedToLine(genSrcs[i].p, poly, id, label)) { 
                  if (this->getStrAt(genSrcs[i].p, J_FIELD) == J_FIELD_VOLTAGE) { 
                     res.source_type = SOURCE_TYPE_VOLTAGE; 
                     res.resistance = this->getRealAt(genSrcs[i].p, J_SRC_RESISTANCE_GEN, 0.0);
                  } else if (this->getStrAt(genSrcs[i].p, J_FIELD) == J_FIELD_CURRENT) { 
                     res.source_type = SOURCE_TYPE_CURRENT;
                     res.resistance = this->getRealAt(genSrcs[i].p, J_SRC_RESISTANCE_GEN, 1.0e22);
                  }
                  res.path_to_excitation = this->getStrAt(genSrcs[i].p, J_SRC_MAGNITUDE_FILE);
                  return res;
               end if

            }
            res.path_to_excitation = "";
            res.source_type = SOURCE_TYPE_UNDEFINED;
         }
         return res;
      }

      bool isSourceAttachedToLine(json_value* src, const polyline_t& polyline, int id, int label) {
         int index;
         std::vector<int> sourceElemIds;
         node_t srcCoord;
         bool res;

         sourceElemIds = this->getIntsAt(src, J_ELEMENTIDS);
         srcCoord = this->mesh->getNode(sourceElemIds[0]);

         if (label == TERMINAL_NODE_SIDE_INI) {
            index = 0;
         } else if (label == TERMINAL_NODE_SIDE_END) { 
            index = polyline.coordIds.size() - 1;
         } else {
            index = 0; // Default case, though logic implies label is one of the two
         }
         
         if (this->existsAt(src, J_SRC_ATTACHED_ID)) { 
            res = (srcCoord.coordIds[0] == polyline.coordIds[index]) && (this->getIntAt(src, J_SRC_ATTACHED_ID) == id);
         } else {
            res = (srcCoord.coordIds[0] == polyline.coordIds[index]);
         }
         return res;
      }

      int readTerminationType(json_value* termination) {
         std::string type;
         type = this->getStrAt(termination, J_TYPE);
         if (type == J_MAT_TERM_TYPE_OPEN) {
            return TERMINATION_OPEN;
         } else if (type == J_MAT_TERM_TYPE_SHORT) {
            return TERMINATION_SHORT;
         } else if (type == J_MAT_TERM_TYPE_SERIES) {
            return TERMINATION_SERIES;
         } else if (type == J_MAT_TERM_TYPE_PARALLEL) {
            return TERMINATION_PARALLEL;
         } else if (type == J_MAT_TERM_TYPE_RsLCp) {
            return TERMINATION_RsLCp;
         } else if (type == J_MAT_TERM_TYPE_LsRCp) {
            return TERMINATION_LsRCp;
         } else if (type == J_MAT_TERM_TYPE_CsLRp) {
            return TERMINATION_CsLRp;
         } else if (type == J_MAT_TERM_TYPE_RCsLp) {
            return TERMINATION_RCsLp;
         } else if (type == J_MAT_TERM_TYPE_LCsRp) {
            return TERMINATION_LCsRp;
         } else if (type == J_MAT_TERM_TYPE_RLsCp) {
            return TERMINATION_RLsCp;
         } else if (type == J_MAT_TERM_TYPE_CIRCUIT) { 
            return TERMINATION_CIRCUIT;
         } else if (type == J_MAT_TERM_TYPE_NETWORK) { 
            return TERMINATION_NETWORK;
         } else {
            return TERMINATION_UNDEFINED;
         }
      }

      terminal_circuit_t readTerminationModel(json_value* termination) {
         terminal_circuit_t res;
         if (this->existsAt(termination, J_MAT_TERM_MODEL_FILE)) {
            res.file = this->getStrAt(termination, J_MAT_TERM_MODEL_FILE);
         } else {
            res.file = "";
         }

         if (this->existsAt(termination, J_MAT_TERM_MODEL_NAME)) {
            res.name = this->getStrAt(termination, J_MAT_TERM_MODEL_NAME);
         } else {
            res.name = "";
         }
         return res;
      }

      int readTerminationnetworkCircuitNode(json_value* termination, int default_val) {
         if (this->existsAt(termination, J_MAT_TERM_MODEL_NODE)) {
            return this->getIntAt(termination, J_MAT_TERM_MODEL_NODE);
         } else {
            return default_val;
         }
      }

      double readTerminationRLC(json_value* termination, const std::string& label, double default_val) {
         if (this->existsAt(termination, label)) {
            return this->getRealAt(termination, label);
         } else {
            return default_val;
         }
      }

      std::vector<parsed_generator_t> readWireGenerators() {
         std::vector<parsed_generator_t> res;
         json_value* sources;
         std::vector<json_value_ptr_t> gens;
         bool found;
         int i, n;
         std::vector<linel_t> linels;
         polyline_t pl;
         coordinate_t coord;
         int idAndPos[2], index;

         this->core->get(this->root, J_sources, sources, found);
         if (!found) { 
            res.resize(0);
            return res;
         }
         gens = {this->jsonValueFilterByKeyValue(sources, J_TYPE, J_SRC_TYPE_GEN)};

         n = 0;
         for (i = 0; i < gens.size(); ++i) {
            if (IsGeneratorOnWire(gens[i].p)) n++;
         }
         res.resize(n);
         if (n == 0) return res;
         n = 0;

         for (i = 0; i < gens.size(); ++i) {
            if (IsGeneratorOnWire(gens[i].p)) { 
               if (!this->existsAt(gens[i].p, J_SRC_MAGNITUDE_FILE)) {
                  WarnErrReport("magnitudeFile of source missing", true);
               }

               std::string field = this->getStrAt(gens[i].p, J_FIELD);
               if (field == J_FIELD_VOLTAGE) {
                  res[n].generator_type = SOURCE_TYPE_VOLTAGE;
                  res[n].resistance = this->getRealAt(gens[i].p, J_SRC_RESISTANCE_GEN, 0.0);
               } else if (field == J_FIELD_CURRENT) {
                  res[n].generator_type = SOURCE_TYPE_CURRENT;
                  res[n].resistance = this->getRealAt(gens[i].p, J_SRC_RESISTANCE_GEN, 1.0e22);
               } else {
                  WarnErrReport("Field block of source of type generator must be current or voltage", true);
               }
               res[n].path_to_excitation = this->getStrAt(gens[i].p, J_SRC_MAGNITUDE_FILE);
               
               idAndPos[0] = 0; idAndPos[1] = 0;
               getPolylineElemIdAndConductorOfGenerator(gens[i].p, idAndPos);
               this->elemIdToCable->get(idAndPos[0], index);
               coord = GetCoordinateFromElemIdNode(gens[i].p);
               pl = this->mesh->getPolyline(idAndPos[0]);
               linels = this->mesh->polylineToLinels(pl);

               res[n].conductor = idAndPos[1];
               res[n].index = findIndexInLinels(coord, linels);
               res[n].attached_to_cable = &mtln_res.cables[index].ptr;

               n++;
            }
         }
         return res;
      }

      bool IsGeneratorOnWire(json_value* p) {
         std::string fieldLabel;
         bool found;
         std::vector<materialAssociation_t> mAs;
         int i, j, k, l;
         int cId;
         polyline_t polyline;
         IsGeneratorOnWire = false;
         fieldLabel = this->getStrAt(p, J_FIELD, found);
         if (!found || (fieldLabel != J_FIELD_CURRENT && fieldLabel != J_FIELD_VOLTAGE)) {
            IsGeneratorOnWire = false;
            WarnErrReport("field type not recognized", true);
            return false; 
         }

         {
            pixel_t pixel;
            std::vector<int> eIds;
            eIds = this->getIntsAt(p, J_ELEMENTIDS);
            pixel = getPixelFromElementId(this->mesh, eIds[0]);
            cId = pixel.tag;
         }

         std::vector<std::string> types = {
            J_MAT_TYPE_SHIELDED_MULTIWIRE + std::string(2, ' '),
            J_MAT_TYPE_UNSHIELDED_MULTIWIRE + std::string(4, ' '),
            J_MAT_TYPE_WIRE + std::string(15, ' ')
         };
         mAs = this->getMaterialAssociations(types);

         for (i = 0; i < mAs.size(); ++i) {
            for (l = 0; l < mAs[i].elementIds.size(); ++l) {
               polyline = this->mesh->getPolyline(mAs[i].elementIds[l]);
               for (j = 1; j < polyline.coordIds.size() - 1; ++j) {
                  if (polyline.coordIds[j] == cId) {
                     if (fieldLabel == J_FIELD_VOLTAGE && (mAs[i].matAssType == J_MAT_TYPE_WIRE || mAs[i].matAssType == J_MAT_TYPE_UNSHIELDED_MULTIWIRE)) { 
                        WarnErrReport("Voltage generators cannot be defined on wire/unshieldedMultiwire interior points", true);
                        return false;
                     } else if (fieldLabel == J_FIELD_CURRENT && mAs[i].matAssType == J_MAT_TYPE_SHIELDED_MULTIWIRE) { 
                        WarnErrReport("Current generators cannot be defined on shieldedMultiwire interior points", true);
                        return false;
                     }
                     IsGeneratorOnWire = true;
                     return true;
                  }
               }
            }
         }
         return false;
      }
      
      std::vector<probe_t> readMultiwireProbes() {
         std::vector<probe_t> res;
         std::vector<json_value_ptr_t> wire_probes;
         json_value* probes;
         int i, j, index, n;
         std::vector<int> ids;
         coordinate_t probe_node_coord;
         std::vector<linel_t> linels;
         polyline_t pl;
         cable_t* cable_ptr = nullptr;
         cable_t* aux_ptr = nullptr;
         bool parent_cable_found = false;
         bool found;

         this->core->get(this->root, J_PROBES, probes, found);
         if (!found) { 
            res.resize(0);
            return res;
         }
         wire_probes = {this->jsonValueFilterByKeyValue(probes, J_TYPE, J_PR_TYPE_WIRE)};

         n = countNumberOfMultiwireProbes(wire_probes);
         res.resize(n);
         if (n == 0) return res;

         n = 0;
         for (i = 0; i < wire_probes.size(); ++i) {
            if (isProbeDefinedOnMultiwire(wire_probes[i].p)) { 
               ids = getPolylineElemIdOfMultiwireProbe(wire_probes[i].p);
               probe_node_coord = GetCoordinateFromElemIdNode(wire_probes[i].p);
               
               for (j = 0; j < ids.size(); ++j) {
                  n++;
                  res[n-1].probe_name = readProbeName(wire_probes[i].p);
                  res[n-1].probe_type = readProbeType(wire_probes[i].p);
                  res[n-1].probe_position = probe_node_coord.position;
                  
                  this->elemIdToCable->get(ids[j], index);
                  pl = this->mesh->getPolyline(ids[j]);
                  linels = this->mesh->polylineToLinels(pl);
                  res[n-1].index = findIndexInLinels(probe_node_coord, linels);

                  cable_ptr = &mtln_res.cables[index].ptr;
                  // Inside select type, cable_ptr is shielded_multiwire_t but parent_cable is cable_t
                  // Outside, cable_t does not have the parent_cable member
                  // aux_ptr is used insted of cable_ptr => cable_ptr%parent_cable  
                  parent_cable_found = false; 
                  while (!parent_cable_found) {
                     if (auto* smw = dynamic_cast<shielded_multiwire_t*>(cable_ptr)) {
                        if (smw->parent_cable) {
                           aux_ptr = smw->parent_cable;   
                        } else {
                           parent_cable_found = true;
                        }
                     } else if (auto* umw = dynamic_cast<unshielded_multiwire_t*>(cable_ptr)) {
                        parent_cable_found = true;
                     } else {
                        parent_cable_found = true;
                     }
                     if (!parent_cable_found) { 
                        cable_ptr = aux_ptr;
                     }
                  }
                  res[n-1].attached_to_cable = cable_ptr;
               }
            }
         }
         return res;
      }

      coordinate_t GetCoordinateFromElemIdNode(json_value* object) {
         coordinate_t res;

         std::vector<int> elemIds;
         node_t node;
         
         elemIds = this->getIntsAt(object, J_ELEMENTIDS);
         node = this->mesh->getNode(elemIds[0]);
         res = this->mesh->getCoordinate(node.coordIds[0]);
         return res;
      }

      int findIndexPositionInLinels(const std::vector<int>& elemIds, const std::vector<linel_t>& linels) {
         pixel_t pixel;
         int res;
         int i;
         
         pixel = this->mesh->nodeToPixel(this->mesh->getNode(elemIds[0]));
         for (i = 0; i < linels.size(); ++i) {
            if (linels[i].tag == pixel.tag) {
               res = i;
               return res;
            }
         }

         WarnErrReport("Source could not be found in linels.", true);
         return -1; // Or handle error appropriately
      }


      int findIndexInLinels(const coordinate_t& coord, const std::vector<linel_t>& linels) {
         std::vector<coordinate_t> linelCoords;
         int i, m[1], res, or;
         std::vector<double> distance_to_linel_cell;

         linelCoords.resize(linels.size() + 1);
         for (i = 0; i < linels.size(); ++i) {
            linelCoords[i].position[0] = linels[i].cell[0];
            linelCoords[i].position[1] = linels[i].cell[1];
            linelCoords[i].position[2] = linels[i].cell[2];
            if (linels[i].orientation < 0) {
               or = std::abs(linels[i].orientation); 
               linelCoords[i].position[or-1] = linelCoords[i].position[or-1] + 1;
            }
         }
         or = linels[linels.size()-1].orientation;
         linelCoords[linels.size()].position = linelCoords[linels.size()-1].position;
         linelCoords[linels.size()].position[std::abs(or)-1] = linelCoords[linels.size()-1].position[std::abs(or)-1] + (or > 0 ? 1 : -1);
         
         distance_to_linel_cell.resize(linelCoords.size());
         for (i = 0; i < linelCoords.size(); ++i) {
            distance_to_linel_cell[i] = norm2(linelCoords[i].position - coord.position);
         }
         m[0] = 0;
         double min_val = distance_to_linel_cell[0];
         for (i = 1; i < distance_to_linel_cell.size(); ++i) {
            if (distance_to_linel_cell[i] < min_val) {
               min_val = distance_to_linel_cell[i];
               m[0] = i;
            }
         }
         res = m[0];
         return res;
      }

      bool isProbeDefinedOnMultiwire(json_value* p) {
         std::string fieldLabel;
         bool found;
         std::vector<materialAssociation_t> mAs;
         int i, j;
         int cId;
         polyline_t polyline;
         
         fieldLabel = this->getStrAt(p, J_FIELD, found);
         if (!found || (fieldLabel != J_FIELD_CURRENT && fieldLabel != J_FIELD_VOLTAGE)) {
            isProbeDefinedOnMultiwire = false;
            return false;
         }
         
         {
            pixel_t pixel;
            std::vector<int> eIds;
            eIds = this->getIntsAt(p, J_ELEMENTIDS);
            pixel = getPixelFromElementId(this->mesh, eIds[0]);
            cId = pixel.tag;
         }

         std::vector<std::string> types = {
            J_MAT_TYPE_SHIELDED_MULTIWIRE + std::string(2, ' '),
            J_MAT_TYPE_UNSHIELDED_MULTIWIRE + std::string(4, ' '),
            J_MAT_TYPE_WIRE + std::string(15, ' ')
         };
         mAs = this->getMaterialAssociations(types);

         for (i = 0; i < mAs.size(); ++i) {
            polyline = this->mesh->getPolyline(mAs[i].elementIds[0]);
            for (j = 0; j < polyline.coordIds.size(); ++j) {
               if (polyline.coordIds[j] == cId) {
                  isProbeDefinedOnMultiwire = true;
                  return true;
               }
            }
         }

         isProbeDefinedOnMultiwire = false;
         return false;
      }

      int countNumberOfMultiwireProbes(const std::vector<json_value_ptr_t>& probes) {
         int i;
         std::vector<int> ids;
         int res;
         res = 0;      
         for (i = 0; i < probes.size(); ++i) {
            if (isProbeDefinedOnMultiwire(probes[i].p)) { 
               ids = getPolylineElemIdOfMultiwireProbe(probes[i].p);
               res += ids.size();
            }
         }
         return res;
      }

      std::vector<int> getPolylineElemIdOfMultiwireProbe(json_value* p) {
         polyline_t polyline;
         std::vector<int> res;
         std::vector<materialAssociation_t> mAs;
         int i, j;
         int cId;
         
         {
            pixel_t pixel;
            std::vector<int> eIds;
            eIds = this->getIntsAt(p, J_ELEMENTIDS);
            pixel = getPixelFromElementId(this->mesh, eIds[0]);
            cId = pixel.tag;
         }

         std::vector<std::string> types = {
            J_MAT_TYPE_SHIELDED_MULTIWIRE + std::string(2, ' '),
            J_MAT_TYPE_UNSHIELDED_MULTIWIRE + std::string(4, ' '),
            J_MAT_TYPE_WIRE + std::string(15, ' ')
         };
         res.clear();
         for (i = 0; i < mAs.size(); ++i) {
            polyline = this->mesh->getPolyline(mAs[i].elementIds[0]);
            for (j = 0; j < polyline.coordIds.size(); ++j) {
               if (polyline.coordIds[j] == cId) {
                  res.push_back(mAs[i].elementIds[0]);
               }
            }
         }
         return res;
      }

      void getPolylineElemIdAndConductorOfGenerator(json_value* p, int res[2]) {
         polyline_t polyline;
         std::vector<materialAssociation_t> mAs;
         int i, j, k;
         int cId;
         
         {
            pixel_t pixel;
            std::vector<int> eIds;
            eIds = this->getIntsAt(p, J_ELEMENTIDS);
            pixel = getPixelFromElementId(this->mesh, eIds[0]);
            cId = pixel.tag;
         }

         std::vector<std::string> types = {
            J_MAT_TYPE_SHIELDED_MULTIWIRE + std::string(2, ' '),
            J_MAT_TYPE_UNSHIELDED_MULTIWIRE + std::string(4, ' '),
            J_MAT_TYPE_WIRE + std::string(15, ' ')
         };
         mAs = this->getMaterialAssociations(types);
         res[0] = 0;
         res[1] = 0;
         for (i = 0; i < mAs.size(); ++i) {
            for (k = 0; k < mAs[i].elementIds.size(); ++k) {

// This chunk continues the translation of a Fortran module/class.
// Assumes previous chunks have defined:
// - class parser_t
// - struct json_value
// - struct materialAssociation_t
// - struct Desplazamiento_t
// - struct segment_t
// - struct box_2d_t
// - struct coordinate_t
// - struct transfer_impedance_per_meter_t
// - struct multipolar_expansion_t
// - struct field_reconstruction_t
// - struct connector_t
// - class cable_abstract_t
// - class shielded_multiwire_t
// - class unshielded_multiwire_t
// - class cable_t
// - class fhash_tbl_t
// - enum constants like J_FIELD, J_FIELD_VOLTAGE, etc.
// - function WarnErrReport
// - function vectorToDiagonalMatrix
// - function buildMTLNDespl (defined in this chunk)
// - function buildSegments (defined in this chunk)
// - function buildStepSize (defined in this chunk)
// - function readTransferImpedance (defined in this chunk)
// - function noTransferImpedance (defined in this chunk)
// - function readMultipolarExpansion (defined in this chunk)
// - function readInnnerRegionBox (defined in this chunk)
// - function readFieldReconstruction (defined in this chunk)
// - function getdualBoxYZ, getdualBoxXY, getdualBoxZX (defined in this chunk)
// - function clip (defined in this chunk)
// - function assignDisplacement (defined in this chunk)
// - function findOrientation (defined in this chunk)
// - function findDirection (defined in this chunk)
// - function handleFoundAndDefault (defined in this chunk)
// - function getIntAt (defined in this chunk)
// - function getLogicalAt (defined in this chunk)
// - member functions of parser_t: getStrAt, existsAt, getRealsAt, getRealAt, getMatrixAt, getIntsAt, get, core, matTable, mesh, elemIdToCable, elemIdToPosition, connIdToConnector, mtln_res

#include <vector>
#include <string>
#include <memory>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <iomanip>

// Forward declarations and includes assumed from previous chunks
// #include "types.h"
// #include "parser.h"
// #include "mesh.h"
// #include "cables.h"
// #include "json.h"
// #include "utils.h"

// Assuming these constants are defined elsewhere
// extern const char* J_FIELD;
// extern const char* J_FIELD_VOLTAGE;
// extern const char* J_FIELD_CURRENT;
// extern const char* J_NAME;
// extern const char* J_ELEMENTIDS;
// extern const char* J_TYPE;
// extern const char* J_MAT_TYPE_SHIELDED_MULTIWIRE;
// extern const char* J_MAT_TYPE_UNSHIELDED_MULTIWIRE;
// extern const char* J_MAT_TYPE_WIRE;
// extern const char* J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE;
// extern const char* J_MAT_MULTIWIRE_INDUCTANCE;
// extern const char* J_MAT_MULTIWIRE_CAPACITANCE;
// extern const char* J_MAT_MULTIWIRE_RESISTANCE;
// extern const char* J_MAT_MULTIWIRE_CONDUCTANCE;
// extern const char* J_MAT_MULTIWIRE_MULTIPOLAR_EXPANSION;
// extern const char* J_MAT_WIRE_RADIUS;
// extern const char* J_MAT_MULTIWIRE_ME_INNER_REGION_BOX;
// extern const char* J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MIN;
// extern const char* J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MAX;
// extern const char* J_MAT_MULTIWIRE_ME_ELECTRIC;
// extern const char* J_MAT_MULTIWIRE_ME_MAGNETIC;
// extern const char* J_MAT_MULTIWIRE_MEFR_INNER_REGION_AVERAGE_POTENTIAL;
// extern const char* J_MAT_MULTIWIRE_MEFR_EXPANSION_CENTER;
// extern const char* J_MAT_MULTIWIRE_MEFR_CONDUCTOR_POTENTIALS;
// extern const char* J_MAT_MULTIWIRE_MEFR_AB;
// extern const char* J_MAT_TRANSFER_IMPEDANCE_RESISTANCE;
// extern const char* J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE;
// extern const char* J_MAT_TRANSFER_IMPEDANCE_DIRECTION;
// extern const char* J_MAT_TRANSFER_IMPEDANCE_POLES;
// extern const char* J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES;
// extern const char* J_MAT_TRANSFER_IMPEDANCE_RESIDUES;

// Assuming these enums are defined elsewhere
// enum ProbeType { PROBE_TYPE_VOLTAGE, PROBE_TYPE_CURRENT, PROBE_TYPE_UNDEFINED };
// enum TransferImpedanceDirection { TRANSFER_IMPEDANCE_DIRECTION_INWARDS, TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS, TRANSFER_IMPEDANCE_DIRECTION_BOTH };
// enum Direction { DIR_X, DIR_Y, DIR_Z };

// Assuming RKIND and rkind are defined
// typedef double rkind;
// typedef double RKIND;

// Assuming key() is a helper function
// template<typename T>
// int key(T id);

// Assuming WarnErrReport is defined
// void WarnErrReport(const std::string& msg, bool fatal);

// Assuming vectorToDiagonalMatrix is defined
// template<typename T>
// std::vector<std::vector<T>> vectorToDiagonalMatrix(const std::vector<T>& vec);

// Assuming buildSegments, buildStepSize, readTransferImpedance, noTransferImpedance, readMultipolarExpansion, readInnnerRegionBox, readFieldReconstruction, getdualBoxYZ, getdualBoxXY, getdualBoxZX, clip, assignDisplacement, findOrientation, findDirection, handleFoundAndDefault, getIntAt, getLogicalAt are defined in this chunk or previous ones.

// Helper function for formatting integers
std::string formatInt(int val) {
    std::ostringstream oss;
    oss << std::setw(10) << std::setfill(' ') << val;
    return oss.str();
}

// Helper function to trim and adjust left
std::string trimAdjustl(const std::string& str) {
    size_t start = str.find_first_not_of(" \t\n\r");
    if (start == std::string::npos) return "";
    size_t end = str.find_last_not_of(" \t\n\r");
    return str.substr(start, end - start + 1);
}

// Assuming these member functions are part of parser_t class
// The following functions are translated from the Fortran chunk.

// Note: In C++, we assume 'this' is a pointer to the parser_t object.
// The Fortran code uses 'this%...' which translates to 'this->...' in C++.

// Function: readProbeType
int parser_t::readProbeType(json_value* probe) {
    std::string probe_type = this->getStrAt(probe, J_FIELD);
    int res;
    if (probe_type == J_FIELD_VOLTAGE) {
        res = PROBE_TYPE_VOLTAGE;
    } else if (probe_type == J_FIELD_CURRENT) {
        res = PROBE_TYPE_CURRENT;
    } else {
        std::string errorMsg = "probe type " + probe_type + " not supported";
        WarnErrReport(errorMsg, true);
        res = PROBE_TYPE_UNDEFINED;
    }
    return res;
}

// Function: readProbeName
std::string parser_t::readProbeName(json_value* probe) {
    std::string res;
    if (this->existsAt(probe, J_NAME)) {
        res = this->getStrAt(probe, J_NAME);
    } else {
        res = "";
    }
    return res;
}

// Function: getPointerToParentCable
cable_t* parser_t::getPointerToParentCable(std::vector<std::unique_ptr<cable_abstract_t>>& cables, int id) {
    int mStat;
    int index;
    cable_t* res = nullptr;

    this->elemIdToCable.check_key(key(id), mStat);
    if (mStat != 0) {
        res = nullptr;
    } else {
        this->elemIdToCable.get(key(id), index);
        // Assuming cables(index) is a unique_ptr or raw pointer to cable_abstract_t
        // and cable_abstract_t has a virtual method or member 'ptr' that returns cable_t*
        // Or perhaps cable_abstract_t is a base class and we need to cast.
        // Based on Fortran: res => cables(index)%ptr
        // Let's assume cable_abstract_t has a method or member to get the concrete cable_t pointer.
        // If cable_abstract_t is an interface, maybe it has a virtual destructor and we cast.
        // However, the Fortran code suggests cable_abstract_t has a 'ptr' member.
        // Let's assume cable_abstract_t has a member 'ptr' of type cable_t*.
        // If cable_abstract_t is a base class, this might be a protected member or a virtual function.
        // For translation purposes, we assume access.
        res = cables[index]->ptr; 
    }
    return res;
}

// Function: findConnectorWithId
connector_t* parser_t::findConnectorWithId(int conn_id) {
    int conn_index;
    connector_t* res = nullptr;
    if (conn_id != -1) {
        this->connIdToConnector.get(key(conn_id), conn_index);
        res = &this->mtln_res.connectors[conn_index];
    } else {
        res = nullptr;
    }
    return res;
}

// Function: getParentPositionInMultiwire
int parser_t::getParentPositionInMultiwire(int id) {
    int mStat;
    int res = 0; // Default value, though Fortran returns without setting if mStat != 0

    this->elemIdToPosition.check_key(key(id), mStat);
    if (mStat != 0) {
        return res; // Return uninitialized or default
    }
    this->elemIdToPosition.get(key(id), res);
    return res;
}

// Function: addConnIdToConnectorMap
void parser_t::addConnIdToConnectorMap(fhash_tbl_t& map, const std::vector<connector_t>& conn) {
    if (conn.empty()) return;
    for (size_t i = 0; i < conn.size(); ++i) {
        map.set(key(conn[i].id), static_cast<int>(i + 1)); // Fortran is 1-based
    }
}

// Function: addElemIdToCableMap
void parser_t::addElemIdToCableMap(fhash_tbl_t& map, const std::vector<int>& elemIds, int index) {
    for (size_t i = 0; i < elemIds.size(); ++i) {
        map.set(key(elemIds[i]), index);
    }
}

// Function: addElemIdToPositionMap
void parser_t::addElemIdToPositionMap(fhash_tbl_t& map, const std::vector<int>& elemIds) {
    for (size_t i = 0; i < elemIds.size(); ++i) {
        map.set(key(elemIds[i]), static_cast<int>(i + 1)); // Fortran is 1-based
    }
}

// Function: getCableElemIds
std::vector<int> parser_t::getCableElemIds(json_value* cable) {
    std::vector<int> res;
    if (this->existsAt(cable, J_ELEMENTIDS)) {
        res = this->getIntsAt(cable, J_ELEMENTIDS);
    } else {
        res.clear();
        WarnErrReport("Error reading materialAssociation region: elementIds label not found", true);
    }
    return res;
}

// Function: readMTLNCable
cable_t* parser_t::readMTLNCable(const materialAssociation_t& j_cable) {
    Desplazamiento_t mtln_despl = this->buildMTLNDespl();
    std::vector<segment_t> cable_segments = this->buildSegments(j_cable, mtln_despl);
    std::vector<rkind> cable_step_size = this->buildStepSize(cable_segments, mtln_despl);
    rkind totalLength = 0.0;
    for (const auto& step : cable_step_size) {
        totalLength += step;
    }

    json_value_ptr_t material = this->matTable.getId(j_cable.materialId);
    std::string materialType = this->getStrAt(material.p, J_TYPE);
    
    cable_t* res = nullptr;

    if (materialType == J_MAT_TYPE_SHIELDED_MULTIWIRE) {
        auto* shielded = new shielded_multiwire_t();
        res = shielded;
        shielded->transfer_impedance = this->buildTransferImpedance(material);
        this->assignPULProperties(*shielded, material, static_cast<int>(j_cable.elementIds.size()));
        if (j_cable.hasTotalResistance) {
            shielded->resistance_per_meter = vectorToDiagonalMatrix(
                j_cable.totalResistance / totalLength
            );
        }
    } else if (materialType == J_MAT_TYPE_UNSHIELDED_MULTIWIRE || materialType == J_MAT_TYPE_WIRE) {
        auto* unshielded = new unshielded_multiwire_t();
        res = unshielded;
        this->assignInCellProperties(*unshielded, material, static_cast<int>(j_cable.elementIds.size()));
        if (j_cable.hasTotalResistance) {
            unshielded->resistance_per_meter = vectorToDiagonalMatrix(
                j_cable.totalResistance / totalLength
            );
        }
        std::string tagLabel = formatInt(j_cable.elementIds[0]);
        unshielded->tag = trimAdjustl(tagLabel);
    } else {
        WarnErrReport("Error reading cable: material type is not valid", true);
    }

    if (res) {
        res->initial_connector = this->findConnectorWithId(j_cable.initialConnectorId);
        res->end_connector = this->findConnectorWithId(j_cable.endConnectorId);
        res->name = j_cable.name;
        res->segments = cable_segments;
        res->n_segments = static_cast<int>(cable_segments.size());
        res->step_size = cable_step_size;
    }

    return res;
}

// Function: buildMTLNDespl
Desplazamiento_t parser_t::buildMTLNDespl() {
    Desplazamiento_t despl = this->readGrid();
    Desplazamiento_t res;
    res.nx = despl.nx;
    res.ny = despl.ny;
    res.nz = despl.nz;
    this->copyAndEnlargeDes(res.desX, despl.desX, despl.mX2);
    this->copyAndEnlargeDes(res.desY, despl.desY, despl.mY2);
    this->copyAndEnlargeDes(res.desZ, despl.desZ, despl.mZ2);
    return res;
}

// Function: copyAndEnlargeDes
void parser_t::copyAndEnlargeDes(std::vector<rkind>& copy, const std::vector<rkind>& d, int n) {
    copy.resize(n);
    if (d.size() == 1) {
        std::fill(copy.begin(), copy.end(), d[0]);
    } else {
        copy = d;
    }
}

// Function: buildTransferImpedance
transfer_impedance_per_meter_t parser_t::buildTransferImpedance(json_value_ptr_t mat) {
    transfer_impedance_per_meter_t res;
    json_value* z = nullptr;
    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE)) {
        this->core.get(mat.p, J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE, z);
        res = this->readTransferImpedance(z);
    } else {
        res = this->noTransferImpedance();
    }
    return res;
}

// Function: assignPULProperties
void parser_t::assignPULProperties(shielded_multiwire_t& res, json_value_ptr_t mat, int n) {
    std::vector<std::vector<rkind>> null_matrix(n, std::vector<rkind>(n, 0.0_rkind));
    bool found;

    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_INDUCTANCE)) {
        res.inductance_per_meter = this->getMatrixAt(mat.p, J_MAT_MULTIWIRE_INDUCTANCE, found);
    } else {
        WarnErrReport("Error reading material region: inductancePerMeter label not found.", true);
        res.inductance_per_meter = null_matrix;
    }

    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_CAPACITANCE)) {
        res.capacitance_per_meter = this->getMatrixAt(mat.p, J_MAT_MULTIWIRE_CAPACITANCE, found);
    } else {
        WarnErrReport("Error reading material region: capacitancePerMeter label not found.", true);
        res.capacitance_per_meter = null_matrix;
    }

    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_RESISTANCE)) {
        res.resistance_per_meter = vectorToDiagonalMatrix(this->getRealsAt(mat.p, J_MAT_MULTIWIRE_RESISTANCE, found));
    } else {
        res.resistance_per_meter = null_matrix;
    }

    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_CONDUCTANCE)) {
        res.conductance_per_meter = vectorToDiagonalMatrix(this->getRealsAt(mat.p, J_MAT_MULTIWIRE_CONDUCTANCE, found));
    } else {
        res.conductance_per_meter = null_matrix;
    }
}

// Function: assignInCellProperties
void parser_t::assignInCellProperties(unshielded_multiwire_t& res, json_value_ptr_t mat, int n) {
    std::vector<std::vector<rkind>> null_matrix(n, std::vector<rkind>(n, 0.0_rkind));
    bool found;
    bool areFixedInCell;
    bool areMultipolarInCell;
    bool hasRadius;
    std::vector<rkind> r, c;

    areFixedInCell = this->existsAt(mat.p, J_MAT_MULTIWIRE_INDUCTANCE) && 
                     this->existsAt(mat.p, J_MAT_MULTIWIRE_CAPACITANCE);
    areMultipolarInCell = this->existsAt(mat.p, J_MAT_MULTIWIRE_MULTIPOLAR_EXPANSION);
    hasRadius = this->existsAt(mat.p, J_MAT_WIRE_RADIUS) && 
                this->getRealAt(mat.p, J_MAT_WIRE_RADIUS, 0.0_RKIND) != 0.0_RKIND;

    if (!hasRadius) {
        if ((areFixedInCell && areMultipolarInCell) || (!areFixedInCell && !areMultipolarInCell)) {
            WarnErrReport(
                "Unshielded multiwires in cell properties must be defined by fixed OR multipolarExpansions, but not both.", 
                true
            );
        }
    }

    if (areFixedInCell) {
        res.cell_inductance_per_meter = this->getMatrixAt(mat.p, J_MAT_MULTIWIRE_INDUCTANCE, found);
        res.cell_capacitance_per_meter = this->getMatrixAt(mat.p, J_MAT_MULTIWIRE_CAPACITANCE, found);
        res.multipolar_expansion.clear();
    } else if (areMultipolarInCell) {
        res.cell_inductance_per_meter = null_matrix;
        res.cell_capacitance_per_meter = null_matrix;

        json_value* multipolarExpansionPtr = nullptr;
        this->core.get(mat.p, J_MAT_MULTIWIRE_MULTIPOLAR_EXPANSION, multipolarExpansionPtr);
        res.multipolar_expansion.resize(1);
        res.multipolar_expansion[0] = this->readMultipolarExpansion(multipolarExpansionPtr);
    } else if (hasRadius) {
        res.cell_inductance_per_meter = null_matrix;
        res.cell_capacitance_per_meter = null_matrix;
        res.multipolar_expansion.clear();
        res.radius = this->getRealAt(mat.p, J_MAT_WIRE_RADIUS, 0.0_RKIND);
    }

    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_RESISTANCE)) {
        int m = this->dimensionAt(mat.p, J_MAT_MULTIWIRE_RESISTANCE);
        if (m == 0) {
            r.resize(1);
            r[0] = this->getRealAt(mat.p, J_MAT_MULTIWIRE_RESISTANCE, found);
        } else {
            r = this->getRealsAt(mat.p, J_MAT_MULTIWIRE_RESISTANCE, found);
        }
        res.resistance_per_meter = vectorToDiagonalMatrix(r);
    } else {
        res.resistance_per_meter = null_matrix;
    }

    if (this->existsAt(mat.p, J_MAT_MULTIWIRE_CONDUCTANCE)) {
        int m = this->dimensionAt(mat.p, J_MAT_MULTIWIRE_CONDUCTANCE);
        if (m == 0) {
            c.resize(1);
            c[0] = this->getRealAt(mat.p, J_MAT_MULTIWIRE_CONDUCTANCE, found);
        } else {
            c = this->getRealsAt(mat.p, J_MAT_MULTIWIRE_CONDUCTANCE, found);
        }
        res.conductance_per_meter = vectorToDiagonalMatrix(c);
    } else {
        res.conductance_per_meter = null_matrix;
    }
}

// Function: readMultipolarExpansion
multipolar_expansion_t parser_t::readMultipolarExpansion(json_value* multipolarExpansionPtr) {
    multipolar_expansion_t res;
    json_value* jvPtr = nullptr;
    bool found;

    this->core.get(multipolarExpansionPtr, J_MAT_MULTIWIRE_ME_INNER_REGION_BOX, jvPtr, found);
    if (!found) {
        WarnErrReport("Error reading multipolar expansion: innerRegionBox label not found", true);
    }
    res.inner_region = this->readInnnerRegionBox(jvPtr);

    this->core.get(multipolarExpansionPtr, J_MAT_MULTIWIRE_ME_ELECTRIC, jvPtr, found);
    if (!found) {
        WarnErrReport("Error reading multipolar expansion electric reconstruction not found", true);
    }
    res.electric = this->readFieldReconstruction(jvPtr);

    this->core.get(multipolarExpansionPtr, J_MAT_MULTIWIRE_ME_MAGNETIC, jvPtr, found);
    if (!found) {
        WarnErrReport("Error reading multipolar expansion magnetic reconstruction not found", true);
    }
    res.magnetic = this->readFieldReconstruction(jvPtr);

    return res;
}

// Function: readInnnerRegionBox
box_2d_t parser_t::readInnnerRegionBox(json_value* ptr) {
    box_2d_t inner_region;
    inner_region.min = this->getRealsAt(ptr, J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MIN);
    inner_region.max = this->getRealsAt(ptr, J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MAX);
    return inner_region;
}

// Function: readFieldReconstruction
std::vector<field_reconstruction_t> parser_t::readFieldReconstruction(json_value* ptr) {
    std::vector<field_reconstruction_t> res;
    int count = this->core.count(ptr);
    res.resize(count);

    for (int j = 0; j < count; ++j) {
        json_value* frPtr = nullptr;
        this->core.get_child(ptr, j, frPtr);

        res[j].inner_region_average_potential = this->getRealAt(frPtr, J_MAT_MULTIWIRE_MEFR_INNER_REGION_AVERAGE_POTENTIAL);
        res[j].expansion_center = this->getRealsAt(frPtr, J_MAT_MULTIWIRE_MEFR_EXPANSION_CENTER);
        res[j].conductor_potentials = this->getRealsAt(frPtr, J_MAT_MULTIWIRE_MEFR_CONDUCTOR_POTENTIALS);

        json_value* absPtr = nullptr;
        bool found;
        this->core.get(frPtr, J_MAT_MULTIWIRE_MEFR_AB, absPtr, found);
        if (!found) {
            WarnErrReport("Error reading multipolar expansion: ab label not found", true);
        }
        
        int abCount = this->core.count(absPtr);
        res[j].ab.resize(abCount);
        for (int i = 0; i < abCount; ++i) {
            json_value* abPtr = nullptr;
            this->core.get_child(absPtr, i, abPtr);
            
            double read_value;
            this->core.get(abPtr, "(1)", read_value);
            res[j].ab[i].a = read_value;
            this->core.get(abPtr, "(2)", read_value);
            res[j].ab[i].b = read_value;
        }
    }
    return res;
}

// Function: buildSegments
std::vector<segment_t> parser_t::buildSegments(const materialAssociation_t& j_cable, const Desplazamiento_t& despl) {
    std::vector<int> elemIds = j_cable.elementIds;
    polyline_t temp = this->mesh.getPolyline(elemIds[0]);
    std::vector<linel_t> linels = this->mesh.polylineToLinels(temp);
    
    std::vector<segment_t> res(linels.size());
    int prevOr = 0;

    for (size_t i = 0; i < linels.size(); ++i) {
        res[i].x = linels[i].cell[0];
        res[i].y = linels[i].cell[1];
        res[i].z = linels[i].cell[2];
        res[i].orientation = linels[i].orientation;

        if (prevOr == std::abs(res[i].orientation)) {
            res[i].dualBox = res[i-1].dualBox;
            res[i].d1 = res[i-1].d1;
            res[i].d2 = res[i-1].d2;
        } else {
            int absOr = std::abs(res[i].orientation);
            if (absOr == DIR_X) {
                res[i].dualBox = this->getdualBoxYZ(res[i], despl);
                res[i].d1 = despl.desY[res[i].y - 1];
                res[i].d2 = despl.desZ[res[i].z - 1];
            } else if (absOr == DIR_Y) {
                res[i].dualBox = this->getdualBoxZX(res[i], despl);
                res[i].d1 = despl.desZ[res[i].z - 1];
                res[i].d2 = despl.desX[res[i].x - 1];
            } else if (absOr == DIR_Z) {
                res[i].dualBox = this->getdualBoxXY(res[i], despl);
                res[i].d1 = despl.desX[res[i].x - 1];
                res[i].d2 = despl.desY[res[i].y - 1];
            }
        }
        prevOr = std::abs(res[i].orientation);
    }
    return res;
}

// Function: clip
int parser_t::clip(int i, int lo, int hi) {
    return std::max(lo, std::min(i, hi));
}

// Function: getdualBoxYZ
box_2d_t parser_t::getdualBoxYZ(const segment_t& segment, const Desplazamiento_t& despl) {
    box_2d_t res;
    int y0 = this->clip(segment.y - 1, 0, static_cast<int>(despl.desY.size()) - 1);
    int y1 = this->clip(segment.y, 0, static_cast<int>(despl.desY.size()) - 1);
    int z0 = this->clip(segment.z - 1, 0, static_cast<int>(despl.desZ.size()) - 1);
    int z1 = this->clip(segment.z, 0, static_cast<int>(despl.desZ.size()) - 1);

    res.min = {-0.5 * despl.desY[y0], -0.5 * despl.desZ[z0]};
    res.max = { 0.5 * despl.desY[y1],  0.5 * despl.desZ[z1]};
    return res;
}

// Function: getdualBoxXY
box_2d_t parser_t::getdualBoxXY(const segment_t& segment, const Desplazamiento_t& despl) {
    box_2d_t res;
    int x0 = this->clip(segment.x - 1, 0, static_cast<int>(despl.desX.size()) - 1);
    int x1 = this->clip(segment.x, 0, static_cast<int>(despl.desX.size()) - 1);
    int y0 = this->clip(segment.y - 1, 0, static_cast<int>(despl.desY.size()) - 1);
    int y1 = this->clip(segment.y, 0, static_cast<int>(despl.desY.size()) - 1);

    res.min = {-0.5 * despl.desX[x0], -0.5 * despl.desY[y0]};
    res.max = { 0.5 * despl.desX[x1],  0.5 * despl.desY[y1]};
    return res;
}

// Function: getdualBoxZX
box_2d_t parser_t::getdualBoxZX(const segment_t& segment, const Desplazamiento_t& despl) {
    box_2d_t res;
    int z0 = this->clip(segment.z - 1, 0, static_cast<int>(despl.desZ.size()) - 1);
    int z1 = this->clip(segment.z, 0, static_cast<int>(despl.desZ.size()) - 1);
    int x0 = this->clip(segment.x - 1, 0, static_cast<int>(despl.desX.size()) - 1);
    int x1 = this->clip(segment.x, 0, static_cast<int>(despl.desX.size()) - 1);

    res.min = {-0.5 * despl.desZ[z0], -0.5 * despl.desX[x0]};
    res.max = { 0.5 * despl.desZ[z1],  0.5 * despl.desX[x1]};
    return res;
}

// Function: buildStepSize
std::vector<rkind> parser_t::buildStepSize(const std::vector<segment_t>& segments, const Desplazamiento_t& despl) {
    std::vector<rkind> res(segments.size());
    for (size_t i = 0; i < segments.size(); ++i) {
        int or_val = std::abs(segments[i].orientation);
        if (or_val == DIR_X) {
            res[i] = despl.desX[segments[i].x];
        } else if (or_val == DIR_Y) {
            res[i] = despl.desY[segments[i].y];
        } else if (or_val == DIR_Z) {
            res[i] = despl.desZ[segments[i].z];
        }
    }
    return res;
}

// Function: readTransferImpedance
transfer_impedance_per_meter_t parser_t::readTransferImpedance(json_value* z) {
    transfer_impedance_per_meter_t res;
    std::string direction;

    if (this->existsAt(z, J_MAT_TRANSFER_IMPEDANCE_RESISTANCE)) {
        res.resistive_term = this->getRealAt(z, J_MAT_TRANSFER_IMPEDANCE_RESISTANCE);
    }

    if (this->existsAt(z, J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE)) {
        res.inductive_term = this->getRealAt(z, J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE);
    }

    if (this->existsAt(z, J_MAT_TRANSFER_IMPEDANCE_DIRECTION)) {
        direction = trimAdjustl(this->getStrAt(z, J_MAT_TRANSFER_IMPEDANCE_DIRECTION));
    } else {
        WarnErrReport("Error reading material: direction of transferImpedancePerMeter missing", true);
    }

    if (direction == "inwards") {
        res.direction = TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
    } else if (direction == "outwards") {
        res.direction = TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS;
    } else if (direction == "both") {
        res.direction = TRANSFER_IMPEDANCE_DIRECTION_BOTH;
    }

    // Block scope for poles and residues
    {
        json_value* poles = nullptr;
        json_value* residues = nullptr;
        json_value* p = nullptr;
        json_value* r = nullptr;
        int n = 0;
        double read_value;

        if (this->existsAt(z, J_MAT_TRANSFER_IMPEDANCE_POLES)) {
            n = this->getIntAt(z, J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES);
            res.poles.resize(n);
            res.residues.resize(n);
            
            this->core.get(z, J_MAT_TRANSFER_IMPEDANCE_POLES, poles);
            this->core.get(z, J_MAT_TRANSFER_IMPEDANCE_RESIDUES, residues);

            for (int i = 0; i < n; ++i) {
                this->core.get_child(poles, i, p);
                this->core.get(p, "(1)", read_value);
                res.poles[i].re = read_value;
                this->core.get(p, "(2)", read_value);
                res.poles[i].im = read_value;

                this->core.get_child(residues, i, r);
                this->core.get(r, "(1)", read_value);
                res.residues[i].re = read_value;
                this->core.get(r, "(2)", read_value);
                res.residues[i].im = read_value;
            }
        } else {
            res.poles.clear();
            res.residues.clear();
        }
    }

    return res;
}

// Function: noTransferImpedance
transfer_impedance_per_meter_t parser_t::noTransferImpedance() {
    transfer_impedance_per_meter_t res;
    res.poles.clear();
    res.residues.clear();
    return res;
}

// Function: assignDisplacement
std::vector<rkind> parser_t::assignDisplacement(const Desplazamiento_t& desp, int axis) {
    std::vector<rkind> res;
    if (axis == 1) {
        res = desp.desX;
    } else if (axis == 2) {
        res = desp.desY;
    } else if (axis == 3) {
        res = desp.desZ;
    }
    return res;
}

// Function: findOrientation
int parser_t::findOrientation(const coordinate_t& coordDiference) {
    int res = 0;
    for (int i = 0; i < 3; ++i) {
        if (coordDiference.position[i] != 0) {
            res = coordDiference.position[i] / std::abs(coordDiference.position[i]);
        }
    }
    return res;
}

// Function: findDirection
int parser_t::findDirection(const coordinate_t& coordDiference) {
    int res = 0;
    for (int i = 0; i < 3; ++i) {
        if (coordDiference.position[i] != 0) {
            res = i + 1; // Fortran is 1-based
        }
    }
    return res;
}

// Function: handleFoundAndDefault
void parser_t::handleFoundAndDefault(const std::string& path, bool found, bool defaultPresent) {
    if (!found && !defaultPresent) {
        std::string errorMsg = "ERROR expecting a value at: " + path;
        WarnErrReport(errorMsg, true);
    }
}

// Function: getLogicalAt
bool parser_t::getLogicalAt(json_value* place, const std::string& path, bool* found, bool defaultVal) {
    bool res;
    bool localFound;
    
    // Assuming core.get has this signature
    this->core.get(place, path, res, localFound, defaultVal);
    
    if (found) {
        *found = localFound;
    } else {
        this->handleFoundAndDefault(path, localFound, true); // present(default) is true if default is provided
    }
    return res;
}

// Function: getIntAt
int parser_t::getIntAt(json_value* place, const std::string& path, bool* found, int defaultVal) {
    int res;
    bool localFound;
    
    // Assuming core.get has this signature
    this->core.get(place, path, res, localFound, defaultVal);
    
    if (found) {
        *found = localFound;
    } else {
        this->handleFoundAndDefault(path, localFound, true); // present(default) is true if default is provided
    }
    return res;
}

// This file is a continuation of the parser_t class implementation.
// Includes from previous chunks should be present here.
// #include "parser_t.hpp"
// #include "json_value.hpp"
// #include "cell_region_t.hpp"
// #include "cell_interval_t.hpp"
// #include "WarnErrReport.hpp"
// #include <vector>
// #include <string>
// #include <optional>
// #include <memory>

// Assuming RKIND, JSON_CK, BUFSIZE, J_ELEMENTIDS, CELL_TYPE_VOXEL are defined in headers

namespace parser_ns { // Or whatever namespace the module maps to

    // ... (Previous methods of parser_t class) ...

    // Helper function declaration (assumed to be defined elsewhere or in header)
    void handleFoundAndDefault(const std::string& path, bool localFound, bool hasDefault);

    // Note: The Fortran code uses member functions. In C++, these are methods of parser_t.
    // The 'this' pointer is implicit.

    // getAt function (returns double, optional found, optional default)
    double parser_t::getAt(const json_value& place, const std::string& path, 
                           std::optional<bool>& found, std::optional<double> defaultVal) {
        double res;
        bool localFound;
        this->core->get(place, path, res, localFound, defaultVal);
        if (found.has_value()) {
            found.value() = localFound;
        } else {
            handleFoundAndDefault(path, localFound, defaultVal.has_value());
        }
        return res;
    }

    // getAt function (returns vector<int>, optional found)
    std::vector<int> parser_t::getAt(const json_value& place, const std::string& path, 
                                     std::optional<bool>& found) {
        std::vector<int> res;
        bool localFound;
        this->core->get(place, path, res, localFound);
        if (found.has_value()) {
            found.value() = localFound;
        } else {
            handleFoundAndDefault(path, localFound, false);
        }
        return res;
    }

    // getAt function (returns double, optional found, optional default)
    double parser_t::getAt(const json_value& place, const std::string& path, 
                           std::optional<bool>& found, std::optional<double> defaultVal) {
        double res;
        bool localFound;
        this->core->get(place, path, res, localFound, defaultVal);
        if (found.has_value()) {
            found.value() = localFound;
        } else {
            handleFoundAndDefault(path, localFound, defaultVal.has_value());
        }
        return res;
    }

    // getAt function (returns vector<double>, optional found)
    std::vector<double> parser_t::getAt(const json_value& place, const std::string& path, 
                                        std::optional<bool>& found) {
        std::vector<double> res;
        bool localFound;
        this->core->get(place, path, res, localFound); // Note: Fortran had 'localfound' (lowercase f)
        if (found.has_value()) {
            found.value() = localFound;
        } else {
            handleFoundAndDefault(path, localFound, false);
        }
        return res;
    }

    // getMatrixAt function
    std::vector<std::vector<double>> parser_t::getMatrixAt(const json_value& place, const std::string& path, 
                                                           std::optional<bool>& found) {
        std::vector<std::vector<double>> res;
        bool localFound;
        json_value* matrix = nullptr;
        
        this->core->get(place, path, matrix, localFound);
        if (found.has_value()) {
            found.value() = localFound;
        } else {
            handleFoundAndDefault(path, localFound, false);
        }
        
        int vartype = 0;
        int nr = 0;
        this->core->info(matrix, vartype, nr);
        
        res.resize(nr, std::vector<double>(nr));
        
        for (int i = 1; i <= nr; ++i) {
            json_value* row = nullptr;
            this->core->get_child(matrix, i, row);
            std::vector<double> res_row;
            this->core->get(row, res_row);
            // Fortran is 1-based, C++ is 0-based. Adjust indices.
            for (int j = 0; j < static_cast<int>(res_row.size()); ++j) {
                res[i-1][j] = res_row[j];
            }
        }
        // Note: Fortran comment "need to check if not found" is ignored as it's not code.
        return res;
    }

    // getStrAt function
    std::string parser_t::getStrAt(const json_value& place, const std::string& path, 
                                   std::optional<bool>& found, std::optional<std::string> defaultVal) {
        std::string res;
        bool localFound;
        this->core->get(place, path, res, localFound, defaultVal);
        if (found.has_value()) {
            found.value() = localFound;
        } else {
            handleFoundAndDefault(path, localFound, defaultVal.has_value());
        }
        return res;
    }

    // existsAt function
    bool parser_t::existsAt(const json_value& place, const std::string& path) {
        bool found = false;
        this->core->info(place, path, found);
        return found;
    }

    // dimensionAt function
    int parser_t::dimensionAt(const json_value& place, const std::string& path) {
        int n_children = 0;
        this->core->info(place, path, n_children);
        return n_children;
    }

    // jsonValueFilterByKeyValues function
    std::vector<json_value_ptr_t> parser_t::jsonValueFilterByKeyValues(const json_value& srcs, 
                                                                       const std::string& key, 
                                                                       const std::vector<std::string>& values) {
        std::vector<json_value_ptr_t> res;
        res.reserve(values.size() * 10); // Heuristic reserve
        
        for (const auto& val : values) {
            std::vector<json_value_ptr_t> foundEntries = this->jsonValueFilterByKeyValue(srcs, key, val);
            if (!foundEntries.empty()) {
                res.insert(res.end(), foundEntries.begin(), foundEntries.end());
            }
        }
        return res;
    }

    // jsonValueFilterByKeyValue function
    std::vector<json_value_ptr_t> parser_t::jsonValueFilterByKeyValue(const json_value& place, 
                                                                      const std::string& key, 
                                                                      const std::string& value) {
        std::vector<json_value_ptr_t> res;
        int n = 0;
        int count = this->core->count(place);
        
        // Precounting
        for (int i = 1; i <= count; ++i) {
            json_value* src = nullptr;
            this->core->get_child(place, i, src);
            
            bool found = false;
            std::string typeStr = this->getStrAt(*src, key, found);
            
            if (!found) {
                std::string errorMsg = "Key: " + key + " not found while doing value filter.";
                WarnErrReport(errorMsg, true);
            }
            
            if (found && typeStr == value) {
                n++;
            }
        }
        
        res.resize(n);
        int j = 0;
        for (int i = 1; i <= count; ++i) {
            json_value* src = nullptr;
            this->core->get_child(place, i, src);
            
            bool found = false;
            std::string typeStr = this->getStrAt(*src, key, found);
            
            if (found && typeStr == value) {
                res[j].p = src;
                j++;
            }
        }
        return res;
    }

    // getSingleVolumeInElementsIds function
    std::vector<cell_interval_t> parser_t::getSingleVolumeInElementsIds(const json_value& pw) {
        bool found = false;
        std::vector<int> elemIds = this->getIntsAt(pw, J_ELEMENTIDS, found);
        
        if (!found) {
            WarnErrReport("Error reading single volume elementIds label not found.", true);
            return std::vector<cell_interval_t>(0);
        }
        
        if (elemIds.empty()) {
            WarnErrReport("Entity elementIds must not be empty.", true);
            return std::vector<cell_interval_t>(0);
        }
        
        if (elemIds.size() != 1) {
            WarnErrReport("Entity must contain a single elementId.", true);
            return std::vector<cell_interval_t>(0);
        }
        
        bool cellFound = false;
        cell_region_t cellRegion = this->mesh.getCellRegion(elemIds[0], cellFound);
        
        if (!cellFound) {
            std::string errorMsg = "Entity elementId " + std::to_string(elemIds[0]) + " not found.";
            WarnErrReport(errorMsg, true);
            return std::vector<cell_interval_t>(0);
        }
        
        std::vector<cell_interval_t> res = cellRegion.getIntervalsOfType(CELL_TYPE_VOXEL);
        
        if (res.size() != 1) {
            std::string errorMsg = "Entity must contain a single cell region defining a volume.";
            WarnErrReport(errorMsg, true);
            return std::vector<cell_interval_t>(0);
        }
        
        return res;
    }

} // namespace parser_ns