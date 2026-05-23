#include "smbjson_m.h"
#include "NFDETypes_extension_m.h"
#include <iostream>

namespace smbjson {

    // ---- JSON accessor implementations ----

    bool parser_t::getLogicalAt(const jmod::json_value* val, const std::string& key, bool default_val, bool* foundOut) {
        jmod::json_value* ptr = nullptr;
        bool found = false;
        core->get(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->data.get<bool>();
        return default_val;
    }

    int parser_t::getIntAt(const jmod::json_value* val, const std::string& key, int default_val, bool* foundOut) {
        jmod::json_value* ptr = nullptr;
        bool found = false;
        core->get(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->data.get<int>();
        return default_val;
    }

    std::vector<int> parser_t::getIntsAt(const jmod::json_value* val, const std::string& key, bool* foundOut) {
        std::vector<int> res;
        jmod::json_value* arr = nullptr;
        bool found = false;
        core->get(val, key, arr, found);
        if (foundOut) *foundOut = found;
        if (found && arr) {
            int n = core->count(arr);
            res.resize(n);
            for (int i = 0; i < n; ++i) {
                jmod::json_value* child = nullptr;
                core->get_child(arr, i + 1, child);
                if (child) res[i] = child->data.get<int>();
            }
        }
        return res;
    }

    double parser_t::getRealAt(const jmod::json_value* val, const std::string& key, double default_val, bool* foundOut) {
        jmod::json_value* ptr = nullptr;
        bool found = false;
        core->get(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->data.get<double>();
        return default_val;
    }

    std::vector<double> parser_t::getRealsAt(const jmod::json_value* val, const std::string& key, bool* foundOut) {
        std::vector<double> res;
        jmod::json_value* arr = nullptr;
        bool found = false;
        core->get(val, key, arr, found);
        if (foundOut) *foundOut = found;
        if (found && arr) {
            int n = core->count(arr);
            res.resize(n);
            for (int i = 0; i < n; ++i) {
                jmod::json_value* child = nullptr;
                core->get_child(arr, i + 1, child);
                if (child) res[i] = child->data.get<double>();
            }
        }
        return res;
    }

    std::vector<std::vector<double>> parser_t::getMatrixAt(const jmod::json_value* val, const std::string& key, bool* foundOut) {
        std::vector<std::vector<double>> res;
        jmod::json_value* arr = nullptr;
        bool found = false;
        core->get(val, key, arr, found);
        if (foundOut) *foundOut = found;
        if (found && arr) {
            int n = core->count(arr);
            for (int i = 0; i < n; ++i) {
                jmod::json_value* child = nullptr;
                core->get_child(arr, i + 1, child);
                res.push_back(getRealsAt(child, ""));
            }
        }
        return res;
    }

    std::string parser_t::getStrAt(const jmod::json_value* val, const std::string& key, const std::string& default_val, bool* foundOut) {
        jmod::json_value* ptr = nullptr;
        bool found = false;
        core->get(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->data.get<std::string>();
        return default_val;
    }

    bool parser_t::existsAt(const jmod::json_value* val, const std::string& key) {
        jmod::json_value* ptr = nullptr;
        bool found = false;
        core->get(val, key, ptr, found);
        return found;
    }

    int parser_t::dimensionAt(const jmod::json_value* val, const std::string& key) {
        jmod::json_value* ptr = nullptr;
        bool found = false;
        core->get(val, key, ptr, found);
        if (found) return core->count(ptr);
        return 0;
    }

    parser_t::domain_t parser_t::getDomain(const jmod::json_value* place, const std::string& path) {
        domain_t res;
        res.tstart = getRealAt(place, path + ".tstart", 0.0);
        res.tstop = getRealAt(place, path + ".tstop", 0.0);
        res.tstep = getRealAt(place, path + ".tstep", 0.0);
        res.fstart = getRealAt(place, path + ".fstart", 0.0);
        res.fstop = getRealAt(place, path + ".fstop", 0.0);
        res.fstep = getRealAt(place, path + ".fstep", 0.0);
        res.filename = getStrAt(place, path + ".filename", "");
        res.type1 = getIntAt(place, path + ".type1", NFDE::NP_T1_PLAIN);
        res.type2 = getIntAt(place, path + ".type2", NFDE::NP_T2_TIME);
        res.isLogarithmicFrequencySpacing = getLogicalAt(place, path + ".isLogarithmicFrequencySpacing", false);
        return res;
    }

    // ---- Missing method implementations ----

    std::string parser_t::buildTagName(int matId, int elementId) {
        return std::to_string(matId) + "@" + std::to_string(elementId);
    }

    std::vector<parser_t::materialAssociation_t> parser_t::getMaterialAssociations(
        const std::vector<std::string>& materialTypes,
        const std::vector<std::string>& elementLabels) {
        std::vector<materialAssociation_t> res;
        jmod::json_value* allMatAss = nullptr;
        bool found = false;
        core->get(root, jlbl::J_MATERIAL_ASSOCIATIONS, allMatAss, found);
        if (!found) return res;

        int nMatAss = core->count(allMatAss);
        for (int i = 0; i < nMatAss; ++i) {
            jmod::json_value* mAPtr = nullptr;
            core->get_child(allMatAss, i + 1, mAPtr);

            materialAssociation_t mA = parseMaterialAssociation(mAPtr);
            jmod::json_value_ptr_t mat = matTable.getId(mA.materialId);
            if (!mat) continue;
            std::string matType = getStrAt(mat, jlbl::J_TYPE, "");

            bool typeMatch = false;
            for (auto& t : materialTypes) {
                if (matType == t) { typeMatch = true; break; }
            }
            if (!typeMatch) continue;

            // Element label filter (simplified, not used for wire/cable)
            bool labelMatch = true;
            if (!elementLabels.empty()) {
                labelMatch = false;
                for (int eid : mA.elementIds) {
                    bool cRf = false;
                    Cell::cell_region_t cR = mesh.getCellRegion(eid, cRf);
                    if (!cRf) continue;
                    jmod::json_value_ptr_t elm = elementTable.getId(eid);
                    if (!elm) continue;
                    std::string elemType = getStrAt(elm, jlbl::J_TYPE, "");
                    std::string elemSubtype = getStrAt(elm, jlbl::J_SUBTYPE, "");
                    for (auto& el : elementLabels) {
                        bool negative = (el[0] == '-');
                        std::string target = negative ? el.substr(1) : el;
                        bool matches = (elemType == target || elemSubtype == target);
                        if (negative) { labelMatch = labelMatch && !matches; }
                        else { labelMatch = labelMatch || matches; }
                    }
                }
            }
            if (!labelMatch) continue;

            mA.matAssType = matType;
            res.push_back(mA);
        }
        return res;
    }

    parser_t::materialAssociation_t parser_t::parseMaterialAssociation(const jmod::json_value* matAss) {
        materialAssociation_t mA;
        bool found = false;

        mA.materialId = getIntAt(matAss, jlbl::J_MATERIAL_ID, 0);
        if (found) mA.materialId = mA.materialId;
        mA.elementIds = getIntsAt(matAss, jlbl::J_ELEMENTIDS);
        mA.name = getStrAt(matAss, jlbl::J_NAME, "");
        mA.initialTerminalId = getIntAt(matAss, jlbl::J_MAT_ASS_CAB_INI_TERM_ID, -1);
        mA.endTerminalId = getIntAt(matAss, jlbl::J_MAT_ASS_CAB_END_TERM_ID, -1);
        mA.initialConnectorId = getIntAt(matAss, jlbl::J_MAT_ASS_CAB_INI_CONN_ID, -1);
        mA.endConnectorId = getIntAt(matAss, jlbl::J_MAT_ASS_CAB_END_CONN_ID, -1);
        mA.containedWithinElementId = getIntAt(matAss, jlbl::J_MAT_ASS_CAB_CONTAINED_WITHIN_ID, -1);
        mA.hasTotalResistance = existsAt(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE);
        if (mA.hasTotalResistance) {
            int dim = dimensionAt(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE);
            if (dim == 0) {
                mA.totalResistance = {getRealAt(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE, 0.0)};
            } else {
                mA.totalResistance = getRealsAt(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE);
            }
        }
        return mA;
    }

    std::vector<Cell::cell_interval_t> parser_t::getSingleVolumeInElementsIds(const jmod::json_value* pw) {
        std::vector<Cell::cell_interval_t> res;
        std::vector<int> elemIds = getIntsAt(pw, jlbl::J_ELEMENTIDS);
        if (elemIds.size() >= 1) {
            bool cRf = false;
            Cell::cell_region_t cR = mesh.getCellRegion(int(elemIds[0]), cRf);
            if (cRf) {
                for (auto& iv : cR.intervals) {
                    if (iv.getType() == Cell::CELL_TYPE_VOXEL) res.push_back(iv);
                }
            }
        }
        return res;
    }

    std::vector<jmod::json_value*> parser_t::jsonValueFilterByKeyValue(
        const jmod::json_value* place, const std::string& key, const std::string& value) {
        return jsonValueFilterByKeyValues(place, key, {value});
    }

    std::vector<jmod::json_value*> parser_t::jsonValueFilterByKeyValues(
        const jmod::json_value* place, const std::string& key, const std::vector<std::string>& values) {
        std::vector<jmod::json_value*> res;
        if (!place) return res;
        int n = core->count(place);
        for (int i = 0; i < n; ++i) {
            jmod::json_value* child = nullptr;
            core->get_child(place, i + 1, child);
            std::string v = getStrAt(child, key, "");
            for (auto& target : values) {
                if (v == target) { res.push_back(child); break; }
            }
        }
        return res;
    }

    int parser_t::labelToBoundaryPlace(const std::string& str) {
        if (str == jlbl::J_BND_XL) return NFDE::F_XL - 1;
        if (str == jlbl::J_BND_XU) return NFDE::F_XU - 1;
        if (str == jlbl::J_BND_YL) return NFDE::F_YL - 1;
        if (str == jlbl::J_BND_YU) return NFDE::F_YU - 1;
        if (str == jlbl::J_BND_ZL) return NFDE::F_ZL - 1;
        if (str == jlbl::J_BND_ZU) return NFDE::F_ZU - 1;
        return -1;
    }

    int parser_t::labelToBoundaryType(const std::string& str) {
        if (str == jlbl::J_BND_TYPE_MUR) return NFDE::F_MUR;
        if (str == jlbl::J_BND_TYPE_PEC) return NFDE::F_PEC;
        if (str == jlbl::J_BND_TYPE_PMC) return NFDE::F_PMC;
        if (str == jlbl::J_BND_TYPE_PML) return NFDE::F_PML;
        if (str == jlbl::J_BND_TYPE_PERIODIC) return NFDE::F_PER;
    }

    void parser_t::readThinWires(NFDE::ThinWires_t& res, NFDE::MasSondas_t& sonda) {
        // Stub: thin wires parsing deferred (MTLN feature)
    }

    parser_t::parser_t(const std::string& filename) {
        this->filename = filename;
        
        jsonfile = new jmod::json_file();
        jsonfile->initialize();
        if (jsonfile->failed()) {
            Report::WarnErrReport("Failed to initialize JSONfile", true);
            return;
        }

        jsonfile->load(filename);
        if (jsonfile->failed()) {
            Report::WarnErrReport("Failed to load JSON file: " + filename, true);
            return;
        }

        core = new jmod::json_core();
        jsonfile->get_core(*core);
        root = new jmod::json_value();
        jsonfile->get(".", *root);

        isInitialized = true;
    }

    NFDE::Parseador_t parser_t::readProblemDescription() {
        NFDE::Parseador_t res;
        
        mesh = readMesh();
        matTable = IdTable::IdChildTable_t::ctor(*core, *root, jlbl::J_MATERIALS);
        elementTable = IdTable::IdChildTable_t::ctor(*core, *root, std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_ELEMENTS));
        
        NFDETypes_extension_m::initializeProblemDescription(res);
        
        res.switches = readAdditionalArguments();
        
        // Basics
        res.general = new NFDE::NFDEGeneral_t(readGeneral());
        res.matriz = new NFDE::MatrizMedios_t(readMediaMatrix());
        res.despl = new NFDE::Desplazamiento_t(readGrid());
        res.front = new NFDE::Frontera_t(readBoundary());
        
        // Materials
        readBackgroundMaterial(*res.Mats);
        res.pecRegs = new NFDE::PECRegions_t(readPECRegions());
        res.pmcRegs = new NFDE::PECRegions_t(readPMCRegions());
        res.DielRegs = new NFDE::DielectricRegions_t(readDielectricRegions());
        res.LossyThinSurfs = new NFDE::LossyThinSurfaces_t(readLossyThinSurfaces());
        
        // Sources
        res.plnSrc = new NFDE::PlaneWaves_t(readPlanewaves());
        res.nodSrc = new NFDE::NodSource_t(readNodalSources());
        
        // Probes
        res.oldSONDA = new NFDE::Sondas_t(readProbes());
        res.Sonda = new NFDE::MasSondas_t(readMoreProbes());
        res.BloquePrb = new NFDE::BloqueProbes_t(readBlockProbes());
        res.VolPrb = new NFDE::VolProbes_t(readVolumicProbes());
        
        // Conformal elements
        res.conformalRegs = new NFDE::ConformalPECRegions_t(readConformalRegions());

        // Thin elements
#ifdef CompileWithMTLN 
        res.mtln = new NFDE::mtln_t(readMTLN());
#else
        readThinWires(*res.tWires, *res.Sonda);
#endif
        res.tSlots = new NFDE::ThinSlots_t(readThinSlots());

        return res;
    }

    Mesh::mesh_t parser_t::readMesh() {
        Mesh::mesh_t res;
        addCoordinates(res);
        addElements(res);
        return res;
    }

    void parser_t::addCoordinates(Mesh::mesh_t& mesh) {
        jmod::json_value* jcs = nullptr;
        jmod::json_value* jc = nullptr;
        bool found = false;
        
        core->get(root, std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_COORDINATES), jcs, found);
        if (found) {
            int numberOfCoordinates = core->count(jcs);
            mesh.allocateCoordinates(50 * numberOfCoordinates);
            for (int i = 1; i <= numberOfCoordinates; ++i) {
                core->get_child(jcs, i, jc);
                int id = getIntAt(jc, jlbl::J_ID, 0);
                std::vector<double> pos = getRealsAt(jc, jlbl::J_COORDINATE_POS);
                Mesh::coordinate_t c;
                std::copy(pos.begin(), pos.end(), c.position);
                mesh.addCoordinate(id, c);
            }
        }
    }

    void parser_t::addElements(Mesh::mesh_t& mesh) {
        std::string elementType;
        jmod::json_value* jes = nullptr;
        jmod::json_value* je = nullptr;
        bool found = false;
        
        core->get(root, std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_ELEMENTS), jes, found);
        int numberOfElements = core->count(jes);
        mesh.allocateElements(50 * numberOfElements);
            
        if (found) {
            for (int i = 1; i <= numberOfElements; ++i) {
                core->get_child(jes, i, je);
                int id = getIntAt(je, jlbl::J_ID, 0);
                elementType = getStrAt(je, jlbl::J_TYPE);
                
                if (elementType == jlbl::J_ELEM_TYPE_NODE) {
                    std::vector<int> coordIds = getIntsAt(je, jlbl::J_COORDINATE_IDS);
                    Mesh::node_t node;
                    node.coordIds = coordIds;
                    mesh.addElement(id, node);
                } else if (elementType == jlbl::J_ELEM_TYPE_POLYLINE) {
                    std::vector<int> coordIds = getIntsAt(je, jlbl::J_COORDINATE_IDS);
                    Mesh::polyline_t polyline;
                    polyline.coordIds = coordIds;
                    mesh.addElement(id, polyline);
                } else if (elementType == jlbl::J_ELEM_TYPE_CELL) {
                    bool isConformal = false;
                    jmod::json_value* triangles = nullptr;
                    core->get(je, jlbl::J_CONF_VOLUME_TRIANGLES, triangles, isConformal);
                    
                    if (!isConformal) {
                        Cell::cell_region_t cR;
                        cR.intervals = readCellIntervals(je, jlbl::J_CELL_INTERVALS);
                        mesh.addCellRegion(id, cR);
                    } else {
                        Mesh::conformal_region_t cV;
                        cV.triangles = readTriangles(je, jlbl::J_CONF_VOLUME_TRIANGLES);
                        for (int k = 1; k <= cV.triangles.size(); ++k) {
                            for (int j = 1; j <= 3; ++j) {
                                bool crFound = false;
                                Mesh::coordinate_t c = mesh.getCoordinate(cV.triangles[k-1].vertices[j-1].id, crFound);
                                cV.triangles[k-1].vertices[j-1].position[0] = c.position[0];
                                cV.triangles[k-1].vertices[j-1].position[1] = c.position[1];
                                cV.triangles[k-1].vertices[j-1].position[2] = c.position[2];
                            }
                        }
                        cV.intervals = readCellIntervals(je, jlbl::J_CELL_INTERVALS);
                        std::string subtype = getStrAt(je, jlbl::J_SUBTYPE);

                        if (subtype == jlbl::J_CONF_SUBTYPE_VOLUME) cV.type = Mesh::REGION_TYPE_VOLUME;
                        if (subtype == jlbl::J_CONF_SUBTYPE_SURFACE) cV.type = Mesh::REGION_TYPE_SURFACE;

                        mesh.addConformalRegion(id, cV);
                    }
                } else {
                    Report::WarnErrReport("Invalid element type", false);
                }
            }
        }
    }

    std::vector<Cell::cell_interval_t> parser_t::readCellIntervals(const jmod::json_value* place, const std::string& path) {
        jmod::json_value* intervalsPlace = nullptr;
        jmod::json_value* interval = nullptr;
        bool containsInterval = false;
        
        core->get(place, path, intervalsPlace, containsInterval);
        if (!containsInterval) {
            return std::vector<Cell::cell_interval_t>();
        }
        int nIntervals = core->count(intervalsPlace);
        std::vector<Cell::cell_interval_t> res(nIntervals);
        for (int i = 1; i <= nIntervals; ++i) {
            core->get_child(intervalsPlace, i, interval);
            std::vector<double> cellIni = getRealsAt(interval, "(1)");
            std::vector<double> cellEnd = getRealsAt(interval, "(2)");
            res[i-1].ini.cell[0] = cellIni[0];
            res[i-1].ini.cell[1] = cellIni[1];
            res[i-1].ini.cell[2] = cellIni[2];
            res[i-1].end.cell[0] = cellEnd[0];
            res[i-1].end.cell[1] = cellEnd[1];
            res[i-1].end.cell[2] = cellEnd[2];
        }
        return res;
    }

    std::vector<Conf::triangle_t> parser_t::readTriangles(const jmod::json_value* place, const std::string& path) {
        jmod::json_value* triangles = nullptr;
        jmod::json_value* triangle_ptr = nullptr;
        bool containsTriangles = false;
        
        core->get(place, path, triangles, containsTriangles);
        if (!containsTriangles) {
            return std::vector<Conf::triangle_t>();
        }
        int nTriangles = core->count(triangles);
        std::vector<Conf::triangle_t> res(nTriangles);
        for (int i = 1; i <= nTriangles; ++i) {
            core->get_child(triangles, i, triangle_ptr);
            std::vector<double> triangle = getRealsAt(triangle_ptr, ""); // Assuming get returns vector for array
            // Note: json-fortran get to vector might need specific handling, assuming helper exists
            for (int j = 1; j <= 3; ++j) {
                res[i-1].vertices[j-1].id = static_cast<int>(triangle[j-1]);
            }
        }
        return res;
    }

    std::string parser_t::readAdditionalArguments() {
        return getStrAt(root, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_ADDITIONAL_ARGUMENTS), " ");
    }

    NFDE::NFDEGeneral_t parser_t::readGeneral() {
        NFDE::NFDEGeneral_t res;
        res.dt = getRealAt(root, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_TIME_STEP), 0.0);
        if (res.dt < 0) Report::WarnErrReport("timStep cannot be negative", true);
        res.nmax = getRealAt(root, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_NUMBER_OF_STEPS), 0.0);
        if (res.nmax <= 0) Report::WarnErrReport("numberOfSteps has to be positive", true);
        res.mtlnProblem = getLogicalAt(root, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_MTLN_PROBLEM), false);
        return res;
    }

    NFDE::MatrizMedios_t parser_t::readMediaMatrix() {
        NFDE::MatrizMedios_t res;
        std::string P = std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_GRID) + "." + jlbl::J_GRID_NUMBER_OF_CELLS;
        res.totalX = getIntAt(root, P + ".(0)", 0) + 1;
        res.totalY = getIntAt(root, P + ".(1)", 0) + 1;
        res.totalZ = getIntAt(root, P + ".(2)", 0) + 1;
        return res;
    }

    NFDE::Desplazamiento_t parser_t::readGrid() {
        NFDE::Desplazamiento_t res;
        std::string P = std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_GRID);
        
        int nX = getIntAt(root, std::string(P) + "." + std::string(jlbl::J_GRID_NUMBER_OF_CELLS) + ".(0)", 0);
        int nY = getIntAt(root, std::string(P) + "." + std::string(jlbl::J_GRID_NUMBER_OF_CELLS) + ".(1)", 0);
        int nZ = getIntAt(root, std::string(P) + "." + std::string(jlbl::J_GRID_NUMBER_OF_CELLS) + ".(2)", 0);

        res.nX = nX;
        res.nY = nY;
        res.nZ = nZ;

        // Helper lambda for assignDes logic
        auto assignDes = [&](const std::string& path, std::vector<double>& dest, int& n) {
            std::vector<double> vec = getRealsAt(root, path);
            if (vec.empty()) {
                Report::WarnErrReport("Error reading grid: steps not found at path: " + path, true);
            }
            if (vec.size() != 1 && vec.size() != static_cast<size_t>(n)) {
                Report::WarnErrReport("Error reading grid: steps must be arrays of size 1 (for regular grids, false) || size equal to the number of cells.", true);
            }

            if (vec.size() == 1) {
                n = 1;
                dest = vec;
            } else {
                dest = vec;
            }
        };

        assignDes(std::string(P) + "." + std::string(jlbl::J_GRID_STEPS) + ".x", res.desX, res.nX);
        assignDes(std::string(P) + "." + std::string(jlbl::J_GRID_STEPS) + ".y", res.desY, res.nY);
        assignDes(std::string(P) + "." + std::string(jlbl::J_GRID_STEPS) + ".z", res.desZ, res.nZ);

        res.originx = getRealAt(root, std::string(P) + "." + std::string(jlbl::J_GRID_ORIGIN) + ".(0)", 0.0);
        res.originy = getRealAt(root, std::string(P) + "." + std::string(jlbl::J_GRID_ORIGIN) + ".(1)", 0.0);
        res.originz = getRealAt(root, std::string(P) + "." + std::string(jlbl::J_GRID_ORIGIN) + ".(2)", 0.0);

        res.mx1 = 0;
        res.my1 = 0;
        res.mz1 = 0;
        res.mx2 = nX;
        res.my2 = nY;
        res.mz2 = nZ;

        return res;
    }

    NFDE::Frontera_t parser_t::readBoundary() {
        NFDE::Frontera_t res;
        std::string bdrType;
        jmod::json_value* bdrs = nullptr;
        bool found = false;
        
        core->get(root, jlbl::J_BOUNDARY, bdrs, found);
        if (!found) {
            Report::WarnErrReport("Error reading boundary: " + std::string(jlbl::J_BOUNDARY) + " not found.", true);
        }
        
        {
            bool foundLocal = false;
            bdrType = getStrAt(bdrs, std::string(jlbl::J_BND_ALL) + "." + std::string(jlbl::J_TYPE), "", &foundLocal);
            if (foundLocal) {
                // Assuming labelToBoundaryType returns int and res.tipoFrontera is an array/vector
                // This part depends heavily on specific struct definitions not fully provided
                // Placeholder logic
                for(int i=0; i<6; ++i) {
                     res.tipoFrontera[i] = labelToBoundaryType(bdrType);
                }
                bool allPML = true;
                 for(int _pf=0; _pf<6; ++_pf) { if(res.tipoFrontera[_pf] != NFDE::F_PML) allPML = false; }
                 if (allPML) {
                    // res.propiedadesPML = readPMLProperties(jlbl::J_BOUNDARY "." jlbl::J_BND_ALL);
                }
                return res;
            }
        }
         
        {
            std::vector<std::string> placeLabels = {jlbl::J_BND_XL, jlbl::J_BND_XU, jlbl::J_BND_YL, jlbl::J_BND_YU, jlbl::J_BND_ZL, jlbl::J_BND_ZU};
            for (int i = 0; i < 6; ++i) {
                bool foundLocal = false;
                bdrType = getStrAt(bdrs, std::string(placeLabels[i]) + "." + std::string(jlbl::J_TYPE), "", &foundLocal);
                if (!foundLocal) {
                    Report::WarnErrReport("ERROR reading boundary: " + placeLabels[i] + " || " + jlbl::J_BND_ALL + " not found.", true);
                }
                int j = labelToBoundaryPlace(placeLabels[i]);
                res.tipoFrontera[j] = labelToBoundaryType(bdrType);
                if (res.tipoFrontera[j] == NFDE::F_PML) {
                    // res.propiedadesPML[j] = readPMLProperties(jlbl::J_BOUNDARY "." placeLabels[i]);
                }
            }
        }

        return res;
    }

    void parser_t::readBackgroundMaterial(NFDE::Materials_t& mats) {
        bool found = false;
        double val = getRealAt(root, std::string(jlbl::J_BACKGROUND) + "." + std::string(jlbl::J_BKG_ABS_PERMITTIVITY), 0.0, &found);
        if (found) mats.Mats[0].eps = val;

        found = false;
        val = getRealAt(root, std::string(jlbl::J_BACKGROUND) + "." + std::string(jlbl::J_BKG_ABS_PERMEABILITY), 0.0, &found);
        if (found) mats.Mats[0].mu = val;
    }

    NFDE::PECRegions_t parser_t::readPECRegions() {
        return buildPECPMCRegions(jlbl::J_MAT_TYPE_PEC);
    }

    NFDE::PECRegions_t parser_t::readPMCRegions() {
        return buildPECPMCRegions(jlbl::J_MAT_TYPE_PMC);
    }

    NFDE::PECRegions_t parser_t::buildPECPMCRegions(const std::string& matType) {
        NFDE::PECRegions_t res;
        std::vector<materialAssociation_t> mAs = getMaterialAssociations(
            {matType}, 
            {std::string("-") + jlbl::J_CONF_SUBTYPE_SURFACE, std::string(jlbl::J_ELEM_TYPE_CELL) + "    ", std::string("-") + jlbl::J_CONF_SUBTYPE_VOLUME + " "}
        );
        
        if (mAs.empty()) { 
            std::vector<NFDE::coords_t> emptyCoords;
            // appendRegion logic simplified for empty case
            res.nLins = 0; res.nLins_max = 0;
            res.nSurfs = 0; res.nSurfs_max = 0;
            res.nVols = 0; res.nVols_max = 0;
            return res;
        }
      
        for (size_t i = 0; i < mAs.size(); ++i) {
            std::vector<NFDE::coords_t> cs;
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_LINEL);
            // appendRegion(res.Lins, res.nLins, res.nLins_max, cs);
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_SURFEL);
            // appendRegion(res.Surfs, res.nSurfs, res.nSurfs_max, cs);
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_VOXEL);
            // appendRegion(res.Vols, res.nVols, res.nVols_max, cs);
            // deallocate(cs);
        }
        return res;
    }

    void parser_t::matAssToCoords(const materialAssociation_t& mA, std::vector<NFDE::coords_t>& res, int cellType) {
        std::vector<NFDE::coords_t> newCoords;
        int nCs = 0;
        
        // Precount
        for (size_t e = 0; e < mA.elementIds.size(); ++e) {
            bool cRf; Cell::cell_region_t cR = mesh.getCellRegion(int(mA.elementIds[e]), cRf);
            // newCoords = cellRegionToCoords(cR, cellType);
            // nCs += newCoords.size();
        }

        // Fills coords
        int jIni = 0;
        res.resize(nCs);
        for (size_t e = 0; e < mA.elementIds.size(); ++e) {
            bool cRf; Cell::cell_region_t cR = mesh.getCellRegion(int(mA.elementIds[e]), cRf);
            std::string tagName = buildTagName(mA.materialId, mA.elementIds[e]);
            // newCoords = cellRegionToCoords(cR, cellType, tag=tagName);
            if (newCoords.empty()) continue;
            int jEnd = jIni + newCoords.size() - 1;
            // res[jIni:jEnd] = newCoords;
            jIni = jEnd + 1; 
        }
    }

    NFDE::ConformalPECRegions_t parser_t::readConformalRegions() {
        NFDE::ConformalPECRegions_t res;
        std::vector<materialAssociation_t> mAs = getMaterialAssociations(
            {jlbl::J_MAT_TYPE_PEC}, 
            {jlbl::J_CONF_SUBTYPE_VOLUME, jlbl::J_CONF_SUBTYPE_SURFACE}
        );

        for (size_t i = 0; i < mAs.size(); ++i) {
            for (size_t j = 0; j < mAs[i].elementIds.size(); ++j) {
                bool found = false;
                Mesh::conformal_region_t cR = mesh.getConformalRegion(mAs[i].elementIds[j], found);
                if (found) { 
                    std::string tagName = buildTagName(mAs[i].materialId, mAs[i].elementIds[j]);
                    if (cR.type == Mesh::REGION_TYPE_VOLUME) { 
                        // appendRegion(res.volumes, cR, tagName);
                    }
                    if (cR.type == Mesh::REGION_TYPE_SURFACE) { 
                        // appendRegion(res.surfaces, cR, tagName);
                    }
                }
            }
        }
        return res;
    }

    NFDE::DielectricRegions_t parser_t::readDielectricRegions() {
        NFDE::DielectricRegions_t res;
        
        fillDielectricsOfCellType(res.Vols, Cell::CELL_TYPE_VOXEL);
        fillDielectricsOfCellType(res.Surfs, Cell::CELL_TYPE_SURFEL);
        fillDielectricsOfCellType(res.Lins, Cell::CELL_TYPE_LINEL);
        
        res.nVols = res.Vols.size();
        res.nSurfs = res.Surfs.size();
        res.nLins = res.Lins.size();

        res.nVols_max = res.nVols;
        res.nSurfs_max = res.nSurfs;
        res.nLins_max = res.nLins;
        return res;
    }

    void parser_t::fillDielectricsOfCellType(std::vector<NFDE::Dielectric_t>& res, int cellType) {
        std::vector<materialAssociation_t> mAs = getMaterialAssociations(
            {std::string(jlbl::J_MAT_TYPE_ISOTROPIC), std::string(jlbl::J_MAT_TYPE_LUMPED)}, {}
        );
        
        if (mAs.empty()) {
            res.clear();
            return;
        }

        int nDielectrics = 0;
        for (size_t i = 0; i < mAs.size(); ++i) {           
            if (containsCellRegionsWithType(mAs[i], cellType)) {
                nDielectrics++;
            } 
        }

        res.resize(nDielectrics);
        
        if (nDielectrics == 0) return;

        int j = 0;
        mAs = getMaterialAssociations({jlbl::J_MAT_TYPE_ISOTROPIC}, {});
        for (size_t i = 0; i < mAs.size(); ++i) {       
            if (!containsCellRegionsWithType(mAs[i], cellType)) continue;
            res[j] = readDielectric(mAs[i], cellType);
            j++;
        }

        mAs = getMaterialAssociations({jlbl::J_MAT_TYPE_LUMPED}, {});
        for (size_t i = 0; i < mAs.size(); ++i) {
            if (!containsCellRegionsWithType(mAs[i], cellType)) continue;
            res[j] = readLumped(mAs[i], cellType);
            j++;
        }
    }

    NFDE::Dielectric_t parser_t::readDielectric(const materialAssociation_t& mA, int cellType) {
        NFDE::Dielectric_t res;
        res.c1P.clear();
        res.n_C1P = 0;
        std::vector<NFDE::coords_t> c2p;
        matAssToCoords(mA, c2p, cellType);
        res.n_C2P = c2p.size();
        
        IdTable::json_value_ptr_t matPtr = matTable.getId(mA.materialId);
        res.sigma  = getRealAt(matPtr, jlbl::J_MAT_ELECTRIC_CONDUCTIVITY, 0.0);
        res.sigmam = getRealAt(matPtr, jlbl::J_MAT_MAGNETIC_CONDUCTIVITY, 0.0);
        res.eps    = getRealAt(matPtr, jlbl::J_MAT_REL_PERMITTIVITY, 1.0) * NFDE::EPSILON_VACUUM;
        res.mu     = getRealAt(matPtr, jlbl::J_MAT_REL_PERMEABILITY, 1.0) * NFDE::MU_VACUUM;
        return res;
    }

    NFDE::Dielectric_t parser_t::readLumped(const materialAssociation_t& mA, int cellType) {
        NFDE::Dielectric_t res;
        res.c1P.clear();
        res.n_C1P = 0;
        std::vector<NFDE::coords_t> c2p;
        matAssToCoords(mA, c2p, cellType);
        res.n_C2P = c2p.size();
        
        IdTable::json_value_ptr_t matPtr = matTable.getId(mA.materialId);
        
        std::string model = getStrAt(matPtr, jlbl::J_MAT_LUMPED_MODEL);
        if (model.empty()) {
            Report::WarnErrReport("ERROR reading lumped material: " + std::to_string(mA.materialId) + " model not found.", true);
        }
        
        res.orient = 1;
        res.DiodOri = 1;
        res.eps = NFDE::EPSILON_VACUUM;
        res.mu = NFDE::MU_VACUUM;

        if (model == jlbl::J_MAT_LUMPED_MODEL_RESISTOR) {
            res.resistor = true;
            res.R = getRealAt(matPtr, jlbl::J_MAT_LUMPED_RESISTANCE, 0.0);
            res.Rtime_on = getRealAt(matPtr, jlbl::J_MAT_LUMPED_STARTING_TIME, 0.0);
            res.Rtime_off = getRealAt(matPtr, jlbl::J_MAT_LUMPED_END_TIME, 1.0);
            if (res.Rtime_on < 0 || res.Rtime_off < 0) {
                Report::WarnErrReport("ERROR reading lumped material: starting || end time is negative", true);
            }
        } else if (model == jlbl::J_MAT_LUMPED_MODEL_INDUCTOR) {
            res.inductor = true;
            res.L = getRealAt(matPtr, jlbl::J_MAT_LUMPED_INDUCTANCE, 0.0);
            res.R = getRealAt(matPtr, jlbl::J_MAT_LUMPED_RESISTANCE, 0.0);
        } else if (model == jlbl::J_MAT_LUMPED_MODEL_CAPACITOR) {
            res.capacitor = true;
            res.C = getRealAt(matPtr, jlbl::J_MAT_LUMPED_CAPACITANCE, 0.0);
            res.R = getRealAt(matPtr, jlbl::J_MAT_LUMPED_RESISTANCE, 0.0);
        } else {
            Report::WarnErrReport("ERROR reading lumped material: invalid model.", true);
        }
        return res;
    }

    bool parser_t::containsCellRegionsWithType(const materialAssociation_t& mA, int cellType) {
        for (size_t e = 0; e < mA.elementIds.size(); ++e) {
            bool cRf; Cell::cell_region_t cR = mesh.getCellRegion(int(mA.elementIds[e]), cRf);
            // if (cellRegionToCoords(cR, cellType).size() != 0) return true;
        }
        return false;
    }

    NFDE::LossyThinSurfaces_t parser_t::readLossyThinSurfaces() {
        NFDE::LossyThinSurfaces_t res;
        std::vector<materialAssociation_t> mAs = getMaterialAssociations({jlbl::J_MAT_TYPE_MULTILAYERED_SURFACE}, {});
        
        int nLossySurfaces = 0;
        for (size_t i = 0; i < mAs.size(); ++i) {
            std::vector<NFDE::coords_t> cs;
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_SURFEL);
            if (!cs.empty()) nLossySurfaces++;
        }

        if (nLossySurfaces == 0) {
            return emptyLossyThinSurfaces();
        }

        res.cs.resize(nLossySurfaces);
        res.length = nLossySurfaces;
        res.length_max = nLossySurfaces;

        int k = 0;
        for (size_t i = 0; i < mAs.size(); ++i) {
            std::vector<NFDE::coords_t> cs;
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_SURFEL);
            if (cs.empty()) continue;
            res.cs[k] = readLossyThinSurface(mAs[i]);
            k++;
        }

        for (size_t i = 0; i < nLossySurfaces; ++i) {
            if (res.nC_max < res.cs[i].c.size()) {
                res.nC_max = res.cs[i].c.size();
            }
        }
        return res;
    }

    NFDE::LossyThinSurface_t parser_t::readLossyThinSurface(const materialAssociation_t& mA) {
        NFDE::LossyThinSurface_t res;
        matAssToCoords(mA, res.c, Cell::CELL_TYPE_SURFEL);
        res.nc = res.c.size();

        IdTable::json_value_ptr_t mat = matTable.getId(mA.materialId);
        res.files = getStrAt(mat, jlbl::J_NAME, " ");
        
        jmod::json_value* layers = nullptr;
        bool tmpFound; core->get(mat, std::string(jlbl::J_MAT_MULTILAYERED_SURF_LAYERS), layers, tmpFound);

        res.numcapas = core->count(layers);
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

        for (int i = 0; i < res.numcapas; ++i) {
            jmod::json_value* layer = nullptr;
            core->get_child(layers, i + 1, layer);
            res.sigma[i] = getRealAt(layer, jlbl::J_MAT_ELECTRIC_CONDUCTIVITY, 0.0);
            res.sigmam[i] = getRealAt(layer, jlbl::J_MAT_MAGNETIC_CONDUCTIVITY, 0.0);
            bool hasAbsPermittivity = false;
            res.eps[i] = getRealAt(layer, jlbl::J_MAT_ABS_PERMITTIVITY, hasAbsPermittivity);
            if (!hasAbsPermittivity) {
                res.eps[i] = getRealAt(layer, jlbl::J_MAT_REL_PERMITTIVITY, 1.0) * NFDE::EPSILON_VACUUM;
            }
            bool hasAbsPermeability = false;
            res.mu[i] = getRealAt(layer, jlbl::J_MAT_ABS_PERMEABILITY, hasAbsPermeability);
            if (!hasAbsPermeability) {
                res.mu[i] = getRealAt(layer, jlbl::J_MAT_REL_PERMEABILITY, 1.0) * NFDE::MU_VACUUM;
            }
            bool found = false;
            res.thk[i] = getRealAt(layer, jlbl::J_MAT_MULTILAYERED_SURF_THICKNESS, found);
            if (!found) {
                Report::WarnErrReport("ERROR reading lossy thin surface: jlbl::J_MAT_MULTILAYERED_SURF_THICKNESS in layer not found.", true);
            }
            res.sigma_devia[i] = 0.0;
            res.eps_devia[i] = 0.0;
            res.mu_devia[i] = 0.0;
            res.sigmam_devia[i] = 0.0;
            res.thk_devia[i] = 0.0;
        }
        return res;
    }

    NFDE::LossyThinSurfaces_t parser_t::emptyLossyThinSurfaces() {
        NFDE::LossyThinSurfaces_t res;
        res.cs.clear();
        res.length = 0;
        res.length_max = 0;
        res.nC_max = 0;
        return res;
    }

    NFDE::PlaneWaves_t parser_t::readPlanewaves() {
        NFDE::PlaneWaves_t res;
        jmod::json_value* sources = nullptr;
        bool found = false;

        core->get(root, jlbl::J_SOURCES, sources, found);
        
        if (!found) {
            res.collection.clear();
            res.nc = 0;
            res.nC_max = 0;
            return res;
        }

        std::vector<IdTable::json_value_ptr_t> pws = jsonValueFilterByKeyValue(sources, jlbl::J_TYPE, jlbl::J_SRC_TYPE_PW);

        res.collection.resize(pws.size());
        for (size_t i = 0; i < pws.size(); ++i) {
            res.collection[i] = readPlanewave(pws[i]);
        }
        res.nc = pws.size();
        res.nC_max = pws.size();

        return res;
    }

    NFDE::PlaneWave_t parser_t::readPlanewave(const jmod::json_value* pw) {
        NFDE::PlaneWave_t res;
        res.nombre_fichero = getStrAt(pw, jlbl::J_SRC_MAGNITUDE_FILE);
        res.atributo = "LOCKED";
        res.theta = getRealAt(pw, std::string(jlbl::J_SRC_PW_DIRECTION) + "." + std::string(jlbl::J_SRC_PW_THETA), 0.0);
        res.phi = getRealAt(pw, std::string(jlbl::J_SRC_PW_DIRECTION) + "." + std::string(jlbl::J_SRC_PW_PHI), 0.0);
        res.alpha = getRealAt(pw, std::string(jlbl::J_SRC_PW_POLARIZATION) + "." + std::string(jlbl::J_SRC_PW_THETA), 0.0);
        res.beta = getRealAt(pw, std::string(jlbl::J_SRC_PW_POLARIZATION) + "." + std::string(jlbl::J_SRC_PW_PHI), 0.0);

        {
            std::vector<Cell::cell_interval_t> cellIntervals = getSingleVolumeInElementsIds(pw);
            if (cellIntervals.empty()) return res;
            std::vector<NFDE::coords_t> nfdeCoords = Pt::cellIntervalsToCoords(cellIntervals);
            if (!nfdeCoords.empty()) {
                res.coor1[0] = nfdeCoords[0].Xi;
                res.coor1[1] = nfdeCoords[0].Yi;
                res.coor1[2] = nfdeCoords[0].Zi;
                res.coor2[0] = nfdeCoords[0].Xe;
                res.coor2[1] = nfdeCoords[0].Ye;
                res.coor2[2] = nfdeCoords[0].Ze;
            }
        }

        res.isRC = false;
        res.numModes = 1;
        res.INCERTMAX = 0.0;
        return res;
    }

    NFDE::NodSource_t parser_t::readNodalSources() {
        NFDE::NodSource_t res;
        jmod::json_value* sources = nullptr;
        bool found = false;

        core->get(root, jlbl::J_SOURCES, sources, found);
        if (!found) {
            res.NodalSource.clear();
            return res;
        }

        std::vector<IdTable::json_value_ptr_t> nodSrcs = jsonValueFilterByKeyValues(sources, jlbl::J_TYPE, {jlbl::J_SRC_TYPE_NS});
        if (nodSrcs.empty()) {
            res.NodalSource.clear();
            return res;
        }

        res.NodalSource.resize(nodSrcs.size());
        res.n_nodSrc = nodSrcs.size();
        res.n_nodSrc_max = nodSrcs.size();
        for (size_t i = 0; i < nodSrcs.size(); ++i) {
            res.NodalSource[i] = readField(nodSrcs[i]);
        }
        for (size_t i = 0; i < nodSrcs.size(); ++i) {
            res.n_C1P_max = std::max(res.n_C1P_max, res.NodalSource[i].n_C1P);
            res.n_C2P_max = std::max(res.n_C2P_max, res.NodalSource[i].n_C2P);
        }
        if (res.n_nodSrc > 0) {
            if (res.n_C1P_max == 0) res.n_C1P_max = 1;
            if (res.n_C2P_max == 0) res.n_C2P_max = 1;
        }
        return res;
    }

    NFDE::Curr_Field_Src_t parser_t::readField(const jmod::json_value* jns) {
        NFDE::Curr_Field_Src_t res;
        std::string field = getStrAt(jns, jlbl::J_FIELD, jlbl::J_FIELD_CURRENT);
        
        if (field == jlbl::J_FIELD_CURRENT) {
            res.isElec = true;
        } else {
            Report::WarnErrReport("Error reading current field source. Field label not recognized.", false);
        }
        
        std::string hardness = getStrAt(jns, jlbl::J_SRC_NS_HARDNESS, jlbl::J_SRC_NS_HARDNESS_SOFT);
        if (hardness == jlbl::J_SRC_NS_HARDNESS_SOFT) {
            res.isHard = false;
        } else if (hardness == jlbl::J_SRC_NS_HARDNESS_HARD) {
            res.isHard = true;
        } else {
            Report::WarnErrReport("Error reading current field source. Hardness label not recognized.", false);
        }
        
        res.isInitialValue = false;
        res.nombre = getStrAt(jns, jlbl::J_SRC_MAGNITUDE_FILE);

        std::vector<int> elementIds = getIntsAt(jns, jlbl::J_ELEMENTIDS);
        std::vector<Cell::cell_region_t> cRs = mesh.getCellRegions(elementIds);
        std::vector<Cell::cell_interval_t> intervals = Pt::getIntervalsInCellRegions(cRs, Cell::CELL_TYPE_LINEL);
        std::vector<NFDE::coords_t> cs = Pt::cellIntervalsToCoords(intervals);
        std::vector<NFDE::coords_scaled_t> scaledCoords = Pt::coordsToScaledCoords(cs);

        int cnt_c1p = 0, cnt_c2p = 0;
        for (const auto& c : scaledCoords) {
            if (c.Xi == c.Xe && c.Yi == c.Ye && c.Zi == c.Ze) cnt_c1p++;
            else cnt_c2p++;
        }
        if (cnt_c1p > 0) res.c1P.resize(cnt_c1p);
        if (cnt_c2p > 0) res.c2P.resize(cnt_c2p);
        cnt_c1p = 0; cnt_c2p = 0;
        for (const auto& c : scaledCoords) {
            if (c.Xi == c.Xe && c.Yi == c.Ye && c.Zi == c.Ze) {
                res.c1P[cnt_c1p++] = c;
            } else {
                res.c2P[cnt_c2p++] = c;
            }
        }
        res.n_C1P = cnt_c1p;
        res.n_C2P = cnt_c2p;

        return res;
    }

    NFDE::Sondas_t parser_t::readProbes() {
        NFDE::Sondas_t res;
        jmod::json_value* allProbes = nullptr;
        bool found = false;

        core->get(root, jlbl::J_PROBES, allProbes, found);
        if (!found) {
            res.probes.clear();
            res.n_probes = 0;
            res.n_probes_max = 0;
            return res;
        }

        std::vector<std::string> validTypes = {jlbl::J_PR_TYPE_FARFIELD};
        std::vector<IdTable::json_value_ptr_t> ps = jsonValueFilterByKeyValues(allProbes, jlbl::J_TYPE, validTypes);

        res.n_probes = ps.size();
        res.n_probes_max = ps.size();
        res.probes.resize(ps.size());
        for (size_t i = 0; i < ps.size(); ++i) {
            res.probes[i] = readFarFieldProbe(ps[i]);
        }
        return res;
    }

    NFDE::abstractSonda_t parser_t::readFarFieldProbe(const jmod::json_value* p) {
        NFDE::abstractSonda_t res;
        res.n_FarField = 1;
        res.n_FarField_max = 1;
        res.FarField.resize(1);
        
        NFDE::Sonda_t* ff = &res.FarField[0].probe;
        ff->grname = " ";
        ff->outputrequest = getStrAt(p, jlbl::J_NAME);

        domain_t domain = getDomain(p, jlbl::J_PR_DOMAIN);
        if (domain.type2 != NFDE::NP_T2_FREQ) {
            Report::WarnErrReport("Only frequency domain is accepted for far field probes.", false);
        }
        ff->tstart = 0.0;
        ff->tstop = 0.0;
        ff->tstep = 0.0;
        ff->fstart = domain.fstart;
        ff->fstop = domain.fstop;
        ff->fstep = domain.fstep;

        {
            bool transferFunctionFound = false;
            std::string fn = getStrAt(p, std::string(jlbl::J_PR_DOMAIN) + jlbl::J_PR_DOMAIN_MAGNITUDE_FILE, "");
            if (!transferFunctionFound) {
                jmod::json_value* sources = nullptr;
                bool sourcesFound = false;
                core->get(root, jlbl::J_SOURCES, sources, sourcesFound);
                if (sourcesFound) {
                    if (core->count(sources) == 1) {
                        jmod::json_value* src = nullptr;
                        core->get_child(sources, 1, src);
                        fn = getStrAt(src, jlbl::J_SRC_MAGNITUDE_FILE, "");
                    }
                }
            }

            if (transferFunctionFound) {
                ff->FileNormalize = fn;
            } else {
                ff->FileNormalize = " ";
            }
        }

        if (domain.isLogarithmicFrequencySpacing) {
            appendLogSufix(ff->outputrequest);
        }

        {
            // nfdeCoords = Pt::cellIntervalsToCoords(getSingleVolumeInElementsIds(p));
            // ff->n_cord = 2;
            // ff->i.resize(2);
            // ff->j.resize(2);
            // ff->k.resize(2);
            // ff->node.clear();
            // ff->i[0] = nfdeCoords[0].Xi;
            // ff->i[1] = nfdeCoords[0].Xe;
            // ff->j[0] = nfdeCoords[0].Yi;
            // ff->j[1] = nfdeCoords[0].Ye;
            // ff->k[0] = nfdeCoords[0].Zi;
            // ff->k[1] = nfdeCoords[0].Ze;
        }

        {
            // readDirection(p, jlbl::J_PR_FAR_FIELD_PHI, ff->phistart, ff->phistop, ff->phistep);
            // readDirection(p, jlbl::J_PR_FAR_FIELD_THETA, ff->thetastart, ff->thetastop, ff->thetastep);
        }
        return res;
    }

    void parser_t::readDirection(const jmod::json_value* p, const std::string& label, double& initial, double& final, double& step) {
        jmod::json_value* dir = nullptr;
        bool found = false;
        core->get(p, label, dir, found);
        if (!found) {
            Report::WarnErrReport("Error reading far field probe. Direction label not found.", true);
        }
        initial = getRealAt(dir, jlbl::J_PR_FAR_FIELD_DIR_INITIAL, 0.0);
        final = getRealAt(dir, jlbl::J_PR_FAR_FIELD_DIR_FINAL, 0.0);
        step = getRealAt(dir, jlbl::J_PR_FAR_FIELD_DIR_STEP, 0.0);
    }

    NFDE::MasSondas_t parser_t::readMoreProbes() {
        NFDE::MasSondas_t res;
        jmod::json_value* allProbes = nullptr;
        bool found = false;

        core->get(root, jlbl::J_PROBES, allProbes, found);
        if (!found) {
            res.collection.clear();
            res.length = 0;
            res.length_max = 0;
            res.len_cor_max = 0;
            return res;
        }

        std::vector<std::string> validTypes = {jlbl::J_PR_TYPE_POINT, jlbl::J_PR_TYPE_LINE};
        std::vector<IdTable::json_value_ptr_t> ps = jsonValueFilterByKeyValues(allProbes, jlbl::J_TYPE, validTypes);
        
        int filtered_size = 0;
        for (size_t i = 0; i < ps.size(); ++i) {
            if (isMoreProbe(ps[i])) { 
                filtered_size++;
            }
        }

        int n = 0;
        res.collection.resize(filtered_size);
        for (size_t i = 0; i < ps.size(); ++i) {
            if (isMoreProbe(ps[i])) { 
                std::string probeLbl = getStrAt(ps[i], jlbl::J_TYPE, jlbl::J_FIELD_ELECTRIC);
                if (probeLbl == jlbl::J_PR_TYPE_POINT) { 
                    res.collection[n] = readPointProbe(ps[i]);
                    n++;
                } else if (probeLbl == jlbl::J_PR_TYPE_LINE) {
                    res.collection[n] = readLineProbe(ps[i]);
                    n++;
                }
            }
        }

        res.length = res.collection.size();
        res.length_max = res.collection.size();
        for (size_t i = 0; i < res.collection.size(); ++i) {
            if (res.collection[i].cordinates.size() > res.len_cor_max) {
                res.len_cor_max = res.collection[i].cordinates.size();
            }
        }
        return res;
    }

    bool parser_t::isMoreProbe(const jmod::json_value* p) {
        return isPointProbe(p) || isLineProbe(p);
    }

    bool parser_t::isLineProbe(const jmod::json_value* p) {
        return getStrAt(p, jlbl::J_TYPE) == jlbl::J_PR_TYPE_LINE;
    }

    bool parser_t::isPointProbe(const jmod::json_value* p) {
        bool found = false;
        std::string typeLabel = getStrAt(p, jlbl::J_TYPE, "", &found);
        if (!found) {
            return false;
        }
        if (typeLabel != jlbl::J_PR_TYPE_POINT) {
            return false;
        }

        std::string fieldLabel = getStrAt(p, jlbl::J_FIELD, jlbl::J_FIELD_ELECTRIC);
        return (fieldLabel == jlbl::J_FIELD_ELECTRIC || fieldLabel == jlbl::J_FIELD_MAGNETIC);
    }

    NFDE::MasSonda_t parser_t::readLineProbe(const jmod::json_value* p) {
        NFDE::MasSonda_t res;
        bool nameFound = false;
        std::string outputName = getStrAt(p, jlbl::J_NAME, "", &nameFound);
        if (!nameFound) {
            Report::WarnErrReport("ERROR: name entry not found for probe.", false);
        }
        res.outputrequest = outputName;

        setDomain(res, getDomain(p, jlbl::J_PR_DOMAIN));

        bool elementIdsFound = false;
        std::vector<int> elemIds = getIntsAt(p, jlbl::J_ELEMENTIDS, &elementIdsFound);
        if (!elementIdsFound) {
            Report::WarnErrReport("ERROR: element ids entry not found for probe.", false);
        }
        if (elemIds.size() != 1) {
            Report::WarnErrReport("ERROR: point probe must contain a single element id.", true);
        }

        bool pf; Mesh::polyline_t polyline = mesh.getPolyline(int(elemIds[0]), pf);
        std::vector<Mesh::linel_t> linels = mesh.polylineToLinels(polyline);
        res.cordinates.resize(linels.size());
        for (size_t i = 0; i < linels.size(); ++i) {
            res.cordinates[i].Xi = linels[i].cell[0];
            res.cordinates[i].Yi = linels[i].cell[1];
            res.cordinates[i].Zi = linels[i].cell[2];
            int orientation = std::abs(linels[i].orientation);
            if (orientation == 1) res.cordinates[i].Xe = linels[i].cell[0] + 1;
            else if (orientation == 2) res.cordinates[i].Ye = linels[i].cell[1] + 1;
            else if (orientation == 3) res.cordinates[i].Ze = linels[i].cell[2] + 1;
            res.cordinates[i].Or = ((linels[i].orientation) > 0 ? 1 : -1) * NFDE::NP_COR_LINE;
            res.cordinates[i].tag = outputName;
        }

        res.len_cor = 1;
        return res;
    }

    NFDE::MasSonda_t parser_t::readPointProbe(const jmod::json_value* p) {
        NFDE::MasSonda_t res;
        bool nameFound = false;
        std::string outputName = getStrAt(p, jlbl::J_NAME, "", &nameFound);
        if (!nameFound) {
            Report::WarnErrReport("Point probes must define a name.", false);
        }
        res.outputrequest = outputName;

        setDomain(res, getDomain(p, jlbl::J_PR_DOMAIN));

        bool elementIdsFound = false;
        std::vector<int> elemIds = getIntsAt(p, jlbl::J_ELEMENTIDS, &elementIdsFound);
        if (!elementIdsFound) {
            Report::WarnErrReport("Element ids entry not found for probe.", false);
        }
        if (elemIds.size() != 1) {
            Report::WarnErrReport("Point probe must contain a single element id.", true);
        }

        Mesh::pixel_t pixel = Pt::getPixelFromElementId(mesh, elemIds[0]);
        
        bool typeLabelFound = false;
        std::string typeLabel = getStrAt(p, jlbl::J_TYPE, "", &typeLabelFound);
        if (!typeLabelFound) {
            Report::WarnErrReport("Point probe type label not found.", true);
        }
        
        if (typeLabel == jlbl::J_PR_TYPE_POINT) {
            jmod::json_value* dirLabelPtr = nullptr;
            bool dirLabelsFound = false;
            core->get(p, jlbl::J_PR_POINT_DIRECTIONS, dirLabelPtr, dirLabelsFound);
            std::vector<char> dirLabels;
            if (dirLabelsFound) {
                dirLabels = buildDirLabels(dirLabelPtr);
            } else {
                dirLabels = {jlbl::J_DIR_X[0], jlbl::J_DIR_Y[0], jlbl::J_DIR_Z[0]};
            }
            
            bool fieldLabelFound = false;
            std::string fieldLabel = getStrAt(p, jlbl::J_FIELD, jlbl::J_FIELD_ELECTRIC);
            res.cordinates.resize(dirLabels.size());
            for (size_t j = 0; j < dirLabels.size(); ++j) {
                res.cordinates[j].tag = outputName;
                res.cordinates[j].Xi = static_cast<int>(pixel.cell[0]);
                res.cordinates[j].Yi = static_cast<int>(pixel.cell[1]);
                res.cordinates[j].Zi = static_cast<int>(pixel.cell[2]);
                res.cordinates[j].Or = strToFieldType(fieldLabel, dirLabels[j]);
            }
        }

        res.len_cor = res.cordinates.size();
        return res;
    }

    std::vector<char> parser_t::buildDirLabels(const jmod::json_value* dirLabelsPtr) {
        std::vector<char> res(core->count(dirLabelsPtr));
        for (int i = 0; i < core->count(dirLabelsPtr); ++i) {
            jmod::json_value* child = nullptr;
            core->get_child(dirLabelsPtr, i + 1, child);
            std::string str = getStrAt(child, "");
            res[i] = str[0];
        }
        return res;
    }

    void parser_t::setDomain(NFDE::MasSonda_t& res, const domain_t& domain) {
        res.tstart = domain.tstart;
        res.tstep = domain.tstep;
        res.tstop = domain.tstop;
        res.fstart = domain.fstart;
        res.fstep = domain.fstep;
        res.fstop = domain.fstop;
        if (!domain.filename.empty()) {
            res.filename = domain.filename;
        } else {
            res.filename = " ";
        }
        res.type1 = domain.type1;
        res.type2 = domain.type2;

        if (domain.isLogarithmicFrequencySpacing) {
            appendLogSufix(res.outputrequest);
        }
    }

    int parser_t::strToFieldType(const std::string& fieldLabel, char dirLabel) {
        if (fieldLabel == jlbl::J_FIELD_ELECTRIC) {
            if (dirLabel == jlbl::J_DIR_X[0]) return NFDE::NP_COR_EX;
            if (dirLabel == jlbl::J_DIR_Y[0]) return NFDE::NP_COR_EY;
            if (dirLabel == jlbl::J_DIR_Z[0]) return NFDE::NP_COR_EZ;
            Report::WarnErrReport("Invalid dir label for electric field probe.", true);
        } else if (fieldLabel == jlbl::J_FIELD_MAGNETIC) {
            if (dirLabel == jlbl::J_DIR_X[0]) return NFDE::NP_COR_HX;
            if (dirLabel == jlbl::J_DIR_Y[0]) return NFDE::NP_COR_HY;
            if (dirLabel == jlbl::J_DIR_Z[0]) return NFDE::NP_COR_HZ;
            Report::WarnErrReport("Invalid dir label for magnetic field probe.", true);
        } else if (fieldLabel == jlbl::J_FIELD_CURRENT) {
            return NFDE::NP_COR_WIRECURRENT;
        } else if (fieldLabel == jlbl::J_FIELD_VOLTAGE) {
            return NFDE::NP_COR_DDP;
        } else if (fieldLabel == jlbl::J_FIELD_CHARGE) {
            return NFDE::NP_COR_CHARGE;
        } else {
            Report::WarnErrReport("Invalid field label for point/wire probe.", true);
        }
        return 0;
    }

    NFDE::BloqueProbes_t parser_t::readBlockProbes() {
        NFDE::BloqueProbes_t res;
        std::vector<IdTable::json_value_ptr_t> bps;
        jmod::json_value* probes = nullptr;
        bool found = false;

        core->get(root, jlbl::J_PROBES, probes, found);
        if (!found) {
            res.bp.clear();
            return res;
        }

        bps = jsonValueFilterByKeyValues(probes, jlbl::J_TYPE, {jlbl::J_PR_TYPE_BULK_CURRENT});
        if (bps.empty()) {
            res.bp.clear();
            return res;
        }

        res.n_bp = bps.size();
        res.n_bp_max = bps.size();
        res.bp.resize(bps.size());
        for (size_t i = 0; i < bps.size(); ++i) {
            res.bp[i] = readBlockProbe(bps[i]);
        }
        return res;
    }

    NFDE::BloqueProbe_t parser_t::readBlockProbe(const jmod::json_value* bp) {
        NFDE::BloqueProbe_t res;
        std::vector<int> elemIds = getIntsAt(bp, jlbl::J_ELEMENTIDS);
        std::vector<Cell::cell_region_t> cRs = mesh.getCellRegions(elemIds);
        
        if (cRs.size() != 1) {
            Report::WarnErrReport("Bulk current probe must be defined by a single cell region.", false);
        }

        if (cRs[0].intervals.size() != 1) {
            Report::WarnErrReport("Bulk current probe must be defined by a single cell interval.", false);
        }
        
        std::vector<NFDE::coords_t> cs = Pt::cellIntervalsToCoords(cRs[0].intervals);

        res.i1 = cs[0].Xi;
        res.i2 = cs[0].Xe;
        res.j1 = cs[0].Yi;
        res.j2 = cs[0].Ye;
        res.k1 = cs[0].Zi;
        res.k2 = cs[0].Ze;
        res.nml = std::abs(cs[0].Or);
        
        if (res.nml == 0) {
            std::string direction = getStrAt(bp, jlbl::J_DIR);
            if (direction == jlbl::J_DIR_X) res.nml = 1;
            else if (direction == jlbl::J_DIR_Y) res.nml = 2;
            else if (direction == jlbl::J_DIR_Z) res.nml = 3;
            else Report::WarnErrReport("Null direction detected for bulk probe. Check definition", true);
        }

        res.outputrequest = getStrAt(bp, jlbl::J_NAME);
        setDomainBlock(res, getDomain(bp, jlbl::J_PR_DOMAIN));

        res.skip = 1;
        res.tag = getStrAt(bp, jlbl::J_NAME, " ");
        res.t = NFDE::BcELECT;
        return res;
    }

    void parser_t::setDomainBlock(NFDE::BloqueProbe_t& res, const domain_t& domain) {
        res.tstart = domain.tstart;
        res.tstep = domain.tstep;
        res.tstop = domain.tstop;
        res.fstart = domain.fstart;
        res.fstep = domain.fstep;
        res.fstop = domain.fstop;
        if (!domain.filename.empty()) {
            res.FileNormalize = domain.filename;
        } else {
            res.FileNormalize = " ";
        }
        res.type2 = domain.type2;

        if (domain.isLogarithmicFrequencySpacing) {
            appendLogSufix(res.outputrequest);
        }
    }

    NFDE::VolProbes_t parser_t::readVolumicProbes() {
        NFDE::VolProbes_t res;
        std::vector<IdTable::json_value_ptr_t> ps;
        jmod::json_value* probes = nullptr;
        bool found = false;

        core->get(root, jlbl::J_PROBES, probes, found);
        if (!found) {
            return buildNoVolProbes();
        }

        ps = jsonValueFilterByKeyValues(probes, jlbl::J_TYPE, {jlbl::J_PR_TYPE_MOVIE});
        if (ps.empty()) {
            return buildNoVolProbes();
        }

        res.length = ps.size();
        res.length_max = ps.size();
        res.len_cor_max = 2 * ps.size();
        res.collection.resize(ps.size());
        for (size_t i = 0; i < ps.size(); ++i) {
            res.collection[i] = readVolProbe(ps[i]);
        }
        return res;
    }

    NFDE::VolProbes_t parser_t::buildNoVolProbes() {
        NFDE::VolProbes_t res;
        res.collection.clear();
        res.length = 0;
        res.length_max = 0;
        res.len_cor_max = 0;
        return res;
    }

    NFDE::VolProbe_t parser_t::readVolProbe(const jmod::json_value* p) {
        NFDE::VolProbe_t res;
        bool found = false;
        std::vector<int> elemIds = getIntsAt(p, jlbl::J_ELEMENTIDS, &found);
        std::vector<Cell::cell_region_t> cRs = mesh.getCellRegions(elemIds);
        
        if (cRs.size() != 1) {
            Report::WarnErrReport("Movie probe must be defined over a single cell region.", false);
        }

        if (cRs[0].intervals.size() != 1) {
            Report::WarnErrReport("Movie probe must be defined by a single cell interval.", false);
        }
        
        std::vector<NFDE::coords_t> cs = Pt::cellIntervalsToCoords(cRs[0].intervals);

        std::string fieldType = getStrAt(p, jlbl::J_FIELD, jlbl::J_FIELD_ELECTRIC);
        jmod::json_value* compsPtr = nullptr;
        bool componentsFound = false;
        core->get(p, jlbl::J_PR_MOVIE_COMPONENT, compsPtr, componentsFound);
        
        res.cordinates.resize(1);
        std::string component;
        if (componentsFound) {
            component = getStrAt(compsPtr, "");
            res.cordinates[0] = cs[0];
            res.cordinates[0].Or = buildVolProbeType(fieldType, component);
        } else {
            component = jlbl::J_DIR_M;
            res.cordinates[0].Or = buildVolProbeType(fieldType, component);
        }
        res.len_cor = res.cordinates.size();
        
        res.outputrequest = getStrAt(p, jlbl::J_NAME, " ");
        setDomainVol(res, getDomain(p, jlbl::J_PR_DOMAIN));
        return res;
    }

    int parser_t::buildVolProbeType(const std::string& fieldType, const std::string& component) {
        if (fieldType == jlbl::J_FIELD_ELECTRIC) {
            if (component == jlbl::J_DIR_X) return NFDE::iExC;
            if (component == jlbl::J_DIR_Y) return NFDE::iEyC;
            if (component == jlbl::J_DIR_Z) return NFDE::iEzC;
            if (component == jlbl::J_DIR_M) return NFDE::iMEC;
        } else if (fieldType == jlbl::J_FIELD_MAGNETIC) {
            if (component == jlbl::J_DIR_X) return NFDE::iHxC;
            if (component == jlbl::J_DIR_Y) return NFDE::iHyC;
            if (component == jlbl::J_DIR_Z) return NFDE::iHzC;
            if (component == jlbl::J_DIR_M) return NFDE::iMHC;
        } else if (fieldType == jlbl::J_FIELD_CURRENT_DENSITY) {
            if (component == jlbl::J_DIR_X) return NFDE::iCurX;
            if (component == jlbl::J_DIR_Y) return NFDE::iCurY;
            if (component == jlbl::J_DIR_Z) return NFDE::iCurZ;
            if (component == jlbl::J_DIR_M) return NFDE::iCur;
        } else {
            Report::WarnErrReport("Invalid field type for movie probe.", true);
        }
        return 0;
    }

    void parser_t::setDomainVol(NFDE::VolProbe_t& res, const domain_t& domain) {
        res.tstart = domain.tstart;
        res.tstep = domain.tstep;
        res.tstop = domain.tstop;
        res.fstart = domain.fstart;
        res.fstep = domain.fstep;
        res.fstop = domain.fstop;
        if (!domain.filename.empty()) {
            res.filename = domain.filename;
        } else {
            res.filename = " ";
        }
        res.type2 = domain.type2;

        if (domain.isLogarithmicFrequencySpacing) {
            appendLogSufix(res.outputrequest);
        }
    }

    void parser_t::appendLogSufix(std::string& fn) {
        fn = fn + "_log_";
    }

    NFDE::ThinSlots_t parser_t::readThinSlots() {
        NFDE::ThinSlots_t res;
        std::vector<materialAssociation_t> mAs = getMaterialAssociations({std::string(jlbl::J_MAT_TYPE_SLOT)}, {});
        
        if (mAs.empty()) {
            res.tg.clear();
            return res;
        }

        res.n_tg = mAs.size();
        res.tg.resize(res.n_tg);
        for (size_t i = 0; i < mAs.size(); ++i) {
            res.tg[i] = readThinSlot(mAs[i]);
        }
        return res;
    }

    NFDE::ThinSlot_t parser_t::readThinSlot(const materialAssociation_t& mA) {
        NFDE::ThinSlot_t res;
        IdTable::json_value_ptr_t mat = matTable.getId(mA.materialId);
        bool found = false;
        res.width = getRealAt(mat, jlbl::J_MAT_THINSLOT_WIDTH, found);
        if (!found) {
            Report::WarnErrReport("Missing thin slot width for material " + std::to_string(mA.materialId), true);
        }
        std::vector<NFDE::coords_t> coords;
        matAssToCoords(mA, coords, Cell::CELL_TYPE_LINEL);
        coordsToThinSlotComp(coords, res.tgc);
        return res;
    }

    void parser_t::coordsToThinSlotComp(const std::vector<NFDE::coords_t>& cs,
                                         std::vector<NFDE::ThinSlotComp_t>& tc) {
        tc.resize(cs.size());
        for (size_t i = 0; i < cs.size(); ++i) {
            tc[i] = buildBaseThinSlotComponent(cs[i]);
        }
    }

    NFDE::ThinSlotComp_t parser_t::buildBaseThinSlotComponent(const NFDE::coords_t& cs) {
        NFDE::ThinSlotComp_t res;
        res.i = cs.Xi;
        res.j = cs.Yi;
        res.K = cs.Zi;
        res.Or = cs.Or;
        return res;
    }

} // namespace smbjson
