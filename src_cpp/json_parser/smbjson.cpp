#include "smbjson_m.h"
#include "NFDETypes_extension_m.h"
#include <iostream>
#include <cmath>
#include <fstream>

namespace smbjson {

    namespace {
        bool parseOneBasedArrayIndex(const std::string& token, int& idx) {
            if (token.size() < 3 || token.front() != '(' || token.back() != ')') return false;
            try {
                idx = std::stoi(token.substr(1, token.size() - 2));
            } catch (...) {
                return false;
            }
            return idx >= 1;
        }

        const nlohmann::json* findPathNode(const nlohmann::json* root, const std::string& path) {
            if (root == nullptr) return nullptr;
            if (path.empty()) return root;

            const nlohmann::json* current = root;
            size_t pos = 0;
            while (pos <= path.size()) {
                const size_t dot = path.find('.', pos);
                const std::string token =
                    path.substr(pos, dot == std::string::npos ? std::string::npos : dot - pos);

                if (!token.empty()) {
                    int idx = 0;
                    if (parseOneBasedArrayIndex(token, idx)) {
                        if (!current->is_array() || idx > static_cast<int>(current->size())) {
                            return nullptr;
                        }
                        current = &(*current)[static_cast<size_t>(idx - 1)];
                    } else {
                        if (!current->is_object()) return nullptr;
                        auto it = current->find(token);
                        if (it == current->end()) return nullptr;
                        current = &(*it);
                    }
                }

                if (dot == std::string::npos) break;
                pos = dot + 1;
            }
            return current;
        }

        void jsonGet(const nlohmann::json* val, const std::string& key,
                     const nlohmann::json*& out, bool& found) {
            out = findPathNode(val, key);
            found = (out != nullptr);
        }

        int jsonCount(const nlohmann::json* val) {
            if (val == nullptr) return 0;
            if (val->is_array() || val->is_object()) {
                return static_cast<int>(val->size());
            }
            return 0;
        }

        bool jsonGetChild(const nlohmann::json* val, int oneBasedIndex,
                          const nlohmann::json*& out) {
            out = nullptr;
            if (val == nullptr || !val->is_array()) return false;
            if (oneBasedIndex < 1 || oneBasedIndex > static_cast<int>(val->size())) {
                return false;
            }
            out = &(*val)[static_cast<size_t>(oneBasedIndex - 1)];
            return true;
        }

        void appendRegion(std::vector<NFDE::coords_t>& dest, int& n, int& nMax,
                          const std::vector<NFDE::coords_t>& cs) {
            if (dest.empty()) {
                dest = cs;
                n = static_cast<int>(cs.size());
                nMax = n;
                return;
            }
            std::vector<NFDE::coords_t> aux = dest;
            dest.resize(aux.size() + cs.size());
            for (size_t i = 0; i < aux.size(); ++i) dest[i] = aux[i];
            for (size_t i = 0; i < cs.size(); ++i) dest[aux.size() + i] = cs[i];
            n = static_cast<int>(dest.size());
            nMax = n;
        }
    }

    // ---- JSON accessor implementations ----

    bool parser_t::getLogicalAt(const nlohmann::json* val, const std::string& key, bool default_val, bool* foundOut) {
        const nlohmann::json* ptr = nullptr;
        bool found = false;
        jsonGet(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->get<bool>();
        return default_val;
    }

    int parser_t::getIntAt(const nlohmann::json* val, const std::string& key, int default_val, bool* foundOut) {
        const nlohmann::json* ptr = nullptr;
        bool found = false;
        jsonGet(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->get<int>();
        return default_val;
    }

    std::vector<int> parser_t::getIntsAt(const nlohmann::json* val, const std::string& key, bool* foundOut) {
        std::vector<int> res;
        const nlohmann::json* arr = nullptr;
        bool found = false;
        jsonGet(val, key, arr, found);
        if (foundOut) *foundOut = found;
        if (found && arr) {
            int n = jsonCount(arr);
            res.resize(n);
            for (int i = 0; i < n; ++i) {
                const nlohmann::json* child = nullptr;
                jsonGetChild(arr, i + 1, child);
                if (child) res[i] = child->get<int>();
            }
        }
        return res;
    }

    double parser_t::getRealAt(const nlohmann::json* val, const std::string& key, double default_val, bool* foundOut) {
        const nlohmann::json* ptr = nullptr;
        bool found = false;
        jsonGet(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->get<double>();
        return default_val;
    }

    std::vector<double> parser_t::getRealsAt(const nlohmann::json* val, const std::string& key, bool* foundOut) {
        std::vector<double> res;
        const nlohmann::json* arr = nullptr;
        bool found = false;
        jsonGet(val, key, arr, found);
        if (foundOut) *foundOut = found;
        if (found && arr) {
            int n = jsonCount(arr);
            res.resize(n);
            for (int i = 0; i < n; ++i) {
                const nlohmann::json* child = nullptr;
                jsonGetChild(arr, i + 1, child);
                if (child) res[i] = child->get<double>();
            }
        }
        return res;
    }

    std::vector<std::vector<double>> parser_t::getMatrixAt(const nlohmann::json* val, const std::string& key, bool* foundOut) {
        std::vector<std::vector<double>> res;
        const nlohmann::json* arr = nullptr;
        bool found = false;
        jsonGet(val, key, arr, found);
        if (foundOut) *foundOut = found;
        if (found && arr) {
            int n = jsonCount(arr);
            for (int i = 0; i < n; ++i) {
                const nlohmann::json* child = nullptr;
                jsonGetChild(arr, i + 1, child);
                res.push_back(getRealsAt(child, ""));
            }
        }
        return res;
    }

    std::string parser_t::getStrAt(const nlohmann::json* val, const std::string& key, const std::string& default_val, bool* foundOut) {
        const nlohmann::json* ptr = nullptr;
        bool found = false;
        jsonGet(val, key, ptr, found);
        if (foundOut) *foundOut = found;
        if (found && ptr) return ptr->get<std::string>();
        return default_val;
    }

    parser_t::domain_t parser_t::getDomain(const nlohmann::json* place, const std::string& path) {
        domain_t res;
        const nlohmann::json* domain = nullptr;
        bool found = false;
        jsonGet(place, path, domain, found);
        if (!found) {
            res.filename = " ";
            return res;
        }

        bool transferFunctionFound = false;
        std::string fn = getStrAt(domain, jlbl::J_PR_DOMAIN_MAGNITUDE_FILE, " ", &transferFunctionFound);
        if (transferFunctionFound) {
            while (!fn.empty() && (fn.front() == ' ' || fn.front() == '\t')) fn.erase(fn.begin());
            while (!fn.empty() && (fn.back() == ' ' || fn.back() == '\t')) fn.pop_back();
            res.filename = fn;
            res.hasTransferFunction = true;
        }

        res.type1 = NFDE::NP_T1_PLAIN;
        std::string domainType = getStrAt(domain, jlbl::J_PR_DOMAIN_TYPE, jlbl::J_PR_DOMAIN_TYPE_TIME);
        res.type2 = getNPDomainType(domainType, transferFunctionFound);

        res.tstart = getRealAt(domain, jlbl::J_PR_DOMAIN_TIME_START, 0.0);
        res.tstop = getRealAt(domain, jlbl::J_PR_DOMAIN_TIME_STOP, 0.0);
        res.tstep = getRealAt(domain, jlbl::J_PR_DOMAIN_TIME_STEP, 0.0);
        res.fstart = getRealAt(domain, jlbl::J_PR_DOMAIN_FREQ_START, 0.0);
        res.fstop = getRealAt(domain, jlbl::J_PR_DOMAIN_FREQ_STOP, 0.0);

        int numberOfFrequencies = getIntAt(domain, jlbl::J_PR_DOMAIN_FREQ_NUMBER, 0);
        if (numberOfFrequencies == 0) {
            res.fstep = 0.0;
        } else {
            res.fstep = (res.fstop - res.fstart) / numberOfFrequencies;
        }

        std::string freqSpacing = getStrAt(domain, jlbl::J_PR_DOMAIN_FREQ_SPACING,
                                           jlbl::J_PR_DOMAIN_FREQ_SPACING_LINEAR);
        res.isLogarithmicFrequencySpacing =
            (freqSpacing == jlbl::J_PR_DOMAIN_FREQ_SPACING_LOGARITHMIC);
        return res;
    }

    int parser_t::getNPDomainType(const std::string& typeLabel, bool hasTransferFunction) {
        bool isTime = false;
        bool isFrequency = false;
        if (typeLabel == jlbl::J_PR_DOMAIN_TYPE_TIME) {
            isTime = true;
        } else if (typeLabel == jlbl::J_PR_DOMAIN_TYPE_FREQ) {
            isFrequency = true;
        } else if (typeLabel == jlbl::J_PR_DOMAIN_TYPE_TIMEFREQ) {
            isTime = true;
            isFrequency = true;
        }

        if (isTime && !isFrequency && !hasTransferFunction) return NFDE::NP_T2_TIME;
        if (!isTime && isFrequency && !hasTransferFunction) return NFDE::NP_T2_FREQ;
        if (!isTime && !isFrequency && hasTransferFunction) return NFDE::NP_T2_TRANSFER;
        if (isTime && isFrequency && !hasTransferFunction) return NFDE::NP_T2_TIMEFREQ;
        if (isTime && !isFrequency && hasTransferFunction) return NFDE::NP_T2_TIMETRANSF;
        if (!isTime && isFrequency && hasTransferFunction) return NFDE::NP_T2_FREQTRANSF;
        if (isTime && isFrequency && hasTransferFunction) return NFDE::NP_T2_TIMEFRECTRANSF;
        Report::WarnErrReport("Invalid domain type for probe.", true);
        return NFDE::NP_T2_TIME;
    }

    // ---- Missing method implementations ----

    std::string parser_t::adaptName(const std::string& str) {
        std::string res = str;
        while (!res.empty() && (res.front() == ' ' || res.front() == '\t')) res.erase(res.begin());
        while (!res.empty() && (res.back() == ' ' || res.back() == '\t')) res.pop_back();
        for (char& c : res) {
            if (c == ' ') c = '_';
        }
        return res;
    }

    void parser_t::checkIsValidName(const std::string& str) {
        if (str.find('@') != std::string::npos) {
            Report::WarnErrReport("ERROR in name: " + str + " contains invalid character @", true);
        }
    }

    std::string parser_t::buildTagName(int matId, int elementId) {
        std::string matName;
        {
            const nlohmann::json* mat = id_map_m::findById(matTable, matId);
            bool found = false;
            matName = getStrAt(mat, jlbl::J_NAME, "", &found);
            if (!found) matName = "material" + std::to_string(matId);
            matName = adaptName(matName);
        }
        std::string layerName;
        {
            const nlohmann::json* elem = id_map_m::findById(elementTable, elementId);
            bool found = false;
            layerName = getStrAt(elem, jlbl::J_NAME, "", &found);
            if (!found) layerName = "layer" + std::to_string(elementId);
            layerName = adaptName(layerName);
        }
        checkIsValidName(matName);
        checkIsValidName(layerName);
        return matName + "@" + layerName;
    }

    std::vector<parser_t::materialAssociation_t> parser_t::getMaterialAssociations(
        const std::vector<std::string>& materialTypes,
        const std::vector<std::string>& elementLabels) {
        std::vector<materialAssociation_t> res;
        const nlohmann::json* allMatAss = nullptr;
        bool found = false;
        jsonGet(&rootJson, jlbl::J_MATERIAL_ASSOCIATIONS, allMatAss, found);
        if (!found) return res;

        int nMatAss = jsonCount(allMatAss);
        for (int i = 0; i < nMatAss; ++i) {
            const nlohmann::json* mAPtr = nullptr;
            jsonGetChild(allMatAss, i + 1, mAPtr);

            materialAssociation_t mA = parseMaterialAssociation(mAPtr);
            const nlohmann::json* mat = id_map_m::findById(matTable, mA.materialId);
            if (!mat) continue;
            std::string matType = getStrAt(mat, jlbl::J_TYPE, "");

            bool typeMatch = false;
            for (auto& t : materialTypes) {
                if (matType == t) { typeMatch = true; break; }
            }
            if (!typeMatch) continue;

            // Element label filter (matches Fortran isAssociatedWithElementLabel)
            bool labelMatch = true;
            if (!elementLabels.empty()) {
                labelMatch = false;
                for (int eid : mA.elementIds) {
                    const nlohmann::json* elm = id_map_m::findById(elementTable, eid);
                    if (!elm) continue;
                    std::string elemType = getStrAt(elm, jlbl::J_TYPE, "");
                    std::string elemSubtype = getStrAt(elm, jlbl::J_SUBTYPE, "");
                    for (const auto& el : elementLabels) {
                        bool negative = (!el.empty() && el[0] == '-');
                        std::string target = negative ? el.substr(1) : el;
                        while (!target.empty() && (target.front() == ' ' || target.front() == '\t'))
                            target.erase(target.begin());
                        while (!target.empty() && (target.back() == ' ' || target.back() == '\t'))
                            target.pop_back();
                        bool matches = (elemType == target || elemSubtype == target);
                        if (negative) {
                            labelMatch = labelMatch && !matches;
                        } else {
                            labelMatch = labelMatch || matches;
                        }
                    }
                }
            }
            if (!labelMatch) continue;

            mA.matAssType = matType;
            res.push_back(mA);
        }
        return res;
    }

    parser_t::materialAssociation_t parser_t::parseMaterialAssociation(const nlohmann::json* matAss) {
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
        const nlohmann::json* totalResistance =
            findPathNode(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE);
        mA.hasTotalResistance = (totalResistance != nullptr);
        if (mA.hasTotalResistance) {
            int dim = jsonCount(totalResistance);
            if (dim == 0) {
                mA.totalResistance = {getRealAt(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE, 0.0)};
            } else {
                mA.totalResistance = getRealsAt(matAss, jlbl::J_MAT_ASS_TOTAL_RESISTANCE);
            }
        }
        return mA;
    }

    std::vector<Cell::cell_interval_t> parser_t::getSingleVolumeInElementsIds(const nlohmann::json* pw) {
        std::vector<Cell::cell_interval_t> res;
        bool found = false;
        std::vector<int> elemIds = getIntsAt(pw, jlbl::J_ELEMENTIDS, &found);
        if (!found) {
            Report::WarnErrReport("Error reading single volume elementIds label not found.", true);
            return res;
        }
        if (elemIds.empty()) {
            Report::WarnErrReport("Entity elementIds must not be empty.", true);
            return res;
        }
        if (elemIds.size() != 1) {
            Report::WarnErrReport("Entity must contain a single elementId.", true);
            return res;
        }
        bool cRf = false;
        Cell::cell_region_t cR = mesh.getCellRegion(int(elemIds[0]), cRf);
        if (!cRf) {
            Report::WarnErrReport("Entity elementId " + std::to_string(elemIds[0]) + " not found.", true);
            return res;
        }
        for (auto& iv : cR.intervals) {
            if (iv.getType() == Cell::CELL_TYPE_VOXEL) res.push_back(iv);
        }
        if (res.size() != 1) {
            Report::WarnErrReport("Entity must contain a single cell region defining a volume.", true);
        }
        return res;
    }

    std::vector<const nlohmann::json*> parser_t::jsonValueFilterByKeyValue(
        const nlohmann::json* place, const std::string& key, const std::string& value) {
        return jsonValueFilterByKeyValues(place, key, {value});
    }

    std::vector<const nlohmann::json*> parser_t::jsonValueFilterByKeyValues(
        const nlohmann::json* place, const std::string& key, const std::vector<std::string>& values) {
        std::vector<const nlohmann::json*> res;
        if (!place) return res;
        int n = jsonCount(place);
        for (int i = 0; i < n; ++i) {
            const nlohmann::json* child = nullptr;
            jsonGetChild(place, i + 1, child);
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
        return NFDE::F_MUR;
    }

    void parser_t::readThinWires(NFDE::ThinWires_t& res, NFDE::MasSondas_t& sonda) {
        auto mwires = getMaterialAssociations({
            std::string(jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE) + "  ",
            jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE
        }, {});
        if (!mwires.empty()) {
            Report::WarnErrReport(
                "ERROR: shieldedMultiwires and unshieldedMultiwires can only be defined if compiled with MTLN",
                true);
        }

        wireMAs_ = getMaterialAssociations({jlbl::J_MAT_TYPE_WIRE}, {});
        wireRes_ = &res;
        wireNGlobal_ = 0;
        wireNNodes_ = 0;
        wireNodeCoordIds_.assign(std::max<size_t>(1, 2 * wireMAs_.size()), 0);
        wireNodeNodeIdx_.assign(std::max<size_t>(1, 2 * wireMAs_.size()), 0);

        int nTw = 0;
        for (const auto& mA : wireMAs_) {
            if (isThinWire(mA)) nTw++;
        }
        res.tw.resize(nTw);
        res.n_tw = nTw;
        res.n_tw_max = nTw;

        int j = 0;
        for (const auto& mA : wireMAs_) {
            if (isThinWire(mA)) {
                res.tw[j] = readThinWire(mA);
                j++;
            }
        }

        const nlohmann::json* allProbes = nullptr;
        bool probesFound = false;
        jsonGet(&rootJson, jlbl::J_PROBES, allProbes, probesFound);
        if (probesFound) {
            auto wireProbePs = jsonValueFilterByKeyValue(allProbes, jlbl::J_TYPE, jlbl::J_PR_TYPE_WIRE);
            int nWireProbes = static_cast<int>(wireProbePs.size());
            if (nWireProbes > 0) {
                int nExisting = sonda.length;
                std::vector<NFDE::MasSonda_t> newCollection(nExisting + nWireProbes);
                for (int k = 0; k < nExisting; ++k) newCollection[k] = sonda.collection[k];
                for (int k = 0; k < nWireProbes; ++k) {
                    newCollection[nExisting + k] = readWireProbe(wireProbePs[k]);
                }
                sonda.collection = std::move(newCollection);
                sonda.length = nExisting + nWireProbes;
                sonda.length_max = nExisting + nWireProbes;
            }
        }

        for (const auto& probe : sonda.collection) {
            if (static_cast<int>(probe.cordinates.size()) > sonda.len_cor_max) {
                sonda.len_cor_max = static_cast<int>(probe.cordinates.size());
            }
        }

        wireRes_ = nullptr;
    }

    parser_t::parser_t(const std::string& filename) {
        this->filename = filename;
        std::ifstream input(filename);
        if (!input.is_open()) {
            Report::WarnErrReport("Failed to load JSON file: " + filename, true);
            return;
        }
        try {
            input >> rootJson;
        } catch (...) {
            Report::WarnErrReport("Failed to parse JSON file: " + filename, true);
            return;
        }
        isInitialized = true;
    }

    NFDE::Parseador_t parser_t::readProblemDescription() {
        NFDE::Parseador_t res;
        
        mesh = readMesh();
        matTable = id_map_m::buildIdMap(rootJson, jlbl::J_MATERIALS);
        elementTable = id_map_m::buildIdMap(
            rootJson, std::string(jlbl::J_MESH) + "." +
                          std::string(jlbl::J_ELEMENTS));
        
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
        const nlohmann::json* jcs = nullptr;
        const nlohmann::json* jc = nullptr;
        bool found = false;
        
        jsonGet(&rootJson, std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_COORDINATES), jcs, found);
        if (found) {
            int numberOfCoordinates = jsonCount(jcs);
            mesh.allocateCoordinates(50 * numberOfCoordinates);
            for (int i = 1; i <= numberOfCoordinates; ++i) {
                jsonGetChild(jcs, i, jc);
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
        const nlohmann::json* jes = nullptr;
        const nlohmann::json* je = nullptr;
        bool found = false;
        
        jsonGet(&rootJson, std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_ELEMENTS), jes, found);
        int numberOfElements = jsonCount(jes);
        mesh.allocateElements(50 * numberOfElements);
            
        if (found) {
            for (int i = 1; i <= numberOfElements; ++i) {
                jsonGetChild(jes, i, je);
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
                    const nlohmann::json* triangles = nullptr;
                    jsonGet(je, jlbl::J_CONF_VOLUME_TRIANGLES, triangles, isConformal);
                    
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

    std::vector<Cell::cell_interval_t> parser_t::readCellIntervals(const nlohmann::json* place, const std::string& path) {
        const nlohmann::json* intervalsPlace = nullptr;
        const nlohmann::json* interval = nullptr;
        bool containsInterval = false;
        
        jsonGet(place, path, intervalsPlace, containsInterval);
        if (!containsInterval) {
            return std::vector<Cell::cell_interval_t>();
        }
        int nIntervals = jsonCount(intervalsPlace);
        std::vector<Cell::cell_interval_t> res(nIntervals);
        for (int i = 1; i <= nIntervals; ++i) {
            jsonGetChild(intervalsPlace, i, interval);
            if (!interval) continue;
            std::vector<double> cellIni = getRealsAt(interval, "(1)");
            std::vector<double> cellEnd = getRealsAt(interval, "(2)");
            if (cellIni.size() < 3 || cellEnd.size() < 3) continue;
            res[i-1].ini.cell[0] = cellIni[0];
            res[i-1].ini.cell[1] = cellIni[1];
            res[i-1].ini.cell[2] = cellIni[2];
            res[i-1].end.cell[0] = cellEnd[0];
            res[i-1].end.cell[1] = cellEnd[1];
            res[i-1].end.cell[2] = cellEnd[2];
        }
        return res;
    }

    std::vector<Conf::triangle_t> parser_t::readTriangles(const nlohmann::json* place, const std::string& path) {
        const nlohmann::json* triangles = nullptr;
        const nlohmann::json* triangle_ptr = nullptr;
        bool containsTriangles = false;
        
        jsonGet(place, path, triangles, containsTriangles);
        if (!containsTriangles) {
            return std::vector<Conf::triangle_t>();
        }
        int nTriangles = jsonCount(triangles);
        std::vector<Conf::triangle_t> res(nTriangles);
        for (int i = 1; i <= nTriangles; ++i) {
            jsonGetChild(triangles, i, triangle_ptr);
            std::vector<double> triangle = getRealsAt(triangle_ptr, ""); // Assuming get returns vector for array
            // Note: json-fortran get to vector might need specific handling, assuming helper exists
            for (int j = 1; j <= 3; ++j) {
                res[i-1].vertices[j-1].id = static_cast<int>(triangle[j-1]);
            }
        }
        return res;
    }

    std::string parser_t::readAdditionalArguments() {
        return getStrAt(&rootJson, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_ADDITIONAL_ARGUMENTS), " ");
    }

    NFDE::NFDEGeneral_t parser_t::readGeneral() {
        NFDE::NFDEGeneral_t res;
        res.dt = getRealAt(&rootJson, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_TIME_STEP), 0.0);
        if (res.dt < 0) Report::WarnErrReport("timStep cannot be negative", true);
        res.nmax = getRealAt(&rootJson, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_NUMBER_OF_STEPS), 0.0);
        if (res.nmax <= 0) Report::WarnErrReport("numberOfSteps has to be positive", true);
        res.mtlnProblem = getLogicalAt(&rootJson, std::string(jlbl::J_GENERAL) + "." + std::string(jlbl::J_GEN_MTLN_PROBLEM), false);
        return res;
    }

    NFDE::MatrizMedios_t parser_t::readMediaMatrix() {
        NFDE::MatrizMedios_t res;
        std::string P = std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_GRID) + "." + jlbl::J_GRID_NUMBER_OF_CELLS;
        res.totalX = getIntAt(&rootJson, P + ".(1)", 0) + 1;
        res.totalY = getIntAt(&rootJson, P + ".(2)", 0) + 1;
        res.totalZ = getIntAt(&rootJson, P + ".(3)", 0) + 1;
        return res;
    }

    NFDE::Desplazamiento_t parser_t::readGrid() {
        NFDE::Desplazamiento_t res;
        std::string P = std::string(jlbl::J_MESH) + "." + std::string(jlbl::J_GRID);
        
        int nX = getIntAt(&rootJson, std::string(P) + "." + std::string(jlbl::J_GRID_NUMBER_OF_CELLS) + ".(1)", 0);
        int nY = getIntAt(&rootJson, std::string(P) + "." + std::string(jlbl::J_GRID_NUMBER_OF_CELLS) + ".(2)", 0);
        int nZ = getIntAt(&rootJson, std::string(P) + "." + std::string(jlbl::J_GRID_NUMBER_OF_CELLS) + ".(3)", 0);

        res.nX = nX;
        res.nY = nY;
        res.nZ = nZ;

        // Helper lambda for assignDes logic
        auto assignDes = [&](const std::string& path, std::vector<double>& dest, int& n) {
            std::vector<double> vec = getRealsAt(&rootJson, path);
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

        res.originx = getRealAt(&rootJson, std::string(P) + "." + std::string(jlbl::J_GRID_ORIGIN) + ".(1)", 0.0);
        res.originy = getRealAt(&rootJson, std::string(P) + "." + std::string(jlbl::J_GRID_ORIGIN) + ".(2)", 0.0);
        res.originz = getRealAt(&rootJson, std::string(P) + "." + std::string(jlbl::J_GRID_ORIGIN) + ".(3)", 0.0);

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
        const nlohmann::json* bdrs = nullptr;
        bool found = false;
        
        jsonGet(&rootJson, jlbl::J_BOUNDARY, bdrs, found);
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
                    for (int _pf = 0; _pf < 6; ++_pf) {
                        res.propiedadesPML[_pf] = readPMLProperties(
                            std::string(jlbl::J_BOUNDARY) + "." + std::string(jlbl::J_BND_ALL));
                    }
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
                    res.propiedadesPML[j] = readPMLProperties(
                        std::string(jlbl::J_BOUNDARY) + "." + placeLabels[i]);
                }
            }
        }

        return res;
    }

    void parser_t::readBackgroundMaterial(NFDE::Materials_t& mats) {
        bool found = false;
        double val = getRealAt(&rootJson, std::string(jlbl::J_BACKGROUND) + "." + std::string(jlbl::J_BKG_ABS_PERMITTIVITY), 0.0, &found);
        if (found) mats.Mats[0].eps = val;

        found = false;
        val = getRealAt(&rootJson, std::string(jlbl::J_BACKGROUND) + "." + std::string(jlbl::J_BKG_ABS_PERMEABILITY), 0.0, &found);
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
            appendRegion(res.Lins, res.nLins, res.nLins_max, emptyCoords);
            appendRegion(res.Surfs, res.nSurfs, res.nSurfs_max, emptyCoords);
            appendRegion(res.Vols, res.nVols, res.nVols_max, emptyCoords);
            return res;
        }

        for (size_t i = 0; i < mAs.size(); ++i) {
            std::vector<NFDE::coords_t> cs;
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_LINEL);
            appendRegion(res.Lins, res.nLins, res.nLins_max, cs);
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_SURFEL);
            appendRegion(res.Surfs, res.nSurfs, res.nSurfs_max, cs);
            matAssToCoords(mAs[i], cs, Cell::CELL_TYPE_VOXEL);
            appendRegion(res.Vols, res.nVols, res.nVols_max, cs);
        }
        return res;
    }

    void parser_t::matAssToCoords(const materialAssociation_t& mA, std::vector<NFDE::coords_t>& res, int cellType) {
        int nCs = 0;
        for (size_t e = 0; e < mA.elementIds.size(); ++e) {
            bool cRf = false;
            Cell::cell_region_t cR = mesh.getCellRegion(int(mA.elementIds[e]), cRf);
            if (cRf) {
                nCs += static_cast<int>(Pt::cellRegionToCoords(cR, cellType).size());
            }
        }

        res.resize(nCs);
        int jIni = 0;
        for (size_t e = 0; e < mA.elementIds.size(); ++e) {
            bool cRf = false;
            Cell::cell_region_t cR = mesh.getCellRegion(int(mA.elementIds[e]), cRf);
            if (!cRf) continue;
            std::string tagName = buildTagName(mA.materialId, mA.elementIds[e]);
            auto newCoords = Pt::cellRegionToCoords(cR, cellType, tagName);
            if (newCoords.empty()) continue;
            for (size_t k = 0; k < newCoords.size(); ++k) {
                res[jIni + static_cast<int>(k)] = newCoords[k];
            }
            jIni += static_cast<int>(newCoords.size());
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
        res.c2P = c2p;
        res.n_C2P = static_cast<int32_t>(c2p.size());
        
        const nlohmann::json* matPtr = id_map_m::findById(matTable, mA.materialId);
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
        res.c2P = c2p;
        res.n_C2P = static_cast<int32_t>(c2p.size());
        
        const nlohmann::json* matPtr = id_map_m::findById(matTable, mA.materialId);
        
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
            bool cRf = false;
            Cell::cell_region_t cR = mesh.getCellRegion(int(mA.elementIds[e]), cRf);
            if (cRf && !Pt::cellRegionToCoords(cR, cellType).empty()) return true;
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

        const nlohmann::json* mat = id_map_m::findById(matTable, mA.materialId);
        res.files = getStrAt(mat, jlbl::J_NAME, " ");
        
        const nlohmann::json* layers = nullptr;
        bool tmpFound; jsonGet(mat, std::string(jlbl::J_MAT_MULTILAYERED_SURF_LAYERS), layers, tmpFound);

        res.numcapas = jsonCount(layers);
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
            const nlohmann::json* layer = nullptr;
            jsonGetChild(layers, i + 1, layer);
            res.sigma[i] = getRealAt(layer, jlbl::J_MAT_ELECTRIC_CONDUCTIVITY, 0.0);
            res.sigmam[i] = getRealAt(layer, jlbl::J_MAT_MAGNETIC_CONDUCTIVITY, 0.0);
            bool hasAbsPermittivity = false;
            res.eps[i] = getRealAt(layer, jlbl::J_MAT_ABS_PERMITTIVITY, 0.0, &hasAbsPermittivity);
            if (!hasAbsPermittivity) {
                res.eps[i] = getRealAt(layer, jlbl::J_MAT_REL_PERMITTIVITY, 1.0) * NFDE::EPSILON_VACUUM;
            }
            bool hasAbsPermeability = false;
            res.mu[i] = getRealAt(layer, jlbl::J_MAT_ABS_PERMEABILITY, 0.0, &hasAbsPermeability);
            if (!hasAbsPermeability) {
                res.mu[i] = getRealAt(layer, jlbl::J_MAT_REL_PERMEABILITY, 1.0) * NFDE::MU_VACUUM;
            }
            bool found = false;
            res.thk[i] = getRealAt(layer, jlbl::J_MAT_MULTILAYERED_SURF_THICKNESS, 0.0, &found);
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
        const nlohmann::json* sources = nullptr;
        bool found = false;

        jsonGet(&rootJson, jlbl::J_SOURCES, sources, found);
        
        if (!found) {
            res.collection.clear();
            res.nc = 0;
            res.nC_max = 0;
            return res;
        }

        std::vector<const nlohmann::json*> pws = jsonValueFilterByKeyValue(sources, jlbl::J_TYPE, jlbl::J_SRC_TYPE_PW);

        res.collection.resize(pws.size());
        for (size_t i = 0; i < pws.size(); ++i) {
            res.collection[i] = readPlanewave(pws[i]);
        }
        res.nc = pws.size();
        res.nC_max = pws.size();

        return res;
    }

    NFDE::PlaneWave_t parser_t::readPlanewave(const nlohmann::json* pw) {
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
        const nlohmann::json* sources = nullptr;
        bool found = false;

        jsonGet(&rootJson, jlbl::J_SOURCES, sources, found);
        if (!found) {
            res.NodalSource.clear();
            return res;
        }

        std::vector<const nlohmann::json*> nodSrcs = jsonValueFilterByKeyValues(sources, jlbl::J_TYPE, {jlbl::J_SRC_TYPE_NS});
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

    NFDE::Curr_Field_Src_t parser_t::readField(const nlohmann::json* jns) {
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
        const nlohmann::json* allProbes = nullptr;
        bool found = false;

        jsonGet(&rootJson, jlbl::J_PROBES, allProbes, found);
        if (!found) {
            res.probes.clear();
            res.n_probes = 0;
            res.n_probes_max = 0;
            return res;
        }

        std::vector<std::string> validTypes = {jlbl::J_PR_TYPE_FARFIELD};
        std::vector<const nlohmann::json*> ps = jsonValueFilterByKeyValues(allProbes, jlbl::J_TYPE, validTypes);

        res.n_probes = ps.size();
        res.n_probes_max = ps.size();
        res.probes.resize(ps.size());
        for (size_t i = 0; i < ps.size(); ++i) {
            res.probes[i] = readFarFieldProbe(ps[i]);
        }
        return res;
    }

    NFDE::abstractSonda_t parser_t::readFarFieldProbe(const nlohmann::json* p) {
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

        if (domain.hasTransferFunction) {
            ff->FileNormalize = domain.filename;
        } else {
            bool transferFunctionFound = false;
            std::string fn = getStrAt(p, std::string(jlbl::J_PR_DOMAIN) + "." + jlbl::J_PR_DOMAIN_MAGNITUDE_FILE,
                                      " ", &transferFunctionFound);
            if (!transferFunctionFound) {
                const nlohmann::json* sources = nullptr;
                bool sourcesFound = false;
                jsonGet(&rootJson, jlbl::J_SOURCES, sources, sourcesFound);
                if (sourcesFound && jsonCount(sources) == 1) {
                    const nlohmann::json* src = nullptr;
                    jsonGetChild(sources, 1, src);
                    fn = getStrAt(src, jlbl::J_SRC_MAGNITUDE_FILE, " ", &transferFunctionFound);
                }
            }
            if (transferFunctionFound) {
                while (!fn.empty() && (fn.front() == ' ' || fn.front() == '\t')) fn.erase(fn.begin());
                while (!fn.empty() && (fn.back() == ' ' || fn.back() == '\t')) fn.pop_back();
                ff->FileNormalize = fn;
            } else {
                ff->FileNormalize = " ";
            }
        }

        if (domain.isLogarithmicFrequencySpacing) {
            appendLogSufix(ff->outputrequest);
        }

        {
            std::vector<NFDE::coords_t> nfdeCoords =
                Pt::cellIntervalsToCoords(getSingleVolumeInElementsIds(p));
            ff->n_cord = 2;
            ff->n_cord_max = 2;
            ff->i = {nfdeCoords[0].Xi, nfdeCoords[0].Xe};
            ff->j = {nfdeCoords[0].Yi, nfdeCoords[0].Ye};
            ff->K = {nfdeCoords[0].Zi, nfdeCoords[0].Ze};
            ff->node.clear();
        }

        readDirection(p, jlbl::J_PR_FAR_FIELD_PHI, ff->phistart, ff->phistop, ff->phistep);
        readDirection(p, jlbl::J_PR_FAR_FIELD_THETA, ff->thetastart, ff->thetastop, ff->thetastep);
        return res;
    }

    void parser_t::readDirection(const nlohmann::json* p, const std::string& label, double& initial, double& final, double& step) {
        const nlohmann::json* dir = nullptr;
        bool found = false;
        jsonGet(p, label, dir, found);
        if (!found) {
            Report::WarnErrReport("Error reading far field probe. Direction label not found.", true);
        }
        initial = getRealAt(dir, jlbl::J_PR_FAR_FIELD_DIR_INITIAL, 0.0);
        final = getRealAt(dir, jlbl::J_PR_FAR_FIELD_DIR_FINAL, 0.0);
        step = getRealAt(dir, jlbl::J_PR_FAR_FIELD_DIR_STEP, 0.0);
    }

    NFDE::MasSondas_t parser_t::readMoreProbes() {
        NFDE::MasSondas_t res;
        const nlohmann::json* allProbes = nullptr;
        bool found = false;

        jsonGet(&rootJson, jlbl::J_PROBES, allProbes, found);
        if (!found) {
            res.collection.clear();
            res.length = 0;
            res.length_max = 0;
            res.len_cor_max = 0;
            return res;
        }

        std::vector<std::string> validTypes = {jlbl::J_PR_TYPE_POINT, jlbl::J_PR_TYPE_LINE};
        std::vector<const nlohmann::json*> ps = jsonValueFilterByKeyValues(allProbes, jlbl::J_TYPE, validTypes);
        
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

    bool parser_t::isMoreProbe(const nlohmann::json* p) {
        return isPointProbe(p) || isLineProbe(p);
    }

    bool parser_t::isLineProbe(const nlohmann::json* p) {
        return getStrAt(p, jlbl::J_TYPE) == jlbl::J_PR_TYPE_LINE;
    }

    bool parser_t::isPointProbe(const nlohmann::json* p) {
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

    NFDE::MasSonda_t parser_t::readLineProbe(const nlohmann::json* p) {
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

    NFDE::MasSonda_t parser_t::readPointProbe(const nlohmann::json* p) {
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
            const nlohmann::json* dirLabelPtr = nullptr;
            bool dirLabelsFound = false;
            jsonGet(p, jlbl::J_PR_POINT_DIRECTIONS, dirLabelPtr, dirLabelsFound);
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

    std::vector<char> parser_t::buildDirLabels(const nlohmann::json* dirLabelsPtr) {
        std::vector<char> res(jsonCount(dirLabelsPtr));
        for (int i = 0; i < jsonCount(dirLabelsPtr); ++i) {
            const nlohmann::json* child = nullptr;
            jsonGetChild(dirLabelsPtr, i + 1, child);
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
        if (domain.hasTransferFunction) {
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
        std::vector<const nlohmann::json*> bps;
        const nlohmann::json* probes = nullptr;
        bool found = false;

        jsonGet(&rootJson, jlbl::J_PROBES, probes, found);
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

    NFDE::BloqueProbe_t parser_t::readBlockProbe(const nlohmann::json* bp) {
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
        if (domain.hasTransferFunction) {
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
        std::vector<const nlohmann::json*> ps;
        const nlohmann::json* probes = nullptr;
        bool found = false;

        jsonGet(&rootJson, jlbl::J_PROBES, probes, found);
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

    NFDE::VolProbe_t parser_t::readVolProbe(const nlohmann::json* p) {
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
        const nlohmann::json* compsPtr = nullptr;
        bool componentsFound = false;
        jsonGet(p, jlbl::J_PR_MOVIE_COMPONENT, compsPtr, componentsFound);
        
        res.cordinates.resize(1);
        res.cordinates[0] = cs[0];
        std::string component;
        if (componentsFound) {
            component = getStrAt(compsPtr, "");
        } else {
            component = jlbl::J_DIR_M;
        }
        res.cordinates[0].Or = buildVolProbeType(fieldType, component);
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
        if (domain.hasTransferFunction) {
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
        const nlohmann::json* mat = id_map_m::findById(matTable, mA.materialId);
        bool found = false;
        res.width = getRealAt(mat, jlbl::J_MAT_THINSLOT_WIDTH, 0.0, &found);
        if (!found) {
            Report::WarnErrReport("Missing thin slot width for material " + std::to_string(mA.materialId), true);
        }
        std::vector<NFDE::coords_t> coords;
        matAssToCoords(mA, coords, Cell::CELL_TYPE_LINEL);
        coordsToThinSlotComp(coords, res.tgc);
        res.n_tgc = static_cast<int32_t>(res.tgc.size());
        return res;
    }

    void parser_t::coordsToThinSlotComp(const std::vector<NFDE::coords_t>& cs,
                                         std::vector<NFDE::ThinSlotComp_t>& tc) {
        int nTgc = 0;
        for (const auto& c : cs) {
            nTgc += (c.Xe - c.Xi + 1) * (c.Ye - c.Yi + 1) * (c.Ze - c.Zi + 1);
        }
        tc.resize(nTgc);
        int j = 0;
        for (const auto& c : cs) {
            switch (std::abs(c.Or)) {
            case NFDE::iEx:
                for (int k = 1; k <= c.Xe - c.Xi + 1; ++k) {
                    tc[j] = buildBaseThinSlotComponent(c);
                    tc[j].i = c.Xi + k - 1;
                    ++j;
                }
                break;
            case NFDE::iEy:
                for (int k = 1; k <= c.Ye - c.Yi + 1; ++k) {
                    tc[j] = buildBaseThinSlotComponent(c);
                    tc[j].j = c.Yi + k - 1;
                    ++j;
                }
                break;
            case NFDE::iEz:
                for (int k = 1; k <= c.Ze - c.Zi + 1; ++k) {
                    tc[j] = buildBaseThinSlotComponent(c);
                    tc[j].K = c.Zi + k - 1;
                    ++j;
                }
                break;
            default:
                break;
            }
        }
    }

    NFDE::ThinSlotComp_t parser_t::buildBaseThinSlotComponent(const NFDE::coords_t& cs) {
        NFDE::ThinSlotComp_t res;
        res.i = cs.Xi;
        res.j = cs.Yi;
        res.K = cs.Zi;
        res.dir = std::abs(cs.Or);
        res.tag = cs.tag;
        return res;
    }

    NFDE::FronteraPML_t parser_t::readPMLProperties(const std::string& path) {
        NFDE::FronteraPML_t res;
        res.numCapas = getIntAt(&rootJson, path + "." + std::string(jlbl::J_BND_PML_LAYERS), 8);
        res.orden = getRealAt(&rootJson, path + "." + std::string(jlbl::J_BND_PML_ORDER), 2.0);
        res.refl = getRealAt(&rootJson, path + "." + std::string(jlbl::J_BND_PML_REFLECTION), 0.001);
        return res;
    }

    int parser_t::strToTerminationType(const std::string& label) {
        if (label == jlbl::J_MAT_TERM_TYPE_OPEN) return NFDE::MATERIAL_CONS;
        if (label == jlbl::J_MAT_TERM_TYPE_SERIES) return NFDE::SERIES_CONS;
        if (label == jlbl::J_MAT_TERM_TYPE_SHORT) return NFDE::MATERIAL_CONS;
        return NFDE::MATERIAL_CONS;
    }

    parser_t::thinwiretermination_t parser_t::readThinWireTermination(const nlohmann::json* terminal) {
        thinwiretermination_t res;
        const nlohmann::json* tms = nullptr;
        bool found = false;
        jsonGet(terminal, jlbl::J_MAT_TERM_TERMINATIONS, tms, found);
        if (!found) {
            Report::WarnErrReport("Error reading wire terminal. terminations not found.", true);
        }
        if (jsonCount(tms) != 1) {
            Report::WarnErrReport("Only terminals with a single termination are allowed for a wire.", true);
        }
        const nlohmann::json* tm = nullptr;
        jsonGetChild(tms, 1, tm);
        bool labelFound = false;
        std::string label = getStrAt(tm, jlbl::J_TYPE, "", &labelFound);
        res.terminationType = strToTerminationType(label);
        if (!labelFound) {
            Report::WarnErrReport("Error reading wire terminal. termination must specify a type.", true);
        }
        if (label == jlbl::J_MAT_TERM_TYPE_OPEN || label == jlbl::J_MAT_TERM_TYPE_SHORT) {
            res.r = 0.0; res.l = 0.0; res.c = 0.0;
        } else if (label == jlbl::J_MAT_TERM_TYPE_SERIES) {
            res.r = getRealAt(tm, jlbl::J_MAT_TERM_RESISTANCE, 0.0);
            res.l = getRealAt(tm, jlbl::J_MAT_TERM_INDUCTANCE, 0.0);
            res.c = getRealAt(tm, jlbl::J_MAT_TERM_CAPACITANCE, 1e22);
        } else {
            Report::WarnErrReport(
                "Error reading wire terminal. Holland wires only support open, short, and series terminations",
                true);
        }
        return res;
    }

    bool parser_t::isThinWire(const materialAssociation_t& mA) {
        if (mA.elementIds.size() != 1) {
            Report::WarnErrReport("Thin wires must be defined by a single element id.", true);
        }
        bool found = false;
        auto pl = mesh.getPolyline(mA.elementIds[0], found);
        if (!found || !mesh.arePolylineSegmentsStructured(pl)) {
            Report::WarnErrReport("Thin wires must be defined by a structured polyline.", true);
        }
        return true;
    }

    int parser_t::getOrAssignNodeIndex(int coordId) {
        for (int k = 0; k < wireNNodes_; ++k) {
            if (wireNodeCoordIds_[k] == coordId) {
                return wireNodeNodeIdx_[k];
            }
        }
        wireNGlobal_++;
        wireNNodes_++;
        if (wireNNodes_ > static_cast<int>(wireNodeCoordIds_.size())) {
            wireNodeCoordIds_.resize(wireNNodes_);
            wireNodeNodeIdx_.resize(wireNNodes_);
        }
        wireNodeCoordIds_[wireNNodes_ - 1] = coordId;
        wireNodeNodeIdx_[wireNNodes_ - 1] = wireNGlobal_;
        return wireNGlobal_;
    }

    int parser_t::orientFieldFromGenerator(const std::vector<Mesh::linel_t>& linels, int position) {
        int orient = linels[position - 1].orientation;
        if (position == 1) {
            return (orient >= 0) ? 1 : -1;
        }
        if (position == static_cast<int>(linels.size())) {
            return (orient >= 0) ? -1 : 1;
        }
        return (orient >= 0) ? 1 : -1;
    }

    int parser_t::findSourcePositionInLinels(const std::vector<int>& srcElemIds,
                                              const std::vector<Mesh::linel_t>& linels) {
        bool found = false;
        auto node = mesh.getNode(srcElemIds[0], found);
        if (!found) {
            Report::WarnErrReport("Source could not be found in linels.", true);
        }
        auto pixel = mesh.nodeToPixel(node);
        for (size_t i = 0; i < linels.size(); ++i) {
            if (linels[i].tag == pixel.tag) {
                return static_cast<int>(i + 1);
            }
        }
        Report::WarnErrReport("Source could not be found in linels.", true);
        return 0;
    }

    std::vector<parser_t::generator_description_t> parser_t::readGeneratorOnThinWire(
        const std::vector<Mesh::linel_t>& linels, const std::vector<int>& plineElemIds) {
        std::vector<generator_description_t> res(linels.size());
        for (auto& g : res) {
            g.srcfile = "None";
            g.srctype = "None";
            g.multiplier = 0.0;
        }

        const nlohmann::json* sources = nullptr;
        bool found = false;
        jsonGet(&rootJson, jlbl::J_SOURCES, sources, found);
        if (!found) return res;

        auto genSrcs = jsonValueFilterByKeyValues(sources, jlbl::J_TYPE, {jlbl::J_SRC_TYPE_GEN});
        if (genSrcs.empty()) return res;

        for (auto* genPtr : genSrcs) {
            std::vector<int> sourceElemIds = getIntsAt(genPtr, jlbl::J_ELEMENTIDS);
            bool nodeFound = false;
            auto srcNode = mesh.getNode(sourceElemIds[0], nodeFound);
            bool plFound = false;
            auto polylineCoords = mesh.getPolyline(plineElemIds[0], plFound);
            if (!nodeFound || !plFound) continue;
            bool inPolyline = false;
            for (int cid : polylineCoords.coordIds) {
                if (!srcNode.coordIds.empty() && cid == srcNode.coordIds[0]) {
                    inPolyline = true;
                    break;
                }
            }
            if (!inPolyline) continue;

            int position = findSourcePositionInLinels(sourceElemIds, linels);
            if (findPathNode(genPtr, jlbl::J_SRC_MAGNITUDE_FILE) == nullptr) {
                Report::WarnErrReport("magnitudeFile of source missing", true);
            }
            std::string field = getStrAt(genPtr, jlbl::J_FIELD);
            if (field == jlbl::J_FIELD_VOLTAGE) {
                res[position - 1].srctype = "VOLT";
            } else if (field == jlbl::J_FIELD_CURRENT) {
                res[position - 1].srctype = "CURR";
            } else {
                Report::WarnErrReport(
                    "Field block of source of type generator must be current or voltage", true);
            }
            res[position - 1].srcfile = getStrAt(genPtr, jlbl::J_SRC_MAGNITUDE_FILE);
            res[position - 1].multiplier = static_cast<double>(orientFieldFromGenerator(linels, position));
        }
        return res;
    }

    NFDE::ThinWire_t parser_t::readThinWire(const materialAssociation_t& cable) {
        NFDE::ThinWire_t res;
        const nlohmann::json* mat = id_map_m::findById(matTable, cable.materialId);
        res.rad = getRealAt(mat, jlbl::J_MAT_WIRE_RADIUS, 0.0);
        res.res = getRealAt(mat, jlbl::J_MAT_WIRE_RESISTANCE, 0.0);
        res.ind = getRealAt(mat, jlbl::J_MAT_WIRE_INDUCTANCE, 0.0);
        res.dispfile = "";
        res.dispfile_LeftEnd = "";
        res.dispfile_RightEnd = "";

        {
            const nlohmann::json* terminal =
                id_map_m::findById(matTable, cable.initialTerminalId);
            auto term = readThinWireTermination(terminal);
            res.tl = term.terminationType;
            res.R_LeftEnd = term.r;
            res.L_LeftEnd = term.l;
            res.C_LeftEnd = term.c;
            res.dispfile_LeftEnd = "";
        }
        {
            const nlohmann::json* terminal =
                id_map_m::findById(matTable, cable.endTerminalId);
            auto term = readThinWireTermination(terminal);
            res.tr = term.terminationType;
            res.R_RightEnd = term.r;
            res.L_RightEnd = term.l;
            res.C_RightEnd = term.c;
            res.dispfile_RightEnd = "";
        }

        bool plFound = false;
        auto polyline = mesh.getPolyline(cable.elementIds[0], plFound);
        auto linels = mesh.polylineToLinels(polyline);

        if (cable.hasTotalResistance) {
            auto despl = readGrid();
            double totalLength = 0.0;
            for (const auto& linel : linels) {
                double stepSize = 0.0;
                int dir = std::abs(linel.orientation);
                if (dir == Mesh::DIR_X) {
                    stepSize = (despl.desX.size() == 1) ? despl.desX[0] : despl.desX[linel.cell[0]];
                } else if (dir == Mesh::DIR_Y) {
                    stepSize = (despl.desY.size() == 1) ? despl.desY[0] : despl.desY[linel.cell[1]];
                } else if (dir == Mesh::DIR_Z) {
                    stepSize = (despl.desZ.size() == 1) ? despl.desZ[0] : despl.desZ[linel.cell[2]];
                }
                totalLength += stepSize;
            }
            res.res = cable.totalResistance[0] / totalLength;
        }

        std::string tagLabel = buildTagName(cable.materialId, cable.elementIds[0]);
        auto genDesc = readGeneratorOnThinWire(linels, cable.elementIds);

        res.n_twc = static_cast<int32_t>(linels.size());
        res.n_twc_max = res.n_twc;
        res.twc.resize(linels.size());
        res.LeftEnd = getOrAssignNodeIndex(polyline.coordIds.front());
        res.RightEnd = getOrAssignNodeIndex(polyline.coordIds.back());

        int baseGlobal = wireNGlobal_;
        for (size_t i = 0; i < linels.size(); ++i) {
            res.twc[i].srcfile = genDesc[i].srcfile;
            res.twc[i].srctype = genDesc[i].srctype;
            res.twc[i].m = genDesc[i].multiplier;
            res.twc[i].i = linels[i].cell[0];
            res.twc[i].j = linels[i].cell[1];
            res.twc[i].K = linels[i].cell[2];
            res.twc[i].d = std::abs(linels[i].orientation);
            res.twc[i].nd = baseGlobal + static_cast<int>(i) + 1;
            res.twc[i].tag = tagLabel;
        }
        wireNGlobal_ = baseGlobal + static_cast<int>(linels.size());
        return res;
    }

    void parser_t::setDomainOfWireProbe(NFDE::MasSonda_t& res, const domain_t& domain) {
        res.tstart = domain.tstart;
        res.tstep = domain.tstep;
        res.tstop = domain.tstop;
        res.fstart = domain.fstart;
        res.fstep = domain.fstep;
        res.fstop = domain.fstop;
        res.filename = domain.filename.empty() ? " " : domain.filename;
        res.type1 = domain.type1;
        res.type2 = domain.type2;
        if (domain.isLogarithmicFrequencySpacing) {
            appendLogSufix(res.outputrequest);
        }
    }

    int parser_t::getSegmentNdWhichMatchesCoord(int coordId, const Mesh::coordinate_t& probe_coord) {
        for (size_t k = 0; k < wireMAs_.size(); ++k) {
            bool plFound = false;
            auto polyline = mesh.getPolyline(wireMAs_[k].elementIds[0], plFound);
            if (!plFound) continue;
            for (int cid : polyline.coordIds) {
                if (cid != coordId) continue;
                auto linels = mesh.polylineToLinels(polyline);
                int n_linels = static_cast<int>(linels.size());
                std::vector<Mesh::coordinate_t> linelCoords(n_linels + 1);
                for (int i = 0; i < n_linels; ++i) {
                    linelCoords[i].position[0] = linels[i].cell[0];
                    linelCoords[i].position[1] = linels[i].cell[1];
                    linelCoords[i].position[2] = linels[i].cell[2];
                    if (linels[i].orientation < 0) {
                        int orDir = std::abs(linels[i].orientation);
                        linelCoords[i].position[orDir - 1] += 1.0;
                    }
                }
                int orDir = linels[n_linels - 1].orientation;
                linelCoords[n_linels] = linelCoords[n_linels - 1];
                linelCoords[n_linels].position[std::abs(orDir) - 1] +=
                    (orDir > 0) ? 1.0 : -1.0;

                int bestIdx = 0;
                double bestDist = 1e30;
                for (int i = 0; i <= n_linels; ++i) {
                    double dx = linelCoords[i].position[0] - probe_coord.position[0];
                    double dy = linelCoords[i].position[1] - probe_coord.position[1];
                    double dz = linelCoords[i].position[2] - probe_coord.position[2];
                    double dist = std::sqrt(dx * dx + dy * dy + dz * dz);
                    if (dist < bestDist) {
                        bestDist = dist;
                        bestIdx = i;
                    }
                }
                int local_idx = std::min(bestIdx + 1, n_linels);
                if (wireRes_ && k < wireRes_->tw.size()) {
                    return wireRes_->tw[k].twc[local_idx - 1].nd;
                }
                return 0;
            }
        }
        return 0;
    }

    NFDE::MasSonda_t parser_t::readWireProbe(const nlohmann::json* p) {
        NFDE::MasSonda_t res;
        bool nameFound = false;
        std::string outputName = getStrAt(p, jlbl::J_NAME, "", &nameFound);
        if (!nameFound) {
            Report::WarnErrReport("Wire probes must define a name.", true);
        }
        res.outputrequest = outputName;
        setDomainOfWireProbe(res, getDomain(p, jlbl::J_PR_DOMAIN));

        bool elementIdsFound = false;
        std::vector<int> elemIds = getIntsAt(p, jlbl::J_ELEMENTIDS, &elementIdsFound);
        if (!elementIdsFound) {
            Report::WarnErrReport("Element ids entry not found for wire probe.", true);
        }
        if (elemIds.size() != 1) {
            Report::WarnErrReport("Wire probe must contain a single element id.", true);
        }

        bool nodeFound = false;
        auto node = mesh.getNode(elemIds[0], nodeFound);
        bool coordFound = false;
        auto probe_coord = mesh.getCoordinate(node.coordIds[0], coordFound);
        std::string fieldLabel = getStrAt(p, jlbl::J_FIELD, jlbl::J_FIELD_VOLTAGE);

        NFDE::coords_t cord;
        cord.tag = outputName;
        cord.Xi = getSegmentNdWhichMatchesCoord(node.coordIds[0], probe_coord);
        cord.Yi = 0;
        cord.Zi = 0;
        if (fieldLabel == jlbl::J_FIELD_CURRENT) {
            cord.Or = NFDE::NP_COR_WIRECURRENT;
        } else if (fieldLabel == jlbl::J_FIELD_VOLTAGE) {
            cord.Or = NFDE::NP_COR_DDP;
        } else {
            Report::WarnErrReport("Invalid field label for wire probe.", true);
        }
        res.cordinates = {cord};
        res.len_cor = 1;
        return res;
    }

} // namespace smbjson
