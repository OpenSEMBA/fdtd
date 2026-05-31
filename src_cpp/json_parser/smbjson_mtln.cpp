#ifdef CompileWithMTLN

#include "smbjson_m.h"
#include "parser_tools_m.h"
#include "cells_m.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <sstream>
#include <unordered_map>

namespace smbjson {

namespace Mtln = mtln_types_m;
namespace Cell = cells_m;
namespace Pt = parser_tools_m;

struct aux_node_t {
    Mtln::terminal_node_t node;
    int cId = 0;
    Mesh::coordinate_t relPos{};
};

class parser_t::MtlnReader {
public:
    MtlnReader(parser_t& parser, Mtln::mtln_t& mtlnRes);

    void read() {
        auto cables = p.getMaterialAssociations(
            {jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_WIRE},
            {});

        mtln_res.connectors = readConnectors();
        addConnIdToConnectorMap(connIdToConnector, mtln_res.connectors);

        if (cables.empty()) {
            mtln_res.time_step = 0.0;
            mtln_res.number_of_steps = 0;
            return;
        }

        {
            auto unshielded = p.getMaterialAssociations(
                {jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE, jlbl::J_MAT_TYPE_WIRE}, {});
            mtln_res.n_unsh = static_cast<int>(unshielded.size());
            auto shielded = p.getMaterialAssociations({jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE}, {});
            mtln_res.n_sh = static_cast<int>(shielded.size());
        }

        if (p.rootJson.contains(jlbl::J_GENERAL) && p.rootJson[jlbl::J_GENERAL].is_object()) {
            const auto& general = p.rootJson[jlbl::J_GENERAL];
            mtln_res.time_step = jsonReal(general, jlbl::J_GEN_TIME_STEP, 0.0);
            mtln_res.number_of_steps = jsonInt(general, jlbl::J_GEN_NUMBER_OF_STEPS, 0);
        } else {
            mtln_res.time_step = 0.0;
            mtln_res.number_of_steps = 0;
        }

        mtln_res.cables.resize(cables.size());
        for (size_t i = 0; i < cables.size(); ++i) {
            auto cable = readMTLNCable(cables[i]);
            stopOnRepeatedName(*cable, mtln_res.cables, static_cast<int>(i));
            mtln_res.cables[i].ptr = std::move(cable);
            addElemIdToCableMap(elemIdToCable, cables[i].elementIds, static_cast<int>(i) + 1);
            addElemIdToPositionMap(elemIdToPosition, cables[i].elementIds);
        }

        for (size_t i = 0; i < cables.size(); ++i) {
            auto* sh = dynamic_cast<Mtln::shielded_multiwire_t*>(mtln_res.cables[i].ptr.get());
            if (sh) {
                sh->parent_cable = assignParentCable(cables[i]);
                sh->conductor_in_parent = assignConductorInParent(cables[i]);
            }
        }

        mtln_res.wireGenerators = readWireGenerators();
        mtln_res.probes = readMultiwireProbes();
        mtln_res.networks = buildNetworks();
    }

private:
    smbjson::parser_t& p;
    Mtln::mtln_t& mtln_res;
    std::unordered_map<int, int> elemIdToCable;
    std::unordered_map<int, int> elemIdToPosition;
    std::unordered_map<int, int> connIdToConnector;

    static std::string jsonString(const nlohmann::json& obj, const char* key, const std::string& fallback = "") {
        if (!obj.is_object() || !obj.contains(key)) return fallback;
        const auto& v = obj[key];
        return v.is_string() ? v.get<std::string>() : fallback;
    }

    static int jsonInt(const nlohmann::json& obj, const char* key, int fallback = 0) {
        if (!obj.is_object() || !obj.contains(key)) return fallback;
        const auto& v = obj[key];
        return v.is_number_integer() ? v.get<int>() : fallback;
    }

    static double jsonReal(const nlohmann::json& obj, const char* key, double fallback = 0.0) {
        if (!obj.is_object() || !obj.contains(key)) return fallback;
        const auto& v = obj[key];
        return v.is_number() ? v.get<double>() : fallback;
    }

    static std::vector<double> jsonRealVector(const nlohmann::json& obj, const char* key) {
        std::vector<double> out;
        if (!obj.is_object() || !obj.contains(key)) return out;
        const auto& arr = obj[key];
        if (!arr.is_array()) return out;
        out.reserve(arr.size());
        for (const auto& v : arr) {
            if (v.is_number()) out.push_back(v.get<double>());
        }
        return out;
    }

    static std::vector<int> jsonIntVector(const nlohmann::json& obj, const char* key) {
        std::vector<int> out;
        if (!obj.is_object() || !obj.contains(key)) return out;
        const auto& arr = obj[key];
        if (!arr.is_array()) return out;
        out.reserve(arr.size());
        for (const auto& v : arr) {
            if (v.is_number_integer()) out.push_back(v.get<int>());
        }
        return out;
    }

    static std::vector<std::vector<double>> jsonMatrix(const nlohmann::json& obj, const char* key) {
        std::vector<std::vector<double>> out;
        if (!obj.is_object() || !obj.contains(key)) return out;
        const auto& arr = obj[key];
        if (!arr.is_array()) return out;
        out.reserve(arr.size());
        for (const auto& row : arr) {
            std::vector<double> rowVals;
            if (row.is_array()) {
                rowVals.reserve(row.size());
                for (const auto& v : row) {
                    if (v.is_number()) rowVals.push_back(v.get<double>());
                }
            }
            out.push_back(std::move(rowVals));
        }
        return out;
    }

    static int jsonDimension(const nlohmann::json& obj, const char* key) {
        if (!obj.is_object() || !obj.contains(key)) return -1;
        const auto& v = obj[key];
        if (!v.is_array()) return 0;
        if (v.empty()) return 1;
        return v.front().is_array() ? 2 : 1;
    }

    const nlohmann::json* materialsArray() const {
        if (!p.rootJson.is_object() || !p.rootJson.contains(jlbl::J_MATERIALS)) return nullptr;
        const auto& mats = p.rootJson[jlbl::J_MATERIALS];
        return mats.is_array() ? &mats : nullptr;
    }

    const nlohmann::json* sourcesArray() const {
        if (!p.rootJson.is_object() || !p.rootJson.contains(jlbl::J_SOURCES)) return nullptr;
        const auto& src = p.rootJson[jlbl::J_SOURCES];
        return src.is_array() ? &src : nullptr;
    }

    const nlohmann::json* probesArray() const {
        if (!p.rootJson.is_object() || !p.rootJson.contains(jlbl::J_PROBES)) return nullptr;
        const auto& pr = p.rootJson[jlbl::J_PROBES];
        return pr.is_array() ? &pr : nullptr;
    }

    const nlohmann::json* materialAssociationsArray() const {
        if (!p.rootJson.is_object() || !p.rootJson.contains(jlbl::J_MATERIAL_ASSOCIATIONS)) return nullptr;
        const auto& ass = p.rootJson[jlbl::J_MATERIAL_ASSOCIATIONS];
        return ass.is_array() ? &ass : nullptr;
    }

    const nlohmann::json* materialById(int materialId) const {
        const auto* mats = materialsArray();
        if (!mats) return nullptr;
        for (const auto& mat : *mats) {
            if (jsonInt(mat, jlbl::J_ID, 0) == materialId) {
                return &mat;
            }
        }
        return nullptr;
    }

    std::string materialTypeById(int materialId) const {
        const auto* mat = materialById(materialId);
        return mat ? jsonString(*mat, jlbl::J_TYPE) : std::string();
    }

    static void addConnIdToConnectorMap(
        std::unordered_map<int, int>& map, const std::vector<Mtln::connector_t>& conn) {
        for (size_t i = 0; i < conn.size(); ++i) {
            map[conn[i].id] = static_cast<int>(i) + 1;
        }
    }

    static void addElemIdToCableMap(
        std::unordered_map<int, int>& map, const std::vector<int>& elemIds, int index) {
        for (int id : elemIds) {
            map[id] = index;
        }
    }

    static void addElemIdToPositionMap(
        std::unordered_map<int, int>& map, const std::vector<int>& elemIds) {
        for (size_t i = 0; i < elemIds.size(); ++i) {
            map[elemIds[i]] = static_cast<int>(i) + 1;
        }
    }

    static void stopOnRepeatedName(
        const Mtln::cable_t& cable,
        const std::vector<Mtln::cable_abstract_t>& cables,
        int n) {
        for (int i = 0; i < n; ++i) {
            if (cables[static_cast<size_t>(i)].ptr && cables[static_cast<size_t>(i)].ptr->name == cable.name) {
                Report::WarnErrReport("Cable name " + cable.name + " has already been used", true);
            }
        }
    }

    std::vector<Mtln::connector_t> readConnectors() {
        const auto* mats = materialsArray();
        if (!mats) {
            return {};
        }

        std::vector<Mtln::connector_t> res;
        for (const auto& mat : *mats) {
            if (jsonString(mat, jlbl::J_TYPE) != jlbl::J_MAT_TYPE_CONNECTOR) {
                continue;
            }
            Mtln::connector_t connector;
            connector.id = jsonInt(mat, jlbl::J_ID, 0);
            connector.resistances = jsonRealVector(mat, jlbl::J_MAT_CONN_RESISTANCES);
            if (mat.contains(jlbl::J_MAT_CONN_TRANSFER_IMPEDANCES) &&
                mat[jlbl::J_MAT_CONN_TRANSFER_IMPEDANCES].is_array()) {
                for (const auto& z : mat[jlbl::J_MAT_CONN_TRANSFER_IMPEDANCES]) {
                    connector.transfer_impedances_per_meter.push_back(readTransferImpedance(z));
                }
            }
            res.push_back(std::move(connector));
        }
        return res;
    }

    std::vector<Mtln::terminal_network_t> buildNetworks() {
        std::vector<aux_node_t> aux_nodes;
        std::vector<Mesh::coordinate_t> networks_coordinates;
        auto cables = p.getMaterialAssociations({jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE}, {});
        auto cables2 = p.getMaterialAssociations({jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE}, {});
        auto cables3 = p.getMaterialAssociations({jlbl::J_MAT_TYPE_WIRE}, {});
        std::vector<smbjson::parser_t::materialAssociation_t> allCables;
        allCables.insert(allCables.end(), cables.begin(), cables.end());
        allCables.insert(allCables.end(), cables2.begin(), cables2.end());
        allCables.insert(allCables.end(), cables3.begin(), cables3.end());

        for (const auto& cable : allCables) {
            const std::string cableType = materialTypeById(cable.materialId);
            bool isShieldedCable = (cableType == jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE);
            const nlohmann::json* terminations_ini = getTerminationsOnSide(cable.initialTerminalId);
            const nlohmann::json* terminations_end = getTerminationsOnSide(cable.endTerminalId);
            for (size_t j = 0; j < cable.elementIds.size(); ++j) {
                int elemId = cable.elementIds[j];
                int conductorIndex = static_cast<int>(j) + 1;
                auto nodeIni = buildNode(terminations_ini, Mtln::TERMINAL_NODE_SIDE_INI, conductorIndex, elemId, isShieldedCable);
                aux_nodes.push_back(nodeIni);
                auto nodeEnd = buildNode(terminations_end, Mtln::TERMINAL_NODE_SIDE_END, conductorIndex, elemId, isShieldedCable);
                aux_nodes.push_back(nodeEnd);
                updateListOfNetworksCoordinates(networks_coordinates, elemId);
            }
        }

        std::vector<Mtln::terminal_network_t> res(networks_coordinates.size());
        for (size_t i = 0; i < networks_coordinates.size(); ++i) {
            res[i] = buildNetwork(networks_coordinates[i], aux_nodes, static_cast<int>(i) + 1);
        }
        return res;
    }

    Mtln::terminal_network_t buildNetwork(
        const Mesh::coordinate_t& network_coordinate,
        const std::vector<aux_node_t>& aux_nodes,
        int network_index) {
        Mtln::terminal_network_t res;
        auto network_nodes = filterNetworkNodesByCoordinate(aux_nodes, network_coordinate);
        auto node_ids = buildListOfNodeIds(network_nodes);
        auto network_circuits = buildNetworkCircuits(network_nodes, node_ids, network_index);
        for (int node_id : node_ids) {
            res.add_connection(buildConnection(node_id, network_nodes, network_circuits));
        }
        return res;
    }

    std::vector<Mtln::network_circuit_t> buildNetworkCircuits(
        const std::vector<aux_node_t>& nodes,
        const std::vector<int>& node_ids,
        int network_index) {
        auto subckt_filtered_nodes = filterNetworkNodesByNetworkCircuit(nodes);
        std::vector<Mtln::network_circuit_t> res;
        for (int node_id : node_ids) {
            auto id_filtered_nodes = filterNetworkNodesById(subckt_filtered_nodes, node_id);
            if (!id_filtered_nodes.empty()) {
                Mtln::network_circuit_t nc;
                nc.nodeId = id_filtered_nodes.front().cId;
                nc.model_name = id_filtered_nodes.front().node.termination.model.name;
                nc.model_file = id_filtered_nodes.front().node.termination.model.file;
                nc.circuit_name = "subckt_" + nc.model_file + "_" + std::to_string(network_index);
                nc.number_of_nodes = readNumberOfNodes(nc.model_file, nc.model_name);
                if (nc.number_of_nodes == 0) {
                    Report::WarnErrReport("Problem in network model. No ports detected", true);
                }
                res.push_back(nc);
            }
        }
        return res;
    }

    static int readNumberOfNodes(const std::string& model_file, const std::string& model_name) {
        std::ifstream in(model_file);
        if (!in) {
            return 0;
        }
        std::string line;
        while (std::getline(in, line)) {
            size_t start = line.find_first_not_of(" \t");
            if (start == std::string::npos) {
                continue;
            }
            line = line.substr(start);
            if (line.empty() || line[0] == '*') {
                continue;
            }
            std::vector<std::string> words;
            Pt::splitLineIntoWords(line, words);
            if (words.size() >= 2 && Pt::to_upper(words[0]) == ".SUBCKT" && words[1] == model_name) {
                return static_cast<int>(words.size()) - 2;
            }
        }
        return 0;
    }

    static std::vector<int> buildListOfNodeIds(const std::vector<aux_node_t>& network_nodes) {
        std::vector<int> res;
        for (const auto& node : network_nodes) {
            if (std::find(res.begin(), res.end(), node.cId) == res.end()) {
                res.push_back(node.cId);
            }
        }
        return res;
    }

    static std::vector<aux_node_t> filterNetworkNodesByCoordinate(
        const std::vector<aux_node_t>& aux_nodes, const Mesh::coordinate_t& network_coordinate) {
        std::vector<aux_node_t> res;
        for (const auto& node : aux_nodes) {
            if (node.relPos == network_coordinate) {
                res.push_back(node);
            }
        }
        return res;
    }

    static std::vector<aux_node_t> filterNetworkNodesById(
        const std::vector<aux_node_t>& aux_nodes, int cId) {
        std::vector<aux_node_t> res;
        for (const auto& node : aux_nodes) {
            if (node.cId == cId) {
                res.push_back(node);
            }
        }
        return res;
    }

    static std::vector<aux_node_t> filterNetworkNodesByNetworkCircuit(const std::vector<aux_node_t>& aux_nodes) {
        std::vector<aux_node_t> res;
        for (const auto& node : aux_nodes) {
            if (node.node.termination.termination_type == Mtln::TERMINATION_NETWORK) {
                res.push_back(node);
            }
        }
        return res;
    }

    static Mtln::terminal_connection_t buildConnection(
        int node_id,
        const std::vector<aux_node_t>& network_nodes,
        const std::vector<Mtln::network_circuit_t>& network_circuits) {
        Mtln::terminal_connection_t res;
        for (const auto& node : network_nodes) {
            if (node.cId == node_id) {
                res.add_node(node.node);
            }
        }
        for (const auto& circuit : network_circuits) {
            if (circuit.nodeId == node_id) {
                res.network_circuit = circuit;
            }
        }
        return res;
    }

    void updateListOfNetworksCoordinates(std::vector<Mesh::coordinate_t>& coordinates, int conductor_index) {
        bool plFound = false;
        auto polyline = p.mesh.getPolyline(conductor_index, plFound);
        if (!plFound || polyline.coordIds.empty()) {
            return;
        }
        bool cFound = false;
        auto coord_ini = p.mesh.getCoordinate(polyline.coordIds.front(), cFound);
        auto coord_end = p.mesh.getCoordinate(polyline.coordIds.back(), cFound);

        bool found_ini = false;
        bool found_end = false;
        for (const auto& c : coordinates) {
            if (c == coord_ini) {
                found_ini = true;
            }
            if (c == coord_end) {
                found_end = true;
            }
        }
        if (!found_ini) {
            coordinates.push_back(coord_ini);
        }
        if (!found_end) {
            coordinates.push_back(coord_end);
        }
    }

    const nlohmann::json* getTerminationsOnSide(int terminationId) {
        if (terminationId == -1) {
            Report::WarnErrReport("Error: missing terminal on cable side", true);
            return nullptr;
        }
        const auto* terminal = materialById(terminationId);
        if (!terminal || !terminal->contains(jlbl::J_MAT_TERM_TERMINATIONS) ||
            !(*terminal)[jlbl::J_MAT_TERM_TERMINATIONS].is_array()) {
            Report::WarnErrReport("Error: missing terminations on terminal", true);
            return nullptr;
        }
        return &(*terminal)[jlbl::J_MAT_TERM_TERMINATIONS];
    }

    aux_node_t buildNode(
        const nlohmann::json* termination_list,
        int label,
        int index,
        int id,
        bool isShieldedCable) {
        aux_node_t res;
        if (!termination_list || !termination_list->is_array()) {
            return res;
        }
        if (index < 1 || static_cast<size_t>(index) > termination_list->size()) {
            return res;
        }
        const auto& termination = (*termination_list)[static_cast<size_t>(index - 1)];
        res.node.termination.termination_type = readTerminationType(termination);
        res.node.termination.capacitance = readTerminationRLC(termination, jlbl::J_MAT_TERM_CAPACITANCE, 1e22);
        res.node.termination.resistance = readTerminationRLC(termination, jlbl::J_MAT_TERM_RESISTANCE, 0.0);
        res.node.termination.inductance = readTerminationRLC(termination, jlbl::J_MAT_TERM_INDUCTANCE, 0.0);
        res.node.termination.source = readGeneratorOnTermination(id, label);
        res.node.termination.model = readTerminationModel(termination);
        res.node.termination.networkCircuitNode = readTerminationNetworkCircuitNode(termination, -1);
        res.node.side = label;
        res.node.conductor_in_cable = index;

        auto it = elemIdToCable.find(id);
        if (it != elemIdToCable.end()) {
            int cable_index = it->second;
            res.node.belongs_to_cable = mtln_res.cables[static_cast<size_t>(cable_index - 1)].ptr.get();
            bool plFound = false;
            auto polyline = p.mesh.getPolyline(id, plFound);
            if (plFound && !polyline.coordIds.empty()) {
                if (label == Mtln::TERMINAL_NODE_SIDE_INI) {
                    res.cId = polyline.coordIds.front();
                    res.relPos = p.mesh.getCoordinate(polyline.coordIds.front(), plFound);
                } else if (label == Mtln::TERMINAL_NODE_SIDE_END) {
                    res.cId = polyline.coordIds.back();
                    res.relPos = p.mesh.getCoordinate(polyline.coordIds.back(), plFound);
                }
                if (res.node.termination.termination_type == Mtln::TERMINATION_SHORT && !isShieldedCable) {
                    if (!terminalTouchesAnyEntity(res.cId, res.relPos, id)) {
                        res.node.termination.termination_type = Mtln::TERMINATION_OPEN;
                        std::string warningMsg =
                            "MTLN terminal on cable " + res.node.belongs_to_cable->name +
                            " (conductor " + std::to_string(index) + ", side " + sideToStr(label) +
                            ") is short but not touching any wire or non-vacuum material. Treating as open.";
                        Report::WarnErrReport(warningMsg, false);
                    }
                }
            }
        }
        return res;
    }

    bool terminalTouchesAnyEntity(int cId, const Mesh::coordinate_t& relPos, int ownElemId) {
        return touchesOtherWire(cId, ownElemId) || touchesNonVacuumMaterial(cId, relPos);
    }

    bool touchesOtherWire(int cId, int ownElemId) {
        auto wireMAs = p.getMaterialAssociations(
            {jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_WIRE},
            {});
        for (const auto& mA : wireMAs) {
            for (int elemId : mA.elementIds) {
                if (elemId == ownElemId) {
                    continue;
                }
                bool found = false;
                auto pl = p.mesh.getPolyline(elemId, found);
                if (found) {
                    for (int cid : pl.coordIds) {
                        if (cid == cId) {
                            return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    bool touchesNonVacuumMaterial(int cId, const Mesh::coordinate_t& relPos) {
        int ix = static_cast<int>(std::lround(relPos.position[0]));
        int iy = static_cast<int>(std::lround(relPos.position[1]));
        int iz = static_cast<int>(std::lround(relPos.position[2]));

        const auto* allMatAss = materialAssociationsArray();
        if (!allMatAss) {
            return false;
        }

        for (const auto& mA : *allMatAss) {
            const int materialId = jsonInt(mA, jlbl::J_MATERIAL_ID, 0);
            const std::string matType = materialTypeById(materialId);
            if (matType == jlbl::J_MAT_TYPE_WIRE || matType == jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE ||
                matType == jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE || matType == jlbl::J_MAT_TYPE_TERMINAL ||
                matType == jlbl::J_MAT_TYPE_CONNECTOR) {
                continue;
            }
            const auto* mat = materialById(materialId);
            if (!mat) {
                continue;
            }
            if (matType == jlbl::J_MAT_TYPE_ISOTROPIC && isVacuumIsotropic(*mat)) {
                continue;
            }
            for (int elemId : jsonIntVector(mA, jlbl::J_ELEMENTIDS)) {
                if (elementTouchesCoordinate(elemId, cId, ix, iy, iz)) {
                    return true;
                }
            }
        }
        return false;
    }

    bool elementTouchesCoordinate(int elemId, int cId, int ix, int iy, int iz) {
        bool found = false;
        auto node = p.mesh.getNode(elemId, found);
        if (found) {
            for (int cid : node.coordIds) {
                if (cid == cId) {
                    return true;
                }
            }
            return false;
        }
        auto pl = p.mesh.getPolyline(elemId, found);
        if (found) {
            for (int cid : pl.coordIds) {
                if (cid == cId) {
                    return true;
                }
            }
            return false;
        }
        auto cr = p.mesh.getCellRegion(elemId, found);
        if (found) {
            for (const auto& interval : cr.intervals) {
                if (intervalContainsNode(interval, ix, iy, iz)) {
                    return true;
                }
            }
        }
        return false;
    }

    static bool intervalContainsNode(const Cell::cell_interval_t& interval, int ix, int iy, int iz) {
        int ax = std::min(interval.ini.cell[0], interval.end.cell[0]);
        int bx = std::max(interval.ini.cell[0], interval.end.cell[0]);
        int ay = std::min(interval.ini.cell[1], interval.end.cell[1]);
        int by = std::max(interval.ini.cell[1], interval.end.cell[1]);
        int az = std::min(interval.ini.cell[2], interval.end.cell[2]);
        int bz = std::max(interval.ini.cell[2], interval.end.cell[2]);
        return ix >= ax && ix <= bx && iy >= ay && iy <= by && iz >= az && iz <= bz;
    }

    bool isVacuumIsotropic(const nlohmann::json& mat) {
        constexpr double tol = 1.0e-12;
        double relEps = jsonReal(mat, jlbl::J_MAT_REL_PERMITTIVITY, 1.0);
        double relMu = jsonReal(mat, jlbl::J_MAT_REL_PERMEABILITY, 1.0);
        double sigmaE = jsonReal(mat, jlbl::J_MAT_ELECTRIC_CONDUCTIVITY, 0.0);
        double sigmaM = jsonReal(mat, jlbl::J_MAT_MAGNETIC_CONDUCTIVITY, 0.0);
        double absEps = jsonReal(mat, jlbl::J_MAT_ABS_PERMITTIVITY, relEps * NFDE::EPSILON_VACUUM);
        double absMu = jsonReal(mat, jlbl::J_MAT_ABS_PERMEABILITY, relMu * NFDE::MU_VACUUM);
        return std::abs(relEps - 1.0) <= tol && std::abs(relMu - 1.0) <= tol &&
               std::abs(absEps - NFDE::EPSILON_VACUUM) <= std::max(tol, tol * NFDE::EPSILON_VACUUM) &&
               std::abs(absMu - NFDE::MU_VACUUM) <= std::max(tol, tol * NFDE::MU_VACUUM) &&
               std::abs(sigmaE) <= tol && std::abs(sigmaM) <= tol;
    }

    static std::string sideToStr(int side) {
        if (side == Mtln::TERMINAL_NODE_SIDE_INI) {
            return "initial";
        }
        if (side == Mtln::TERMINAL_NODE_SIDE_END) {
            return "end";
        }
        return "undefined";
    }

    Mtln::node_source_t readGeneratorOnTermination(int id, int label) {
        Mtln::node_source_t res;
        if (!p.rootJson.contains(jlbl::J_SOURCES) || !p.rootJson[jlbl::J_SOURCES].is_array()) {
            return res;
        }

        bool plFound = false;
        auto poly = p.mesh.getPolyline(id, plFound);
        if (!plFound) {
            return res;
        }
        for (const auto& gen : p.rootJson[jlbl::J_SOURCES]) {
            if (jsonString(gen, jlbl::J_TYPE) != jlbl::J_SRC_TYPE_GEN) {
                continue;
            }
            if (!gen.contains(jlbl::J_SRC_MAGNITUDE_FILE) || !gen[jlbl::J_SRC_MAGNITUDE_FILE].is_string()) {
                Report::WarnErrReport("magnitudeFile of source missing", true);
                return res;
            }
            if (!gen.contains(jlbl::J_FIELD) || !gen[jlbl::J_FIELD].is_string()) {
                Report::WarnErrReport("Type of generator is ambigous", true);
                return res;
            }
            std::string field = gen[jlbl::J_FIELD].get<std::string>();
            if (field != jlbl::J_FIELD_VOLTAGE && field != jlbl::J_FIELD_CURRENT) {
                Report::WarnErrReport("Only voltage and current generators are supported", true);
                return res;
            }
            if (isSourceAttachedToLine(gen, poly, id, label)) {
                if (field == jlbl::J_FIELD_VOLTAGE) {
                    res.source_type = Mtln::SOURCE_TYPE_VOLTAGE;
                    res.resistance = jsonReal(gen, jlbl::J_SRC_RESISTANCE_GEN, 0.0);
                } else {
                    res.source_type = Mtln::SOURCE_TYPE_CURRENT;
                    res.resistance = jsonReal(gen, jlbl::J_SRC_RESISTANCE_GEN, 1.0e22);
                }
                res.path_to_excitation = gen[jlbl::J_SRC_MAGNITUDE_FILE].get<std::string>();
                return res;
            }
        }
        return res;
    }

    bool isSourceAttachedToLine(
        const nlohmann::json& src, const Mesh::polyline_t& polyline, int id, int label) {
        if (!src.contains(jlbl::J_ELEMENTIDS) || !src[jlbl::J_ELEMENTIDS].is_array() ||
            src[jlbl::J_ELEMENTIDS].empty() || !src[jlbl::J_ELEMENTIDS][0].is_number_integer()) {
            return false;
        }
        std::vector<int> sourceElemIds;
        for (const auto& v : src[jlbl::J_ELEMENTIDS]) {
            if (v.is_number_integer()) sourceElemIds.push_back(v.get<int>());
        }
        if (sourceElemIds.empty()) return false;
        bool nodeFound = false;
        auto srcCoord = p.mesh.getNode(sourceElemIds[0], nodeFound);
        if (!nodeFound || srcCoord.coordIds.empty() || polyline.coordIds.empty()) return false;
        int index = (label == Mtln::TERMINAL_NODE_SIDE_INI) ? 0 : static_cast<int>(polyline.coordIds.size()) - 1;
        if (src.contains(jlbl::J_SRC_ATTACHED_ID) && src[jlbl::J_SRC_ATTACHED_ID].is_number_integer()) {
            return !srcCoord.coordIds.empty() && srcCoord.coordIds[0] == polyline.coordIds[static_cast<size_t>(index)] &&
                   src[jlbl::J_SRC_ATTACHED_ID].get<int>() == id;
        }
        return !srcCoord.coordIds.empty() && srcCoord.coordIds[0] == polyline.coordIds[static_cast<size_t>(index)];
    }

    static int readTerminationType(const nlohmann::json& termination) {
        std::string type = jsonString(termination, jlbl::J_TYPE);
        if (type == jlbl::J_MAT_TERM_TYPE_OPEN) return Mtln::TERMINATION_OPEN;
        if (type == jlbl::J_MAT_TERM_TYPE_SHORT) return Mtln::TERMINATION_SHORT;
        if (type == jlbl::J_MAT_TERM_TYPE_SERIES) return Mtln::TERMINATION_SERIES;
        if (type == jlbl::J_MAT_TERM_TYPE_PARALLEL) return Mtln::TERMINATION_PARALLEL;
        if (type == jlbl::J_MAT_TERM_TYPE_RsLCp) return Mtln::TERMINATION_RsLCp;
        if (type == jlbl::J_MAT_TERM_TYPE_LsRCp) return Mtln::TERMINATION_LsRCp;
        if (type == jlbl::J_MAT_TERM_TYPE_CsLRp) return Mtln::TERMINATION_CsLRp;
        if (type == jlbl::J_MAT_TERM_TYPE_RCsLp) return Mtln::TERMINATION_RCsLp;
        if (type == jlbl::J_MAT_TERM_TYPE_LCsRp) return Mtln::TERMINATION_LCsRp;
        if (type == jlbl::J_MAT_TERM_TYPE_RLsCp) return Mtln::TERMINATION_RLsCp;
        if (type == jlbl::J_MAT_TERM_TYPE_CIRCUIT) return Mtln::TERMINATION_CIRCUIT;
        if (type == jlbl::J_MAT_TERM_TYPE_NETWORK) return Mtln::TERMINATION_NETWORK;
        return Mtln::TERMINATION_UNDEFINED;
    }

    static Mtln::terminal_circuit_t readTerminationModel(const nlohmann::json& termination) {
        Mtln::terminal_circuit_t res;
        if (termination.contains(jlbl::J_MAT_TERM_MODEL_FILE) &&
            termination[jlbl::J_MAT_TERM_MODEL_FILE].is_string()) {
            res.file = termination[jlbl::J_MAT_TERM_MODEL_FILE].get<std::string>();
        }
        if (termination.contains(jlbl::J_MAT_TERM_MODEL_NAME) &&
            termination[jlbl::J_MAT_TERM_MODEL_NAME].is_string()) {
            res.name = termination[jlbl::J_MAT_TERM_MODEL_NAME].get<std::string>();
        }
        return res;
    }

    static int readTerminationNetworkCircuitNode(const nlohmann::json& termination, int defaultVal) {
        if (termination.contains(jlbl::J_MAT_TERM_MODEL_NODE) &&
            termination[jlbl::J_MAT_TERM_MODEL_NODE].is_number_integer()) {
            return termination[jlbl::J_MAT_TERM_MODEL_NODE].get<int>();
        }
        return defaultVal;
    }

    static double readTerminationRLC(
        const nlohmann::json& termination, const std::string& label, double defaultVal) {
        if (termination.contains(label) && termination[label].is_number()) {
            return termination[label].get<double>();
        }
        return defaultVal;
    }

    double readTerminationRLC(const nlohmann::json& termination, const char* label, double defaultVal) {
        return readTerminationRLC(termination, std::string(label), defaultVal);
    }

    std::vector<Mtln::parsed_generator_t> readWireGenerators() {
        const auto* sources = sourcesArray();
        if (!sources) {
            return {};
        }
        std::vector<Mtln::parsed_generator_t> res;
        for (const auto& gen : *sources) {
            if (jsonString(gen, jlbl::J_TYPE) != jlbl::J_SRC_TYPE_GEN) {
                continue;
            }
            if (!isGeneratorOnWire(gen)) {
                continue;
            }
            if (!gen.contains(jlbl::J_SRC_MAGNITUDE_FILE) || !gen[jlbl::J_SRC_MAGNITUDE_FILE].is_string()) {
                Report::WarnErrReport("magnitudeFile of source missing", true);
            }
            Mtln::parsed_generator_t g;
            std::string field = jsonString(gen, jlbl::J_FIELD);
            if (field == jlbl::J_FIELD_VOLTAGE) {
                g.generator_type = Mtln::SOURCE_TYPE_VOLTAGE;
                g.resistance = jsonReal(gen, jlbl::J_SRC_RESISTANCE_GEN, 0.0);
            } else if (field == jlbl::J_FIELD_CURRENT) {
                g.generator_type = Mtln::SOURCE_TYPE_CURRENT;
                g.resistance = jsonReal(gen, jlbl::J_SRC_RESISTANCE_GEN, 1.0e22);
            } else {
                Report::WarnErrReport(
                    "Field block of source of type generator must be current or voltage", true);
            }
            g.path_to_excitation = jsonString(gen, jlbl::J_SRC_MAGNITUDE_FILE);
            auto idAndPos = getPolylineElemIdAndConductorOfGenerator(gen);
            auto it = elemIdToCable.find(idAndPos.first);
            if (it == elemIdToCable.end()) {
                continue;
            }
            int index = it->second;
            auto coord = getCoordinateFromElemIdNode(gen);
            bool plFound = false;
            auto pl = p.mesh.getPolyline(idAndPos.first, plFound);
            if (!plFound) {
                continue;
            }
            auto linels = p.mesh.polylineToLinels(pl);
            g.conductor = idAndPos.second;
            g.index = findIndexInLinels(coord, linels);
            g.attached_to_cable = mtln_res.cables[static_cast<size_t>(index - 1)].ptr.get();
            res.push_back(g);
        }
        return res;
    }

    bool isGeneratorOnWire(const nlohmann::json& src) {
        std::string fieldLabel = jsonString(src, jlbl::J_FIELD);
        if (fieldLabel.empty() ||
            (fieldLabel != jlbl::J_FIELD_CURRENT && fieldLabel != jlbl::J_FIELD_VOLTAGE)) {
            Report::WarnErrReport("field type not recognized", true);
            return false;
        }
        auto eIds = jsonIntVector(src, jlbl::J_ELEMENTIDS);
        if (eIds.empty()) return false;
        auto pixel = Pt::getPixelFromElementId(p.mesh, eIds[0]);
        int cId = pixel.tag;

        auto mAs = p.getMaterialAssociations(
            {jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_WIRE},
            {});
        for (const auto& mA : mAs) {
            for (int elemId : mA.elementIds) {
                bool plFound = false;
                auto polyline = p.mesh.getPolyline(elemId, plFound);
                if (!plFound || polyline.coordIds.size() < 2) {
                    continue;
                }
                // Match Fortran IsGeneratorOnWire: only interior polyline nodes (j=2..n-1).
                for (size_t j = 1; j + 1 < polyline.coordIds.size(); ++j) {
                    if (polyline.coordIds[j] != cId) {
                        continue;
                    }
                    const bool interior = true;
                    if (interior) {
                        if (fieldLabel == jlbl::J_FIELD_VOLTAGE &&
                            (mA.matAssType == jlbl::J_MAT_TYPE_WIRE ||
                             mA.matAssType == jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE)) {
                            Report::WarnErrReport(
                                "Voltage generators cannot be defined on wire/unshieldedMultiwire interior points",
                                true);
                        } else if (fieldLabel == jlbl::J_FIELD_CURRENT &&
                                   mA.matAssType == jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE) {
                            Report::WarnErrReport(
                                "Current generators cannot be defined on shieldedMultiwire interior points", true);
                        }
                    }
                    return true;
                }
            }
        }
        return false;
    }

    std::vector<Mtln::probe_t> readMultiwireProbes() {
        const auto* probesRoot = probesArray();
        if (!probesRoot) {
            return {};
        }
        std::vector<const nlohmann::json*> wire_probes;
        for (const auto& probe : *probesRoot) {
            if (jsonString(probe, jlbl::J_TYPE) == jlbl::J_PR_TYPE_WIRE) {
                wire_probes.push_back(&probe);
            }
        }
        int n = countNumberOfMultiwireProbes(wire_probes);
        std::vector<Mtln::probe_t> res;
        res.reserve(static_cast<size_t>(n));
        for (const auto* probePtr : wire_probes) {
            if (!probePtr || !isProbeDefinedOnMultiwire(*probePtr)) {
                continue;
            }
            auto ids = getPolylineElemIdOfMultiwireProbe(*probePtr);
            auto probe_node_coord = getCoordinateFromElemIdNode(*probePtr);
            for (int elemId : ids) {
                Mtln::probe_t probe;
                probe.probe_name = readProbeName(*probePtr);
                probe.probe_type = readProbeType(*probePtr);
                probe.probe_position = {
                    probe_node_coord.position[0],
                    probe_node_coord.position[1],
                    probe_node_coord.position[2],
                };
                auto it = elemIdToCable.find(elemId);
                if (it == elemIdToCable.end()) {
                    continue;
                }
                int index = it->second;
                bool plFound = false;
                auto pl = p.mesh.getPolyline(elemId, plFound);
                auto linels = p.mesh.polylineToLinels(pl);
                probe.index = findIndexInLinels(probe_node_coord, linels);
                Mtln::cable_t* cable_ptr = mtln_res.cables[static_cast<size_t>(index - 1)].ptr.get();
                bool parent_cable_found = false;
                while (!parent_cable_found) {
                    auto* shielded = dynamic_cast<Mtln::shielded_multiwire_t*>(cable_ptr);
                    if (shielded) {
                        if (shielded->parent_cable) {
                            cable_ptr = shielded->parent_cable;
                        } else {
                            parent_cable_found = true;
                        }
                    } else {
                        parent_cable_found = true;
                    }
                }
                probe.attached_to_cable = cable_ptr;
                res.push_back(probe);
            }
        }
        return res;
    }

    Mesh::coordinate_t getCoordinateFromElemIdNode(const nlohmann::json& object) {
        auto elemIds = jsonIntVector(object, jlbl::J_ELEMENTIDS);
        bool found = false;
        if (elemIds.empty()) {
            return Mesh::coordinate_t{};
        }
        auto node = p.mesh.getNode(elemIds[0], found);
        if (!found || node.coordIds.empty()) {
            return Mesh::coordinate_t{};
        }
        return p.mesh.getCoordinate(node.coordIds[0], found);
    }

    static int findIndexInLinels(const Mesh::coordinate_t& coord, const std::vector<Mesh::linel_t>& linels) {
        std::vector<Mesh::coordinate_t> linelCoords(linels.size() + 1);
        for (size_t i = 0; i < linels.size(); ++i) {
            linelCoords[i].position[0] = linels[i].cell[0];
            linelCoords[i].position[1] = linels[i].cell[1];
            linelCoords[i].position[2] = linels[i].cell[2];
            int orAbs = std::abs(linels[i].orientation);
            if (orAbs > 0 && orAbs <= 3) {
                int idx = orAbs - 1;
                if (linels[i].orientation < 0) {
                    linelCoords[i].position[idx] += 1.0;
                }
            }
        }
        if (!linels.empty()) {
            int tailOr = linels.back().orientation;
            linelCoords[linels.size()] = linelCoords[linels.size() - 1];
            int tailOrAbs = std::abs(tailOr);
            if (tailOrAbs > 0 && tailOrAbs <= 3) {
                int idx = tailOrAbs - 1;
                linelCoords[linels.size()].position[idx] +=
                    (tailOr > 0 ? 1.0 : -1.0);
            }
        }
        int best = 1;
        double bestDist = 1e300;
        for (size_t i = 0; i < linelCoords.size(); ++i) {
            double dist = 0.0;
            for (int d = 0; d < 3; ++d) {
                double diff = linelCoords[i].position[d] - coord.position[d];
                dist += diff * diff;
            }
            if (dist < bestDist) {
                bestDist = dist;
                best = static_cast<int>(i) + 1;
            }
        }
        return best;
    }

    bool isProbeDefinedOnMultiwire(const nlohmann::json& probePtr) {
        std::string fieldLabel = jsonString(probePtr, jlbl::J_FIELD);
        if (fieldLabel.empty() ||
            (fieldLabel != jlbl::J_FIELD_CURRENT && fieldLabel != jlbl::J_FIELD_VOLTAGE)) {
            return false;
        }
        auto eIds = jsonIntVector(probePtr, jlbl::J_ELEMENTIDS);
        if (eIds.empty()) {
            return false;
        }
        auto pixel = Pt::getPixelFromElementId(p.mesh, eIds[0]);
        int cId = pixel.tag;
        auto mAs = p.getMaterialAssociations(
            {jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_WIRE},
            {});
        for (const auto& mA : mAs) {
            if (mA.elementIds.empty()) {
                continue;
            }
            bool plFound = false;
            auto polyline = p.mesh.getPolyline(mA.elementIds[0], plFound);
            if (!plFound) {
                continue;
            }
            for (int cid : polyline.coordIds) {
                if (cid == cId) {
                    return true;
                }
            }
        }
        return false;
    }

    int countNumberOfMultiwireProbes(const std::vector<const nlohmann::json*>& probes) {
        int count = 0;
        for (const auto* probePtr : probes) {
            if (probePtr && isProbeDefinedOnMultiwire(*probePtr)) {
                count += static_cast<int>(getPolylineElemIdOfMultiwireProbe(*probePtr).size());
            }
        }
        return count;
    }

    std::vector<int> getPolylineElemIdOfMultiwireProbe(const nlohmann::json& probePtr) {
        auto eIds = jsonIntVector(probePtr, jlbl::J_ELEMENTIDS);
        if (eIds.empty()) {
            return {};
        }
        auto pixel = Pt::getPixelFromElementId(p.mesh, eIds[0]);
        int cId = pixel.tag;
        std::vector<int> res;
        auto mAs = p.getMaterialAssociations(
            {jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_WIRE},
            {});
        for (const auto& mA : mAs) {
            if (mA.elementIds.empty()) {
                continue;
            }
            bool plFound = false;
            auto polyline = p.mesh.getPolyline(mA.elementIds[0], plFound);
            if (!plFound) {
                continue;
            }
            for (int cid : polyline.coordIds) {
                if (cid == cId) {
                    res.push_back(mA.elementIds[0]);
                }
            }
        }
        return res;
    }

    std::pair<int, int> getPolylineElemIdAndConductorOfGenerator(const nlohmann::json& src) {
        auto eIds = jsonIntVector(src, jlbl::J_ELEMENTIDS);
        if (eIds.empty()) {
            Report::WarnErrReport(
                "Generator does not define elementIds on wire, unshielded multiwire or shielded multiwire", true);
            return {0, 0};
        }
        auto pixel = Pt::getPixelFromElementId(p.mesh, eIds[0]);
        int cId = pixel.tag;
        std::pair<int, int> res{0, 0};
        auto mAs = p.getMaterialAssociations(
            {jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE,
             jlbl::J_MAT_TYPE_WIRE},
            {});
        for (const auto& mA : mAs) {
            for (size_t k = 0; k < mA.elementIds.size(); ++k) {
                bool plFound = false;
                auto polyline = p.mesh.getPolyline(mA.elementIds[k], plFound);
                if (!plFound || polyline.coordIds.size() < 2) {
                    continue;
                }
                for (size_t j = 0; j < polyline.coordIds.size(); ++j) {
                    if (polyline.coordIds[j] == cId) {
                        res.first = mA.elementIds[k];
                        res.second = static_cast<int>(k) + 1;
                        break;
                    }
                }
                if (res.first != 0) {
                    break;
                }
            }
            if (res.first != 0) {
                break;
            }
        }
        if (res.first == 0 && res.second == 0) {
            Report::WarnErrReport(
                "Generator does not belong to any wire, unshielded multiwire or shielded multiwire", true);
        }
        return res;
    }

    int readProbeType(const nlohmann::json& probePtr) {
        std::string probe_type = jsonString(probePtr, jlbl::J_FIELD);
        if (probe_type == jlbl::J_FIELD_VOLTAGE) {
            return Mtln::PROBE_TYPE_VOLTAGE;
        }
        if (probe_type == jlbl::J_FIELD_CURRENT) {
            return Mtln::PROBE_TYPE_CURRENT;
        }
        Report::WarnErrReport("probe type " + probe_type + " not supported", true);
        return Mtln::PROBE_TYPE_UNDEFINED;
    }

    std::string readProbeName(const nlohmann::json& probePtr) {
        if (probePtr.contains(jlbl::J_NAME) && probePtr[jlbl::J_NAME].is_string()) {
            return probePtr[jlbl::J_NAME].get<std::string>();
        }
        return "";
    }

    Mtln::cable_t* assignParentCable(const parser_t::materialAssociation_t& cable) {
        const std::string matType = materialTypeById(cable.materialId);
        if (matType == jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE) {
            if (cable.containedWithinElementId == -1) {
                return nullptr;
            }
            return getPointerToParentCable(cable.containedWithinElementId);
        }
        if (matType == jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE || matType == jlbl::J_MAT_TYPE_WIRE) {
            return nullptr;
        }
        Report::WarnErrReport("ERROR: Material type not recognized", true);
        return nullptr;
    }

    int assignConductorInParent(const parser_t::materialAssociation_t& cable) {
        const std::string matType = materialTypeById(cable.materialId);
        if (matType == jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE) {
            if (cable.containedWithinElementId == -1) {
                return 0;
            }
            return getParentPositionInMultiwire(cable.containedWithinElementId);
        }
        if (matType == jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE || matType == jlbl::J_MAT_TYPE_WIRE) {
            return 0;
        }
        Report::WarnErrReport("ERROR: Material type not recognized", true);
        return 0;
    }

    Mtln::cable_t* getPointerToParentCable(int id) {
        auto it = elemIdToCable.find(id);
        if (it == elemIdToCable.end()) {
            return nullptr;
        }
        return mtln_res.cables[static_cast<size_t>(it->second - 1)].ptr.get();
    }

    Mtln::connector_t* findConnectorWithId(int conn_id) {
        if (conn_id == -1) {
            return nullptr;
        }
        auto it = connIdToConnector.find(conn_id);
        if (it == connIdToConnector.end()) {
            return nullptr;
        }
        return &mtln_res.connectors[static_cast<size_t>(it->second - 1)];
    }

    int getParentPositionInMultiwire(int id) {
        auto it = elemIdToPosition.find(id);
        if (it == elemIdToPosition.end()) {
            return 0;
        }
        return it->second;
    }

    std::unique_ptr<Mtln::cable_t> readMTLNCable(const parser_t::materialAssociation_t& j_cable) {
        const nlohmann::json* material = materialById(j_cable.materialId);
        if (!material) {
            Report::WarnErrReport("Error reading cable: material not found", true);
            return std::make_unique<Mtln::unshielded_multiwire_t>();
        }
        std::string materialType = jsonString(*material, jlbl::J_TYPE);
        auto mtln_despl = buildMTLNDespl();
        auto cable_segments = buildSegments(j_cable, mtln_despl);
        auto cable_step_size = buildStepSize(cable_segments, mtln_despl);
        double totalLength = 0.0;
        for (double s : cable_step_size) {
            totalLength += s;
        }

        std::unique_ptr<Mtln::cable_t> res;
        if (materialType == jlbl::J_MAT_TYPE_SHIELDED_MULTIWIRE) {
            auto sh = std::make_unique<Mtln::shielded_multiwire_t>();
            sh->transfer_impedance = buildTransferImpedance(*material);
            assignPULProperties(*sh, *material, static_cast<int>(j_cable.elementIds.size()));
            if (j_cable.hasTotalResistance) {
                sh->resistance_per_meter =
                    Pt::vectorToDiagonalMatrix(scaleVector(j_cable.totalResistance, 1.0 / totalLength));
            }
            res = std::move(sh);
        } else if (materialType == jlbl::J_MAT_TYPE_UNSHIELDED_MULTIWIRE ||
                   materialType == jlbl::J_MAT_TYPE_WIRE) {
            auto unsh = std::make_unique<Mtln::unshielded_multiwire_t>();
            assignInCellProperties(*unsh, *material, static_cast<int>(j_cable.elementIds.size()));
            if (j_cable.hasTotalResistance) {
                unsh->resistance_per_meter =
                    Pt::vectorToDiagonalMatrix(scaleVector(j_cable.totalResistance, 1.0 / totalLength));
            }
            if (!j_cable.elementIds.empty()) {
                unsh->tag = std::to_string(j_cable.elementIds[0]);
            }
            res = std::move(unsh);
        } else {
            Report::WarnErrReport("Error reading cable: material type is not valid", true);
            res = std::make_unique<Mtln::unshielded_multiwire_t>();
        }

        res->initial_connector = findConnectorWithId(j_cable.initialConnectorId);
        res->end_connector = findConnectorWithId(j_cable.endConnectorId);
        res->name = j_cable.name;
        res->segments = std::move(cable_segments);
        res->n_segments = static_cast<int>(res->segments.size());
        res->step_size = std::move(cable_step_size);
        return res;
    }

    static std::vector<double> scaleVector(const std::vector<double>& v, double factor) {
        std::vector<double> out(v.size());
        for (size_t i = 0; i < v.size(); ++i) {
            out[i] = v[i] * factor;
        }
        return out;
    }

    NFDE::Desplazamiento_t buildMTLNDespl() {
        auto despl = p.readGrid();
        NFDE::Desplazamiento_t res;
        res.nX = despl.nX;
        res.nY = despl.nY;
        res.nZ = despl.nZ;
        res.mx1 = despl.mx1;
        res.my1 = despl.my1;
        res.mz1 = despl.mz1;
        res.mx2 = despl.mx2;
        res.my2 = despl.my2;
        res.mz2 = despl.mz2;
        copyAndEnlargeDes(res.desX, despl.desX, despl.mx2);
        copyAndEnlargeDes(res.desY, despl.desY, despl.my2);
        copyAndEnlargeDes(res.desZ, despl.desZ, despl.mz2);
        return res;
    }

    static void copyAndEnlargeDes(std::vector<double>& copy, const std::vector<double>& d, int n) {
        if (d.size() == 1) {
            copy.assign(static_cast<size_t>(n), d[0]);
        } else {
            copy = d;
        }
    }

    Mtln::transfer_impedance_per_meter_t buildTransferImpedance(const nlohmann::json& mat) {
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE)) {
            return readTransferImpedance(mat[jlbl::J_MAT_MULTIWIRE_TRANSFER_IMPEDANCE]);
        }
        return noTransferImpedance();
    }

    void assignPULProperties(Mtln::shielded_multiwire_t& res, const nlohmann::json& mat, int n) {
        std::vector<std::vector<double>> null_matrix(n, std::vector<double>(n, 0.0));
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_INDUCTANCE)) {
            res.inductance_per_meter = jsonMatrix(mat, jlbl::J_MAT_MULTIWIRE_INDUCTANCE);
        } else {
            Report::WarnErrReport("Error reading material region: inductancePerMeter label not found.", true);
            res.inductance_per_meter = null_matrix;
        }
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_CAPACITANCE)) {
            res.capacitance_per_meter = jsonMatrix(mat, jlbl::J_MAT_MULTIWIRE_CAPACITANCE);
        } else {
            Report::WarnErrReport("Error reading material region: capacitancePerMeter label not found.", true);
            res.capacitance_per_meter = null_matrix;
        }
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_RESISTANCE)) {
            res.resistance_per_meter =
                Pt::vectorToDiagonalMatrix(jsonRealVector(mat, jlbl::J_MAT_MULTIWIRE_RESISTANCE));
        } else {
            res.resistance_per_meter = null_matrix;
        }
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_CONDUCTANCE)) {
            res.conductance_per_meter =
                Pt::vectorToDiagonalMatrix(jsonRealVector(mat, jlbl::J_MAT_MULTIWIRE_CONDUCTANCE));
        } else {
            res.conductance_per_meter = null_matrix;
        }
    }

    void assignInCellProperties(Mtln::unshielded_multiwire_t& res, const nlohmann::json& mat, int n) {
        std::vector<std::vector<double>> null_matrix(n, std::vector<double>(n, 0.0));
        bool areFixedInCell = mat.contains(jlbl::J_MAT_MULTIWIRE_INDUCTANCE) &&
                              mat.contains(jlbl::J_MAT_MULTIWIRE_CAPACITANCE);
        bool areMultipolarInCell = mat.contains(jlbl::J_MAT_MULTIWIRE_MULTIPOLAR_EXPANSION);
        bool hasRadius = mat.contains(jlbl::J_MAT_WIRE_RADIUS) &&
                         jsonReal(mat, jlbl::J_MAT_WIRE_RADIUS, 0.0) != 0.0;
        if (!hasRadius) {
            if ((areFixedInCell && areMultipolarInCell) || (!areFixedInCell && !areMultipolarInCell)) {
                Report::WarnErrReport(
                    "Unshielded multiwires in cell properties must be defined by fixed OR multipolarExpansions, but not both.",
                    true);
            }
        }
        if (areFixedInCell) {
            res.cell_inductance_per_meter = jsonMatrix(mat, jlbl::J_MAT_MULTIWIRE_INDUCTANCE);
            res.cell_capacitance_per_meter = jsonMatrix(mat, jlbl::J_MAT_MULTIWIRE_CAPACITANCE);
            res.multipolar_expansion.clear();
        } else if (areMultipolarInCell) {
            res.cell_inductance_per_meter = null_matrix;
            res.cell_capacitance_per_meter = null_matrix;
            res.multipolar_expansion.resize(1);
            res.multipolar_expansion[0] = readMultipolarExpansion(mat[jlbl::J_MAT_MULTIWIRE_MULTIPOLAR_EXPANSION]);
        } else if (hasRadius) {
            res.cell_inductance_per_meter = null_matrix;
            res.cell_capacitance_per_meter = null_matrix;
            res.multipolar_expansion.clear();
            res.radius = jsonReal(mat, jlbl::J_MAT_WIRE_RADIUS, 0.0);
        }
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_RESISTANCE)) {
            int m = jsonDimension(mat, jlbl::J_MAT_MULTIWIRE_RESISTANCE);
            std::vector<double> r;
            if (m == 0) {
                r = {jsonReal(mat, jlbl::J_MAT_MULTIWIRE_RESISTANCE, 0.0)};
            } else {
                r = jsonRealVector(mat, jlbl::J_MAT_MULTIWIRE_RESISTANCE);
            }
            res.resistance_per_meter = Pt::vectorToDiagonalMatrix(r);
        } else {
            res.resistance_per_meter = null_matrix;
        }
        if (mat.contains(jlbl::J_MAT_MULTIWIRE_CONDUCTANCE)) {
            int m = jsonDimension(mat, jlbl::J_MAT_MULTIWIRE_CONDUCTANCE);
            std::vector<double> c;
            if (m == 0) {
                c = {jsonReal(mat, jlbl::J_MAT_MULTIWIRE_CONDUCTANCE, 0.0)};
            } else {
                c = jsonRealVector(mat, jlbl::J_MAT_MULTIWIRE_CONDUCTANCE);
            }
            res.conductance_per_meter = Pt::vectorToDiagonalMatrix(c);
        } else {
            res.conductance_per_meter = null_matrix;
        }
    }

    Mtln::multipolar_expansion_t readMultipolarExpansion(const nlohmann::json& multipolarExpansionPtr) {
        Mtln::multipolar_expansion_t res;
        if (!multipolarExpansionPtr.contains(jlbl::J_MAT_MULTIWIRE_ME_INNER_REGION_BOX)) {
            Report::WarnErrReport("Error reading multipolar expansion: innerRegionBox label not found", true);
        }
        res.inner_region = readInnerRegionBox(multipolarExpansionPtr[jlbl::J_MAT_MULTIWIRE_ME_INNER_REGION_BOX]);
        if (!multipolarExpansionPtr.contains(jlbl::J_MAT_MULTIWIRE_ME_ELECTRIC)) {
            Report::WarnErrReport("Error reading multipolar expansion electric reconstruction not found", true);
        }
        res.electric = readFieldReconstruction(multipolarExpansionPtr[jlbl::J_MAT_MULTIWIRE_ME_ELECTRIC]);
        if (!multipolarExpansionPtr.contains(jlbl::J_MAT_MULTIWIRE_ME_MAGNETIC)) {
            Report::WarnErrReport("Error reading multipolar expansion magnetic reconstruction not found", true);
        }
        res.magnetic = readFieldReconstruction(multipolarExpansionPtr[jlbl::J_MAT_MULTIWIRE_ME_MAGNETIC]);
        return res;
    }

    Mtln::box_2d_t readInnerRegionBox(const nlohmann::json& ptr) {
        Mtln::box_2d_t inner_region;
        inner_region.min = jsonRealVector(ptr, jlbl::J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MIN);
        inner_region.max = jsonRealVector(ptr, jlbl::J_MAT_MULTIWIRE_ME_INNER_REGION_BOX_MAX);
        return inner_region;
    }

    std::vector<Mtln::field_reconstruction_t> readFieldReconstruction(const nlohmann::json& ptr) {
        if (!ptr.is_array()) {
            return {};
        }
        int count = static_cast<int>(ptr.size());
        std::vector<Mtln::field_reconstruction_t> res(count);
        for (int j = 1; j <= count; ++j) {
            const auto& frPtr = ptr[static_cast<size_t>(j - 1)];
            res[static_cast<size_t>(j - 1)].inner_region_average_potential =
                jsonReal(frPtr, jlbl::J_MAT_MULTIWIRE_MEFR_INNER_REGION_AVERAGE_POTENTIAL, 0.0);
            res[static_cast<size_t>(j - 1)].expansion_center =
                jsonRealVector(frPtr, jlbl::J_MAT_MULTIWIRE_MEFR_EXPANSION_CENTER);
            res[static_cast<size_t>(j - 1)].conductor_potentials =
                jsonRealVector(frPtr, jlbl::J_MAT_MULTIWIRE_MEFR_CONDUCTOR_POTENTIALS);
            if (!frPtr.contains(jlbl::J_MAT_MULTIWIRE_MEFR_AB) ||
                !frPtr[jlbl::J_MAT_MULTIWIRE_MEFR_AB].is_array()) {
                Report::WarnErrReport("Error reading multipolar expansion: ab label not found", true);
                continue;
            }
            const auto& absPtr = frPtr[jlbl::J_MAT_MULTIWIRE_MEFR_AB];
            int abCount = static_cast<int>(absPtr.size());
            res[static_cast<size_t>(j - 1)].ab.resize(abCount);
            for (int i = 1; i <= abCount; ++i) {
                const auto& abPtr = absPtr[static_cast<size_t>(i - 1)];
                if (abPtr.is_array() && abPtr.size() >= 2 &&
                    abPtr[0].is_number() && abPtr[1].is_number()) {
                    res[static_cast<size_t>(j - 1)].ab[static_cast<size_t>(i - 1)].a = abPtr[0].get<double>();
                    res[static_cast<size_t>(j - 1)].ab[static_cast<size_t>(i - 1)].b = abPtr[1].get<double>();
                }
            }
        }
        return res;
    }

    std::vector<Mtln::segment_t> buildSegments(
        const parser_t::materialAssociation_t& j_cable, const NFDE::Desplazamiento_t& despl) {
        if (j_cable.elementIds.empty()) {
            return {};
        }
        bool plFound = false;
        auto temp = p.mesh.getPolyline(j_cable.elementIds[0], plFound);
        auto linels = p.mesh.polylineToLinels(temp);
        std::vector<Mtln::segment_t> res(linels.size());
        int prevOr = 0;
        for (size_t i = 0; i < linels.size(); ++i) {
            res[i].x = linels[i].cell[0];
            res[i].y = linels[i].cell[1];
            res[i].z = linels[i].cell[2];
            res[i].orientation = linels[i].orientation;
            if (i > 0 && prevOr == std::abs(res[i].orientation)) {
                res[i].dualBox = res[i - 1].dualBox;
                res[i].d1 = res[i - 1].d1;
                res[i].d2 = res[i - 1].d2;
            } else {
                switch (std::abs(res[i].orientation)) {
                case Cell::DIR_X:
                    res[i].dualBox = getDualBoxYZ(res[i], despl);
                    res[i].d1 = despl.desY[clip(res[i].y - 1, 0, static_cast<int>(despl.desY.size()) - 1)];
                    res[i].d2 = despl.desZ[clip(res[i].z - 1, 0, static_cast<int>(despl.desZ.size()) - 1)];
                    break;
                case Cell::DIR_Y:
                    res[i].dualBox = getDualBoxZX(res[i], despl);
                    res[i].d1 = despl.desZ[clip(res[i].z - 1, 0, static_cast<int>(despl.desZ.size()) - 1)];
                    res[i].d2 = despl.desX[clip(res[i].x - 1, 0, static_cast<int>(despl.desX.size()) - 1)];
                    break;
                case Cell::DIR_Z:
                    res[i].dualBox = getDualBoxXY(res[i], despl);
                    res[i].d1 = despl.desX[clip(res[i].x - 1, 0, static_cast<int>(despl.desX.size()) - 1)];
                    res[i].d2 = despl.desY[clip(res[i].y - 1, 0, static_cast<int>(despl.desY.size()) - 1)];
                    break;
                default:
                    break;
                }
            }
            prevOr = std::abs(res[i].orientation);
        }
        return res;
    }

    static int clip(int i, int lo, int hi) {
        return std::max(lo, std::min(i, hi));
    }

    static Mtln::box_2d_t getDualBoxYZ(const Mtln::segment_t& segment, const NFDE::Desplazamiento_t& despl) {
        Mtln::box_2d_t res;
        int y0 = clip(segment.y - 1, 0, static_cast<int>(despl.desY.size()) - 1);
        int y1 = clip(segment.y, 0, static_cast<int>(despl.desY.size()) - 1);
        int z0 = clip(segment.z - 1, 0, static_cast<int>(despl.desZ.size()) - 1);
        int z1 = clip(segment.z, 0, static_cast<int>(despl.desZ.size()) - 1);
        res.min = {-0.5 * despl.desY[static_cast<size_t>(y0)], -0.5 * despl.desZ[static_cast<size_t>(z0)]};
        res.max = {0.5 * despl.desY[static_cast<size_t>(y1)], 0.5 * despl.desZ[static_cast<size_t>(z1)]};
        return res;
    }

    static Mtln::box_2d_t getDualBoxXY(const Mtln::segment_t& segment, const NFDE::Desplazamiento_t& despl) {
        Mtln::box_2d_t res;
        int x0 = clip(segment.x - 1, 0, static_cast<int>(despl.desX.size()) - 1);
        int x1 = clip(segment.x, 0, static_cast<int>(despl.desX.size()) - 1);
        int y0 = clip(segment.y - 1, 0, static_cast<int>(despl.desY.size()) - 1);
        int y1 = clip(segment.y, 0, static_cast<int>(despl.desY.size()) - 1);
        res.min = {-0.5 * despl.desX[static_cast<size_t>(x0)], -0.5 * despl.desY[static_cast<size_t>(y0)]};
        res.max = {0.5 * despl.desX[static_cast<size_t>(x1)], 0.5 * despl.desY[static_cast<size_t>(y1)]};
        return res;
    }

    static Mtln::box_2d_t getDualBoxZX(const Mtln::segment_t& segment, const NFDE::Desplazamiento_t& despl) {
        Mtln::box_2d_t res;
        int z0 = clip(segment.z - 1, 0, static_cast<int>(despl.desZ.size()) - 1);
        int z1 = clip(segment.z, 0, static_cast<int>(despl.desZ.size()) - 1);
        int x0 = clip(segment.x - 1, 0, static_cast<int>(despl.desX.size()) - 1);
        int x1 = clip(segment.x, 0, static_cast<int>(despl.desX.size()) - 1);
        res.min = {-0.5 * despl.desZ[static_cast<size_t>(z0)], -0.5 * despl.desX[static_cast<size_t>(x0)]};
        res.max = {0.5 * despl.desZ[static_cast<size_t>(z1)], 0.5 * despl.desX[static_cast<size_t>(x1)]};
        return res;
    }

    static std::vector<double> buildStepSize(
        const std::vector<Mtln::segment_t>& segments, const NFDE::Desplazamiento_t& despl) {
        std::vector<double> res(segments.size());
        for (size_t i = 0; i < segments.size(); ++i) {
            int orAbs = std::abs(segments[i].orientation);
            switch (orAbs) {
            case Cell::DIR_X:
                res[i] = despl.desX[clip(segments[i].x, 0, static_cast<int>(despl.desX.size()) - 1)];
                break;
            case Cell::DIR_Y:
                res[i] = despl.desY[clip(segments[i].y, 0, static_cast<int>(despl.desY.size()) - 1)];
                break;
            case Cell::DIR_Z:
                res[i] = despl.desZ[clip(segments[i].z, 0, static_cast<int>(despl.desZ.size()) - 1)];
                break;
            default:
                res[i] = 0.0;
                break;
            }
        }
        return res;
    }

    Mtln::transfer_impedance_per_meter_t readTransferImpedance(const nlohmann::json& z) {
        Mtln::transfer_impedance_per_meter_t res;
        if (z.contains(jlbl::J_MAT_TRANSFER_IMPEDANCE_RESISTANCE) &&
            z[jlbl::J_MAT_TRANSFER_IMPEDANCE_RESISTANCE].is_number()) {
            res.resistive_term = z[jlbl::J_MAT_TRANSFER_IMPEDANCE_RESISTANCE].get<double>();
        }
        if (z.contains(jlbl::J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE) &&
            z[jlbl::J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE].is_number()) {
            res.inductive_term = z[jlbl::J_MAT_TRANSFER_IMPEDANCE_INDUCTANCE].get<double>();
        }
        std::string direction;
        if (z.contains(jlbl::J_MAT_TRANSFER_IMPEDANCE_DIRECTION) &&
            z[jlbl::J_MAT_TRANSFER_IMPEDANCE_DIRECTION].is_string()) {
            direction = z[jlbl::J_MAT_TRANSFER_IMPEDANCE_DIRECTION].get<std::string>();
        } else {
            Report::WarnErrReport(
                "Error reading material: direction of transferImpedancePerMeter missing", true);
        }
        if (direction == "inwards") {
            res.direction = Mtln::TRANSFER_IMPEDANCE_DIRECTION_INWARDS;
        } else if (direction == "outwards") {
            res.direction = Mtln::TRANSFER_IMPEDANCE_DIRECTION_OUTWARDS;
        } else if (direction == "both") {
            res.direction = Mtln::TRANSFER_IMPEDANCE_DIRECTION_BOTH;
        }
        if (z.contains(jlbl::J_MAT_TRANSFER_IMPEDANCE_POLES) &&
            z[jlbl::J_MAT_TRANSFER_IMPEDANCE_POLES].is_array()) {
            int n = 0;
            if (z.contains(jlbl::J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES) &&
                z[jlbl::J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES].is_number_integer()) {
                n = z[jlbl::J_MAT_TRANSFER_IMPEDANCE_NUMBER_POLES].get<int>();
            }
            res.poles.resize(n);
            res.residues.resize(n);
            const auto& poles = z[jlbl::J_MAT_TRANSFER_IMPEDANCE_POLES];
            const auto* residuesPtr = z.contains(jlbl::J_MAT_TRANSFER_IMPEDANCE_RESIDUES)
                                          ? &z[jlbl::J_MAT_TRANSFER_IMPEDANCE_RESIDUES]
                                          : nullptr;
            if (!residuesPtr || !residuesPtr->is_array()) {
                return res;
            }
            const auto& residues = *residuesPtr;
            for (int i = 1; i <= n; ++i) {
                const size_t idx = static_cast<size_t>(i - 1);
                if (idx < poles.size() && poles[idx].is_array() && poles[idx].size() >= 2 &&
                    poles[idx][0].is_number() && poles[idx][1].is_number()) {
                    res.poles[static_cast<size_t>(i - 1)] = {
                        poles[idx][0].get<double>(), poles[idx][1].get<double>()};
                }
                if (idx < residues.size() && residues[idx].is_array() && residues[idx].size() >= 2 &&
                    residues[idx][0].is_number() && residues[idx][1].is_number()) {
                    res.residues[static_cast<size_t>(i - 1)] = {
                        residues[idx][0].get<double>(), residues[idx][1].get<double>()};
                }
            }
        }
        return res;
    }

    static Mtln::transfer_impedance_per_meter_t noTransferImpedance() {
        Mtln::transfer_impedance_per_meter_t res;
        return res;
    }
};

parser_t::MtlnReader::MtlnReader(parser_t& parser, Mtln::mtln_t& mtlnRes)
    : p(parser), mtln_res(mtlnRes) {}

mtln_types_m::mtln_t parser_t::readMTLN() {
    mtln_types_m::mtln_t res;
    MtlnReader reader(*this, res);
    reader.read();
    return std::move(res);
}

} // namespace smbjson

#endif // CompileWithMTLN
