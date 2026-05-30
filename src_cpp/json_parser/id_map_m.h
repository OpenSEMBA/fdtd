#ifndef ID_MAP_M_H
#define ID_MAP_M_H

#include <string>
#include <unordered_map>

#include <nlohmann/json.hpp>

#include "smbjson_labels_m.h"

namespace id_map_m {

using id_map_t = std::unordered_map<int, const nlohmann::json*>;

inline const nlohmann::json* getPathNode(const nlohmann::json* root,
                                         const std::string& path) {
    if (root == nullptr) return nullptr;
    if (path.empty()) return root;

    const nlohmann::json* current = root;
    size_t pos = 0;
    while (pos <= path.size()) {
        const size_t dot = path.find('.', pos);
        const std::string token = path.substr(
            pos, dot == std::string::npos ? std::string::npos : dot - pos);
        if (!token.empty()) {
            if (token.front() == '(' && token.back() == ')') {
                const std::string idxStr = token.substr(1, token.size() - 2);
                int idx = 0;
                try {
                    idx = std::stoi(idxStr);
                } catch (...) {
                    return nullptr;
                }
                if (!current->is_array() || idx < 1 ||
                    idx > static_cast<int>(current->size())) {
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

inline id_map_t buildIdMap(const nlohmann::json& root, const std::string& path) {
    id_map_t res;
    const nlohmann::json* entries = getPathNode(&root, path);
    if (entries == nullptr || !entries->is_array()) return res;

    const int numberOfEntries = static_cast<int>(entries->size());
    for (int i = 0; i < numberOfEntries; ++i) {
        const nlohmann::json* entry = &(*entries)[static_cast<size_t>(i)];
        if (entry == nullptr || !entry->is_object()) continue;
        if (!entry->contains(smbjson_labels_m::J_ID)) continue;

        int id = 0;
        (*entry)[smbjson_labels_m::J_ID].get_to(id);
        res[id] = entry;
    }
    return res;
}

inline const nlohmann::json* findById(const id_map_t& map, int id) {
    auto it = map.find(id);
    if (it == map.end()) return nullptr;
    return it->second;
}

inline bool containsId(const id_map_t& map, int id) {
    return map.find(id) != map.end();
}

}  // namespace id_map_m

#endif  // ID_MAP_M_H
