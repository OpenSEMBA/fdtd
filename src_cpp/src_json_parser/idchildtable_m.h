#ifndef IDCHILDTABLE_M_H
#define IDCHILDTABLE_M_H

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include "json_module.h"
#include "smbjson_labels_m.h"
#include "Report_m.h"

namespace idchildtable_m {

    using json_value_ptr_t = json_module::json_value*;

    class IdChildTable_t {
    private:
        std::unordered_map<int, json_value_ptr_t> idToChilds;

    public:
        static IdChildTable_t ctor(json_module::json_core& core, json_module::json_value& root, const std::string& path);

        int totalSize();
        int checkId(int id);
        json_value_ptr_t getId(int id);
    };

    inline IdChildTable_t IdChildTable_t::ctor(json_module::json_core& core, json_module::json_value& root, const std::string& path) {
        IdChildTable_t res;
        json_module::json_value* jentries = nullptr;
        bool found = false;

        core.get(root, path, jentries, found);
        if (!found) return res;

        int numberOfEntries = core.count(jentries);
        for (int i = 1; i <= numberOfEntries; ++i) {
            json_module::json_value* jentry = nullptr;
            core.get_child(jentries, i, jentry);
            if (jentry == nullptr) continue;

            int id = 0;
            bool idFound = false;
            jentry->data[smbjson_labels_m::J_ID].get_to(id);

            res.idToChilds[id] = jentry;
        }
        return res;
    }

    inline int IdChildTable_t::totalSize() {
        return static_cast<int>(idToChilds.size());
    }

    inline int IdChildTable_t::checkId(int id) {
        return (idToChilds.find(id) == idToChilds.end()) ? 1 : 0;
    }

    inline json_value_ptr_t IdChildTable_t::getId(int id) {
        auto it = idToChilds.find(id);
        if (it != idToChilds.end()) {
            return it->second;
        }
        return nullptr;
    }

} // namespace idchildtable_m

#endif // IDCHILDTABLE_M_H
