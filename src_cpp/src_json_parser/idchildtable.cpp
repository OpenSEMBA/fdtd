#include <unordered_map>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
namespace json_module { struct json_core {}; }
namespace fhash { struct fhash_tbl_t {}; }
namespace idchildtable_m {
    struct my_key_t { std::string name; int index; };
    struct IdChildTable_t {
        std::unordered_map<std::string, int> keys;
        IdChildTable_t() {}
        IdChildTable_t(const json_module::json_core*, const std::string&) {}
        bool hasKey(const std::string& k) const { return keys.count(k) > 0; }
        int getKeyIndex(const std::string& k) const { auto it = keys.find(k); return it != keys.end() ? it->second : -1; }
        int size() const { return keys.size(); }
        std::string getName(int i) const { for(auto& p : keys) { if(p.second == i) return p.first; } return ""; }
        void addKey(const std::string&, int) {}
    };
}
