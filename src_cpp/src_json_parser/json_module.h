#ifndef JSON_MODULE_H
#define JSON_MODULE_H

#include <nlohmann/json.hpp>
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <stdexcept>

// Lightweight wrapper around nlohmann::json that mimics the json-fortran
// API used by smbjson.cpp, but with nlohmann idioms under the hood.

namespace json_module {

    // json_value — thin wrapper around nlohmann::json
    class json_value {
    public:
        nlohmann::json data;
    };

    // json_core — provides the typed get() interface
    class json_core {
    public:
        json_core() = default;

        // Navigate to a path and return the json_value at that path
        void get(const json_value& root, const std::string& path, json_value*& out, bool& found) const {
            out = nullptr;
            found = false;
            try {
                auto* j = &root.data;
                std::string current;
                for (char c : path) {
                    if (c == '.') {
                        if (!current.empty()) {
                            if (current[0] == '(') {
                                int idx = 0;
                                try { idx = std::stoi(current); } catch (...) {}
                                if (j->is_array() && idx >= 0 && idx < (int)j->size()) {
                                    j = &(*j)[idx];
                                } else { return; }
                            } else {
                                if (j->contains(current)) {
                                    j = &(*j)[current];
                                } else { return; }
                            }
                        }
                        current.clear();
                    } else {
                        current += c;
                    }
                }
                if (!current.empty()) {
                    if (current[0] == '(') {
                        int idx = 0;
                        try { idx = std::stoi(current); } catch (...) {}
                        if (j->is_array() && idx >= 0 && idx < (int)j->size()) {
                            j = &(*j)[idx];
                        } else { return; }
                    } else {
                        if (j->contains(current)) {
                            j = &(*j)[current];
                        } else { return; }
                    }
                }
                out = new json_value();
                out->data = *j;
                found = true;
            } catch (...) {}
        }

        void get(const json_value& root, const std::string& path, json_value*& out) const {
            bool found = false;
            get(root, path, out, found);
            if (!found) out = nullptr;
        }

        int count(const json_value* val) const {
            if (val == nullptr) return 0;
            if (val->data.is_array() || val->data.is_object()) return (int)val->data.size();
            return 0;
        }

        bool get_child(const json_value* parent, int index, json_value*& child) const {
            child = nullptr;
            if (parent == nullptr) return false;
            if (!parent->data.is_array()) return false;
            if (index < 1 || index > (int)parent->data.size()) return false;
            child = new json_value();
            child->data = parent->data[index - 1];
            return true;
        }
    };

    // json_file — file loading
    class json_file {
    public:
        json_file() : _failed(false), _initialized(false) {}

        void initialize() {
            _root = nlohmann::json::object();
            _initialized = true;
            _failed = false;
        }

        bool failed() const { return _failed; }

        void load(const std::string& filename) {
            if (!_initialized) { _failed = true; return; }
            std::ifstream f(filename);
            if (!f.is_open()) { _failed = true; return; }
            try { f >> _root; _failed = false; }
            catch (...) { _failed = true; }
        }

        void get_core(json_core& core) { core = json_core(); }

        void get(const std::string& /*path*/, json_value& out) const {
            out.data = _root;
        }

        const nlohmann::json& get_root() const { return _root; }

    private:
        nlohmann::json _root;
        bool _failed;
        bool _initialized;
    };

} // namespace json_module

#endif // JSON_MODULE_H
