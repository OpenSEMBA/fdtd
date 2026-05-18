#include <string>
#include <vector>
#include <iostream>
#include <cctype>
#include <algorithm>

// Forward declarations for external types used in this module
// These would typically come from the included modules: mesh_m, cells_m, json_module, json_kinds, NFDETypes_m

// Placeholder for types from mesh_m, cells_m, etc.
// In a real translation, these would be actual struct/class definitions from the respective headers.

// Placeholder for json_module types
struct json_value {
    // Placeholder
};

// Placeholder for json_kinds
using RKIND = double;
using rkind = double;

// Placeholder for NFDETypes_m types
struct coordinate_t {
    // Placeholder
};

struct terminal_node_t {
    // Placeholder
};

struct cell_region_t {
    std::vector<int> intervals; // Placeholder for cell_interval_t vector
    int getType() const { return 0; } // Placeholder
    std::vector<int> getIntervalsOfType(int type) const { return {}; } // Placeholder
};

struct cell_interval_t {
    struct {
        int cell(int dir) const { return 0; }
    } ini;
    struct {
        int cell(int dir) const { return 0; }
    } end;
    int getOrientation() const { return 0; }
};

struct coords_t {
    int Xi, Xe, Yi, Ye, Zi, Ze;
    int Or;
    std::string tag;
};

struct coords_scaled_t {
    int Xi, Xe, Yi, Ye, Zi, Ze;
    int Or;
    std::string tag;
    double xc, yc, zc;
};

struct pixel_t {
    // Placeholder
};

struct mesh_t {
    struct node_t {
        // Placeholder
    };
    node_t getNode(int id, bool& found) { found = false; return node_t(); }
    pixel_t nodeToPixel(const node_t& node) { return pixel_t(); }
};

struct node_t {
    // Placeholder
};

// Constants
constexpr int BUFSIZE = 256;
constexpr int iEx = 1;
constexpr int iEy = 2;
constexpr int iEz = 3;
constexpr int DIR_X = 1;
constexpr int DIR_Y = 2;
constexpr int DIR_Z = 3;

// Type for JSON value pointer
struct json_value_ptr_t {
    json_value* p;
};

// Auxiliary node type (conditional on CompileWithMTLN)
#ifdef CompileWithMTLN
struct aux_node_t {
    terminal_node_t node;
    int cId;
    coordinate_t relPos;
};
#endif

namespace parser_tools_m {

    // Helper function to count intervals of a specific type
    int countIntervalsOfType(const std::vector<int>& intervals, int cellType) {
        int count = 0;
        for (int interval : intervals) {
            if (interval == cellType) {
                count++;
            }
        }
        return count;
    }

    std::vector<int> getIntervalsInCellRegions(const std::vector<cell_region_t>& cellRegions, int cellType = -1) {
        int numberOfIntervals = 0;
        for (size_t i = 0; i < cellRegions.size(); ++i) {
            if (cellType != -1) {
                numberOfIntervals += countIntervalsOfType(cellRegions[i].intervals, cellType);
            } else {
                numberOfIntervals += cellRegions[i].intervals.size();
            }
        }

        std::vector<int> intervals;
        intervals.reserve(numberOfIntervals);
        int copiedIntervals = 0;
        for (size_t i = 0; i < cellRegions.size(); ++i) {
            std::vector<int> intervalsInRegion;
            if (cellType != -1) {
                intervalsInRegion = cellRegions[i].getIntervalsOfType(cellType);
            } else {
                intervalsInRegion = cellRegions[i].intervals;
            }
            for (size_t j = 0; j < intervalsInRegion.size(); ++j) {
                intervals.push_back(intervalsInRegion[j]);
                copiedIntervals++;
            }
        }
        return intervals;
    }

    std::vector<coords_t> cellRegionToCoords(const cell_region_t& cellRegion, int cellType = -1, const std::string& tag = "") {
        std::vector<int> intervals = getIntervalsInCellRegions({cellRegion}, cellType);
        std::vector<coords_t> cs;
        if (!tag.empty()) {
            cs = cellIntervalsToCoords(intervals, tag);
        } else {
            cs = cellIntervalsToCoords(intervals);
        }
        return cs;
    }

    std::vector<coords_scaled_t> coordsToScaledCoords(const std::vector<coords_t>& cs) {
        std::vector<coords_scaled_t> res(cs.size());
        for (size_t i = 0; i < cs.size(); ++i) {
            res[i].Xi = cs[i].Xi;
            res[i].Xe = cs[i].Xe;
            res[i].Yi = cs[i].Yi;
            res[i].Ye = cs[i].Ye;
            res[i].Zi = cs[i].Zi;
            res[i].Ze = cs[i].Ze;
            res[i].Or = cs[i].Or;
            res[i].tag = cs[i].tag;

            res[i].xc = 0.0;
            res[i].yc = 0.0;
            res[i].zc = 0.0;

            switch (cs[i].Or) {
                case iEx:
                    res[i].xc = 1.0;
                    break;
                case -iEx:
                    res[i].xc = -1.0;
                    break;
                case iEy:
                    res[i].yc = 1.0;
                    break;
                case -iEy:
                    res[i].yc = -1.0;
                    break;
                case iEz:
                    res[i].zc = 1.0;
                    break;
                case -iEz:
                    res[i].zc = -1.0;
                    break;
            }
        }
        return res;
    }

    void cellRegionsToScaledCoords(std::vector<coords_scaled_t>& res, const std::vector<cell_region_t>& cellRegions, const std::string& tag = "") {
        std::vector<int> intervals = getIntervalsInCellRegions(cellRegions, CELL_TYPE_LINEL);
        std::vector<coords_t> cs;
        if (!tag.empty()) {
            cs = cellIntervalsToCoords(intervals, tag);
        } else {
            cs = cellIntervalsToCoords(intervals);
        }
        std::vector<coords_scaled_t> scaledCoords = coordsToScaledCoords(cs);
        res = scaledCoords;
    }

    std::vector<coords_t> cellIntervalsToCoords(const std::vector<int>& ivls, const std::string& tag = "") {
        std::vector<coords_t> res(ivls.size());
        for (size_t i = 0; i < ivls.size(); ++i) {
            // Placeholder: In real code, ivls would be cell_interval_t objects
            // Here we assume ivls is a vector of indices or simplified data
            // For the sake of translation, we'll assume a helper function exists
            // or that cell_interval_t can be accessed via index.
            // Since the original code uses ivls(i)%getOrientation(), we need cell_interval_t objects.
            // Let's assume ivls is actually std::vector<cell_interval_t> in the real scenario.
            // But the signature says type(cell_interval_t), dimension(:), intent(in) :: ivls
            // So let's redefine the function signature to match the intent.
        }
        return res;
    }

    // Redefining cellIntervalsToCoords to take cell_interval_t objects
    std::vector<coords_t> cellIntervalsToCoords(const std::vector<cell_interval_t>& ivls, const std::string& tag = "") {
        std::vector<coords_t> res(ivls.size());
        for (size_t i = 0; i < ivls.size(); ++i) {
            res[i].Or = ivls[i].getOrientation();
            convertInterval(res[i].Xi, res[i].Xe, ivls[i], DIR_X);
            convertInterval(res[i].Yi, res[i].Ye, ivls[i], DIR_Y);
            convertInterval(res[i].Zi, res[i].Ze, ivls[i], DIR_Z);
            res[i].tag = tag;
        }
        return res;
    }

    void convertInterval(int& xi, int& xe, const cell_interval_t& interval, int dir) {
        int a = interval.ini.cell(dir);
        int b = interval.end.cell(dir);
        if (a < b) {
            xi = a;
            xe = b - 1;
        } else if (a == b) {
            xi = a;
            xe = b;
        } else {
            xi = b;
            xe = a - 1;
        }
    }

    pixel_t getPixelFromElementId(const mesh_t& mesh, int id) {
        bool nodeFound = false;
        mesh_t::node_t node = mesh.getNode(id, nodeFound);
        pixel_t res;
        if (nodeFound) {
            res = mesh.nodeToPixel(node);
        } else {
            std::cerr << "Error converting pixel. Node not found." << std::endl;
            // In Fortran, stop would terminate the program. In C++, we might throw or exit.
            exit(1);
        }
        return res;
    }

    std::vector<std::vector<double>> vectorToDiagonalMatrix(const std::vector<double>& vector) {
        int n = vector.size();
        std::vector<std::vector<double>> res(n, std::vector<double>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            res[i][i] = vector[i];
        }
        return res;
    }

    std::vector<std::vector<double>> scalarToMatrix(double scalar) {
        std::vector<std::vector<double>> res(1, std::vector<double>(1, 0.0));
        res[0][0] = scalar;
        return res;
    }

    std::vector<std::string> splitLineIntoWords(const std::string& line) {
        std::vector<std::string> words;
        int lenstr = line.length();
        int i = 0;
        int nwords = 0;

        // First pass: count words
        while (i < lenstr) {
            while (i < lenstr && (line[i] == ' ' || line[i] == '\t')) {
                i++;
            }
            if (i > lenstr) break;
            nwords++;
            while (i < lenstr && (line[i] != ' ' && line[i] != '\t')) {
                i++;
            }
        }

        if (nwords == 0) {
            return words;
        }

        words.resize(nwords);
        i = 0;
        int start = 0;
        int wordIdx = 0;
        while (i < lenstr) {
            while (i < lenstr && (line[i] == ' ' || line[i] == '\t')) {
                i++;
            }
            if (i > lenstr) break;
            start = i;
            while (i < lenstr && (line[i] != ' ' && line[i] != '\t')) {
                i++;
            }
            int wlen = i - start;
            words[wordIdx++] = line.substr(start, wlen);
        }
        return words;
    }

    std::string to_upper(const std::string& str) {
        std::string res = str;
        for (size_t i = 0; i < res.length(); ++i) {
            if (res[i] >= 'a' && res[i] <= 'z') {
                res[i] = res[i] - 32; // ASCII adjustment
            }
        }
        return res;
    }

} // namespace parser_tools_m