#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cctype>

// Forward declarations and includes based on Fortran use statements
// Assuming these types exist in the respective headers
#include "mesh_m.h"
#include "cells_m.h"
#include "json_module.h"
#include "json_kinds.h"
#include "NFDETypes_m.h"

// Constants
constexpr int BUFSIZE = 256; // Assumed size, typically defined in kinds or similar
using rkind = double;
using RKIND = double;

// Enumerations assumed from context
enum Orientation {
    iEx = 1,
    iEy = 2,
    iEz = 3,
    DIR_X = 1,
    DIR_Y = 2,
    DIR_Z = 3
};

namespace parser_tools_m {

    // Helper struct for JSON value pointer if needed, though often handled by json_module
    struct json_value_ptr_t {
        json_value* p = nullptr;
    };

#ifdef CompileWithMTLN
    struct aux_node_t {
        terminal_node_t node;
        int cId;
        coordinate_t relPos;
    };
#endif

    // Function: getIntervalsInCellRegions
    std::vector<cell_interval_t> getIntervalsInCellRegions(
        const std::vector<cell_region_t>& cellRegions,
        int cellType = -1) // -1 indicates not present
    {
        int numberOfIntervals = 0;
        for (size_t i = 0; i < cellRegions.size(); ++i) {
            if (cellType != -1) {
                // Assuming cell_region_t has a method or member to get intervals and check type
                // Fortran: count(cellRegions(i)%intervals%getType() == cellType)
                // C++ equivalent depends on implementation. Assuming a helper or loop.
                // For translation purposes, we assume cell_region_t provides access to intervals
                // and cell_interval_t has getType().
                const auto& intervals = cellRegions[i].getIntervals(); // Placeholder for access
                for (const auto& interval : intervals) {
                    if (interval.getType() == cellType) {
                        numberOfIntervals++;
                    }
                }
            } else {
                numberOfIntervals += static_cast<int>(cellRegions[i].getIntervals().size());
            }
        }

        std::vector<cell_interval_t> intervals(numberOfIntervals);
        int copiedIntervals = 0;
        for (size_t i = 0; i < cellRegions.size(); ++i) {
            std::vector<cell_interval_t> intervalsInRegion;
            if (cellType != -1) {
                // Assuming getIntervalsOfType exists
                intervalsInRegion = cellRegions[i].getIntervalsOfType(cellType);
            } else {
                intervalsInRegion = cellRegions[i].getIntervals();
            }
            
            for (const auto& interval : intervalsInRegion) {
                if (copiedIntervals < numberOfIntervals) {
                    intervals[copiedIntervals] = interval;
                }
                copiedIntervals++;
            }
        }
        return intervals;
    }

    // Function: cellRegionToCoords
    std::vector<coords_t> cellRegionToCoords(
        const cell_region_t& cellRegion,
        int cellType = -1,
        const std::string& tag = "")
    {
        std::vector<cell_interval_t> intervals = getIntervalsInCellRegions({cellRegion}, cellType);
        std::vector<coords_t> cs;
        if (!tag.empty()) {
            cs = cellIntervalsToCoords(intervals, tag);
        } else {
            cs = cellIntervalsToCoords(intervals);
        }
        return cs;
    }

    // Function: coordsToScaledCoords
    std::vector<coords_scaled_t> coordsToScaledCoords(
        const std::vector<coords_t>& cs)
    {
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
                default:
                    break;
            }
        }
        return res;
    }

    // Subroutine: cellRegionsToScaledCoords
    void cellRegionsToScaledCoords(
        std::vector<coords_scaled_t>& res,
        const std::vector<cell_region_t>& cellRegions,
        const std::string& tag = "")
    {
        std::vector<cell_interval_t> intervals = getIntervalsInCellRegions(cellRegions, CELL_TYPE_LINEL);
        std::vector<coords_t> cs;
        if (!tag.empty()) {
            cs = cellIntervalsToCoords(intervals, tag);
        } else {
            cs = cellIntervalsToCoords(intervals);
        }
        std::vector<coords_scaled_t> scaledCoords = coordsToScaledCoords(cs);
        res = scaledCoords;
    }

    // Function: cellIntervalsToCoords
    std::vector<coords_t> cellIntervalsToCoords(
        const std::vector<cell_interval_t>& ivls,
        const std::string& tag = "")
    {
        std::vector<coords_t> res(ivls.size());
        for (size_t i = 0; i < ivls.size(); ++i) {
            res[i].Or = ivls[i].getOrientation();
            
            int xi, xe;
            convertInterval(xi, xe, ivls[i], DIR_X);
            res[i].Xi = xi;
            res[i].Xe = xe;

            int yi, ye;
            convertInterval(yi, ye, ivls[i], DIR_Y);
            res[i].Yi = yi;
            res[i].Ye = ye;

            int zi, ze;
            convertInterval(zi, ze, ivls[i], DIR_Z);
            res[i].Zi = zi;
            res[i].Ze = ze;

            if (!tag.empty()) {
                res[i].tag = tag;
            } else {
                res[i].tag = "";
            }
        }
        return res;
    }

    // Helper subroutine for cellIntervalsToCoords
    void convertInterval(
        int& xi,
        int& xe,
        const cell_interval_t& interval,
        int dir)
    {
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

    // Function: getPixelFromElementId
    pixel_t getPixelFromElementId(
        const mesh_t& mesh,
        int id)
    {
        pixel_t res;
        bool nodeFound = false;
        node_t node = mesh.getNode(id, nodeFound);
        
        if (nodeFound) {
            res = mesh.nodeToPixel(node);
        } else {
            std::cerr << "Error converting pixel. Node not found." << std::endl;
            // In Fortran 'stop' terminates the program. In C++, we might throw or exit.
            // Exiting to match behavior closely.
            exit(1);
        }
        return res;
    }

    // Function: vectorToDiagonalMatrix
    std::vector<std::vector<rkind>> vectorToDiagonalMatrix(
        const std::vector<rkind>& vector)
    {
        int n = static_cast<int>(vector.size());
        std::vector<std::vector<rkind>> res(n, std::vector<rkind>(n, 0.0));
        for (int i = 0; i < n; ++i) {
            res[i][i] = vector[i];
        }
        return res;
    }

    // Function: scalarToMatrix
    std::vector<std::vector<rkind>> scalarToMatrix(
        rkind scalar)
    {
        std::vector<std::vector<rkind>> res(1, std::vector<rkind>(1, 0.0));
        res[0][0] = scalar;
        return res;
    }

    // Subroutine: splitLineIntoWords
    void splitLineIntoWords(
        const std::string& line,
        std::vector<std::string>& words)
    {
        int lenstr = static_cast<int>(line.length());
        int nwords = 0;
        int i = 0;

        // First pass: count words
        while (i < lenstr) {
            while (i < lenstr && (line[i] == ' ' || line[i] == '\t')) {
                i++;
            }
            if (i >= lenstr) break;
            nwords++;
            while (i < lenstr && (line[i] != ' ' && line[i] != '\t')) {
                i++;
            }
        }

        if (nwords == 0) {
            words.clear();
            return;
        }

        words.resize(nwords);
        i = 0;
        int start = 0;
        nwords = 0;
        while (i < lenstr) {
            while (i < lenstr && (line[i] == ' ' || line[i] == '\t')) {
                i++;
            }
            if (i >= lenstr) break;
            start = i;
            while (i < lenstr && (line[i] != ' ' && line[i] != '\t')) {
                i++;
            }
            int wlen = i - start;
            nwords++;
            words[nwords - 1] = line.substr(start, wlen);
        }
    }

    // Function: to_upper
    std::string to_upper(
        const std::string& str)
    {
        std::string res = str;
        for (size_t i = 0; i < res.length(); ++i) {
            if (res[i] >= 'a' && res[i] <= 'z') {
                res[i] = static_cast<char>(static_cast<int>(res[i]) - 32);
            }
        }
        return res;
    }

} // namespace parser_tools_m