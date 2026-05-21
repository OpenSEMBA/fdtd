#ifndef PARSER_TOOLS_M_H
#define PARSER_TOOLS_M_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <cstdlib>

// Forward declarations
namespace json_module {
    class json_core;
    class json_value;
}

#include "mesh_m.h"
#include "cells_m.h"

namespace parser_tools_m {

    constexpr int PARSER_TOOLS_BUFSIZE = 256;
    using rkind = double;

    constexpr int iEx = 1;
    constexpr int iEy = 2;
    constexpr int iEz = 3;

    // Function declarations
    std::vector<cells_m::cell_interval_t> getIntervalsInCellRegions(
        const std::vector<cells_m::cell_region_t>& cellRegions,
        int cellType = -1);

    mesh_m::pixel_t getPixelFromElementId(
        const mesh_m::mesh_t& mesh, int id);

    std::vector<std::vector<rkind>> vectorToDiagonalMatrix(
        const std::vector<rkind>& vector);

    std::vector<std::vector<rkind>> scalarToMatrix(rkind scalar);

    void splitLineIntoWords(const std::string& line, std::vector<std::string>& words);

    std::string to_upper(const std::string& str);

    // Inline implementations

    inline std::vector<cells_m::cell_interval_t> getIntervalsInCellRegions(
        const std::vector<cells_m::cell_region_t>& cellRegions, int cellType) {
        int numberOfIntervals = 0;
        for (const auto& region : cellRegions) {
            if (cellType != -1) {
                for (const auto& iv : region.intervals) {
                    if (iv.getType() == cellType) numberOfIntervals++;
                }
            } else {
                numberOfIntervals += static_cast<int>(region.intervals.size());
            }
        }
        std::vector<cells_m::cell_interval_t> result(numberOfIntervals);
        int idx = 0;
        for (const auto& region : cellRegions) {
            auto ivs = (cellType != -1) ? region.getIntervalsOfType(cellType) : region.intervals;
            for (const auto& iv : ivs) { if (idx < numberOfIntervals) result[idx++] = iv; }
        }
        return result;
    }

    inline mesh_m::pixel_t getPixelFromElementId(const mesh_m::mesh_t& mesh, int id) {
        mesh_m::pixel_t res{};
        bool found = false;
        auto node = mesh.getNode(id, found);
        if (found) res = mesh.nodeToPixel(node);
        else {
            std::cerr << "Error converting pixel. Node not found." << std::endl;
            std::exit(1);
        }
        return res;
    }

    inline std::vector<std::vector<rkind>> vectorToDiagonalMatrix(const std::vector<rkind>& v) {
        int n = static_cast<int>(v.size());
        std::vector<std::vector<rkind>> res(n, std::vector<rkind>(n, 0.0));
        for (int i = 0; i < n; ++i) res[i][i] = v[i];
        return res;
    }

    inline std::vector<std::vector<rkind>> scalarToMatrix(rkind scalar) {
        return std::vector<std::vector<rkind>>(1, std::vector<rkind>(1, scalar));
    }

    inline void splitLineIntoWords(const std::string& line, std::vector<std::string>& words) {
        int lenstr = static_cast<int>(line.length());
        int nwords = 0, i = 0;
        while (i < lenstr) {
            while (i < lenstr && (line[i] == ' ' || line[i] == '\t')) i++;
            if (i >= lenstr) break;
            nwords++;
            while (i < lenstr && line[i] != ' ' && line[i] != '\t') i++;
        }
        if (nwords == 0) { words.clear(); return; }
        words.resize(nwords);
        i = 0; int start = 0, n = 0;
        while (i < lenstr) {
            while (i < lenstr && (line[i] == ' ' || line[i] == '\t')) i++;
            if (i >= lenstr) break;
            start = i;
            while (i < lenstr && line[i] != ' ' && line[i] != '\t') i++;
            int wlen = i - start;
            n++;
            words[n - 1] = line.substr(start, wlen);
        }
    }

    inline std::string to_upper(const std::string& str) {
        std::string res = str;
        for (size_t i = 0; i < res.length(); ++i) {
            if (res[i] >= 'a' && res[i] <= 'z') {
                res[i] = static_cast<char>(static_cast<int>(res[i]) - 32);
            }
        }
        return res;
    }

} // namespace parser_tools_m

#endif // PARSER_TOOLS_M_H
