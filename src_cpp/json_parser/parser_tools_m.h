#ifndef PARSER_TOOLS_M_H
#define PARSER_TOOLS_M_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <cstdlib>

#include "mesh_m.h"
#include "cells_m.h"
#include "../main/nfde_types.h"

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


    inline std::vector<NFDETypes_m::coords_t> cellIntervalsToCoords(
        const std::vector<cells_m::cell_interval_t>& ivls, const std::string& tag = "") {
        std::vector<NFDETypes_m::coords_t> res;
        res.resize(ivls.size());
        for (size_t i = 0; i < ivls.size(); ++i) {
            res[i].Or = ivls[i].getOrientation();
            int dir;
            for (dir = 1; dir <= 3; ++dir) {
                if (ivls[i].ini.cell[dir-1] == ivls[i].end.cell[dir-1]) {
                    switch(dir) {
                        case 1: res[i].Xi = ivls[i].ini.cell[0]; res[i].Xe = ivls[i].end.cell[0]; break;
                        case 2: res[i].Yi = ivls[i].ini.cell[1]; res[i].Ye = ivls[i].end.cell[1]; break;
                        case 3: res[i].Zi = ivls[i].ini.cell[2]; res[i].Ze = ivls[i].end.cell[2]; break;
                    }
                } else {
                    switch(dir) {
                        case 1: res[i].Xi = std::min(ivls[i].ini.cell[0], ivls[i].end.cell[0]);
                                res[i].Xe = std::max(ivls[i].ini.cell[0], ivls[i].end.cell[0]) - 1; break;
                        case 2: res[i].Yi = std::min(ivls[i].ini.cell[1], ivls[i].end.cell[1]);
                                res[i].Ye = std::max(ivls[i].ini.cell[1], ivls[i].end.cell[1]) - 1; break;
                        case 3: res[i].Zi = std::min(ivls[i].ini.cell[2], ivls[i].end.cell[2]);
                                res[i].Ze = std::max(ivls[i].ini.cell[2], ivls[i].end.cell[2]) - 1; break;
                    }
                }
            }
            if (!tag.empty()) res[i].tag = tag;
        }
        return res;
    }

    inline std::vector<NFDETypes_m::coords_t> cellRegionToCoords(
        const cells_m::cell_region_t& cellRegion, int cellType,
        const std::string& tag = "") {
        auto intervals = getIntervalsInCellRegions({cellRegion}, cellType);
        return cellIntervalsToCoords(intervals, tag);
    }

    inline std::vector<NFDETypes_m::coords_scaled_t> coordsToScaledCoords(
        const std::vector<NFDETypes_m::coords_t>& cs) {
        std::vector<NFDETypes_m::coords_scaled_t> res;
        res.resize(cs.size());
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
            int or_abs = std::abs(cs[i].Or);
            if (cs[i].Or == iEx) res[i].xc = 1.0;
            else if (cs[i].Or == -iEx) res[i].xc = -1.0;
            else if (cs[i].Or == iEy) res[i].yc = 1.0;
            else if (cs[i].Or == -iEy) res[i].yc = -1.0;
            else if (cs[i].Or == iEz) res[i].zc = 1.0;
            else if (cs[i].Or == -iEz) res[i].zc = -1.0;
        }
        return res;
    }

} // namespace parser_tools_m

#endif // PARSER_TOOLS_M_H
