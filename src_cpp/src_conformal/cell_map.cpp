#include <any>
#include <vector>

#include "conformal_types.h"
#include "fhash_m.h"
#include "cell_map_m.h"

namespace cell_map_m {

namespace detail {

bool isNewSide(const std::vector<side_t>& sides, const side_t& side) {
    for (const auto& s : sides) {
        if (s.isEquiv(side)) {
            return false;
        }
    }
    return true;
}

std::vector<int> cellKey3(const std::vector<int>& k) {
    std::vector<int> out(3, 0);
    for (std::size_t i = 0; i < k.size() && i < 3; ++i) {
        out[i] = k[i];
    }
    return out;
}

} // namespace detail

bool cell_hasKey(const cell_map_t& tbl, const std::vector<int>& k) {
    int stat = 0;
    tbl.check_key(key(detail::cellKey3(k)), stat);
    return stat == 0;
}

bool side_hasKey(const side_tris_map_t& tbl, const std::vector<int>& k) {
    int stat = 0;
    tbl.check_key(key(k), stat);
    return stat == 0;
}

bool cell_ratio_hasKey(const cell_ratios_map_t& tbl, const std::vector<int>& k) {
    int stat = 0;
    tbl.check_key(key(detail::cellKey3(k)), stat);
    return stat == 0;
}

bool cell_map_t::hasKey(const std::vector<int>& k) const {
    return cell_hasKey(*this, k);
}

bool side_tris_map_t::hasKey(const std::vector<int>& k) const {
    return side_hasKey(*this, k);
}

bool cell_ratios_map_t::hasKey(const std::vector<int>& k) const {
    return cell_ratio_hasKey(*this, k);
}

void buildSideToTrisMap(side_triangle_map_t& res, const std::vector<triangle_t>& triangles) {
    if (res.keys.empty()) {
        res.keys.clear();
    }
    for (const auto& tri : triangles) {
        const std::vector<side_t> sides = tri.getSides();
        for (const auto& side : sides) {
            if (side.isOnAnyEdge()) {
                res.addTriangleToSide(side, tri);
            }
        }
    }
}

void addTriangleToSideImpl(side_triangle_map_t& tbl, const side_t& side, const triangle_t& triangle) {
    std::vector<int> side_dir(4);
    const std::vector<int> cell = side.getCell();
    side_dir[0] = cell[0];
    side_dir[1] = cell[1];
    side_dir[2] = cell[2];
    side_dir[3] = side.getEdge();

    if (tbl.hasKey(side_dir)) {
        std::any alloc_list;
        tbl.get_raw(key(side_dir), alloc_list);
        auto elems = std::any_cast<element_set_t>(alloc_list);
        elems.triangles.push_back(triangle);
        tbl.set(key(side_dir), elems);
    } else {
        element_set_t elems;
        elems.triangles.push_back(triangle);
        tbl.set(key(side_dir), elems);

        side_dir_t entry;
        entry.side_dir = side_dir;
        tbl.keys.push_back(entry);
    }
}

void side_triangle_map_t::addTriangleToSide(const side_t& side, const triangle_t& triangle) {
    addTriangleToSideImpl(*this, side, triangle);
}

void buildSideMap(side_tris_map_t& res, const std::vector<triangle_t>& triangles) {
    side_triangle_map_t tri_map;
    buildSideToTrisMap(tri_map, triangles);

    const std::vector<side_dir_t> keys = tri_map.keys;
    for (const auto& entry : keys) {
        element_set_t elems;
        elems.triangles = tri_map.getTrianglesFromSide(entry.side_dir);
        res.set(key(entry.side_dir), elems);
    }
    res.keys = keys;
}

void buildCellMap(cell_map_t& res, const ConformalPECElements_t& volume) {
    triangle_map_t tri_map;
    interval_map_t interval_map;
    side_map_t side_map;
    side_map_t side_map_on;

    buildMapOfTrisOnFaces(tri_map, volume.triangles);
    buildMapOfIntervals(interval_map, volume.intervals);
    buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, volume.triangles);
    buildMapOfSidesOnEdgeFromTrisOnFaces(side_map_on, volume.triangles);

    std::vector<cell_t> keys = mergeKeys(tri_map.keys, side_map.keys);
    keys = mergeKeys(keys, side_map_on.keys);
    keys = mergeKeys(keys, interval_map.keys);

    for (const auto& entry : keys) {
        element_set_t elems;
        elems.triangles = tri_map.getTrianglesInCell(entry.cell);
        elems.sides = side_map.getSidesInCell(entry.cell);
        elems.sides_on = side_map_on.getOnSidesInCell(entry.cell);
        elems.intervals = interval_map.getIntervalsInCell(entry.cell);
        res.set(key(entry.cell), elems);
    }
    res.keys = keys;
}

std::vector<cell_t> mergeKeys(const std::vector<cell_t>& tri_keys, const std::vector<cell_t>& side_keys) {
    std::vector<cell_t> res;
    if (tri_keys.empty()) {
        res.clear();
    } else {
        res = tri_keys;
    }
    for (const auto& sk : side_keys) {
        if (isNewKey(tri_keys, sk)) {
            addKey(res, sk);
        }
    }
    return res;
}

void addKey(std::vector<cell_t>& keys, const cell_t& new_key) {
    keys.push_back(new_key);
}

bool isNewKey(const std::vector<cell_t>& keys, const cell_t& k) {
    for (const auto& entry : keys) {
        if (entry.cell == k.cell) {
            return false;
        }
    }
    return true;
}

void buildMapOfTrisOnFaces(triangle_map_t& res, const std::vector<triangle_t>& triangles) {
    if (res.keys.empty()) {
        res.keys.clear();
    }
    for (const auto& tri : triangles) {
        if (tri.isOnAnyFace()) {
            res.addTriangle(tri);
        }
    }
}

void buildMapOfIntervals(interval_map_t& res, const std::vector<interval_t>& intervals) {
    if (res.keys.empty()) {
        res.keys.clear();
    }
    for (const auto& interval : intervals) {
        res.addInterval(interval);
    }
}

void buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map_t& res,
                                                   const std::vector<triangle_t>& triangles) {
    if (res.keys.empty()) {
        res.keys.clear();
    }
    for (const auto& tri : triangles) {
        if (!tri.isOnAnyFace()) {
            const std::vector<side_t> sides = tri.getSides();
            for (const auto& side : sides) {
                if (side.isOnAnyFace() || side.isOnAnyEdge()) {
                    res.addSide(side);
                }
            }
        }
    }
}

void buildMapOfSidesOnEdgeFromTrisOnFaces(side_map_t& res, const std::vector<triangle_t>& triangles) {
    if (res.keys.empty()) {
        res.keys.clear();
    }
    for (const auto& tri : triangles) {
        if (tri.isOnAnyFace()) {
            const std::vector<side_t> sides = tri.getSides();
            for (const auto& side : sides) {
                if (side.isOnAnyEdge()) {
                    res.addSideOn(side);
                }
            }
        }
    }
}

void addTriangleImpl(triangle_map_t& tbl, const triangle_t& triangle) {
    const std::vector<int> cell = triangle.getCell();
    if (tbl.hasKey(cell)) {
        std::any alloc_list;
        tbl.get_raw(key(cell), alloc_list);
        auto elems = std::any_cast<element_set_t>(alloc_list);
        elems.triangles.push_back(triangle);
        tbl.set(key(cell), elems);
    } else {
        element_set_t elems;
        elems.triangles.push_back(triangle);
        tbl.set(key(cell), elems);

        cell_t entry;
        entry.cell = cell;
        tbl.keys.push_back(entry);
    }
}

void triangle_map_t::addTriangle(const triangle_t& triangle) {
    addTriangleImpl(*this, triangle);
}

void addIntervalImpl(interval_map_t& tbl, const interval_t& interval) {
    side_t aux;
    aux.init.position = {
        static_cast<double>(interval.ini.cell[0]),
        static_cast<double>(interval.ini.cell[1]),
        static_cast<double>(interval.ini.cell[2]),
    };
    aux.end.position = {
        static_cast<double>(interval.end.cell[0]),
        static_cast<double>(interval.end.cell[1]),
        static_cast<double>(interval.end.cell[2]),
    };
    const std::vector<int> cell = aux.getCell();

    if (tbl.hasKey(cell)) {
        std::any alloc_list;
        tbl.get_raw(key(cell), alloc_list);
        auto elems = std::any_cast<element_set_t>(alloc_list);
        elems.intervals.push_back(interval);
        tbl.set(key(cell), elems);
    } else {
        element_set_t elems;
        elems.intervals.push_back(interval);
        tbl.set(key(cell), elems);

        cell_t entry;
        entry.cell = aux.getCell();
        tbl.keys.push_back(entry);
    }
}

void interval_map_t::addInterval(const interval_t& interval) {
    addIntervalImpl(*this, interval);
}

void addSideImpl(side_map_t& tbl, const side_t& side) {
    const std::vector<int> cell = side.getCell();
    if (tbl.hasKey(cell)) {
        std::any alloc_list;
        tbl.get_raw(key(cell), alloc_list);
        auto elems = std::any_cast<element_set_t>(alloc_list);
        if (detail::isNewSide(elems.sides, side)) {
            elems.sides.push_back(side);
            tbl.set(key(cell), elems);
        }
    } else {
        element_set_t elems;
        elems.sides.push_back(side);
        tbl.set(key(cell), elems);

        cell_t entry;
        entry.cell = cell;
        tbl.keys.push_back(entry);
    }
}

void side_map_t::addSide(const side_t& side) {
    addSideImpl(*this, side);
}

void addSideOnImpl(side_map_t& tbl, const side_t& side) {
    const std::vector<int> cell = side.getCell();
    if (tbl.hasKey(cell)) {
        std::any alloc_list;
        tbl.get_raw(key(cell), alloc_list);
        auto elems = std::any_cast<element_set_t>(alloc_list);
        elems.sides_on.push_back(side);
        tbl.set(key(cell), elems);
    } else {
        element_set_t elems;
        elems.sides_on.push_back(side);
        tbl.set(key(cell), elems);

        cell_t entry;
        entry.cell = cell;
        tbl.keys.push_back(entry);
    }
}

void side_map_t::addSideOn(const side_t& side) {
    addSideOnImpl(*this, side);
}

std::vector<triangle_t> getTrianglesInCellImpl(const cell_map_t& tbl, const std::vector<int>& k) {
    if (tbl.hasKey(k)) {
        std::any alloc_list;
        tbl.get_raw(key(detail::cellKey3(k)), alloc_list);
        return std::any_cast<element_set_t>(alloc_list).triangles;
    }
    return {};
}

std::vector<interval_t> getIntervalsInCellImpl(const cell_map_t& tbl, const std::vector<int>& k) {
    if (tbl.hasKey(k)) {
        std::any alloc_list;
        tbl.get_raw(key(detail::cellKey3(k)), alloc_list);
        return std::any_cast<element_set_t>(alloc_list).intervals;
    }
    return {};
}

std::vector<triangle_t> getTrianglesFromSideImpl(const side_tris_map_t& tbl, const std::vector<int>& k) {
    if (tbl.hasKey(k)) {
        std::any alloc_list;
        tbl.get_raw(key(k), alloc_list);
        return std::any_cast<element_set_t>(alloc_list).triangles;
    }
    return {};
}

std::vector<side_t> getSidesInCellImpl(const cell_map_t& tbl, const std::vector<int>& k) {
    if (tbl.hasKey(k)) {
        std::any alloc_list;
        tbl.get_raw(key(detail::cellKey3(k)), alloc_list);
        return std::any_cast<element_set_t>(alloc_list).sides;
    }
    return {};
}

std::vector<side_t> getOnSidesInCellImpl(const cell_map_t& tbl, const std::vector<int>& k) {
    if (tbl.hasKey(k)) {
        std::any alloc_list;
        tbl.get_raw(key(detail::cellKey3(k)), alloc_list);
        return std::any_cast<element_set_t>(alloc_list).sides_on;
    }
    return {};
}

std::vector<triangle_t> cell_map_t::getTrianglesInCell(const std::vector<int>& k) const {
    return getTrianglesInCellImpl(*this, k);
}

std::vector<side_t> cell_map_t::getSidesInCell(const std::vector<int>& k) const {
    return getSidesInCellImpl(*this, k);
}

std::vector<side_t> cell_map_t::getOnSidesInCell(const std::vector<int>& k) const {
    return getOnSidesInCellImpl(*this, k);
}

std::vector<interval_t> cell_map_t::getIntervalsInCell(const std::vector<int>& k) const {
    return getIntervalsInCellImpl(*this, k);
}

std::vector<triangle_t> side_tris_map_t::getTrianglesFromSide(const std::vector<int>& k) const {
    return getTrianglesFromSideImpl(*this, k);
}

void addFaceRatioImpl(cell_ratios_map_t& tbl, const cell_t& cell, int direction, rkind ratio) {
    if (direction < conformal_types_m::FACE_X || direction > conformal_types_m::FACE_Z) {
        return;
    }
    if (tbl.hasKey(cell.cell)) {
        std::any alloc_list;
        tbl.get_raw(key(cell.cell), alloc_list);
        auto ratios = std::any_cast<cell_ratios_t>(alloc_list);
        ratios.area[static_cast<std::size_t>(direction - 1)] = ratio;
        tbl.set(key(cell.cell), ratios);
    } else {
        cell_ratios_t ratios;
        ratios.area[static_cast<std::size_t>(direction - 1)] = ratio;
        tbl.set(key(cell.cell), ratios);
        tbl.keys.push_back(cell);
    }
}

void cell_ratios_map_t::addFaceRatio(const cell_t& cell, int direction, rkind ratio) {
    addFaceRatioImpl(*this, cell, direction, ratio);
}

void addEdgeRatioImpl(cell_ratios_map_t& tbl, const cell_t& cell, int direction, rkind ratio) {
    if (direction < conformal_types_m::EDGE_X || direction > conformal_types_m::EDGE_Z) {
        return;
    }
    if (tbl.hasKey(cell.cell)) {
        std::any alloc_list;
        tbl.get_raw(key(cell.cell), alloc_list);
        auto ratios = std::any_cast<cell_ratios_t>(alloc_list);
        ratios.length[static_cast<std::size_t>(direction - 1)] = ratio;
        tbl.set(key(cell.cell), ratios);
    } else {
        cell_ratios_t ratios;
        ratios.length[static_cast<std::size_t>(direction - 1)] = ratio;
        tbl.set(key(cell.cell), ratios);
        tbl.keys.push_back(cell);
    }
}

void cell_ratios_map_t::addEdgeRatio(const cell_t& cell, int direction, rkind ratio) {
    addEdgeRatioImpl(*this, cell, direction, ratio);
}

cell_ratios_t getCellRatiosInCellImpl(const cell_ratios_map_t& tbl, const std::vector<int>& k) {
    cell_ratios_t res;
    if (tbl.hasKey(k)) {
        std::any alloc_list;
        tbl.get_raw(key(detail::cellKey3(k)), alloc_list);
        res = std::any_cast<cell_ratios_t>(alloc_list);
    }
    return res;
}

cell_ratios_t cell_ratios_map_t::getCellRatiosInCell(const std::vector<int>& k) const {
    return getCellRatiosInCellImpl(*this, k);
}

} // namespace cell_map_m
