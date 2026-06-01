#ifndef CELL_MAP_M_H
#define CELL_MAP_M_H

#include <string>
#include <vector>

#include "conformal_types.h"
#include "fhash_m.h"

namespace cell_map_m {

using conformal_types_m::interval_t;
using conformal_types_m::side_t;
using conformal_types_m::triangle_t;

using rkind = double;
using fhash_tbl_t = fhash_m::fhash_tbl_t;
using fhash_key_t = fhash_m::fhash_key_t;

inline fhash_key_t key(const std::vector<int>& k) { return fhash_m::key(k); }
inline fhash_key_t key(int k) { return fhash_m::key(k); }
inline fhash_key_t key(const std::string& k) { return fhash_m::key(k); }

struct ConformalPECElements_t {
    std::vector<triangle_t> triangles;
    std::vector<interval_t> intervals;
    std::string tag;
};

struct element_set_t {
    std::vector<triangle_t> triangles;
    std::vector<side_t> sides;
    std::vector<side_t> sides_on;
    std::vector<interval_t> intervals;
};

struct cell_t {
    std::vector<int> cell;
};

struct cell_ratios_t {
    std::vector<rkind> area = {1, 1, 1};
    std::vector<rkind> length = {1, 1, 1};
};

class cell_ratios_map_t : public fhash_tbl_t {
public:
    std::vector<cell_t> keys;

    bool hasKey(const std::vector<int>& k) const;
    void addFaceRatio(const cell_t& cell, int direction, rkind ratio);
    void addEdgeRatio(const cell_t& cell, int direction, rkind ratio);
    cell_ratios_t getCellRatiosInCell(const std::vector<int>& k) const;
};

class cell_map_t : public fhash_tbl_t {
public:
    std::vector<cell_t> keys;

    bool hasKey(const std::vector<int>& k) const;
    std::vector<triangle_t> getTrianglesInCell(const std::vector<int>& k) const;
    std::vector<side_t> getSidesInCell(const std::vector<int>& k) const;
    std::vector<side_t> getOnSidesInCell(const std::vector<int>& k) const;
    std::vector<interval_t> getIntervalsInCell(const std::vector<int>& k) const;
};

class triangle_map_t : public cell_map_t {
public:
    void addTriangle(const triangle_t& triangle);
};

class interval_map_t : public cell_map_t {
public:
    void addInterval(const interval_t& interval);
};

class side_map_t : public cell_map_t {
public:
    void addSide(const side_t& side);
    void addSideOn(const side_t& side);
};

struct side_dir_t {
    std::vector<int> side_dir;
};

class side_tris_map_t : public fhash_tbl_t {
public:
    std::vector<side_dir_t> keys;

    bool hasKey(const std::vector<int>& k) const;
    std::vector<triangle_t> getTrianglesFromSide(const std::vector<int>& k) const;
};

class side_triangle_map_t : public side_tris_map_t {
public:
    void addTriangleToSide(const side_t& side, const triangle_t& triangle);
};

void buildSideToTrisMap(side_triangle_map_t& res, const std::vector<triangle_t>& triangles);
void buildSideMap(side_tris_map_t& res, const std::vector<triangle_t>& triangles);
void buildCellMap(cell_map_t& res, const ConformalPECElements_t& volume);

std::vector<cell_t> mergeKeys(const std::vector<cell_t>& tri_keys, const std::vector<cell_t>& side_keys);
void addKey(std::vector<cell_t>& keys, const cell_t& new_key);
bool isNewKey(const std::vector<cell_t>& keys, const cell_t& k);

void buildMapOfTrisOnFaces(triangle_map_t& res, const std::vector<triangle_t>& triangles);
void buildMapOfIntervals(interval_map_t& res, const std::vector<interval_t>& intervals);
void buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map_t& res,
                                                   const std::vector<triangle_t>& triangles);
void buildMapOfSidesOnEdgeFromTrisOnFaces(side_map_t& res,
                                          const std::vector<triangle_t>& triangles);

} // namespace cell_map_m

#endif // CELL_MAP_M_H
