#include <vector>
#include <memory>
#include <algorithm>
#include <iostream>

// Forward declarations and includes for external modules/types
// Assuming these are defined in other headers based on the use statements

namespace geometry_m {
    // Placeholder for geometry types
    struct triangle_t {
        std::vector<double> getSides() const; // Returns vector of side_t
        bool isOnAnyFace() const;
        bool isOnAnyEdge() const;
        std::vector<int> getCell() const;
    };

    struct side_t {
        bool isOnAnyEdge() const;
        bool isOnAnyFace() const;
        std::vector<int> getCell() const;
        int getEdge() const;
        struct {
            std::vector<int> position;
        } ini, end;
        void init; // Placeholder for method
    };

    struct interval_t {
        struct {
            std::vector<int> cell;
        } ini, end;
    };

    enum Face { FACE_X, FACE_Y, FACE_Z };

    bool isNewSide(const std::vector<side_t>& sides, const side_t& side) {
        for (const auto& s : sides) {
            // Assuming side_t has equality operator or comparison logic
            if (s == side) return false;
        }
        return true;
    }
}

namespace fhash {
    // Placeholder for hash table
    class fhash_tbl_t {
    public:
        std::vector<std::vector<int>> keys; // Simplified key storage
        bool hasKey(const std::vector<int>& k) const;
        void check_key(const std::vector<int>& k, int& stat) const;
        void get_raw(const std::vector<int>& k, std::shared_ptr<void>& alloc_list) const;
        void set(const std::vector<int>& k, std::shared_ptr<void> value);
    };

    std::vector<int> key(const std::vector<int>& k) {
        return k;
    }
}

namespace NFDETypes_m {
    using rkind = double;

    struct ConformalPECElements_t {
        std::vector<geometry_m::triangle_t> triangles;
        std::vector<geometry_m::interval_t> intervals;
    };
}

namespace cell_map_m {

    using triangle_t = geometry_m::triangle_t;
    using side_t = geometry_m::side_t;
    using interval_t = geometry_m::interval_t;
    using rkind = NFDETypes_m::rkind;
    using ConformalPECElements_t = NFDETypes_m::ConformalPECElements_t;
    using fhash_tbl_t = fhash::fhash_tbl_t;
    using key = fhash::key;

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

        bool hasKey(const std::vector<int>& k) const {
            return cell_ratio_hasKey(*this, k);
        }

        void addFaceRatio(const cell_t& cell, int direction, rkind ratio) {
            ::addFaceRatio(*this, cell, direction, ratio);
        }

        void addEdgeRatio(const cell_t& cell, int direction, rkind ratio) {
            ::addEdgeRatio(*this, cell, direction, ratio);
        }

        cell_ratios_t getCellRatiosInCell(const std::vector<int>& k) const {
            return getCellRatiosInCell(*this, k);
        }
    };

    class cell_map_t : public fhash_tbl_t {
    public:
        std::vector<cell_t> keys;

        bool hasKey(const std::vector<int>& k) const {
            return cell_hasKey(*this, k);
        }

        std::vector<triangle_t> getTrianglesInCell(const std::vector<int>& k) const {
            return getTrianglesInCell(*this, k);
        }

        std::vector<side_t> getSidesInCell(const std::vector<int>& k) const {
            return getSidesInCell(*this, k);
        }

        std::vector<side_t> getOnSidesInCell(const std::vector<int>& k) const {
            return getOnSidesInCell(*this, k);
        }

        std::vector<interval_t> getIntervalsInCell(const std::vector<int>& k) const {
            return getIntervalsInCell(*this, k);
        }
    };

    class triangle_map_t : public cell_map_t {
    public:
        void addTriangle(const triangle_t& triangle) {
            ::addTriangle(*this, triangle);
        }
    };

    class interval_map_t : public cell_map_t {
    public:
        void addInterval(const interval_t& interval) {
            ::addInterval(*this, interval);
        }
    };

    class side_map_t : public cell_map_t {
    public:
        void addSide(const side_t& side) {
            ::addSide(*this, side);
        }

        void addSideOn(const side_t& side) {
            ::addSideOn(*this, side);
        }
    };

    struct side_dir_t {
        std::vector<int> side_dir;
    };

    class side_tris_map_t : public fhash_tbl_t {
    public:
        std::vector<side_dir_t> keys;

        bool hasKey(const std::vector<int>& k) const {
            return side_hasKey(*this, k);
        }

        std::vector<triangle_t> getTrianglesFromSide(const std::vector<int>& k) const {
            return getTrianglesFromSide(*this, k);
        }
    };

    class side_triangle_map_t : public side_tris_map_t {
    public:
        void addTriangleToSide(const side_t& side, const triangle_t& triangle) {
            ::addTriangleToSide(*this, side, triangle);
        }
    };

    // Helper functions and subroutines

    void buildSideToTrisMap(side_triangle_map_t& res, const std::vector<triangle_t>& triangles) {
        if (res.keys.empty()) {
            res.keys = {};
        }
        for (size_t i = 0; i < triangles.size(); ++i) {
            std::vector<side_t> sides = triangles[i].getSides();
            for (int j = 0; j < 3; ++j) {
                if (sides[j].isOnAnyEdge()) {
                    res.addTriangleToSide(sides[j], triangles[i]);
                }
            }
        }
    }

    void addTriangleToSide(side_triangle_map_t& this_obj, const side_t& side, const triangle_t& triangle) {
        std::vector<int> side_dir(4);
        side_dir[0] = side.getCell()[0];
        side_dir[1] = side.getCell()[1];
        side_dir[2] = side.getCell()[2];
        side_dir[3] = side.getEdge();

        if (this_obj.hasKey(side_dir)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(side_dir), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                element_set_t aux_list_copy;
                aux_list_copy.triangles.resize(aux_list->triangles.size() + 1);
                for (size_t i = 0; i < aux_list->triangles.size(); ++i) {
                    aux_list_copy.triangles[i] = aux_list->triangles[i];
                }
                aux_list_copy.triangles[aux_list->triangles.size()] = triangle;
                
                // Update the stored object
                aux_list->triangles = aux_list_copy.triangles;
                this_obj.set(key(side_dir), std::shared_ptr<void>(aux_list));
            }
        } else {
            element_set_t aux_list;
            aux_list.triangles.resize(1);
            aux_list.triangles[0] = triangle;
            this_obj.set(key(side_dir), std::make_shared<element_set_t>(aux_list));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({side_dir});
            this_obj.keys = aux_keys;
        }
    }

    void buildSideMap(side_tris_map_t& res, const std::vector<triangle_t>& triangles) {
        side_triangle_map_t tri_map;
        buildSideToTrisMap(tri_map, triangles);
        
        std::vector<side_dir_t> keys = tri_map.keys;
        for (size_t i = 0; i < keys.size(); ++i) {
            element_set_t elems;
            elems.triangles = tri_map.getTrianglesFromSide(keys[i].side_dir);
            res.set(key(keys[i].side_dir), std::make_shared<element_set_t>(elems));
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

        for (size_t i = 0; i < keys.size(); ++i) {
            element_set_t elems;
            elems.triangles = tri_map.getTrianglesInCell(keys[i].cell);
            elems.sides = side_map.getSidesInCell(keys[i].cell);
            elems.sides_on = side_map_on.getOnSidesInCell(keys[i].cell);
            elems.intervals = interval_map.getIntervalsInCell(keys[i].cell);
            res.set(key(keys[i].cell), std::make_shared<element_set_t>(elems));
        }
        res.keys = keys;
    }

    std::vector<cell_t> mergeKeys(const std::vector<cell_t>& tri_keys, const std::vector<cell_t>& side_keys) {
        std::vector<cell_t> res;
        if (tri_keys.empty()) {
            res = {};
        } else {
            res = tri_keys;
        }
        
        for (const auto& sk : side_keys) {
            if (isNewKey(res, sk)) {
                addKey(res, sk);
            }
        }
        return res;
    }

    void addKey(std::vector<cell_t>& keys, const cell_t& new_key) {
        std::vector<cell_t> aux = keys;
        aux.push_back(new_key);
        keys = aux;
    }

    bool isNewKey(const std::vector<cell_t>& keys, const cell_t& k) {
        for (const auto& key : keys) {
            if (key.cell == k.cell) {
                return false;
            }
        }
        return true;
    }

    void buildMapOfTrisOnFaces(triangle_map_t& res, const std::vector<triangle_t>& triangles) {
        if (res.keys.empty()) {
            res.keys = {};
        }
        for (const auto& tri : triangles) {
            if (tri.isOnAnyFace()) {
                res.addTriangle(tri);
            }
        }
    }

    void buildMapOfIntervals(interval_map_t& res, const std::vector<interval_t>& intervals) {
        if (res.keys.empty()) {
            res.keys = {};
        }
        for (const auto& interval : intervals) {
            res.addInterval(interval);
        }
    }

    void buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map_t& res, const std::vector<triangle_t>& triangles) {
        if (res.keys.empty()) {
            res.keys = {};
        }
        for (const auto& tri : triangles) {
            if (!tri.isOnAnyFace()) {
                std::vector<side_t> sides = tri.getSides();
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
            res.keys = {};
        }
        for (const auto& tri : triangles) {
            if (tri.isOnAnyFace()) {
                std::vector<side_t> sides = tri.getSides();
                for (const auto& side : sides) {
                    if (side.isOnAnyEdge()) {
                        res.addSideOn(side);
                    }
                }
            }
        }
    }

    bool cell_hasKey(const cell_map_t& this_obj, const std::vector<int>& k) {
        int stat = 0;
        this_obj.check_key(key(k), stat);
        return stat == 0;
    }

    bool side_hasKey(const side_tris_map_t& this_obj, const std::vector<int>& k) {
        int stat = 0;
        this_obj.check_key(key(k), stat);
        return stat == 0;
    }

    bool cell_ratio_hasKey(const cell_ratios_map_t& this_obj, const std::vector<int>& k) {
        int stat = 0;
        this_obj.check_key(key(k), stat);
        return stat == 0;
    }

    void addTriangle(triangle_map_t& this_obj, const triangle_t& triangle) {
        std::vector<int> cell = triangle.getCell();
        if (this_obj.hasKey(cell)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                element_set_t aux_list_copy;
                aux_list_copy.triangles.resize(aux_list->triangles.size() + 1);
                for (size_t i = 0; i < aux_list->triangles.size(); ++i) {
                    aux_list_copy.triangles[i] = aux_list->triangles[i];
                }
                aux_list_copy.triangles[aux_list->triangles.size()] = triangle;
                
                aux_list->triangles = aux_list_copy.triangles;
                this_obj.set(key(cell), std::shared_ptr<void>(aux_list));
            }
        } else {
            element_set_t aux_list;
            aux_list.triangles.resize(1);
            aux_list.triangles[0] = triangle;
            this_obj.set(key(cell), std::make_shared<element_set_t>(aux_list));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell});
            this_obj.keys = aux_keys;
        }
    }

    void addInterval(interval_map_t& this_obj, const interval_t& interval) {
        side_t aux;
        aux.ini.position = interval.ini.cell;
        aux.end.position = interval.end.cell;
        std::vector<int> cell = aux.getCell();

        if (this_obj.hasKey(cell)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                element_set_t aux_list_copy;
                aux_list_copy.intervals.resize(aux_list->intervals.size() + 1);
                for (size_t i = 0; i < aux_list->intervals.size(); ++i) {
                    aux_list_copy.intervals[i] = aux_list->intervals[i];
                }
                aux_list_copy.intervals[aux_list->intervals.size()] = interval;
                
                aux_list->intervals = aux_list_copy.intervals;
                this_obj.set(key(cell), std::shared_ptr<void>(aux_list));
            }
        } else {
            element_set_t aux_list;
            aux_list.intervals.resize(1);
            aux_list.intervals[0] = interval;
            this_obj.set(key(cell), std::make_shared<element_set_t>(aux_list));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({aux.getCell()});
            this_obj.keys = aux_keys;
        }
    }

    void addSide(side_map_t& this_obj, const side_t& side) {
        std::vector<int> cell = side.getCell();
        if (this_obj.hasKey(cell)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                if (isNewSide(aux_list->sides, side)) {
                    element_set_t aux_list_copy;
                    aux_list_copy.sides.resize(aux_list->sides.size() + 1);
                    for (size_t i = 0; i < aux_list->sides.size(); ++i) {
                        aux_list_copy.sides[i] = aux_list->sides[i];
                    }
                    aux_list_copy.sides[aux_list->sides.size()] = side;
                    
                    aux_list->sides = aux_list_copy.sides;
                    this_obj.set(key(cell), std::shared_ptr<void>(aux_list));
                }
            }
        } else {
            element_set_t aux_list;
            aux_list.sides.resize(1);
            aux_list.sides[0] = side;
            this_obj.set(key(cell), std::make_shared<element_set_t>(aux_list));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell});
            this_obj.keys = aux_keys;
        }
    }

    void addSideOn(side_map_t& this_obj, const side_t& side) {
        std::vector<int> cell = side.getCell();
        if (this_obj.hasKey(cell)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                element_set_t aux_list_copy;
                aux_list_copy.sides_on.resize(aux_list->sides_on.size() + 1);
                for (size_t i = 0; i < aux_list->sides_on.size(); ++i) {
                    aux_list_copy.sides_on[i] = aux_list->sides_on[i];
                }
                aux_list_copy.sides_on[aux_list->sides_on.size()] = side;
                
                aux_list->sides_on = aux_list_copy.sides_on;
                this_obj.set(key(cell), std::shared_ptr<void>(aux_list));
            }
        } else {
            element_set_t aux_list;
            aux_list.sides_on.resize(1);
            aux_list.sides_on[0] = side;
            this_obj.set(key(cell), std::make_shared<element_set_t>(aux_list));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell});
            this_obj.keys = aux_keys;
        }
    }

    std::vector<triangle_t> getTrianglesInCell(const cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<triangle_t> res;
        if (this_obj.hasKey(k)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(k), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                res = aux_list->triangles;
            }
        } else {
            res = {};
        }
        return res;
    }

    std::vector<interval_t> getIntervalsInCell(const cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<interval_t> res;
        if (this_obj.hasKey(k)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(k), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                res = aux_list->intervals;
            }
        } else {
            res = {};
        }
        return res;
    }

    std::vector<triangle_t> getTrianglesFromSide(const side_tris_map_t& this_obj, const std::vector<int>& k) {
        std::vector<triangle_t> res;
        if (this_obj.hasKey(k)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(k), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                res = aux_list->triangles;
            }
        } else {
            res = {};
        }
        return res;
    }

    std::vector<side_t> getSidesInCell(const cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<side_t> res;
        if (this_obj.hasKey(k)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(k), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                res = aux_list->sides;
            }
        } else {
            res = {};
        }
        return res;
    }

    std::vector<side_t> getOnSidesInCell(const cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<side_t> res;
        if (this_obj.hasKey(k)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(k), alloc_list);
            if (auto* aux_list = dynamic_cast<element_set_t*>(alloc_list.get())) {
                res = aux_list->sides_on;
            }
        } else {
            res = {};
        }
        return res;
    }

    void addFaceRatio(cell_ratios_map_t& this_obj, const cell_t& cell, int direction, rkind ratio) {
        if (this_obj.hasKey(cell.cell)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(cell.cell), alloc_list);
            if (auto* aux_list = dynamic_cast<cell_ratios_t*>(alloc_list.get())) {
                aux_list->area[direction] = ratio;
                this_obj.set(key(cell.cell), std::shared_ptr<void>(aux_list));
            }
        } else {
            cell_ratios_t aux_cell_ratio;
            aux_cell_ratio.area[direction] = ratio;
            this_obj.set(key(cell.cell), std::make_shared<cell_ratios_t>(aux_cell_ratio));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back(cell);
            this_obj.keys = aux_keys;
        }
    }

    void addEdgeRatio(cell_ratios_map_t& this_obj, const cell_t& cell, int direction, rkind ratio) {
        if (this_obj.hasKey(cell.cell)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(cell.cell), alloc_list);
            if (auto* aux_list = dynamic_cast<cell_ratios_t*>(alloc_list.get())) {
                aux_list->length[direction] = ratio;
                this_obj.set(key(cell.cell), std::shared_ptr<void>(aux_list));
            }
        } else {
            cell_ratios_t aux_cell_ratio;
            aux_cell_ratio.length[direction] = ratio;
            this_obj.set(key(cell.cell), std::make_shared<cell_ratios_t>(aux_cell_ratio));

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back(cell);
            this_obj.keys = aux_keys;
        }
    }

    cell_ratios_t getCellRatiosInCell(const cell_ratios_map_t& this_obj, const std::vector<int>& k) {
        cell_ratios_t res;
        if (this_obj.hasKey(k)) {
            std::shared_ptr<void> alloc_list;
            this_obj.get_raw(key(k), alloc_list);
            if (auto* aux_list = dynamic_cast<cell_ratios_t*>(alloc_list.get())) {
                res = *aux_list;
            }
        }
        return res;
    }

}