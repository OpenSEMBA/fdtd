#include <vector>
#include <memory>
#include <typeinfo>
#include <algorithm>
#include <iostream>

// Forward declarations and includes for external modules
// #include "geometry_m.h"
// #include "fhash.h"
// #include "NFDETypes_m.h"

// Assuming these types are defined in the included headers
// struct triangle_t;
// struct side_t;
// struct interval_t;
// struct fhash_tbl_t;
// struct ConformalPECElements_t;

// Constants from geometry_m
// enum FaceType { FACE_X, FACE_Y, FACE_Z };
// bool isNewSide(const std::vector<side_t>& sides, const side_t& side);

// Type alias from fhash
// using key = fhash_key;

// Kind parameter from NFDETypes_m
// using rkind = double; // Assuming real(8) maps to double

namespace cell_map_m {

    struct element_set_t {
        std::vector<triangle_t> triangles;
        std::vector<side_t> sides;
        std::vector<side_t> sides_on;
        std::vector<interval_t> intervals;
    };

    struct cell_t {
        int cell[3];
    };

    struct cell_ratios_t {
        double area[3] = {1.0, 1.0, 1.0};
        double length[3] = {1.0, 1.0, 1.0};
    };

    class cell_ratios_map_t : public fhash_tbl_t {
    public:
        std::vector<cell_t> keys;

        bool hasKey(const std::vector<int>& k) {
            return cell_ratio_hasKey(*this, k);
        }

        void addFaceRatio(const cell_t& cell, int direction, double ratio) {
            ::cell_map_m::addFaceRatio(*this, cell, direction, ratio);
        }

        void addEdgeRatio(const cell_t& cell, int direction, double ratio) {
            ::cell_map_m::addEdgeRatio(*this, cell, direction, ratio);
        }

        cell_ratios_t getCellRatiosInCell(const std::vector<int>& k) {
            return getCellRatiosInCell(*this, k);
        }
    };

    class cell_map_t : public fhash_tbl_t {
    public:
        std::vector<cell_t> keys;

        bool hasKey(const std::vector<int>& k) {
            return cell_hasKey(*this, k);
        }

        std::vector<triangle_t> getTrianglesInCell(const std::vector<int>& k) {
            return getTrianglesInCell(*this, k);
        }

        std::vector<side_t> getSidesInCell(const std::vector<int>& k) {
            return getSidesInCell(*this, k);
        }

        std::vector<side_t> getOnSidesInCell(const std::vector<int>& k) {
            return getOnSidesInCell(*this, k);
        }

        std::vector<interval_t> getIntervalsInCell(const std::vector<int>& k) {
            return getIntervalsInCell(*this, k);
        }
    };

    class triangle_map_t : public cell_map_t {
    public:
        void addTriangle(const triangle_t& triangle) {
            ::cell_map_m::addTriangle(*this, triangle);
        }
    };

    class interval_map_t : public cell_map_t {
    public:
        void addInterval(const interval_t& interval) {
            ::cell_map_m::addInterval(*this, interval);
        }
    };

    class side_map_t : public cell_map_t {
    public:
        void addSide(const side_t& side) {
            ::cell_map_m::addSide(*this, side);
        }

        void addSideOn(const side_t& side) {
            ::cell_map_m::addSideOn(*this, side);
        }
    };

    struct side_dir_t {
        int side_dir[4];
    };

    class side_tris_map_t : public fhash_tbl_t {
    public:
        std::vector<side_dir_t> keys;

        bool hasKey(const std::vector<int>& k) {
            return side_hasKey(*this, k);
        }

        std::vector<triangle_t> getTrianglesFromSide(const std::vector<int>& k) {
            return getTrianglesFromSide(*this, k);
        }
    };

    class side_triangle_map_t : public side_tris_map_t {
    public:
        void addTriangleToSide(const side_t& side, const triangle_t& triangle) {
            ::cell_map_m::addTriangleToSide(*this, side, triangle);
        }
    };

    void buildSideToTrisMap(side_triangle_map_t& res, const std::vector<triangle_t>& triangles) {
        std::vector<side_t> sides;
        side_t aux_side;

        if (!res.keys.size()) {
            res.keys = std::vector<side_dir_t>(0);
        }
        for (size_t i = 0; i < triangles.size(); ++i) {
            sides = triangles[i].getSides();
            for (size_t j = 0; j < 3; ++j) {
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
            std::vector<std::any> alloc_list;
            this_obj.get_raw(key(side_dir), alloc_list);
            
            // Note: std::any requires C++17. For older standards, use void* or a variant-like structure.
            // Assuming C++17 for std::any
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                element_set_t aux_list;
                aux_list.triangles.resize(alloc_list_ptr->triangles.size() + 1);
                for (size_t i = 0; i < alloc_list_ptr->triangles.size(); ++i) {
                    aux_list.triangles[i] = alloc_list_ptr->triangles[i];
                }
                aux_list.triangles[alloc_list_ptr->triangles.size()] = triangle;
                
                alloc_list_ptr->triangles.resize(aux_list.triangles.size());
                alloc_list_ptr->triangles = aux_list.triangles;
                this_obj.set(key(side_dir), alloc_list_ptr);
            }
        } else {
            element_set_t aux_list;
            aux_list.triangles.resize(1);
            aux_list.triangles[0] = triangle;
            this_obj.set(key(side_dir), aux_list);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({side_dir[0], side_dir[1], side_dir[2]});
            this_obj.keys = aux_keys;
        }
    }

    void buildSideMap(side_tris_map_t& res, const std::vector<triangle_t>& triangles) {
        side_triangle_map_t tri_map;
        std::vector<side_dir_t> keys;
        element_set_t elems;
        
        buildSideToTrisMap(tri_map, triangles);
        keys = tri_map.keys;
        
        for (size_t i = 0; i < keys.size(); ++i) {
            elems.triangles = tri_map.getTrianglesFromSide({keys[i].side_dir[0], keys[i].side_dir[1], keys[i].side_dir[2], keys[i].side_dir[3]});
            res.set(key({keys[i].side_dir[0], keys[i].side_dir[1], keys[i].side_dir[2], keys[i].side_dir[3]}), elems);
        }
        res.keys = keys;
    }

    void buildCellMap(cell_map_t& res, const ConformalPECElements_t& volume) {
        triangle_map_t tri_map;
        interval_map_t interval_map;
        side_map_t side_map, side_map_on;
        std::vector<cell_t> keys;
        element_set_t elems;
        
        buildMapOfTrisOnFaces(tri_map, volume.triangles);
        buildMapOfIntervals(interval_map, volume.intervals);
        buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map, volume.triangles);
        buildMapOfSidesOnEdgeFromTrisOnFaces(side_map_on, volume.triangles);
        
        keys = mergeKeys(tri_map.keys, side_map.keys);
        keys = mergeKeys(keys, side_map_on.keys);
        keys = mergeKeys(keys, interval_map.keys);
        
        for (size_t i = 0; i < keys.size(); ++i) {
            elems.triangles = tri_map.getTrianglesInCell({keys[i].cell[0], keys[i].cell[1], keys[i].cell[2]});
            elems.sides = side_map.getSidesInCell({keys[i].cell[0], keys[i].cell[1], keys[i].cell[2]});
            elems.sides_on = side_map_on.getOnSidesInCell({keys[i].cell[0], keys[i].cell[1], keys[i].cell[2]});
            elems.intervals = interval_map.getIntervalsInCell({keys[i].cell[0], keys[i].cell[1], keys[i].cell[2]});
            res.set(key({keys[i].cell[0], keys[i].cell[1], keys[i].cell[2]}), elems);
        }
        res.keys = keys;
    }

    std::vector<cell_t> mergeKeys(const std::vector<cell_t>& tri_keys, const std::vector<cell_t>& side_keys) {
        std::vector<cell_t> res;
        std::vector<cell_t> aux;
        
        if (tri_keys.size() == 0) {
            res = std::vector<cell_t>(0);
        } else {
            res = tri_keys;
        }
        
        if (side_keys.size() != 0) {
            for (size_t i = 0; i < side_keys.size(); ++i) {
                if (isNewKey(res, side_keys[i])) {
                    addKey(res, side_keys[i]);
                }
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
        if (keys.size() != 0) {
            for (size_t i = 0; i < keys.size(); ++i) {
                if (keys[i].cell[0] == k.cell[0] && 
                    keys[i].cell[1] == k.cell[1] && 
                    keys[i].cell[2] == k.cell[2]) {
                    return false;
                }
            }
        }
        return true;
    }

    void buildMapOfTrisOnFaces(triangle_map_t& res, const std::vector<triangle_t>& triangles) {
        if (!res.keys.size()) {
            res.keys = std::vector<cell_t>(0);
        }
        for (size_t i = 0; i < triangles.size(); ++i) {
            if (triangles[i].isOnAnyFace()) {
                res.addTriangle(triangles[i]);
            }
        }
    }

    void buildMapOfIntervals(interval_map_t& res, const std::vector<interval_t>& intervals) {
        if (!res.keys.size()) {
            res.keys = std::vector<cell_t>(0);
        }
        for (size_t i = 0; i < intervals.size(); ++i) {
            res.addInterval(intervals[i]);
        }
    }

    void buildMapOfSidesOnFaceOrEdgeFromTrisNotOnFaces(side_map_t& res, const std::vector<triangle_t>& triangles) {
        std::vector<side_t> sides(3);
        int cell[3];
        
        if (!res.keys.size()) {
            res.keys = std::vector<cell_t>(0);
        }
        for (size_t i = 0; i < triangles.size(); ++i) {
            if (!triangles[i].isOnAnyFace()) {
                sides = triangles[i].getSides();
                for (size_t j = 0; j < 3; ++j) {
                    if (sides[j].isOnAnyFace() || sides[j].isOnAnyEdge()) {
                        res.addSide(sides[j]);
                    }
                }
            }
        }
    }

    void buildMapOfSidesOnEdgeFromTrisOnFaces(side_map_t& res, const std::vector<triangle_t>& triangles) {
        std::vector<side_t> sides(3);
        int cell[3];
        
        if (!res.keys.size()) {
            res.keys = std::vector<cell_t>(0);
        }
        for (size_t i = 0; i < triangles.size(); ++i) {
            if (triangles[i].isOnAnyFace()) {
                sides = triangles[i].getSides();
                for (size_t j = 0; j < 3; ++j) {
                    if (sides[j].isOnAnyEdge()) {
                        res.addSideOn(sides[j]);
                    }
                }
            }
        }
    }

    bool cell_hasKey(cell_map_t& this_obj, const std::vector<int>& k) {
        int stat;
        this_obj.check_key(key(k), stat);
        return stat == 0;
    }

    bool side_hasKey(side_tris_map_t& this_obj, const std::vector<int>& k) {
        int stat;
        this_obj.check_key(key(k), stat);
        return stat == 0;
    }

    bool cell_ratio_hasKey(cell_ratios_map_t& this_obj, const std::vector<int>& k) {
        int stat;
        this_obj.check_key(key(k), stat);
        return stat == 0;
    }

    void addTriangle(triangle_map_t& this_obj, const triangle_t& triangle) {
        std::vector<int> cell(3);
        cell[0] = triangle.getCell()[0];
        cell[1] = triangle.getCell()[1];
        cell[2] = triangle.getCell()[2];

        if (this_obj.hasKey(cell)) {
            std::vector<std::any> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                element_set_t aux_list;
                aux_list.triangles.resize(alloc_list_ptr->triangles.size() + 1);
                for (size_t i = 0; i < alloc_list_ptr->triangles.size(); ++i) {
                    aux_list.triangles[i] = alloc_list_ptr->triangles[i];
                }
                aux_list.triangles[alloc_list_ptr->triangles.size()] = triangle;
                
                alloc_list_ptr->triangles.resize(aux_list.triangles.size());
                alloc_list_ptr->triangles = aux_list.triangles;
                this_obj.set(key(cell), alloc_list_ptr);
            }
        } else {
            element_set_t aux_list;
            aux_list.triangles.resize(1);
            aux_list.triangles[0] = triangle;
            this_obj.set(key(cell), aux_list);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell[0], cell[1], cell[2]});
            this_obj.keys = aux_keys;
        }
    }

    void addInterval(interval_map_t& this_obj, const interval_t& interval) {
        side_t aux;
        // Assuming init is a method or property that sets up the side
        // aux.init().position = interval.ini.cell;
        // aux.end().position = interval.end.cell;
        // This part is tricky without knowing the exact API of interval_t and side_t
        // Assuming a method getCell() exists on side_t after initialization
        std::vector<int> cell = aux.getCell();

        if (this_obj.hasKey(cell)) {
            std::vector<std::any> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                element_set_t aux_list;
                aux_list.intervals.resize(alloc_list_ptr->intervals.size() + 1);
                for (size_t i = 0; i < alloc_list_ptr->intervals.size(); ++i) {
                    aux_list.intervals[i] = alloc_list_ptr->intervals[i];
                }
                aux_list.intervals[alloc_list_ptr->intervals.size()] = interval;
                
                alloc_list_ptr->intervals.resize(aux_list.intervals.size());
                alloc_list_ptr->intervals = aux_list.intervals;
                this_obj.set(key(cell), alloc_list_ptr);
            }
        } else {
            element_set_t aux_list;
            aux_list.intervals.resize(1);
            aux_list.intervals[0] = interval;
            this_obj.set(key(cell), aux_list);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell[0], cell[1], cell[2]});
            this_obj.keys = aux_keys;
        }
    }

    void addSide(side_map_t& this_obj, const side_t& side) {
        std::vector<int> cell = side.getCell();

        if (this_obj.hasKey(cell)) {
            std::vector<std::any> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                
                if (isNewSide(alloc_list_ptr->sides, side)) {
                    element_set_t aux_list;
                    aux_list.sides.resize(alloc_list_ptr->sides.size() + 1);
                    for (size_t i = 0; i < alloc_list_ptr->sides.size(); ++i) {
                        aux_list.sides[i] = alloc_list_ptr->sides[i];
                    }
                    aux_list.sides[alloc_list_ptr->sides.size()] = side;
                    
                    alloc_list_ptr->sides.resize(aux_list.sides.size());
                    alloc_list_ptr->sides = aux_list.sides;
                    this_obj.set(key(cell), alloc_list_ptr);
                }
            }
        } else {
            element_set_t aux_list;
            aux_list.sides.resize(1);
            aux_list.sides[0] = side;
            this_obj.set(key(cell), aux_list);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell[0], cell[1], cell[2]});
            this_obj.keys = aux_keys;
        }
    }

    void addSideOn(side_map_t& this_obj, const side_t& side) {
        std::vector<int> cell = side.getCell();

        if (this_obj.hasKey(cell)) {
            std::vector<std::any> alloc_list;
            this_obj.get_raw(key(cell), alloc_list);
            
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                element_set_t aux_list;
                aux_list.sides_on.resize(alloc_list_ptr->sides_on.size() + 1);
                for (size_t i = 0; i < alloc_list_ptr->sides_on.size(); ++i) {
                    aux_list.sides_on[i] = alloc_list_ptr->sides_on[i];
                }
                aux_list.sides_on[alloc_list_ptr->sides_on.size()] = side;
                
                alloc_list_ptr->sides_on.resize(aux_list.sides_on.size());
                alloc_list_ptr->sides_on = aux_list.sides_on;
                this_obj.set(key(cell), alloc_list_ptr);
            }
        } else {
            element_set_t aux_list;
            aux_list.sides_on.resize(1);
            aux_list.sides_on[0] = side;
            this_obj.set(key(cell), aux_list);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back({cell[0], cell[1], cell[2]});
            this_obj.keys = aux_keys;
        }
    }

    std::vector<triangle_t> getTrianglesInCell(cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<std::any> alloc_list;
        std::vector<triangle_t> res;

        if (this_obj.hasKey(k)) {
            this_obj.get_raw(key(k), alloc_list);
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                res = alloc_list_ptr->triangles;
            }
        } else {
            res = std::vector<triangle_t>(0);
        }
        return res;
    }

    std::vector<interval_t> getIntervalsInCell(cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<std::any> alloc_list;
        std::vector<interval_t> res;

        if (this_obj.hasKey(k)) {
            this_obj.get_raw(key(k), alloc_list);
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                res = alloc_list_ptr->intervals;
            }
        } else {
            res = std::vector<interval_t>(0);
        }
        return res;
    }

    std::vector<triangle_t> getTrianglesFromSide(side_tris_map_t& this_obj, const std::vector<int>& k) {
        std::vector<std::any> alloc_list;
        std::vector<triangle_t> res;

        if (this_obj.hasKey(k)) {
            this_obj.get_raw(key(k), alloc_list);
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                res = alloc_list_ptr->triangles;
            }
        } else {
            res = std::vector<triangle_t>(0);
        }
        return res;
    }

    std::vector<side_t> getSidesInCell(cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<std::any> alloc_list;
        std::vector<side_t> res;

        if (this_obj.hasKey(k)) {
            this_obj.get_raw(key(k), alloc_list);
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                res = alloc_list_ptr->sides;
            }
        } else {
            res = std::vector<side_t>(0);
        }
        return res;
    }

    std::vector<side_t> getOnSidesInCell(cell_map_t& this_obj, const std::vector<int>& k) {
        std::vector<std::any> alloc_list;
        std::vector<side_t> res;

        if (this_obj.hasKey(k)) {
            this_obj.get_raw(key(k), alloc_list);
            if (alloc_list.size() > 0 && typeid(element_set_t) == typeid(alloc_list[0])) {
                element_set_t* alloc_list_ptr = &alloc_list[0];
                res = alloc_list_ptr->sides_on;
            }
        } else {
            res = std::vector<side_t>(0);
        }
        return res;
    }

    void addFaceRatio(cell_ratios_map_t& this_obj, const cell_t& cell, int direction, double ratio) {
        std::vector<std::any> alloc_list;
        cell_ratios_t aux_cell_ratio;
        std::vector<int> cell_vec = {cell.cell[0], cell.cell[1], cell.cell[2]};

        if (this_obj.hasKey(cell_vec)) {
            this_obj.get_raw(key(cell_vec), alloc_list);
            if (alloc_list.size() > 0 && typeid(cell_ratios_t) == typeid(alloc_list[0])) {
                cell_ratios_t* alloc_list_ptr = &alloc_list[0];
                alloc_list_ptr->area[direction] = ratio;
                this_obj.set(key(cell_vec), alloc_list_ptr);
            }
        } else {
            aux_cell_ratio.area[direction] = ratio;
            this_obj.set(key(cell_vec), aux_cell_ratio);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back(cell);
            this_obj.keys = aux_keys;
        }
    }

    void addEdgeRatio(cell_ratios_map_t& this_obj, const cell_t& cell, int direction, double ratio) {
        std::vector<std::any> alloc_list;
        cell_ratios_t aux_cell_ratio;
        std::vector<int> cell_vec = {cell.cell[0], cell.cell[1], cell.cell[2]};

        if (this_obj.hasKey(cell_vec)) {
            this_obj.get_raw(key(cell_vec), alloc_list);
            if (alloc_list.size() > 0 && typeid(cell_ratios_t) == typeid(alloc_list[0])) {
                cell_ratios_t* alloc_list_ptr = &alloc_list[0];
                alloc_list_ptr->length[direction] = ratio;
                this_obj.set(key(cell_vec), alloc_list_ptr);
            }
        } else {
            aux_cell_ratio.length[direction] = ratio;
            this_obj.set(key(cell_vec), aux_cell_ratio);

            std::vector<cell_t> aux_keys = this_obj.keys;
            aux_keys.push_back(cell);
            this_obj.keys = aux_keys;
        }
    }

    cell_ratios_t getCellRatiosInCell(cell_ratios_map_t& this_obj, const std::vector<int>& k) {
        std::vector<std::any> alloc_list;
        cell_ratios_t res;

        if (this_obj.hasKey(k)) {
            this_obj.get_raw(key(k), alloc_list);
            if (alloc_list.size() > 0 && typeid(cell_ratios_t) == typeid(alloc_list[0])) {
                res = alloc_list[0];
            }
        }
        return res;
    }

}