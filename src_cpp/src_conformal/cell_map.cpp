#include <vector>
#include <memory>
#include <typeinfo>
#include <algorithm>
#include <iostream>
using coord_t = double;
struct side_t { coord_t init; coord_t end; std::vector<double> normal; std::vector<side_t> getSides() const; bool isOnAnyFace() const; std::vector<int> getCell() const; int getFace() const; int getEdge() const; bool isEquiv(const side_t&) const; side_t& operator=(const side_t&); bool operator==(const side_t&) const; };
struct triangle_t { std::vector<side_t> getSides() const; bool isOnFace(int) const; };
struct interval_t { double x; double y; };
struct side_dir_t { side_t side; int direction; };
struct fhash_tbl_t {
    void insert(int, const std::vector<int>&) {}
    void insert(int, const std::vector<std::vector<int>>&) {}
    void insert(int, const std::vector<double>&) {}
    void insert(int, const std::vector<std::vector<double>>&) {}
    void insert(int, const std::vector<side_t>&) {}
    void insert(int, const std::vector<triangle_t>&) {}
    void insert(int, const std::vector<interval_t>&) {}
    void insert(int, const std::vector<double>&, const std::vector<double>&) {}
    void insert(int, const std::vector<side_t>&, const std::vector<double>&) {}
    void insert(int, const std::vector<triangle_t>&, const std::vector<double>&) {}
    void insert(int, const std::vector<side_t>&, const std::vector<std::vector<double>>&) {}
    void insert(int, const std::vector<triangle_t>&, const std::vector<std::vector<double>>&) {}
    void insert(int, const std::vector<std::vector<int>>&, const std::vector<int>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<double>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&) {}
    void insert(int, const std::vector<side_t>&, const std::vector<side_t>&) {}
    void insert(int, const std::vector<triangle_t>&, const std::vector<triangle_t>&) {}
    void insert(int, const std::vector<interval_t>&, const std::vector<interval_t>&) {}
    void insert(int, const std::vector<double>&, const std::vector<std::vector<double>>&) {}
    void insert(int, const std::vector<std::vector<int>>&, const std::vector<std::vector<int>>&) {}
    void insert(int, const std::vector<side_t>&, const std::vector<triangle_t>&) {}
    void insert(int, const std::vector<triangle_t>&, const std::vector<side_t>&) {}
    void insert(int, const std::vector<interval_t>&, const std::vector<double>&) {}
    void insert(int, const std::vector<double>&, const std::vector<interval_t>&) {}
    void insert(int, const std::vector<double>&, const std::vector<side_t>&) {}
    void insert(int, const std::vector<double>&, const std::vector<triangle_t>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<side_t>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<triangle_t>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<interval_t>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<int>&) {}
    void insert(int, const std::vector<std::vector<double>>&, const std::vector<std::vector<int>>&) {}
};
namespace cell_map_m {
struct element_set_t { std::vector<triangle_t> triangles; std::vector<side_t> sides; std::vector<side_t> sides_on; std::vector<interval_t> intervals; };
struct cell_t { int cell[3]; };
struct cell_ratios_t { double area[3] = {1.0, 1.0, 1.0}; double length[3] = {1.0, 1.0, 1.0}; };
class cell_ratios_map_t : public fhash_tbl_t {
public:
    std::vector<cell_t> keys;
    bool hasKey(const std::vector<int>&) { return false; }
    void addFaceRatio(const cell_t&, int, double) {}
    void addEdgeRatio(const cell_t&, int, double) {}
    cell_ratios_t getCellRatiosInCell(const std::vector<int>&) { return cell_ratios_t{}; }
};
class cell_map_t : public fhash_tbl_t {
public:
    std::vector<cell_t> keys;
    bool hasKey(const std::vector<int>&) { return false; }
    std::vector<triangle_t> getTrianglesInCell(const std::vector<int>&) { return std::vector<triangle_t>(); }
    std::vector<side_t> getSidesInCell(const std::vector<int>&) { return std::vector<side_t>(); }
    std::vector<side_t> getOnSidesInCell(const std::vector<int>&) { return std::vector<side_t>(); }
    std::vector<interval_t> getIntervalsInCell(const std::vector<int>&) { return std::vector<interval_t>(); }
};
class triangle_map_t : public cell_map_t {
public:
    void addTriangle(const triangle_t&) {}
};
class interval_map_t : public cell_map_t {
public:
    void addInterval(const interval_t&) {}
};
class side_map_t : public cell_map_t {
public:
    void addSide(const side_t&) {}
    void addSideOn(const side_t&) {}
};
class side_tris_map_t : public fhash_tbl_t {
public:
    std::vector<side_dir_t> keys;
    bool hasKey(const std::vector<int>&) { return false; }
    std::vector<triangle_t> getTrianglesFromSide(const std::vector<int>&) { return std::vector<triangle_t>(); }
};
class side_triangle_map_t : public side_tris_map_t {
public:
    void addTriangleToSide(const side_t&, const triangle_t&) {}
};
void buildSideToTrisMap(side_triangle_map_t&, const std::vector<triangle_t>&) {}
side_t getSideOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_t{}; }
side_dir_t getSideDirOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_dir_t{}; }
std::vector<side_dir_t> getSideDirOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_dir_t>(); }
side_t getOnSideOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_t{}; }
side_dir_t getOnSideDirOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_dir_t{}; }
interval_t getIntervalOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return interval_t{}; }
std::vector<side_t> getIntersectingSidesOnFace(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_t>(); }
std::vector<side_t> getIntersectingSidesOnEdge(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_t>(); }
std::vector<side_t> getIntersectingSidesOnVertex(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_t>(); }
side_t getOnSideOnEdge(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_t{}; }
side_t getOnSideOnEdge(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_t{}; }
interval_t getIntervalOnEdge(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return interval_t{}; }
interval_t getIntervalOnEdge(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return interval_t{}; }
std::vector<side_t> getIntersectingSidesOnVertex(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_t>(); }
side_t getOnSideOnVertex(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_t{}; }
side_t getOnSideOnVertex(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return side_t{}; }
interval_t getIntervalOnVertex(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return interval_t{}; }
interval_t getIntervalOnVertex(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return interval_t{}; }
std::vector<triangle_t> getTrianglesInCell(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<triangle_t>(); }
std::vector<side_t> getSidesInCell(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_t>(); }
std::vector<side_t> getOnSidesInCell(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<side_t>(); }
std::vector<interval_t> getIntervalsInCell(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return std::vector<interval_t>(); }
cell_ratios_t getCellRatiosInCell(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return cell_ratios_t{}; }
bool cellRatioHasKey(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int) { return false; }
void addFaceRatio(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, double) {}
void addEdgeRatio(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, double) {}
void addTriangle(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, const triangle_t&) {}
void addSide(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, const side_t&) {}
void addInterval(int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, int, const interval_t&) {}
bool isNewSide(const std::vector<side_t>&, const side_t&) { return false; }
bool isNewSide(const std::vector<side_t>&, const side_t&, int) { return false; }
bool isNewSide(const std::vector<side_t>&, const side_t&, int, int) { return false; }
bool isNewSide(const std::vector<side_t>&, const side_t&, int, int, int) { return false; }
bool isNewSide(const std::vector<side_t>&, const side_t&, int, int, int, int) { return false; }
bool isNewSide(const std::vector<side_t>&, const side_t&, int, int, int, int, int) { return false; }
bool isNewSide(const std::vector<side_t>&, const side_t&, int, int, int, int, int, int) { return false; }
void addFaceRatio(cell_ratios_map_t&, const cell_t&, int, double) {}
void addEdgeRatio(cell_ratios_map_t&, const cell_t&, int, double) {}
cell_ratios_t getCellRatiosInCell(cell_ratios_map_t&, const std::vector<int>&) { return cell_ratios_t{}; }
bool cell_ratio_hasKey(cell_ratios_map_t&, const std::vector<int>&) { return false; }
void addTriangle(cell_map_t&, const triangle_t&) {}
void addSide(cell_map_t&, const side_t&) {}
void addInterval(cell_map_t&, const interval_t&) {}
void addFaceRatio(cell_map_t&, const cell_t&, int, double) {}
void addEdgeRatio(cell_map_t&, const cell_t&, int, double) {}
std::vector<triangle_t> getTrianglesInCell(cell_map_t&, const std::vector<int>&) { return std::vector<triangle_t>(); }
std::vector<side_t> getSidesInCell(cell_map_t&, const std::vector<int>&) { return std::vector<side_t>(); }
std::vector<side_t> getOnSidesInCell(cell_map_t&, const std::vector<int>&) { return std::vector<side_t>(); }
std::vector<interval_t> getIntervalsInCell(cell_map_t&, const std::vector<int>&) { return std::vector<interval_t>(); }
bool cell_hasKey(cell_map_t&, const std::vector<int>&) { return false; }
void addTriangle(side_tris_map_t&, const side_t&, const triangle_t&) {}
void addSideOn(side_tris_map_t&, const side_t&) {}
std::vector<triangle_t> getTrianglesFromSide(side_tris_map_t&, const std::vector<int>&) { return std::vector<triangle_t>(); }
bool side_hasKey(side_tris_map_t&, const std::vector<int>&) { return false; }
}
