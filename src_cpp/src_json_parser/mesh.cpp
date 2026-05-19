#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <algorithm>
struct triangle_t { std::vector<int> vertices; std::vector<double> normal; };
struct cell_interval_t { int cell[3]; double interval; };
struct side_t { double init; double end; std::vector<double> normal; };
struct triangle_t2 { std::vector<int> vertices; };
struct interval_t { double x; double y; };
struct side_dir_t { side_t side; int direction; };
struct fhash_tbl_t { void insert(int, const std::vector<int>&) {} void insert(int, const std::vector<std::vector<int>>&) {} void insert(int, const std::vector<double>&) {} void insert(int, const std::vector<std::vector<double>>&) {} };
namespace mesh_m {
    void addCoordinates(std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&) {}
    void addElements(std::vector<std::vector<int>>&, const std::vector<std::vector<int>>&) {}
    std::vector<std::vector<double>> readCoordinates(const std::string&) { return {}; }
    std::vector<std::vector<int>> readElements(const std::string&) { return {}; }
    std::vector<triangle_t> readTriangles(const std::string&) { return {}; }
    std::vector<cell_interval_t> readCellIntervals(const std::string&) { return {}; }
}
