#include <vector>
#include <cmath>
using coord_t = double;
struct side_t { coord_t init; coord_t end; std::vector<double> normal; std::vector<side_t> getSides() const { return {}; } bool isOnAnyFace() const { return false; } std::vector<int> getCell() const { return {}; } int getFace() const { return 0; } int getEdge() const { return 0; } bool isEquiv(const side_t&) const { return false; } side_t& operator=(const side_t&) { return *this; } bool operator==(const side_t&) const { return false; } };
struct triangle_t { std::vector<side_t> getSides() const { return {}; } bool isOnFace(int) const { return false; } };
struct interval_t { double x; double y; };
struct side_dir_t { side_t side; int direction; };
const int NOT_ON_FACE = -1;
const int FACE_X = 1, FACE_Y = 2, FACE_Z = 3;
const int EDGE_X = 1, EDGE_Y = 2, EDGE_Z = 3;
std::vector<double> cross(const std::vector<double>& v1, const std::vector<double>& v2) { std::vector<double> res(3); res[0]=v1[1]*v2[2]-v1[2]*v2[1]; res[1]=-(v1[0]*v2[2]-v1[2]*v2[0]); res[2]=v1[0]*v2[1]-v1[1]*v2[0]; return res; }
int abs_int(int x) { return std::abs(x); }
int floor_int(double x) { return static_cast<int>(std::floor(x)); }
double contourArea(const std::vector<side_t>&) { return 0.0; }
double getArea(const triangle_t&) { return 0.0; }
