#include <vector>
#include <string>
#include <cmath>
#include <iostream>
struct triangle_t { std::vector<int> vertices; std::vector<double> normal; };
struct cell_interval_t { int cell[3]; double interval; };
std::vector<std::vector<double>> cellIntervalsToCoords(const std::vector<cell_interval_t>&, int, int, int, int, int, int) { return {}; }
double convertInterval(double, double, double) { return 0.0; }
const int CELL_TYPE_LINEL = 1, CELL_TYPE_QUADL = 2, CELL_TYPE_TETL = 3, CELL_TYPE_HEXL = 4;
