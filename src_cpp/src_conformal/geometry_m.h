#ifndef GEOMETRY_M_H
#define GEOMETRY_M_H

#include <vector>
#include "conformal_types.h"

namespace geometry_m {

using conformal_types_m::coord_t;
using conformal_types_m::side_t;
using conformal_types_m::triangle_t;

double getArea(const triangle_t& triangle);
std::vector<int> findContourCell(const std::vector<side_t>& contour);
int findContourFace(const std::vector<side_t>& contour);
double contourArea(const std::vector<side_t>& contour, int orientation = conformal_types_m::NOT_ON_FACE);

std::vector<side_t> getSidesOnFace(const std::vector<side_t>& sides, int face);
std::vector<side_t> getSidesOnEdge(const std::vector<side_t>& sides, int edge);
std::vector<side_t> getPathOnFace(const std::vector<side_t>& sides);
std::vector<side_t> buildSidesContour(const std::vector<side_t>& sides);
std::vector<side_t> buildVertexToVertexContour(const std::vector<side_t>& inner_path);
std::vector<side_t> buildVertexToSideContour(const std::vector<side_t>& inner_path);
std::vector<side_t> buildSideToVertexContour(const std::vector<side_t>& inner_path);
std::vector<side_t> buildSideToSideContour(const std::vector<side_t>& inner_path);

std::vector<triangle_t> getTrianglesOffFaces(const std::vector<triangle_t>& triangles);
std::vector<triangle_t> getTrianglesOnFace(const std::vector<triangle_t>& tris, int face);

} // namespace geometry_m

#endif
