#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <array>
#include <memory>

// Forward declarations and includes for external modules/types
// These would typically be in their respective headers
// #include "geometry_m.h"
// #include "cell_map_m.h"
// #include "NFDETypes_m.h"

// Assuming NFDETypes_m provides these types and rkind
// If not, they need to be defined or included from NFDETypes_m.h
// For the sake of translation, we assume the following definitions exist or are included:
// using rkind = double; // Assuming rkind is double precision
// struct ConformalPECRegions_t;
// struct ConformalPECElements_t;
// struct ConformalMedia_t;
// struct edge_t;
// struct face_t;
// struct conformal_face_media_t;
// struct conformal_edge_media_t;
// struct side_tris_map_t;
// struct cell_map_t;
// struct side_t;
// struct triangle_t;
// struct interval_t;
// struct coord_t;
// struct cell_ratios_map_t;
// struct cell_ratios_t;

// Constants from geometry_m or similar
// #define FACE_X 0
// #define FACE_Y 1
// #define FACE_Z 2
// #define EDGE_X 0
// #define EDGE_Y 1
// #define EDGE_Z 2
// #define NOT_ON_EDGE -1

namespace conformal_m {

   constexpr double EDGE_RATIO_EQ_TOLERANCE = 1e-5;
   constexpr double FACE_RATIO_EQ_TOLERANCE = 1e-3;

   // Helper functions that are likely defined in other modules or need to be stubbed/declared
   // Since the prompt asks to translate THIS file, we assume these functions are available
   // in the linked environment or other headers. We will declare them here if they are
   // part of the interface or defined elsewhere. For a complete translation, these would
   // be in their respective .h files.

   // Declarations for functions called but not defined in this module
   void buildCellMap(cell_map_t& cell_map, const ConformalPECElements_t& volume);
   void fillElements(const cell_map_t& cell_map, std::vector<face_t>& faces, std::vector<edge_t>& edges);
   void addRatio(std::vector<double>& ratios, double ratio);
   bool isNewRatio(const std::vector<double>& ratios, double ratio, double tol);
   std::vector<conformal_edge_media_t> addEdgeMedia(const std::vector<edge_t>& edges, const std::vector<double>& edge_ratios);
   std::vector<conformal_face_media_t> addFaceMedia(const std::vector<face_t>& faces, const std::vector<double>& face_ratios);
   std::vector<edge_t> filterEdgesByMedia(const std::vector<edge_t>& edges, double ratio);
   std::vector<face_t> filterFacesByMedia(const std::vector<face_t>& faces, double ratio);
   void buildSideMap(side_tris_map_t& map, const std::vector<triangle_t>& triangles);
   ConformalMedia_t buildConformalVolume(const ConformalPECElements_t& volume);
   
   // Helper functions used in computeTimeStepScalingFactor
   // These are methods of cell_ratios_map_t, so they are declared within the struct if it's a class
   // or as free functions if it's a struct with methods. Assuming struct with methods.

   // Helper functions used in fillElements
   std::vector<side_t> getSidesInCell(const cell_map_t& cell_map, const std::array<int, 3>& cell);
   std::vector<triangle_t> getTrianglesInCell(const cell_map_t& cell_map, const std::array<int, 3>& cell);
   std::vector<interval_t> getIntervalsInCell(const cell_map_t& cell_map, const std::array<int, 3>& cell);
   std::vector<side_t> getSidesOnFace(const std::vector<side_t>& sides, int face);
   std::vector<side_t> findLargestContour(const std::vector<side_t>& sides);
   void fillFaceFromContour(const std::vector<side_t>& contour, std::vector<face_t>& faces);
   void fillEdgesFromContour(const std::vector<side_t>& contour, std::vector<edge_t>& edges);
   std::vector<triangle_t> getTrianglesOnFace(const std::vector<triangle_t>& tris, int face);
   void fillFullFaces(const std::vector<triangle_t>& tris_on_face, std::vector<face_t>& faces, std::vector<edge_t>& edges);
   void fillIntervals(const std::vector<interval_t>& intervals, std::vector<edge_t>& edges, std::vector<face_t>& faces);
   std::vector<side_t> getSidesOnEdge(const std::vector<side_t>& sides, int edge);
   void fillEdges(const std::vector<side_t>& sides, std::vector<edge_t>& edges);
   double getArea(const triangle_t& tri);
   void addFace(std::vector<face_t>& faces, const std::array<int, 3>& cell, int face, double ratio);
   bool isNewEdge(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, double ratio);
   void addEdge(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, const side_t& side);
   bool isEdgeFilled(const std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge);
   void fillSmallerRatio(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, const side_t& side);
   void reduceEdgeRatio(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, const side_t& side);
   void fillEdgesFromInterval(std::vector<edge_t>& edges, const interval_t& interval);
   void fillFaceFromInterval(std::vector<face_t>& faces, const interval_t& interval);
   std::array<int, 3> findContourCell(const std::vector<side_t>& contour);
   int findContourFace(const std::vector<side_t>& contour);
   double contourArea(const std::vector<side_t>& contour);
   std::vector<side_t> buildSidesContour(const std::vector<side_t>& sides);
   std::vector<side_t> buildSidesFromCellInterval(const interval_t& interval);

   // eq_ratio is defined in this module
   bool eq_ratio(double a, double b, double tol) {
       return std::abs(a - b) < tol;
   }

   std::vector<side_tris_map_t> buildSideMaps(const ConformalPECRegions_t& regions) {
       std::vector<side_tris_map_t> res(regions.volumes.size());
       for (size_t i = 0; i < regions.volumes.size(); ++i) {
           buildSideMap(res[i], regions.volumes[i].triangles);
       }
       return res;
   }

   std::vector<ConformalMedia_t> buildConformalMedia(const ConformalPECRegions_t& regions) {
       std::vector<ConformalMedia_t> res(regions.volumes.size());
       for (size_t i = 0; i < regions.volumes.size(); ++i) {
           res[i] = buildConformalVolume(regions.volumes[i]);
       }
       return res;
   }

   ConformalMedia_t buildConformalVolume(const ConformalPECElements_t& volume) {
       ConformalMedia_t res;
       cell_map_t cell_map;
       std::vector<double> edge_ratios, face_ratios;
       std::vector<edge_t> edges, filtered_edges;
       std::vector<face_t> faces, filtered_faces;

       buildCellMap(cell_map, volume);
       fillElements(cell_map, faces, edges);
       addNewRatios(edges, faces, edge_ratios, face_ratios);
       
       // Note: addEdgeMedia and addFaceMedia return vectors of pointers in Fortran
       // In C++, we need to handle memory management. 
       // Assuming conformal_edge_media_t and conformal_face_media_t have dynamic arrays
       // The Fortran code uses pointers for res. In C++, we can use std::vector of structs
       // if the structs contain vectors, or std::vector of unique_ptr if they contain raw pointers.
       // Given the complexity, let's assume the structs have std::vector members for edges/faces
       // or we manage the pointers carefully.
       // For simplicity and safety, let's assume the structs have std::vector members.
       // If they have raw pointers, we need to delete them.
       
       // Let's assume the following structure for conformal_edge_media_t and conformal_face_media_t
       // based on usage:
       // struct conformal_edge_media_t {
       //     std::vector<edge_t> edges;
       //     double ratio;
       //     int n_elements;
       // };
       // struct conformal_face_media_t {
       //     std::vector<face_t> faces;
       //     double ratio;
       //     int n_elements;
       // };

       // If the original Fortran uses pointers, we might need to use std::vector<std::unique_ptr<...>>
       // However, without the definition of these types, we'll assume they are value types with internal vectors.

       auto edge_media = addEdgeMedia(edges, edge_ratios);
       auto face_media = addFaceMedia(faces, face_ratios);

       res.edge_media = edge_media; // Assuming edge_media is a vector of conformal_edge_media_t
       res.face_media = face_media; // Assuming face_media is a vector of conformal_face_media_t

       res.n_edges_media = res.edge_media.size();
       res.n_faces_media = res.face_media.size();

       res.time_step_scale_factor = computeTimeStepScalingFactor(res.edge_media, res.face_media);
       res.tag = volume.tag;
       return res;
   }

   void addNewRatios(const std::vector<edge_t>& edges, const std::vector<face_t>& faces, std::vector<double>& edge_ratios, std::vector<double>& face_ratios) {
       edge_ratios.clear();
       face_ratios.clear();
       for (size_t i = 0; i < edges.size(); ++i) {
           if (isNewRatio(edge_ratios, edges[i].ratio, EDGE_RATIO_EQ_TOLERANCE)) {
               addRatio(edge_ratios, edges[i].ratio);
           }
       }
       for (size_t i = 0; i < faces.size(); ++i) {
           if (isNewRatio(face_ratios, faces[i].ratio, FACE_RATIO_EQ_TOLERANCE)) {
               addRatio(face_ratios, faces[i].ratio);
           }
       }
   }

   std::vector<conformal_edge_media_t> addEdgeMedia(const std::vector<edge_t>& edges, const std::vector<double>& edge_ratios) {
       std::vector<conformal_edge_media_t> res(edge_ratios.size());
       for (size_t i = 0; i < edge_ratios.size(); ++i) {
           auto filtered_edges = filterEdgesByMedia(edges, edge_ratios[i]);
           res[i].edges = filtered_edges;
           res[i].ratio = edge_ratios[i];
           res[i].n_elements = filtered_edges.size();
       }
       return res;
   }

   std::vector<conformal_face_media_t> addFaceMedia(const std::vector<face_t>& faces, const std::vector<double>& face_ratios) {
       std::vector<conformal_face_media_t> res(face_ratios.size());
       for (size_t i = 0; i < face_ratios.size(); ++i) {
           auto filtered_faces = filterFacesByMedia(faces, face_ratios[i]);
           res[i].faces = filtered_faces;
           res[i].ratio = face_ratios[i];
           res[i].n_elements = filtered_faces.size();
       }
       return res;
   }

   double computeTimeStepScalingFactor(const std::vector<conformal_face_media_t>& faces_media, const std::vector<conformal_edge_media_t>& edges_media) {
       double res = 1.0;
       cell_ratios_map_t cell_ratio_map;
       
       // Initialize cell_ratio_map%keys if it's a vector
       if (cell_ratio_map.keys.empty()) {
           cell_ratio_map.keys = std::vector<std::array<int, 3>>();
       }

       for (size_t i = 0; i < faces_media.size(); ++i) {
           for (size_t j = 0; j < faces_media[i].faces.size(); ++j) {
               cell_ratio_map.addFaceRatio(faces_media[i].faces[j].cell, faces_media[i].faces[j].direction, faces_media[i].faces[j].ratio);
           }
       }
       for (size_t i = 0; i < edges_media.size(); ++i) {
           for (size_t j = 0; j < edges_media[i].edges.size(); ++j) {
               cell_ratio_map.addEdgeRatio(edges_media[i].edges[j].cell, edges_media[i].edges[j].direction, edges_media[i].edges[j].ratio);
           }
       }

       double l_ratio = 0.0;
       for (size_t i = 0; i < cell_ratio_map.keys.size(); ++i) {
           auto cell = cell_ratio_map.keys[i].cell;
           auto aux_cell = cell;
           auto cell_ratio_info = cell_ratio_map.getCellRatiosInCell(cell);
           
           for (int j = FACE_X; j <= FACE_Z; ++j) {
               double area = cell_ratio_info.area[j];
               int idx1 = (j % 3) + 1;
               int idx2 = (j + 1) % 3 + 1;
               
               l_ratio = std::max(cell_ratio_info.length[idx1], cell_ratio_info.length[idx2]);
               
               aux_cell[idx1] = aux_cell[idx1] + 1;
               if (cell_ratio_map.hasKey(aux_cell)) {
                   cell_ratio_info = cell_ratio_map.getCellRatiosInCell(aux_cell);
                   l_ratio = std::max(l_ratio, cell_ratio_info.length[idx2]);
               } else {
                   l_ratio = 1.0;
               }
               
               aux_cell = cell;
               aux_cell[idx2] = aux_cell[idx2] + 1;
               if (cell_ratio_map.hasKey(aux_cell)) {
                   cell_ratio_info = cell_ratio_map.getCellRatiosInCell(aux_cell);
                   l_ratio = std::max(l_ratio, cell_ratio_info.length[idx1]);
               } else {
                   l_ratio = 1.0;
               }
               
               if (area != 0.0 && l_ratio != 0.0) {
                   res = std::min(res, std::sqrt(area / l_ratio));
               }
           }
       }
       return res;
   }

   void fillElements(const cell_map_t& cell_map, std::vector<face_t>& faces, std::vector<edge_t>& edges) {
       faces.clear();
       edges.clear();
       
       for (size_t i = 0; i < cell_map.keys.size(); ++i) {
           auto cell = cell_map.keys[i].cell;
           auto sides = cell_map.getSidesInCell(cell);
           auto tris = cell_map.getTrianglesInCell(cell);
           auto intervals = cell_map.getIntervalsInCell(cell);
           
           for (int face = FACE_X; face <= FACE_Z; ++face) {
               auto sides_on_face = getSidesOnFace(sides, face);
               if (!sides_on_face.empty()) {
                   auto contour = findLargestContour(sides_on_face);
                   fillFaceFromContour(contour, faces);
                   fillEdgesFromContour(contour, edges);
               }
               auto tris_on_face = getTrianglesOnFace(tris, face);
               if (!tris_on_face.empty()) {
                   fillFullFaces(tris_on_face, faces, edges);
               }
           }
           fillIntervals(intervals, edges, faces);
       }

       for (size_t i = 0; i < cell_map.keys.size(); ++i) {
           auto cell = cell_map.keys[i].cell;
           auto sides = cell_map.getSidesInCell(cell);
           
           for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
               auto sides_on_edge = getSidesOnEdge(sides, edge);
               fillEdges(sides_on_edge, edges);
           }
           
           auto on_sides = cell_map.getOnSidesInCell(cell);
           for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
               auto sides_on_edge = getSidesOnEdge(on_sides, edge);
               fillEdges(sides_on_edge, edges);
           }
       }
   }

   std::vector<side_t> buildSidesFromCellInterval(const interval_t& interval) {
       std::vector<side_t> res(4);
       side_t aux;
       aux.init.position = interval.ini.cell;
       aux.end.position = interval.end.cell;
       int face = aux.getFace();
       std::array<int, 3> cs1 = aux.getCell();
       std::array<int, 3> cs2, cs3, cs4;
       
       switch(face) {
           case FACE_X:
               cs2 = cs1; cs2[1] += 1;
               cs3 = cs1; cs3[1] += 1; cs3[2] += 1;
               cs4 = cs1; cs4[2] += 1;
               break;
           case FACE_Y:
               cs2 = cs1; cs2[2] += 1;
               cs3 = cs1; cs3[0] += 1; cs3[2] += 1;
               cs4 = cs1; cs4[0] += 1;
               break;
           case FACE_Z:
               cs2 = cs1; cs2[0] += 1;
               cs3 = cs1; cs3[0] += 1; cs3[1] += 1;
               cs4 = cs1; cs4[1] += 1;
               break;
       }
       
       res[0].init.position = cs1;
       res[0].end.position = cs2;
       res[1].init.position = cs2;
       res[1].end.position = cs3;
       res[2].init.position = cs3;
       res[2].end.position = cs4;
       res[3].init.position = cs4;
       res[3].end.position = cs1;
       
       return res;
   }

   void fillIntervals(const std::vector<interval_t>& intervals, std::vector<edge_t>& edges, std::vector<face_t>& faces) {
       for (size_t i = 0; i < intervals.size(); ++i) {
           fillEdgesFromInterval(edges, intervals[i]);
           fillFaceFromInterval(faces, intervals[i]);
       }
   }

   void fillFullFaces(const std::vector<triangle_t>& tris_on_face, std::vector<face_t>& faces, std::vector<edge_t>& edges) {
       double area = 0.0;
       double ratio = 0.0;
       for (size_t j = 0; j < tris_on_face.size(); ++j) {
           area += getArea(tris_on_face[j]);
       }
       if (std::abs(area - 1.0) < 1e-4) {
           auto cell = tris_on_face[0].getCell();
           int face = tris_on_face[0].getFace();
           addFace(faces, cell, face, ratio);
           for (size_t k = 0; k < tris_on_face.size(); ++k) {
               auto tri_sides = tris_on_face[k].getSides();
               for (size_t s = 0; s < 3; ++s) {
                   if (tri_sides[s].isOnAnyEdge()) {
                       auto cell_s = tri_sides[s].getCell();
                       int edge = tri_sides[s].getEdge();
                       if (isNewEdge(edges, cell_s, edge, ratio)) {
                           addEdge(edges, cell_s, edge, tri_sides[s]);
                       }
                   }
               }
           }
       }
   }

   void fillEdgesFromContour(const std::vector<side_t>& contour, std::vector<edge_t>& edges) {
       for (size_t i = 0; i < contour.size(); ++i) {
           int edge = contour[i].getEdge();
           auto cell = contour[i].getCell();
           if (edge != NOT_ON_EDGE) {
               if (isEdgeFilled(edges, cell, edge)) {
                   fillSmallerRatio(edges, cell, edge, contour[i]);
               } else {
                   addEdge(edges, cell, edge, contour[i]);
               }
           }
       }
   }

   void fillEdges(const std::vector<side_t>& sides, std::vector<edge_t>& edges) {
       for (size_t i = 0; i < sides.size(); ++i) {
           int edge = sides[i].getEdge();
           auto cell = sides[i].getCell();
           if (edge != NOT_ON_EDGE) {
               if (isEdgeFilled(edges, cell, edge)) {
                   reduceEdgeRatio(edges, cell, edge, sides[i]);
               } else {
                   addEdge(edges, cell, edge, sides[i]);
               }
           }
       }
   }

   void fillEdgesFromInterval(std::vector<edge_t>& edges, const interval_t& interval) {
       auto sides = buildSidesFromCellInterval(interval);
       for (size_t i = 0; i < sides.size(); ++i) {
           int edge = sides[i].getEdge();
           if (edge != NOT_ON_EDGE) {
               if (isEdgeFilled(edges, sides[i].getCell(), edge)) {
                   fillSmallerRatio(edges, sides[i].getCell(), edge, sides[i]);
               } else {
                   addEdge(edges, sides[i].getCell(), edge, sides[i]);
               }
           }
       }
   }

   void fillFaceFromInterval(std::vector<face_t>& faces, const interval_t& interval) {
       side_t aux;
       aux.init.position = interval.ini.cell;
       aux.end.position = interval.end.cell;
       double ratio = 0.0;
       addFace(faces, aux.getCell(), aux.getFace(), ratio);
   }

   void fillFaceFromContour(const std::vector<side_t>& contour, std::vector<face_t>& faces) {
       double area;
       int face;
       auto cell = findContourCell(contour);
       face = findContourFace(contour);
       if (!contour.empty()) {
           area = 1.0 - contourArea(contour);
           addFace(faces, cell, face, area);
       }
   }

   std::vector<side_t> findLargestContour(const std::vector<side_t>& sides) {
       std::vector<side_t> res;
       std::vector<side_t> aux_contour;
       std::vector<side_t> aux_side(1);
       double area = 0;
       double contour_area;
       
       for (size_t i = 0; i < sides.size(); ++i) {
           aux_side[0] = sides[i];
           aux_contour = buildSidesContour(aux_side);
           contour_area = contourArea(aux_contour);
           if (contour_area > area) {
               res = aux_contour;
               area = contour_area;
           }
       }
       return res;
   }

   bool isNewEdge(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, double ratio) {
       for (size_t i = 0; i < edges.size(); ++i) {
           if (edges[i].cell == cell && edges[i].ratio == ratio && edges[i].direction == edge) {
               return false;
           }
       }
       return true;
   }

   bool isEdgeFilled(const std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge) {
       for (size_t i = 0; i < edges.size(); ++i) {
           if (edges[i].cell == cell && edges[i].direction == edge) {
               return true;
           }
       }
       return false;
   }

   void reduceEdgeRatio(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, const side_t& side) {
       for (size_t i = 0; i < edges.size(); ++i) {
           if (edges[i].cell == cell && edges[i].direction == edge) {
               if (edges[i].material_coords[0] != std::min(side.init.position[edge], side.end.position[edge]) &&
                   edges[i].material_coords[1] != std::max(side.init.position[edge], side.end.position[edge]) &&
                   edges[i].ratio != 0) {
                   edges[i].ratio -= side.length();
               }
               return;
           }
       }
   }

   void fillSmallerRatio(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, const side_t& side) {
       double new_ratio = 1.0 - side.length();
       for (size_t i = 0; i < edges.size(); ++i) {
           if (edges[i].cell == cell && edges[i].direction == edge) {
               if (new_ratio < edges[i].ratio) {
                   edges[i].ratio = new_ratio;
               }
               return;
           }
       }
   }

   void addEdge(std::vector<edge_t>& edges, const std::array<int, 3>& cell, int edge, const side_t& side) {
       double ratio = 1.0 - side.length();
       std::vector<edge_t> aux(edges.size() + 1);
       for (size_t i = 0; i < edges.size(); ++i) {
           aux[i] = edges[i];
       }
       edge_t new_edge;
       new_edge.cell = cell;
       new_edge.ratio = ratio;
       new_edge.direction = edge;
       new_edge.material_coords[0] = std::min(side.init.position[edge], side.end.position[edge]);
       new_edge.material_coords[1] = std::max(side.init.position[edge], side.end.position[edge]);
       aux[edges.size()] = new_edge;
       edges = aux;
   }

   void addFace(std::vector<face_t>& faces, const std::array<int, 3>& cell, int face, double ratio) {
       std::vector<face_t> aux(faces.size() + 1);
       for (size_t i = 0; i < faces.size(); ++i) {
           aux[i] = faces[i];
       }
       face_t new_face;
       new_face.cell = cell;
       new_face.ratio = ratio;
       new_face.direction = face;
       aux[faces.size()] = new_face;
       faces = aux;
   }

   void addRatio(std::vector<double>& ratios, double ratio) {
       if (ratios.empty()) {
           ratios.push_back(ratio);
       } else {
           ratios.push_back(ratio);
       }
   }

   bool isNewRatio(const std::vector<double>& ratios, double ratio, double tol) {
       for (size_t i = 0; i < ratios.size(); ++i) {
           if (eq_ratio(ratios[i], ratio, tol)) {
               return false;
           }
       }
       return true;
   }

   std::vector<edge_t> filterEdgesByMedia(const std::vector<edge_t>& edges, double ratio) {
       int n = 0;
       for (size_t i = 0; i < edges.size(); ++i) {
           if (eq_ratio(edges[i].ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) {
               n++;
           }
       }
       std::vector<edge_t> res(n);
       n = 0;
       for (size_t i = 0; i < edges.size(); ++i) {
           if (eq_ratio(edges[i].ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) {
               res[n] = edges[i];
               n++;
           }
       }
       return res;
   }

   std::vector<face_t> filterFacesByMedia(const std::vector<face_t>& faces, double ratio) {
       int n = 0;
       for (size_t i = 0; i < faces.size(); ++i) {
           if (eq_ratio(faces[i].ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) {
               n++;
           }
       }
       std::vector<face_t> res(n);
       n = 0;
       for (size_t i = 0; i < faces.size(); ++i) {
           if (eq_ratio(faces[i].ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) {
               res[n] = faces[i];
               n++;
           }
       }
       return res;
   }

} // namespace conformal_m