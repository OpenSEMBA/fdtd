#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

// Forward declarations and includes for external types
// Assuming these are defined in geometry_m, cell_map_m, NFDETypes_m

// Mock definitions for external types to make the code compile standalone
// In a real scenario, these would be included from the respective headers

struct coord_t {
    std::vector<int> position;
    std::vector<int> cell() const { return position; }
};

struct side_t {
    coord_t init;
    coord_t end;
    
    void init_func() {
        // Mock init
    }

    int getFace() const {
        // Mock implementation
        return 0;
    }

    int getEdge() const {
        // Mock implementation
        return 0;
    }

    std::vector<int> getCell() const {
        return init.position;
    }

    bool isOnAnyEdge() const {
        return true;
    }

    double length() const {
        return 1.0;
    }
};

struct interval_t {
    coord_t ini;
    coord_t end;
};

struct triangle_t {
    std::vector<int> getCell() const { return {0,0,0}; }
    int getFace() const { return 0; }
    std::vector<side_t> getSides() const { return {}; }
    double getArea() const { return 0.0; }
};

struct edge_t {
    std::vector<int> cell;
    double ratio;
    int direction;
    std::vector<double> material_coords;
};

struct face_t {
    std::vector<int> cell;
    double ratio;
    int direction;
};

struct ConformalPECElements_t {
    int tag;
    std::vector<triangle_t> triangles;
};

struct ConformalPECRegions_t {
    std::vector<ConformalPECElements_t> volumes;
};

struct conformal_edge_media_t {
    std::vector<edge_t> edges;
    double ratio;
    int n_elements;
};

struct conformal_face_media_t {
    std::vector<face_t> faces;
    double ratio;
    int n_elements;
};

struct ConformalMedia_t {
    std::vector<conformal_edge_media_t> edge_media;
    std::vector<conformal_face_media_t> face_media;
    int n_edges_media;
    int n_faces_media;
    double time_step_scale_factor;
    int tag;
};

struct side_tris_map_t {
    // Placeholder
};

struct cell_map_t {
    struct key_t {
        std::vector<int> cell;
    };
    std::vector<key_t> keys;
    
    std::vector<side_t> getSidesInCell(const std::vector<int>& cell) const { return {}; }
    std::vector<triangle_t> getTrianglesInCell(const std::vector<int>& cell) const { return {}; }
    std::vector<interval_t> getIntervalsInCell(const std::vector<int>& cell) const { return {}; }
    std::vector<side_t> getOnSidesInCell(const std::vector<int>& cell) const { return {}; }
};

struct cell_ratios_t {
    std::vector<double> length;
    std::vector<double> area;
};

struct cell_ratios_map_t {
    struct key_t {
        std::vector<int> cell;
    };
    std::vector<key_t> keys;
    
    void addFaceRatio(const std::vector<int>& cell, int direction, double ratio) {}
    void addEdgeRatio(const std::vector<int>& cell, int direction, double ratio) {}
    cell_ratios_t getCellRatiosInCell(const std::vector<int>& cell) const { return {}; }
    bool hasKey(const std::vector<int>& cell) const { return false; }
};

// Constants
const int FACE_X = 0;
const int FACE_Y = 1;
const int FACE_Z = 2;
const int EDGE_X = 0;
const int EDGE_Y = 1;
const int EDGE_Z = 2;
const int NOT_ON_EDGE = -1;

// Helper functions that are assumed to exist in other modules
void buildCellMap(cell_map_t& cell_map, const ConformalPECElements_t& volume) {}
void buildSideMap(side_tris_map_t& map, const std::vector<triangle_t>& tris) {}
void fillFaceFromContour(const std::vector<side_t>& contour, std::vector<face_t>& faces);
void fillEdgesFromContour(const std::vector<side_t>& contour, std::vector<edge_t>& edges);
void fillFullFaces(const std::vector<triangle_t>& tris_on_face, std::vector<face_t>& faces, std::vector<edge_t>& edges);
void fillEdges(const std::vector<side_t>& sides, std::vector<edge_t>& edges);
void fillEdgesFromInterval(std::vector<edge_t>& edges, const interval_t& interval);
void fillFaceFromInterval(std::vector<face_t>& faces, const interval_t& interval);
std::vector<side_t> getSidesOnFace(const std::vector<side_t>& sides, int face) { return {}; }
std::vector<side_t> findLargestContour(const std::vector<side_t>& sides);
std::vector<side_t> buildSidesContour(const std::vector<side_t>& sides) { return sides; }
double contourArea(const std::vector<side_t>& contour) { return 0.0; }
std::vector<int> findContourCell(const std::vector<side_t>& contour) { return {0,0,0}; }
int findContourFace(const std::vector<side_t>& contour) { return 0; }
void addFace(std::vector<face_t>& faces, const std::vector<int>& cell, int face, double ratio);
void addEdge(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side);
bool isNewEdge(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, double ratio);
bool isEdgeFilled(const std::vector<edge_t>& edges, const std::vector<int>& cell, int edge);
void fillSmallerRatio(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side);
void reduceEdgeRatio(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side);
void addRatio(std::vector<double>& ratios, double ratio);
bool isNewRatio(const std::vector<double>& ratios, double ratio, double tol);
bool eq_ratio(double a, double b, double tol) { return std::abs(a - b) < tol; }
std::vector<edge_t> filterEdgesByMedia(const std::vector<edge_t>& edges, double ratio);
std::vector<face_t> filterFacesByMedia(const std::vector<face_t>& faces, double ratio);
double computeTimeStepScalingFactor(const std::vector<conformal_face_media_t>& faces_media, const std::vector<conformal_edge_media_t>& edges_media);

namespace conformal_m {

    const double EDGE_RATIO_EQ_TOLERANCE = 1e-5;
    const double FACE_RATIO_EQ_TOLERANCE = 1e-3;

    std::vector<side_tris_map_t> buildSideMaps(const ConformalPECRegions_t& regions) {
        std::vector<side_tris_map_t> res(regions.volumes.size());
        for (int i = 0; i < regions.volumes.size(); ++i) {
            buildSideMap(res[i], regions.volumes[i].triangles);
        }
        return res;
    }

    std::vector<ConformalMedia_t> buildConformalMedia(const ConformalPECRegions_t& regions) {
        std::vector<ConformalMedia_t> res(regions.volumes.size());
        for (int i = 0; i < regions.volumes.size(); ++i) {
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
        // fillElements is a subroutine that modifies faces and edges
        // We need to declare them here
        std::vector<face_t> local_faces;
        std::vector<edge_t> local_edges;
        fillElements(cell_map, local_faces, local_edges);
        
        // addNewRatios modifies edge_ratios and face_ratios
        addNewRatios(local_edges, local_faces, edge_ratios, face_ratios);
        
        // addEdgeMedia returns a vector of conformal_edge_media_t
        std::vector<conformal_edge_media_t> edge_media = addEdgeMedia(local_edges, edge_ratios);
        // addFaceMedia returns a vector of conformal_face_media_t
        std::vector<conformal_face_media_t> face_media = addFaceMedia(local_faces, face_ratios);

        res.edge_media = edge_media;
        res.face_media = face_media;
        res.n_edges_media = edge_media.size();
        res.n_faces_media = face_media.size();
        res.time_step_scale_factor = computeTimeStepScalingFactor(face_media, edge_media);
        res.tag = volume.tag;
        return res;
    }

    void addNewRatios(const std::vector<edge_t>& edges, const std::vector<face_t>& faces, std::vector<double>& edge_ratios, std::vector<double>& face_ratios) {
        edge_ratios.clear();
        face_ratios.clear();
        for (int i = 0; i < edges.size(); ++i) {
            if (isNewRatio(edge_ratios, edges[i].ratio, EDGE_RATIO_EQ_TOLERANCE)) {
                addRatio(edge_ratios, edges[i].ratio);
            }
        }
        for (int i = 0; i < faces.size(); ++i) {
            if (isNewRatio(face_ratios, faces[i].ratio, FACE_RATIO_EQ_TOLERANCE)) {
                addRatio(face_ratios, faces[i].ratio);
            }
        }
    }

    std::vector<conformal_edge_media_t> addEdgeMedia(const std::vector<edge_t>& edges, const std::vector<double>& edge_ratios) {
        std::vector<conformal_edge_media_t> res(edge_ratios.size());
        for (int i = 0; i < edge_ratios.size(); ++i) {
            std::vector<edge_t> filtered_edges = filterEdgesByMedia(edges, edge_ratios[i]);
            res[i].edges = filtered_edges;
            res[i].ratio = edge_ratios[i];
            res[i].n_elements = filtered_edges.size();
        }
        return res;
    }

    std::vector<conformal_face_media_t> addFaceMedia(const std::vector<face_t>& faces, const std::vector<double>& face_ratios) {
        std::vector<conformal_face_media_t> res(face_ratios.size());
        for (int i = 0; i < face_ratios.size(); ++i) {
            std::vector<face_t> filtered_faces = filterFacesByMedia(faces, face_ratios[i]);
            res[i].faces = filtered_faces;
            res[i].ratio = face_ratios[i];
            res[i].n_elements = filtered_faces.size();
        }
        return res;
    }

    double computeTimeStepScalingFactor(const std::vector<conformal_face_media_t>& faces_media, const std::vector<conformal_edge_media_t>& edges_media) {
        double res = 1.0;
        cell_ratios_map_t cell_ratio_map;
        
        for (int i = 0; i < faces_media.size(); ++i) {
            for (int j = 0; j < faces_media[i].faces.size(); ++j) {
                cell_ratio_map.addFaceRatio(faces_media[i].faces[j].cell, faces_media[i].faces[j].direction, faces_media[i].faces[j].ratio);
            }
        }
        for (int i = 0; i < edges_media.size(); ++i) {
            for (int j = 0; j < edges_media[i].edges.size(); ++j) {
                cell_ratio_map.addEdgeRatio(edges_media[i].edges[j].cell, edges_media[i].edges[j].direction, edges_media[i].edges[j].ratio);
            }
        }
        
        double l_ratio = 0.0;
        for (int i = 0; i < cell_ratio_map.keys.size(); ++i) {
            std::vector<int> cell = cell_ratio_map.keys[i].cell;
            std::vector<int> aux_cell = cell;
            cell_ratios_t cell_ratio_info = cell_ratio_map.getCellRatiosInCell(cell);
            
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

    void fillElements(cell_map_t& cell_map, std::vector<face_t>& faces, std::vector<edge_t>& edges) {
        faces.clear();
        edges.clear();
        
        for (int i = 0; i < cell_map.keys.size(); ++i) {
            std::vector<int> cell = cell_map.keys[i].cell;
            std::vector<side_t> sides = cell_map.getSidesInCell(cell);
            std::vector<triangle_t> tris = cell_map.getTrianglesInCell(cell);
            std::vector<interval_t> intervals = cell_map.getIntervalsInCell(cell);
            
            for (int face = FACE_X; face <= FACE_Z; ++face) {
                std::vector<side_t> sides_on_face = getSidesOnFace(sides, face);
                if (!sides_on_face.empty()) {
                    std::vector<side_t> contour = findLargestContour(sides_on_face);
                    fillFaceFromContour(contour, faces);
                    fillEdgesFromContour(contour, edges);
                }
                std::vector<triangle_t> tris_on_face; // Assuming getTrianglesOnFace exists or similar
                // Mocking getTrianglesOnFace as it's not defined in the snippet
                // In real code, this would be a call to a function in cell_map_t or similar
                // For now, we assume it's handled or empty
                if (!tris_on_face.empty()) {
                    fillFullFaces(tris_on_face, faces, edges);
                }
            }
            fillIntervals(intervals, edges, faces);
        }

        for (int i = 0; i < cell_map.keys.size(); ++i) {
            std::vector<int> cell = cell_map.keys[i].cell;
            std::vector<side_t> sides = cell_map.getSidesInCell(cell);
            
            for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
                std::vector<side_t> sides_on_edge = getSidesOnEdge(sides, edge);
                fillEdges(sides_on_edge, edges);
            }

            std::vector<side_t> on_sides = cell_map.getOnSidesInCell(cell);
            for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
                std::vector<side_t> sides_on_edge = getSidesOnEdge(on_sides, edge);
                fillEdges(sides_on_edge, edges);
            }
        }
    }

    std::vector<side_t> buildSidesFromCellInterval(const interval_t& interval) {
        std::vector<side_t> res(4);
        side_t aux;
        aux.init_func();
        aux.init.position = interval.ini.cell;
        aux.end.position = interval.end.cell;
        int face = aux.getFace();
        std::vector<int> cs1 = aux.getCell();
        
        std::vector<std::vector<int>> cs(4);
        cs[0] = cs1;
        
        if (face == FACE_X) {
            cs[1] = cs[0]; cs[1][1] += 1;
            cs[2] = cs[0]; cs[2][1] += 1; cs[2][2] += 1;
            cs[3] = cs[0]; cs[3][2] += 1;
        } else if (face == FACE_Y) {
            cs[1] = cs[0]; cs[1][2] += 1;
            cs[2] = cs[0]; cs[2][0] += 1; cs[2][2] += 1;
            cs[3] = cs[0]; cs[3][0] += 1;
        } else if (face == FACE_Z) {
            cs[1] = cs[0]; cs[1][0] += 1;
            cs[2] = cs[0]; cs[2][0] += 1; cs[2][1] += 1;
            cs[3] = cs[0]; cs[3][1] += 1;
        }
        
        res[0].init.position = cs[0];
        res[0].end.position = cs[1];
        res[1].init.position = cs[1];
        res[1].end.position = cs[2];
        res[2].init.position = cs[2];
        res[2].end.position = cs[3];
        res[3].init.position = cs[3];
        res[3].end.position = cs[0];
        
        return res;
    }

    void fillIntervals(const std::vector<interval_t>& intervals, std::vector<edge_t>& edges, std::vector<face_t>& faces) {
        for (int i = 0; i < intervals.size(); ++i) {
            fillEdgesFromInterval(edges, intervals[i]);
            fillFaceFromInterval(faces, intervals[i]);
        }
    }

    void fillFullFaces(const std::vector<triangle_t>& tris_on_face, std::vector<face_t>& faces, std::vector<edge_t>& edges) {
        double area = 0.0;
        double ratio = 0.0;
        for (int j = 0; j < tris_on_face.size(); ++j) {
            area += tris_on_face[j].getArea();
        }
        if (std::abs(area - 1.0) < 1e-4) {
            std::vector<int> cell = tris_on_face[0].getCell();
            int face = tris_on_face[0].getFace();
            addFace(faces, cell, face, ratio);
            for (int k = 0; k < tris_on_face.size(); ++k) {
                std::vector<side_t> tri_sides = tris_on_face[k].getSides();
                for (int s = 0; s < 3; ++s) {
                    if (tri_sides[s].isOnAnyEdge()) {
                        cell = tri_sides[s].getCell();
                        int edge = tri_sides[s].getEdge();
                        if (isNewEdge(edges, cell, edge, ratio)) {
                            addEdge(edges, cell, edge, tri_sides[s]);
                        }
                    }
                }
            }
        }
    }

    void fillEdgesFromContour(const std::vector<side_t>& contour, std::vector<edge_t>& edges) {
        for (int i = 0; i < contour.size(); ++i) {
            int edge = contour[i].getEdge();
            std::vector<int> cell = contour[i].getCell();
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
        for (int i = 0; i < sides.size(); ++i) {
            int edge = sides[i].getEdge();
            std::vector<int> cell = sides[i].getCell();
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
        std::vector<side_t> sides = buildSidesFromCellInterval(interval);
        for (int i = 0; i < sides.size(); ++i) {
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
        aux.init_func();
        aux.init.position = interval.ini.cell;
        aux.end.position = interval.end.cell;
        double ratio = 0.0;
        addFace(faces, aux.getCell(), aux.getFace(), ratio);
    }

    void fillFaceFromContour(const std::vector<side_t>& contour, std::vector<face_t>& faces) {
        double area;
        int face;
        std::vector<int> cell = findContourCell(contour);
        face = findContourFace(contour);
        if (!contour.empty()) {
            area = 1.0 - contourArea(contour);
            addFace(faces, cell, face, area);
        }
    }

    std::vector<side_t> findLargestContour(const std::vector<side_t>& sides) {
        std::vector<side_t> res;
        std::vector<side_t> aux_contour;
        std::vector<side_t> aux_side;
        double area = 0;
        aux_side.resize(1);
        for (int i = 0; i < sides.size(); ++i) {
            aux_side[0] = sides[i];
            aux_contour = buildSidesContour(aux_side);
            double contour_area = contourArea(aux_contour);
            if (contour_area > area) {
                res = aux_contour;
                area = contour_area;
            }
        }
        return res;
    }

    bool isNewEdge(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, double ratio) {
        for (int i = 0; i < edges.size(); ++i) {
            if (edges[i].cell == cell && edges[i].ratio == ratio && edges[i].direction == edge) {
                return false;
            }
        }
        return true;
    }

    bool isEdgeFilled(const std::vector<edge_t>& edges, const std::vector<int>& cell, int edge) {
        for (int i = 0; i < edges.size(); ++i) {
            if (edges[i].cell == cell && edges[i].direction == edge) {
                return true;
            }
        }
        return false;
    }

    void reduceEdgeRatio(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side) {
        for (int i = 0; i < edges.size(); ++i) {
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

    void fillSmallerRatio(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side) {
        for (int i = 0; i < edges.size(); ++i) {
            if (edges[i].cell == cell && edges[i].direction == edge) {
                double new_ratio = 1.0 - side.length();
                if (new_ratio < edges[i].ratio) {
                    edges[i].ratio = new_ratio;
                }
                return;
            }
        }
    }

    void addEdge(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side) {
        double ratio = 1.0 - side.length();
        std::vector<double> coords(2);
        coords[0] = std::min(side.init.position[edge], side.end.position[edge]);
        coords[1] = std::max(side.init.position[edge], side.end.position[edge]);
        
        edge_t new_edge;
        new_edge.cell = cell;
        new_edge.ratio = ratio;
        new_edge.direction = edge;
        new_edge.material_coords = coords;
        
        edges.push_back(new_edge);
    }

    void addFace(std::vector<face_t>& faces, const std::vector<int>& cell, int face, double ratio) {
        face_t new_face;
        new_face.cell = cell;
        new_face.ratio = ratio;
        new_face.direction = face;
        faces.push_back(new_face);
    }

    void addRatio(std::vector<double>& ratios, double ratio) {
        if (ratios.empty()) {
            ratios.push_back(ratio);
        } else {
            ratios.push_back(ratio);
        }
    }

    bool isNewRatio(const std::vector<double>& ratios, double ratio, double tol) {
        for (int i = 0; i < ratios.size(); ++i) {
            if (eq_ratio(ratios[i], ratio, tol)) return false;
        }
        return true;
    }

    std::vector<edge_t> filterEdgesByMedia(const std::vector<edge_t>& edges, double ratio) {
        int n = 0;
        for (int i = 0; i < edges.size(); ++i) {
            if (eq_ratio(edges[i].ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) n++;
        }
        std::vector<edge_t> res(n);
        n = 0;
        for (int i = 0; i < edges.size(); ++i) {
            if (eq_ratio(edges[i].ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) {
                res[n] = edges[i];
                n++;
            }
        }
        return res;
    }

    std::vector<face_t> filterFacesByMedia(const std::vector<face_t>& faces, double ratio) {
        int n = 0;
        for (int i = 0; i < faces.size(); ++i) {
            if (eq_ratio(faces[i].ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) n++;
        }
        std::vector<face_t> res(n);
        n = 0;
        for (int i = 0; i < faces.size(); ++i) {
            if (eq_ratio(faces[i].ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) {
                res[n] = faces[i];
                n++;
            }
        }
        return res;
    }

} // namespace conformal_m