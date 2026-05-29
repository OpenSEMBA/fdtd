#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "conformal_m.h"
#include "conformal_types.h"
#include "geometry_m.h"

namespace conformal_m {

namespace {

using conformal_types_m::EDGE_X;
using conformal_types_m::EDGE_Y;
using conformal_types_m::EDGE_Z;
using conformal_types_m::FACE_X;
using conformal_types_m::FACE_Y;
using conformal_types_m::FACE_Z;
using conformal_types_m::NOT_ON_EDGE;
using conformal_types_m::interval_t;
using conformal_types_m::side_t;
using conformal_types_m::triangle_t;

constexpr double EDGE_RATIO_EQ_TOLERANCE = 1e-5;
constexpr double FACE_RATIO_EQ_TOLERANCE = 1e-3;

bool eq_ratio(double a, double b, double tol) {
    return std::abs(a - b) < tol;
}

bool cellsEqual(const int32_t cell[3], const std::vector<int>& other) {
    return cell[0] == other[0] && cell[1] == other[1] && cell[2] == other[2];
}

void copyCell(const std::vector<int>& from, int32_t to[3]) {
    for (int i = 0; i < 3; ++i) {
        to[i] = static_cast<int32_t>(from[i]);
    }
}

bool isNewRatio(const std::vector<double>& ratios, double ratio, double tol) {
    for (double r : ratios) {
        if (eq_ratio(r, ratio, tol)) {
            return false;
        }
    }
    return true;
}

void addRatio(std::vector<double>& ratios, double ratio) {
    ratios.push_back(ratio);
}

std::vector<edge_t> filterEdgesByMedia(const std::vector<edge_t>& edges, double ratio) {
    std::vector<edge_t> res;
    for (const auto& e : edges) {
        if (eq_ratio(e.ratio, ratio, EDGE_RATIO_EQ_TOLERANCE)) {
            res.push_back(e);
        }
    }
    return res;
}

std::vector<face_t> filterFacesByMedia(const std::vector<face_t>& faces, double ratio) {
    std::vector<face_t> res;
    for (const auto& f : faces) {
        if (eq_ratio(f.ratio, ratio, FACE_RATIO_EQ_TOLERANCE)) {
            res.push_back(f);
        }
    }
    return res;
}

bool isNewEdge(const std::vector<edge_t>& edges,
               const std::vector<int>& cell,
               int edge,
               double ratio) {
    for (const auto& e : edges) {
        if (cellsEqual(e.cell, cell) && eq_ratio(e.ratio, ratio, EDGE_RATIO_EQ_TOLERANCE) &&
            e.direction == edge) {
            return false;
        }
    }
    return true;
}

bool isEdgeFilled(const std::vector<edge_t>& edges, const std::vector<int>& cell, int edge) {
    for (const auto& e : edges) {
        if (cellsEqual(e.cell, cell) && e.direction == edge) {
            return true;
        }
    }
    return false;
}

void addEdge(std::vector<edge_t>& edges, const std::vector<int>& cell, int edge, const side_t& side) {
    if (edge < EDGE_X || edge > EDGE_Z) {
        return;
    }
    const double ratio = 1.0 - side.length();
    edge_t new_edge;
    copyCell(cell, new_edge.cell);
    new_edge.ratio = ratio;
    new_edge.direction = edge;
    new_edge.material_coords[0] =
        std::min(side.init.position[static_cast<std::size_t>(edge - 1)],
                 side.end.position[static_cast<std::size_t>(edge - 1)]);
    new_edge.material_coords[1] =
        std::max(side.init.position[static_cast<std::size_t>(edge - 1)],
                 side.end.position[static_cast<std::size_t>(edge - 1)]);
    edges.push_back(new_edge);
}

void addFace(std::vector<face_t>& faces, const std::vector<int>& cell, int face, double ratio) {
    if (face < FACE_X || face > FACE_Z) {
        return;
    }
    face_t new_face;
    copyCell(cell, new_face.cell);
    new_face.ratio = ratio;
    new_face.direction = face;
    faces.push_back(new_face);
}

void reduceEdgeRatio(std::vector<edge_t>& edges,
                     const std::vector<int>& cell,
                     int edge,
                     const side_t& side) {
    const int idx = edge - 1;
    for (auto& e : edges) {
        if (cellsEqual(e.cell, cell) && e.direction == edge) {
            if (e.material_coords[0] !=
                    std::min(side.init.position[static_cast<std::size_t>(idx)],
                             side.end.position[static_cast<std::size_t>(idx)]) &&
                e.material_coords[1] !=
                    std::max(side.init.position[static_cast<std::size_t>(idx)],
                             side.end.position[static_cast<std::size_t>(idx)]) &&
                e.ratio != 0.0) {
                e.ratio -= side.length();
            }
            return;
        }
    }
}

void fillSmallerRatio(std::vector<edge_t>& edges,
                      const std::vector<int>& cell,
                      int edge,
                      const side_t& side) {
    for (auto& e : edges) {
        if (cellsEqual(e.cell, cell) && e.direction == edge) {
            const double new_ratio = 1.0 - side.length();
            if (new_ratio < e.ratio) {
                e.ratio = new_ratio;
            }
            return;
        }
    }
}

void fillEdgesFromContour(const std::vector<side_t>& contour, std::vector<edge_t>& edges) {
    for (const auto& s : contour) {
        const int edge = s.getEdge();
        const std::vector<int> cell = s.getCell();
        if (edge != NOT_ON_EDGE) {
            if (isEdgeFilled(edges, cell, edge)) {
                fillSmallerRatio(edges, cell, edge, s);
            } else {
                addEdge(edges, cell, edge, s);
            }
        }
    }
}

void fillEdges(const std::vector<side_t>& sides, std::vector<edge_t>& edges) {
    for (const auto& s : sides) {
        const int edge = s.getEdge();
        const std::vector<int> cell = s.getCell();
        if (edge != NOT_ON_EDGE) {
            if (isEdgeFilled(edges, cell, edge)) {
                reduceEdgeRatio(edges, cell, edge, s);
            } else {
                addEdge(edges, cell, edge, s);
            }
        }
    }
}

std::vector<side_t> buildSidesFromCellInterval(const interval_t& interval) {
    side_t aux;
    aux.init.position = {static_cast<double>(interval.ini.cell[0]),
                         static_cast<double>(interval.ini.cell[1]),
                         static_cast<double>(interval.ini.cell[2])};
    aux.end.position = {static_cast<double>(interval.end.cell[0]),
                        static_cast<double>(interval.end.cell[1]),
                        static_cast<double>(interval.end.cell[2])};
    const int face = aux.getFace();
    std::vector<int> cs0 = aux.getCell();
    std::vector<std::vector<int>> cs(4, cs0);

    if (face == FACE_X) {
        cs[1][1] += 1;
        cs[2][1] += 1;
        cs[2][2] += 1;
        cs[3][2] += 1;
    } else if (face == FACE_Y) {
        cs[1][2] += 1;
        cs[2][0] += 1;
        cs[2][2] += 1;
        cs[3][0] += 1;
    } else if (face == FACE_Z) {
        cs[1][0] += 1;
        cs[2][0] += 1;
        cs[2][1] += 1;
        cs[3][1] += 1;
    }

    std::vector<side_t> res(4);
    for (int i = 0; i < 4; ++i) {
        res[i].init.position = {static_cast<double>(cs[i][0]),
                                static_cast<double>(cs[i][1]),
                                static_cast<double>(cs[i][2])};
    }
    res[0].end.position = res[1].init.position;
    res[1].end.position = res[2].init.position;
    res[2].end.position = res[3].init.position;
    res[3].end.position = res[0].init.position;
    return res;
}

void fillEdgesFromInterval(std::vector<edge_t>& edges, const interval_t& interval) {
    const std::vector<side_t> sides = buildSidesFromCellInterval(interval);
    for (const auto& s : sides) {
        const int edge = s.getEdge();
        if (edge != NOT_ON_EDGE) {
            const std::vector<int> cell = s.getCell();
            if (isEdgeFilled(edges, cell, edge)) {
                fillSmallerRatio(edges, cell, edge, s);
            } else {
                addEdge(edges, cell, edge, s);
            }
        }
    }
}

void fillFaceFromInterval(std::vector<face_t>& faces, const interval_t& interval) {
    side_t aux;
    aux.init.position = {static_cast<double>(interval.ini.cell[0]),
                         static_cast<double>(interval.ini.cell[1]),
                         static_cast<double>(interval.ini.cell[2])};
    aux.end.position = {static_cast<double>(interval.end.cell[0]),
                        static_cast<double>(interval.end.cell[1]),
                        static_cast<double>(interval.end.cell[2])};
    const double ratio = 0.0;
    addFace(faces, aux.getCell(), aux.getFace(), ratio);
}

void fillFaceFromContour(const std::vector<side_t>& contour, std::vector<face_t>& faces) {
    if (contour.empty()) {
        return;
    }
    const std::vector<int> cell = geometry_m::findContourCell(contour);
    const int face = geometry_m::findContourFace(contour);
    const double area = 1.0 - geometry_m::contourArea(contour);
    addFace(faces, cell, face, area);
}

std::vector<side_t> findLargestContour(const std::vector<side_t>& sides) {
    std::vector<side_t> res;
    double area = 0.0;
    for (const auto& s : sides) {
        const std::vector<side_t> aux_contour = geometry_m::buildSidesContour({s});
        const double contour_area = geometry_m::contourArea(aux_contour);
        if (contour_area > area) {
            res = aux_contour;
            area = contour_area;
        }
    }
    return res;
}

void fillFullFaces(const std::vector<triangle_t>& tris_on_face,
                   std::vector<face_t>& faces,
                   std::vector<edge_t>& edges) {
    double area = 0.0;
    const double ratio = 0.0;
    for (const auto& tri : tris_on_face) {
        area += geometry_m::getArea(tri);
    }
    if (std::abs(area - 1.0) < 1e-4) {
        const std::vector<int> cell = tris_on_face[0].getCell();
        const int face = tris_on_face[0].getFace();
        addFace(faces, cell, face, ratio);
        for (const auto& tri : tris_on_face) {
            const std::vector<side_t> tri_sides = tri.getSides();
            for (const auto& s : tri_sides) {
                if (s.isOnAnyEdge()) {
                    const std::vector<int> scell = s.getCell();
                    const int edge = s.getEdge();
                    if (isNewEdge(edges, scell, edge, ratio)) {
                        addEdge(edges, scell, edge, s);
                    }
                }
            }
        }
    }
}

void fillIntervals(const std::vector<interval_t>& intervals,
                   std::vector<edge_t>& edges,
                   std::vector<face_t>& faces) {
    for (const auto& interval : intervals) {
        fillEdgesFromInterval(edges, interval);
        fillFaceFromInterval(faces, interval);
    }
}

void fillElements(const cell_map_m::cell_map_t& cell_map,
                  std::vector<face_t>& faces,
                  std::vector<edge_t>& edges) {
    faces.clear();
    edges.clear();

    for (const auto& key_entry : cell_map.keys) {
        const std::vector<int> cell = key_entry.cell;
        std::vector<side_t> sides = cell_map.getSidesInCell(cell);
        const std::vector<triangle_t> tris = cell_map.getTrianglesInCell(cell);
        const std::vector<interval_t> intervals = cell_map.getIntervalsInCell(cell);

        for (int face = FACE_X; face <= FACE_Z; ++face) {
            const std::vector<side_t> sides_on_face = geometry_m::getSidesOnFace(sides, face);
            if (!sides_on_face.empty()) {
                const std::vector<side_t> contour = findLargestContour(sides_on_face);
                fillFaceFromContour(contour, faces);
                fillEdgesFromContour(contour, edges);
            }
            const std::vector<triangle_t> tris_on_face = geometry_m::getTrianglesOnFace(tris, face);
            if (!tris_on_face.empty()) {
                fillFullFaces(tris_on_face, faces, edges);
            }
        }
        fillIntervals(intervals, edges, faces);
    }

    for (const auto& key_entry : cell_map.keys) {
        const std::vector<int> cell = key_entry.cell;
        std::vector<side_t> sides = cell_map.getSidesInCell(cell);

        for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
            fillEdges(geometry_m::getSidesOnEdge(sides, edge), edges);
        }

        sides = cell_map.getOnSidesInCell(cell);
        for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
            fillEdges(geometry_m::getSidesOnEdge(sides, edge), edges);
        }
    }
}

void addNewRatios(const std::vector<edge_t>& edges,
                  const std::vector<face_t>& faces,
                  std::vector<double>& edge_ratios,
                  std::vector<double>& face_ratios) {
    edge_ratios.clear();
    face_ratios.clear();
    for (const auto& e : edges) {
        if (isNewRatio(edge_ratios, e.ratio, EDGE_RATIO_EQ_TOLERANCE)) {
            addRatio(edge_ratios, e.ratio);
        }
    }
    for (const auto& f : faces) {
        if (isNewRatio(face_ratios, f.ratio, FACE_RATIO_EQ_TOLERANCE)) {
            addRatio(face_ratios, f.ratio);
        }
    }
}

std::vector<conformal_edge_media_t> addEdgeMedia(const std::vector<edge_t>& edges,
                                                 const std::vector<double>& edge_ratios) {
    std::vector<conformal_edge_media_t> res(edge_ratios.size());
    for (std::size_t i = 0; i < edge_ratios.size(); ++i) {
        res[i].edges = filterEdgesByMedia(edges, edge_ratios[i]);
        res[i].ratio = edge_ratios[i];
        res[i].n_elements = static_cast<int32_t>(res[i].edges.size());
    }
    return res;
}

std::vector<conformal_face_media_t> addFaceMedia(const std::vector<face_t>& faces,
                                                 const std::vector<double>& face_ratios) {
    std::vector<conformal_face_media_t> res(face_ratios.size());
    for (std::size_t i = 0; i < face_ratios.size(); ++i) {
        res[i].faces = filterFacesByMedia(faces, face_ratios[i]);
        res[i].ratio = face_ratios[i];
        res[i].n_elements = static_cast<int32_t>(res[i].faces.size());
    }
    return res;
}

double computeTimeStepScalingFactor(const std::vector<conformal_face_media_t>& faces_media,
                                    const std::vector<conformal_edge_media_t>& edges_media) {
    double res = 1.0;
    cell_map_m::cell_ratios_map_t cell_ratio_map;

    for (const auto& fm : faces_media) {
        for (const auto& f : fm.faces) {
            cell_map_m::cell_t c;
            c.cell = {f.cell[0], f.cell[1], f.cell[2]};
            cell_ratio_map.addFaceRatio(c, f.direction, f.ratio);
        }
    }
    for (const auto& em : edges_media) {
        for (const auto& e : em.edges) {
            cell_map_m::cell_t c;
            c.cell = {e.cell[0], e.cell[1], e.cell[2]};
            cell_ratio_map.addEdgeRatio(c, e.direction, e.ratio);
        }
    }

    for (const auto& key_entry : cell_ratio_map.keys) {
        std::vector<int> cell = key_entry.cell;
        std::vector<int> aux_cell = cell;
        cell_map_m::cell_ratios_t cell_ratio_info = cell_ratio_map.getCellRatiosInCell(cell);

        for (int j = FACE_X; j <= FACE_Z; ++j) {
            const double area = cell_ratio_info.area[static_cast<std::size_t>(j - 1)];
            const int idx1 = j % 3;
            const int idx2 = (j + 1) % 3;
            double l_ratio =
                std::max(cell_ratio_info.length[static_cast<std::size_t>(idx1)],
                         cell_ratio_info.length[static_cast<std::size_t>(idx2)]);

            aux_cell[idx1] = aux_cell[idx1] + 1;
            if (cell_ratio_map.hasKey(aux_cell)) {
                cell_ratio_info = cell_ratio_map.getCellRatiosInCell(aux_cell);
                l_ratio = std::max(l_ratio, cell_ratio_info.length[static_cast<std::size_t>(idx2)]);
            } else {
                l_ratio = 1.0;
            }

            aux_cell = cell;
            aux_cell[idx2] = aux_cell[idx2] + 1;
            if (cell_ratio_map.hasKey(aux_cell)) {
                cell_ratio_info = cell_ratio_map.getCellRatiosInCell(aux_cell);
                l_ratio = std::max(l_ratio, cell_ratio_info.length[static_cast<std::size_t>(idx1)]);
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

ConformalMedia_t buildConformalVolume(const cell_map_m::ConformalPECElements_t& volume) {
    ConformalMedia_t res;
    cell_map_m::cell_map_t cell_map;
    std::vector<double> edge_ratios;
    std::vector<double> face_ratios;
    std::vector<edge_t> edges;
    std::vector<face_t> faces;

    cell_map_m::buildCellMap(cell_map, volume);
    fillElements(cell_map, faces, edges);
    addNewRatios(edges, faces, edge_ratios, face_ratios);

    res.edge_media = addEdgeMedia(edges, edge_ratios);
    res.face_media = addFaceMedia(faces, face_ratios);
    res.n_edges_media = static_cast<int32_t>(res.edge_media.size());
    res.n_faces_media = static_cast<int32_t>(res.face_media.size());
    res.time_step_scale_factor =
        computeTimeStepScalingFactor(res.face_media, res.edge_media);
    res.tag = volume.tag;
    return res;
}

} // namespace

std::vector<ConformalMedia_t> buildConformalMedia(const ConformalPECRegions_t& regions) {
    std::vector<ConformalMedia_t> res(regions.volumes.size());
    for (std::size_t i = 0; i < regions.volumes.size(); ++i) {
        res[i] = buildConformalVolume(regions.volumes[i]);
    }
    return res;
}

} // namespace conformal_m
