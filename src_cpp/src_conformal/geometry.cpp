#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream>

// Forward declarations for types defined in conformal_types_m
// Assuming conformal_types_m defines: triangle_t, side_t, coord_t
// And constants: NOT_ON_FACE, FACE_X, FACE_Y, FACE_Z, EDGE_X, EDGE_Y, EDGE_Z

// Placeholder for types from conformal_types_m
// In a real translation, these would be included from a header file
struct coord_t {
    std::vector<double> position; // Assuming 3D coordinates
    // Methods assumed to exist based on usage
    bool isOnVertex() const;
    bool isOnAnyFace() const;
    int getCell() const; // Returns integer array/vector
    int getFace() const;
    int getEdge() const;
    // Comparison operators assumed
    bool operator==(const coord_t& other) const;
};

struct side_t {
    coord_t init;
    coord_t end;
    std::vector<double> normal; // Assuming 3D normal vector

    // Methods assumed to exist based on usage
    std::vector<double> getSides() const; // Returns vector of side_t? Or positions?
    // Based on getArea usage: sides = triangle.getSides(); contourArea(sides);
    // getSides likely returns a vector of side_t or similar structure representing edges.
    // Let's assume it returns std::vector<side_t>
    std::vector<side_t> getSides() const;

    bool isOnAnyFace() const;
    std::vector<int> getCell() const; // Returns integer array/vector
    int getFace() const;
    int getEdge() const;
    bool isEquiv(const side_t& other) const;
    // Assignment operator assumed
    side_t& operator=(const side_t& other);
    // Comparison operators assumed
    bool operator==(const side_t& other) const;
};

struct triangle_t {
    // Methods assumed to exist based on usage
    std::vector<side_t> getSides() const;
    bool isOnFace(int face) const;
};

// Constants assumed from conformal_types_m or similar
const int NOT_ON_FACE = -1;
const int FACE_X = 1;
const int FACE_Y = 2;
const int FACE_Z = 3;
const int EDGE_X = 1;
const int EDGE_Y = 2;
const int EDGE_Z = 3;

// Helper function for cross product
std::vector<double> cross(const std::vector<double>& v1, const std::vector<double>& v2) {
    std::vector<double> res(3);
    res[0] = v1[1]*v2[2] - v1[2]*v2[1];
    res[1] = -(v1[0]*v2[2] - v1[2]*v2[0]);
    res[2] = v1[0]*v2[1] - v1[1]*v2[0];
    return res;
}

// Helper function for absolute value
int abs_int(int x) {
    return std::abs(x);
}

// Helper function for floor
int floor_int(double x) {
    return static_cast<int>(std::floor(x));
}

// Helper function for all() comparison for vectors
bool all_equal(const std::vector<double>& v1, const std::vector<double>& v2) {
    if (v1.size() != v2.size()) return false;
    for (size_t i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i]) return false;
    }
    return true;
}

bool all_equal(const std::vector<int>& v1, const std::vector<int>& v2) {
    if (v1.size() != v2.size()) return false;
    for (size_t i = 0; i < v1.size(); ++i) {
        if (v1[i] != v2[i]) return false;
    }
    return true;
}

namespace geometry_m {

    double getArea(const triangle_t& triangle) {
        std::vector<side_t> sides = triangle.getSides();
        return contourArea(sides);
    }

    std::vector<double> cross_vec(const std::vector<double>& v1, const std::vector<double>& v2) {
        std::vector<double> res(3);
        res[0] = v1[1]*v2[2] - v1[2]*v2[1];
        res[1] = -(v1[0]*v2[2] - v1[2]*v2[0]);
        res[2] = v1[0]*v2[1] - v1[1]*v2[0];
        return res;
    }

    side_t mergeSides(const std::vector<side_t>& sides, int edge) {
        std::vector<side_t> sides_copy = sides;
        double c;
        side_t res;
        for (int i = 0; i < static_cast<int>(sides_copy.size()); ++i) {
            if (sides_copy[i].init.position[edge] > sides_copy[i].end.position[edge]) {
                c = sides_copy[i].init.position[edge];
                sides_copy[i].init.position[edge] = sides_copy[i].end.position[edge];
                sides_copy[i].end.position[edge] = c;
            }
        }
        res = sides_copy[0];
        for (int i = 1; i < static_cast<int>(sides_copy.size()); ++i) {
            if (sides_copy[i].init.position[edge] < res.init.position[edge]) {
                res.init.position[edge] = sides_copy[i].init.position[edge];
            }
            if (sides_copy[i].end.position[edge] > res.end.position[edge]) {
                res.end.position[edge] = sides_copy[i].end.position[edge];
            }
        }
        return res;
    }

    std::vector<int> findContourCell(const std::vector<side_t>& contour) {
        std::vector<int> res(3);
        for (int i = 0; i < static_cast<int>(contour.size()); ++i) {
            if (contour[i].isOnAnyFace()) {
                res = contour[i].getCell();
            }
        }
        return res;
    }

    int findContourFace(const std::vector<side_t>& contour) {
        int res = NOT_ON_FACE;
        for (int i = 0; i < static_cast<int>(contour.size()); ++i) {
            if (contour[i].isOnAnyFace()) {
                res = contour[i].getFace();
            }
        }
        if (res == NOT_ON_FACE) {
            throw std::runtime_error("Contour face could not be identified");
        }
        return res;
    }

    std::vector<side_t> buildCellSideSet(const std::vector<side_t>& sides, const std::vector<side_t>& on_sides) {
        std::vector<side_t> sides_on_face;
        std::vector<side_t> contour;
        std::vector<side_t> aux;
        std::vector<side_t> res;
        int face, edge;
        res.clear(); // Represents allocate(res(0))

        for (face = FACE_X; face <= FACE_Z; ++face) {
            sides_on_face = getSidesOnFace(sides, face);
            contour = buildSidesContour(sides_on_face);
            addNewSides(res, contour);
        }
        for (edge = EDGE_X; edge <= EDGE_Z; ++edge) {
            aux = getSidesOnEdge(sides, edge);
            addNewSides(res, aux);
            aux = getSidesOnEdge(on_sides, edge);
            addNewSides(res, aux);
        }
        return res;
    }

    void addNewSides(std::vector<side_t>& sides, const std::vector<side_t>& new_sides) {
        for (int i = 0; i < static_cast<int>(new_sides.size()); ++i) {
            if (isNewSide(sides, new_sides[i])) {
                addNewSide(sides, new_sides[i]);
            }
        }
    }

    bool isNewSide(const std::vector<side_t>& sides, const side_t& side) {
        bool result = true;
        for (int i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isEquiv(side)) {
                result = false;
            }
        }
        return result;
    }

    void addNewSide(std::vector<side_t>& sides, const side_t& side) {
        if (sides.empty()) {
            sides.clear();
            sides.push_back(side);
        } else {
            std::vector<side_t> aux(sides.size() + 1);
            for (size_t i = 0; i < sides.size(); ++i) {
                aux[i] = sides[i];
            }
            aux[sides.size()] = side;
            sides.clear();
            sides.resize(aux.size());
            sides = aux;
        }
    }

    std::vector<side_t> buildSidesContour(const std::vector<side_t>& sides) {
        std::vector<side_t> inner_path;
        std::vector<side_t> res;
        coord_t init, end;

        if (sides.empty()) {
            res.clear();
        } else {
            inner_path = getPathOnFace(sides);
            init = inner_path[0].init;
            end = inner_path[inner_path.size() - 1].end;
            if (init.isOnVertex() && end.isOnVertex()) {
                res = buildVertexToVertexContour(inner_path);
            } else if (init.isOnVertex() && !end.isOnVertex()) {
                res = buildVertexToSideContour(inner_path);
            } else if (!init.isOnVertex() && end.isOnVertex()) {
                res = buildSideToVertexContour(inner_path);
            } else {
                res = buildSideToSideContour(inner_path);
            }
        }
        return res;
    }

    std::vector<side_t> buildVertexToVertexContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        int mid_corner_idx;
        std::vector<std::vector<double>> corners(3, std::vector<double>(4));

        corners = buildCorners(inner_path[0], inner_path[0].getFace());

        res.resize(inner_path.size() + 2);
        for (size_t i = 0; i < inner_path.size(); ++i) {
            res[i] = inner_path[i];
        }
        mid_corner_idx = cornerIndex(corners, inner_path[inner_path.size() - 1].end.position) % 4 + 1;
        res[inner_path.size() + 1] = buildSide(inner_path[inner_path.size() - 1].end.position, corners[0][mid_corner_idx], corners[1][mid_corner_idx], corners[2][mid_corner_idx]);
        res[inner_path.size() + 2] = buildSide(std::vector<double>{corners[0][mid_corner_idx], corners[1][mid_corner_idx], corners[2][mid_corner_idx]}, inner_path[0].init.position);
        return res;
    }

    std::vector<side_t> buildVertexToSideContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        coord_t init, end;
        int i, idx;
        std::vector<std::vector<double>> corners(3, std::vector<double>(4));
        side_t cell_side;

        init = inner_path[0].init;
        end = inner_path[inner_path.size() - 1].end;
        corners = buildCorners(inner_path[0], inner_path[0].getFace());

        res.resize(inner_path.size());
        for (size_t i = 0; i < inner_path.size(); ++i) {
            res[i] = inner_path[i];
        }
        for (i = 0; i < 4; ++i) {
            cell_side.init.position = std::vector<double>{corners[0][i], corners[1][i], corners[2][i]};
            cell_side.end.position = std::vector<double>{corners[0][(i+1)%4], corners[1][(i+1)%4], corners[2][(i+1)%4]};
            if (all_equal(cell_side.getCell(), std::vector<int>{floor_int(end.position[0]), floor_int(end.position[1]), floor_int(end.position[2])}) &&
                (cell_side.getEdge() == end.getEdge())) {
                idx = i;
                break;
            }
        }
        addSide(res, buildSide(end.position, std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]}));
        while (!all_equal(std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]}, init.position)) {
            addSide(res, buildSide(std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]}, std::vector<double>{corners[0][(idx+2)%4], corners[1][(idx+2)%4], corners[2][(idx+2)%4]}));
            idx = idx + 1;
        }
        return res;
    }

    std::vector<side_t> buildSideToVertexContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        coord_t init, end;
        int idx;
        std::vector<std::vector<double>> corners(3, std::vector<double>(4));
        side_t cell_side;

        init = inner_path[0].init;
        end = inner_path[inner_path.size() - 1].end;
        corners = buildCorners(inner_path[0], inner_path[0].getFace());

        res.resize(inner_path.size());
        for (size_t i = 0; i < inner_path.size(); ++i) {
            res[i] = inner_path[i];
        }

        idx = cornerIndex(corners, end.position);

        cell_side.init.position = std::vector<double>{corners[0][idx], corners[1][idx], corners[2][idx]};
        cell_side.end.position = std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]};
        while (!((all_equal(cell_side.getCell(), std::vector<int>{floor_int(init.position[0]), floor_int(init.position[1]), floor_int(init.position[2])})) &&
                 (cell_side.getEdge() == init.getEdge()))) {
            addSide(res, buildSide(cell_side.init.position, cell_side.end.position));
            cell_side.init.position = std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]};
            cell_side.end.position = std::vector<double>{corners[0][(idx+2)%4], corners[1][(idx+2)%4], corners[2][(idx+2)%4]};
            idx = idx + 1;
        }
        addSide(res, buildSide(cell_side.init.position, init.position));

        return res;
    }

    std::vector<side_t> buildSideToSideContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        side_t cell_side;
        int i, idx_i, idx_e, idx;
        std::vector<std::vector<double>> corners(3, std::vector<double>(4));
        coord_t init, end;

        init = inner_path[0].init;
        end = inner_path[inner_path.size() - 1].end;

        corners = buildCorners(inner_path[0], inner_path[0].getFace());
        res.resize(inner_path.size());
        for (size_t i = 0; i < inner_path.size(); ++i) {
            res[i] = inner_path[i];
        }
        for (i = 0; i < 4; ++i) {
            cell_side.init.position = std::vector<double>{corners[0][i], corners[1][i], corners[2][i]};
            cell_side.end.position = std::vector<double>{corners[0][(i+1)%4], corners[1][(i+1)%4], corners[2][(i+1)%4]};
            if (all_equal(cell_side.getCell(), std::vector<int>{floor_int(init.position[0]), floor_int(init.position[1]), floor_int(init.position[2])}) &&
                (cell_side.getEdge() == init.getEdge())) {
                idx_i = i;
            }
            if (all_equal(cell_side.getCell(), std::vector<int>{floor_int(end.position[0]), floor_int(end.position[1]), floor_int(end.position[2])}) &&
                (cell_side.getEdge() == end.getEdge())) {
                idx_e = i;
            }
        }
        idx = (idx_e % 4) + 1;
        addSide(res, buildSide(end.position, std::vector<double>{corners[0][idx], corners[1][idx], corners[2][idx]}));

        cell_side.init.position = std::vector<double>{corners[0][idx], corners[1][idx], corners[2][idx]};
        cell_side.end.position = std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]};
        while (!((all_equal(cell_side.getCell(), std::vector<int>{floor_int(init.position[0]), floor_int(init.position[1]), floor_int(init.position[2])})) &&
                 (cell_side.getEdge() == init.getEdge()))) {
            addSide(res, buildSide(cell_side.init.position, cell_side.end.position));
            cell_side.init.position = std::vector<double>{corners[0][(idx+1)%4], corners[1][(idx+1)%4], corners[2][(idx+1)%4]};
            cell_side.end.position = std::vector<double>{corners[0][(idx+2)%4], corners[1][(idx+2)%4], corners[2][(idx+2)%4]};
            idx = idx + 1;
        }
        addSide(res, buildSide(cell_side.init.position, init.position));

        return res;
    }

    void addSide(std::vector<side_t>& sides, const side_t& side) {
        std::vector<side_t> aux(sides.size() + 1);
        for (size_t i = 0; i < sides.size(); ++i) {
            aux[i] = sides[i];
        }
        aux[sides.size()] = side;
        sides.clear();
        sides.resize(aux.size());
        sides = aux;
    }

    side_t buildSide(const std::vector<double>& c1, const std::vector<double>& c2) {
        side_t res;
        res.init.position = c1;
        res.end.position = c2;
        return res;
    }

    int cornerIndex(const std::vector<std::vector<double>>& corners, const std::vector<double>& vertex) {
        int result = 0;
        for (int i = 0; i < 4; ++i) {
            if (all_equal(vertex, std::vector<double>{corners[0][i], corners[1][i], corners[2][i]})) {
                result = i + 1;
                break;
            }
        }
        return result;
    }

    std::vector<std::vector<int>> buildCorners(const side_t& side, int face) {
        std::vector<std::vector<int>> res(3, std::vector<int>(4));
        std::vector<int> cell = side.getCell();
        std::vector<int> aux(3);

        if (face == FACE_X) {
            res[0][0] = cell[0]; res[1][0] = cell[1]; res[2][0] = cell[2];
            res[0][1] = cell[0]; res[1][1] = cell[1] + 1; res[2][1] = cell[2];
            res[0][2] = cell[0]; res[1][2] = cell[1] + 1; res[2][2] = cell[2] + 1;
            res[0][3] = cell[0]; res[1][3] = cell[1]; res[2][3] = cell[2] + 1;
        } else if (face == FACE_Y) {
            res[0][0] = cell[0]; res[1][0] = cell[1]; res[2][0] = cell[2];
            res[0][3] = cell[0] + 1; res[1][3] = cell[1]; res[2][3] = cell[2];
            res[0][2] = cell[0] + 1; res[1][2] = cell[1]; res[2][2] = cell[2] + 1;
            res[0][1] = cell[0]; res[1][1] = cell[1]; res[2][1] = cell[2] + 1;
        } else if (face == FACE_Z) {
            res[0][0] = cell[0]; res[1][0] = cell[1]; res[2][0] = cell[2];
            res[0][1] = cell[0] + 1; res[1][1] = cell[1]; res[2][1] = cell[2];
            res[0][2] = cell[0] + 1; res[1][2] = cell[1] + 1; res[2][2] = cell[2];
            res[0][3] = cell[0]; res[1][3] = cell[1] + 1; res[2][3] = cell[2];
        }

        if (isClockwise(side, face)) {
            aux = res[1];
            res[1] = res[3];
            res[3] = aux;
        }
        return res;
    }

    bool isClockwise(const side_t& side, int face) {
        bool result = true;
        std::vector<double> x_prod(3), diff(3);
        diff[0] = side.end.position[0] - side.init.position[0];
        diff[1] = side.end.position[1] - side.init.position[1];
        diff[2] = side.end.position[2] - side.init.position[2];
        x_prod = cross_vec(diff, side.normal);
        if (x_prod[face - 1] < 0) result = false;
        return result;
    }

    double contourArea(const std::vector<side_t>& contour, int orientation) {
        std::vector<side_t> aux_contour;
        double res = 0;
        int face, i, dir1, dir2;

        face = orientation;
        aux_contour = contour;
        if (isClockwise(aux_contour[0], face)) {
            for (i = 0; i < static_cast<int>(aux_contour.size()); ++i) {
                aux_contour[aux_contour.size() - 1 - i].init = aux_contour[i].end;
                aux_contour[aux_contour.size() - 1 - i].end = aux_contour[i].init;
            }
        }

        dir1 = (face % 3) + 1;
        dir2 = ((face + 1) % 3) + 1;
        res = 0;
        for (i = 0; i < static_cast<int>(aux_contour.size()); ++i) {
            res += aux_contour[i].init.position[dir1 - 1] * aux_contour[i].end.position[dir2 - 1] -
                   aux_contour[i].end.position[dir1 - 1] * aux_contour[i].init.position[dir2 - 1];
        }
        res = 0.5 * res;
        return res;
    }

    std::vector<side_t> getPathOnFace(const std::vector<side_t>& sides) {
        std::vector<side_t> res;
        int i, n;
        res.resize(sides.size());
        n = 0;
        while (n < static_cast<int>(res.size())) {
            for (i = 0; i < static_cast<int>(sides.size()); ++i) {
                if (n == 0) {
                    if (!sides[i].init.isOnAnyFace()) {
                        n = n + 1;
                        res[n - 1] = sides[i];
                    }
                } else if (n != 0 && all_equal(sides[i].init.position, res[n - 1].end.position)) {
                    n = n + 1;
                    res[n - 1] = sides[i];
                }
            }
        }
        return res;
    }

    int getCellDistance(const std::vector<int>& ref_cell, const std::vector<int>& cell, int edge) {
        int res = 0;
        int i;
        for (i = 0; i < 3; ++i) {
            res += (i + 1) * abs_int(ref_cell[i] - cell[i]);
        }
        if (edge == EDGE_X) {
            if (res == 2) res = 1;
            if (res == 3) res = 2;
            if (res == 5) res = 3;
        } else if (edge == EDGE_Y) {
            if (res == 3) res = 2;
            if (res == 4) res = 3;
        }
        res = res + 1;
        return res;
    }

    std::vector<triangle_t> getTrianglesOffFaces(const std::vector<triangle_t>& triangles) {
        std::vector<triangle_t> res;
        int i, j, n;
        n = 0;
        for (i = 0; i < static_cast<int>(triangles.size()); ++i) {
            if (!(triangles[i].isOnFace(1) || triangles[i].isOnFace(2) || triangles[i].isOnFace(3))) {
                n = n + 1;
            }
        }
        res.resize(n);
        n = 0;
        for (i = 0; i < static_cast<int>(triangles.size()); ++i) {
            if (!(triangles[i].isOnFace(1) || triangles[i].isOnFace(2) || triangles[i].isOnFace(3))) {
                n = n + 1;
                res[n - 1] = triangles[i]; // Note: Original code had 'triangles(j)' which is likely a typo for 'triangles(i)'
            }
        }
        return res;
    }

    std::vector<triangle_t> getTrianglesOnFace(const std::vector<triangle_t>& tris, int face) {
        std::vector<triangle_t> res;
        int i, n;
        n = 0;
        for (i = 0; i < static_cast<int>(tris.size()); ++i) {
            if (tris[i].isOnFace(face)) {
                n = n + 1;
            }
        }
        res.resize(n);
        n = 0;
        for (i = 0; i < static_cast<int>(tris.size()); ++i) {
            if (tris[i].isOnFace(face)) {
                n = n + 1;
                res[n - 1] = tris[i];
            }
        }
        return res;
    }

    std::vector<side_t> getSidesOnFace(const std::vector<side_t>& sides, int face) {
        std::vector<side_t> res;
        int i, n;
        n = 0;
        for (i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isOnFace(face)) {
                n = n + 1;
            }
        }
        res.resize(n);
        n = 0;
        for (i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isOnFace(face)) {
                n = n + 1;
                res[n - 1] = sides[i];
            }
        }
        return res;
    }

    std::vector<side_t> getSidesOnEdge(const std::vector<side_t>& sides, int edge) {
        std::vector<side_t> res;
        int i, n;
        n = 0;
        for (i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isOnEdge(edge)) {
                n = n + 1;
            }
        }
        res.resize(n);
        n = 0;
        for (i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isOnEdge(edge)) {
                n = n + 1;
                res[n - 1] = sides[i];
            }
        }
        return res;
    }

    std::vector<side_t> getSidesOnAdjacentEdges(const std::vector<side_t>& sides, int face) {
        std::vector<side_t> res;
        int i, n;
        n = 0;
        for (i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isOnEdge((face % 3) + 1) ||
                sides[i].isOnEdge(((face + 1) % 3) + 1)) {
                n = n + 1;
            }
        }
        res.resize(n);
        n = 0;
        for (i = 0; i < static_cast<int>(sides.size()); ++i) {
            if (sides[i].isOnEdge((face % 3) + 1) ||
                sides[i].isOnEdge(((face + 1) % 3) + 1)) {
                n = n + 1;
                res[n - 1] = sides[i];
            }
        }
        return res;
    }

} // namespace geometry_m