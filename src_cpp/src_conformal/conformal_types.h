#ifndef CONFORMAL_TYPES_H
#define CONFORMAL_TYPES_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>

namespace conformal_types_m {

    constexpr int FACE_X = 1;
    constexpr int FACE_Y = 2;
    constexpr int FACE_Z = 3;
    constexpr int NOT_ON_FACE = -1;

    constexpr int EDGE_X = 1;
    constexpr int EDGE_Y = 2;
    constexpr int EDGE_Z = 3;
    constexpr int NOT_ON_EDGE = -1;

    constexpr double POS_TOL = 1e-4;

    struct coord_t {
        std::vector<double> position;
        int id = -1;

        coord_t() : position(3, 0.0), id(-1) {}
        coord_t(const std::vector<double>& pos, int i = -1) : position(pos), id(i) {}

        bool isOnVertex() const;
        bool isOnEdge(int edge) const;
        bool isOnAnyEdge() const;
        bool isOnFace(int face) const;
        bool isOnAnyFace() const;
        int getEdge() const;
    };

    struct side_t {
        coord_t init;
        coord_t end;
        std::vector<double> normal;

        side_t() : normal(3, 0.0) {}

        int getEdge() const;
        std::vector<int> getCell() const;
        std::vector<side_t> getSides() const; // Note: Original Fortran 'getSides' is a method of triangle_t, but listed here in side_t? 
                                              // Looking at Fortran: 'procedure :: getSides' is inside triangle_t. 
                                              // However, the prompt asks to preserve names. 
                                              // In Fortran, 'getSides' is a method of triangle_t. 
                                              // I will place it in triangle_t as per Fortran structure.
        bool isOnEdge(int edge) const;
        bool isInCell(const std::vector<int>& cell) const;
        bool isOnAnyEdge() const;
        bool isOnFace(int face) const;
        bool isOnAnyFace() const;
        double length() const;
        bool isEquiv(const side_t& side) const;
        int getFace() const;
    };

    struct side_list_t {
        std::vector<side_t> sides;
    };

    struct triangle_t {
        std::vector<coord_t> vertices;

        triangle_t() : vertices(3) {}

        std::vector<double> getNormal() const;
        int getFace() const;
        std::vector<side_t> getSides() const;
        std::vector<int> getCell() const;
        bool isOnFace(int face) const;
        bool isOnAnyFace() const;
        
    private:
        std::vector<double> centroid() const;
    };

    struct point_t {
        std::vector<int> cell;

        point_t() : cell(3, 0) {}
    };

    struct interval_t {
        point_t ini;
        point_t end;
    };

    // Implementation of coord_t methods

    bool coord_t::isOnVertex() const {
        std::vector<double> delta(3);
        for (int i = 0; i < 3; ++i) {
            delta[i] = position[i] - std::floor(position[i]);
        }
        return (delta[0] == 0.0) && (delta[1] == 0.0) && (delta[2] == 0.0);
    }

    bool coord_t::isOnEdge(int edge) const {
        std::vector<double> delta(3);
        for (int i = 0; i < 3; ++i) {
            delta[i] = position[i] - std::floor(position[i]);
        }
        // Fortran: delta(edge) /= 0 .and. delta(mod(edge,3) + 1) == 0 .and. delta(mod(edge+1,3) + 1) == 0
        // Note: Fortran arrays are 1-based. C++ vectors are 0-based.
        // edge is 1, 2, or 3.
        // delta(edge) -> delta[edge-1]
        // mod(edge,3) + 1: 
        // if edge=1, mod(1,3)=1, +1=2 -> index 1
        // if edge=2, mod(2,3)=2, +1=3 -> index 2
        // if edge=3, mod(3,3)=0, +1=1 -> index 0
        // mod(edge+1,3) + 1:
        // if edge=1, mod(2,3)=2, +1=3 -> index 2
        // if edge=2, mod(3,3)=0, +1=1 -> index 0
        // if edge=3, mod(4,3)=1, +1=2 -> index 1
        
        int idx1 = edge - 1;
        int idx2 = (edge % 3) + 1 - 1; // Convert 1-based result of mod+1 to 0-based index
        int idx3 = ((edge + 1) % 3) + 1 - 1;

        // Fortran logic: delta(edge) > POS_TOL AND delta(idx2) < POS_TOL AND delta(idx3) < POS_TOL
        return (delta[idx1] > POS_TOL) && (delta[idx2] < POS_TOL) && (delta[idx3] < POS_TOL);
    }

    bool coord_t::isOnAnyEdge() const {
        return isOnEdge(EDGE_X) || isOnEdge(EDGE_Y) || isOnEdge(EDGE_Z);
    }

    bool coord_t::isOnFace(int face) const {
        std::vector<double> delta(3);
        for (int i = 0; i < 3; ++i) {
            delta[i] = position[i] - std::floor(position[i]);
        }
        // Fortran: delta(face) < POS_TOL .and. delta(mod(face,3) + 1) > POS_TOL .and. delta(mod(face+1,3) + 1) > POS_TOL
        int idx1 = face - 1;
        int idx2 = (face % 3) + 1 - 1;
        int idx3 = ((face + 1) % 3) + 1 - 1;

        return (delta[idx1] < POS_TOL) && (delta[idx2] > POS_TOL) && (delta[idx3] > POS_TOL);
    }

    bool coord_t::isOnAnyFace() const {
        return isOnFace(FACE_X) || isOnFace(FACE_Y) || isOnFace(FACE_Z);
    }

    int coord_t::getEdge() const {
        int res = NOT_ON_EDGE;
        if (!isOnVertex()) {
            for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
                if (isOnEdge(edge)) {
                    res = edge;
                }
            }
        }
        return res;
    }

    // Implementation of side_t methods

    int side_t::getEdge() const {
        int res = NOT_ON_EDGE;
        for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
            if (isOnEdge(edge)) {
                res = edge;
            }
        }
        return res;
    }

    std::vector<int> side_t::getCell() const {
        std::vector<double> c(3);
        for (int i = 0; i < 3; ++i) {
            c[i] = 0.5 * (init.position[i] + end.position[i]);
        }
        std::vector<int> res(3);
        for (int i = 0; i < 3; ++i) {
            res[i] = static_cast<int>(std::floor(c[i]));
        }
        return res;
    }

    bool side_t::isInCell(const std::vector<int>& cell) const {
        std::vector<int> myCell = getCell();
        return (myCell[0] == cell[0]) && (myCell[1] == cell[1]) && (myCell[2] == cell[2]);
    }

    bool side_t::isOnEdge(int edge) const {
        coord_t c(std::vector<double>(3));
        for (int i = 0; i < 3; ++i) {
            c.position[i] = 0.5 * (end.position[i] + init.position[i]);
        }
        return c.isOnEdge(edge);
    }

    int side_t::getFace() const {
        int res = NOT_ON_FACE;
        for (int face = FACE_X; face <= FACE_Z; ++face) {
            if (isOnFace(face)) {
                res = face;
            }
        }
        return res;
    }

    bool side_t::isOnFace(int face) const {
        coord_t mean;
        for (int i = 0; i < 3; ++i) {
            mean.position[i] = 0.5 * (init.position[i] + end.position[i]);
        }
        return mean.getEdge() == NOT_ON_EDGE && mean.isOnFace(face);
    }

    bool side_t::isOnAnyFace() const {
        return isOnFace(FACE_X) || isOnFace(FACE_Y) || isOnFace(FACE_Z);
    }

    bool side_t::isOnAnyEdge() const {
        return isOnEdge(EDGE_X) || isOnEdge(EDGE_Y) || isOnEdge(EDGE_Z);
    }

    double side_t::length() const {
        double sum = 0.0;
        for (int i = 0; i < 3; ++i) {
            double diff = init.position[i] - end.position[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
    }

    bool side_t::isEquiv(const side_t& side) const {
        bool eq = true;
        bool eq_inv = true;
        for (int i = 0; i < 3; ++i) {
            eq = eq && (std::abs(init.position[i] - side.init.position[i]) < 0.01) && 
                       (std::abs(end.position[i] - side.end.position[i]) < 0.01);
            eq_inv = eq_inv && (std::abs(init.position[i] - side.end.position[i]) < 0.01) && 
                           (std::abs(end.position[i] - side.init.position[i]) < 0.01);
        }
        return eq || eq_inv;
    }

    // Implementation of triangle_t methods

    std::vector<double> triangle_t::getNormal() const {
        std::vector<double> v1(3);
        std::vector<double> v2(3);
        for (int i = 0; i < 3; ++i) {
            v1[i] = vertices[1].position[i] - vertices[0].position[i];
            v2[i] = vertices[2].position[i] - vertices[1].position[i];
        }
        std::vector<double> res(3);
        res[0] = v1[1]*v2[2] - v1[2]*v2[1];
        res[1] = -(v1[0]*v2[2] - v1[2]*v2[0]);
        res[2] = v1[0]*v2[1] - v1[1]*v2[0];
        
        double n = std::sqrt(res[0]*res[0] + res[1]*res[1] + res[2]*res[2]);
        if (n > 0) {
            res[0] /= n;
            res[1] /= n;
            res[2] /= n;
        }
        return res;
    }

    int triangle_t::getFace() const {
        int res = NOT_ON_FACE;
        for (int face = FACE_X; face <= FACE_Z; ++face) {
            if (isOnFace(face)) {
                res = face;
            }
        }
        return res;
    }

    std::vector<side_t> triangle_t::getSides() const {
        std::vector<side_t> res(3);
        for (int i = 0; i < 3; ++i) {
            res[i].init.position = vertices[i].position;
            int nextIdx = (i % 3) + 1; // 0->1, 1->2, 2->0
            res[i].end.position = vertices[nextIdx].position;
            res[i].init.id = vertices[i].id;
            res[i].end.id = vertices[nextIdx].id;
            res[i].normal = getNormal();
        }
        return res;
    }

    std::vector<int> triangle_t::getCell() const {
        std::vector<int> res(3);
        // Fortran: min(a,b,c)
        res[0] = std::min({vertices[0].position[0], vertices[1].position[0], vertices[2].position[0]});
        res[1] = std::min({vertices[0].position[1], vertices[1].position[1], vertices[2].position[1]});
        res[2] = std::min({vertices[0].position[2], vertices[1].position[2], vertices[2].position[2]});
        
        // Apply floor
        for (int i = 0; i < 3; ++i) {
            res[i] = static_cast<int>(std::floor(static_cast<double>(res[i])));
        }
        return res;
    }

    bool triangle_t::isOnFace(int face) const {
        coord_t c(std::vector<double>(3));
        c.position = centroid();
        return c.isOnFace(face);
    }

    bool triangle_t::isOnAnyFace() const {
        return isOnFace(FACE_X) || isOnFace(FACE_Y) || isOnFace(FACE_Z);
    }

    std::vector<double> triangle_t::centroid() const {
        std::vector<double> res(3, 0.0);
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                res[j] += vertices[i].position[j];
            }
        }
        for (int i = 0; i < 3; ++i) {
            res[i] /= 3.0;
        }
        return res;
    }

} // namespace conformal_types_m

#endif // CONFORMAL_TYPES_H