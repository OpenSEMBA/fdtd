#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iostream>

#include "conformal_types.h"

namespace geometry_m {

    using namespace conformal_types_m;

    bool positionsEqualTol(const std::vector<double>& a, const std::vector<double>& b) {
        return a.size() == 3u && b.size() == 3u &&
               std::abs(a[0] - b[0]) < POS_TOL && std::abs(a[1] - b[1]) < POS_TOL &&
               std::abs(a[2] - b[2]) < POS_TOL;
    }

    double contourArea(const std::vector<side_t>& contour, int orientation = NOT_ON_FACE);
    std::vector<side_t> getSidesOnFace(const std::vector<side_t>& sides, int face);
    std::vector<side_t> buildSidesContour(const std::vector<side_t>& sides);
    void addNewSides(std::vector<side_t>& sides, const std::vector<side_t>& new_sides);
    bool isNewSide(const std::vector<side_t>& sides, const side_t& side);
    void addNewSide(std::vector<side_t>& sides, const side_t& side);
    std::vector<side_t> getSidesOnEdge(const std::vector<side_t>& sides, int edge);
    std::vector<side_t> getPathOnFace(const std::vector<side_t>& sides);
    std::vector<side_t> buildVertexToVertexContour(const std::vector<side_t>& inner_path);
    std::vector<side_t> buildVertexToSideContour(const std::vector<side_t>& inner_path);
    std::vector<side_t> buildSideToVertexContour(const std::vector<side_t>& inner_path);
    std::vector<side_t> buildSideToSideContour(const std::vector<side_t>& inner_path);
    std::vector<std::vector<int>> buildCorners(const side_t& side, int face);
    int cornerIndex(const std::vector<std::vector<int>>& corners, const std::vector<double>& vertex);
    side_t buildSide(const std::vector<double>& c1, const std::vector<double>& c2);
    void addSide(std::vector<side_t>& sides, const side_t& side);
    bool isClockwise(const side_t& side, int face);

    double getArea(const triangle_t& triangle) {
        std::vector<side_t> sides = triangle.getSides();
        return contourArea(sides);
    }

    std::vector<double> cross(const std::vector<double>& v1, const std::vector<double>& v2) {
        std::vector<double> res(3);
        res[0] = v1[1]*v2[2]-v1[2]*v2[1];
        res[1] = -(v1[0]*v2[2]-v1[2]*v2[0]);
        res[2] = v1[0]*v2[1]-v1[1]*v2[0];
        return res;
    }

    side_t mergeSides(const std::vector<side_t>& sides, int edge) {
        std::vector<side_t> sides_copy = sides;
        double c;
        side_t res;
        for (int i = 0; i < (int)sides_copy.size(); ++i) {
            if (sides_copy[i].init.position[edge - 1] > sides_copy[i].end.position[edge - 1]) { 
                c = sides_copy[i].init.position[edge - 1];
                sides_copy[i].init.position[edge - 1] = sides_copy[i].end.position[edge - 1];
                sides_copy[i].end.position[edge - 1] = c;
            }
        }
        res = sides_copy[0];
        for (int i = 1; i < (int)sides_copy.size(); ++i) {
            if (sides_copy[i].init.position[edge - 1] < res.init.position[edge - 1]) { 
                res.init.position[edge - 1] = sides_copy[i].init.position[edge - 1];
            }
            if (sides_copy[i].end.position[edge - 1] > res.end.position[edge - 1]) { 
                res.end.position[edge - 1] = sides_copy[i].end.position[edge - 1];
            }
        }
        return res;
    }

    std::vector<int> findContourCell(const std::vector<side_t>& contour) {
        std::vector<int> res(3);
        for (int i = 0; i < (int)contour.size(); ++i) {
            if (contour[i].isOnAnyFace()) { 
                res = contour[i].getCell();
            }
        }
        return res;
    }

    int findContourFace(const std::vector<side_t>& contour) {
        int res = NOT_ON_FACE;
        for (int i = 0; i < (int)contour.size(); ++i) {
            if (contour[i].isOnAnyFace()) { 
                res = contour[i].getFace();
            }
        }
        if (res == NOT_ON_FACE) throw std::runtime_error("Contour face could not be identified");
        return res;
    }

    std::vector<side_t> buildCellSideSet(const std::vector<side_t>& sides, const std::vector<side_t>& on_sides) {
        std::vector<side_t> res;
        for (int face = FACE_X; face <= FACE_Z; ++face) {
            std::vector<side_t> sides_on_face = getSidesOnFace(sides, face);
            std::vector<side_t> contour = buildSidesContour(sides_on_face);
            // addNewSides modifies res in place, so we pass it by reference
            addNewSides(res, contour);
        }
        for (int edge = EDGE_X; edge <= EDGE_Z; ++edge) {
            std::vector<side_t> aux = getSidesOnEdge(sides, edge);
            addNewSides(res, aux);
            aux = getSidesOnEdge(on_sides, edge);
            addNewSides(res, aux);
        }
        return res;
    } 

    void addNewSides(std::vector<side_t>& sides, const std::vector<side_t>& new_sides) {
        for (int i = 0; i < (int)new_sides.size(); ++i) {
            if (isNewSide(sides, new_sides[i])) {
                addNewSide(sides, new_sides[i]);
            }
        }
    }

    bool isNewSide(const std::vector<side_t>& sides, const side_t& side) {
        bool result = true;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isEquiv(side)) {
                result = false;
            }
        }
        return result;
    }    

    void addNewSide(std::vector<side_t>& sides, const side_t& side) {
        if (sides.size() == 0) { 
            sides.clear();
            sides.push_back(side);
        } else { 
            std::vector<side_t> aux = sides;
            aux.push_back(side);
            sides = aux;
        }
    }

    std::vector<side_t> buildSidesContour(const std::vector<side_t>& sides) {
        std::vector<side_t> res;
        if (sides.size() == 0) { 
            res.clear();
        } else {
            std::vector<side_t> inner_path = getPathOnFace(sides);
            if (inner_path.empty()) {
                return res;
            }
            coord_t init = inner_path[0].init;
            coord_t end = inner_path[inner_path.size()-1].end;
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
        std::vector<std::vector<int>> corners;

        corners = buildCorners(inner_path[0], inner_path[0].getFace());
        
        res.resize(inner_path.size() + 2);
        for(size_t i=0; i<inner_path.size(); ++i) res[i] = inner_path[i];
        
        std::vector<double> vertex(3);
        vertex[0] = inner_path[inner_path.size()-1].end.position[0];
        vertex[1] = inner_path[inner_path.size()-1].end.position[1];
        vertex[2] = inner_path[inner_path.size()-1].end.position[2];
        
        int corner_idx = cornerIndex(corners, vertex);
        if (corner_idx < 0) {
            corner_idx = 1;
        }
        // Fortran: mod(cornerIndex, 4) + 1 (1-based column) -> 0-based index corner_idx % 4
        mid_corner_idx = corner_idx % 4;
        
        // buildSide takes two real(3) vectors
        std::vector<double> c1(3);
        c1[0] = inner_path[inner_path.size()-1].end.position[0];
        c1[1] = inner_path[inner_path.size()-1].end.position[1];
        c1[2] = inner_path[inner_path.size()-1].end.position[2];
        
        std::vector<double> c2(3);
        c2[0] = static_cast<double>(corners[0][mid_corner_idx]);
        c2[1] = static_cast<double>(corners[1][mid_corner_idx]);
        c2[2] = static_cast<double>(corners[2][mid_corner_idx]);
        
        res[inner_path.size()] = buildSide(c1, c2);
        
        c1[0] = static_cast<double>(corners[0][mid_corner_idx]);
        c1[1] = static_cast<double>(corners[1][mid_corner_idx]);
        c1[2] = static_cast<double>(corners[2][mid_corner_idx]);
        
        c2[0] = inner_path[0].init.position[0];
        c2[1] = inner_path[0].init.position[1];
        c2[2] = inner_path[0].init.position[2];
        
        res[inner_path.size() + 1] = buildSide(c1, c2);
        
        return res;
    }
    
    std::vector<side_t> buildVertexToSideContour(const std::vector<side_t>& inner_path) {
        coord_t init = inner_path[0].init;
        coord_t end = inner_path[inner_path.size()-1].end;
        std::vector<std::vector<int>> corners = buildCorners(inner_path[0], inner_path[0].getFace());
        int idx = 0;

        std::vector<side_t> res = inner_path;

        side_t cell_side;
        for (int i = 1; i <= 4; ++i) {
            const int col0 = i - 1;
            const int col1 = i % 4;
            cell_side.init.position[0] = static_cast<double>(corners[0][col0]);
            cell_side.init.position[1] = static_cast<double>(corners[1][col0]);
            cell_side.init.position[2] = static_cast<double>(corners[2][col0]);
            cell_side.end.position[0] = static_cast<double>(corners[0][col1]);
            cell_side.end.position[1] = static_cast<double>(corners[1][col1]);
            cell_side.end.position[2] = static_cast<double>(corners[2][col1]);

            const std::vector<int> cell = cell_side.getCell();
            const std::vector<int> floor_end = {
                static_cast<int>(std::floor(end.position[0])),
                static_cast<int>(std::floor(end.position[1])),
                static_cast<int>(std::floor(end.position[2])),
            };

            if (cell == floor_end && cell_side.getEdge() == end.getEdge()) {
                idx = i;
                break;
            }
        }

        auto cornerPos = [&](int col0) {
            return std::vector<double>{
                static_cast<double>(corners[0][col0]),
                static_cast<double>(corners[1][col0]),
                static_cast<double>(corners[2][col0]),
            };
        };

        addSide(res, buildSide(
            std::vector<double>{end.position[0], end.position[1], end.position[2]},
            cornerPos(idx % 4)));

        while (cornerPos(idx % 4)[0] != init.position[0] ||
               cornerPos(idx % 4)[1] != init.position[1] ||
               cornerPos(idx % 4)[2] != init.position[2]) {
            const int col0 = idx % 4;
            const int col1 = (idx + 1) % 4;
            addSide(res, buildSide(cornerPos(col0), cornerPos(col1)));
            ++idx;
        }

        return res;
    }

    std::vector<side_t> buildSideToVertexContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        coord_t init = inner_path[0].init;
        coord_t end = inner_path[inner_path.size()-1].end;
        std::vector<std::vector<int>> corners = buildCorners(inner_path[0], inner_path[0].getFace());
        int idx = cornerIndex(corners, std::vector<double>{end.position[0], end.position[1], end.position[2]}) - 1;
        if (idx < 0) {
            idx = 0;
        }

        res = inner_path;

        side_t cell_side;
        cell_side.init.position[0] = static_cast<double>(corners[0][idx]);
        cell_side.init.position[1] = static_cast<double>(corners[1][idx]);
        cell_side.init.position[2] = static_cast<double>(corners[2][idx]);
        
        int next_idx = (idx + 1) % 4;
        cell_side.end.position[0] = static_cast<double>(corners[0][next_idx]);
        cell_side.end.position[1] = static_cast<double>(corners[1][next_idx]);
        cell_side.end.position[2] = static_cast<double>(corners[2][next_idx]);
        
        std::vector<int> floor_init(3);
        floor_init[0] = static_cast<int>(std::floor(init.position[0]));
        floor_init[1] = static_cast<int>(std::floor(init.position[1]));
        floor_init[2] = static_cast<int>(std::floor(init.position[2]));

        for (int step = 0; step < 4; ++step) {
            std::vector<int> cell = cell_side.getCell();
            if (cell == floor_init && cell_side.getEdge() == init.getEdge()) {
                break;
            }
            
            std::vector<double> c1(3);
            c1[0] = cell_side.init.position[0];
            c1[1] = cell_side.init.position[1];
            c1[2] = cell_side.init.position[2];
            
            std::vector<double> c2(3);
            c2[0] = cell_side.end.position[0];
            c2[1] = cell_side.end.position[1];
            c2[2] = cell_side.end.position[2];
            
            addSide(res, buildSide(c1, c2));
            
            idx = next_idx;
            next_idx = (idx + 1) % 4;
            cell_side.init.position[0] = static_cast<double>(corners[0][idx]);
            cell_side.init.position[1] = static_cast<double>(corners[1][idx]);
            cell_side.init.position[2] = static_cast<double>(corners[2][idx]);
            cell_side.end.position[0] = static_cast<double>(corners[0][next_idx]);
            cell_side.end.position[1] = static_cast<double>(corners[1][next_idx]);
            cell_side.end.position[2] = static_cast<double>(corners[2][next_idx]);
        }
        
        std::vector<double> c1(3);
        c1[0] = cell_side.init.position[0];
        c1[1] = cell_side.init.position[1];
        c1[2] = cell_side.init.position[2];
        
        std::vector<double> c2(3);
        c2[0] = init.position[0];
        c2[1] = init.position[1];
        c2[2] = init.position[2];
        
        addSide(res, buildSide(c1, c2));

        return res;
    }
    
    std::vector<side_t> buildSideToSideContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        side_t cell_side;
        int idx_i = -1, idx_e = -1;
        std::vector<std::vector<int>> corners;
        coord_t init = inner_path[0].init;
        coord_t end = inner_path[inner_path.size()-1].end;

        corners = buildCorners(inner_path[0], inner_path[0].getFace());
        res = inner_path;
        
        for (int i = 0; i < 4; ++i) {
            cell_side.init.position[0] = static_cast<double>(corners[0][i]);
            cell_side.init.position[1] = static_cast<double>(corners[1][i]);
            cell_side.init.position[2] = static_cast<double>(corners[2][i]);
            int next_i = (i + 1) % 4;
            cell_side.end.position[0] = static_cast<double>(corners[0][next_i]);
            cell_side.end.position[1] = static_cast<double>(corners[1][next_i]);
            cell_side.end.position[2] = static_cast<double>(corners[2][next_i]);
            
            std::vector<int> cell = cell_side.getCell();
            std::vector<int> floor_init(3);
            floor_init[0] = static_cast<int>(std::floor(init.position[0]));
            floor_init[1] = static_cast<int>(std::floor(init.position[1]));
            floor_init[2] = static_cast<int>(std::floor(init.position[2]));
            
            if (cell == floor_init && cell_side.getEdge() == init.getEdge()) { 
                idx_i = i;
            }
            
            std::vector<int> floor_end(3);
            floor_end[0] = static_cast<int>(std::floor(end.position[0]));
            floor_end[1] = static_cast<int>(std::floor(end.position[1]));
            floor_end[2] = static_cast<int>(std::floor(end.position[2]));
            
            if (cell == floor_end && cell_side.getEdge() == end.getEdge()) { 
                idx_e = i;
            }
        }
        
        int idx = (idx_e + 1) % 4;
        
        std::vector<double> c1(3);
        c1[0] = end.position[0];
        c1[1] = end.position[1];
        c1[2] = end.position[2];
        
        std::vector<double> c2(3);
        c2[0] = static_cast<double>(corners[0][idx]);
        c2[1] = static_cast<double>(corners[1][idx]);
        c2[2] = static_cast<double>(corners[2][idx]);
        
        addSide(res, buildSide(c1, c2));

        cell_side.init.position[0] = static_cast<double>(corners[0][idx]);
        cell_side.init.position[1] = static_cast<double>(corners[1][idx]);
        cell_side.init.position[2] = static_cast<double>(corners[2][idx]);
        
        int next_idx = (idx + 1) % 4;
        cell_side.end.position[0] = static_cast<double>(corners[0][next_idx]);
        cell_side.end.position[1] = static_cast<double>(corners[1][next_idx]);
        cell_side.end.position[2] = static_cast<double>(corners[2][next_idx]);
        
        std::vector<int> floor_init(3);
        floor_init[0] = static_cast<int>(std::floor(init.position[0]));
        floor_init[1] = static_cast<int>(std::floor(init.position[1]));
        floor_init[2] = static_cast<int>(std::floor(init.position[2]));

        for (int step = 0; step < 4; ++step) {
            std::vector<int> cell = cell_side.getCell();
            if (cell == floor_init && cell_side.getEdge() == init.getEdge()) {
                break;
            }
            
            c1[0] = cell_side.init.position[0];
            c1[1] = cell_side.init.position[1];
            c1[2] = cell_side.init.position[2];
            
            c2[0] = cell_side.end.position[0];
            c2[1] = cell_side.end.position[1];
            c2[2] = cell_side.end.position[2];
            
            addSide(res, buildSide(c1, c2));
            
            idx = next_idx;
            next_idx = (idx + 1) % 4;
            cell_side.init.position[0] = static_cast<double>(corners[0][idx]);
            cell_side.init.position[1] = static_cast<double>(corners[1][idx]);
            cell_side.init.position[2] = static_cast<double>(corners[2][idx]);
            cell_side.end.position[0] = static_cast<double>(corners[0][next_idx]);
            cell_side.end.position[1] = static_cast<double>(corners[1][next_idx]);
            cell_side.end.position[2] = static_cast<double>(corners[2][next_idx]);
        }
        
        c1[0] = cell_side.init.position[0];
        c1[1] = cell_side.init.position[1];
        c1[2] = cell_side.init.position[2];
        
        c2[0] = init.position[0];
        c2[1] = init.position[1];
        c2[2] = init.position[2];
        
        addSide(res, buildSide(c1, c2));

        return res;
    }

    void addSide(std::vector<side_t>& sides, const side_t& side) { 
        std::vector<side_t> aux = sides;
        aux.push_back(side);
        sides = aux;        
    }

    side_t buildSide(const std::vector<double>& c1, const std::vector<double>& c2) {
        side_t res;
        res.init.position[0] = c1[0];
        res.init.position[1] = c1[1];
        res.init.position[2] = c1[2];
        res.end.position[0] = c2[0];
        res.end.position[1] = c2[1];
        res.end.position[2] = c2[2];
        return res;
    }

    int cornerIndex(const std::vector<std::vector<int>>& corners, const std::vector<double>& vertex) {
        for (int i = 0; i < 4; ++i) {
            if (vertex[0] == static_cast<double>(corners[0][i]) &&
                vertex[1] == static_cast<double>(corners[1][i]) &&
                vertex[2] == static_cast<double>(corners[2][i])) {
                return i + 1;
            }
        }
        return -1;
    }

    std::vector<std::vector<int>> buildCorners(const side_t& side, int face) {
        std::vector<std::vector<int>> res(3, std::vector<int>(4));
        std::vector<int> cell = side.getCell();

        if (face == FACE_X) { 
            res[0][0] = cell[0]; res[1][0] = cell[1]; res[2][0] = cell[2];
            res[0][1] = cell[0]; res[1][1] = cell[1]+1; res[2][1] = cell[2];
            res[0][2] = cell[0]; res[1][2] = cell[1]+1; res[2][2] = cell[2]+1;
            res[0][3] = cell[0]; res[1][3] = cell[1]; res[2][3] = cell[2]+1;
        } else if (face == FACE_Y) { 
            res[0][0] = cell[0]; res[1][0] = cell[1]; res[2][0] = cell[2];
            res[0][3] = cell[0]+1; res[1][3] = cell[1]; res[2][3] = cell[2];
            res[0][2] = cell[0]+1; res[1][2] = cell[1]; res[2][2] = cell[2]+1;
            res[0][1] = cell[0]; res[1][1] = cell[1]; res[2][1] = cell[2]+1;
        } else if (face == FACE_Z) { 
            res[0][0] = cell[0]; res[1][0] = cell[1]; res[2][0] = cell[2];
            res[0][1] = cell[0]+1; res[1][1] = cell[1]; res[2][1] = cell[2];
            res[0][2] = cell[0]+1; res[1][2] = cell[1]+1; res[2][2] = cell[2];
            res[0][3] = cell[0]; res[1][3] = cell[1]+1; res[2][3] = cell[2];
        }

        if (isClockwise(side, face)) {
            for (int r = 0; r < 3; ++r) {
                std::swap(res[r][1], res[r][3]);
            }
        }
        
        return res;
    }

    bool isClockwise(const side_t& side, int face) {
        bool result = true;
        std::vector<double> diff(3);
        diff[0] = side.end.position[0] - side.init.position[0];
        diff[1] = side.end.position[1] - side.init.position[1];
        diff[2] = side.end.position[2] - side.init.position[2];
        
        std::vector<double> x_prod = cross(diff, side.normal);
        
        if (x_prod[face - 1] < 0) result = false;
        return result;
    }

    double contourArea(const std::vector<side_t>& contour, int orientation) {
        std::vector<side_t> aux_contour = contour;
        int face;
        if (orientation != 0) { // present(orientation) check is tricky in C++. 
            // The Fortran code: if (present(orientation)) then face = orientation else ...
            // In C++, we can't easily check if an optional argument was passed if it's a value.
            // However, the function signature in Fortran is: function contourArea(contour, orientation) result(res)
            // integer, optional :: orientation
            // If we call it with one argument, orientation is not present.
            // In C++, we can use a flag or a special value.
            // Let's assume the caller passes a valid face index or NOT_ON_FACE if not present.
            // But NOT_ON_FACE is -1.
            // Let's change the signature to take an optional int via a pointer or reference, or use a default value that indicates "not present".
            // Since I can't change the signature arbitrarily without breaking the "preserve names" rule if it implies interface compatibility,
            // I will assume the C++ caller handles the optional nature.
            // For this translation, I'll assume orientation is passed. If it's NOT_ON_FACE, we treat it as not present.
            if (orientation == NOT_ON_FACE) {
                 // Not present
                 for (int i = 0; i < (int)contour.size(); ++i) {
                    face = contour[i].getFace();
                    if (face != NOT_ON_FACE) break;
                 }
            } else {
                face = orientation;
            }
        } else {
             for (int i = 0; i < (int)contour.size(); ++i) {
                face = contour[i].getFace();
                if (face != NOT_ON_FACE) break;
             }
        }
        
        // Re-implementing the logic properly for C++
        // The Fortran code:
        // if (present(orientation)) then 
        //     face = orientation
        // else
        //     do i = 1, size(contour)
        //         face = contour(i)%getFace()
        //         if (face /= NOT_ON_FACE) exit
        //     end do
        // end if
        
        // Let's assume the C++ function is called with a second argument.
        // If the user wants to simulate "optional", they might pass a specific value.
        // I will add a boolean flag to the C++ function to indicate presence.
        // But the rule says "Preserve ALL names". So I must keep the signature.
        // I will assume that if orientation is NOT_ON_FACE, it is considered "not present".
        
        if (orientation == NOT_ON_FACE) {
             for (int i = 0; i < (int)contour.size(); ++i) {
                face = contour[i].getFace();
                if (face != NOT_ON_FACE) break;
             }
        } else {
            face = orientation;
        }

        if (isClockwise(contour[0], face)) { 
            for (int i = 0; i < (int)contour.size(); ++i) {
                aux_contour[contour.size() - 1 - i].init = contour[i].end;
                aux_contour[contour.size() - 1 - i].end = contour[i].init;
            }
        }

        int dir1 = face % 3 + 1; // 1-based index
        int dir2 = (face + 1) % 3 + 1; // 1-based index
        
        double res = 0;
        for (int i = 0; i < (int)aux_contour.size(); ++i) {
           res += aux_contour[i].init.position[dir1 - 1] * aux_contour[i].end.position[dir2 - 1] - 
                       aux_contour[i].end.position[dir1 - 1] * aux_contour[i].init.position[dir2 - 1];
        }
        res = 0.5 * res;
        return res;
    }
  

    std::vector<side_t> getPathOnFace(const std::vector<side_t>& sides) {
        if (sides.empty()) {
            return {};
        }
        std::vector<side_t> res(sides.size());
        int n = 0;
        while (n < static_cast<int>(sides.size())) {
            const int prev_n = n;
            for (int i = 0; i < static_cast<int>(sides.size()); ++i) {
                if (n == 0) {
                    if (!sides[i].init.isOnAnyFace()) {
                        n++;
                        res[n - 1] = sides[i];
                    }
                } else if (positionsEqualTol(sides[i].init.position, res[n - 1].end.position)) {
                    n++;
                    res[n - 1] = sides[i];
                }
            }
            if (n == prev_n) {
                break;
            }
        }
        res.resize(static_cast<std::size_t>(n));
        return res;
    }

    int getCellDistance(const std::vector<int>& ref_cell, const std::vector<int>& cell, int edge) {
        int res = 0;
        for (int i = 0; i < 3; ++i) {
            res += (i+1)*std::abs(ref_cell[i]-cell[i]);
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
        int n = 0;
        for (int i = 0; i < (int)triangles.size(); ++i) {
           if (!(triangles[i].isOnFace(1) || triangles[i].isOnFace(2) || triangles[i].isOnFace(3))) n++;
        }
        res.resize(n);
        n = 0;
        for (int i = 0; i < (int)triangles.size(); ++i) {
           if (!(triangles[i].isOnFace(1) || triangles[i].isOnFace(2) || triangles[i].isOnFace(3))) { 
              n++;
              res[n-1] = triangles[i]; // Note: Fortran code had 'triangles(j)' which is likely a typo for 'triangles(i)'
           }
        }
        return res;
    }
 
    std::vector<triangle_t> getTrianglesOnFace(const std::vector<triangle_t>& tris, int face) {
        std::vector<triangle_t> res;
        int n = 0;
        for (int i = 0; i < (int)tris.size(); ++i) {
            if (tris[i].isOnFace(face)) n++;
        }
        res.resize(n);
        n = 0;
        for (int i = 0; i < (int)tris.size(); ++i) {
            if (tris[i].isOnFace(face)) { 
                n++;
                res[n-1] = tris[i];
            }
        }
        return res;
    }

    std::vector<side_t> getSidesOnFace(const std::vector<side_t>& sides, int face) {
        std::vector<side_t> res;
        int n = 0;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isOnFace(face)) {
                n++;
            }
        }
        res.resize(n);
        n = 0;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isOnFace(face)) {
                n++;
                res[n-1] = sides[i];
            }
        }
        return res;
    }
 
    std::vector<side_t> getSidesOnEdge(const std::vector<side_t>& sides, int edge) {
        std::vector<side_t> res;
        int n = 0;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isOnEdge(edge)) {
                n++;
            }
        }
        res.resize(n);
        n = 0;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isOnEdge(edge)) {
                n++;
                res[n-1] = sides[i];
            }
        }
        return res;
    }
 
    std::vector<side_t> getSidesOnAdjacentEdges(const std::vector<side_t>& sides, int face) {
        std::vector<side_t> res;
        int n = 0;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isOnEdge(face % 3 + 1) || 
                sides[i].isOnEdge((face + 1) % 3 + 1)) {
                n++;
            }
        }
        res.resize(n);
        n = 0;
        for (int i = 0; i < (int)sides.size(); ++i) {
            if (sides[i].isOnEdge(face % 3 + 1) || 
                sides[i].isOnEdge((face + 1) % 3 + 1)) {
                n++;
                res[n-1] = sides[i];
            }
        }
        return res;
    }

}