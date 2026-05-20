#include <vector>
#include <cmath>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <iostream>

// Assuming these types are defined in conformal_types_m
// We need to forward declare or include them. 
// Since we don't have the content of conformal_types_m, we assume standard definitions 
// or that they are available in the global namespace or another included header.
// For the sake of translation, we assume:
// - triangle_t has a method getSides() returning std::vector<side_t>
// - side_t has members init and end of type coord_t
// - coord_t has a method position(int) returning double, and isOnVertex(), isOnAnyFace(), getCell(), getFace(), getEdge(), isEquiv(side_t), normal
// - NOT_ON_FACE, FACE_X, FACE_Y, FACE_Z, EDGE_X, EDGE_Y, EDGE_Z are integer constants

// Placeholder definitions to make the code compile if conformal_types_m is not provided.
// In a real scenario, these would come from the actual header.
struct coord_t {
    double position(int dim) const { return pos[dim]; }
    void position(int dim, double val) { pos[dim] = val; }
    bool isOnVertex() const { return false; } // Placeholder
    bool isOnAnyFace() const { return false; } // Placeholder
    std::vector<int> getCell() const { return {}; } // Placeholder
    int getFace() const { return 0; } // Placeholder
    int getEdge() const { return 0; } // Placeholder
    bool isEquiv(const struct side_t& other) const { return false; } // Placeholder
    double normal[3];
    double pos[3];
};

struct side_t {
    coord_t init;
    coord_t end;
    std::vector<side_t> getSides() const { return {}; } // Placeholder for triangle_t method
    bool isOnAnyFace() const { return false; }
    std::vector<int> getCell() const { return {}; }
    int getFace() const { return 0; }
    int getEdge() const { return 0; }
    bool isEquiv(const side_t& other) const { return false; }
};

struct triangle_t {
    std::vector<side_t> getSides() const { return {}; }
    bool isOnFace(int face) const { return false; }
};

// Constants
const int NOT_ON_FACE = -1;
const int FACE_X = 0;
const int FACE_Y = 1;
const int FACE_Z = 2;
const int EDGE_X = 0;
const int EDGE_Y = 1;
const int EDGE_Z = 2;

namespace geometry_m {

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
            if (sides_copy[i].init.position(edge) > sides_copy[i].end.position(edge)) { 
                c = sides_copy[i].init.position(edge);
                sides_copy[i].init.position(edge, sides_copy[i].end.position(edge));
                sides_copy[i].end.position(edge, c);
            }
        }
        res = sides_copy[0];
        for (int i = 1; i < (int)sides_copy.size(); ++i) {
            if (sides_copy[i].init.position(edge) < res.init.position(edge)) { 
                res.init.position(edge, sides_copy[i].init.position(edge));
            }
            if (sides_copy[i].end.position(edge) > res.end.position(edge)) { 
                res.end.position(edge, sides_copy[i].end.position(edge));
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
            if (sides[i].init.isEquiv(side)) { // Assuming isEquiv is on coord_t or side_t. Based on context, likely side_t has isEquiv or coord_t does. 
                // The Fortran code says sides(i)%isEquiv(side). side_t doesn't have isEquiv in my placeholder.
                // Let's assume side_t has isEquiv.
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
        std::vector<std::vector<double>> corners(3, std::vector<double>(4));

        corners = buildCorners(inner_path[0], inner_path[0].getFace());
        
        res.resize(inner_path.size() + 2);
        for(size_t i=0; i<inner_path.size(); ++i) res[i] = inner_path[i];
        
        // cornerIndex returns 1-based index in Fortran, so we need to adjust for 0-based C++
        // mod(cornerIndex(...), 4) + 1 gives 1..4. 
        // In C++, we want 0..3. So (cornerIndex(...) - 1) % 4.
        int ci = cornerIndex(corners, inner_path[inner_path.size()-1].end.position); // Assuming position is accessible
        // Wait, inner_path[inner_path.size()-1].end is a coord_t. 
        // Let's assume coord_t has a method to get position as vector or we pass the vector.
        // The Fortran code: cornerIndex(corners, inner_path(size(inner_path))%end%position)
        // position is a function returning a real. But cornerIndex takes a real(3).
        // So we need to extract the 3 coordinates.
        std::vector<double> vertex(3);
        vertex[0] = inner_path[inner_path.size()-1].end.position(0);
        vertex[1] = inner_path[inner_path.size()-1].end.position(1);
        vertex[2] = inner_path[inner_path.size()-1].end.position(2);
        
        mid_corner_idx = (cornerIndex(corners, vertex) - 1) % 4; // 0-based index for corners(:, mid_corner_idx)
        
        // buildSide takes two real(3) vectors
        std::vector<double> c1(3);
        c1[0] = inner_path[inner_path.size()-1].end.position(0);
        c1[1] = inner_path[inner_path.size()-1].end.position(1);
        c1[2] = inner_path[inner_path.size()-1].end.position(2);
        
        std::vector<double> c2(3);
        c2[0] = corners[0][mid_corner_idx];
        c2[1] = corners[1][mid_corner_idx];
        c2[2] = corners[2][mid_corner_idx];
        
        res[inner_path.size()] = buildSide(c1, c2);
        
        c1[0] = corners[0][mid_corner_idx];
        c1[1] = corners[1][mid_corner_idx];
        c1[2] = corners[2][mid_corner_idx];
        
        c2[0] = inner_path[0].init.position(0);
        c2[1] = inner_path[0].init.position(1);
        c2[2] = inner_path[0].init.position(2);
        
        res[inner_path.size() + 1] = buildSide(c1, c2);
        
        return res;
    }
    
    std::vector<side_t> buildVertexToSideContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        coord_t init = inner_path[0].init;
        coord_t end = inner_path[inner_path.size()-1].end;
        std::vector<std::vector<double>> corners = buildCorners(inner_path[0], inner_path[0].getFace());
        int idx = -1;

        res = inner_path;
        
        side_t cell_side;
        for (int i = 0; i < 4; ++i) {
            cell_side.init.position(0, corners[0][i]);
            cell_side.init.position(1, corners[1][i]);
            cell_side.init.position(2, corners[2][i]);
            
            int next_i = (i + 1) % 4;
            cell_side.end.position(0, corners[0][next_i]);
            cell_side.end.position(1, corners[1][next_i]);
            cell_side.end.position(2, corners[2][next_i]);
            
            std::vector<int> cell = cell_side.getCell();
            std::vector<int> floor_end(3);
            floor_end[0] = std::floor(end.position(0));
            floor_end[1] = std::floor(end.position(1));
            floor_end[2] = std::floor(end.position(2));
            
            if (cell == floor_end && (cell_side.getEdge() == end.getEdge())) { 
                idx = i;
                break;
            }
        }
        
        std::vector<double> c1(3);
        c1[0] = end.position(0);
        c1[1] = end.position(1);
        c1[2] = end.position(2);
        
        std::vector<double> c2(3);
        int next_idx = (idx + 1) % 4;
        c2[0] = corners[0][next_idx];
        c2[1] = corners[1][next_idx];
        c2[2] = corners[2][next_idx];
        
        addSide(res, buildSide(c1, c2));
        
        while (!(cell_side.getCell() == floor_end && cell_side.getEdge() == end.getEdge())) {
             // Re-evaluate cell_side for the loop condition? 
             // The Fortran code updates cell_side inside the loop.
             // Let's replicate the logic.
             
             // Current cell_side is from the previous iteration or initial setup?
             // In Fortran, cell_side is updated at the end of the loop.
             // The condition checks the CURRENT cell_side.
             // But wait, the loop condition is checked BEFORE the body.
             // The first check uses the cell_side defined in the for-loop? No, that was local to the for-loop scope in C++ if declared there.
             // In Fortran, cell_side is declared outside.
             
             // Let's restart the logic carefully.
             // We need to track the current corner index.
             int curr_idx = next_idx; // Start from the one we just added side to? 
             // No, the next side to add is from corners(:, mod(idx,4)+1) to corners(:, mod(idx+1,4)+1).
             // We already added side from end to corners(:, mod(idx,4)+1).
             // Now we need to add sides along the edge until we hit init.
             
             // Let's use a variable for the current corner index.
             int current_corner_idx = next_idx;
             
             while (true) {
                 std::vector<double> start_c(3);
                 start_c[0] = corners[0][current_corner_idx];
                 start_c[1] = corners[1][current_corner_idx];
                 start_c[2] = corners[2][current_corner_idx];
                 
                 int next_corner_idx = (current_corner_idx + 1) % 4;
                 std::vector<double> end_c(3);
                 end_c[0] = corners[0][next_corner_idx];
                 end_c[1] = corners[1][next_corner_idx];
                 end_c[2] = corners[2][next_corner_idx];
                 
                 cell_side.init.position(0, start_c[0]);
                 cell_side.init.position(1, start_c[1]);
                 cell_side.init.position(2, start_c[2]);
                 cell_side.end.position(0, end_c[0]);
                 cell_side.end.position(1, end_c[1]);
                 cell_side.end.position(2, end_c[2]);
                 
                 std::vector<int> cell = cell_side.getCell();
                 std::vector<int> floor_init(3);
                 floor_init[0] = std::floor(init.position(0));
                 floor_init[1] = std::floor(init.position(1));
                 floor_init[2] = std::floor(init.position(2));
                 
                 if (cell == floor_init && cell_side.getEdge() == init.getEdge()) {
                     break;
                 }
                 
                 addSide(res, buildSide(start_c, end_c));
                 current_corner_idx = next_corner_idx;
             }
             break; // Only one while loop needed
        }
        
        // The above logic is getting messy. Let's rewrite buildVertexToSideContour more cleanly.
        res.clear();
        res = inner_path;
        
        init = inner_path[0].init;
        end = inner_path[inner_path.size()-1].end;
        corners = buildCorners(inner_path[0], inner_path[0].getFace());
        
        int idx = -1;
        side_t temp_side;
        for (int i = 0; i < 4; ++i) {
            temp_side.init.position(0, corners[0][i]);
            temp_side.init.position(1, corners[1][i]);
            temp_side.init.position(2, corners[2][i]);
            int next_i = (i + 1) % 4;
            temp_side.end.position(0, corners[0][next_i]);
            temp_side.end.position(1, corners[1][next_i]);
            temp_side.end.position(2, corners[2][next_i]);
            
            std::vector<int> cell = temp_side.getCell();
            std::vector<int> floor_end(3);
            floor_end[0] = std::floor(end.position(0));
            floor_end[1] = std::floor(end.position(1));
            floor_end[2] = std::floor(end.position(2));
            
            if (cell == floor_end && (temp_side.getEdge() == end.getEdge())) { 
                idx = i;
                break;
            }
        }
        
        std::vector<double> c1(3);
        c1[0] = end.position(0);
        c1[1] = end.position(1);
        c1[2] = end.position(2);
        
        std::vector<double> c2(3);
        int next_idx = (idx + 1) % 4;
        c2[0] = corners[0][next_idx];
        c2[1] = corners[1][next_idx];
        c2[2] = corners[2][next_idx];
        
        addSide(res, buildSide(c1, c2));
        
        int curr = next_idx;
        while (true) {
            std::vector<double> start_c(3);
            start_c[0] = corners[0][curr];
            start_c[1] = corners[1][curr];
            start_c[2] = corners[2][curr];
            
            int next_curr = (curr + 1) % 4;
            std::vector<double> end_c(3);
            end_c[0] = corners[0][next_curr];
            end_c[1] = corners[1][next_curr];
            end_c[2] = corners[2][next_curr];
            
            temp_side.init.position(0, start_c[0]);
            temp_side.init.position(1, start_c[1]);
            temp_side.init.position(2, start_c[2]);
            temp_side.end.position(0, end_c[0]);
            temp_side.end.position(1, end_c[1]);
            temp_side.end.position(2, end_c[2]);
            
            std::vector<int> cell = temp_side.getCell();
            std::vector<int> floor_init(3);
            floor_init[0] = std::floor(init.position(0));
            floor_init[1] = std::floor(init.position(1));
            floor_init[2] = std::floor(init.position(2));
            
            if (cell == floor_init && temp_side.getEdge() == init.getEdge()) {
                break;
            }
            
            addSide(res, buildSide(start_c, end_c));
            curr = next_curr;
        }
        
        return res;
    }

    std::vector<side_t> buildSideToVertexContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        coord_t init = inner_path[0].init;
        coord_t end = inner_path[inner_path.size()-1].end;
        std::vector<std::vector<double>> corners = buildCorners(inner_path[0], inner_path[0].getFace());
        int idx = cornerIndex(corners, std::vector<double>{end.position(0), end.position(1), end.position(2)}) - 1; // 0-based

        res = inner_path;

        side_t cell_side;
        cell_side.init.position(0, corners[0][idx]);
        cell_side.init.position(1, corners[1][idx]);
        cell_side.init.position(2, corners[2][idx]);
        
        int next_idx = (idx + 1) % 4;
        cell_side.end.position(0, corners[0][next_idx]);
        cell_side.end.position(1, corners[1][next_idx]);
        cell_side.end.position(2, corners[2][next_idx]);
        
        std::vector<int> floor_init(3);
        floor_init[0] = std::floor(init.position(0));
        floor_init[1] = std::floor(init.position(1));
        floor_init[2] = std::floor(init.position(2));

        while (true) {
            std::vector<int> cell = cell_side.getCell();
            if (cell == floor_init && cell_side.getEdge() == init.getEdge()) {
                break;
            }
            
            std::vector<double> c1(3);
            c1[0] = cell_side.init.position(0);
            c1[1] = cell_side.init.position(1);
            c1[2] = cell_side.init.position(2);
            
            std::vector<double> c2(3);
            c2[0] = cell_side.end.position(0);
            c2[1] = cell_side.end.position(1);
            c2[2] = cell_side.end.position(2);
            
            addSide(res, buildSide(c1, c2));
            
            idx = next_idx;
            next_idx = (idx + 1) % 4;
            cell_side.init.position(0, corners[0][idx]);
            cell_side.init.position(1, corners[1][idx]);
            cell_side.init.position(2, corners[2][idx]);
            cell_side.end.position(0, corners[0][next_idx]);
            cell_side.end.position(1, corners[1][next_idx]);
            cell_side.end.position(2, corners[2][next_idx]);
        }
        
        std::vector<double> c1(3);
        c1[0] = cell_side.init.position(0);
        c1[1] = cell_side.init.position(1);
        c1[2] = cell_side.init.position(2);
        
        std::vector<double> c2(3);
        c2[0] = init.position(0);
        c2[1] = init.position(1);
        c2[2] = init.position(2);
        
        addSide(res, buildSide(c1, c2));

        return res;
    }
    
    std::vector<side_t> buildSideToSideContour(const std::vector<side_t>& inner_path) {
        std::vector<side_t> res;
        side_t cell_side;
        int idx_i = -1, idx_e = -1;
        std::vector<std::vector<double>> corners;
        coord_t init = inner_path[0].init;
        coord_t end = inner_path[inner_path.size()-1].end;

        corners = buildCorners(inner_path[0], inner_path[0].getFace());
        res = inner_path;
        
        for (int i = 0; i < 4; ++i) {
            cell_side.init.position(0, corners[0][i]);
            cell_side.init.position(1, corners[1][i]);
            cell_side.init.position(2, corners[2][i]);
            int next_i = (i + 1) % 4;
            cell_side.end.position(0, corners[0][next_i]);
            cell_side.end.position(1, corners[1][next_i]);
            cell_side.end.position(2, corners[2][next_i]);
            
            std::vector<int> cell = cell_side.getCell();
            std::vector<int> floor_init(3);
            floor_init[0] = std::floor(init.position(0));
            floor_init[1] = std::floor(init.position(1));
            floor_init[2] = std::floor(init.position(2));
            
            if (cell == floor_init && cell_side.getEdge() == init.getEdge()) { 
                idx_i = i;
            }
            
            std::vector<int> floor_end(3);
            floor_end[0] = std::floor(end.position(0));
            floor_end[1] = std::floor(end.position(1));
            floor_end[2] = std::floor(end.position(2));
            
            if (cell == floor_end && cell_side.getEdge() == end.getEdge()) { 
                idx_e = i;
            }
        }
        
        int idx = (idx_e + 1) % 4;
        
        std::vector<double> c1(3);
        c1[0] = end.position(0);
        c1[1] = end.position(1);
        c1[2] = end.position(2);
        
        std::vector<double> c2(3);
        c2[0] = corners[0][idx];
        c2[1] = corners[1][idx];
        c2[2] = corners[2][idx];
        
        addSide(res, buildSide(c1, c2));

        cell_side.init.position(0, corners[0][idx]);
        cell_side.init.position(1, corners[1][idx]);
        cell_side.init.position(2, corners[2][idx]);
        
        int next_idx = (idx + 1) % 4;
        cell_side.end.position(0, corners[0][next_idx]);
        cell_side.end.position(1, corners[1][next_idx]);
        cell_side.end.position(2, corners[2][next_idx]);
        
        std::vector<int> floor_init(3);
        floor_init[0] = std::floor(init.position(0));
        floor_init[1] = std::floor(init.position(1));
        floor_init[2] = std::floor(init.position(2));

        while (true) {
            std::vector<int> cell = cell_side.getCell();
            if (cell == floor_init && cell_side.getEdge() == init.getEdge()) {
                break;
            }
            
            std::vector<double> c1(3);
            c1[0] = cell_side.init.position(0);
            c1[1] = cell_side.init.position(1);
            c1[2] = cell_side.init.position(2);
            
            std::vector<double> c2(3);
            c2[0] = cell_side.end.position(0);
            c2[1] = cell_side.end.position(1);
            c2[2] = cell_side.end.position(2);
            
            addSide(res, buildSide(c1, c2));
            
            idx = next_idx;
            next_idx = (idx + 1) % 4;
            cell_side.init.position(0, corners[0][idx]);
            cell_side.init.position(1, corners[1][idx]);
            cell_side.init.position(2, corners[2][idx]);
            cell_side.end.position(0, corners[0][next_idx]);
            cell_side.end.position(1, corners[1][next_idx]);
            cell_side.end.position(2, corners[2][next_idx]);
        }
        
        std::vector<double> c1(3);
        c1[0] = cell_side.init.position(0);
        c1[1] = cell_side.init.position(1);
        c1[2] = cell_side.init.position(2);
        
        std::vector<double> c2(3);
        c2[0] = init.position(0);
        c2[1] = init.position(1);
        c2[2] = init.position(2);
        
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
        res.init.position(0, c1[0]);
        res.init.position(1, c1[1]);
        res.init.position(2, c1[2]);
        res.end.position(0, c2[0]);
        res.end.position(1, c2[1]);
        res.end.position(2, c2[2]);
        return res;
    }

    int cornerIndex(const std::vector<std::vector<double>>& corners, const std::vector<double>& vertex) {
        for (int i = 0; i < 4; ++i) {
            if (vertex[0] == corners[0][i] && vertex[1] == corners[1][i] && vertex[2] == corners[2][i]) { 
                return i + 1; // 1-based index
            }
        }
        return -1; // Should not happen
    }

    std::vector<std::vector<int>> buildCorners(const side_t& side, int face) {
        std::vector<std::vector<int>> res(3, std::vector<int>(4));
        std::vector<int> cell = side.getCell();
        std::vector<int> aux(3);

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
            aux = res[1];
            res[1] = res[3];
            res[3] = aux;
        }
        
        return res;
    }

    bool isClockwise(const side_t& side, int face) {
        bool result = true;
        std::vector<double> diff(3);
        diff[0] = side.end.position(0) - side.init.position(0);
        diff[1] = side.end.position(1) - side.init.position(1);
        diff[2] = side.end.position(2) - side.init.position(2);
        
        std::vector<double> x_prod = cross(diff, std::vector<double>{side.init.normal[0], side.init.normal[1], side.init.normal[2]});
        // Note: side.normal is an array in coord_t? Or side has a normal member?
        // In Fortran: side%normal. In my placeholder, coord_t has normal.
        // side.init is coord_t, side.end is coord_t.
        // Let's assume side has a normal member or we use init.normal.
        // The Fortran code uses side%normal.
        // If side_t doesn't have normal, we might need to add it.
        // For now, assuming side.init.normal is the normal.
        
        if (x_prod[face] < 0) result = false;
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
           res += aux_contour[i].init.position(dir1) * aux_contour[i].end.position(dir2) - 
                       aux_contour[i].end.position(dir1) * aux_contour[i].init.position(dir2);
        }
        res = 0.5 * res;
        return res;
    }
  

    std::vector<side_t> getPathOnFace(const std::vector<side_t>& sides) {
        std::vector<side_t> res(sides.size());
        int n = 0;
        while (n < (int)res.size()) {
           for (int i = 0; i < (int)sides.size(); ++i) { 
              if (n == 0) { 
                if(!sides[i].init.isOnAnyFace()) { 
                 n = n + 1;
                 res[n-1] = sides[i];
                }
              } else if (n != 0 && sides[i].init.position(0) == res[n-1].end.position(0) &&
                         sides[i].init.position(1) == res[n-1].end.position(1) &&
                         sides[i].init.position(2) == res[n-1].end.position(2)) { 
                 n = n + 1;
                 res[n-1] = sides[i];
              }
           }
        }
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