#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <typeinfo>
#include <memory>
#include <stdexcept>

// Forward declarations and includes for external dependencies
// #include "fhash.h" // Assuming fhash_tbl_t is defined here
// #include "cells_m.h" // Assuming cell_region_t, cell_interval_t, linel_t, pixel_t are defined here
// #include "geometry_m.h" // Assuming triangle_t is defined here

// Constants from Fortran
#ifdef CompileWithSMBJSON

namespace mesh_m {

    // Constants
    constexpr int MAX_LINE = 256;
    constexpr int REGION_TYPE_VOLUME = 3;
    constexpr int REGION_TYPE_SURFACE = 2;

    // Assuming DIR_X, DIR_Z are defined elsewhere, typically 0 and 2 for x, y, z indices
    // If not, they might need to be defined or passed. 
    // Based on context, they are likely integer constants.
    // Let's assume they are defined in a common header or here if not present.
    // Since the prompt says preserve names, I will assume they are available in scope or defined elsewhere.
    // For compilation to succeed in a vacuum, I'll add dummy definitions if needed, 
    // but strictly following rules, I should not invent names. 
    // However, 'DIR_X' and 'DIR_Z' are used in loops. 
    // I will assume they are integers defined elsewhere. 
    // If this code is standalone, it won't compile without them. 
    // I will add a comment.
    // extern int DIR_X;
    // extern int DIR_Z;

    // Base Element Type
    struct element_t {
        std::vector<int> coordIds;
    };

    // Node Type
    struct node_t : public element_t {
        // coordIds must be size 1.
    };

    // Polyline Type
    struct polyline_t : public element_t {
        // coordIds must be size >1.
    };

    // Coordinate Type
    struct coordinate_t {
        std::vector<double> position; // real, dimension(3)

        coordinate_t() : position(3, 0.0) {}
        
        // Operator -
        coordinate_t operator-(const coordinate_t& b) const {
            coordinate_t res;
            for (int i = 0; i < 3; ++i) {
                res.position[i] = position[i] - b.position[i];
            }
            return res;
        }

        // Operator ==
        bool operator==(const coordinate_t& b) const {
            return (position[0] == b.position[0]) && 
                   (position[1] == b.position[1]) && 
                   (position[2] == b.position[2]);
        }
    };

    // Conformal Region Type
    struct conformal_region_t {
        std::vector<triangle_t> triangles;
        std::vector<cell_interval_t> intervals;
        int type;
    };

    // Mesh Type
    class mesh_t {
    private:
        fhash_tbl_t coordinates; // Map of CoordinateIds to relative coordinates.
        fhash_tbl_t elements;    // Map of ElementIds to elements/cellsRegions.    

    public:
        // Constructor/Destructor if needed
        mesh_t() {}
        ~mesh_t() {}

        // Methods
        void allocateCoordinates(int buck) {
            coordinates.allocate(buck);
        }

        void allocateElements(int buck) {
            elements.allocate(buck);
        }

        void printHashInfo() {
            std::cout << "Coordinates hash info:" << std::endl;
            int num_buckets, num_items, num_collisions, max_depth;
            coordinates.stats(num_buckets, num_items, num_collisions, max_depth);
            std::cout << "  Number of buckets allocated: " << num_buckets << std::endl;
            std::cout << "  Number of key-value pairs stored: " << num_items << std::endl;
            std::cout << "  Total number of hash-collisions: " << num_collisions << std::endl;
            std::cout << "  The worst case bucket depth is " << max_depth << std::endl;
            
            std::cout << "Elements hash info:" << std::endl;
            elements.stats(num_buckets, num_items, num_collisions, max_depth);
            std::cout << "  Number of buckets allocated: " << num_buckets << std::endl;
            std::cout << "  Number of key-value pairs stored: " << num_items << std::endl;
            std::cout << "  Total number of hash-collisions: " << num_collisions << std::endl;
            std::cout << "  The worst case bucket depth is " << max_depth << std::endl;
        }

        int checkCoordinateId(int id) {
            int stat;
            // Assuming key() creates a key object from int
            coordinates.check_key(key(id), stat);
            return stat;
        }
        
        int checkElementId(int id) {
            int stat;
            elements.check_key(key(id), stat);
            return stat;
        }

        void addCoordinate(int id, const coordinate_t& coordinate) {
            coordinates.set(key(id), coordinate);
        }

        void addElement(int id, const element_t& e) {
            elements.set(key(id), e);
        }

        void addCellRegion(int id, const cell_region_t& e) {
            elements.set(key(id), e);
        }

        void addConformalRegion(int id, const conformal_region_t& e) {
            elements.set(key(id), e);
        }

        coordinate_t getCoordinate(int id, bool* found = nullptr) {
            coordinate_t res;
            int stat;
            std::unique_ptr<coordinate_t> d; // Using unique_ptr to simulate allocatable class(*)

            if (found) *found = false;

            // get_raw needs to retrieve the value into d
            // Assuming fhash_tbl_t has a method that returns a void* or similar, 
            // or we need a template helper. 
            // Since C++ doesn't have class(*), we assume the hash table stores void* or std::any.
            // For this translation, we assume get_raw handles the casting or we use a generic pointer.
            // Let's assume get_raw takes a reference to a void* or similar.
            // However, to keep it simple and close to Fortran semantics:
            // We'll assume a helper or that the hash table is typed.
            // Given the complexity of fhash_tbl_t, I will assume it has a method:
            // void get_raw(const Key& k, void*& ptr, int& stat)
            
            void* raw_ptr = nullptr;
            coordinates.get_raw(key(id), raw_ptr, stat);
            
            if (stat != 0) return res;

            // Dynamic cast or type check
            // Since we don't know the exact type stored without RTTI or a tag,
            // we assume the hash table stores the correct type or we check typeid.
            // For simplicity, assuming it stores coordinate_t pointers.
            coordinate_t* c_ptr = dynamic_cast<coordinate_t*>(raw_ptr);
            if (c_ptr) {
                res = *c_ptr;
                if (found) *found = true;
            }
            
            return res;
        }

        node_t getNode(int id, bool* found = nullptr) {
            node_t res;
            int status;
            void* raw_ptr = nullptr;

            if (found) *found = false;
            
            elements.get_raw(key(id), raw_ptr, status);
            if (status != 0) return res;

            node_t* n_ptr = dynamic_cast<node_t*>(raw_ptr);
            if (n_ptr) {
                res = *n_ptr;
                if (found) *found = true;
            }
            
            return res;
        }

        polyline_t getPolyline(int id, bool* found = nullptr) {
            polyline_t res;
            int stat;
            void* raw_ptr = nullptr;

            if (found) *found = false;
            
            elements.get_raw(key(id), raw_ptr, stat);
            if (stat != 0) return res;

            polyline_t* p_ptr = dynamic_cast<polyline_t*>(raw_ptr);
            if (p_ptr) {
                res = *p_ptr;
                if (found) *found = true;
            }
            
            return res;
        }

        cell_region_t getCellRegion(int id, bool* found = nullptr) {
            cell_region_t res;
            int stat;
            void* raw_ptr = nullptr;

            if (found) *found = false;
            
            elements.get_raw(key(id), raw_ptr, stat);
            if (stat != 0) return res;

            // Check for cell_region_t first
            cell_region_t* cr_ptr = dynamic_cast<cell_region_t*>(raw_ptr);
            if (cr_ptr) {
                res = *cr_ptr;
                if (found) *found = true;
                return res;
            }
            
            // Check for conformal_region_t
            conformal_region_t* cfr_ptr = dynamic_cast<conformal_region_t*>(raw_ptr);
            if (cfr_ptr) {
                if (cfr_ptr->intervals.empty()) return res;
                res.intervals = cfr_ptr->intervals;
                if (found) *found = true;
            }
            
            return res;
        }

        std::vector<cell_region_t> getCellRegions(const std::vector<int>& ids) {
            std::vector<cell_region_t> res;
            cell_region_t cR;
            bool found;
            int numberOfCellRegions = 0;

            // Precounts
            for (int i = 0; i < ids.size(); ++i) {
                cR = getCellRegion(ids[i], &found);
                if (found) {
                    numberOfCellRegions++;
                }
            }
            
            res.resize(numberOfCellRegions);
            int j = 0;
            for (int i = 0; i < ids.size(); ++i) {
                cR = getCellRegion(ids[i], &found);
                if (found) {
                    res[j] = cR;
                    j++;
                }
            }

            return res;
        }

        conformal_region_t getConformalRegion(int id, bool* found = nullptr) {
            conformal_region_t res;
            int stat;
            void* raw_ptr = nullptr;

            if (found) *found = false;
            
            elements.get_raw(key(id), raw_ptr, stat);
            if (stat != 0) return res;

            // Check for cell_region_t
            cell_region_t* cr_ptr = dynamic_cast<cell_region_t*>(raw_ptr);
            if (cr_ptr) {
                return res; // Return empty/default
            }
            
            // Check for conformal_region_t
            conformal_region_t* cfr_ptr = dynamic_cast<conformal_region_t*>(raw_ptr);
            if (cfr_ptr) {
                res = *cfr_ptr;
                if (found) *found = true;
            }
            
            return res;
        }

        std::vector<conformal_region_t> getConformalRegions(const std::vector<int>& ids) {
            std::vector<conformal_region_t> res;
            conformal_region_t cR;
            bool found;
            int numberOfConformalRegions = 0;

            // Precounts
            for (int i = 0; i < ids.size(); ++i) {
                cR = getConformalRegion(ids[i], &found);
                if (found) {
                    numberOfConformalRegions++;
                }
            }
            
            res.resize(numberOfConformalRegions);
            int j = 0;
            for (int i = 0; i < ids.size(); ++i) {
                cR = getConformalRegion(ids[i], &found);
                if (found) {
                    res[j] = cR;
                    j++;
                }
            }

            return res;
        }

        int countPolylineSegments(const polyline_t& pl) {
            int res = 0;
            for (int i = 0; i < pl.coordIds.size() - 1; ++i) {
                coordinate_t iC = getCoordinate(pl.coordIds[i]);
                coordinate_t eC = getCoordinate(pl.coordIds[i+1]);
                cell_interval_t interval;
                // Assuming int() cast for position
                for(int k=0; k<3; ++k) {
                    interval.ini.cell[k] = static_cast<int>(iC.position[k]);
                    interval.end.cell[k] = static_cast<int>(eC.position[k]);
                }
                res += interval.getSize();
            }
            return res;
        }

        bool arePolylineSegmentsStructured(const polyline_t& pl) {
            bool res = true;
            for (int i = 0; i < pl.coordIds.size() - 1; ++i) {
                coordinate_t iC = getCoordinate(pl.coordIds[i]);
                coordinate_t eC = getCoordinate(pl.coordIds[i+1]);
                
                // Check if position is integer
                bool is_integer_i = true;
                bool is_integer_e = true;
                for(int k=0; k<3; ++k) {
                    if (std::floor(iC.position[k]) != iC.position[k]) is_integer_i = false;
                    if (std::floor(eC.position[k]) != eC.position[k]) is_integer_e = false;
                }
                
                if (!is_integer_i || !is_integer_e) {
                    res = false;
                    return res;
                } 

                int numberOfVaryingDirections = 0;
                // Assuming DIR_X and DIR_Z are 0 and 2
                for (int d = 0; d < 3; ++d) { // Iterate all dimensions, assuming DIR_X=0, DIR_Z=2
                    if (iC.position[d] != eC.position[d]) {
                        numberOfVaryingDirections++;
                    }
                }
                if (numberOfVaryingDirections > 1) {
                    res = false;
                    return res;
                }
            }

            res = true;
            return res;
        }

        std::vector<linel_t> polylineToLinels(const polyline_t& pl) {
            std::vector<linel_t> res;
            
            if (!arePolylineSegmentsStructured(pl)) {
                res.resize(0);
                return res;
            }

            int num_segments = countPolylineSegments(pl);
            res.resize(num_segments);
            if (res.empty()) return res;

            int lastSegment = 0; // 0-indexed
            for (int i = 0; i < pl.coordIds.size() - 1; ++i) {
                coordinate_t iC = getCoordinate(pl.coordIds[i]);
                coordinate_t eC = getCoordinate(pl.coordIds[i+1]);
                cell_interval_t interval;
                for(int k=0; k<3; ++k) {
                    interval.ini.cell[k] = static_cast<int>(iC.position[k]);
                    interval.end.cell[k] = static_cast<int>(eC.position[k]);
                }
                
                bool positions_differ = false;
                for(int k=0; k<3; ++k) {
                    if (iC.position[k] != eC.position[k]) {
                        positions_differ = true;
                        break;
                    }
                }

                if (positions_differ) {
                    std::vector<int> segment(3);
                    int size = interval.getSize();
                    for(int k=0; k<3; ++k) {
                        segment[k] = (interval.end.cell[k] - interval.ini.cell[k]) / size;
                    }
                    
                    res[lastSegment].tag = pl.coordIds[i];
                    for (int j = 0; j < size; ++j) {
                        coordinate_t mC;
                        for(int k=0; k<3; ++k) {
                            mC.position[k] = iC.position[k] + segment[k] * (static_cast<double>(j-1) + 0.5);
                        }
                        for(int k=0; k<3; ++k) {
                            res[lastSegment].cell[k] = static_cast<int>(std::floor(mC.position[k]));
                        }
                        res[lastSegment].orientation = interval.getOrientation();
                        lastSegment++;
                    }
                }
            }

            res[0].tag = pl.coordIds[0];
            res[lastSegment-1].tag = pl.coordIds[pl.coordIds.size()-1];
            
            return res;
        }

        pixel_t nodeToPixel(const node_t& node) {
            pixel_t res;
            coordinate_t c;
            bool coordFound = false;

            // Get first coordinate
            c = getCoordinate(node.coordIds[0], &coordFound);
            if (!coordFound) {
                std::cerr << "ERROR: converting node to pixel. Coordinate not found." << std::endl;
                return res;
            }
            for(int k=0; k<3; ++k) {
                res.cell[k] = c.position[k];
            }
            res.tag = node.coordIds[0];
            return res;
        }
    };

} // namespace mesh_m

#endif