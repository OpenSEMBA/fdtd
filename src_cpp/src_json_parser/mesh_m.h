#ifndef MESH_M_H
#define MESH_M_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <memory>
#include <unordered_map>
#include <iostream>

#include "cells_m.h"
#include "conformal_types.h"

namespace mesh_m {

    // Constants (mesh uses 0-based directions)
    constexpr int DIR_X = 0;
    constexpr int DIR_Y = 1;
    constexpr int DIR_Z = 2;

    struct element_t {
        std::vector<int> coordIds;
    };

    struct node_t : public element_t {
        // coordIds must be size 1
    };

    struct polyline_t : public element_t {
        // coordIds must be size > 1
    };

    struct coordinate_t {
        double position[3];

        coordinate_t operator-(const coordinate_t& b) const {
            coordinate_t res;
            for (int i = 0; i < 3; ++i) res.position[i] = this->position[i] - b.position[i];
            return res;
        }

        bool operator==(const coordinate_t& b) const {
            for (int i = 0; i < 3; ++i) {
                if (this->position[i] != b.position[i]) return false;
            }
            return true;
        }
    };

    struct pixel_t {
        double cell[3];
        int tag;
    };

    struct linel_t {
        int tag;
        double cell[3];
        int orientation;
    };

    struct conformal_region_t {
        std::vector<conformal_types_m::triangle_t> triangles;
        std::vector<cells_m::cell_interval_t> intervals;
        int type;
    };

    // Polymorphic base for element storage
    struct mesh_item_base {
        virtual ~mesh_item_base() = default;
        virtual mesh_item_base* clone() const = 0;
    };

    template <typename T>
    struct mesh_item : public mesh_item_base {
        T value;
        mesh_item(const T& v) : value(v) {}
        mesh_item_base* clone() const override { return new mesh_item<T>(value); }
    };

    class mesh_t {
    private:
        std::unordered_map<int, coordinate_t> coordinates;
        std::unordered_map<int, std::shared_ptr<mesh_item_base>> elements;

    public:
        inline void allocateCoordinates(int) {}
        inline void allocateElements(int) {}

        inline void addCoordinate(int id, const coordinate_t& coord) {
            coordinates[id] = coord;
        }

        coordinate_t getCoordinate(int id, bool& found) {
            coordinate_t res{};
            found = false;
            auto it = coordinates.find(id);
            if (it != coordinates.end()) {
                res = it->second;
                found = true;
            }
            return res;
        }

        inline void checkCoordinateId(int id, int& stat) {
            stat = (coordinates.find(id) == coordinates.end()) ? 1 : 0;
        }

        void addElement(int id, const element_t& e) {
            elements[id] = std::make_shared<mesh_item<element_t>>(e);
        }

        void addElement(int id, const node_t& e) {
            elements[id] = std::make_shared<mesh_item<node_t>>(e);
        }

        void addElement(int id, const polyline_t& e) {
            elements[id] = std::make_shared<mesh_item<polyline_t>>(e);
        }

        void addCellRegion(int id, const cells_m::cell_region_t& e) {
            elements[id] = std::make_shared<mesh_item<cells_m::cell_region_t>>(e);
        }

        void addConformalRegion(int id, const conformal_region_t& e) {
            elements[id] = std::make_shared<mesh_item<conformal_region_t>>(e);
        }

        node_t getNode(int id, bool& found) {
            node_t res{};
            found = false;
            auto it = elements.find(id);
            if (it != elements.end()) {
                auto* item = dynamic_cast<mesh_item<node_t>*>(it->second.get());
                if (item) {
                    res = item->value;
                    found = true;
                }
            }
            return res;
        }

        polyline_t getPolyline(int id, bool& found) {
            polyline_t res{};
            found = false;
            auto it = elements.find(id);
            if (it != elements.end()) {
                auto* item = dynamic_cast<mesh_item<polyline_t>*>(it->second.get());
                if (item) {
                    res = item->value;
                    found = true;
                }
            }
            return res;
        }

        cells_m::cell_region_t getCellRegion(int id, bool& found) {
            cells_m::cell_region_t res;
            found = false;
            auto it = elements.find(id);
            if (it != elements.end()) {
                auto* item = dynamic_cast<mesh_item<cells_m::cell_region_t>*>(it->second.get());
                if (item) {
                    res = item->value;
                    found = true;
                } else {
                    auto* item2 = dynamic_cast<mesh_item<conformal_region_t>*>(it->second.get());
                    if (item2 && !item2->value.intervals.empty()) {
                        res.intervals = item2->value.intervals;
                        found = true;
                    }
                }
            }
            return res;
        }

        conformal_region_t getConformalRegion(int id, bool& found) {
            conformal_region_t res;
            found = false;
            auto it = elements.find(id);
            if (it != elements.end()) {
                auto* item = dynamic_cast<mesh_item<conformal_region_t>*>(it->second.get());
                if (item) {
                    res = item->value;
                    found = true;
                }
            }
            return res;
        }

        std::vector<cells_m::cell_region_t> getCellRegions(const std::vector<int>& ids) {
            std::vector<cells_m::cell_region_t> res;
            int count = 0;
            for (int id : ids) { bool found = false; auto cr = getCellRegion(id, found); if (found) count++; }
            res.resize(count);
            int j = 0;
            for (int id : ids) { bool found = false; auto cr = getCellRegion(id, found); if (found) { res[j++] = cr; } }
            return res;
        }

        std::vector<conformal_region_t> getConformalRegions(const std::vector<int>& ids) {
            std::vector<conformal_region_t> res;
            int count = 0;
            for (int id : ids) { bool found = false; auto cr = getConformalRegion(id, found); if (found) count++; }
            res.resize(count);
            int j = 0;
            for (int id : ids) { bool found = false; auto cr = getConformalRegion(id, found); if (found) { res[j++] = cr; } }
            return res;
        }

        int countPolylineSegments(const polyline_t& pl) {
            int res = 0;
            if (pl.coordIds.size() <= 1) return 0;
            for (size_t i = 0; i < pl.coordIds.size() - 1; ++i) {
                bool f1=false, f2=false;
                coordinate_t iC = getCoordinate(pl.coordIds[i], f1);
                coordinate_t eC = getCoordinate(pl.coordIds[i+1], f2);
                cells_m::cell_interval_t interval;
                for(int k=0;k<3;++k) { interval.ini.cell[k]=(int)std::floor(iC.position[k]); interval.end.cell[k]=(int)std::floor(eC.position[k]); }
                res += interval.getSize();
            }
            return res;
        }

        bool arePolylineSegmentsStructured(const polyline_t& pl) {
            if (pl.coordIds.size() <= 1) return false;
            for (size_t i = 0; i < pl.coordIds.size() - 1; ++i) {
                bool f1=false, f2=false;
                coordinate_t iC = getCoordinate(pl.coordIds[i], f1);
                coordinate_t eC = getCoordinate(pl.coordIds[i+1], f2);
                if (!f1 || !f2) return false;
                for(int k=0;k<3;++k) {
                    if (std::floor(iC.position[k]) != iC.position[k]) return false;
                    if (std::floor(eC.position[k]) != eC.position[k]) return false;
                }
                int nvd = 0;
                for (int d=0;d<3;++d) if (iC.position[d] != eC.position[d]) nvd++;
                if (nvd > 1) return false;
            }
            return true;
        }

        std::vector<linel_t> polylineToLinels(const polyline_t& pl) {
            if (!arePolylineSegmentsStructured(pl)) return std::vector<linel_t>();
            int count = countPolylineSegments(pl);
            std::vector<linel_t> res(count);
            if (count == 0) return res;
            int lastSegment = 0;
            for (size_t i = 0; i < pl.coordIds.size() - 1; ++i) {
                bool f1=false, f2=false;
                coordinate_t iC = getCoordinate(pl.coordIds[i], f1);
                coordinate_t eC = getCoordinate(pl.coordIds[i+1], f2);
                cells_m::cell_interval_t interval;
                for(int k=0;k<3;++k) { interval.ini.cell[k]=(int)std::floor(iC.position[k]); interval.end.cell[k]=(int)std::floor(eC.position[k]); }
                if (iC.position[0]!=eC.position[0] || iC.position[1]!=eC.position[1] || iC.position[2]!=eC.position[2]) {
                    int segment[3];
                    int size=interval.getSize();
                    for(int k=0;k<3;++k) segment[k]=(interval.end.cell[k]-interval.ini.cell[k])/size;
                    res[lastSegment].tag=pl.coordIds[i];
                    for(int j=1;j<=size;++j) {
                        coordinate_t mC;
                        for(int k=0;k<3;++k) mC.position[k]=iC.position[k]+segment[k]*(double(j-1)+0.5);
                        res[lastSegment].cell[0]=(int)std::floor(mC.position[0]);
                        res[lastSegment].cell[1]=(int)std::floor(mC.position[1]);
                        res[lastSegment].cell[2]=(int)std::floor(mC.position[2]);
                        res[lastSegment].orientation=cells_m::DIR_NULL;
                        int diff;
                        for(int d=0;d<3;++d) {
                            diff=(int)eC.position[d]-(int)iC.position[d];
                            if(diff>0){res[lastSegment].orientation=d+1;break;}
                            if(diff<0){res[lastSegment].orientation=-(d+1);break;}
                        }
                        lastSegment++;
                    }
                }
            }
            if (!res.empty()) {
                res[0].tag=pl.coordIds[0];
                res[lastSegment-1].tag=pl.coordIds[pl.coordIds.size()-1];
            }
            return res;
        }

        pixel_t nodeToPixel(const node_t& node) {
            pixel_t res{};
            bool cf=false;
            coordinate_t c=getCoordinate(node.coordIds[0], cf);
            if(!cf) return res;
            res.cell[0]=c.position[0]; res.cell[1]=c.position[1]; res.cell[2]=c.position[2];
            res.tag=node.coordIds[0];
            return res;
        }

        inline void printHashInfo() {
            // No-op
        }
    };

} // namespace mesh_m

#endif // MESH_M_H
