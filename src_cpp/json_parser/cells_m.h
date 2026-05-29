#ifndef CELLS_M_H
#define CELLS_M_H

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

namespace cells_m {

    constexpr int DIR_NULL = 0;
    constexpr int DIR_X    = 1;
    constexpr int DIR_Y    = 2;
    constexpr int DIR_Z    = 3;

    constexpr int NO_TAG   = -1;

    constexpr int CELL_TYPE_PIXEL  = 0;
    constexpr int CELL_TYPE_LINEL  = 1;
    constexpr int CELL_TYPE_SURFEL = 2;
    constexpr int CELL_TYPE_VOXEL  = 3;

    struct cell_t {
        int cell[3];
    };

    struct pixel_t : public cell_t {
        int tag = NO_TAG;
    };

    struct linel_t : public cell_t {
        int orientation = DIR_NULL;
        int tag = NO_TAG;
    };

    struct surfel_t : public cell_t {
        int orientation = DIR_NULL;
        int tag = NO_TAG;
    };

    struct voxel_t : public cell_t {
        std::string tag;
    };

    class cell_interval_t {
    public:
        cell_t ini;
        cell_t end;

        int getType() const;
        int getOrientation() const;
        int getSize() const;

    private:
        int varyingDirections() const;
    };

    class cell_region_t {
    public:
        std::vector<cell_interval_t> intervals;

        std::vector<pixel_t> toPixels() const;
        std::vector<cell_interval_t> getIntervalsOfType(int cellType) const;
    };

    // Forward declarations — implemented inline below
    bool operator==(const linel_t& a, const linel_t& b);
    bool operator==(const pixel_t& a, const pixel_t& b);

    // ========== Inline implementations ==========

    inline int cell_interval_t::varyingDirections() const {
        int res = 0;
        for (int i = DIR_X; i <= DIR_Z; ++i) {
            if ((end.cell[i - 1] - ini.cell[i - 1]) != 0) {
                res += 1;
            }
        }
        return res;
    }

    inline int cell_interval_t::getType() const {
        return varyingDirections();
    }

    inline int cell_interval_t::getOrientation() const {
        int res = DIR_NULL;
        int type = getType();

        if (type == CELL_TYPE_LINEL) {
            int diff;
            for (int i = DIR_X; i <= DIR_Z; ++i) {
                diff = end.cell[i - 1] - ini.cell[i - 1];
                if (diff > 0) { res = i; return res; }
                else if (diff < 0) { res = -i; return res; }
            }
        } else if (type == CELL_TYPE_SURFEL) {
            int diff[3];
            for (int i = 0; i < 3; ++i) {
                diff[i] = end.cell[i] - ini.cell[i];
            }
            int idx = 0;
            for (int i = DIR_X; i <= DIR_Z; ++i) {
                if (diff[i - 1] == 0) { idx = i; }
            }
            res = idx;
            auto fortran_mod = [](int a, int b) -> int {
                int r = a % b;
                if (r < 0) r += b;
                return r;
            };
            int idx1 = fortran_mod(res, 3) + 1;
            int idx2 = fortran_mod(res + 1, 3) + 1;
            if (diff[idx1 - 1] < 0 && diff[idx2 - 1] < 0) {
                res = -res;
            }
        } else {
            res = DIR_NULL;
        }
        return res;
    }

    inline int cell_interval_t::getSize() const {
        int res = 1;
        int diff[3];
        for (int i = 0; i < 3; ++i) {
            diff[i] = std::abs(end.cell[i] - ini.cell[i]);
        }
        for (int i = DIR_X; i <= DIR_Z; ++i) {
            if (diff[i - 1] != 0) {
                res = res * diff[i - 1];
            }
        }
        if (diff[0] == 0 && diff[1] == 0 && diff[2] == 0) {
            res = 0;
        }
        return res;
    }

    inline std::vector<cell_interval_t> cell_region_t::getIntervalsOfType(int cellType) const {
        int count = 0;
        for (const auto& interval : intervals) {
            if (interval.getType() == cellType) {
                count++;
            }
        }
        std::vector<cell_interval_t> res;
        res.reserve(count);
        int j = 0;
        for (const auto& interval : intervals) {
            if (interval.getType() == cellType) {
                res.push_back(interval);
                j++;
            }
        }
        return res;
    }

    inline std::vector<pixel_t> cell_region_t::toPixels() const {
        std::vector<cell_interval_t> pixelIntervals = getIntervalsOfType(CELL_TYPE_PIXEL);
        std::vector<pixel_t> res;
        res.reserve(pixelIntervals.size());
        for (const auto& interval : pixelIntervals) {
            pixel_t p;
            p.cell[0] = interval.ini.cell[0];
            p.cell[1] = interval.ini.cell[1];
            p.cell[2] = interval.ini.cell[2];
            p.tag = NO_TAG;
            res.push_back(p);
        }
        return res;
    }

    inline bool operator==(const pixel_t& a, const pixel_t& b) {
        return (a.cell[0] == b.cell[0]) &&
               (a.cell[1] == b.cell[1]) &&
               (a.cell[2] == b.cell[2]) &&
               (a.tag == b.tag);
    }

    inline bool operator==(const linel_t& a, const linel_t& b) {
        return (a.cell[0] == b.cell[0]) &&
               (a.cell[1] == b.cell[1]) &&
               (a.cell[2] == b.cell[2]) &&
               (a.tag == b.tag) &&
               (a.orientation == b.orientation);
    }

} // namespace cells_m

#endif // CELLS_M_H
