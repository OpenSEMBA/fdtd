#ifndef PREPROCESS_TAGS_H
#define PREPROCESS_TAGS_H

#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>

#include "nfde_types.h"

namespace Preprocess_m {

struct tagtype_t {
    int32_t numertags = 0;
    std::vector<std::string> tag;
};

inline std::string trim_copy(const std::string& s) {
    const auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    const auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

inline int32_t searchtag(const tagtype_t& tagtype, const std::string& tag) {
    const std::string needle = trim_copy(tag);
    for (int32_t i = 0; i < tagtype.numertags; ++i) {
        if (trim_copy(tagtype.tag[static_cast<size_t>(i)]) == needle) {
            return i + 1;
        }
    }
    return -1;
}

inline std::string dielectric_tag_at(const NFDETypes_m::Dielectric_t& component,
                                     const std::string& coord_type,
                                     int32_t idx) {
    const size_t i = static_cast<size_t>(idx - 1);
    if (coord_type == "c1P") {
        return trim_copy(component.c1P.at(i).tag);
    }
    return trim_copy(component.c2P.at(i).tag);
}

inline void checkDielectricTagForDuplicate(const NFDETypes_m::Dielectric_t& component,
                                           const std::vector<NFDETypes_m::Dielectric_t>& prev_components,
                                           int32_t n_prev,
                                           int32_t idx,
                                           const std::string& coord_type,
                                           int32_t& numertag,
                                           tagtype_t& tagtype,
                                           int32_t precounting,
                                           const std::string& error_msg) {
    const std::string tagToCheck = dielectric_tag_at(component, coord_type, idx);
    bool foundDuplicate = false;

    if (idx > 1) {
        for (int32_t k = 1; k < idx; ++k) {
            if (tagToCheck == dielectric_tag_at(component, coord_type, k)) {
                foundDuplicate = true;
                break;
            }
        }
    }

    if (!foundDuplicate) {
        for (int32_t m = 0; m < n_prev; ++m) {
            const auto& prev = prev_components[static_cast<size_t>(m)];
            for (int32_t k = 1; k <= prev.n_C1P; ++k) {
                if (tagToCheck == dielectric_tag_at(prev, "c1P", k)) {
                    throw std::runtime_error(error_msg + " Duplicate tag found: " + tagToCheck);
                }
            }
            for (int32_t k = 1; k <= prev.n_C2P; ++k) {
                if (tagToCheck == dielectric_tag_at(prev, "c2P", k)) {
                    throw std::runtime_error(error_msg + " Duplicate tag found: " + tagToCheck);
                }
            }
        }
    }

    if (foundDuplicate) {
        --numertag;
    } else if (precounting == 1) {
        if (static_cast<size_t>(numertag) > tagtype.tag.size()) {
            tagtype.tag.resize(static_cast<size_t>(numertag));
        }
        tagtype.tag[static_cast<size_t>(numertag - 1)] = tagToCheck;
    }
}

inline void checkLossyTagForDuplicate(const NFDETypes_m::LossyThinSurface_t& component,
                                      const std::vector<NFDETypes_m::LossyThinSurface_t>& prev_components,
                                      int32_t n_prev,
                                      int32_t idx,
                                      int32_t& numertag,
                                      tagtype_t& tagtype,
                                      int32_t precounting) {
    const std::string tagToCheck = trim_copy(component.c[static_cast<size_t>(idx - 1)].tag);
    bool foundDuplicate = false;

    if (idx > 1) {
        for (int32_t k = 1; k < idx; ++k) {
            if (tagToCheck == trim_copy(component.c[static_cast<size_t>(k - 1)].tag)) {
                foundDuplicate = true;
                break;
            }
        }
    }

    if (!foundDuplicate && n_prev > 0) {
        for (int32_t m = 0; m < n_prev; ++m) {
            const auto& prev = prev_components[static_cast<size_t>(m)];
            for (int32_t k = 1; k <= prev.nc; ++k) {
                if (tagToCheck == trim_copy(prev.c[static_cast<size_t>(k - 1)].tag)) {
                    foundDuplicate = true;
                }
            }
        }
    }

    if (foundDuplicate) {
        --numertag;
    } else if (precounting == 1) {
        if (static_cast<size_t>(numertag) > tagtype.tag.size()) {
            tagtype.tag.resize(static_cast<size_t>(numertag));
        }
        tagtype.tag[static_cast<size_t>(numertag - 1)] = tagToCheck;
    }
}

} // namespace Preprocess_m

#endif
