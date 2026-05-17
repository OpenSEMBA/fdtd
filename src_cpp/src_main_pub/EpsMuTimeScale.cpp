#include <cmath>
#include <limits>
#include <string>

// Assuming FDETYPES_m provides RKind definition. 
// In a real scenario, this would be an include of the header defining RKind.
// For translation purposes, we assume RKind is a type alias for double.
#ifndef RKind
#define RKind double
#endif

namespace EpsMuTimeScale_m {

    struct EpsMuTimeScale_input_parameters_t {
        RKind tini;
        RKind tend;
        RKind alpha_max;
        bool electric;
        bool magnetic;
        bool are_there;

        // Default constructor to initialize members
        EpsMuTimeScale_input_parameters_t() : tini(0.0), tend(0.0), alpha_max(0.0), 
                                              electric(false), magnetic(false), are_there(false) {}

        // Method corresponding to new_input_
        void init0() {
            new_input_(*this);
        }

        // Method corresponding to get_slope_
        RKind get_slope() const {
            return get_slope_(*this);
        }

        // Method corresponding to checkError_
        int checkError() const {
            return checkError_(*this);
        }
    };

    // Private helper functions

    void new_input_(EpsMuTimeScale_input_parameters_t& this_obj) {
        this_obj.alpha_max = 1.0_RKind;
        this_obj.tini = std::numeric_limits<RKind>::max(); // 1e20_Rkind approximated by max or specific large value
        this_obj.tend = std::numeric_limits<RKind>::max(); // 1e20_Rkind approximated by max or specific large value
        this_obj.electric = false;
        this_obj.magnetic = false;
        this_obj.are_there = false;
    }

    RKind get_slope_(const EpsMuTimeScale_input_parameters_t& this_obj) {
        return (this_obj.alpha_max - 1.0_RKind) / (this_obj.tend - this_obj.tini);
    }

    int checkError_(const EpsMuTimeScale_input_parameters_t& this_obj) {
        int res = 0;
        if (this_obj.alpha_max <= 0.0 || this_obj.tini < 0.0) {
            res = -1;
        }
        return res;
    }

} // namespace EpsMuTimeScale_m