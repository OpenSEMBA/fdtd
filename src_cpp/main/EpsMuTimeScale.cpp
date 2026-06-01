#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

// Assuming FDETYPES_m defines RKind as double
// If not, we define it here for completeness
#ifndef RKind
#define RKind double
#endif

namespace EpsMuTimeScale_m {

    class EpsMuTimeScale_input_parameters_t {
    public:
        double tini;
        double tend;
        double alpha_max;
        bool electric;
        bool magnetic;
        bool are_there;

        EpsMuTimeScale_input_parameters_t() {
            init0();
        }

        double get_slope() const {
            return (alpha_max - 1.0) / (tend - tini);
        }

        void init0() {
            alpha_max = 1.0;
            tini = 1e20;
            tend = 1e20;
            electric = false;
            magnetic = false;
            are_there = false;
        }

        int checkError() const {
            int res = 0;
            if (alpha_max <= 0.0 || tini < 0.0) {
                res = -1;
            }
            return res;
        }
    };

} // namespace EpsMuTimeScale_m