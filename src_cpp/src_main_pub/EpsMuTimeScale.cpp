#include <vector>
#include <string>
#include <iostream>
namespace EpsMuTimeScale_m {
struct EpsMuTimeScale_input_parameters_t {
    int checkError() const { return 0; }
};
void epsMuTimeScale(const std::vector<double>&, const std::vector<double>&, std::vector<double>&, const EpsMuTimeScale_input_parameters_t&) {}
}
