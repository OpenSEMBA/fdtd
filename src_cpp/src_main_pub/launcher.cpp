#include "semba_fdtd.h"

#include <cstdlib>
#include <iostream>
#include <string>

extern "C" {
    struct semba_fdtd_t;
    semba_fdtd_t* create_semba_fdtd();
    void destroy_semba_fdtd(semba_fdtd_t* p);
    void semba_fdtd_init(semba_fdtd_t* p, const char* flags);
    void semba_fdtd_launch(semba_fdtd_t* p);
    void semba_fdtd_end(semba_fdtd_t* p, const char* case_name);
}

int main(int argc, char** argv) {
    std::string flags;
    for (int i = 1; i < argc; ++i) {
        if (i > 1) flags += ' ';
        flags += argv[i];
    }

    semba_fdtd_t* semba = create_semba_fdtd();
    try {
        semba_fdtd_init(semba, flags.c_str());
        semba_fdtd_launch(semba);
        const std::string input_file = SEMBA_FDTD_m::resolveInputFileFromFlags(flags);
        const std::string case_name = SEMBA_FDTD_m::extractCaseNameFromInput(input_file);
        semba_fdtd_end(semba, case_name.c_str());
    } catch (const std::exception& ex) {
        std::cerr << ex.what() << std::endl;
        destroy_semba_fdtd(semba);
        return 1;
    }
    destroy_semba_fdtd(semba);
    return 0;
}
