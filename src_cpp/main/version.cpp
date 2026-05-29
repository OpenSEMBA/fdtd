// version_m.h
// This file contains the version information module converted to a C++ namespace.

#ifndef VERSION_M_H
#define VERSION_M_H

#include <string>

namespace version_m {
    constexpr const char* program_name = "semba-fdtd";
    constexpr const char* compilation_date = __DATE__ " " __TIME__;
    constexpr const char* git_commit = "158c6f3c";
    constexpr const char* compiler_id = "GNU 13.3.0";
    constexpr const char* cmake_build_type = "";
    constexpr const char* compilation_flags = "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch";
    constexpr const char* compilation_flags_debug = "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch -g -O0 -fno-inline -fcheck=all -fbacktrace";
    constexpr const char* compilation_flags_release = "-fopenmp -ffree-form -ffree-line-length-none -fdec -fallow-argument-mismatch -Ofast";
}

#endif // VERSION_M_H