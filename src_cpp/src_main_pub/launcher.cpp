#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <cstdlib>
#include <algorithm>
#include "semba_fdtd.cpp"

int main(int argc, char* argv[]) {
    std::string input_file = "input.fdtd.json";
    for (int i = 1; i < argc; i++) {
        if (std::string(argv[i]) == "-i" && i + 1 < argc) {
            input_file = argv[i+1];
            i++;
        }
    }
    
    // Extract case name from filename
    std::string caseName = input_file;
    size_t pos = caseName.find_last_of("/\\");
    if (pos != std::string::npos) caseName = caseName.substr(pos + 1);
    pos = caseName.find(".fdtd");
    if (pos != std::string::npos) caseName = caseName.substr(0, pos);
    else {
        pos = caseName.find(".json");
        if (pos != std::string::npos) caseName = caseName.substr(0, pos);
    }
    
    std::cout << "semba-fdtd C++ translator (first iteration)" << std::endl;
    std::cout << "Input: " << input_file << std::endl;
    std::cout << "WARNING: This is a direct Fortran-to-C++ translation." << std::endl;
    std::cout << "Many functions are placeholder stubs requiring manual implementation." << std::endl;
    
    SEMBA_FDTD_m::semba_fdtd_t semba;
    semba.init(input_file);
    semba.launch();
    semba.end(caseName);
    
    return 0;
}
