#ifndef REPORT_M_H
#define REPORT_M_H

#include <string>
#include <iostream>
#include <stdexcept>

namespace Report_m {

    // WarnErrReport: prints message to stderr, throws on fatal
    inline void WarnErrReport(const std::string& msg, bool fatal) {
        std::cerr << "Error: " << msg << std::endl;
        if (fatal) {
            throw std::runtime_error("Fatal error in smbjson: " + msg);
        }
    }

} // namespace Report_m

#endif // REPORT_M_H
