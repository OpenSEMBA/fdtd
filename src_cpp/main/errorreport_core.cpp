// Minimal Report/error helpers for semba-fdtd-cpp (Fortran: errorreport.F90).
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

#ifndef BUFSIZE
#define BUFSIZE 256
#endif

struct SGGFDTDINFO_t {
    int NumMedia = 0;
};

struct tiempo_t {
    double segundos = 0.0;
    char hora[BUFSIZE] = {};
    char fecha[BUFSIZE] = {};
};

const std::string SEPARADOR = "========================================";
const std::string separador = "========================================";

void openclosedelete(const std::string& ficherin) {
    int my_iostat = 0;
    int retries = 0;
retry_open:
    if (my_iostat != 0) {
        std::cout << '.' << std::flush;
    }
    if (retries++ > 100) {
        return;
    }

    std::string filename = ficherin;
    const size_t start = filename.find_first_not_of(' ');
    if (start == std::string::npos) {
        filename.clear();
    } else {
        filename = filename.substr(start);
    }

    {
        std::ofstream file(filename, std::ios::out);
        if (!file.is_open()) {
            my_iostat = 1;
            goto retry_open;
        }
        file << "!END" << std::endl;
        file.close();
    }
    std::remove(filename.c_str());
}

void openclose(const std::string& ficherin) {
    int my_iostat = 0;
retry_open:
    if (my_iostat != 0) {
        std::cout << '.' << std::flush;
    }

    std::string filename = ficherin;
    const size_t start = filename.find_first_not_of(' ');
    if (start == std::string::npos) {
        filename.clear();
    } else {
        filename = filename.substr(start);
    }

    std::ofstream file(filename, std::ios::out);
    if (!file.is_open()) {
        my_iostat = 1;
        goto retry_open;
    }
    file << "!END" << std::endl;
    file.close();
}

void erasesignalingfiles(bool simu_devia) {
    if (simu_devia) {
        return;
    }
    const char* names[] = {
        "stop", "stopflushing", "flush", "flushdata", "unpack", "postprocess", "flushxdmf", "flushvtk", "snap",
        "stop_only", "stopflushing_only", "flush_only", "flushdata_only",
        "stop_dontwritevtk", "stop_only_dontwritevtk", "stopflushing_dontwritevtk",
        "stopflushing_only_dontwritevtk", "flush_dontwritevtk", "flush_only_dontwritevtk",
        "unpack", "postprocess", "flushxdmf", "flushvtk",
    };
    for (const char* name : names) {
        openclosedelete(name);
    }
}

void WarnErrReport(const std::string& msg, bool error = false) {
    std::cerr << (error ? "Error: " : "Warning: ") << msg << std::endl;
    if (error) {
        std::exit(1);
    }
}

void reportmedia(SGGFDTDINFO_t&) {}

namespace Report_m {

bool printea = false;

void print11(int layoutnumber, const std::string& message, bool force) {
    if (printea || force) {
        std::cout << "[Rank " << layoutnumber << "] " << message << std::endl;
    }
}

void stoponerror(int layoutnumber, int num_procs, const std::string& message, bool /*fatal*/) {
    char whoami_buf[BUFSIZE];
    std::snprintf(whoami_buf, BUFSIZE, "(%d/%d) ", layoutnumber + 1, num_procs);
    print11(layoutnumber, std::string(whoami_buf) + " ERROR: " + message, true);
    openclosedelete("running");
    openclosedelete("pause");
    openclosedelete("relaunch");
    std::exit(1);
}

void StopOnError(int layoutnumber, int num_procs, const std::string& message, bool calledfrommain) {
    stoponerror(layoutnumber, num_procs, message, calledfrommain);
}

void get_secnds(tiempo_t* t) {
    if (t == nullptr) {
        return;
    }
    const std::time_t now = std::time(nullptr);
    t->segundos = static_cast<double>(now);
    if (std::tm* ltime = std::localtime(&now)) {
        std::strftime(t->hora, BUFSIZE, "%H:%M:%S", ltime);
        std::strftime(t->fecha, BUFSIZE, "%Y-%m-%d", ltime);
    }
}

void get_secnds(void* time_out2) {
    get_secnds(static_cast<tiempo_t*>(time_out2));
}

void print_basic_help(void*) {}
void print_help(void*) {}
void print_credits(void*) {}
void removeintraspaces(std::string&) {}
void buscaswitchficheroinput(void*) {}
void default_flags(void*) {}

void erasesignalingfiles() {
    ::erasesignalingfiles(false);
}

void openclosedelete_global(const std::string& filename) {
    openclosedelete(filename);
}

} // namespace Report_m
