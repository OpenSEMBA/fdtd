#include "semba_fdtd.h"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#ifdef CompileWithMPI
#include <mpi.h>
#endif

extern "C" {
    struct semba_fdtd_t;
    semba_fdtd_t* create_semba_fdtd();
    void destroy_semba_fdtd(semba_fdtd_t* p);
    void semba_fdtd_init(semba_fdtd_t* p, const char* flags);
    void semba_fdtd_launch(semba_fdtd_t* p);
    void semba_fdtd_end(semba_fdtd_t* p, const char* case_name);
}

namespace {

void write_stderr_line(const char* text) {
    if (text == nullptr) return;
    std::fputs(text, stderr);
    std::fputc('\n', stderr);
    std::fflush(stderr);
}

#ifdef CompileWithMPI
[[noreturn]] void abort_mpi_world(int code) {
    MPI_Abort(MPI_COMM_WORLD, code);
    std::abort();
}
#endif

} // namespace

int main(int argc, char** argv) {
#ifdef CompileWithMPI
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) {
        MPI_Init(&argc, &argv);
    }
#endif

    if (argc > 0 && argv[0] != nullptr) {
        setenv("SEMBA_FDTD_BINARY_PATH", argv[0], 1);
    }

    std::string flags;
    for (int i = 1; i < argc; ++i) {
        if (i > 1) flags += ' ';
        flags += argv[i];
    }

    semba_fdtd_t* semba = create_semba_fdtd();
    try {
        semba_fdtd_init(semba, flags.c_str());
        semba_fdtd_launch(semba);
        semba_fdtd_end(semba, "");
    } catch (const std::exception& ex) {
        if (std::string(ex.what()) == "__SEMBA_FDTD_STOP_1__") {
            write_stderr_line("STOP 1");
        } else {
            write_stderr_line(ex.what());
        }
        destroy_semba_fdtd(semba);
#ifdef CompileWithMPI
        int mpi_finalized = 0;
        int mpi_initialized_now = 0;
        MPI_Initialized(&mpi_initialized_now);
        MPI_Finalized(&mpi_finalized);
        if (mpi_initialized_now && !mpi_finalized) {
            abort_mpi_world(1);
        }
#endif
        return 1;
    } catch (...) {
        write_stderr_line("Unhandled non-standard exception");
        destroy_semba_fdtd(semba);
#ifdef CompileWithMPI
        int mpi_finalized = 0;
        int mpi_initialized_now = 0;
        MPI_Initialized(&mpi_initialized_now);
        MPI_Finalized(&mpi_finalized);
        if (mpi_initialized_now && !mpi_finalized) {
            abort_mpi_world(1);
        }
#endif
        return 1;
    }
    destroy_semba_fdtd(semba);
#ifdef CompileWithMPI
    int mpi_finalized = 0;
    MPI_Finalized(&mpi_finalized);
    if (!mpi_finalized && !mpi_initialized) {
        MPI_Finalize();
    }
#endif
    return 0;
}
