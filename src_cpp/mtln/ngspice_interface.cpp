#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <functional>

// Forward declarations for C bindings if needed, 
// but since these are external C functions, we declare them here.
// Note: In a real scenario, you might include the actual C header file.

namespace ngspice_interface_m {

    // Corresponds to Fortran type, bind(c) :: vectorInfo_t
    struct vectorInfo_t {
        void* vName;          // type(c_ptr)
        int32_t vType;        // integer(c_int)
        int16_t vFlags;       // integer(c_short)
        void* vRealData;      // type(c_ptr) !real
        void* vCompData;      // type(c_ptr) !ngcomplex
        int32_t vLength;      // integer(c_int)
    };

    // Interface declarations for external C functions

    // subroutine command(input) bind (C, name = "command")
    // character(kind=c_char), dimension(*), intent(in) :: input
    // We assume input is a null-terminated string or a pointer to char array.
    // Since Fortran character(*) with intent(in) usually maps to a C string,
    // we use const char*.
    extern "C" void command(const char* input);

    // subroutine circ(input) bind(C, name = "circ")
    // type(c_ptr), intent(in) :: input(*)
    // This is an array of pointers. In C++, this is typically void** or const void**.
    // The size is not specified in the interface, so we might need a count or assume null-terminated array of pointers.
    // However, without a count, it's hard to translate directly to a safe C++ function.
    // Assuming it's an array of c_ptr (void*) terminated by a null pointer or passed with a count in a real scenario.
    // For strict translation, we'll use a vector of void* or raw pointer array.
    // Let's assume it takes an array of void* pointers.
    extern "C" void circ(void** input);

    // subroutine start() bind (C, name = "start")
    extern "C" void start();

    // type(c_ptr) function get_vector_info(name) bind (C, name="get_vector_info")
    // character(kind=c_char), dimension(*), intent(in) :: name
    extern "C" void* get_vector_info(const char* name);

    // type(c_ptr) function get_all_plots() bind (C, name="get_all_plots")
    extern "C" void* get_all_plots();

    // integer(c_int) function has_error() bind (C, name="has_error")
    extern "C" int32_t has_error();

} // namespace ngspice_interface_m