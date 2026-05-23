#ifndef NGSPICE_INTERFACE_M_H
#define NGSPICE_INTERFACE_M_H

#include <cstdint>

extern "C" {
void command(const char* input);
void circ(char** input);
void start();
void* get_vector_info(const char* name);
char** get_all_plots();
int32_t has_error();
}

#endif
