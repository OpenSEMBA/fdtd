#include <string>
#include <vector>
#include <iostream>
#include <cstdlib>
struct getargs_t {
    int narg = 0;
    std::vector<std::string> arg;
};
void getargs(int, char**, getargs_t&) {}
bool switch_present(const std::string&, const getargs_t&) { return false; }
std::string get_switch_value(const std::string&, const getargs_t&) { return ""; }
