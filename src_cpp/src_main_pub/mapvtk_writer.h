#ifndef MAPVTK_WRITER_H
#define MAPVTK_WRITER_H

#include <nlohmann/json.hpp>
#include <string>

namespace mapvtk {

bool flagsContainMapVtk(const std::string& input_flags);

// Match probeOutputPrefix / Fortran fichin stems: "case.fdtd__MAP_..." not "case.fdtd.fdtd__MAP_...".
std::string mapOutputStem(const std::string& case_name);

void writeMapVtkFromJson(const std::string& case_name, const nlohmann::json& root);
void writeCurrentMapVtkFromJson(const std::string& case_name, const nlohmann::json& root);

} // namespace mapvtk

#endif
