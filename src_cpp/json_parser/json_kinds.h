#ifndef JSON_KINDS_H
#define JSON_KINDS_H

#include <cstdint>

// json-fortran kind definitions (from external/json-fortran/src/json_kinds.F90)
// Default build: REAL64 (8 bytes), INT32 (4 bytes)

constexpr int json_RK = 8;    // real64
constexpr int json_IK = 4;    // int32
constexpr int json_CDK = 1;   // DEFAULT character kind
constexpr int json_LK = 4;    // default logical kind
constexpr int json_CK = 1;    // character kind (same as CDK when not using UCS4)
constexpr int real32 = 4;     // 4-byte real

#endif // JSON_KINDS_H
