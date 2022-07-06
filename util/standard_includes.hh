#pragma once

#include <cstdlib>
#include <cassert>
#include <cstdint>

#include <iostream>
#include <iomanip>

#include <vector>

// typedefs (taken from Guido)
using byte_t = unsigned char;

using int_vt   = std::vector<int> ;
using int32_vt = std::vector<int32_t>;

using uint     = unsigned int;
using uint_vt  = std::vector<uint>;
using uint_vvt = std::vector<uint_vt>;

using uint32_vt = std::vector<uint32_t>;

using int64_vt  = std::vector<int64_t>;
using uint64_vt = std::vector<uint64_t>;

using double_vt = std::vector<double>;

using constcharptr_t = const char*;

using bool_vt = std::vector<bool>;

using string_vt = std::vector<std::string>;
