#pragma once

#include <iostream>

template <typename Tattr>
struct tuple_2_t { Tattr k, a; };

template <typename Tattr>
struct tuple_3_t { Tattr k, a, b; };

// [k, a]
using tuple_uint32_2_t = tuple_2_t<uint32_t>;
using tuple_uint64_2_t = tuple_2_t<uint64_t>;

// [k, a, b]
using tuple_uint32_3_t = tuple_3_t<uint32_t>;
using tuple_uint64_3_t = tuple_3_t<uint64_t>;

// output operators
inline std::ostream&
operator<<(std::ostream& os, const tuple_uint32_2_t& t) {
  os << "[" << t.k << "," << t.a << "]";
  return os;
}

inline std::ostream&
operator<<(std::ostream& os, const tuple_uint32_3_t& t) {
  os << "[" << t.k << "," << t.a << "]";
  return os;
}
