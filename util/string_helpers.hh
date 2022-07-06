#pragma once

#include "standard_includes.hh"
#include <cctype>

namespace df::infra {

// trim from start (in place)
inline void
ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
    return !std::isspace(ch);
  }));
}

// trim from end (in place)
inline void
rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
    return !std::isspace(ch);
  }).base(), s.end());
}

inline std::string&
to_lower(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [] (unsigned char c) { return std::tolower(c); });
  return s;
}

inline std::string&
to_upper(std::string &s) {
  std::transform(s.begin(), s.end(), s.begin(),
                 [] (unsigned char c) { return std::toupper(c); });
  return s;
}

} // namespace df::infra
