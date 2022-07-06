#pragma once

#include <chrono>
#include <cstdint>
#include <stdexcept>
#include <ratio>

#include "standard_includes.hh"
#include "debugging_helpers.hh"
#include "math.hh"

namespace df::infra {

// time units with floating-point tick representation, i.e. with a fractional part
using milliseconds = std::chrono::duration<double, std::milli>;
using microseconds = std::chrono::duration<double, std::micro>;
using seconds      = std::chrono::duration<double>;
using minutes      = std::chrono::duration<double, std::ratio<60>>;
using hours        = std::chrono::duration<double, std::ratio<3600>>;


/*
 * to_string methods for std::ratio and std::chrono::duration
 *
 * TODO: consider re-implementing like this:
 * https://codereview.stackexchange.com/questions/185477/string-conversion-and-stream-insertion-of-stdratio-and-stdchronoduration
 */

/*
 * Get the string representation of an SI prefix stored as a std::ratio,
 * e.g., Ratio<1, 1000> -> m (milli)
 */
template <typename Ratio>
inline std::string to_si_prefix() {
  if (std::ratio_equal_v<Ratio, std::ratio<1, 1>>) { return ""; }

  //if (std::ratio_equal_v<std::ratio<Num, Denom>, std::yocto>) {
  //}
  //if (std::ratio_equal_v<std::ratio<Num, Denom>, std::zepto>) {
  //}
  if constexpr (std::ratio_equal_v<Ratio, std::atto>)  { return "a"; }
  if constexpr (std::ratio_equal_v<Ratio, std::atto>)  { return "a"; }
  if constexpr (std::ratio_equal_v<Ratio, std::femto>) { return "f"; }
  if constexpr (std::ratio_equal_v<Ratio, std::pico>)  { return "p"; }
  if constexpr (std::ratio_equal_v<Ratio, std::nano>)  { return "n"; }
  if constexpr (std::ratio_equal_v<Ratio, std::micro>) { return "u"; }
  if constexpr (std::ratio_equal_v<Ratio, std::milli>) { return "m"; }
  if constexpr (std::ratio_equal_v<Ratio, std::centi>) { return "c"; }
  if constexpr (std::ratio_equal_v<Ratio, std::deci>)  { return "d"; }
  if constexpr (std::ratio_equal_v<Ratio, std::deca>)  { return "da"; }
  if constexpr (std::ratio_equal_v<Ratio, std::hecto>) { return "h"; }
  if constexpr (std::ratio_equal_v<Ratio, std::kilo>)  { return "k"; }
  if constexpr (std::ratio_equal_v<Ratio, std::mega>)  { return "M"; }
  if constexpr (std::ratio_equal_v<Ratio, std::giga>)  { return "G"; }
  if constexpr (std::ratio_equal_v<Ratio, std::tera>)  { return "T"; }
  if constexpr (std::ratio_equal_v<Ratio, std::peta>)  { return "P"; }
  if constexpr (std::ratio_equal_v<Ratio, std::exa>)   { return "E"; }
  //if (std::ratio_equal_v<std::ratio<Num, Denom>, std::zetta>) {
  //}
  //if (std::ratio_equal_v<std::ratio<Num, Denom>, std::yotta>) {
  //}
  throw std::invalid_argument("Ratio must be a standard typedef from std::ratio");
}

/*
 * Get the string representation of a std::duration (time period),
 * including the unit,
 * e.g., std::chrono::milliseconds(4) -> "4ms"
 *
 * This is only needed until there is compiler support for
 * std::chrono::operator<<(std::ostream&, const std::chrono::duration&)
 * (C++20 feature, see
 * https://en.cppreference.com/w/cpp/chrono/duration/operator_ltlt)
 *
 * TODO: only handles duration units <= seconds (i.e., ns, us, ms, s),
 *       but not e.g., hours...
 */
template <typename Rep, typename Period>
inline std::string to_string(const std::chrono::duration<Rep, Period>& aDuration,
                             const bool aSpace = false) {
  std::string res = std::to_string(aDuration.count());
  res += (aSpace ? " " : "") + to_si_prefix<Period>() + "s";
  return res;
}


}
