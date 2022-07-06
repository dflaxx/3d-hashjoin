#include "standard_includes.hh"

#include <functional>

#include "chrono_helpers.hh"

namespace df::infra {

/*
 * Repeat a function aFunc until its cumulative execution time is at least aMinTime.
 * Return the cumulative execution time and the number of repetitions.
 * aTeardown is a std::function that is executed after each repetition
 * to prepare for the next, e.g. to clear data structures.
 */
inline
std::pair<std::chrono::nanoseconds, size_t>
repeat_mintime(const std::chrono::nanoseconds aMinTime,
                  const std::function<void(void)> aFunc,
                  const std::function<void(void)> aTeardown = [](){},
                  const bool aTeardownAfterLast = false,
                  const size_t aMinRepeat = 1) {
  size_t n = aMinRepeat;
  using namespace std::chrono_literals;
  std::chrono::time_point<std::chrono::steady_clock> t0, t1;
  std::chrono::nanoseconds lTotalTime = 0s;

  for (size_t i = 0; i < n; ++i) {
    t0 = std::chrono::steady_clock::now();
    aFunc();
    t1 = std::chrono::steady_clock::now();
    lTotalTime += (t1 - t0);

    if (i == n-1 && lTotalTime < aMinTime) {
      n *= 2;
    }

    if (i != n-1 || aTeardownAfterLast) { aTeardown(); }
  }

  return {lTotalTime, n};
}

};
