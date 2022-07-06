#pragma once
#include "standard_includes.hh"
#include <cmath>

/*
 * Some math functions
 */

namespace df::infra {

// factorial of an unsigned integer
template <typename Tuint>
Tuint factorial(const Tuint n) {
  Tuint lRes = 1;
  for (Tuint i = 2; i <= n; ++i) {
    Tuint lNewRes = lRes * i;
    if (lRes != 0 && lNewRes / lRes != i) {
      std::cout << __PRETTY_FUNCTION__ << ": Warning. Overflow." << std::endl;
    }
    lRes = lNewRes;
  }
  return lRes;
}

// binomial coefficient "n over k" for unsigned integers n, k
// Source: https://www.codewhoop.com/numbers/binomial-coefficient.html
template <typename Tuint>
Tuint binomial(const Tuint n, Tuint k) {
  assert(n >= k);
  Tuint lRes = 1;  // binom(n, 0) = binom(n, n) = 1
  if (k > (n - k)) {
    k = n - k;     // binom(n, k) = binom(n, n-k)  (symmetry)
  }
  for (size_t i = 0; i < k; ++i) {
    lRes *= (n - i);  // numerator: n * (n-1) * ... * (n-k+1)
    lRes /= (i + 1);  // denominator: 1 * 2 * ... * k
  }
  return lRes;
}

// given insigned integer n >= 0 in base 10, return number of digits needed to represent in base b.
// https://www.geeksforgeeks.org/given-number-n-decimal-base-find-number-digits-base-base-b/
template <typename Tuint>
size_t number_of_digits(const Tuint n, const Tuint b = 10) {
  if (0 == n) { return 1; }
  return static_cast<size_t>(std::floor(std::log(n) / std::log(b)) + 1);
}

// given an unsigned integer n >= 0, check if n is a power of b
// https://codereview.stackexchange.com/a/117204
template <typename Tuint>
bool is_power_of(Tuint n, const Tuint b = 10) {
  while (n >= b && n % b == 0) {
    n /= b;
  }
  return (n == 1);
}

}
