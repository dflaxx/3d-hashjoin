#ifndef INFRA_ZIPF_DISTRIBUTION_HH
#define INFRA_ZIPF_DISTRIBUTION_HH

#include <algorithm>
#include <cmath>
#include <random>
#include "standard_includes.hh"

/*
 * Zipf-like random distribution.
 * returns zipf-distributed number in [1,n]
 *
 * "Rejection-inversion to generate variates from monotone discrete
 * distributions", Wolfgang Hörmann and Gerhard Derflinger
 * ACM TOMACS 6.3 (1996): 169-184
 *
 * Implementation by Stack Overflow user drobilla (https://stackoverflow.com/users/1582476/drobilla):
 * https://stackoverflow.com/questions/9983239/how-to-generate-zipf-distributed-numbers-efficiently
 * https://stackoverflow.com/questions/9983239/how-to-generate-zipf-distributed-numbers-efficiently/44154095#44154095
 */

template<class IntType = unsigned long, class RealType = double>
class zipf_distribution
{
public:
    typedef RealType input_type;
    typedef IntType result_type;

    static_assert(std::numeric_limits<IntType>::is_integer, "");
    static_assert(!std::numeric_limits<RealType>::is_integer, "");

    zipf_distribution(const IntType n=std::numeric_limits<IntType>::max(),
                      const RealType q=1.0)
        : n(n)
        , q(q)
        , H_x1(H(1.5) - 1.0)
        , H_n(H(n + 0.5))
        , dist(H_x1, H_n)
        , pmf_denom(0.0)
    {
        // 
        for (IntType i = 1; i <= n; ++i) {
            pmf_denom = pmf_denom + (1.0 / std::pow((double)i, q));
        }
        pmf_denom = 1.0 / pmf_denom;
    }

    IntType operator()(std::mt19937& rng)
    {
        while (true) {
            const RealType u = dist(rng);
            const RealType x = H_inv(u);
            const IntType  k = clamp<IntType>(std::round(x), 1, n);
            if (u >= H(k + 0.5) - h(k)) {
                return k;
            }
        }
    }

    // <- added by dflachs, 04/2021 --
    // probability mass function (pmf) for value k in [1,n] (or, value with rank k)
    RealType pmf(const IntType k)
    {
        assert(1 <= k && k <= n);
        return (1.0 / std::pow(k, q)) * pmf_denom;
    }

    // cumulative distribution function (cdf) for value k in [1,n] (or, value with rank k)
    RealType cdf(const IntType k)
    {
        assert(1 <= k && k <= n);
        RealType res = 0.0;
        for (IntType i = 1; i <= k; ++i) {
            res += pmf(i);
        }
        return res;
    }
    // -->

private:
    /** Clamp x to [min, max]. */
    // static constexpr T clamp(const T x, const T min, const T max)
    template<typename T>
    static const T clamp(const T x, const T min, const T max)
    {
        return std::max(min, std::min(max, x));
    }

    /** exp(x) - 1 / x */
    static double
    expxm1bx(const double x)
    {
        return (std::abs(x) > epsilon)
            ? std::expm1(x) / x
            : (1.0 + x/2.0 * (1.0 + x/3.0 * (1.0 + x/4.0)));
    }

    /** H(x) = log(x) if q == 1, (x^(1-q) - 1)/(1 - q) otherwise.
     * H(x) is an integral of h(x).
     *
     * Note the numerator is one less than in the paper order to work with all
     * positive q.
     */
    const RealType H(const RealType x)
    {
        const RealType log_x = std::log(x);
        return expxm1bx((1.0 - q) * log_x) * log_x;
    }

    /** log(1 + x) / x */
    static RealType
    log1pxbx(const RealType x)
    {
        return (std::abs(x) > epsilon)
            ? std::log1p(x) / x
            : 1.0 - x * ((1/2.0) - x * ((1/3.0) - x * (1/4.0)));
    }

    /** The inverse function of H(x) */
    const RealType H_inv(const RealType x)
    {
        const RealType t = std::max(-1.0, x * (1.0 - q));
        return std::exp(log1pxbx(t) * x);
    }

    /** That hat function h(x) = 1 / (x ^ q) */
    const RealType h(const RealType x)
    {
        return std::exp(-q * std::log(x));
    }

private:
    // static constexpr RealType epsilon = 1e-8;
    static const RealType epsilon;

private:
    IntType                                  n;     ///< Number of elements
    RealType                                 q;     ///< Exponent
    RealType                                 H_x1;  ///< H(x_1)
    RealType                                 H_n;   ///< H(n)
    std::uniform_real_distribution<RealType> dist;  ///< [H(x_1), H(n)]

    // added by dflachs, 04/2021
    RealType                                 pmf_denom; ///< denominator of the probability mass
                                                        //   function (cf. Wikipedia: Zipf's law)

};


template<class IntType, class RealType>
const RealType zipf_distribution<IntType,RealType>::epsilon = 1e-8;

#endif

