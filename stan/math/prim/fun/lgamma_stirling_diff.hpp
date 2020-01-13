#ifndef STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP
#define STAN_MATH_PRIM_SCAL_FUN_LGAMMA_STIRLING_DIFF_HPP

#include <cmath>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/lgamma.hpp>
#include <stan/math/prim/fun/inv.hpp>
#include <stan/math/prim/fun/square.hpp>
#include <stan/math/prim/fun/lgamma_stirling.hpp>
#include <stan/math/prim/err/check_nonnegative.hpp>

namespace stan {
namespace math {

constexpr double lgamma_stirling_diff_useful = 10;

/**
 * Return the difference between log of the gamma function and it's Stirling
 * approximation.
 * This is useful to stably compute log of ratios of gamma functions with large
 * arguments where the Stirling approximation allows for analytic solution
 * and the (small) differences can be added afterwards.
 * This is for example used in the implementation of lbeta.
 *
 * The function will return correct value for all arguments, but the can add
 * precision only when x >= lgamma_stirling_diff_useful.
 *
   \f[
   \mbox{lgamma_stirling_diff}(x) =
    \log(\Gamma(x)) - \frac{1}{2} \log(2\pi) +
    (x-\frac{1}{2})*\log(x) - x
   \f]

 *
 * @param x value
 * @return Difference between lgamma(x) and it's Stirling approximation.
 * @tparam T Type of  value.
 */

template <typename T>
T lgamma_stirling_diff(const T x) {
  static const char* function = "lgamma_stirling_diff";
  check_nonnegative(function, "argument", x);

  if (value_of(x) < lgamma_stirling_diff_useful) {
    return lgamma(x) - lgamma_stirling(x);
  } else {
    // Using the Stirling series as expressed in formula 5.11.1. at
    // https://dlmf.nist.gov/5.11
    constexpr double stirling_series[] = {
        0.0833333333333333333333333,
        -0.00277777777777777777777778,
        0.000793650793650793650793651,
        -0.000595238095238095238095238,
    };
    constexpr int n_stirling_terms = 3;
    T result(0.0);
    T multiplier = inv(x);
    T inv_x_squared = inv(square(x));
    for (int n = 0; n < n_stirling_terms; n++, multiplier *= inv_x_squared) {
      result += stirling_series[n] * multiplier;
    }
    return result;
  }
}

}  // namespace math
}  // namespace stan

#endif
