#ifndef STAN_MATH_PRIM_PROB_UNIFORM_LCCDF_HPP
#define STAN_MATH_PRIM_PROB_UNIFORM_LCCDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/promote_scalar.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

template <typename T_y, typename T_low, typename T_high>
return_type_t<T_y, T_low, T_high> uniform_lccdf(const T_y& y,
                                                const T_low& alpha,
                                                const T_high& beta) {
  using T_partials_return = partials_return_t<T_y, T_low, T_high>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_alpha_ref = ref_type_if_t<!is_constant<T_low>::value, T_low>;
  using T_beta_ref = ref_type_if_t<!is_constant<T_high>::value, T_high>;
  static const char* function = "uniform_lccdf";
  check_consistent_sizes(function, "Random variable", y,
                         "Lower bound parameter", alpha,
                         "Upper bound parameter", beta);

  T_y_ref y_ref = y;
  T_alpha_ref alpha_ref = alpha;
  T_beta_ref beta_ref = beta;

  const auto& y_col = as_column_vector_or_scalar(y_ref);
  const auto& alpha_col = as_column_vector_or_scalar(alpha_ref);
  const auto& beta_col = as_column_vector_or_scalar(beta_ref);

  const auto& y_arr = as_array_or_scalar(y_col);
  const auto& alpha_arr = as_array_or_scalar(alpha_col);
  const auto& beta_arr = as_array_or_scalar(beta_col);

  ref_type_t<decltype(value_of(y_arr))> y_val = value_of(y_arr);
  ref_type_t<decltype(value_of(alpha_arr))> alpha_val = value_of(alpha_arr);
  ref_type_t<decltype(value_of(beta_arr))> beta_val = value_of(beta_arr);

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Lower bound parameter", alpha_val);
  check_finite(function, "Upper bound parameter", beta_val);
  check_greater(function, "Upper bound parameter", beta_val, alpha_val);

  if (size_zero(y, alpha, beta)) {
    return 0;
  }

  if (sum(promote_scalar<int>(y_val < alpha_val))
      || sum(promote_scalar<int>(beta_val < y_val))) {
    return 0;
  }

  auto ops_partials = operands_and_partials(y_ref, alpha_ref, beta_ref);

  const auto& b_minus_a
      = to_ref_if<!is_constant_all<T_y, T_low, T_high>::value>(beta_val
                                                               - alpha_val);
  const auto& ccdf_log_n
      = to_ref_if<!is_constant_all<T_y, T_low, T_high>::value>(
          1 - (y_val - alpha_val) / b_minus_a);

  T_partials_return ccdf_log = sum(log(ccdf_log_n));

  if (!is_constant_all<T_y>::value) {
    edge<0>(ops_partials).partials_ = inv(-b_minus_a * ccdf_log_n);
  }
  if (!is_constant_all<T_low, T_high>::value) {
    const auto& rep_deriv = to_ref_if<(!is_constant_all<T_low>::value
                                       && !is_constant_all<T_high>::value)>(
        inv(b_minus_a * b_minus_a * ccdf_log_n));
    if (!is_constant_all<T_low>::value) {
      edge<1>(ops_partials).partials_ = (beta_val - y_val) * rep_deriv;
    }
    if (!is_constant_all<T_high>::value) {
      edge<2>(ops_partials).partials_ = (y_val - alpha_val) * rep_deriv;
    }
  }
  return ops_partials.build(ccdf_log);
}

}  // namespace math
}  // namespace stan
#endif
