#ifndef STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP
#define STAN_MATH_PRIM_PROB_NORMAL_LPDF_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/as_column_vector_or_scalar.hpp>
#include <stan/math/prim/fun/as_array_or_scalar.hpp>
#include <stan/math/prim/fun/as_value_column_array_or_scalar.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/log.hpp>
#include <stan/math/prim/fun/max_size.hpp>
#include <stan/math/prim/fun/size.hpp>
#include <stan/math/prim/fun/size_zero.hpp>
#include <stan/math/prim/fun/to_ref.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <stan/math/prim/functor/operands_and_partials.hpp>
#include <cmath>

namespace stan {
namespace math {

enum class ProbReturnType {Scalar, Vector};

template <typename T, typename = void>
struct prob_broadcaster;

template <typename T>
struct prob_broadcaster<T, require_stan_scalar_t<T>> {
  T ret_;
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  prob_broadcaster(EigArr&& x) : ret_(sum(std::forward<EigArr>(x))) {}

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  prob_broadcaster(Scalar&& x) : ret_(x) {}

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator=(EigArr&& x) {
    ret_ = sum(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator=(Scalar x) {
    ret_ = x;
    return *this;
  }

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator+=(EigArr&& x) {
    ret_ += sum(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator+=(Scalar&& x) {
    ret_ += x;
    return *this;
  }

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator-=(EigArr&& x) {
    ret_ -= sum(x);
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator-=(Scalar&& x) {
    ret_ -= x;
    return *this;
  }
  inline auto ret() noexcept {
    static_assert(!is_var<T>::value, "NOOO");
    return ret_;
  }
  template <typename... Types>
  static auto zero(int /* */) {
    return return_type_t<Types...>(0);
  }

};

template <typename T>
struct prob_broadcaster<T, require_eigen_t<T>> {
  T ret_;
  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  prob_broadcaster(EigArr&& x) : ret_(std::forward<EigArr>(x)) {}

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator=(EigArr&& x) {
    ret_ = x;
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator=(Scalar x) {
    ret_ = Eigen::Array<value_type_t<T>, -1, 1>::Constant(x, ret_.size());
    return *this;
  }

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator+=(EigArr&& x) {
    ret_ += x;
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator+=(Scalar&& x) {
    ret_ += x;
    return *this;
  }

  template <typename EigArr, require_eigen_t<EigArr>* = nullptr>
  inline auto operator-=(EigArr&& x) {
    ret_ -= x;
    return *this;
  }

  template <typename Scalar, require_stan_scalar_t<Scalar>* = nullptr>
  inline auto operator-=(Scalar&& x) {
    ret_ -= x;
    return *this;
  }

  inline auto&& ret() noexcept {
    return std::move(ret_);
  }

  template <typename... Types>
  static auto zero(int size) {
    return Eigen::Array<return_type_t<Types...>, -1, 1>::Constant(0, size).eval();
  }

};



template <ProbReturnType ReturnType, typename... Types>
using prob_return_t = prob_broadcaster<std::conditional_t<ReturnType == ProbReturnType::Scalar, return_type_t<Types...>, Eigen::Array<return_type_t<Types...>, -1, 1>>>;

/** \ingroup prob_dists
 * The log of the normal density for the specified scalar(s) given
 * the specified mean(s) and deviation(s). y, mu, or sigma can
 * each be either a scalar or a vector. Any vector inputs
 * must be the same length.
 *
 * <p>The result log probability is defined to be the sum of the
 * log probabilities for each observation/mean/deviation triple.
 *
 * @tparam T_y type of scalar
 * @tparam T_loc type of location parameter
 * @tparam T_scale type of scale parameter
 * @param y (Sequence of) scalar(s).
 * @param mu (Sequence of) location parameter(s)
 * for the normal distribution.
 * @param sigma (Sequence of) scale parameters for the normal distribution.
 * @return The log of the product of the densities.
 * @throw std::domain_error if the scale is not positive.
 */
template <bool propto, ProbReturnType RetType = ProbReturnType::Scalar, typename T_y, typename T_loc, typename T_scale,
          require_all_not_nonscalar_prim_or_rev_kernel_expression_t<
              T_y, T_loc, T_scale>* = nullptr>
inline auto normal_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  using T_partials_return = partials_return_t<T_y, T_loc, T_scale>;
  using T_y_ref = ref_type_if_t<!is_constant<T_y>::value, T_y>;
  using T_mu_ref = ref_type_if_t<!is_constant<T_loc>::value, T_loc>;
  using T_sigma_ref = ref_type_if_t<!is_constant<T_scale>::value, T_scale>;
  static const char* function = "normal_lpdf";
  check_consistent_sizes(function, "Random variable", y, "Location parameter",
                         mu, "Scale parameter", sigma);
  T_y_ref y_ref = y;
  T_mu_ref mu_ref = mu;
  T_sigma_ref sigma_ref = sigma;

  decltype(auto) y_val = to_ref(as_value_column_array_or_scalar(y_ref));
  decltype(auto) mu_val = to_ref(as_value_column_array_or_scalar(mu_ref));
  decltype(auto) sigma_val = to_ref(as_value_column_array_or_scalar(sigma_ref));

  check_not_nan(function, "Random variable", y_val);
  check_finite(function, "Location parameter", mu_val);
  check_positive(function, "Scale parameter", sigma_val);
  using ret_t = prob_return_t<RetType, T_partials_return>;
  const size_t N = max_size(y, mu, sigma);
  if (size_zero(y, mu, sigma)) {
    return ret_t::template zero<T_y, T_loc, T_scale>(N);
  }
  if (!include_summand<propto, T_y, T_loc, T_scale>::value) {
    return ret_t::template zero<T_y, T_loc, T_scale>(N);
  }

  operands_and_partials<T_y_ref, T_mu_ref, T_sigma_ref> ops_partials(
      y_ref, mu_ref, sigma_ref);

  const auto& inv_sigma
      = to_ref_if<!is_constant_all<T_y, T_scale, T_loc>::value>(inv(sigma_val));
  const auto& y_scaled = to_ref((y_val - mu_val) * inv_sigma);
  const auto& y_scaled_sq
      = to_ref_if<!is_constant_all<T_scale>::value>(y_scaled * y_scaled);

  prob_return_t<RetType, T_partials_return> logp = -0.5 * y_scaled_sq;
  if (include_summand<propto>::value) {
    logp += NEG_LOG_SQRT_TWO_PI * N;
  }
  if (include_summand<propto, T_scale>::value) {
    if (RetType == ProbReturnType::Scalar) {
      logp -= sum(log(sigma_val)) * N / size(sigma);
    } else {
      logp -= log(sigma_val);
    }
  }

  if (!is_constant_all<T_y, T_scale, T_loc>::value) {
    auto scaled_diff = to_ref_if<!is_constant_all<T_y>::value
                                     + !is_constant_all<T_scale>::value
                                     + !is_constant_all<T_loc>::value
                                 >= 2>(inv_sigma * y_scaled);
    if (!is_constant_all<T_y>::value) {
      ops_partials.edge1_.partials_ = -scaled_diff;
    }
    if (!is_constant_all<T_scale>::value) {
      ops_partials.edge3_.partials_ = inv_sigma * y_scaled_sq - inv_sigma;
    }
    if (!is_constant_all<T_loc>::value) {
      ops_partials.edge2_.partials_ = std::move(scaled_diff);
    }
  }
  return ops_partials.build(logp.ret());
}

template <typename T_y, typename T_loc, typename T_scale>
inline auto normal_lpdf(const T_y& y,
                                                      const T_loc& mu,
                                                      const T_scale& sigma) {
  return normal_lpdf<false>(y, mu, sigma);
}

}  // namespace math
}  // namespace stan
#endif
