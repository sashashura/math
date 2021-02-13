#ifndef STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP
#define STAN_MATH_PRIM_FUN_UB_CONSTRAIN_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/exp.hpp>
#include <stan/math/prim/fun/identity_constrain.hpp>
#include <stan/math/prim/fun/identity_free.hpp>
#include <stan/math/prim/fun/subtract.hpp>
#include <stan/math/prim/fun/sum.hpp>
#include <stan/math/prim/fun/value_of.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * Return the upper-bounded value for the specified unconstrained
 * matrix and upper bound.
 *
 * <p>The transform is
 *
 * <p>\f$f(x) = U - \exp(x)\f$
 *
 * <p>where \f$U\f$ is the upper bound.
 *
 * @tparam T type of Matrix
 * @tparam U type of upper bound
 * @param[in] x free Matrix.
 * @param[in] ub upper bound
 * @return matrix constrained to have upper bound
 */
 template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr, require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub) {
  if (is_positive_infinity(ub)) {
    return identity_constrain(x, ub);
  } else {
    //check_greater("ub_constrain", "ub", value_of(x), value_of(ub));
    return subtract(ub, exp(x));
  }
}

/**
 * Return the upper-bounded value for the specified unconstrained
 * scalar and upper bound and increment the specified log
 * probability reference with the log absolute Jacobian
 * determinant of the transform.
 *
 * <p>The transform is as specified for
 * <code>ub_constrain(T, double)</code>.  The log absolute Jacobian
 * determinant is
 *
 * <p>\f$ \log | \frac{d}{dx} -\mbox{exp}(x) + U |
 *     = \log | -\mbox{exp}(x) + 0 | = x\f$.
 *
 * @tparam T type of scalar
 * @tparam U type of upper bound
 * @tparam S type of log probability
 * @param[in] x free scalar
 * @param[in] ub upper bound
 * @param[in,out] lp log density
 * @return scalar constrained to have upper bound
 */
template <typename T, typename L, require_all_stan_scalar_t<T, L>* = nullptr, require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  if (is_positive_infinity(ub)) {
    return identity_constrain(x, ub);
  } else {
    //check_greater("ub_constrain", "ub", value_of(x), value_of(ub));
    lp += x;
    return subtract(ub, exp(x));
  }
}


template <typename T, typename L, require_eigen_t<T>* = nullptr,
  require_stan_scalar_t<L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub) {
  return x.unaryExpr([ub](auto&& xx) {
    return ub_constrain(xx, ub);
  });
}

template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub) {
  return x.binaryExpr(ub, [](auto&& xx, auto&& ubb) {
    return ub_constrain(xx, ubb);
  });
}

template <typename T, typename L, require_eigen_t<T>* = nullptr,
  require_stan_scalar_t<L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  return eval(to_ref(x).unaryExpr([ub, &lp](auto&& xx) {
    return ub_constrain(xx, ub, lp);
  }));
}

template <typename T, typename L, require_all_eigen_t<T, L>* = nullptr,
  require_all_not_st_var<T, L>* = nullptr>
inline auto ub_constrain(const T& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  return eval(to_ref(x).binaryExpr(ub, [&lp](auto&& xx, auto&& ubb) {
    return ub_constrain(xx, ubb, lp);
  }));
}

// VEC

template <typename T, typename L,
  require_stan_scalar_t<L>* = nullptr>
inline auto ub_constrain(const std::vector<T>& x, const L& ub) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub);
  }
  return ret;
}

template <typename T, typename L, require_container_t<L>* = nullptr>
inline auto ub_constrain(const std::vector<T>& x, const L& ub) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub[i]);
  }
  return ret;
}

template <typename T, typename L,
  require_stan_scalar_t<L>* = nullptr>
inline auto ub_constrain(const std::vector<T>& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub, lp);
  }
  return ret;
}

template <typename T, typename L, require_container_t<L>* = nullptr>
inline auto ub_constrain(const std::vector<T>& x, const L& ub, std::decay_t<return_type_t<T, L>>& lp) {
  std::vector<promote_scalar_t<return_type_t<T, L>, T>> ret;
  ret.reserve(x.size());
  for (size_t i = 0; i < x.size(); ++i) {
    ret[i] = ub_constrain(x[i], ub[i], lp);
  }
  return ret;
}


}  // namespace math
}  // namespace stan

#endif
