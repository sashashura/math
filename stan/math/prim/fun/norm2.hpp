#ifndef STAN_MATH_PRIM_FUN_NORM2_HPP
#define STAN_MATH_PRIM_FUN_NORM2_HPP

#include <stan/math/prim/meta.hpp>
#include <stan/math/prim/err.hpp>
#include <stan/math/prim/fun/Eigen.hpp>

namespace stan {
namespace math {

/**
 * Returns L2 norm of a vector. For vectors that equals the square-root of the
 * sum of squares of the elements.
 *
 * @tparam T type of the vector (must be derived from \c Eigen::MatrixBase)
 * @param x Vector.
 * @return L2 norm of x.
 */
template <typename T, require_eigen_vt<std::is_arithmetic, T>* = nullptr>
inline double norm2(const T& v) {
  ref_type_t<T> v_ref = v;
  return v_ref.template lpNorm<2>();
}

/**
 * Returns L2 norm of a vector. For vectors that equals the square-root of the
 * sum of squares of the elements.
 *
 * @tparam Container type of the vector (must be derived from \c std::vector)
 * @param v Vector.
 * @return L2 norm of v.
 */
template <typename Container, require_std_vector_t<Container>* = nullptr>
inline auto norm2(const Container& x) {
  return apply_vector_unary<Container>::reduce(
      x, [](const auto& v) { return norm2(v); });
}

}  // namespace math
}  // namespace stan

#endif
