#ifndef STAN_MATH_REV_FUN_LOG2_HPP
#define STAN_MATH_REV_FUN_LOG2_HPP

#include <stan/math/rev/meta.hpp>
#include <stan/math/rev/core.hpp>
#include <stan/math/prim/fun/log2.hpp>
#include <stan/math/prim/fun/constants.hpp>

namespace stan {
namespace math {

/**
 * Returns the base 2 logarithm of the specified variable (C99).
 *
 * See log2() for the double-based version.
 *
 * The derivative is
 *
 * \f$\frac{d}{dx} \log_2 x = \frac{1}{x \log 2}\f$.
 *
   \f[
   \mbox{log2}(x) =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0 \\
     \log_2(x) & \mbox{if } x\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{log2}(x)}{\partial x} =
   \begin{cases}
     \textrm{NaN} & \mbox{if } x < 0 \\
     \frac{1}{x\ln2} & \mbox{if } x\geq 0 \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]
 *
 * @tparam T Arithmetic or a type inheriting from `EigenBase`.
 * @param a The variable.
 * @return Base 2 logarithm of the variable.
 */
inline auto log2(const var& a) {
  return make_callback_var(log2(a.val()), [a](auto& vi) mutable {
    a.adj() += vi.adj() / (LOG_TWO * a.val());
  });
}

template <typename T, require_rev_matrix_t<T>* = nullptr>
inline auto log2(const T& a) {
  auto a_arena = to_arena(a);
  return make_callback_rev_matrix<T>(log2(a_arena.val()), [a_arena](auto&& vi) mutable {
    a_arena.adj().array() += vi.adj().array() / (LOG_TWO * a_arena.val().array());
  });
}

}  // namespace math
}  // namespace stan
#endif
