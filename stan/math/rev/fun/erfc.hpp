#ifndef STAN_MATH_REV_FUN_ERFC_HPP
#define STAN_MATH_REV_FUN_ERFC_HPP

#include <stan/math/rev/core.hpp>
#include <stan/math/rev/meta.hpp>
#include <stan/math/prim/fun/constants.hpp>
#include <stan/math/prim/fun/erfc.hpp>
#include <cmath>

namespace stan {
namespace math {

/**
 * The complementary error function for variables (C99).
 *
 * The derivative is
 *
 * \f$\frac{d}{dx} \mbox{erfc}(x) = - \frac{2}{\sqrt{\pi}} \exp(-x^2)\f$.
 *
 *
   \f[
   \mbox{erfc}(x) =
   \begin{cases}
     \operatorname{erfc}(x) & \mbox{if } -\infty\leq x \leq \infty \\[6pt]
     \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \frac{\partial\, \mbox{erfc}(x)}{\partial x} =
   \begin{cases}
     \frac{\partial\, \operatorname{erfc}(x)}{\partial x} & \mbox{if }
 -\infty\leq x\leq \infty \\[6pt] \textrm{NaN} & \mbox{if } x = \textrm{NaN}
   \end{cases}
   \f]

   \f[
   \operatorname{erfc}(x)=\frac{2}{\sqrt{\pi}}\int_x^\infty e^{-t^2}dt
   \f]

   \f[
   \frac{\partial \, \operatorname{erfc}(x)}{\partial x} = -\frac{2}{\sqrt{\pi}}
 e^{-x^2} \f]
 *
 * @param a The variable.
 * @return Complementary error function applied to the variable.
 */
inline var erfc(const var& a) {
  return make_callback_var(erfc(a.val()), [a](auto& vi) mutable {
    a.adj() -= vi.adj() * TWO_OVER_SQRT_PI * std::exp(-a.val() * a.val());
  });
}

}  // namespace math
}  // namespace stan
#endif
