#include <stan/math/rev/core.hpp>
#include <test/unit/math/rev/util.hpp>
#include <gtest/gtest.h>
#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

/**
 * Test Eigen::VectorXd return types
 */
struct SinFunctor {
  int N_;
  double* x_mem_;
  template <std::size_t size>
  Eigen::VectorXd operator()(const std::array<bool, size>& needs_adj,
                             const Eigen::VectorXd& x) {
    N_ = x.size();
    Eigen::VectorXd out(N_);
    x_mem_
        = stan::math::ChainableStack::instance_->memalloc_.alloc_array<double>(
            N_);

    for (int n = 0; n < N_; ++n) {
      x_mem_[n] = x(n);
      out(n) = sin(x(n));
    }

    return out;
  }

  template <std::size_t size>
  auto multiply_adjoint_jacobian(const std::array<bool, size>& needs_adj,
                                 const Eigen::VectorXd& adj) {
    Eigen::VectorXd out(N_);

    for (int n = 0; n < N_; ++n) {
      out(n) = cos(x_mem_[n]) * adj(n);
    }

    return std::make_tuple(out);
  }
};

TEST(AgradRev, test_vector_sin_multiple_jac_static_matrix) {
  using stan::math::var_value;
  Eigen::Matrix<double, Eigen::Dynamic, 1> x1(1), x2(2);
  x1 << 1.0;
  x2 << 2.0, 1.0;
  stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>> xx1(x1);
  stan::math::var_value<Eigen::Matrix<double, Eigen::Dynamic, 1>> xx2(x2);
  auto y1 = stan::math::adj_jac_apply<SinFunctor>(xx1);
  auto y2 = stan::math::adj_jac_apply<SinFunctor>(xx2);

  y1.grad();
  EXPECT_FLOAT_EQ(xx1.adj()(0), 0.5403023058681398);
  EXPECT_FLOAT_EQ(xx2.adj()(0), 0.0);
  EXPECT_FLOAT_EQ(xx2.adj()(1), 0.0);

  stan::math::set_zero_all_adjoints();
  stan::math::set_zero_all_adjoints();

  y2.grad();
  EXPECT_FLOAT_EQ(xx1.adj()(0), 0.0);
  EXPECT_FLOAT_EQ(xx2.adj()(0), -0.4161468365471424);
  EXPECT_FLOAT_EQ(xx2.adj()(1), 0.5403023058681398);

  stan::math::set_zero_all_adjoints();
  Eigen::Matrix<double, 1, -1> sum_vec(2);
  sum_vec << 1.73, 1.57;
  auto sum_y2 = sum_vec * y2;
  sum_y2.grad();
  EXPECT_FLOAT_EQ(xx1.adj()(0), 0.0);
  EXPECT_FLOAT_EQ(xx2.adj()(0), 1.73 * -0.4161468365471424);
  EXPECT_FLOAT_EQ(xx2.adj()(1), 1.57 * 0.5403023058681398);
}