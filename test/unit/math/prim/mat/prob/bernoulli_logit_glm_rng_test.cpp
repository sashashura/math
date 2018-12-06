#include <stan/math/prim/mat.hpp>
#include <gtest/gtest.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/distributions.hpp>
#include <test/unit/math/prim/scal/prob/util.hpp>
#include <limits>
#include <vector>

TEST(ProbDistributionsBernoulliLogitGlm, NotVectorized) {
  boost::random::mt19937 rng;
  Eigen::Matrix<double, 1, Eigen::Dynamic> x(2);
  x << 3.5, -1.5;
  double alpha = 2.0;
  std::vector<double> beta{2.0, 4.5};
  EXPECT_NO_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, beta, rng));
}

TEST(ProbDistributionsBernoulliLogitGlm, Vectorized) {
  boost::random::mt19937 rng;
  Eigen::MatrixXd x(2, 2);
  x << 3.5, -1.5,
       4.0, 1.2;
  std::vector<double> alpha{2.0, 1.0};
  std::vector<double> beta{2.0, 4.5};
  EXPECT_NO_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, beta, rng));
}

//  We check that the right errors are thrown.
TEST(ProbDistributionsPoissonLogGLM,
     glm_matches_bernoulli_logit_error_checking) {
  boost::random::mt19937 rng;

  int N = 3;
  int M = 2;
  int W = 4;

  Eigen::MatrixXd x = Eigen::MatrixXd::Random(N, M);
  Eigen::MatrixXd xw1 = Eigen::MatrixXd::Random(W, M);
  Eigen::MatrixXd xw2 = Eigen::MatrixXd::Random(N, W);
  Eigen::MatrixXd xw3 = Eigen::MatrixXd::Random(N, M) * NAN;
  Eigen::VectorXd alpha = Eigen::VectorXd::Random(N, 1);
  Eigen::VectorXd alphaw1 = Eigen::VectorXd::Random(W, 1);
  Eigen::VectorXd alphaw2 = Eigen::VectorXd::Random(N, 1) * NAN;
  Eigen::VectorXd beta = Eigen::VectorXd::Random(M, 1);
  Eigen::VectorXd betaw1 = Eigen::VectorXd::Random(W, 1);
  Eigen::VectorXd betaw2 = Eigen::VectorXd::Random(M, 1) * NAN;

  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(xw1, alpha, beta, rng),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(xw2, alpha, beta, rng),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(xw3, alpha, beta, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alphaw1, beta, rng),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alphaw2, beta, rng),
               std::domain_error);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, betaw1, rng),
               std::invalid_argument);
  EXPECT_THROW(stan::math::bernoulli_logit_glm_rng(x, alpha, betaw2, rng),
               std::domain_error);
}

TEST(ProbDistributionsBernoulliLogitGlm, marginalChiSquareGoodnessFitTest) {
  boost::random::mt19937 rng;
  Eigen::MatrixXd x(2, 2);
  x << 3.5, -1.5,
       2.0, -1.2;
  std::vector<double> alpha{2.0, 1.0};
  std::vector<double> beta{2.0, 4.5};

  //  sage: x = matrix([[3.5, -1.5], [2.0, -1.2]])
  //  sage: alpha = matrix([[2.0], [1.0]])
  //  sage: beta = matrix([[2.0], [4.5]])
  //  sage: z = alpha + x * beta
  //  sage: z
  //
  //  [  2.25000000000000]
  //  [-0.399999999999999]
  //  sage: p1 = 1 / (1 + e**(-z[0][0]))
  //  sage: p2 = 1 / (1 + e**(-z[1][0]))
  //  sage: p1
  //  0.904650535100891
  //  sage: p2
  //  0.401312339887548

  double p1 = 0.904650535100891;
  double p2 = 0.401312339887548;

  int N = 10000;

  // First bin is failures, second is successes. Take N samples, take
  // their first component. Should be (1 - p1) * 100%
  // failures. Likewise for the failures bin, and the second
  // component.
  std::vector<double> percentile_cutpoints{0.1, 1.1};
  std::vector<double> ps1{(1 - p1), p1};
  std::vector<double> ps2{(1 - p2), p2};

  std::vector<double> samples1;
  std::vector<double> samples2;
  for (int i = 0; i < N; ++i) {
    std::vector<int> sample =
      stan::math::bernoulli_logit_glm_rng(x, alpha, beta, rng);
    samples1.push_back(sample[0]);
    samples2.push_back(sample[1]);
  }

  assert_matches_cutpoints(samples1, percentile_cutpoints, ps1, 1e-6);
  assert_matches_cutpoints(samples2, percentile_cutpoints, ps2, 1e-6);
}
