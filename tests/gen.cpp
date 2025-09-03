#include "gen.h"

#include <ctime>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "virtest/vir/test.h"

using namespace gen;

TEST(index44) {
  // Call with all valid values and expect passing
  int indices[4][4] = {{0, 1, 3, 6}, {1, 2, 4, 7}, {3, 4, 5, 8}, {6, 7, 8, 9}};
  for (int row = 0; row < 4; row++) {
    for (int column = 0; column < 4; column++) {
      VERIFY(index44(row, column) == indices[row][column]);
    }
  }
}

TEST(index44_index_out_of_range) {
  int valid_index[4] = {0, 1, 2, 3};
  int invalid_index[4] = {-2, -1, 4, 5};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      // first index out of range
      try {
        index44(invalid_index[i], valid_index[j]);
        std::cout << "index44 unexpectedly passed with first argument "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
      // second index out of range
      try {
        index44(valid_index[i], invalid_index[j]);
        std::cout << "index44 unexpectedly passed with second argument "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
      // both indices out of range
      try {
        index44(invalid_index[i], invalid_index[j]);
        std::cout << "index44 unexpectedly passed with both arguments "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
    }
  }
}

TEST(fillBoostMatrix) {
  const double v[3] = {0.3, 0.2, -0.1};
  const double gamma =
      1.0 / sqrt(1.0 - (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
  const double u[4] = {gamma, gamma * v[0], gamma * v[1], gamma * v[2]};

  double boostMatrix[4][4];
  fillBoostMatrix(-u[1] / gamma, -u[2] / gamma, -u[3] / gamma, boostMatrix);

  // Perform boost on u to get back to (1,0,0,0)
  double u_prime[4];
  for (int mu = 0; mu < 4; ++mu) {
    u_prime[mu] = 0.0;
    for (int nu = 0; nu < 4; ++nu) {
      u_prime[mu] += boostMatrix[mu][nu] * u[nu];
    }
  }

  // Check that u' is approximately (1,0,0,0)
  const double tolerance = 1e-9;
  VERIFY(std::abs(u_prime[0] - 1.0) < tolerance);
  VERIFY(std::abs(u_prime[1]) < tolerance);
  VERIFY(std::abs(u_prime[2]) < tolerance);
  VERIFY(std::abs(u_prime[3]) < tolerance);
}

TEST(create_covariant_boost_matrix) {
  // create covariant 4-velocity
  const double v[3] = {0.3, 0.2, -0.1};
  const double gamma =
      1.0 / sqrt(1.0 - (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
  const double u[4] = {gamma, gamma * v[0], gamma * v[1], gamma * v[2]};
  const double u_cov[4] = {gamma, -gamma * v[0], -gamma * v[1], -gamma * v[2]};

  double boost_contravariant[4][4];
  fillBoostMatrix(-u[1] / gamma, -u[2] / gamma, -u[3] / gamma,
                  boost_contravariant);

  FourMatrix boost_covariant =
      create_covariant_boost_matrix(boost_contravariant);

  // Perform boost on u to get back to (1,0,0,0)
  double u_prime[4];
  for (int mu = 0; mu < 4; ++mu) {
    u_prime[mu] = 0.0;
    for (int nu = 0; nu < 4; ++nu) {
      u_prime[mu] += boost_covariant[mu][nu] * u_cov[nu];
    }
  }

  // Check that u' is approximately (1,0,0,0)
  const double tolerance = 1e-9;
  VERIFY(std::abs(u_prime[0] - 1.0) < tolerance);
  VERIFY(std::abs(u_prime[1]) < tolerance);
  VERIFY(std::abs(u_prime[2]) < tolerance);
  VERIFY(std::abs(u_prime[3]) < tolerance);
}
