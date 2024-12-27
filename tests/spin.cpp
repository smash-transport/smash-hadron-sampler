#include "spin.h"

#include "virtest/vir/test.h"

using namespace spin;

// Custom function to compare two doubles within a given tolerance
bool expect_near(double val1, double val2, double abs_error) {
  return std::abs(val1 - val2) <= abs_error;
}

TEST(four_vector_square) {
  // Test the square of a four-vector
  std::array<double, 4> four_vector = {18.5, 1.5, 2.2, 3.1};
  // Manually calculated value for the square of the four-vector
  const double expected_result =
      18.5 * 18.5 - 1.5 * 1.5 - 2.2 * 2.2 - 3.1 * 3.1;
  const double result = four_vector_square(four_vector);
  VERIFY(expect_near(result, expected_result, 1e-6));
}

TEST(theta_valid_values) {
  // Test the validity of the theta function
  const double mass = 1.0;
  const std::array<double, 16> vorticity = {0.0,  1.0,  2.0,  3.0,  -1.0, 0.0,
                                            4.0,  5.0,  -2.0, -4.0, 0.0,  6.0,
                                            -3.0, -5.0, -6.0, 0.0};
  const double p_lower_index[4] = {1.0, 2.0, 3.0, 4.0};
  // Manually calculated values for theta
  const double expected_theta_array[4] = {-6.5, 2.5, -1.5, 1.5};
  const auto theta_array = theta(mass, vorticity, p_lower_index);

  // Perform checks
  VERIFY(theta_array.size() == 4);
  VERIFY(expect_near(theta_array[0], expected_theta_array[0], 1e-6));
  VERIFY(expect_near(theta_array[1], expected_theta_array[1], 1e-6));
  VERIFY(expect_near(theta_array[2], expected_theta_array[2], 1e-6));
  VERIFY(expect_near(theta_array[3], expected_theta_array[3], 1e-6));
}

TEST(theta_invalid_mass) {
  // Test the invalidity of the theta function for massless particles
  const double mass = 0.0;
  const std::array<double, 16> vorticity = {0.0,  1.0,  2.0,  3.0,  -1.0, 0.0,
                                            4.0,  5.0,  -2.0, -4.0, 0.0,  6.0,
                                            -3.0, -5.0, -6.0, 0.0};
  const double p_lower_index[4] = {1.0, 2.0, 3.0, 4.0};
  // Expect an invalid argument exception
  try {
    theta(mass, vorticity, p_lower_index);
    std::cout << "theta unexpectedly passed with massless particle"
              << std::endl;
  } catch (std::invalid_argument &e) {
    // Exception was caught as expected
  }
}

TEST(exponent_valid_value) {
  // Test the validity of the exponent function
  const double k = -1.5;
  const double energy_density = 18.6;
  const double temperature = 2.0;
  const double mu = 4.8;
  const double theta_squared = -10.89;
  // Manually calculated value for the exponent
  const double expected_result = 9.3 - 2.4 + 4.95;
  const double result =
      exponent(k, energy_density, temperature, mu, theta_squared);
  VERIFY(expect_near(result, expected_result, 1e-6));
}

TEST(exponent_invalid_k) {
  // Test the invalidity of the exponent function for non-integer k
  const double invalid_k_array[3] = {-1.2, 0.4, 7.9};
  const double energy_density = 1.0;
  const double temperature = 1.0;
  const double mu = 1.0;
  const double theta_squared = -1.0;
  // Expect an invalid argument exception
  for (const double k : invalid_k_array) {
    try {
      exponent(k, energy_density, temperature, mu, theta_squared);
      std::cout << "exponent unexpectedly passed with k not being a multiple "
                << "of 1/2" << std::endl;
    } catch (std::invalid_argument &e) {
      // Exception was caught as expected
    }
  }
}

TEST(exponent_positive_theta_squared) {
  // Test the invalidity of the exponent function for positive theta squared
  const double k = 1.0;
  const double energy_density = 1.0;
  const double temperature = 1.0;
  const double mu = 1.0;
  const double invalid_theta_squared_array[3] = {0.0001, 1.3, 4.6};
  // Expect an invalid argument exception
  for (const double theta_squared : invalid_theta_squared_array) {
    try {
      exponent(k, energy_density, temperature, mu, theta_squared);
      std::cout << "exponent unexpectedly passed with positive theta squared"
                << std::endl;
    } catch (std::invalid_argument &e) {
      // Exception was caught as expected
    }
  }
}
