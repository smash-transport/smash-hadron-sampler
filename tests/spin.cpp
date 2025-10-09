#include "spin.h"

#include <array>
#include <cmath>
#include <iomanip>
#include <mutex>
#include <random>

#include "gen.h"
#include "smash/pdgcode.h"
#include "virtest/vir/test.h"
#include "vorticity.h"

using namespace spin;

// Ensures the global ParticleType list is initialized only once using
// std::call_once. Avoids the "Type list was already built!" exception when
// multiple tests need it.
void ensure_particletype_initialized() {
  static std::once_flag flag;
  std::call_once(flag, [] {
    smash::ParticleType::create_type_list(
        "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
        "π⁰ 0.1380 0      - 111\n"
        "ρ  0.776  0.149  - 113  213\n"
        "N⁺ 0.938  0      + 2212\n"
        "Δ  1.232  0.117  + 2224 2214 2114 1114\n"
        "Λ  1.116 2.5e-15 + 3122  # stable by default");
  });
}

// Precomputed lookup table for Levi-Civita symbol
int levi_civita(int i, int j, int k, int l) {
  static const std::array<std::array<int, 4>, 24> permutations = {
      {{{0, 1, 2, 3}}, {{0, 1, 3, 2}}, {{0, 2, 1, 3}}, {{0, 2, 3, 1}},
       {{0, 3, 1, 2}}, {{0, 3, 2, 1}}, {{1, 0, 2, 3}}, {{1, 0, 3, 2}},
       {{1, 2, 0, 3}}, {{1, 2, 3, 0}}, {{1, 3, 0, 2}}, {{1, 3, 2, 0}},
       {{2, 0, 1, 3}}, {{2, 0, 3, 1}}, {{2, 1, 0, 3}}, {{2, 1, 3, 0}},
       {{2, 3, 0, 1}}, {{2, 3, 1, 0}}, {{3, 0, 1, 2}}, {{3, 0, 2, 1}},
       {{3, 1, 0, 2}}, {{3, 1, 2, 0}}, {{3, 2, 0, 1}}, {{3, 2, 1, 0}}}};

  static const std::array<int, 24> signs = {+1, -1, -1, +1, +1, -1, -1, +1,
                                            +1, -1, -1, +1, +1, -1, -1, +1,
                                            +1, -1, -1, +1, +1, -1, -1, +1};

  // Check for repeated indices
  if (i == j || i == k || i == l || j == k || j == l || k == l) {
    return 0;
  }

  // Match the input with precomputed permutations
  for (size_t p = 0; p < permutations.size(); ++p) {
    if (permutations[p] == std::array<int, 4>{i, j, k, l}) {
      return signs[p];
    }
  }
  return 0;  // Should never reach here
}

// Custom function to compare two doubles within a given tolerance
bool expect_near(double val1, double val2, double abs_error) {
  return std::abs(val1 - val2) <= abs_error;
}

TEST(four_vector_square) {
  // Create a random number generator
  std::random_device rd;   // Seed source
  std::mt19937 gen(rd());  // Mersenne Twister generator

  // Define the range [-10.0, 10.0]
  std::uniform_real_distribution<double> dist(-10.0, 10.0);

  // Generate random 4-vector components and calculate the square
  for (int i = 0; i < 20; ++i) {
    double v0 = dist(gen);
    double v1 = dist(gen);
    double v2 = dist(gen);
    double v3 = dist(gen);
    std::array<double, 4> four_vector = {v0, v1, v2, v3};
    double expected_result = v0 * v0 - v1 * v1 - v2 * v2 - v3 * v3;
    double result = four_vector_square(four_vector);
    VERIFY(expect_near(result, expected_result, 1e-6));
  }
}

TEST(theta_valid_values) {
  // Test the validity of the theta function with an antisymmetric
  // vorticity tensor
  const double mass = 0.5;
  const std::array<double, 16> vorticity = {0.0,  1.0,  2.0,  3.0,  -1.0, 0.0,
                                            4.0,  5.0,  -2.0, -4.0, 0.0,  6.0,
                                            -3.0, -5.0, -6.0, 0.0};
  const std::array<double, 4> p = {1.0, 2.0, 3.0, 4.0};
  // Manually calculated values for theta
  const double expected_theta_array[4] = {26.0 * hbarC, 14.0 * hbarC,
                                          -14.0 * hbarC, 10.0 * hbarC};
  const auto theta_array = theta(mass, vorticity, p);
  // Perform checks
  VERIFY(theta_array.size() == 4);
  VERIFY(expect_near(theta_array[0], expected_theta_array[0], 1e-9));
  VERIFY(expect_near(theta_array[1], expected_theta_array[1], 1e-9));
  VERIFY(expect_near(theta_array[2], expected_theta_array[2], 1e-9));
  VERIFY(expect_near(theta_array[3], expected_theta_array[3], 1e-9));

  // Test with second set of values
  const double mass2 = 0.9;
  const std::array<double, 16> vorticity2 = {0.0,  1.9,  -5.1, 3.7, -1.9, 0.0,
                                             -6.3, -4.9, 5.1,  6.3, 0.0,  -5.6,
                                             -3.7, 4.9,  5.6,  0.0};
  const std::array<double, 4> p2 = {2.8, -1.3, 3.4, -4.1};
  const double expected_theta_array2[4] = {
      55.3 * hbarC, -26.67777777778 * hbarC, 11.9333333333333 * hbarC,
      -19.411111111111 * hbarC};
  const auto theta_array2 = theta(mass2, vorticity2, p2);
  // Perform checks
  VERIFY(theta_array2.size() == 4);
  VERIFY(expect_near(theta_array2[0], expected_theta_array2[0], 1e-9));
  VERIFY(expect_near(theta_array2[1], expected_theta_array2[1], 1e-9));
  VERIFY(expect_near(theta_array2[2], expected_theta_array2[2], 1e-9));
  VERIFY(expect_near(theta_array2[3], expected_theta_array2[3], 1e-9));
}

TEST(theta_tensor_structure) {
  // This test performs all tensor contractions explicitely to verify that the
  // tensor structure of the theta function is implemented correctly.
  const double mass = 0.55;
  const std::array<double, 4> p = {5.5, 2.2, 3.3, 1.1};
  const std::array<double, 16> vorticity = {0.0,  1.7,  2.1,  3.4,  -1.7, 0.0,
                                            4.4,  5.9,  -2.1, -4.4, 0.0,  6.3,
                                            -3.4, -5.9, -6.3, 0.0};
  const std::array<std::array<int, 4>, 4> g = {
      {{1, 0, 0, 0}, {0, -1, 0, 0}, {0, 0, -1, 0}, {0, 0, 0, -1}}};

  auto at = [](int i, int j) { return (i * 4 + j); };

  // Calculate the theta tensor components from the function
  std::array<double, 4> theta_calculated = theta(mass, vorticity, p);

  // Calculate the theta tensor components manually
  double theta_expected[4] = {0.0, 0.0, 0.0, 0.0};
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = 0; nu < 4; nu++) {
      for (int rho = 0; rho < 4; rho++) {
        for (int sigma = 0; sigma < 4; sigma++) {
          for (int gamma = 0; gamma < 4; gamma++) {
            theta_expected[mu] += hbarC * (-1.0 / (2.0 * mass)) *
                                  levi_civita(mu, nu, rho, sigma) *
                                  vorticity[at(nu, rho)] * g[sigma][gamma] *
                                  p[gamma];
          }
        }
      }
    }
  }
  // Compare the calculated and expected values
  VERIFY(expect_near(theta_calculated[0], theta_expected[0], 1e-9));
  VERIFY(expect_near(theta_calculated[1], theta_expected[1], 1e-9));
  VERIFY(expect_near(theta_calculated[2], theta_expected[2], 1e-9));
  VERIFY(expect_near(theta_calculated[3], theta_expected[3], 1e-9));
}

TEST(theta_invalid_mass) {
  // Test the invalidity of the theta function for massless particles
  const double mass = 0.0;
  const std::array<double, 16> vorticity = {0.0,  1.0,  2.0,  3.0,  -1.0, 0.0,
                                            4.0,  5.0,  -2.0, -4.0, 0.0,  6.0,
                                            -3.0, -5.0, -6.0, 0.0};
  const std::array<double, 4> p = {1.0, 2.0, 3.0, 4.0};
  // Expect an invalid argument exception
  try {
    theta(mass, vorticity, p);
    std::cout << "theta unexpectedly passed with massless particle"
              << std::endl;
  } catch (std::invalid_argument &e) {
    // Exception was caught as expected
  }
}

TEST(exponent_valid_value) {
  // Test the validity of the exponent function
  const double k = -1.5;
  const double energy = 18.6;
  const double temperature = 2.0;
  const double mu = 4.8;
  const double theta_squared = -10.89;
  // Manually calculated value for the exponent
  const double expected_result = 9.3 - 2.4 + 4.95;
  const double result = exponent(k, energy, temperature, mu, theta_squared);
  VERIFY(expect_near(result, expected_result, 1e-6));
}

TEST(exponent_inconsistent_input_values) {
  // We start with a valid value set and then change the values one by one
  // to test that the exponent function does not coincide with the expected
  // result anymore.
  double input_values[5] = {2.5, 7.54, 3.71, 2.33, -22.5625};
  // Manually calculated value for the exponent
  const double expected_result = 2.03234501 - 0.62803235 - 11.875;
  double result = exponent(input_values[0], input_values[1], input_values[2],
                           input_values[3], input_values[4]);

  VERIFY(expect_near(result, expected_result, 1e-6));

  for (int index = 0; index < 5; index++) {
    // Change one value to an invalid value and store the original value
    double tmp = input_values[index];
    input_values[index] = -0.5;
    // Ensure that the result is not equal to the expected result anymore
    result = exponent(input_values[0], input_values[1], input_values[2],
                      input_values[3], input_values[4]);
    VERIFY(!expect_near(result, expected_result, 1e-6));
    // Reset the value to the original value
    input_values[index] = tmp;
  }
}

TEST(exponent_invalid_k) {
  // Test the invalidity of the exponent function for non-integer k
  const double invalid_k_array[3] = {-1.2, 0.4, 7.9};
  const double energy = 1.0;
  const double temperature = 1.0;
  const double mu = 1.0;
  const double theta_squared = -1.0;
  // Expect an invalid argument exception
  for (const double k : invalid_k_array) {
    try {
      exponent(k, energy, temperature, mu, theta_squared);
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
  const double energy = 1.0;
  const double temperature = 1.0;
  const double mu = 1.0;
  const double invalid_theta_squared_array[3] = {0.0001, 1.3, 4.6};
  // Expect an invalid argument exception
  for (const double theta_squared : invalid_theta_squared_array) {
    try {
      exponent(k, energy, temperature, mu, theta_squared);
      std::cout << "exponent unexpectedly passed with positive theta squared"
                << std::endl;
    } catch (std::invalid_argument &e) {
      // Exception was caught as expected
    }
  }
}

TEST(spin_vector_valid_values) {
  // Initialize the particle type list with a pion (spin 0),
  // proton (spin 1/2), and delta plus (spin 3/2)
  ensure_particletype_initialized();

  // Define the PDG code for the particles
  smash::PdgCode pdg_pion = 0x111;
  smash::PdgCode pdg_proton = 0x2212;
  smash::PdgCode pdg_rho = 0x213;
  smash::PdgCode pdg_delta_plus = 0x2214;

  // Create the particle data for the pion, proton, and delta plus
  smash::ParticleData pion =
      smash::ParticleData(smash::ParticleType::find(pdg_pion));
  smash::ParticleData proton =
      smash::ParticleData(smash::ParticleType::find(pdg_proton));
  smash::ParticleData rho =
      smash::ParticleData(smash::ParticleType::find(pdg_rho));
  smash::ParticleData delta_plus =
      smash::ParticleData(smash::ParticleType::find(pdg_delta_plus));

  // Set particle momenta
  const double pion_mom[4] = {
      std::sqrt(pion.pole_mass() * pion.pole_mass() + 0.11 * 0.11 +
                0.22 * 0.22 + 0.33 * 0.33),
      0.11, 0.22, 0.33};
  const std::array<double, 4> proton_mom = {
      std::sqrt(proton.pole_mass() * proton.pole_mass() + 0.11 * 0.11 +
                0.22 * 0.22 + 0.33 * 0.33),
      0.11, 0.22, 0.33};
  const std::array<double, 4> rho_mom = {
      std::sqrt(rho.pole_mass() * rho.pole_mass() + 0.11 * 0.11 + 0.22 * 0.22 +
                0.33 * 0.33),
      0.11, 0.22, 0.33};
  const std::array<double, 4> delta_plus_mom = {
      std::sqrt(delta_plus.pole_mass() * delta_plus.pole_mass() + 0.11 * 0.11 +
                0.22 * 0.22 + 0.33 * 0.33),
      0.11, 0.22, 0.33};
  pion.set_4momentum(pion.pole_mass(), pion_mom[1], pion_mom[2], pion_mom[3]);
  proton.set_4momentum(proton.pole_mass(), proton_mom[1], proton_mom[2],
                       proton_mom[3]);
  rho.set_4momentum(rho.pole_mass(), rho_mom[1], rho_mom[2], rho_mom[3]);
  delta_plus.set_4momentum(delta_plus.pole_mass(), delta_plus_mom[1],
                           delta_plus_mom[2], delta_plus_mom[3]);

  gen::element surf_element;
  // Set all values of the surface element to 1.0 which
  // do not affect the spin sampling
  surf_element.four_position[0] = 1.0;
  surf_element.four_position[1] = 1.0;
  surf_element.four_position[2] = 1.0;
  surf_element.four_position[3] = 1.0;
  for (int i = 0; i < 4; i++) {
    surf_element.u[i] = 1.0;
    surf_element.dsigma[i] = 1.0;
  }
  surf_element.muq = 1.7;
  surf_element.mus = 2.3;
  for (int i = 0; i < 10; i++) {
    surf_element.pi[i] = 1.0;
  }
  surf_element.Pi = 1.0;
  // Set all relevant values for the spin sampling
  surf_element.T = 4.3;
  surf_element.mub = 0.7;
  surf_element.e = 1.7;

  // This test does not test the boost behavior, so we can set the boost
  // velocity to zero. The particle momentum is already given in the
  // fluid rest frame, so no boost is needed
  const smash::ThreeVector boost_velocity(0.0, 0.0, 0.0);

  // Create a Vorticity instance and set the vorticity components in the surface
  auto vorticity = std::make_unique<Vorticity>();
  vorticity->set_vorticity(
      std::array<double, 16>{0.0, 1.9, -5.1, 3.7, -1.9, 0.0, -6.3, -4.9, 5.1,
                             6.3, 0.0, -5.6, -3.7, 4.9, 5.6, 0.0});
  surf_element.vorticity = gen::OptionalVorticity(std::move(vorticity));

  // TEST FOR SPIN O (PION)
  // Newly initialized pion should have a spin vector of nan
  VERIFY(
      std::isnan(pion.spin_vector()[0]) && std::isnan(pion.spin_vector()[1]) &&
      std::isnan(pion.spin_vector()[2]) && std::isnan(pion.spin_vector()[3]));

  // Calculate and set the spin vector for the pion
  calculate_and_set_spin_vector(0, surf_element, &pion, boost_velocity);

  // Perform checks
  VERIFY(expect_near(pion.spin_vector()[0], 0.0, 1e-9));
  VERIFY(expect_near(pion.spin_vector()[1], 0.0, 1e-9));
  VERIFY(expect_near(pion.spin_vector()[2], 0.0, 1e-9));
  VERIFY(expect_near(pion.spin_vector()[3], 0.0, 1e-9));

  // TEST FOR SPIN 1/2 (PROTON)
  // Newly initialized proton should have a spin vector of nan
  VERIFY(std::isnan(proton.spin_vector()[0]) &&
         std::isnan(proton.spin_vector()[1]) &&
         std::isnan(proton.spin_vector()[2]) &&
         std::isnan(proton.spin_vector()[3]));

  // Calculate and set the spin vector for the proton
  calculate_and_set_spin_vector(0, surf_element, &proton, boost_velocity);

  // Calculate the expected values for the spin vector
  std::array<double, 4> theta_vector =
      theta(proton.pole_mass(), (**surf_element.vorticity).get_vorticity(),
            proton_mom);
  double theta_squared = spin::four_vector_square(theta_vector);

  double mu_proton = chemical_potential(&proton, surf_element);

  double denominator_1 =
      1 /
      (std::exp(exponent(-0.5, proton_mom[0], 4.3, mu_proton, theta_squared)) +
       1);
  double denominator_2 =
      1 /
      (std::exp(exponent(0.5, proton_mom[0], 4.3, mu_proton, theta_squared)) +
       1);

  double numerator_1 = -0.5 * denominator_1;
  double numerator_2 = 0.5 * denominator_2;

  double total_weight =
      (1 / std::sqrt(-theta_squared)) *
      ((numerator_1 + numerator_2) / (denominator_1 + denominator_2));

  double expected_spin_vector[4] = {
      theta_vector[0] * total_weight, theta_vector[1] * total_weight,
      theta_vector[2] * total_weight, theta_vector[3] * total_weight};

  // Perform checks
  VERIFY(expect_near(proton.spin_vector()[0], expected_spin_vector[0], 1e-9));
  VERIFY(expect_near(proton.spin_vector()[1], expected_spin_vector[1], 1e-9));
  VERIFY(expect_near(proton.spin_vector()[2], expected_spin_vector[2], 1e-9));
  VERIFY(expect_near(proton.spin_vector()[3], expected_spin_vector[3], 1e-9));

  // TEST FOR SPIN 1 (RHO)
  // Newly initialized rho should have a spin vector of nan
  VERIFY(std::isnan(rho.spin_vector()[0]) && std::isnan(rho.spin_vector()[1]) &&
         std::isnan(rho.spin_vector()[2]) && std::isnan(rho.spin_vector()[3]));

  // Calculate and set the spin vector for the rho
  calculate_and_set_spin_vector(0, surf_element, &rho, boost_velocity);

  // Calculate the expected values for the spin vector
  theta_vector = theta(rho.pole_mass(),
                       (**surf_element.vorticity).get_vorticity(), rho_mom);
  theta_squared = spin::four_vector_square(theta_vector);

  double mu_rho = chemical_potential(&rho, surf_element);

  denominator_1 =
      1 /
      (std::exp(exponent(-1.0, rho_mom[0], 4.3, mu_rho, theta_squared)) - 1);
  denominator_2 =
      1 / (std::exp(exponent(0.0, rho_mom[0], 4.3, mu_rho, theta_squared)) - 1);
  double denominator_3 =
      1 / (std::exp(exponent(1.0, rho_mom[0], 4.3, mu_rho, theta_squared)) - 1);

  numerator_1 = -1.0 * denominator_1;
  numerator_2 = 0.0 * denominator_2;
  double numerator_3 = 1.0 * denominator_3;

  total_weight = (1 / std::sqrt(-theta_squared)) *
                 ((numerator_1 + numerator_2 + numerator_3) /
                  (denominator_1 + denominator_2 + denominator_3));

  for (int i = 0; i < 4; ++i) {
    expected_spin_vector[i] = theta_vector[i] * total_weight;
  }

  // Perform checks
  VERIFY(expect_near(rho.spin_vector()[0], expected_spin_vector[0], 1e-9));
  VERIFY(expect_near(rho.spin_vector()[1], expected_spin_vector[1], 1e-9));
  VERIFY(expect_near(rho.spin_vector()[2], expected_spin_vector[2], 1e-9));
  VERIFY(expect_near(rho.spin_vector()[3], expected_spin_vector[3], 1e-9));

  // TEST FOR SPIN 3/2 (DELTA PLUS)
  // Newly initialized delta plus should have a spin vector of nan
  VERIFY(std::isnan(delta_plus.spin_vector()[0]) &&
         std::isnan(delta_plus.spin_vector()[1]) &&
         std::isnan(delta_plus.spin_vector()[2]) &&
         std::isnan(delta_plus.spin_vector()[3]));

  // Calculate and set the spin vector for the delta plus
  calculate_and_set_spin_vector(0, surf_element, &delta_plus, boost_velocity);

  // Calculate the expected values for the spin vector
  theta_vector =
      theta(delta_plus.pole_mass(), (**surf_element.vorticity).get_vorticity(),
            delta_plus_mom);
  theta_squared = spin::four_vector_square(theta_vector);

  double mu_delta_plus = chemical_potential(&delta_plus, surf_element);

  denominator_1 = 1 / (std::exp(exponent(-1.5, delta_plus_mom[0], 4.3,
                                         mu_delta_plus, theta_squared)) +
                       1);
  denominator_2 = 1 / (std::exp(exponent(-0.5, delta_plus_mom[0], 4.3,
                                         mu_delta_plus, theta_squared)) +
                       1);
  denominator_3 = 1 / (std::exp(exponent(0.5, delta_plus_mom[0], 4.3,
                                         mu_delta_plus, theta_squared)) +
                       1);
  double denominator_4 = 1 / (std::exp(exponent(1.5, delta_plus_mom[0], 4.3,
                                                mu_delta_plus, theta_squared)) +
                              1);

  numerator_1 = -1.5 * denominator_1;
  numerator_2 = -0.5 * denominator_2;
  numerator_3 = 0.5 * denominator_3;
  double numerator_4 = 1.5 * denominator_4;

  total_weight =
      (1 / std::sqrt(-theta_squared)) *
      ((numerator_1 + numerator_2 + numerator_3 + numerator_4) /
       (denominator_1 + denominator_2 + denominator_3 + denominator_4));

  for (int i = 0; i < 4; ++i) {
    expected_spin_vector[i] = theta_vector[i] * total_weight;
  }

  // Perform checks
  VERIFY(
      expect_near(delta_plus.spin_vector()[0], expected_spin_vector[0], 1e-9));
  VERIFY(
      expect_near(delta_plus.spin_vector()[1], expected_spin_vector[1], 1e-9));
  VERIFY(
      expect_near(delta_plus.spin_vector()[2], expected_spin_vector[2], 1e-9));
  VERIFY(
      expect_near(delta_plus.spin_vector()[3], expected_spin_vector[3], 1e-9));
}
