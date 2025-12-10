#include "gen.h"

#include <ctime>
#include <limits>
#include <numeric>
#include <random>
#include <vector>
#include <mutex>

#include "virtest/vir/test.h"

using namespace gen;

// Custom function to compare two doubles within a given tolerance
static bool expect_near(double val1, double val2, double abs_error) {
  return std::abs(val1 - val2) <= abs_error;
}

// Ensures the global ParticleType list is initialized only once using
// std::call_once. Avoids the "Type list was already built!" exception when
// multiple tests need it.
static void ensure_particletype_initialized() {
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

TEST(chemical_potential) {
  // Initialize the particle type list with a pion (spin 0),
  // proton (spin 1/2), and delta plus (spin 3/2)
  ensure_particletype_initialized();

  // Define the PDG code for the particles
  const smash::PdgCode pdg_pion = 0x111;
  const smash::PdgCode pdg_proton = 0x2212;
  const smash::PdgCode pdg_rho = 0x213;
  const smash::PdgCode pdg_lambda = 0x3122;

  // Create the particle data for the pion, proton, and delta plus
  smash::ParticleData pion =
      smash::ParticleData(smash::ParticleType::find(pdg_pion));
  smash::ParticleData proton =
      smash::ParticleData(smash::ParticleType::find(pdg_proton));
  smash::ParticleData rho =
      smash::ParticleData(smash::ParticleType::find(pdg_rho));
  smash::ParticleData lambda =
      smash::ParticleData(smash::ParticleType::find(pdg_lambda));

  // Set all values of the surface element to 1.0 which
  // do not affect the spin sampling
  gen::element surf_element;
  surf_element.mub = 0.9;
  surf_element.muq = 3.1;
  surf_element.mus = 7.5;

  // Calculate expected chemical potentials
  double expected_chemical_potential_pion = 0.0;
  double expected_chemical_potential_proton = 1.0 * 0.9 + 1.0 * 3.1;
  double expected_chemical_potential_rho = 1.0 * 3.1;
  double expected_chemical_potential_lambda = 1.0 * 0.9 - 1.0 * 7.5;

  // Chemical potentials from function
  double chemical_potential_pion = chemical_potential(&pion, surf_element);
  double chemical_potential_proton = chemical_potential(&proton, surf_element);
  double chemical_potential_rho = chemical_potential(&rho, surf_element);
  double chemical_potential_lambda = chemical_potential(&lambda, surf_element);

  // Perform checks
  VERIFY(expect_near(chemical_potential_pion, expected_chemical_potential_pion,
                     1e-12));
  VERIFY(expect_near(chemical_potential_proton,
                     expected_chemical_potential_proton, 1e-12));
  VERIFY(expect_near(chemical_potential_rho, expected_chemical_potential_rho,
                     1e-12));
  VERIFY(expect_near(chemical_potential_lambda,
                     expected_chemical_potential_lambda, 1e-12));
}

TEST(chemical_potential_from_type_overload) {
  // Ensure particle types are available
  ensure_particletype_initialized();

  // PDG codes
  const smash::PdgCode pdg_pion = 0x111;
  const smash::PdgCode pdg_proton = 0x2212;
  const smash::PdgCode pdg_rho = 0x213;
  const smash::PdgCode pdg_lambda = 0x3122;

  // Build ParticleData so we can easily access type() (const ParticleType&)
  smash::ParticleData pion =
      smash::ParticleData(smash::ParticleType::find(pdg_pion));
  smash::ParticleData proton =
      smash::ParticleData(smash::ParticleType::find(pdg_proton));
  smash::ParticleData rho =
      smash::ParticleData(smash::ParticleType::find(pdg_rho));
  smash::ParticleData lambda =
      smash::ParticleData(smash::ParticleType::find(pdg_lambda));

  // Surface element with chemical potentials
  gen::element surf_element;
  surf_element.mub = 0.9;
  surf_element.muq = 3.1;
  surf_element.mus = 7.5;

  // Expected values (B, Q, S factors)
  const double exp_pion = 0.0;          // B=0, Q=0, S=0
  const double exp_proton = 0.9 + 3.1;  // B=1, Q=1, S=0
  const double exp_rho = 3.1;           // B=0, Q=1, S=0
  const double exp_lambda = 0.9 - 7.5;  // B=1, Q=0, S=-1

  // Call the ParticleType& overload via type()
  const double mu_pion_type =
      gen::chemical_potential(pion.type(), surf_element);
  const double mu_proton_type =
      gen::chemical_potential(proton.type(), surf_element);
  const double mu_rho_type = gen::chemical_potential(rho.type(), surf_element);
  const double mu_lambda_type =
      gen::chemical_potential(lambda.type(), surf_element);

  // Independent checks against expected numbers
  VERIFY(expect_near(mu_pion_type, exp_pion, 1e-12));
  VERIFY(expect_near(mu_proton_type, exp_proton, 1e-12));
  VERIFY(expect_near(mu_rho_type, exp_rho, 1e-12));
  VERIFY(expect_near(mu_lambda_type, exp_lambda, 1e-12));

  // Cross-check: type overload matches pointer overload (one is enough, do
  // proton)
  const double mu_proton_ptr = gen::chemical_potential(&proton, surf_element);
  VERIFY(expect_near(mu_proton_ptr, mu_proton_type, 1e-12));
}
