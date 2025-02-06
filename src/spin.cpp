#include "spin.h"

#include <cmath>
#include <vector>

#include "const.h"
#include "gen.h"
#include "smash/particles.h"
#include "vorticity.h"

namespace spin {

/* The following definitions correspond to the conventions used in Eq. (60)
 * of the paper "Exact spin polarization of massive and massless particles in
 * relativistic fluids at global equilibrium" by A. Palermo and F. Becattini
 * (arXiv:2304.02276v2).
 */

// Calculate theta from EQ. 42 in arXiv:2304.02276v2 which is needed to
// calculate the full EQ. 60 afterwards. The momentum needs to be given with
// upper indices (as calculated by the sampler), the function will consider
std::array<double, 4> theta(const double mass,
                            const std::array<double, 16> &vorticity,
                            const double (&p)[4]) {
  if (mass < small_value) {
    throw std::invalid_argument(
        "Theta cannot be calculated for massless particles.");
  }
  // Transform the row and column indices for a 4x4 matrix into a 1d index
  auto at = [](int i, int j) { return (i * 4 + j); };

  // Here, we explicitely performed the contractions with the Levi-Civita tensor
  // assuming that the vorticity tensor is antisymmetric (which canceled the
  // initial factor of 1/2).
  const double theta0 =
      (-1.0 / mass) * (vorticity[at(2, 3)] * p[1] + vorticity[at(3, 1)] * p[2] +
                       vorticity[at(1, 2)] * p[3]);

  const double theta1 =
      (-1.0 / mass) * (vorticity[at(3, 2)] * p[0] + vorticity[at(0, 3)] * p[2] +
                       vorticity[at(2, 0)] * p[3]);

  const double theta2 =
      (-1.0 / mass) * (vorticity[at(1, 3)] * p[0] + vorticity[at(3, 0)] * p[1] +
                       vorticity[at(0, 1)] * p[3]);

  const double theta3 =
      (-1.0 / mass) * (vorticity[at(2, 1)] * p[0] + vorticity[at(0, 2)] * p[1] +
                       vorticity[at(1, 0)] * p[2]);

  // The additional minus sign is due to the position of the indicies of the
  // vorticity tensor and the momentum which all have lower indices. As we only
  // get quantities with upper indices from the sampler, we need to lower these
  // indices with three metric tensors. Calculating the 0 component, the product
  // of these three metric tensors gives a factor of (-1)^3 = -1, while for the
  // spatial components, the product gives a factor of (-1)^2 = 1.
  return {-theta0, theta1, theta2, theta3};
}

// Calculate and set the spin vector for a given particle
void calculate_and_set_spin_vector(const gen::element &freezeout_element,
                                   smash::ParticleData *particle) {
  // Ensure that the optional values in the freezeout element are set
  if (!freezeout_element.vorticity.has_value()) {
    throw std::runtime_error("Vorticity tensor not set in surface element.");
  }
  if (!freezeout_element.e.has_value()) {
    throw std::runtime_error("Energy density not set in surface element.");
  }
  const double tiny_value = 1e-8;
  const int spin = particle->spin();

  if (spin == 0) {
    const smash::FourVector spin_vec(0.0, 0.0, 0.0, 0.0);
    particle->set_spin_vector(spin_vec);
  } else if (spin > 0) {
    const double energy_density = *freezeout_element.e;
    const double temperature = freezeout_element.T;
    const double mu = freezeout_element.mub;
    const std::array<double, 16> vorticity =
        (**freezeout_element.vorticity).get_vorticity();

    double p[4] = {particle->momentum().x0(), particle->momentum().x1(),
                   particle->momentum().x2(), particle->momentum().x3()};

    const std::array<double, 4> theta_array =
        theta(particle->pole_mass(), vorticity, p);

    const double theta_squared = four_vector_square(theta_array);

    if (theta_squared > 0) {
      throw std::runtime_error(
          "theta_squared must be negative for valid spin vector calculation.");
    }

    // If \sqrt{-theta^2} is sufficiently small, we will use the approximation
    // given right after Eq. (61) in arXiv:2304.02276v2. This prevents division
    // by zero for corona cells in which the vorticity tensor is zero.
    if (std::abs(std::sqrt(-theta_squared)) < tiny_value) {
      const double distribution_argument = (energy_density - mu) / temperature;
      const double factor =
          ((spin / 2.) * ((spin / 2.) + 1.)) / 3. +
          (1 + (spin % 2 == 0 ? 1 : -1) *
                   fermi_bose_distribution(spin, distribution_argument));

      smash::FourVector spin_vec(0.0, 0.0, 0.0, 0.0);
      for (int i = 0; i < 4; i++) {
        spin_vec[i] = theta_array[i] * factor;
      }
      // Set the spin vector in the particle data
      particle->set_spin_vector(spin_vec);

    } else {
      double numerator = 0.0;
      double denominator = 0.0;
      // Sum all terms of in the numerator and denominator
      for (int k = -spin; k <= spin; k += 2) {
        double sum_index = k / 2.0;
        double exponential = std::exp(exponent(sum_index, energy_density,
                                               temperature, mu, theta_squared));
        double denominator_term = 1 / (exponential - (spin % 2 == 0 ? 1 : -1));
        double numerator_term = sum_index * denominator_term;

        numerator += numerator_term;
        denominator += denominator_term;
      }

      if (std::abs(denominator) < 1e-8) {
        throw std::runtime_error(
            "Denominator must be sufficiently negative for valid spin vector "
            "calculation.");
      }

      // Calculate the spin vector
      std::array<double, 4> spin_vector;
      for (int i = 0; i < 4; i++) {
        spin_vector[i] = (theta_array[i] * numerator) /
                         (std::sqrt(-theta_squared) * denominator);
      }

      const smash::FourVector spin_vec(spin_vector[0], spin_vector[1],
                                       spin_vector[2], spin_vector[3]);

      // Set the spin vector in the particle data
      particle->set_spin_vector(spin_vec);
    }
  } else {
    throw std::runtime_error("Spin of particle is invalid or unset.");
  }
}

}  // namespace spin