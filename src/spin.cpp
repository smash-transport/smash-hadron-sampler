#include "spin.h"

#include <cmath>
#include <vector>

#include "const.h"
#include "gen.h"
#include "params.h"
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
// upper indices (as calculated by the sampler).
std::array<double, 4> theta(const double mass,
                            const std::array<double, 16> &vorticity,
                            const std::array<double, 4> &p) {
  if (mass < small_value) {
    throw std::invalid_argument(
        "Theta cannot be calculated for massless particles.");
  }
  // Lower index of the momentum (vorticity already has lower indices from
  // vHLLE)
  const std::array<double, 4> p_ = {p[0], -p[1], -p[2], -p[3]};

  // Transform the row and column indices for a 4x4 matrix into a 1d index
  auto at = [](int i, int j) { return (i * 4 + j); };

  // Here, we explicitely performed the contractions with the Levi-Civita tensor
  // assuming that the vorticity tensor is antisymmetric (which canceled the
  // initial factor of 1/2).
  // The factor of hbarC is added to ensure that theta is dimensionless as vHLLE
  // calculates it in 1/fm which needs to be compensated.
  const double theta0 =
      hbarC * (-1.0 / mass) *
      (vorticity[at(2, 3)] * p_[1] + vorticity[at(3, 1)] * p_[2] +
       vorticity[at(1, 2)] * p_[3]);

  const double theta1 =
      hbarC * (-1.0 / mass) *
      (vorticity[at(3, 2)] * p_[0] + vorticity[at(0, 3)] * p_[2] +
       vorticity[at(2, 0)] * p_[3]);

  const double theta2 =
      hbarC * (-1.0 / mass) *
      (vorticity[at(1, 3)] * p_[0] + vorticity[at(3, 0)] * p_[1] +
       vorticity[at(0, 1)] * p_[3]);

  const double theta3 =
      hbarC * (-1.0 / mass) *
      (vorticity[at(2, 1)] * p_[0] + vorticity[at(0, 2)] * p_[1] +
       vorticity[at(1, 0)] * p_[2]);

  return {theta0, theta1, theta2, theta3};
}

void add_entry_to_theta_storage(const int index_event,
                                const smash::ParticleData *particle,
                                const std::array<double, 16> &vorticity) {
  const std::array<double, 4> p = {
      particle->momentum().x0(), particle->momentum().x1(),
      particle->momentum().x2(), particle->momentum().x3()};

  const std::array<double, 4> txyz = {
      particle->position().x0(), particle->position().x1(),
      particle->position().x2(), particle->position().x3()};

  const std::array<double, 4> theta_array =
      theta(particle->pole_mass(), vorticity, p);

  (*gen::thetaStorage)[index_event].push_back({txyz, theta_array});
}

// Calculate and set the spin vector for a given particle
void calculate_and_set_spin_vector(const int index_event,
                                   const gen::element &freezeout_element,
                                   smash::ParticleData *particle) {
  // Ensure that the optional values in the freezeout element are set
  if (!freezeout_element.vorticity.has_value()) {
    throw std::runtime_error("Vorticity tensor not set in surface element.");
  }
  if (!freezeout_element.e.has_value()) {
    throw std::runtime_error("Energy density not set in surface element.");
  }
  // Linearity applies already for values smaller than 1e-3
  const double tiny_value = 1e-3;
  const int spin = particle->spin();

  if (spin == 0) {
    // Store the vorticity vector if vorticity_output is enabled
    if (params::vorticity_output_enabled) {
      const std::array<double, 16> vorticity =
          (**freezeout_element.vorticity).get_vorticity();
      add_entry_to_theta_storage(index_event, particle, vorticity);
    }

    const smash::FourVector spin_vec(0.0, 0.0, 0.0, 0.0);
    particle->set_spin_vector(spin_vec);

  } else if (spin > 0) {
    const double temperature = freezeout_element.T;
    const double mu = gen::chemical_potential(particle, freezeout_element);
    const std::array<double, 16> vorticity =
        (**freezeout_element.vorticity).get_vorticity();

    // Store the vorticity vector if vorticity_output is enabled
    if (params::vorticity_output_enabled) {
      add_entry_to_theta_storage(index_event, particle, vorticity);
    }

    std::array<double, 4> p = {
        particle->momentum().x0(), particle->momentum().x1(),
        particle->momentum().x2(), particle->momentum().x3()};

    const double particle_energy = particle->momentum().x0();

    const std::array<double, 4> theta_array =
        theta(particle->pole_mass(), vorticity, p);
    double theta_squared = four_vector_square(theta_array);

    if (theta_squared > 0) {
      std::cerr << "Warning: theta_squared = " << theta_squared
                << " is positive, setting theta_squared = 0.\n";
      theta_squared = -small_value;
    }

    // If \sqrt{-theta^2} is sufficiently small, we will use the approximation
    // given right after Eq. (61) in arXiv:2304.02276v2. This prevents division
    // by zero for corona cells in which the vorticity tensor is zero.
    if (std::abs(std::sqrt(-theta_squared)) < tiny_value) {
      const double distribution_argument = (particle_energy - mu) / temperature;
      const double factor =
          (((spin / 2.) * ((spin / 2.) + 1.)) / 3.) *
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
        double exponential = std::exp(exponent(sum_index, particle_energy,
                                               temperature, mu, theta_squared));
        double denominator_term = 1 / (exponential - (spin % 2 == 0 ? 1 : -1));
        double numerator_term = sum_index * denominator_term;

        numerator += numerator_term;
        denominator += denominator_term;
      }

      if (denominator == 0.0) {
        throw std::runtime_error(
            "Denominator in spin vector calculation is zero.");
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
