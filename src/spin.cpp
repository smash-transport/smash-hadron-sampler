#include "spin.h"

#include <array>
#include <cmath>
#include <stdexcept>
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

// Calculate theta in EQ. (60) from arXiv:2304.02276v2
static std::array<double, 4> theta(const double mass,
                                   const std::array<double, 16> &vorticity,
                                   const double (&p_lower_index)[4]) {
  if (mass < small_value) {
    throw std::invalid_argument(
        "Theta cannot be calculated for massless particles.");
  }
  // Transform the row and column indices for a 4x4 matrix into a 1d index
  auto at = [](int i, int j) { return (i * 4 + j); };

  // Here, we explicitely performed the contractions with the Levi-Civita tensor
  const double theta0 =
      (-1.0 / mass) * (vorticity[at(2, 3)] * p_lower_index[1] +
                       vorticity[at(3, 1)] * p_lower_index[2] +
                       vorticity[at(1, 2)] * p_lower_index[3]);

  const double theta1 =
      (-1.0 / mass) * (vorticity[at(3, 2)] * p_lower_index[0] +
                       vorticity[at(0, 3)] * p_lower_index[2] +
                       vorticity[at(2, 0)] * p_lower_index[3]);

  const double theta2 =
      (-1.0 / mass) * (vorticity[at(1, 3)] * p_lower_index[0] +
                       vorticity[at(3, 0)] * p_lower_index[1] +
                       vorticity[at(0, 1)] * p_lower_index[3]);

  const double theta3 =
      (-1.0 / mass) * (vorticity[at(2, 1)] * p_lower_index[0] +
                       vorticity[at(0, 2)] * p_lower_index[1] +
                       vorticity[at(1, 0)] * p_lower_index[2]);
  return {theta0, theta1, theta2, theta3};
}

// Calculate the square of a four-vector in Minkowski space
static double four_vector_square(std::array<double, 4> &four_vector) {
  double result = 0.0;
  for (int i = 0; i < 4; i++) {
    result += metric[i] * four_vector[i] * four_vector[i];
  }
  return result;
}

// Calculate the exponent of exp(...) in EQ. (60) from arXiv:2304.02276v2
static double exponent(const int k, const double energy_density,
                       const double temperature, const double mu,
                       const double theta_squared) {
  return (energy_density - mu) / temperature - k * sqrt(-theta_squared);
}

void calculate_and_set_spin_vector(const gen::element &freezeout_element,
                                   int index_event,
                                   smash::ParticleData ***particle_list,
                                   int *npart) {
  // Ensure that the optional values in the freezeout element are set
  if (!freezeout_element.e || !freezeout_element.vorticity) {
    throw std::runtime_error(
        "Energy density or vorticity tensor is not set in the freezeout "
        "element.");
  }
  // As gen::acceptParticle() increases the particle number in the corresponding
  // event by 1, we need to subtract 1 to get the correct index of the particle.
  const int index_particle = npart[index_event] - 1;
  smash::ParticleData *particle = particle_list[index_event][index_particle];
  const int spin = particle->spin();

  if (spin == 0) {
    const smash::FourVector spin_vec(0, 0, 0, 0);
    particle->set_spin_vector(spin_vec);
  } else if (spin != 0) {
    const double energy_density = *freezeout_element.e;
    const double temperature = freezeout_element.T;
    const double mu = freezeout_element.mub;
    const std::array<double, 16> vorticity =
        (**freezeout_element.vorticity).get_vorticity();
    // Momentum with lower index
    double p_[4];
    smash::FourVector four_momentum = particle->momentum();
    for (int i = 0; i < 4; i++) {
      p_[i] = four_momentum[i] * metric[i];
    }
    std::array<double, 4> theta_array =
        theta(particle->pole_mass(), vorticity, p_);
    const double theta_squared = four_vector_square(theta_array);

    // Calculate the numerator and denominator parts with the spin sums
    const int spin_degeneracy = spin + 1;

    std::vector<double> numerator_vector;
    std::vector<double> denominator_vector;

    int array_index = 0;
    // Calculate all terms of the sums in the numerator and denominator
    // separately and store them in arrays
    for (int k = -spin; k <= spin; k++) {
      // Ensure that the loop index is within the array bounds
      if (array_index < 0 || array_index >= spin_degeneracy) {
        throw std::out_of_range("Array index is out of range.");
      }
      double exponential =
          std::exp(exponent(k, energy_density, temperature, mu, theta_squared));
      denominator_vector.push_back(1 /
                                   (exponential - (spin % 2 == 0 ? 1 : -1)));
      numerator_vector.push_back(k * denominator_vector[array_index]);
      array_index++;
    }
    // Performing the sum over k in the numerator and denominator
    const double numerator =
        std::accumulate(numerator_vector.begin(), numerator_vector.end(), 0.0);
    const double denominator = std::accumulate(denominator_vector.begin(),
                                               denominator_vector.end(), 0.0);

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
  } else {
    throw std::runtime_error("Spin of particle is invalid or unset.");
  }
}

}  // namespace spin