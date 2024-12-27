#ifndef INCLUDE_SPIN_H_
#define INCLUDE_SPIN_H_

#include <array>
#include <stdexcept>

#include "const.h"
#include "gen.h"
#include "smash/particles.h"

namespace spin {
/* The following definitions correspond to the conventions used in Eq. (60)
 * of the paper "Exact spin polarization of massive and massless particles in
 * relativistic fluids at global equilibrium" by A. Palermo and F. Becattini
 * (arXiv:2304.02276v2).
 */

// Calculate the square of a four-vector in Minkowski space
inline double four_vector_square(std::array<double, 4> &four_vector) {
  double result = 0.0;
  for (int i = 0; i < 4; i++) {
    result += metric[i] * four_vector[i] * four_vector[i];
  }
  return result;
}

// Calculate theta from EQ. 42 in arXiv:2304.02276v2 which is needed to
// calculate the full EQ. 60 afterwards.
std::array<double, 4> theta(const double mass,
                            const std::array<double, 16> &vorticity,
                            const double (&p_lower_index)[4]);

// Calculate the exponent of exp(...) in EQ. (60) from arXiv:2304.02276v2
inline double exponent(const double k, const double energy_density,
                       const double temperature, const double mu,
                       const double theta_squared) {
  // Check that k is a multiple of 1/2
  if (std::abs(std::fmod(k, 0.5)) > 1e-6) {
    throw std::invalid_argument("k must be a multiple of 1/2.");
  } else if (theta_squared > 0) {
    throw std::invalid_argument("theta^2 must be negative.");
  }
  return (energy_density - mu) / temperature - k * sqrt(-theta_squared);
}

void calculate_and_set_spin_vector(const gen::element &freezeout_element,
                                   int index_event,
                                   smash::ParticleData ***particle_list,
                                   int *npart);

}  // namespace spin
#endif  // INCLUDE_SPIN_H_