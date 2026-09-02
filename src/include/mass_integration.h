
#ifndef INCLUDE_MASS_INTEGRATION_H_
#define INCLUDE_MASS_INTEGRATION_H_

#include <map>

#include "params.h"

#include "smash/particletype.h"
#include "smash/tabulation.h"

namespace MassIntegration {

namespace details {
double spectral_function(double m, const smash::ParticleType& type);
// M^2 A(M) K_2(jM/T)
double integrand(bool is_bulk, int sum_index, double m,  double T,
                 const smash::ParticleType& type);
void reset_tabulated_integrals();
} // namespace MassIntegration::details

// Calculate table of the integral over temperature for a given species and sum index.
smash::Tabulation tabulate_over_temperature(bool is_bulk, int sum_index,
                                            const smash::ParticleType& type);

// Return cached or newly tabulated value for a given temperature.
double fugacity_expansion_coefficient(bool is_bulk, int sum_index,
                                      const smash::ParticleType& type, double T);

double sample_mass(const smash::ParticleType& type);

}  // namespace MassIntegration

#endif // INCLUDE_MASS_INTEGRATION_H_
