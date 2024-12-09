#ifndef INCLUDE_SPIN_H_
#define INCLUDE_SPIN_H_

#include <array>

#include "gen.h"
#include "smash/particles.h"

namespace spin {
/* The following definitions correspond to the conventions used in Eq. (60)
 * of the paper "Exact spin polarization of massive and massless particles in
 * relativistic fluids at global equilibrium" by A. Palermo and F. Becattini
 * (arXiv:2304.02276v2).
 */
void calculate_and_set_spin_vector(const gen::element &freezeout_element,
                                   int index_event,
                                   smash::ParticleData ***particle_list,
                                   int *npart);

}  // namespace spin
#endif  // INCLUDE_SPIN_H_