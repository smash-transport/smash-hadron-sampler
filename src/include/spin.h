#ifndef INCLUDE_SPIN_H_
#define INCLUDE_SPIN_H_

/* The following definitions correspond to the conventions used in Eq. (60)
 * of the paper "Exact spin polarization of massive and massless particles in
 * relativistic fluids at global equilibrium" by A. Palermo and F. Becattini
 * (arXiv:2304.02276v2). 
 */
std::array<double, 4> spin_vector(const gen::element& freezeout_element,
                                  const smash::ParticleData& particle);
#endif  // INCLUDE_SPIN_H_