#ifndef INCLUDE_GEN_H_
#define INCLUDE_GEN_H_

#include <array>
#include <memory>
#include <optional>

#include "smash/particles.h"
#include "smash/pdgcode.h"
#include "smash/setup_particles_decaymodes.h"
#include "vorticity.h"

class TRandom3;
class DatabasePDG2;
class Particle;

// 4×4 matrix alias
using FourMatrix = std::array<std::array<double, 4>, 4>;

namespace gen {

// creating aliases for Vorticity and energy density
using OptionalVorticity = std::optional<std::unique_ptr<Vorticity>>;
using OptionalEnergy = std::optional<double>;

// This struct is only used if vorticity_output_enabled is set to true in the
// config file. It is used to store the vorticity vector and the corresponding
// coordinates for each sampled particle.
struct ThetaStruct {
  std::array<double, 4> coordinates;
  std::array<double, 4> vorticity_vector;
};

// If vorticity_vector == 1 in the config, thetaStorage will be used to store
// the vorticity vector for each sampled particle
extern std::unique_ptr<std::vector<std::vector<ThetaStruct>>> thetaStorage;

struct element {
  double four_position[4];
  double u[4];
  double dsigma[4];
  double T, mub, muq, mus;
  double pi[10];
  double Pi;
  // optional pointer to the energy density.
  OptionalEnergy e = std::nullopt;
  // Optional pointer to the thermal vorticity tensor omega_{mu nu} following
  // the index structure {mu nu} = [{0 0}, {0 1}, {0 2}, {0 3}, {1 0}, ... ]
  OptionalVorticity vorticity = std::nullopt;
};

// typedef std::vector<Particle*> ParticleList ; // TODO in far future
//  data
extern TRandom3 *rnd;
extern smash::ParticleData ***pList;  // particle arrays
extern int *npart;
const int NPartBuf = 30000;  // dimension of particle buffer for each event

// Core implementation on ParticleType to calculate the full chemical potential
// for a given particle and freezeout element
inline double chemical_potential(const smash::ParticleType &type,
                                 const gen::element &freezeout_element) {
  return type.baryon_number() * freezeout_element.mub +
         type.charge() * freezeout_element.muq +
         type.strangeness() * freezeout_element.mus;
}

// Convenience overload for ParticleData*
inline double chemical_potential(const smash::ParticleData *particle,
                                 const gen::element &freezeout_element) {
  // debug check
  assert(particle != nullptr);
  return chemical_potential(particle->type(), freezeout_element);
}

// active Lorentz boost
void fillBoostMatrix(double vx, double vy, double vz, double boostMatrix[4][4]);

// Convert a Lorentz boost matrix from contravariant form (Λ^μ_ν, acts on upper
// indices) to covariant form (Λ_μ^ν, acts on lower indices) using the (+,-,-,-)
// metric.
//
// Index relation:   (Λ_cov)_μ^ν = g_{μα} (Λ_contra)^α_β g^{βν}
// Matrix relation:  Λ_cov = g * Λ_contra * g,  with g = diag(+1,-1,-1,-1)
FourMatrix create_covariant_boost_matrix(
    const double boost_contravariant[4][4]);
void generate();
void load(const char *filename, int N);
double ffthermal(double *x, double *par);
int index44(const int &i, const int &j);
// Allocate memory for the vorticity vector for each sampled particle
void enable_vorticity_storage();

void generate();
smash::ParticleData *acceptParticle(int event,
                                    const smash::ParticleTypePtr &ldef,
                                    smash::FourVector position,
                                    smash::FourVector momentum);
}  // namespace gen

#endif  // INCLUDE_GEN_H_
