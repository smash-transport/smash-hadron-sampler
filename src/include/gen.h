#ifndef INCLUDE_GEN_H_
#define INCLUDE_GEN_H_

#include <array>
#include <memory>
#include <optional>

#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"
#include "vorticity.h"

class TRandom3;
class DatabasePDG2;
class Particle;

int index44(const int &i, const int &j);

namespace gen {
// typedef std::vector<Particle*> ParticleList ; // TODO in far future
//  data
extern TRandom3 *rnd;
extern smash::ParticleData ***pList;  // particle arrays
extern int *npart;
const int NPartBuf = 30000;  // dimension of particle buffer for each event

// creating aliases for Vorticity and energy density
using OptionalVorticity = std::optional<std::unique_ptr<Vorticity>>;
using OptionalEnergy = std::optional<double>;

// Define the structure of the elements of the freeze-out surface
struct element {
  // Milne coorinates
  double tau, x, y, eta;
  // fluid 4-velocity
  double u[4];
  // normal vector to the freezeout cell
  double dsigma[4];
  // temperature and chemical potentials
  double T, mub, muq, mus;
  // shear stress tensor
  double pi[10];
  // bulk pressure
  double Pi;
  // optional pointer to the energy density.
  OptionalEnergy e = std::nullopt;
  // Optional pointer to the thermal vorticity tensor omega_{mu nu} following
  // the index structure {mu nu} = [{0 0}, {0 1}, {0 2}, {0 3}, {1 0}, ... ]
  OptionalVorticity vorticity = std::nullopt;
};

// functions
void load(char *filename, int N);
int generate();
void acceptParticle(int event, const smash::ParticleTypePtr &ldef,
                    smash::FourVector position, smash::FourVector momentum);
}  // namespace gen

#endif  // INCLUDE_GEN_H_
