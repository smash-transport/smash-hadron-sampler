#ifndef INCLUDE_GEN_H_
#define INCLUDE_GEN_H_

#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

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
  // energy density
  double e;
  // shear stress tensor
  double pi[10];
  // bulk pressure
  double Pi;
  // thermal vorticity tensor omega_{mu nu} following the index structure
  // {mu nu} = [{0 0}, {0 1}, {0 2}, {0 3}, {1 0}, ..., {3 3}]
  double vorticity[16];
};

// functions
void load(char *filename, int N);
int generate();
void acceptParticle(int event, const smash::ParticleTypePtr &ldef,
                    smash::FourVector position, smash::FourVector momentum,
                    double vorticity_cell);
/*
 * If spin sampling is enabled, we need the energy ensity of every cell, which
 * is stored in the extended freezeout surface format. This function checks if
 * the extended freezeout surface is used and throws an exception if not.
 */
void ensure_extended_freezeout_is_used(const std::string &filename);
}  // namespace gen

#endif  // INCLUDE_GEN_H_
