#ifndef INCLUDE_GEN_H_
#define INCLUDE_GEN_H_

#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

class TRandom3;
class DatabasePDG2;
class Particle;

namespace gen {
// typedef std::vector<Particle*> ParticleList ; // TODO in far future
//  data
extern TRandom3 *rnd;
extern smash::ParticleData ***pList; // particle arrays
extern int *npart;
const int NPartBuf = 30000; // dimension of particle buffer for each event

// functions
void fillBoostMatrix(double vx, double vy, double vz, double boostMatrix[4][4]);
void generate();
void load(const char *filename, int N);
double ffthermal(double *x, double *par);
int index44(const int &i, const int &j);
} // namespace gen

#endif // INCLUDE_GEN_H_
