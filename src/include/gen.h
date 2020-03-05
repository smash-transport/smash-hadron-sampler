#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

class TRandom3 ;
class DatabasePDG2;
class Particle ;

namespace gen{
//typedef std::vector<Particle*> ParticleList ; // TODO in far future
// data
extern TRandom3 *rnd ;
extern Particle ***pList ; // particle arrays
extern int *npart ;
const int NPartBuf = 15000; // dimension of particle buffer for each event

// functions
void load(char *filename, int N) ;
int generate() ;
}
