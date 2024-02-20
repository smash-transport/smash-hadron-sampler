#ifndef INCLUDE_GEN_H_
#define INCLUDE_GEN_H_

#include "smash/particles.h"
#include "smash/setup_particles_decaymodes.h"

class TRandom3;
class DatabasePDG2;
class Particle;

int index44(const int &i, const int &j);

namespace gen{
//typedef std::vector<Particle*> ParticleList ; // TODO in far future
// data
extern TRandom3 *rnd ;
extern smash::ParticleData ***pList ; // particle arrays
extern int *npart ;
const int NPartBuf = 30000; // dimension of particle buffer for each event

struct MinMax {
  double minimum;
  double maximum;
};


// functions
void load(char *filename, int N, MinMax &min_max_vorticity);
int generate();
void acceptParticle(int event, const smash::ParticleTypePtr &ldef,
                    smash::FourVector position, smash::FourVector momentum);

/*
 * Get the most likely spin projection in one cell in multiples of 1/2, based on
 * it's vorticity
 * \return Favored spin projection in multiples of 1/2
 */
int get_favored_spin_projection_in_cell(double vorticity_cell,
                                        MinMax &min_max_vorticity);

/*
 * Calculate the z projection of the vorticity in the given cell.
 * This is used to sample spin projections of particles.
 */
double get_vorticity_z_projection_in_cell(double (&u)[4],
                                          double (&u_derivatives)[16]);
/*
 * Updates the min and max values in min_max_vorticity with every loop over each
 * element, such that min_max_vorticity stores the absolute minimum and maximum
 * of the z projection of the vorticity over the whole surface in the end
 */
void update_min_max_vorticity_values(int iterator, double vorticity_cell,
                                     MinMax &min_max_vorticity);

/*
 * Sample the spin projection of a particle based on the vorticity (z
 * projection) of the corresponding cell. First the most likely spin state
 * is determined, afterwards the actual spin state is sampled with a
 * Monte-Carlo sampling, with the acceptance range adjusted according to the
 * most likely spin projection determined before.
 */
int get_sampled_spin_projection();
}  // namespace gen

#endif  // INCLUDE_GEN_H_
