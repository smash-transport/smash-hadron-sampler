#ifndef INCLUDE_GEN_H_
#define INCLUDE_GEN_H_

#include <array>
#include <memory>
#include <optional>
#include <string>

#include "smash/particles.h"
#include "smash/pdgcode.h"
#include "smash/setup_particles_decaymodes.h"
#include "vorticity.h"

class TRandom3;
class TF1;
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

// Summary produced by momentum generation diagnostics.
struct MomentumGenerationDiagnostics {
  double expected_density = 0.0;
  double generated_newvisc_density = 0.0;
  double generated_oldvisc_density = 0.0;
  double generated_equilibrium_density = 0.0;
  double acceptance_newvisc = 0.0;
  double acceptance_oldvisc = 0.0;
  double acceptance_equilibrium = 0.0;
  double chi2_newvisc = 0.0;
  double chi2_oldvisc = 0.0;
  double chi2_equilibrium = 0.0;
  int ndf_newvisc = 0;
  int ndf_oldvisc = 0;
  int ndf_equilibrium = 0;
  double pvalue_newvisc = 0.0;
  bool passed_newvisc_test = false;
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
extern element *surf;
extern double *cumulantDensity;
extern double totalDensity;
extern const smash::ParticleTypeList *database;
extern TF1 *fthermal;
const int NPartBuf = 2000000;  // dimension of particle buffer for each event
extern double dvMax, dsigmaMax;
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

double* calculate_particle_densities(element &surf_element, const std::vector<smash::ParticleType> *db);
bool generate_particle(int iel, int ievent, double dvEff);

std::tuple<double, double, double> sample_momentum_newvisc(int iel, double mass, double muf, double stat);
std::tuple<double, double, double> sample_momentum_equilibrium(int iel, double mass, double muf, double stat);
std::tuple<double, double, double> sample_momentum_oldvisc(int iel, double mass, double muf, double stat);
std::tuple<double, double, double> sample_momentum_newvisc_fullrejection(int iel, double mass, double muf, double stat);
std::tuple<double, double, double> sample_momentum(int iel, double mass, double muf, double stat, bool random_angles);

MomentumGenerationDiagnostics run_momentum_generation_diagnostics(
  int num_samples = 10000, bool save_root_output = true,
  const std::string &output_filename = "pion_momentum_test.root");


double W_shear_correction(double momArray[4], double muf, double stat, const element &surf_elem);
double W_bulk_correction(double p, double mass, double muf, double stat, const element &surf_elem);

void generate();
smash::ParticleData *acceptParticle(int event,
                                    const smash::ParticleTypePtr &ldef,
                                    smash::FourVector position,
                                    smash::FourVector momentum);
}  // namespace gen

#endif  // INCLUDE_GEN_H_
