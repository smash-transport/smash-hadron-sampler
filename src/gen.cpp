#include "gen.h"

#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TROOT.h>
#include <TRandom3.h>

#include <cmath>
#include <fstream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>

#include "const.h"
#include "params.h"
#include "spin.h"
#include "vorticity.h"

using namespace std;

// ##########################################################
// #  this version works with arbitrary T/mu distribution   #
// #  on freezeout hypersurface (May'2012)                  #
// #  also, pre-BOOSTED dsigma is used                      #
// ##########################################################

namespace gen {

int Nelem;
int NPART;
double *ntherm, dvMax, dsigmaMax;
TRandom3 *rnd;
smash::ParticleData ***pList;  // particle arrays
std::unique_ptr<std::vector<std::vector<ThetaStruct>>> thetaStorage = nullptr;

element *surf;
int *npart;               // number of generated particles in each event
double *cumulantDensity;  // particle densities (thermal). Seems to be
                          // redundant, but needed for fast generation
double totalDensity;      // sum of all thermal densities

std::vector<smash::PdgCode> species_to_exclude;  // species to exclude from sampling
const smash::ParticleTypeList *database = nullptr;  // Hadron database, initialized in load()
TF1 *fthermal = nullptr;  // Thermal momentum distribution, initialized in generate()
TLorentzVector mom;  // Momentum 4-vector, reused in generate()
const double gmumu[4] = {1., -1., -1., -1.};
int nmaxiter = 0;  // for debug purposes, to check the maximum number of iterations in the momentum generation loop

// active Lorentz boost
void fillBoostMatrix(double vx, double vy, double vz, double boostMatrix[4][4])
// here in boostMatrix [0]=t, [1]=x, [2]=y, [3]=z
{
  const double vv[3] = {vx, vy, vz};
  const double v2 = vx * vx + vy * vy + vz * vz;
  const double gamma = 1.0 / sqrt(1.0 - v2);
  if (std::isinf(gamma) || std::isnan(gamma)) {
    std::cout << "boost vector invalid; exiting\n";
    exit(1);
  }
  boostMatrix[0][0] = gamma;
  boostMatrix[0][1] = boostMatrix[1][0] = vx * gamma;
  boostMatrix[0][2] = boostMatrix[2][0] = vy * gamma;
  boostMatrix[0][3] = boostMatrix[3][0] = vz * gamma;
  if (v2 > 0.0) {
    for (int i = 1; i < 4; i++)
      for (int j = 1; j < 4; j++)
        boostMatrix[i][j] = (gamma - 1.0) * vv[i - 1] * vv[j - 1] / v2;
  } else {
    for (int i = 1; i < 4; i++)
      for (int j = 1; j < 4; j++) boostMatrix[i][j] = 0.0;
  }
  for (int i = 1; i < 4; i++) boostMatrix[i][i] += 1.0;
}

// index44: returns an index of pi^{mu nu} mu,nu component in a plain 1D array
int index44(const int &i, const int &j) {
  if (i > 3 || j > 3 || i < 0 || j < 0) {
    throw std::out_of_range("index44: indices must be in [0, 3]");
  }
  if (j < i)
    return (i * (i + 1)) / 2 + j;
  else
    return (j * (j + 1)) / 2 + i;
}

// ######## load the elements
void load(const char *filename, int N) {
  ROOT::EnableThreadSafety();
  double vEff = 0.0, vEffOld = 0.0, dvEff, dvEffOld;
  int nfail = 0, ncut = 0;
  TLorentzVector dsigma;
  Nelem = N;
  surf = new element[Nelem];

  pList = new smash::ParticleData **[params::number_of_events];
  for (int i = 0; i < params::number_of_events; i++) {
    pList[i] = new smash::ParticleData *[NPartBuf];
  }
  npart = new int[params::number_of_events];

  std::cout << "Reading " << N << " lines from '" << filename << "'\n";
  ifstream fin(filename);
  if (!fin) {
    std::cout << "cannot read file " << filename << std::endl;
    exit(1);
  }
  dvMax = 0.;
  dsigmaMax = 0.;

  // Read the vorticity tensor from file and set it in all surface cells
  if (params::spin_vector_enabled) {
    std::cout << "Setting vorticity tensor in all surface cells from file "
              << params::vorticity_file << std::endl;
    Vorticity::set_vorticity_in_surface_cells(surf, Nelem);
  }

  // ---- reading loop
  string line;
  istringstream instream;
  std::cout << "1?: failbit=" << instream.fail() << std::endl;
  for (int n = 0; n < Nelem; n++) {
    getline(fin, line);
    instream.str(line);
    instream.seekg(0);
    instream.clear();  // does not work with gcc 4.1 otherwise
    instream >> surf[n].four_position[0] >> surf[n].four_position[1] >>
        surf[n].four_position[2] >> surf[n].four_position[3] >>
        surf[n].dsigma[0] >> surf[n].dsigma[1] >> surf[n].dsigma[2] >>
        surf[n].dsigma[3] >> surf[n].u[0] >> surf[n].u[1] >> surf[n].u[2] >>
        surf[n].u[3] >> surf[n].T >> surf[n].mub >> surf[n].muq >> surf[n].mus;
    for (int i = 0; i < 10; i++) instream >> surf[n].pi[i];
    instream >> surf[n].Pi;

    // If spin sampling is enabled, load the energy density from the
    // extended freezeout surface
    if (params::spin_vector_enabled) {
      double tmp_e;
      instream >> tmp_e;
      surf[n].e = tmp_e;  // set energy density
    }

    if (surf[n].muq > 0.12) {
      surf[n].muq = 0.12;  // omit charge ch.pot. for test
      ncut++;
    }
    if (surf[n].muq < -0.12) {
      surf[n].muq = -0.12;  // omit charge ch.pot. for test
      ncut++;
    }
    if (instream.fail()) {
      std::cout << "reading failed at line " << n << "; exiting\n";
      exit(1);
    }
    // calculate in the old way
    dvEffOld =
        surf[n].dsigma[0] * surf[n].u[0] + surf[n].dsigma[1] * surf[n].u[1] +
        surf[n].dsigma[2] * surf[n].u[2] + surf[n].dsigma[3] * surf[n].u[3];
    vEffOld += dvEffOld;
    if (dvEffOld < 0.0) {
      // cout << "!!! dvOld!=dV " << dvEffOld <<"  " << dV << "  "
      // << surf[n].four_position[0] << endl ;
      nfail++;
    }
    // if(nfail==100) exit(1) ;
    //  ---- boost
    dsigma.SetXYZT(-surf[n].dsigma[1], -surf[n].dsigma[2], -surf[n].dsigma[3],
                   surf[n].dsigma[0]);
    dsigma.Boost(-surf[n].u[1] / surf[n].u[0], -surf[n].u[2] / surf[n].u[0],
                 -surf[n].u[3] / surf[n].u[0]);
    // ######################################################################
    // ###     boost surf.dsigma to the fluid rest frame                   ##
    // ######################################################################
    surf[n].dsigma[0] = dsigma.T();
    surf[n].dsigma[1] = -dsigma.X();
    surf[n].dsigma[2] = -dsigma.Y();
    surf[n].dsigma[3] = -dsigma.Z();
    dvEff = surf[n].dsigma[0];
    vEff += dvEff;
    if (dvMax < dvEff) dvMax = dvEff;
    // maximal value of the weight max(W) = max(dsigma_0+|\vec dsigma_i|)   for
    // equilibrium DFs
    if (dsigma.T() + dsigma.Rho() > dsigmaMax)
      dsigmaMax = dsigma.T() + dsigma.Rho();
    // ########################
    // pi^{mu nu} boost to fluid rest frame
    // ########################
    double boostMatrix[4][4];
    if (params::shear_viscosity_enabled) {
      fillBoostMatrix(-surf[n].u[1] / surf[n].u[0],
                      -surf[n].u[2] / surf[n].u[0],
                      -surf[n].u[3] / surf[n].u[0], boostMatrix);
    }

    /* _pi^{μν} (upper,upper) uses the contravariant boostMatrix acting on
    upper indices */
    if (params::shear_viscosity_enabled) {
      double _pi[10];
      for (int i = 0; i < 4; i++)
        for (int j = i; j < 4; j++) {
          _pi[index44(i, j)] = 0.0;
          for (int k = 0; k < 4; k++)
            for (int l = 0; l < 4; l++)
              _pi[index44(i, j)] += surf[n].pi[index44(k, l)] *
                                    boostMatrix[i][k] * boostMatrix[j][l];
        }
      for (int i = 0; i < 10; i++) surf[n].pi[i] = _pi[i];
    }  // end pi boost
  }
  if (params::shear_viscosity_enabled)
    dsigmaMax *= 2.0;  // *2.0: jun17. default: *1.5
  else
    dsigmaMax *= 1.3;

  std::cout << "..done.\n";
  std::cout << "Veff = " << vEff << "  dvMax = " << dvMax << std::endl;
  std::cout << "Veff(old) = " << vEffOld << std::endl;
  std::cout << "failed elements = " << nfail << std::endl;
  std::cout << "mu_cut elements = " << ncut << std::endl;
  // ---- prepare some stuff to calculate thermal densities

  // Initialize the thermal momentum distribution
  fthermal = new TF1("fthermal", ffthermal, 0.0, 10.0, 4);

  // Load SMASH hadron list
  smash::initialize_default_particles_and_decaymodes();
  database = &smash::ParticleType::list_all();

  // Dump list of hadronic states to terminal
  // for (auto& HadronState: database) {
  //   std::cout << HadronState << '\n';
  // }

  // List species that should not be sampled: photon, electron, muon, tau
  // Sigma meson needs to be excluded to generate correct multiplicities
  species_to_exclude = {0x11, -0x11, 0x13, -0x13,
                        0x15, -0x15, 0x22, 0x9000221};

  // NPART = total number of hadron states
  NPART = database->size();
  std::cout << "NPART=" << NPART << std::endl;
  std::cout << "dsigmaMax=" << dsigmaMax << "\n\n";
  cumulantDensity = new double[NPART];
}

double W_bulk_correction(double p, double mass, double muf, double stat, const element &surf_elem) {
  double feq =
      C_Feq / (exp((sqrt(p * p + mass * mass) - muf) / surf_elem.T) - stat);
  return (1.0 + stat * feq) * surf_elem.Pi *
          (mass * mass / (3 * sqrt(p * p + mass * mass)) -
            sqrt(p * p + mass * mass) * (1.0 / 3.0 - params::speed_of_sound_squared)) /
          (15 * (1.0 / 3.0 - params::speed_of_sound_squared) *
            (1.0 / 3.0 - params::speed_of_sound_squared) * surf_elem.T *
            (params::ecrit +
            params::ecrit * params::ratio_pressure_energydensity));
}

double W_shear_correction(double momArray[4], double muf, double stat, const element &surf_elem) {
  double pipp = 0;
  double feq =
      C_Feq / (exp((momArray[0] - muf) / surf_elem.T) - stat);
  for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          pipp += momArray[i] * momArray[j] * gmumu[i] * gmumu[j] *
                  surf_elem.pi[index44(i, j)];
  return (1.0 + stat * feq) * pipp /
          (2. * surf_elem.T * surf_elem.T *
            (params::ecrit +
            params::ecrit * params::ratio_pressure_energydensity));
}

void enable_vorticity_storage() {
  // Allocate memory for the vorticity vector for each sampled particle
  // if spin sampling and vorticity output are enabled in the config
  if (params::spin_vector_enabled && params::vorticity_output_enabled) {
    thetaStorage = std::make_unique<std::vector<std::vector<ThetaStruct>>>(
        params::number_of_events);
  } else if (!params::spin_vector_enabled &&
             params::vorticity_output_enabled) {
    throw std::runtime_error(
        "Vorticity output is enabled but spin vector computation is not. "
        "Enable spin vector computation in the config file by adding "
        " the line 'compute_spin_vector 1'.");
  }
}

double ffthermal(double *x, double *par) {
  double &T = par[0];
  double &mu = par[1];
  double &mass = par[2];
  double &stat = par[3];
  return x[0] * x[0] / (exp((sqrt(x[0] * x[0] + mass * mass) - mu) / T) - stat);
}

double* calculate_particle_densities(int iel, const std::vector<smash::ParticleType> *db) {
    cumulantDensity = new double[db->size()];
    int ip = 0;
    for (auto &particle : *db) {
        double density = 0.;
        const bool exclude_species =
            std::find(species_to_exclude.begin(), species_to_exclude.end(),
                      particle.pdgcode()) != species_to_exclude.end();
        if (exclude_species || !particle.is_hadron() ||
            particle.pdgcode().charmness() != 0) {
            density = 0;
        } else {
            const double mass = particle.mass();
            const double J = particle.spin() * 0.5;
            const double stat = static_cast<int>(round(2. * J)) & 1 ? -1. : 1.;
            const double muf = chemical_potential(particle, surf[iel]);
            for (int i = 1; i < 11; i++)
                density += (2. * J + 1.) * pow(gevtofm, 3) /
                            (2. * pow(TMath::Pi(), 2)) * mass * mass * surf[iel].T *
                            pow(stat, i + 1) *
                            TMath::BesselK(2, i * mass / surf[iel].T) *
                            exp(i * muf / surf[iel].T) / i;
        }
        if (ip > 0)
            cumulantDensity[ip] = cumulantDensity[ip - 1] + density;
        else
            cumulantDensity[ip] = density;
        totalDensity += density;
        ip += 1;
    }
    return cumulantDensity;
}

// Sample momentum with equilibrium distribution (avoids oversampling)
std::tuple<double, double, double> sample_momentum_equilibrium(int iel, double mass, double muf, double stat, bool random_angles = true) {
  double p=0.0, phi=0.0, sinth=0.0, rval=0.0, W=0.0;
  
  fthermal->SetParameters(surf[iel].T, muf, mass, stat);
  do {            // fast momentum generation loop
    p = fthermal->GetRandom();
    if (random_angles) {
    phi = 2.0 * TMath::Pi() * rnd->Rndm();
    sinth = -1.0 + 2.0 * rnd->Rndm();
    }
    mom.SetPxPyPzE(p * sqrt(1.0 - sinth * sinth) * cos(phi),
                    p * sqrt(1.0 - sinth * sinth) * sin(phi), p * sinth,
                    sqrt(p * p + mass * mass));
    W = (surf[iel].dsigma[0] * mom.E() + surf[iel].dsigma[1] * mom.Px() +
          surf[iel].dsigma[2] * mom.Py() +
          surf[iel].dsigma[3] * mom.Pz()) /
        mom.E();
    rval = rnd->Rndm() * dsigmaMax;
    //niter++;
  } while (rval > W);  // end fast momentum generation
  return std::tuple<double, double, double>(p, phi, sinth);
}

std::tuple<double, double, double> sample_momentum(int iel, double mass, double muf, double stat, bool random_angles = true) {
  double p=0.0, phi=0.0, sinth=0.0, W=0.0, WviscFactor=1.0;

  fthermal->SetParameters(surf[iel].T, muf, mass, stat);
  p = fthermal->GetRandom();
  if (random_angles) {
    phi = 2.0 * TMath::Pi() * rnd->Rndm();
    sinth = -1.0 + 2.0 * rnd->Rndm();
  } else {
    // For testing, one might want to sample only in one direction to compare to the distribution function
    // In that case, sample only in x-direction
    phi = 0.0;
    sinth = 0.0;
  }
  mom.SetPxPyPzE(p * sqrt(1.0 - sinth * sinth) * cos(phi),
                  p * sqrt(1.0 - sinth * sinth) * sin(phi), p * sinth,
                  sqrt(p * p + mass * mass));
  W = (surf[iel].dsigma[0] * mom.E() + surf[iel].dsigma[1] * mom.Px() +
        surf[iel].dsigma[2] * mom.Py() +
        surf[iel].dsigma[3] * mom.Pz()) /
      mom.E();

  double momArray[4] = {mom[3], mom[0], mom[1], mom[2]};
  WviscFactor += W_shear_correction(momArray, muf, stat, surf[iel]);
  WviscFactor -= W_bulk_correction(p, mass, muf, stat, surf[iel]);

  if (WviscFactor < 0.1) WviscFactor = 0.1;
  if (WviscFactor > 1.5) WviscFactor = 1.5;  // test, jul17. before: 1.5 // March26: Upper limit needed to avoid large corrections )

  double keep_sigma = W / surf[iel].dsigma[0]; // acceptance probability from the sigma factor
  double keep_viscous = 0.5 * WviscFactor; //
  double acceptance_probability = keep_sigma * keep_viscous;
  double random_number = rnd->Rndm();
  if (random_number > acceptance_probability) {
    return std::tuple<double, double, double>(
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN()
  );  // reject the particle and return nan
  }

  return std::tuple<double, double, double>(p, phi, sinth);
}

bool generate_particle(int iel, int ievent, double dvEff) {
  int isort = 0;
  // SMASH random number [0..1]
  double xsort = rnd->Rndm() * totalDensity;  // throw dice, particle sort
  while (cumulantDensity[isort] < xsort) isort++;
  auto &part = (*database)[isort];
  // By definition, the spin in SMASH is defined as twice the spin of the
  // multiplet, so that it can be stored as an integer. Hence, it needs to
  // be multiplied by 1/2
  const double J = part.spin() * 0.5;
  const double mass = part.mass();
  const double stat = static_cast<int>(round(2. * J)) & 1 ? -1. : 1.;
  // SMASH quantum charges for the hadron state1
  const double muf = chemical_potential(part, surf[iel]);
  if (muf >= mass) {
    std::cout << " ^^ muf = " << muf << "  " << part.pdgcode()
              << std::endl;
  }
  fthermal->SetParameters(surf[iel].T, muf, mass, stat);
  // const double dfMax = part->GetFMax() ;
  double p=0.0, phi=0.0, sinth=0.0;
  if (params::shear_viscosity_enabled || params::bulk_viscosity_enabled) {
    std::tie(p, phi, sinth) = sample_momentum(iel, mass, muf, stat, true);
    if (std::isnan(p) || std::isnan(phi) || std::isnan(sinth)) {
      return false;  // reject the particle
    }
  }
  else {
    std::tie(p, phi, sinth) = sample_momentum_equilibrium(iel, mass, muf, stat, true);
  }
  // Position and boost
  const double x = surf[iel].four_position[1];
  const double y = surf[iel].four_position[2];
  double t = 0, z = 0, vx = 0, vy = 0, vz = 0;
  
  params::deta_dz = std::cbrt(dvEff);
  double smearing_eta_z = params::deta_dz * (-0.5 + rnd->Rndm());
  
  if (params::hydro_coordinate_system == "tau-eta") {
    smearing_eta_z /=
        (surf[iel].four_position[0] * cosh(surf[iel].four_position[3]));
    const double etaF = 0.5 * log((surf[iel].u[0] + surf[iel].u[3]) /
                                  (surf[iel].u[0] - surf[iel].u[3]));
    vx = surf[iel].u[1] / surf[iel].u[0] * cosh(etaF) /
         cosh(etaF + smearing_eta_z);
    vy = surf[iel].u[2] / surf[iel].u[0] * cosh(etaF) /
         cosh(etaF + smearing_eta_z);
    vz = tanh(etaF + smearing_eta_z);
    t = surf[iel].four_position[0] *
        cosh(surf[iel].four_position[3] + smearing_eta_z);
    z = surf[iel].four_position[0] *
        sinh(surf[iel].four_position[3] + smearing_eta_z);
  } else if (params::hydro_coordinate_system == "cartesian") {
    vx = surf[iel].u[1] / surf[iel].u[0];
    vy = surf[iel].u[2] / surf[iel].u[0];
    vz = surf[iel].u[3] / surf[iel].u[0];
    t = surf[iel].four_position[0];
    z = surf[iel].four_position[3] + smearing_eta_z;
  }
  
  mom.Boost(vx, vy, vz);
  smash::FourVector momentum(mom.E(), mom.Px(), mom.Py(), mom.Pz());
  smash::FourVector position(t, x, y, z);
  smash::ParticleData *particle_ptr =
      acceptParticle(ievent, &part, position, momentum);

  if (params::spin_vector_enabled) {
    spin::calculate_and_set_spin_vector(ievent, surf[iel], particle_ptr);
  }
  return true;  // accept the particle and return
}

void generate() {
  ROOT::EnableThreadSafety();
  for (int iev = 0; iev < params::number_of_events; iev++) npart[iev] = 0;
  int ntherm_fail = 0;                                             
  for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
    // ---> thermal densities, for each surface element
    totalDensity = 0.0;
    if (surf[iel].T <= 0.) {
      ntherm_fail++;
      continue;
    }
    // Generate particle densities:
    calculate_particle_densities(iel, database);

    if (totalDensity < 0. || totalDensity > 100.) {
      ntherm_fail++;
      continue;
    }
    // cout<<"thermal densities calculated.\n" ;
    // cout<<cumulantDensity[NPART-1]<<" = "<<totalDensity<<endl ;
    // ---< end thermal densities calc
    double dvEff = 0.0;
    // dvEff = dsigma_mu * u^mu
    dvEff = surf[iel].dsigma[0];
    for (int ievent = 0; ievent < params::number_of_events; ievent++) {
      // ---- number of particles to generate
      int nToGen = 0;
      if (dvEff * totalDensity < 0.01) {
        // SMASH random number [0..1]
        double x = rnd->Rndm();  // throw dice
        if (x < dvEff * totalDensity) nToGen = 1;
      } else {
        // SMASH random number according to Poisson DF
        nToGen = rnd->Poisson(dvEff * totalDensity);
        // Corrections to the number of particles due to viscous effects
        // requires a rejection scheme where more particles are generated
        // than in the equilibrium case, and then rejected with a certain probability 
        // According to "10.1016/j.cpc.2020.107604"
        if (params::bulk_viscosity_enabled || params::shear_viscosity_enabled) {
          nToGen = rnd->Poisson(2.0 * dvEff * totalDensity);
        }
      }
      // ---- we generate a particle!
      for (int ipart = 0; ipart < nToGen; ipart++) {
        generate_particle(iel, ievent, dvEff); 
      }  // end particle generation loop
    }  // events loop
    
    if (Nelem > 50) { //otherwise throws error for small number of elements
      if (iel % (Nelem / 50) == 0) {
        int progress_in_percent = round(iel / (Nelem * 0.01));
        std::printf("[%3i%%] done\t(maxiter: %10i)\n", progress_in_percent,
                    nmaxiter);
        std::fflush(stdout);
      }
    }
  }  // loop over all elements
  std::cout << "\nThermodynamically failed elements: " << ntherm_fail
            << "\n(caused by negative temperatures or if the sum\n"
               "of thermal densities is below 0 or above 100)\n\n";
  delete fthermal;
}

smash::ParticleData *acceptParticle(int ievent,
                                    const smash::ParticleTypePtr &ldef,
                                    smash::FourVector position,
                                    smash::FourVector momentum) {
  int &npart1 = npart[ievent];

  smash::ParticleData *new_particle = new smash::ParticleData(*ldef);
  new_particle->set_4momentum(momentum);
  new_particle->set_4position(position);

  pList[ievent][npart1] = new_particle;
  npart1++;
  if (std::isinf(momentum.x0()) || std::isnan(momentum.x0())) {
    std::cout << "acceptPart nan: known, coord=" << position << std::endl;
    std::exit(1);
  }
  if (npart1 > NPartBuf) {
    std::cerr << "ERROR: Please increase gen::NPartBuf\n";
    std::exit(1);
  }
  return new_particle;
}

// ################### end #################
}  // end namespace gen
