#include "gen.h"

#include <TF1.h>
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

// Convert a Lorentz boost matrix from contravariant form (Λ^μ_ν, acts on upper
// indices) to covariant form (Λ_μ^ν, acts on lower indices) using the (+,-,-,-)
// metric.
//
// Index relation:   (Λ_cov)_μ^ν = g_{μα} (Λ_contra)^α_β g^{βν}
// Matrix relation:  Λ_cov = g * Λ_contra * g,  with g = diag(+1,-1,-1,-1)
FourMatrix create_covariant_boost_matrix(
    const double boost_contravariant[4][4]) {
  // Metric signature
  const double g[4] = {+1.0, -1.0, -1.0, -1.0};
  FourMatrix boost_covariant{};

  for (int mu = 0; mu < 4; ++mu)
    for (int nu = 0; nu < 4; ++nu)
      boost_covariant[mu][nu] = g[mu] * boost_contravariant[mu][nu] * g[nu];

  return boost_covariant;
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
  if (params::spin_sampling_enabled) {
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
    instream.clear(); // does not work with gcc 4.1 otherwise
    instream >> surf[n].four_position[0] >> surf[n].four_position[1] >>
        surf[n].four_position[2] >> surf[n].four_position[3] >>
        surf[n].dsigma[0] >> surf[n].dsigma[1] >> surf[n].dsigma[2] >>
        surf[n].dsigma[3] >> surf[n].u[0] >> surf[n].u[1] >> surf[n].u[2] >>
        surf[n].u[3] >> surf[n].T >> surf[n].mub >> surf[n].muq >> surf[n].mus;
    for (int i = 0; i < 10; i++) instream >> surf[n].pi[i];
    instream >> surf[n].Pi;

    // If spin sampling is enabled, load the energy density from the
    // extended freezeout surface
    if (params::spin_sampling_enabled) {
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
    FourMatrix boostMatrixCov{};
    bool boost_matrix_needed =
        params::shear_viscosity_enabled || params::spin_sampling_enabled;

    if (boost_matrix_needed) {
      fillBoostMatrix(-surf[n].u[1] / surf[n].u[0],
                      -surf[n].u[2] / surf[n].u[0],
                      -surf[n].u[3] / surf[n].u[0], boostMatrix);

      // As boostMatrix is a contravariant boost matrix (acts on upper indices),
      // and vorticity is a tensor with lower indices, we need to convert the
      // boost matrix to a covariant form (acts on lower indices).
      boostMatrixCov = create_covariant_boost_matrix(boostMatrix);
    }

    // _pi^{μν} (upper,upper) uses the contravariant boostMatrix
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

    // boost vorticity tensor to fluid rest frame
    if (params::spin_sampling_enabled) {
      if (!surf[n].vorticity.has_value()) {
        throw std::runtime_error("Error: Vorticity is not set for element " +
                                 std::to_string(n));
      }
      (*surf[n].vorticity)->boost_vorticity_to_fluid_rest_frame(boostMatrixCov);
    }
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

  // Load SMASH hadron list
  smash::initialize_default_particles_and_decaymodes();
  const smash::ParticleTypeList &database = smash::ParticleType::list_all();

  // Dump list of hadronic states to terminal
  // for (auto& HadronState: database) {
  //   std::cout << HadronState << '\n';
  // }

  // NPART = total number of hadron states
  NPART = database.size();
  std::cout << "NPART=" << NPART << std::endl;
  std::cout << "dsigmaMax=" << dsigmaMax << "\n\n";
  cumulantDensity = new double[NPART];
}

void enable_vorticity_storage() {
  // Allocate memory for the vorticity vector for each sampled particle
  // if spin sampling and vorticity output are enabled in the config
  if (params::spin_sampling_enabled && params::vorticity_output_enabled) {
    thetaStorage = std::make_unique<std::vector<std::vector<ThetaStruct>>>(
        params::number_of_events);
  } else if (!params::spin_sampling_enabled &&
             params::vorticity_output_enabled) {
    throw std::runtime_error(
        "Vorticity output is enabled but spin sampling is not. "
        "Enable spin sampling in the config file by adding "
        " the line 'sample_spin 1'.");
  }
}

double ffthermal(double *x, double *par) {
  double &T = par[0];
  double &mu = par[1];
  double &mass = par[2];
  double &stat = par[3];
  return x[0] * x[0] / (exp((sqrt(x[0] * x[0] + mass * mass) - mu) / T) - stat);
}

void generate() {
  ROOT::EnableThreadSafety();
  const double gmumu[4] = {1., -1., -1., -1.};
  TF1 *fthermal = new TF1("fthermal", ffthermal, 0.0, 10.0, 4);
  TLorentzVector mom;
  for (int iev = 0; iev < params::number_of_events; iev++) npart[iev] = 0;
  int nmaxiter = 0;
  int ntherm_fail = 0;

  // List species that should not be sampled: photon, electron, muon, tau
  // Sigma meson needs to be excluded to generate correct multiplicities
  std::vector<smash::PdgCode> species_to_exclude{0x11, -0x11, 0x13, -0x13,
                                                 0x15, -0x15, 0x22, 0x9000221};

  for (int iel = 0; iel < Nelem; iel++) {  // loop over all elements
    // ---> thermal densities, for each surface element
    totalDensity = 0.0;
    if (surf[iel].T <= 0.) {
      ntherm_fail++;
      continue;
    }

    const smash::ParticleTypeList &database = smash::ParticleType::list_all();
    int ip = 0;
    for (auto &particle : database) {
      double density = 0.;
      const bool exclude_species =
          std::find(species_to_exclude.begin(), species_to_exclude.end(),
                    particle.pdgcode()) != species_to_exclude.end();
      if (exclude_species || !particle.is_hadron() ||
          particle.pdgcode().charmness() != 0) {
        density = 0;
      } else {
        const double mass = particle.mass();
        // By definition, the spin in SMASH is defined as twice the spin of the
        // multiplet, so that it can be stored as an integer. Hence, it needs to
        // be multiplied by 1/2
        const double J = particle.spin() * 0.5;
        const double stat = static_cast<int>(round(2. * J)) & 1 ? -1. : 1.;
        // SMASH quantum charges for the hadron state
        const double muf = particle.baryon_number() * surf[iel].mub +
                           particle.strangeness() * surf[iel].mus +
                           particle.charge() * surf[iel].muq;
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

    if (totalDensity < 0. || totalDensity > 100.) {
      ntherm_fail++;
      continue;
    }
    // cout<<"thermal densities calculated.\n" ;
    // cout<<cumulantDensity[NPART-1]<<" = "<<totalDensity<<endl ;
    // ---< end thermal densities calc
    double rval, dvEff = 0., W;
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
      }
      // ---- we generate a particle!
      for (int ipart = 0; ipart < nToGen; ipart++) {
        int isort = 0;
        // SMASH random number [0..1]
        double xsort = rnd->Rndm() * totalDensity;  // throw dice, particle sort
        while (cumulantDensity[isort] < xsort) isort++;
        auto &part = database[isort];
        const double J = part.spin() * 0.5;
        const double mass = part.mass();
        const double stat = static_cast<int>(round(2. * J)) & 1 ? -1. : 1.;
        // SMASH quantum charges for the hadron state
        const double muf = part.baryon_number() * surf[iel].mub +
                           part.strangeness() * surf[iel].mus +
                           part.charge() * surf[iel].muq;
        if (muf >= mass)
          std::cout << " ^^ muf = " << muf << "  " << part.pdgcode()
                    << std::endl;
        fthermal->SetParameters(surf[iel].T, muf, mass, stat);
        // const double dfMax = part->GetFMax() ;
        int niter = 0;  // number of iterations, for debug purposes
        do {            // fast momentum generation loop
          const double p = fthermal->GetRandom();
          const double phi = 2.0 * TMath::Pi() * rnd->Rndm();
          const double sinth = -1.0 + 2.0 * rnd->Rndm();
          mom.SetPxPyPzE(p * sqrt(1.0 - sinth * sinth) * cos(phi),
                         p * sqrt(1.0 - sinth * sinth) * sin(phi), p * sinth,
                         sqrt(p * p + mass * mass));
          W = (surf[iel].dsigma[0] * mom.E() + surf[iel].dsigma[1] * mom.Px() +
               surf[iel].dsigma[2] * mom.Py() +
               surf[iel].dsigma[3] * mom.Pz()) /
              mom.E();
          double WviscFactor = 1.0;
          if (params::shear_viscosity_enabled) {
            const double feq =
                C_Feq /
                (exp((sqrt(p * p + mass * mass) - muf) / surf[iel].T) - stat);
            double pipp = 0;
            double momArray[4] = {mom[3], mom[0], mom[1], mom[2]};
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                pipp += momArray[i] * momArray[j] * gmumu[i] * gmumu[j] *
                        surf[iel].pi[index44(i, j)];
            WviscFactor +=
                (1.0 + stat * feq) * pipp /
                (2. * surf[iel].T * surf[iel].T *
                 (params::ecrit +
                  params::ecrit * params::ratio_pressure_energydensity));
          }
          if (params::bulk_viscosity_enabled) {
            const double feq =
                C_Feq /
                (exp((sqrt(p * p + mass * mass) - muf) / surf[iel].T) - stat);
            WviscFactor -=
                (1.0 + stat * feq) * surf[iel].Pi *
                (mass * mass / (3 * mom.E()) -
                 mom.E() * (1.0 / 3.0 - params::speed_of_sound_squared)) /
                (15 * (1.0 / 3.0 - params::speed_of_sound_squared) *
                 (1.0 / 3.0 - params::speed_of_sound_squared) * surf[iel].T *
                 (params::ecrit +
                  params::ecrit * params::ratio_pressure_energydensity));
          }
          if (WviscFactor < 0.1) WviscFactor = 0.1;
          // test, jul17; before: 0.5
          // if(WviscFactor>1.2) WviscFactor = 1.2 ; //              before: 1.5
          W *= WviscFactor;
          rval = rnd->Rndm() * dsigmaMax;
          niter++;
        } while (rval > W); // end fast momentum generation
        if (niter > nmaxiter)
          nmaxiter = niter;
        const double x = surf[iel].four_position[1];
        const double y = surf[iel].four_position[2];
        double t = 0, z = 0, vx = 0, vy = 0, vz = 0;
        /* The deta_dz is an estimate of spatial extent based on the volume of
         * the respective freezeout hypersurface element which is used as a
         * smearing parameter in eta or z direction (depending on the hydro
         * coordinate system).
         * Note: No smearing in x and y direction implemented at the moment.
         */
        params::deta_dz = std::cbrt(dvEff);
        double smearing_eta_z = params::deta_dz * (-0.5 + rnd->Rndm());
        if (params::hydro_coordinate_system == "tau-eta") {
          smearing_eta_z /=
              (surf[iel].four_position[0] * cosh(surf[iel].four_position[3]));
          // additional random smearing over eta
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
  
        // Calculate and set the spin vector if spin sampling is enabled
        if (params::spin_sampling_enabled) {
          spin::calculate_and_set_spin_vector(ievent, surf[iel], particle_ptr);
        }

      } // coordinate accepted
    }   // events loop
    if (iel % (Nelem / 50) == 0) {
      int progress_in_percent = round(iel / (Nelem * 0.01));
      std::printf("[%3i%%] done\t(maxiter: %10i)\n", progress_in_percent,
                  nmaxiter);
      std::fflush(stdout);
    }
  } // loop over all elements
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
