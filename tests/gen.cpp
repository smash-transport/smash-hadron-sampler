#include "gen.h"
#include "const.h"
#include "params.h"
#include "iostream"
#include <ctime>
#include <limits>
#include <mutex>
#include <numeric>
#include <random>
#include <vector>

#include <TF1.h>
#include <TH1D.h>
#include <TFile.h>
#include <TRandom3.h>

#include <TLorentzVector.h>
#include "virtest/vir/test.h"

using namespace gen;

// Custom function to compare two doubles within a given tolerance
static bool expect_near(double val1, double val2, double abs_error) {
  return std::abs(val1 - val2) <= abs_error;
}

// Ensures the global ParticleType list is initialized only once using
// std::call_once. Avoids the "Type list was already built!" exception when
// multiple tests need it.
static void ensure_particletype_initialized() {
  static std::once_flag flag;
  std::call_once(flag, [] {
    smash::ParticleType::create_type_list(
        "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
        "π⁰ 0.1380 0      - 111\n"
        "ρ  0.776  0.149  - 113  213\n"
        "N⁺ 0.938  0      + 2212\n"
        "Δ  1.232  0.117  + 2224 2214 2114 1114\n"
        "Λ  1.116 2.5e-15 + 3122  # stable by default");
  });
}

static void initialize_test_hypersurface() {
    // Initialize surface element 
  gen::surf = new element[1];
  gen::surf[0].four_position[0] = 1.0;
  gen::surf[0].four_position[1] = 0.0;
  gen::surf[0].four_position[2] = 0.0;
  gen::surf[0].four_position[3] = 0.0;
  gen::surf[0].u[0] = 1.0;
  gen::surf[0].u[1] = .0;
  gen::surf[0].u[2] = .0;
  gen::surf[0].u[3] = .0;
  gen::surf[0].dsigma[0] = 10.0; // Effective spatial volume in rest frame
  gen::surf[0].dsigma[1] = 0.0;
  gen::surf[0].dsigma[2] = 0.0;
  gen::surf[0].dsigma[3] = 0.0;
  gen::surf[0].T = 0.155;
  gen::surf[0].mub = 0.0;
  gen::surf[0].muq = 0.0;
  gen::surf[0].mus = 0.0;
  for (int i = 0; i < 10; i++) gen::surf[0].pi[i] = 0.0;
  gen::surf[0].Pi = 0.0;

  gen::dsigmaMax = gen::surf[0].dsigma[0] + sqrt(gen::surf[0].dsigma[1] * gen::surf[0].dsigma[1] +
                                                 gen::surf[0].dsigma[2] * gen::surf[0].dsigma[2] + 
                                                 gen::surf[0].dsigma[3] * gen::surf[0].dsigma[3]);
}



TEST(index44) {
  // Call with all valid values and expect passing
  int indices[4][4] = {{0, 1, 3, 6}, {1, 2, 4, 7}, {3, 4, 5, 8}, {6, 7, 8, 9}};
  for (int row = 0; row < 4; row++) {
    for (int column = 0; column < 4; column++) {
      VERIFY(index44(row, column) == indices[row][column]);
    }
  }
}

TEST(index44_index_out_of_range) {
  int valid_index[4] = {0, 1, 2, 3};
  int invalid_index[4] = {-2, -1, 4, 5};
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      // first index out of range
      try {
        index44(invalid_index[i], valid_index[j]);
        std::cout << "index44 unexpectedly passed with first argument "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
      // second index out of range
      try {
        index44(valid_index[i], invalid_index[j]);
        std::cout << "index44 unexpectedly passed with second argument "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
      // both indices out of range
      try {
        index44(invalid_index[i], invalid_index[j]);
        std::cout << "index44 unexpectedly passed with both arguments "
                  << "out of range" << std::endl;
      } catch (std::out_of_range &e) {
        // Exception was caught as expected
      }
    }
  }
}

TEST(fillBoostMatrix) {
  const double v[3] = {0.3, 0.2, -0.1};
  const double gamma =
      1.0 / sqrt(1.0 - (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]));
  const double u[4] = {gamma, gamma * v[0], gamma * v[1], gamma * v[2]};

  double boostMatrix[4][4];
  fillBoostMatrix(-u[1] / gamma, -u[2] / gamma, -u[3] / gamma, boostMatrix);

  // Perform boost on u to get back to (1,0,0,0)
  double u_prime[4];
  for (int mu = 0; mu < 4; ++mu) {
    u_prime[mu] = 0.0;
    for (int nu = 0; nu < 4; ++nu) {
      u_prime[mu] += boostMatrix[mu][nu] * u[nu];
    }
  }

  // Check that u' is approximately (1,0,0,0)
  const double tolerance = 1e-9;
  VERIFY(std::abs(u_prime[0] - 1.0) < tolerance);
  VERIFY(std::abs(u_prime[1]) < tolerance);
  VERIFY(std::abs(u_prime[2]) < tolerance);
  VERIFY(std::abs(u_prime[3]) < tolerance);
}

TEST(chemical_potential) {
  // Initialize the particle type list with a pion (spin 0),
  // proton (spin 1/2), and delta plus (spin 3/2)
  ensure_particletype_initialized();

  // Define the PDG code for the particles
  const smash::PdgCode pdg_pion = 0x111;
  const smash::PdgCode pdg_proton = 0x2212;
  const smash::PdgCode pdg_rho = 0x213;
  const smash::PdgCode pdg_lambda = 0x3122;

  // Create the particle data for the pion, proton, and delta plus
  smash::ParticleData pion =
      smash::ParticleData(smash::ParticleType::find(pdg_pion));
  smash::ParticleData proton =
      smash::ParticleData(smash::ParticleType::find(pdg_proton));
  smash::ParticleData rho =
      smash::ParticleData(smash::ParticleType::find(pdg_rho));
  smash::ParticleData lambda =
      smash::ParticleData(smash::ParticleType::find(pdg_lambda));

  // Set all values of the surface element to 1.0 which
  // do not affect the spin sampling
  gen::element surf_element;
  surf_element.mub = 0.9;
  surf_element.muq = 3.1;
  surf_element.mus = 7.5;

  // Calculate expected chemical potentials
  double expected_chemical_potential_pion = 0.0;
  double expected_chemical_potential_proton = 1.0 * 0.9 + 1.0 * 3.1;
  double expected_chemical_potential_rho = 1.0 * 3.1;
  double expected_chemical_potential_lambda = 1.0 * 0.9 - 1.0 * 7.5;

  // Chemical potentials from function
  double chemical_potential_pion = chemical_potential(&pion, surf_element);
  double chemical_potential_proton = chemical_potential(&proton, surf_element);
  double chemical_potential_rho = chemical_potential(&rho, surf_element);
  double chemical_potential_lambda = chemical_potential(&lambda, surf_element);

  // Perform checks
  VERIFY(expect_near(chemical_potential_pion, expected_chemical_potential_pion,
                     1e-12));
  VERIFY(expect_near(chemical_potential_proton,
                     expected_chemical_potential_proton, 1e-12));
  VERIFY(expect_near(chemical_potential_rho, expected_chemical_potential_rho,
                     1e-12));
  VERIFY(expect_near(chemical_potential_lambda,
                     expected_chemical_potential_lambda, 1e-12));
}

TEST(chemical_potential_from_type_overload) {
  // Ensure particle types are available
  ensure_particletype_initialized();

  // PDG codes
  const smash::PdgCode pdg_pion = 0x111;
  const smash::PdgCode pdg_proton = 0x2212;
  const smash::PdgCode pdg_rho = 0x213;
  const smash::PdgCode pdg_lambda = 0x3122;

  // Build ParticleData so we can easily access type() (const ParticleType&)
  smash::ParticleData pion =
      smash::ParticleData(smash::ParticleType::find(pdg_pion));
  smash::ParticleData proton =
      smash::ParticleData(smash::ParticleType::find(pdg_proton));
  smash::ParticleData rho =
      smash::ParticleData(smash::ParticleType::find(pdg_rho));
  smash::ParticleData lambda =
      smash::ParticleData(smash::ParticleType::find(pdg_lambda));

  // Surface element with chemical potentials
  gen::element surf_element;
  surf_element.mub = 0.9;
  surf_element.muq = 3.1;
  surf_element.mus = 7.5;

  // Expected values (B, Q, S factors)
  const double exp_pion = 0.0;          // B=0, Q=0, S=0
  const double exp_proton = 0.9 + 3.1;  // B=1, Q=1, S=0
  const double exp_rho = 3.1;           // B=0, Q=1, S=0
  const double exp_lambda = 0.9 - 7.5;  // B=1, Q=0, S=-1

  // Call the ParticleType& overload via type()
  const double mu_pion_type =
      gen::chemical_potential(pion.type(), surf_element);
  const double mu_proton_type =
      gen::chemical_potential(proton.type(), surf_element);
  const double mu_rho_type = gen::chemical_potential(rho.type(), surf_element);
  const double mu_lambda_type =
      gen::chemical_potential(lambda.type(), surf_element);

  // Independent checks against expected numbers
  VERIFY(expect_near(mu_pion_type, exp_pion, 1e-12));
  VERIFY(expect_near(mu_proton_type, exp_proton, 1e-12));
  VERIFY(expect_near(mu_rho_type, exp_rho, 1e-12));
  VERIFY(expect_near(mu_lambda_type, exp_lambda, 1e-12));

  // Cross-check: type overload matches pointer overload (one is enough, do
  // proton)
  const double mu_proton_ptr = gen::chemical_potential(&proton, surf_element);
  VERIFY(expect_near(mu_proton_ptr, mu_proton_type, 1e-12));
}

// Gaussian quadrature integration over momentum space
// Evaluates: ∫ p·d³σ (f_eq + δf) · p² dp
// Using Gauss-Laguerre quadrature for the momentum integral
double integrate_distribution_gaussian_quadrature(const element &surf_elem,
                                             const smash::ParticleType &particle) {
  // Gauss-Legendre quadrature points and weights for integration over angles
  // Using 3-point rule for cos(theta), which covers theta integration
  const int n_theta = 3;
  const double theta_weights[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};  // Standard 3-point GL weights, sum=2
  const double cos_theta_pts[3] = {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};  // Standard 3-point GL points
  
  // 10-point rule for better accuracy
  const int n_p = 10;

  // Gauss–Laguerre abscissae and weights for ∫₀^∞ e^{-x} f(x) dx
  const double lag_x[10] = {
    0.137793470540, 0.729454549503, 1.808342901740, 3.401433697855,
    5.552496140064, 8.330152746764, 11.843785837900, 16.279257831378,
    21.996585811981, 29.920697012274
  };

  const double lag_w[10] = {
    0.308441115765, 0.401119929155, 0.218068287612, 0.0620874560987,
    0.00950151697518, 0.000753008388588, 0.0000282592334960,
    0.000000424931398496, 0.00000000183956482398,
    0.00000000000991182721961
  };
  
  const double mass = particle.mass();
  const double T = surf_elem.T;
  const double muf = chemical_potential(particle, surf_elem);
  const double J = particle.spin() * 0.5;
  const double stat = static_cast<int>(round(2. * J)) & 1 ? -1. : 1.;
  //const double gmumu[4] = {1., -1., -1., -1.};
  
  double integral = 0.0;
  double p_dot_dsigma = 0.0;
  
  // Transform normalized Legendre points to [0, p_max]
  for (int ip = 0; ip < n_p; ip++) {
    const double x  = lag_x[ip];     // x = p/T
    const double wx = lag_w[ip];

    const double p  = T * x;         // p = T x
    const double dp_factor = T;      // dp = T dx

    // Angular integration: loop over theta and phi points
    for (int ith = 0; ith < n_theta; ith++) {
      double cos_theta = cos_theta_pts[ith];
      double sin_theta = sqrt(1.0 - cos_theta * cos_theta);
      double theta_weight = theta_weights[ith];
      
      // For phi integration: use a simple Gauss-Legendre 2-point rule for demo
      // or average over several points
      const int n_phi = 4;  // Use 4 points for phi integration
      double phi_weight = 2.0 * TMath::Pi() / n_phi;
      
      for (int iphi = 0; iphi < n_phi; iphi++) {
        double phi = 2.0 * TMath::Pi() * iphi / n_phi;
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);
        
        double E = sqrt(p * p + mass * mass);
        
        // Full 4-vector dot product: p·d³σ = p^μ d³σ_μ
        // With metric signature (+, -, -, -) and p^μ = (E, px, py, pz)
        // p·d³σ = E*dsigma[0] - px*dsigma[1] - py*dsigma[2] - pz*dsigma[3]
        double px = p * sin_theta * cos_phi;
        double py = p * sin_theta * sin_phi;
        double pz = p * cos_theta;
        p_dot_dsigma = E * surf_elem.dsigma[0]
                            + px * surf_elem.dsigma[1]
                            + py * surf_elem.dsigma[2]
                            + pz * surf_elem.dsigma[3];
        
        if (p_dot_dsigma <= 0.0) continue;  // Only positive contributions (particles exiting)

        // Fermi-Dirac or Bose-Einstein distribution
        double feq =
        C_Feq / (exp((sqrt(p * p + mass * mass) - muf) / surf_elem.T) - stat);
        // Viscous corrections factor
        double WviscFactor = 1.0;
        // Shear viscosity correction:
        if (params::shear_viscosity_enabled) {
          double momArray[4] = {E, px, py, pz};
          WviscFactor += gen::W_shear_correction(momArray, muf, stat, surf_elem);
        }
        // Bulk viscosity correction:
        if (params::bulk_viscosity_enabled) {
          WviscFactor -= gen::W_bulk_correction(p, mass, muf, stat, surf_elem);
        }

        if (WviscFactor < 0.1) WviscFactor = 0.1;
        if (WviscFactor > 1.5) WviscFactor = 1.5;  // test, jul17. before: 1.5

        double integrand = (p_dot_dsigma / E) * feq * WviscFactor * p * p;
        
        integral += integrand
                    * std::exp(x)
                    * dp_factor
                    * wx
                    * theta_weight
                    * phi_weight;
      }
    }
  }

  integral *= (2. * J + 1.0); // Scale by (2J+1) for spin degeneracy
  return integral;
}

// This checks if the number of hadrons from grand canonical ensemble matches the number calculated from integrating the distribution function over the surface element
TEST(number_of_hadrons_to_generate) {
  // Initialize globals
  gen::rnd = new TRandom3(0);
  
  ensure_particletype_initialized();
  gen::database = &smash::ParticleType::list_all();
  
  initialize_test_hypersurface();
  
  // Calculate particle densities
  gen::totalDensity = 0.0;

  // Cumulative densities for each particle type
  double* cumulantDens = gen::calculate_particle_densities(0, gen::database);
  
  double dvEff = gen::surf[0].dsigma[0];
  
  double expected_multiplicity = 0.0;
  double relative_error = 0.0;
  int idx_part = 0;
  for (const auto& particle : smash::ParticleType::list_all()) {
    expected_multiplicity += integrate_distribution_gaussian_quadrature(gen::surf[0], particle);
    relative_error = std::abs(expected_multiplicity - dvEff * cumulantDens[idx_part]) / expected_multiplicity;
    VERIFY(expect_near(relative_error, 0.0, 1e-3));
    idx_part++;
  }
}

TEST(number_of_hadrons_without_visc_corrections) {
  // Initialize globals
  gen::rnd = new TRandom3(0);
  
  ensure_particletype_initialized();
  gen::database = &smash::ParticleType::list_all();
  
  // Initialize surface element 
  initialize_test_hypersurface();
  
  // Initialize fthermal
  gen::fthermal = new TF1("fthermal", gen::ffthermal, 0.0, 10.0, 4);
  
  // Initialize particle storage arrays
  int numEvents = 10000;
  gen::npart = new int[numEvents];
  gen::pList = new smash::ParticleData**[numEvents];
  for (int i = 0; i < numEvents; i++) {
    gen::npart[i] = 0;
    gen::pList[i] = new smash::ParticleData*[gen::NPartBuf];
  }
  
  // Calculate particle densities (this sets global cumulantDensity and totalDensity)
  gen::totalDensity = 0.0;                                     
  double* cumulantDens = gen::calculate_particle_densities(0, gen::database);
  double dvEff = gen::surf[0].dsigma[0];
  
  params::shear_viscosity_enabled = false;
  params::bulk_viscosity_enabled = false;
  params::ecrit = 0.5;  // GeV/fm³

  double mean_nToGen = dvEff * gen::totalDensity;

  double nparticles = 0;
  double niter = 0;
  for (int ievent = 0; ievent < numEvents; ievent++) {
    // Simple sampling for number of particles to generate: use the mean directly for testing (less noisy than Poisson sampling)
    int nToGen = static_cast<int>(mean_nToGen);  // floor
    double fractional = mean_nToGen - nToGen;
    if (gen::rnd->Rndm() < fractional) nToGen++;  // add 1 with probability = fractional part
    for (int ipart = 0; ipart < nToGen; ipart++) {
      bool accepted = generate_particle(0, ievent, dvEff);
      niter++;
      if (accepted) {
        nparticles++;
      }
    }
  }
  double expected_multiplicity = 0.0;
  double relative_error = 0.0;
  int idx_part = 0;
  for (const auto& particle : smash::ParticleType::list_all()) {
    expected_multiplicity += integrate_distribution_gaussian_quadrature(gen::surf[0], particle);
    relative_error = std::abs(expected_multiplicity - dvEff * cumulantDens[idx_part]) / expected_multiplicity;
    VERIFY(expect_near(relative_error, 0.0, 1e-3));
    idx_part++;
  }
}

TEST(number_of_hadrons_with_visc_corrections) {
  // Initialize globals
  gen::rnd = new TRandom3(0);
  
  ensure_particletype_initialized();
  gen::database = &smash::ParticleType::list_all();
  
  initialize_test_hypersurface();  

  gen::surf[0].pi[gen::index44(0, 3)] = -0.05;
  gen::surf[0].pi[gen::index44(0, 0)] = 0.025;
  gen::surf[0].pi[gen::index44(1, 1)] = 0.01;
  gen::surf[0].pi[gen::index44(2, 2)] = 0.01; 
  // Bulk viscous pressure
  gen::surf[0].Pi = -0.2; 
  // Initialize fthermal
  gen::fthermal = new TF1("fthermal", gen::ffthermal, 0.0, 10.0, 4);
  
  // Initialize particle storage arrays
  int numEvents = 10000;
  gen::npart = new int[numEvents];
  gen::pList = new smash::ParticleData**[numEvents];
  for (int i = 0; i < numEvents; i++) {
    gen::npart[i] = 0;
    gen::pList[i] = new smash::ParticleData*[gen::NPartBuf];
  }

  // Calculate particle densities (this sets global cumulantDensity and totalDensity)
  gen::totalDensity = 0.0;                                               
  double* cumulantDens = gen::calculate_particle_densities(0, gen::database);

  double dvEff = gen::surf[0].dsigma[0];
  
  params::shear_viscosity_enabled = true;
  params::bulk_viscosity_enabled = true;
  params::ecrit = 0.5;  // GeV/fm³
  
  double mean_nToGen = 2.0 * dvEff * gen::totalDensity;
  std::cout << ", mean_to_gen: " << mean_nToGen << "\n";
  //double nToGen = rnd->Poisson(2.0 * dvEff * gen::totalDensity);
  double nparticles = 0;
  double niter = 0;
  for (int ievent = 0; ievent < numEvents; ievent++) {
    // Simple sampling for number of particles to generate: use the mean directly for testing
    int nToGen = static_cast<int>(mean_nToGen);  // floor
    double fractional = mean_nToGen - nToGen;
    if (gen::rnd->Rndm() < fractional) nToGen++;  // add 1 with probability = fractional part
    for (int ipart = 0; ipart < nToGen; ipart++) {
      bool accepted = generate_particle(0, ievent, dvEff);
      niter++;
      if (accepted) {
        nparticles++;
      }
    }
  }
  double expected_total_multiplicity = 0.0;
  for (const auto& particle : smash::ParticleType::list_all()) {
    double delta_N = integrate_distribution_gaussian_quadrature(gen::surf[0], particle);
    expected_total_multiplicity += delta_N;
  }

  double relative_error = std::abs(expected_total_multiplicity - nparticles / numEvents) / expected_total_multiplicity;
  std::cout << "Expected total multiplicity: " << expected_total_multiplicity << "\n";
  std::cout << "Average generated multiplicity: " << nparticles / numEvents << "\n";
  std::cout << "Relative error: " << relative_error << "\n";
  std::cout << "Acceptance rate: " << nparticles / niter << "\n";
  VERIFY(expect_near(relative_error, 0.0, 2e-2));
  delete[] gen::surf;
}

// Test wether the momentum distribution of generated particles matches the expected distribution
TEST(momentum_generation) {
  // Initialize globals
  gen::rnd = new TRandom3(0);
  
  ensure_particletype_initialized();
  gen::database = &smash::ParticleType::list_all();
  
  initialize_test_hypersurface();

  // Set shear viscosity
  gen::surf[0].pi[gen::index44(0, 3)] = -0.1;
  gen::surf[0].pi[gen::index44(0, 0)] = 0.025;
  gen::surf[0].pi[gen::index44(1, 1)] = 0.01;
  gen::surf[0].pi[gen::index44(2, 2)] = 0.01; 
  // Bulk viscous pressure
  gen::surf[0].Pi = 0.01;
  
  // Initialize fthermal
  gen::fthermal = new TF1("fthermal", gen::ffthermal, 0.0, 10.0, 4);
  
  gen::totalDensity = 0.0;                       
  gen::calculate_particle_densities(0, gen::database);
  double dvEff = gen::surf[0].dsigma[0];
  
  params::shear_viscosity_enabled = true;
  params::bulk_viscosity_enabled = true;
  params::ecrit = 0.5;

  // ============= Pion Momentum Distribution Test =============
  int numSamples = 100; 
  double mass = 0.135;   // pion mass
  double muf = 0.0;      // chemical potential
  double stat = 1.0;     // bosons
  int iel = 0;
  double T = gen::surf[iel].T;

  gen::fthermal->SetParameters(T, muf, mass, stat);
  
  TH1D* hNewVisc = new TH1D("hNewVisc", "New Visc Sampling;p [GeV];dN/dp", 50, 0, 3.0);
  
  // Sample momenta
  int acceptedNew = 0;
  for (int i = 0; i < numSamples; i++) {
    auto [p, phi, sinth] = gen::sample_momentum(iel, mass, muf, stat, false);
    if (!std::isnan(p)) {
      hNewVisc->Fill(p);
      acceptedNew++;
    }
  }

  // Theory with viscous corrections: f_eq * p² * (1 + δf_shear - δf_bulk)
  auto theoryFunc = [&](double* x, double* ) -> double {
    double p = x[0];
    double E = sqrt(p*p + mass*mass);
    double f_eq = p*p / (exp((E - muf) / T) - stat);
    
    // Use existing viscous correction functions (px=p, py=pz=0)
    double momArray[4] = {E, p, 0, 0};
    double viscFactor = 1.0;
    viscFactor += gen::W_shear_correction(momArray, muf, stat, gen::surf[iel]);
    viscFactor -= gen::W_bulk_correction(p, mass, muf, stat, gen::surf[iel]);
    if (viscFactor < 0.1) viscFactor = 0.1;
    if (viscFactor > 1.5) viscFactor = 1.5;
    return f_eq * viscFactor;
  };

  TF1* fTheory = new TF1("fTheory", theoryFunc, 0, 3.0, 0);

  // Normalize histograms
  double binWidth = hNewVisc->GetBinWidth(1);
  hNewVisc->Scale(1.0 / (acceptedNew * binWidth));

  // Normalize theory with viscous corrections
  double theoryIntegral = fTheory->Integral(0, 3.0);
  TH1D* hTheory = new TH1D("hTheory", "Theory (visc);p [GeV];dN/dp", 50, 0, 3.0);
  for (int bin = 1; bin <= hTheory->GetNbinsX(); bin++) {
    double p = hTheory->GetBinCenter(bin);
    hTheory->SetBinContent(bin, fTheory->Eval(p) / theoryIntegral);
  }

  // Chi-squared
  double chi2New = 0;
  int nBinsUsed = 0;
  for (int bin = 1; bin <= hNewVisc->GetNbinsX(); bin++) {
    double obs = hNewVisc->GetBinContent(bin);
    double exp = hTheory->GetBinContent(bin);
    double err = hNewVisc->GetBinError(bin);
    if (err > 0 && exp > 0) {
      chi2New += (obs - exp) * (obs - exp) / (err * err);
      nBinsUsed++;
    }
  }
  double pvalue = TMath::Prob(chi2New, nBinsUsed);
  std::cout << "p-value = " << pvalue << std::endl;
  VERIFY(pvalue > 0.01); 

  delete hNewVisc;
  delete hTheory;
  delete fTheory;
}