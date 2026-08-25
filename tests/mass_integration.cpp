#include "gen.h"
#include "const.h"
#include "params.h"
#include "mass_integration.h"

#include <cmath>
#include <iostream>
#include <random>
#include <string>

#include <TF1.h>

#include "smash/decaymodes.h"

#include "virtest/vir/test.h"

// maximum relative error
static double tolerance = 0.01;

static void initialize() {
  smash::ParticleType::create_type_list(
      "# NAME MASS[GEV] WIDTH[GEV] PARITY PDG\n"
        "π 0.1380 0      - 111 211\n"
        "ρ  0.776 0.149  - 113  213\n"
        "f₂ 1.276 0.187 + 225\n"
      );

  smash::DecayModes::load_decaymodes(
      "ρ\n"
      "1. 1 π π\n"
      "f₂\n"
      "0.5 2 π π\n"
      "0.5 2 ρ ρ\n"
    );

}

TEST(initialize) {
  initialize();
}

TEST_CATCH(cannot_call_integrand_with_pole_mass, std::invalid_argument) {
  params::spectral_function_method = SpectralFunction::PoleMass;
  const auto& simple = smash::ParticleType::find(0x113);
  MassIntegration::details::integrand(1, 0.7, 0.15, simple);
}

TEST_CATCH(cannot_call_integrand_with_stable, std::invalid_argument) {
  params::spectral_function_method = SpectralFunction::BreitWigner;
  const auto& stable = smash::ParticleType::find(0x111);
  MassIntegration::details::integrand(1, stable.mass(), 0.15, stable);
}

TEST(spectral_function_normalization) {
  smash::Integrator integrate;
  for (const auto method : {SpectralFunction::BreitWigner, SpectralFunction::FullVacuum}) {
    params::spectral_function_method = method;
    const auto& resonance = smash::ParticleType::find(0x225);
    const double integral = integrate(resonance.min_mass_kinematic(), 100, [&](double m) {
      return MassIntegration::details::spectral_function(m, resonance);
    });
    COMPARE_RELATIVE_ERROR(integral, 1., tolerance);
  }
}

TEST(interpolation_works_stable) { 
  const auto& stable = smash::ParticleType::find(0x111);  
  const double m = stable.mass();
  auto random_value = smash::random::make_uniform_distribution(0.1, 0.16); // not a tabulation point
  for (auto method: {SpectralFunction::PoleMass, SpectralFunction::BreitWigner, SpectralFunction::FullVacuum}) {
    for (size_t i = 1; i < 11; ++i) {
      params::spectral_function_method = method;
      const int sum_index = i;
      const auto table = MassIntegration::tabulate_over_temperature(sum_index, stable);
      double T_test = random_value();
      double value_at_T_test = table.get_value_linear(T_test);
      double expected_value =  m * m * std::cyl_bessel_k(2, sum_index*m/T_test);
      COMPARE_RELATIVE_ERROR(value_at_T_test, expected_value, tolerance);
    }
  }
}

TEST(writeout_densities) {
  gen::surf = new gen::element[1];
  gen::database = &smash::ParticleType::list_all();

  gen::surf[0].four_position[0] = 1.0;
  gen::surf[0].u[0] = 1.0;
  gen::surf[0].dsigma[0] = 1.0;
  gen::surf[0].T = 0.14;
  gen::surf[0].mub = 0.1;
  // Remaining globals are 0 by default
  
  gen::database = &smash::ParticleType::list_all();
  for (const auto method : {SpectralFunction::PoleMass, SpectralFunction::BreitWigner, SpectralFunction::FullVacuum}) {
    gen::totalDensity = 0.0;
    params::spectral_function_method = method;
    std::vector<double> cumulantDens = gen::calculate_partial_densities(0, *gen::database);
    std::cout << "Cumulant densities for method " << static_cast<int>(method) << ":\n";
    for (size_t i = 0; i < cumulantDens.size(); ++i) {
      double density = cumulantDens[i];
      if (density > 0) {
         density -= cumulantDens[i-1];
      }
      std::cout << (*gen::database)[i].name() << ": " << density << "\n";
    }
    MassIntegration::details::reset_tabulated_integrals();
  }

  delete[] gen::surf;
}
