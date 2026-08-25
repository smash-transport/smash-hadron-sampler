#include "mass_integration.h"

#include <cmath>
#include <iostream>
#include <utility>
#include <stdexcept>
#include "smash/integrate.h"

namespace MassIntegration {
namespace {
constexpr int T_npoints = 100;
constexpr double T_min = 0.01;
constexpr double T_max = 0.16;
constexpr double mass_limit = 20.0;
std::map<std::pair<int, int>, smash::Tabulation> tabulated_integrals;
}  // namespace

namespace details {

double spectral_function(double m, const smash::ParticleType& type) {
    if (type.is_stable()) {
        throw std::invalid_argument("Cannot call for spectral function for stable particle");
    }
    switch (params::spectral_function_method) {
        case SpectralFunction::PoleMass:
            throw std::invalid_argument("Cannot call for spectral function when using PoleMass method");
        case SpectralFunction::BreitWigner: {
            double norm = 0.5 + M_1_PI * std::atan(2*(type.mass() - type.min_mass_kinematic())/type.width_at_pole());
            return type.breit_wigner_spectral_function(m) / norm;
        }
        case SpectralFunction::FullVacuum:
            return type.full_spectral_function(m);
        default:
            throw std::invalid_argument("Invalid spectral function method");
    }
}

double integrand(int sum_index, double m, double T, const smash::ParticleType& type) {
    const double K_n = std::cyl_bessel_k(2, sum_index*m/T);
    const double SF = spectral_function(m, type);
    return m * m * SF * K_n;
}

// Meant for testing purposes.
void reset_tabulated_integrals() {
    tabulated_integrals.clear();
}

}  // namespace MassIntegration::details

smash::Tabulation tabulate_over_temperature(int sum_index,
                                            const smash::ParticleType& type) {
    if (params::spectral_function_method == SpectralFunction::PoleMass || type.is_stable()) {
        // For the pole mass case, the integral reduces to a simple evaluation at the pole mass
        return smash::Tabulation(T_min, T_max, T_npoints, [&](double T) {
            double m = type.mass();
            return m * m * std::cyl_bessel_k(2, sum_index*m/T);
        });
    }
    static smash::Integrator integrate;
    const double m_min = type.min_mass_kinematic();
    return smash::Tabulation(T_min, T_max, T_npoints, [&](double T) {
        return integrate(m_min, mass_limit, [&](double m) {
            return details::integrand(sum_index, m, T, type);
        });
    });
}

double fugacity_expansion_coefficient(int sum_index, const smash::ParticleType& type, double T) {
    if (T<T_min) {
        std::cout << "Temperature " << T << " below tabulation minimum of "
                    << T_min << ", returning 0.\n";
        return 0.0;
    }
    if (T>T_max) {
        std::cout << "Temperature " << T << " above tabulation maximum of "
                    << T_max << ", interpolating linearly.\n";
    }
    const auto key = std::make_pair(sum_index, type.pdgcode().get_decimal());
    auto it = tabulated_integrals.find(key);
    if (it == tabulated_integrals.end()) {
        // It does not exist yet, tabulate and store in the map
        smash::Tabulation tabulation = tabulate_over_temperature(sum_index, type);
        tabulated_integrals[key] = tabulation;
        return tabulation.get_value_linear(T);
    }
    return it->second.get_value_linear(T);
}

double sample_mass(const smash::ParticleType& type) {
    if (type.is_stable()) {
        return type.mass();
    }
    switch (params::spectral_function_method) {
        case SpectralFunction::PoleMass:
            return type.mass();
        case SpectralFunction::BreitWigner:
            return type.sample_breit_wigner_spectral_function();
        case SpectralFunction::FullVacuum:
            return type.sample_full_spectral_function();
        default:
            throw std::invalid_argument("Invalid spectral function method");
    }
}

}  // namespace MassIntegration
