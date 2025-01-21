#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

#include <string>

namespace params {
extern std::string surface_file, output_directory;
extern bool bulk_viscosity_enabled, createRootOutput, shear_viscosity_enabled;
extern int NEVENTS;
extern double dx, dy, deta;
extern double ecrit, speed_of_sound_squared, ratio_pressure_energydensity;
// extern double Temp, mu_b, mu_q, mu_s;

// ---- Routines ----
void readParams(const std::string &filename);
void printParameters();
} // namespace params

#endif // INCLUDE_PARAMS_H_
