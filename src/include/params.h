#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

namespace params {
extern char surface_file[255], output_directory[255];
extern bool bulk_viscosity_enabled, createRootOutput, shear_viscosity_enabled;
extern int NEVENTS;
extern double dx, dy, deta;
extern double ecrit, speed_of_sound_squared, ratio_pressure_energydensity;
// extern double Temp, mu_b, mu_q, mu_s;

// ---- Routines ----
void readParams(char *filename);
void printParameters();
} // namespace params

#endif // INCLUDE_PARAMS_H_
