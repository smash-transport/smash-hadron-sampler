#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

namespace params {
extern char sSurface [255], sSpectraDir [255];
extern bool bulk, createRootOutput, shear;
extern int NEVENTS;
extern double dx, dy, deta;
extern double ecrit, cs2, ratio_pressure_energydensity;
//extern double Temp, mu_b, mu_q, mu_s;

// ---- Routines ----
void readParams(char* filename);
void printParameters();
}

#endif // INCLUDE_PARAMS_H_
