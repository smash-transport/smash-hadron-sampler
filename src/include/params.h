#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

namespace params{
extern char sSurface [255], sSpectraDir [255];
extern bool weakContribution ;
extern bool rescatter ;
extern bool shear ;
extern bool bulk ;
//extern double Temp, mu_b, mu_q, mu_s ;
extern int NEVENTS ;
extern double NBINS, QMAX ;
extern double dx, dy, deta ;
extern double ecrit, cs2, ratio_pressure_energydensity ;

// ---- rooutines ----
void readParams(char* filename) ;
void printParameters() ;
}

#endif // INCLUDE_PARAMS_H_
