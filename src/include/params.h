#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

namespace params {
extern char sSurface[255], sSpectraDir[255];
extern bool weakContribution;
extern bool rescatter;
extern bool shear;
extern bool bulk;
// extern double Temp, mu_b, mu_q, mu_s ;
extern int NEVENTS;
extern double NBINS, QMAX;
extern double dx, dy, deta;
extern double ecrit, cs2, ratio_pressure_energydensity;
/**
 * Controls whether spin sampling is used during simulations, in order to set
 * the projections of spin degrees of freedom at freezeout.
 *
 * **Default Value:** False
 */
extern bool is_spin_sampling_on;
/**
 * This parameter influences the probability of sampling a spin projection that
 * aligns with the vorticity by more than 50%.
 *
 * A value of 0.0 indicates no preference (equal probability for any spin
 * projection), while larger values favor spins aligned with the (local)
 * vorticity and negative values anti-alignment.
 * 
 * **Default Value:** False
 */
extern double global_polarization;

// ---- rooutines ----
void readParams(char* filename);
void printParameters();
}  // namespace params

#endif  // INCLUDE_PARAMS_H_
