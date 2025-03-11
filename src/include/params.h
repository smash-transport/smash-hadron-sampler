#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

namespace params {
/**
 * sVorticity is the path to the file containing the 16 components of the
 * thermal vorticity tensor which is given as dbeta.dat by vHLLE.
 */
extern char sSurface[255], sSpectraDir[255], sVorticity[255];
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
extern bool is_spin_sampling_on, vorticity_output_enabled;

/**
 * Helper function to get the directory part of a file path.
 */
std::string getDirectory(const std::string& filePath);

// ---- rooutines ----
void readParams(char* filename);
void printParameters();
}  // namespace params

#endif  // INCLUDE_PARAMS_H_
