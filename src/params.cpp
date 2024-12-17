#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include "params.h"

using namespace std ;

namespace params{

char sSurface [255], sSpectraDir [255];
char sVorticity [255] = "unset";
bool weakContribution ;
bool rescatter ;
bool shear ;
bool bulk ;
bool is_spin_sampling_on;
//double Temp, mu_b, mu_q, mu_s ;
int NEVENTS ;
double NBINS, QMAX ;
double dx, dy, deta ;
double ecrit ;
double cs2=0.15;
double ratio_pressure_energydensity=0.15;
double global_polarization=0.0;

// Helper function to get the directory of a given file path
string getDirectory(const string& filePath) {
    size_t pos = filePath.find_last_of("/\\");
    if (pos != string::npos) {
        return filePath.substr(0, pos + 1);
    }
	// No directory found, return empty
    return "";
}

// ############ reading and processing the parameters
void readParams(char* filename)
{
	char parName [255], parValue [255] ;
	ifstream fin(filename) ;
	if(!fin.is_open()) { cout << "cannot open parameters file " << filename << endl ; exit(1) ; }

	string dbetaFile = "";
	while(fin.good()){
		string line ;
		getline(fin, line) ;
		istringstream sline (line) ;
		sline >> parName >> parValue ;
		if     (strcmp(parName,"surface")==0) strcpy(sSurface, parValue) ;
		else if(strcmp(parName,"dbeta")==0) {
			strcpy(sVorticity, parValue) ;
			dbetaFile = sVorticity;
		}
		else if(strcmp(parName,"spectra_dir")==0) strcpy(sSpectraDir, parValue) ;
		else if(strcmp(parName,"Nbins")==0) NBINS = atoi(parValue) ;
		else if(strcmp(parName,"q_max")==0) QMAX = atof(parValue) ;
		else if(strcmp(parName,"number_of_events")==0) NEVENTS = atoi(parValue) ;
		else if(strcmp(parName,"rescatter")==0) rescatter = atoi(parValue) ;
		else if(strcmp(parName,"weakContribution")==0) weakContribution = atoi(parValue) ;
		else if(strcmp(parName,"shear")==0) shear = atoi(parValue) ;
		else if(strcmp(parName,"bulk")==0) bulk = atoi(parValue) ;
		else if(strcmp(parName,"ecrit")==0) ecrit = atof(parValue) ;
		else if(strcmp(parName,"cs2")==0) cs2 = atof(parValue) ;
		else if(strcmp(parName,"ratio_pressure_energydensity")==0) ratio_pressure_energydensity = atof(parValue) ;
		else if(strcmp(parName, "sample_spin") == 0) is_spin_sampling_on = atoi(parValue);
		else if(parName[0]=='!') cout << "CCC " << sline.str() << endl ;
		else cout << "UUU " << sline.str() << endl ;
	}

	// If no dbeta file was specified, set the directory to the same as the surface file as a default
    if (dbetaFile.empty()) {
    // Construct the path based on the directory of sSurface
    string dir = getDirectory(sSurface);
    string defaultFile = dir + "beta.dat";

    // Safely copy the constructed string to sVorticity
    strncpy(sVorticity, defaultFile.c_str(), sizeof(sVorticity) - 1);
    sVorticity[sizeof(sVorticity) - 1] = '\0'; // Ensure null-termination
}

 deta=0.05 ; dx=dy=0.0 ; // TODO!
}

void printParameters()
{
  cout << "====== parameters ======\n" ;
  cout << "surface = " << sSurface << endl ;
  if(is_spin_sampling_on) cout << "vorticity = " << sVorticity << endl ;
  cout << "spectraDir = " << sSpectraDir << endl ;
  cout << "numberOfEvents = " << NEVENTS << endl ;
  cout << "isRescatter = " << rescatter << endl ;
  cout << "weakContribution = " << weakContribution << endl ;
  cout << "shear_visc_on = " << shear << endl ;
  cout << "bulk_visc_on = " << bulk << endl ;
  cout << "e_critical = " << ecrit << endl ;
  cout << "Nbins = " << NBINS << endl ;
  cout << "q_max = " << QMAX << endl ;
  cout << "cs2 = " << cs2 << endl ;
  cout << "ratio_pressure_energydensity = " << ratio_pressure_energydensity << endl ;
  cout << "sample_spin = " << is_spin_sampling_on << endl ;
  cout << "======= end parameters =======\n" ;
}

}
