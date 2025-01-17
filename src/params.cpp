#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "params.h"

using namespace std;

namespace params {

char surface_file[255], output_directory[255];
bool bulk_viscosity_enabled{false}, createRootOutput{false},
    shear_viscosity_enabled{false};
int NEVENTS;
double dx{0}, dy{0}, deta{0.05};
double ecrit, speed_of_sound_squared{0.15}, ratio_pressure_energydensity{0.15};
// double Temp, mu_b, mu_q, mu_s ;

// ############ reading and processing the parameters

void readParams(char *filename) {
  char parName[255], parValue[255];
  ifstream fin(filename);
  if (!fin.is_open()) {
    cout << "ERROR: Cannot open config file " << filename << endl;
    exit(1);
  }
  while (fin.good()) {
    string line;
    getline(fin, line);
    istringstream sline(line);
    sline >> parName >> parValue;
    if (strcmp(parName, "surface") == 0)
      strcpy(surface_file, parValue);
    else if (strcmp(parName, "output_dir") == 0)
      strcpy(output_directory, parValue);
    else if (strcmp(parName, "number_of_events") == 0)
      NEVENTS = atoi(parValue);
    else if (strcmp(parName, "shear") == 0)
      shear_viscosity_enabled = atoi(parValue);
    else if (strcmp(parName, "bulk") == 0)
      bulk_viscosity_enabled = atoi(parValue);
    else if (strcmp(parName, "ecrit") == 0)
      ecrit = atof(parValue);
    else if (strcmp(parName, "cs2") == 0)
      speed_of_sound_squared = atof(parValue);
    else if (strcmp(parName, "ratio_pressure_energydensity") == 0)
      ratio_pressure_energydensity = atof(parValue);
    else if (strcmp(parName, "createRootOutput") == 0)
      createRootOutput = atoi(parValue);
    else if (parName[0] == '!')
      cout << "CCC " << sline.str() << endl;
    else
      cout << "UUU " << sline.str() << endl;
  }
}

void printParameters() {
  cout << "======= parameters ===========\n";
  cout << "surface = " << surface_file << endl;
  cout << "output_dir = " << output_directory << endl;
  cout << "number_of_events = " << NEVENTS << endl;
  cout << "shear_visc_on = " << shear_viscosity_enabled << endl;
  cout << "bulk_visc_on = " << bulk_viscosity_enabled << endl;
  cout << "ecrit = " << ecrit << endl;
  cout << "cs2 = " << speed_of_sound_squared << endl;
  cout << "ratio_pressure_energydensity = " << ratio_pressure_energydensity
       << endl;
  cout << "createRootOutput = " << createRootOutput << endl;
  cout << "======= end parameters =======\n";
}

} // namespace params
