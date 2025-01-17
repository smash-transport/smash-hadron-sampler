#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "params.h"

using namespace std;

namespace params {

char sSurface[255], sSpectraDir[255];
bool bulk{false}, createRootOutput{false}, shear{false};
int NEVENTS;
double dx{0}, dy{0}, deta{0.05};
double ecrit, cs2{0.15}, ratio_pressure_energydensity{0.15};
// double Temp, mu_b, mu_q, mu_s ;

// ############ reading and processing the parameters

void readParams(char *filename) {
  char parName[255], parValue[255];
  ifstream fin(filename);
  if (!fin.is_open()) {
    cout << "cannot open parameters file " << filename << endl;
    exit(1);
  }
  while (fin.good()) {
    string line;
    getline(fin, line);
    istringstream sline(line);
    sline >> parName >> parValue;
    if (strcmp(parName, "surface") == 0)
      strcpy(sSurface, parValue);
    else if (strcmp(parName, "output_dir") == 0)
      strcpy(sSpectraDir, parValue);
    else if (strcmp(parName, "number_of_events") == 0)
      NEVENTS = atoi(parValue);
    else if (strcmp(parName, "shear") == 0)
      shear = atoi(parValue);
    else if (strcmp(parName, "bulk") == 0)
      bulk = atoi(parValue);
    else if (strcmp(parName, "ecrit") == 0)
      ecrit = atof(parValue);
    else if (strcmp(parName, "cs2") == 0)
      cs2 = atof(parValue);
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
  cout << "surface = " << sSurface << endl;
  cout << "output_dir = " << sSpectraDir << endl;
  cout << "number_of_events = " << NEVENTS << endl;
  cout << "shear_visc_on = " << shear << endl;
  cout << "bulk_visc_on = " << bulk << endl;
  cout << "ecrit = " << ecrit << endl;
  cout << "cs2 = " << cs2 << endl;
  cout << "ratio_pressure_energydensity = " << ratio_pressure_energydensity
       << endl;
  cout << "createRootOutput = " << createRootOutput << endl;
  cout << "======= end parameters =======\n";
}

} // namespace params
