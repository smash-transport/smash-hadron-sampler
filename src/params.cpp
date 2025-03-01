#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

#include "params.h"

using namespace std;

namespace params {

std::string surface_file{"unset"}, output_directory{"unset"},
    hydro_coordinate_system{"milne"};
bool bulk_viscosity_enabled{false}, createRootOutput{false},
    shear_viscosity_enabled{false}, transversal_smearing_enabled{false};
int NEVENTS;
double dx{0}, dy{0}, deta_dz{0};
double ecrit, speed_of_sound_squared{0.15}, ratio_pressure_energydensity{0.15};
// double Temp, mu_b, mu_q, mu_s ;

// ############ reading and processing the parameters

void readParams(const std::string &filename) {
  std::string parName, parValue;
  ifstream fin(filename.c_str());
  if (!fin.is_open()) {
    cout << "ERROR: Cannot open config file " << filename << endl;
    exit(1);
  }
  while (fin.good()) {
    std::string line;
    getline(fin, line);
    istringstream sline(line);
    sline >> parName >> parValue;
    if (parName == "surface_file")
      surface_file = parValue;
    else if (parName == "output_dir")
      output_directory = parValue;
    else if (parName == "number_of_events")
      NEVENTS = std::stoi(parValue);
    else if (parName == "shear")
      shear_viscosity_enabled = std::stoi(parValue);
    else if (parName == "bulk")
      bulk_viscosity_enabled = std::stoi(parValue);
    else if (parName == "ecrit")
      ecrit = std::stod(parValue);
    else if (parName == "cs2")
      speed_of_sound_squared = std::stod(parValue);
    else if (parName == "ratio_pressure_energydensity")
      ratio_pressure_energydensity = std::stod(parValue);
    else if (parName == "createRootOutput")
      createRootOutput = std::stoi(parValue);
    else if (parName == "hydro_coordinate_system") {
      hydro_coordinate_system = parValue;
      for (unsigned int i = 0; i < hydro_coordinate_system.size(); i++) {
        hydro_coordinate_system[i] = tolower(hydro_coordinate_system[i]);
      }
      if (hydro_coordinate_system != "milne" &&
          hydro_coordinate_system != "cartesian") {
        throw std::invalid_argument(
            std::string(
                "Only 'Milne' or 'Cartesian' are supported as coordinate "
                "systems for the hydro evolution. Provided was '") +
            hydro_coordinate_system +
            std::string("'. Please update the config and try again."));
      }
    } else if (parName == "transversal_smearing") {
      transversal_smearing_enabled = std::stoi(parValue);
    } else if (parName[0] == '!')
      cout << "CCC " << sline.str() << endl;
    else
      cout << "UUU " << sline.str() << endl;
  }
}

void printParameters() {
  cout << "======= parameters ===========\n";
  cout << "surface_file = " << surface_file << endl;
  cout << "output_dir = " << output_directory << endl;
  cout << "number_of_events = " << NEVENTS << endl;
  cout << "shear_visc_on = " << shear_viscosity_enabled << endl;
  cout << "bulk_visc_on = " << bulk_viscosity_enabled << endl;
  cout << "ecrit = " << ecrit << endl;
  cout << "cs2 = " << speed_of_sound_squared << endl;
  cout << "ratio_pressure_energydensity = " << ratio_pressure_energydensity
       << endl;
  cout << "createRootOutput = " << createRootOutput << endl;
  cout << "hydro_coordinate_system = " << hydro_coordinate_system << endl;
  cout << "transversal_smearing = " << transversal_smearing_enabled << endl;
  cout << "======= end parameters =======\n";
}

} // namespace params
