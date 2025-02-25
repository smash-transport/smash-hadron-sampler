#include <fstream>
#include <iostream>
#include <sstream>

#include "params.h"

using namespace std;

namespace params {

std::string surface_file{"unset"}, output_directory{"unset"};
bool bulk_viscosity_enabled{false}, create_root_output{false},
    shear_viscosity_enabled{false};
int number_of_events;
double dx{0}, dy{0}, deta{0.05};
double ecrit, speed_of_sound_squared{0.15}, ratio_pressure_energydensity{0.15};
// double Temp, mu_b, mu_q, mu_s ;

/// Read and process configuration file parameters
void read_configuration_file(const std::string &filename) {
  std::string parName, parValue;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) {
    std::cerr << "ERROR: Cannot open config file " << filename << std::endl;
    exit(1);
  }
  while (fin.good()) {
    std::string line;
    getline(fin, line);
    std::istringstream sline(line);
    sline >> parName >> parValue;
    if (parName == "surface_file") {
      surface_file = parValue;
    } else if (parName == "output_dir") {
      output_directory = parValue;
    } else if (parName == "number_of_events") {
      number_of_events = std::stoi(parValue);
    } else if (parName == "shear") {
      shear_viscosity_enabled = std::stoi(parValue);
    } else if (parName == "bulk") {
      bulk_viscosity_enabled = std::stoi(parValue);
    } else if (parName == "ecrit") {
      ecrit = std::stod(parValue);
    } else if (parName == "cs2") {
      speed_of_sound_squared = std::stod(parValue);
    } else if (parName == "ratio_pressure_energydensity") {
      ratio_pressure_energydensity = std::stod(parValue);
    } else if (parName == "create_root_output") {
      create_root_output = std::stoi(parValue);
    } else if (parName[0] == '!') {
      std::cout << "Comment in config: '" << sline.str() << "'" << std::endl;
    } else {
      std::cout << "Unknown config parameter: '" << sline.str() << "'"
                << std::endl;
    }
  }
}

/// Print config parameters to terminal
void print_config_parameters() {
  std::cout << "# ---------------------- parameters used for sampler run "
               "-----------------------"
            << "\n"
            << "surface_file:                 " << surface_file << "\n"
            << "output_dir:                   " << output_directory << "\n"
            << "number_of_events:             " << number_of_events << "\n"
            << "shear:                        " << shear_viscosity_enabled
            << "\n"
            << "bulk:                         " << bulk_viscosity_enabled
            << "\n"
            << "ecrit:                        " << ecrit << "\n"
            << "cs2:                          " << speed_of_sound_squared
            << "\n"
            << "ratio_pressure_energydensity: " << ratio_pressure_energydensity
            << "\n"
            << "create_root_output:           " << create_root_output << "\n"
            << "# -------------------------------------------------------"
               "-----------------------"
            << "\n"
            << "\n";
}

} // namespace params
