#ifndef INCLUDE_PARAMS_H_
#define INCLUDE_PARAMS_H_

#include <string>
#include <vector>

namespace params {
extern std::string surface_file, vorticity_file, output_directory;
extern bool bulk_viscosity_enabled, create_root_output, shear_viscosity_enabled,
    spin_sampling_enabled, vorticity_output_enabled;
extern int number_of_events;
extern double dx, dy, deta;
extern double ecrit, speed_of_sound_squared, ratio_pressure_energydensity;
// extern double Temp, mu_b, mu_q, mu_s;
extern std::vector<std::string> comments_in_config_file,
    unknown_parameters_in_config_file;

/**
 * Helper function to get the directory part of a file path.
 */
std::string getDirectory(const std::string& filePath);

void print_config_parameters();
void print_comments_and_unknown_parameters_of_config_file(
    std::vector<std::string> comments_or_unknowns_in_config_file);
void read_configuration_file(const std::string& filename);
}  // namespace params

#endif  // INCLUDE_PARAMS_H_
