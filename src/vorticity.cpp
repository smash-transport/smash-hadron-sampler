#include "vorticity.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "gen.h"
#include "params.h"

// Check if the vorticity file exists. If not, throw an exception.
void Vorticity::ensure_vorticity_file_exists_and_contains_16_values() {
  std::ifstream file(params::sVorticity);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << params::sVorticity << std::endl;
  }
  // Check if the first line contains exactly 16 values
  std::string line;
  if (std::getline(file, line)) {
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;

    // Split the line into tokens based on whitespace
    while (iss >> token) {
      tokens.push_back(token);
    }

    // Check if the number of tokens is exactly 16
    if (tokens.size() == 16) {
      file.close();
    } else {
      std::cerr << "First line in file " << params::sVorticity
                << " does not contain exactly 16 values." << std::endl;
    }
  } else {
    std::cerr << "File " << params::sVorticity
              << " is empty or could not be read." << std::endl;
  }

  file.close();
}

// Given the complete freezeout surface, this function sets the vorticity tensor
// in all surface cells from the vorticity file.
void Vorticity::set_vorticity_in_all_surface_cells(element* surf, int N) {
  // Open the vorticity file
  std::ifstream file(params::sVorticity);
  if (!file.is_open()) {
    std::cerr << "Error opening file: " << params::sVorticity << std::endl;
  }
  // Read the vorticity file line by line
  std::string line;
  std::istringstream instream;
  for (int n = 0; n < N; n++) {
    if (std::getline(file, line)) {
      // Initialize Vorticity
      surf[n].vorticity = std::make_unique<Vorticity>();

      // Set the elements of the vorticity tensor from the line
      instream.str(line);
      instream.seekg(0);
      instream.clear();
      for (int column = 0; column < 16; column++) {
        instream >> (**surf[n].vorticity)[column];
      }
    } else {
      std::cerr << "Not enough lines in the vorticity file." << std::endl;
      break;
    }
  }
  file.close();
}