#include "vorticity.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gen.h"
#include "params.h"

int Vorticity::num_corona_cells_ = -1;

// Check if the vorticity file exists and has the correct format.
// If not, throw an exception.
void Vorticity::ensure_vorticity_file_exists_and_check_format() {
  std::ifstream file(params::sVorticity);
  if (!file.is_open()) {
    throw std::runtime_error(std::string("Error opening file: ") +
                             params::sVorticity);
  }

  std::string line;
  int comment_lines_read = 0;

  // Read and validate the first three comment lines
  while (comment_lines_read < 3 && std::getline(file, line)) {
    if (line.empty() || line[0] != '#') {
      std::ostringstream error;
      error << "Line " << (comment_lines_read + 1) << " in file "
            << params::sVorticity
            << " is not a valid comment line. Expected a line starting with "
               "'#'.\n"
            << "Found: " << line;
      throw std::runtime_error(error.str());
    }

    comment_lines_read++;
  }

  if (comment_lines_read == 3) {
    // Check if the third line matches the expected format
    const std::string expected_third_line =
        "#  τ  x  y  η  dΣ[0]  dΣ[1]  dΣ[2]  dΣ[3]  u[0]  u[1]  u[2]  u[3]  "
        "T  μB  μQ  μS  ∂₀β₀  ∂₀β₁  ∂₀β₂  ∂₀β₃  ∂₁β₀  ∂₁β₁  ∂₁β₂  ∂₁β₃  "
        "∂₂β₀  ∂₂β₁  ∂₂β₂  ∂₂β₃  ∂₃β₀  ∂₃β₁  ∂₃β₂  ∂₃β₃  ϵ";

    if (line != expected_third_line) {
      std::ostringstream error;
      error << "Third comment line in file " << params::sVorticity
            << " does not match the expected format.\n"
            << "Found: " << line << "\n"
            << "Expected: " << expected_third_line;
      throw std::runtime_error(error.str());
    }
  }

  // Ensure that exactly 3 comment lines were read
  if (comment_lines_read < 3) {
    std::ostringstream error;
    error << "File " << params::sVorticity
          << " does not contain the expected 3 comment lines at the beginning.";
    throw std::runtime_error(error.str());
  }

  // Validate the first data line (fourth line)
  if (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') {
      std::ostringstream error;
      error << "First data line (line 4) in file " << params::sVorticity
            << " is empty or incorrectly marked as a comment.";
      throw std::runtime_error(error.str());
    }

    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;

    while (iss >> token) {
      tokens.push_back(token);
    }

    if (tokens.size() != 33) {
      std::ostringstream error;
      error << "First data line in file " << params::sVorticity
            << " does not contain exactly 33 values. Found " << tokens.size()
            << " values.";
      throw std::runtime_error(error.str());
    }
  } else {
    std::ostringstream error;
    error << "File " << params::sVorticity
          << " does not contain a valid data line.";
    throw std::runtime_error(error.str());
  }

  file.close();
}

// Set the number of corona cells in Vorticity from
// the second comment line of the vorticity file.
void Vorticity::set_number_of_corona_cells() {
  // Open the vorticity file
  std::ifstream file(params::sVorticity);
  if (!file) {
    std::ostringstream error;
    error << "Failed to open vorticity file: " << params::sVorticity;
    throw std::runtime_error(error.str());
  }

  std::string line;
  // Skip the first line
  if (!std::getline(file, line)) {
    throw std::runtime_error(
        "Vorticity file is empty or does not contain enough lines.");
  }

  // Read the second line
  if (!std::getline(file, line)) {
    throw std::runtime_error("Vorticity file does not contain a second line.");
  }

  // Check if the line is a comment
  if (line.empty() || line.front() != '#') {
    throw std::runtime_error(
        "Expected a comment line starting with '#' but found: " + line);
  }

  // Find the last word in the line and attempt to parse it as an integer
  std::istringstream line_stream(line);
  std::string word, last_word;
  while (line_stream >> word) {
    last_word = word;
  }

  // Check if the last word is a valid integer
  if (last_word.empty() ||
      !std::all_of(last_word.begin(), last_word.end(), ::isdigit)) {
    throw std::runtime_error(
        "Failed to parse an integer from the second comment line: " + line);
  }

  // Convert the last word to an integer and set the number of corona cells
  num_corona_cells_ = std::stoi(last_word);

  file.close();
}

// Given the complete freezeout surface, this function sets the vorticity tensor
// in all surface cells from the vorticity file.
void Vorticity::set_vorticity_and_energy_in_surface_cells(gen::element* surf,
                                                          int N) {
  // Read the vorticity file line by line
  std::ifstream file(params::sVorticity);
  std::string line;

  // Skip comment lines (lines starting with '#')
  while (std::getline(file, line)) {
    if (!line.empty() && line[0] != '#') {
      break;  // Found the start of data
    }
  }

  // Ensure the first non-comment line exists
  if (line.empty()) {
    std::cerr << "Error: Vorticity file contains no data after comments."
              << std::endl;
    return;
  }

  // Ensure that the number of corona cells has been set
    if (num_corona_cells_ < 0) {
        std::cerr << "Error: Number of corona cells not set in Vorticity."
                << std::endl;
        return;
    }

  // Process corona cells (vorticity set to 0)
  for (int i = 0; i < num_corona_cells_; i++) {
    if (line.empty() || !std::getline(file, line)) {
      std::cerr
          << "Error: Unexpected end of file while processing corona cells."
          << std::endl;
      return;
    }

    // Initialize vorticity for this corona cell.
    // The constructor already sets all components to 0.
    surf[i].vorticity = std::make_unique<Vorticity>();
  }

  // Process freezeout cells (set the actual vorticity tensor values)
  for (int i = 0; i < N; i++) {
    if (line.empty() || !std::getline(file, line)) {
      std::cerr
          << "Error: Unexpected end of file while processing freezeout cells."
          << std::endl;
      return;
    }

    // Initialize vorticity for this freezeout cell
    surf[num_corona_cells_ + i].vorticity = std::make_unique<Vorticity>();

    // Parse and set the relevant components (16th to 32nd values)
    std::istringstream instream(line);
    double value;
    for (int column = 0; column <= 32; column++) {
      instream >> value;
      if (column >= 16 && column <= 31) {
        (**surf[num_corona_cells_ + i].vorticity)[column - 16] = value;
      } else if (column == 32) {
        // Set the energy density for this cell
        surf[num_corona_cells_ + i].e = value;
      }
    }
  }

  // Check for unexpected extra lines
  if (std::getline(file, line)) {
    std::cerr << "Warning: Extra lines detected in the vorticity file after "
                 "processing all cells."
              << std::endl;
  }

  // Close the file after processing
  file.close();
}
