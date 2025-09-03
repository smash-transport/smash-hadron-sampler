#include "vorticity.h"

#include <fstream>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "gen.h"
#include "params.h"
#include "vorticity.h"

int Vorticity::num_corona_cells_ = -1;

// Check if the vorticity file exists and has the correct format.
// If not, throw an exception.
void Vorticity::ensure_vorticity_file_exists_and_check_format() {
  std::ifstream file(params::vorticity_file);
  if (!file.is_open()) {
    throw std::runtime_error(std::string("Error opening file: ") +
                             params::surface_file);
  }

  std::string line;
  int comment_lines_read = 0;

  // Read and validate the first three comment lines
  while (comment_lines_read < 3 && std::getline(file, line)) {
    // Check if the line is a comment line
    if (line.empty() || line[0] != '#') {
      std::ostringstream error;
      error << "Line " << (comment_lines_read + 1) << " in file "
            << params::vorticity_file
            << " is not a valid comment line. Expected a line starting with "
               "'#'.\n"
            << "Found: " << line;
      throw std::invalid_argument(error.str());
    }

    // Check the second comment line for the expected format
    if (comment_lines_read == 1) {
      std::regex regex("# Number of corona cells: \\d+");
      if (!std::regex_match(line, regex)) {
        std::ostringstream error;
        error << "Second comment line in file " << params::vorticity_file
              << " does not match the expected format.\n"
              << "Found: " << line << "\n"
              << "Expected: '# Number of corona cells: <integer>'";
        throw std::domain_error(error.str());
      }
    }

    // Check if the third line matches the expected format
    if (comment_lines_read == 2) {
      const std::string expected_third_line =
          "#  τ  x  y  η  dΣ[0]  dΣ[1]  dΣ[2]  dΣ[3]  u[0]  u[1]  u[2]  u[3]  "
          "T  μB  μQ  μS  ∂₀β₀  ∂₀β₁  ∂₀β₂  ∂₀β₃  ∂₁β₀  ∂₁β₁  ∂₁β₂  ∂₁β₃  "
          "∂₂β₀  ∂₂β₁  ∂₂β₂  ∂₂β₃  ∂₃β₀  ∂₃β₁  ∂₃β₂  ∂₃β₃  ϵ";

      if (line != expected_third_line) {
        std::ostringstream error;
        error << "Third comment line in file " << params::vorticity_file
              << " does not match the expected format.\n"
              << "Found: " << line << "\n"
              << "Expected: " << expected_third_line;
        throw std::domain_error(error.str());
      }
    }

    comment_lines_read++;
  }

  // Ensure that exactly 3 comment lines were read
  if (comment_lines_read < 3) {
    std::ostringstream error;
    error << "File " << params::vorticity_file
          << " does not contain the expected 3 comment lines at the beginning.";
    throw std::invalid_argument(error.str());
  }

  // Validate the first data line (fourth line)
  if (std::getline(file, line)) {
    if (line.empty() || line[0] == '#') {
      std::ostringstream error;
      error << "First data line (line 4) in file " << params::vorticity_file
            << " is empty or a comment line. Only 3 comment lines are expected";
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
      error << "First data line in file " << params::vorticity_file
            << " does not contain exactly 33 values. Found " << tokens.size()
            << " values.";
      throw std::length_error(error.str());
    }
  } else {
    std::ostringstream error;
    error << "File " << params::vorticity_file
          << " does not contain a valid data line.";
    throw std::runtime_error(error.str());
  }

  file.close();
}

void Vorticity::ensure_extended_freezeout_is_given() {
  std::ifstream surface_file(params::surface_file);
  if (!surface_file.is_open()) {
    throw std::runtime_error(std::string("Error opening file: ") +
                             params::vorticity_file);
  }
  std::string line;

  while (std::getline(surface_file, line)) {
    if (line[0] != '#') {
      // Even though there are no comment lines in the vhlle output, we include
      // it in case of potential header lines in the future.
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::string token;

      while (iss >> token) {
        tokens.push_back(token);
      }
      // Check for 29 columns in the freezeout file which corresponds to the
      // extended output format
      if (tokens.size() != 29) {
        std::ostringstream error;
        error
            << "First data line in file " << params::surface_file
            << " does not contain exactly 29 values. Found " << tokens.size()
            << " values. Either the extended format was not used or the number "
            << "of columns in the extended freezeout has changed.";
        throw std::runtime_error(error.str());
      }
      break;
    } else {
      continue;
    }
  }
}

// Set the number of corona cells in Vorticity from
// the second comment line of the vorticity file.
void Vorticity::set_number_of_corona_cells() {
  // Open the vorticity file
  std::ifstream file(params::vorticity_file);
  if (!file) {
    std::ostringstream error;
    error << "Failed to open vorticity file: " << params::vorticity_file;
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
void Vorticity::set_vorticity_in_surface_cells(gen::element* surf, int N) {
  // Read the vorticity file line by line
  std::ifstream file(params::vorticity_file);
  std::string line;

  const int num_freezeout_cells = N - num_corona_cells_;
  int num_comment_lines = 0;

  // Ensure that the number of corona cells has been set
  if (num_corona_cells_ < 0) {
    throw std::runtime_error(
        "Error: Number of corona cells not set in Vorticity.");
  }

  // Process corona cells (vorticity set to 0)
  for (int i = 0; i < num_corona_cells_; i++) {
    // Initialize vorticity for this corona cell.
    // The constructor already sets all components to 0.
    surf[i].vorticity = std::make_unique<Vorticity>();

    // TODO: Set energy density in corona cells
  }

  // Skip comment lines in vorticity file (lines starting with '#')
  while (std::getline(file, line)) {
    if (line.empty()) {
      throw std::runtime_error(
          "Error: Unexpected empty line encountered in head of the vorticity "
          "file.");
    } else if (line[0] != '#') {
      break;  // Found the start of data
    }
    num_comment_lines++;
  }

  // Process freezeout cells (set the actual vorticity tensor values)
  for (int i = 0; i < num_freezeout_cells; i++) {
    // Check if the current line is valid
    if (line.empty()) {
      const int current_line = i + 1 + num_comment_lines;
      throw std::runtime_error(
          "Error: Encountered an empty line in the file \"" +
          std::string(params::vorticity_file) + "\" at line " +
          std::to_string(current_line));
    }

    // Initialize vorticity for this freezeout cell
    surf[num_corona_cells_ + i].vorticity = std::make_unique<Vorticity>();

    // Parse and set the relevant components (16th to 32nd values)
    std::istringstream instream(line);
    double value;
    for (int column = 0; column <= 31; column++) {
      instream >> value;
      if (column >= 16 && column <= 31) {
        (**surf[num_corona_cells_ + i].vorticity)[column - 16] = value;
      }
    }

    // Read the next line for the next iteration, if not the last iteration
    if (i < num_freezeout_cells - 1) {
      if (!std::getline(file, line)) {
        // Next line has already been read, therefore the iterator is i+2
        const int current_line = i + 2 + num_comment_lines;
        throw std::runtime_error("Error: Unexpected end in the file \"" +
                                 std::string(params::vorticity_file) +
                                 "\" at line " + std::to_string(current_line));
      }
    }
  }

  // Check for unexpected extra lines
  if (std::getline(file, line)) {
    throw std::runtime_error(
        "Warning: Extra lines detected in the vorticity file after processing "
        "all cells.");
  }

  // Close the file after processing
  file.close();
}

// Boost the vorticity tensor to the fluid rest frame with the given boost
// matrix
void Vorticity::boost_vorticity_to_fluid_rest_frame(
    const FourMatrix& boostMatrix) {
  // Transform the row and column indices for a 4x4 matrix into a 1d index
  auto at = [](int i, int j) { return (i * 4 + j); };

  // Create a copy of the vorticity tensor
  std::array<double, 16> vorticity_copy = vorticity_;

  // Boost the vorticity tensor to the fluid rest frame
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      double vorticity_ij = 0.0;
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          vorticity_ij +=
               boostMatrix[i][k] * boostMatrix[j][l] * vorticity_copy[at(k, l)];
        }
      }
      vorticity_[at(i, j)] = vorticity_ij;
    }
  }
}
