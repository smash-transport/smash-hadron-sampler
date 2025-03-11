#include "vorticity.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

#include "params.h"
#include "virtest/vir/test.h"

namespace test_utils {

namespace fs = std::filesystem;

// Create a temporary directory and return its path
fs::path create_temp_directory(const std::string& dir_name = "test_temp") {
  fs::path tempDir = fs::current_path() / dir_name;
  fs::create_directory(tempDir);

  if (!fs::exists(tempDir) || !fs::is_directory(tempDir)) {
    throw std::runtime_error("Failed to create the temporary directory: " +
                             tempDir.string());
  }

  return tempDir;
}

// Write content to a file
void write_to_file(const fs::path& file_path, const std::string& content) {
  // Open the file in append mode (std::ios::app)
  std::ofstream file(file_path, std::ios::app);
  if (!file) {
    throw std::ios_base::failure("Failed to open file for writing: " +
                                 file_path.string());
  }
  // Write content followed by a newline
  file << content << std::endl;
  file.close();
}

void set_params(const fs::path& configFilePath) {
  std::string configFilePathStr = configFilePath.string();
  char* configFilePathCStr = const_cast<char*>(configFilePathStr.c_str());
  params::readParams(configFilePathCStr);
}

// Cleanup the temporary directory and all its contents
void cleanup_temp_directory(const fs::path& tempDir) {
  try {
    fs::remove_all(tempDir);
  } catch (const std::exception& e) {
    throw std::runtime_error("Cleanup failed for directory: " +
                             tempDir.string() + " with error: " + e.what());
  }
}

}  // namespace test_utils

// Custom function to compare two doubles within a given tolerance
bool expect_near(double val1, double val2, double abs_error) {
  return std::abs(val1 - val2) <= abs_error;
}

TEST(vorticity_default_constructor) {
  Vorticity vorticity;
  for (int i = 0; i < 16; ++i) {
    VERIFY(expect_near(vorticity[i], 0.0, 1e-9));
  }
}

TEST(vorticity_set_vorticity) {
  Vorticity vorticity;
  std::array<double, 16> new_vorticity = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,
                                          7.0,  8.0,  9.0,  10.0, 11.0, 12.0,
                                          13.0, 14.0, 15.0, 16.0};
  vorticity.set_vorticity(new_vorticity);
  for (int i = 0; i < 16; ++i) {
    VERIFY(expect_near(vorticity[i], new_vorticity[i], 1e-9));
  }
}

TEST(vorticity_get_vorticity) {
  Vorticity vorticity;
  std::array<double, 16> new_vorticity = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,
                                          7.0,  8.0,  9.0,  10.0, 11.0, 12.0,
                                          13.0, 14.0, 15.0, 16.0};
  vorticity.set_vorticity(new_vorticity);
  std::array<double, 16> vorticity_array = vorticity.get_vorticity();
  for (size_t i = 0; i < 16; ++i) {
    VERIFY(expect_near(vorticity_array[i], new_vorticity[i], 1e-9));
  }
}

TEST(vorticity_at) {
  Vorticity vorticity;
  std::array<double, 16> new_vorticity = {1.0,  2.0,  3.0,  4.0,  5.0,  6.0,
                                          7.0,  8.0,  9.0,  10.0, 11.0, 12.0,
                                          13.0, 14.0, 15.0, 16.0};
  vorticity.set_vorticity(new_vorticity);
  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j < 4; ++j) {
      VERIFY(expect_near(vorticity.at(i, j), new_vorticity[i * 4 + j], 1e-9));
    }
  }
}

TEST(vorticity_file_does_not_exist) {
  namespace fs = std::filesystem;

  // Create temporary directory
  fs::path tempDir = test_utils::create_temp_directory();

  // Create paths for the config file and the vorticity file
  fs::path configFilePath = tempDir / "config_temp";
  fs::path betaFilePath = tempDir / "beta.dat";

  // Write config file
  test_utils::write_to_file(configFilePath,
                            "dbeta " + betaFilePath.string() + "\n");

  // Pass config to params
  test_utils::set_params(configFilePath);

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with non-existing beta.dat file");
  } catch (const std::runtime_error& e) {
    // Error was caught as expected
  }

  // Cleanup
  try {
    test_utils::cleanup_temp_directory(tempDir);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Cleanup failed: ") + e.what());
  }
}

TEST(vorticity_file_too_few_comment_lines) {
  namespace fs = std::filesystem;

  // Create temporary directory
  fs::path tempDir = test_utils::create_temp_directory();

  // Create paths for the config file and the vorticity file
  fs::path configFilePath = tempDir / "config_temp";
  fs::path betaFilePath = tempDir / "beta.dat";
  fs::path freezeoutFilePath = tempDir / "freezeout.dat";

  // Write config file (Here we use the keyword "surface" instead of "dbeta",
  // as the sampler looks for the beta.dat file in the same directory as the
  // surface file as default. This is ensured here)
  test_utils::write_to_file(configFilePath,
                            "surface " + freezeoutFilePath.string() + "\n");

  // Write beta.dat file with one comment line and expect error
  test_utils::write_to_file(betaFilePath, "# Header Line 1");

  // Pass config to params
  test_utils::set_params(configFilePath);

  // Ensure the beta.dat file exists and was set from default
  VERIFY(fs::exists(params::sVorticity));

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with too few comment lines");
  } catch (const std::invalid_argument& e) {
    // Error was caught as expected
  }

  // Add a second comment line to the beta.dat file and expect error
  test_utils::write_to_file(betaFilePath, "# Number of corona cells: 1");

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with too few comment lines");
  } catch (const std::invalid_argument& e) {
    // Error was caught as expected
  }

  // Add a third comment line not matching the expected line and expect error
  test_utils::write_to_file(betaFilePath, "# Header Line 3");

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with incorrect third comment line");
  } catch (const std::domain_error& e) {
    // Error was caught as expected
  }

  // Cleanup
  try {
    test_utils::cleanup_temp_directory(tempDir);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Cleanup failed: ") + e.what());
  }
}

TEST(vorticity_file_invalid_second_header_line) {
  namespace fs = std::filesystem;

  // Create temporary directory
  fs::path tempDir = test_utils::create_temp_directory();

  // Create paths for the config file and the vorticity file
  fs::path configFilePath = tempDir / "config_temp";
  fs::path betaFilePath = tempDir / "beta.dat";

  // Write config file
  test_utils::write_to_file(configFilePath,
                            "dbeta " + betaFilePath.string() + "\n");

  // Write beta.dat file with invalid second comment line and expect error
  test_utils::write_to_file(betaFilePath, "# Header Line 1");
  test_utils::write_to_file(betaFilePath, "# Header Line");

  // Pass config to params
  test_utils::set_params(configFilePath);

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with invalid second comment line");
  } catch (const std::domain_error& e) {
    // Error was caught as expected
  }

  // Delete beta file and write a valid second comment line
  fs::remove(betaFilePath);
  test_utils::write_to_file(betaFilePath, "# Header Line 1");
  test_utils::write_to_file(betaFilePath, "# Number of corona cells: 1.3");

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with invalid second comment line");
  } catch (const std::domain_error& e) {
    // Error was caught as expected
  }

  // Cleanup
  try {
    test_utils::cleanup_temp_directory(tempDir);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Cleanup failed: ") + e.what());
  }
}

TEST(vorticity_file_too_many_comment_lines) {
  namespace fs = std::filesystem;

  // Create temporary directory
  fs::path tempDir = test_utils::create_temp_directory();

  // Create paths for the config file and the vorticity file
  fs::path configFilePath = tempDir / "config_temp";
  fs::path betaFilePath = tempDir / "beta.dat";

  // Write config file
  test_utils::write_to_file(configFilePath,
                            "dbeta " + betaFilePath.string() + "\n");

  // Write beta.dat file with four comment lines and expect error
  test_utils::write_to_file(betaFilePath, "# Header Line 1");
  test_utils::write_to_file(betaFilePath, "# Number of corona cells: 2");
  test_utils::write_to_file(
      betaFilePath,
      "# τ  x  y  η  dΣ[0]  dΣ[1]  dΣ[2]  dΣ[3]  u[0]  u[1]  u[2]  u[3]  "
      "T  μB  μQ  μS  ∂₀β₀  ∂₀β₁  ∂₀β₂  ∂₀β₃  ∂₁β₀  ∂₁β₁  ∂₁β₂  ∂₁β₃  "
      "∂₂β₀  ∂₂β₁  ∂₂β₂  ∂₂β₃  ∂₃β₀  ∂₃β₁  ∂₃β₂  ∂₃β₃  ϵ");
  test_utils::write_to_file(betaFilePath, "# Header Line 4");

  // Pass config to params
  test_utils::set_params(configFilePath);

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly succeeded "
        "with too many comment lines");
  } catch (const std::runtime_error& e) {
    // Error was caught as expected
  }

  // Cleanup
  try {
    test_utils::cleanup_temp_directory(tempDir);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Cleanup failed: ") + e.what());
  }
}

TEST(vorticity_file_valid) {
  namespace fs = std::filesystem;

  // Create temporary directory
  fs::path tempDir = test_utils::create_temp_directory();

  // Create paths for the config file and the vorticity file
  fs::path configFilePath = tempDir / "config_temp";
  fs::path betaFilePath = tempDir / "beta.dat";

  // Write config file
  test_utils::write_to_file(configFilePath,
                            "dbeta " + betaFilePath.string() + "\n");

  // Write beta.dat file with the expected comment lines
  test_utils::write_to_file(betaFilePath, "# Header Line 1");
  test_utils::write_to_file(betaFilePath, "# Number of corona cells: 2");
  test_utils::write_to_file(
      betaFilePath,
      "# τ  x  y  η  dΣ[0]  dΣ[1]  dΣ[2]  dΣ[3]  u[0]  u[1]  u[2]  u[3]  "
      "T  μB  μQ  μS  ∂₀β₀  ∂₀β₁  ∂₀β₂  ∂₀β₃  ∂₁β₀  ∂₁β₁  ∂₁β₂  ∂₁β₃  "
      "∂₂β₀  ∂₂β₁  ∂₂β₂  ∂₂β₃  ∂₃β₀  ∂₃β₁  ∂₃β₂  ∂₃β₃  ϵ");
  test_utils::write_to_file(betaFilePath,
                            "0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 "
                            "8.0 9.0 10.0 11.0 12.0 13.0 14.0 "
                            "15.0 16.0 17.0 18.0 19.0 20.0 21.0 "
                            "22.0 23.0 24.0 25.0 26.0 27.0 28.0 "
                            "29.0 30.0 31.0 32.0 \n");

  // Pass config to params
  test_utils::set_params(configFilePath);

  try {
    Vorticity::ensure_vorticity_file_exists_and_check_format();
  } catch (const std::exception& e) {
    throw std::runtime_error(
        "ensure_vorticity_file_exists_and_check_format unexpectedly failed "
        "with valid beta.dat file");
  }

  // Cleanup
  try {
    test_utils::cleanup_temp_directory(tempDir);
  } catch (const std::exception& e) {
    throw std::runtime_error(std::string("Cleanup failed: ") + e.what());
  }
}
