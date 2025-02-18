#include <TFile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <string>

#include "build_metadata.h"
#include "gen.h"
#include "oscaroutput.h"
#include "params.h"
#include "tree.h"

using namespace std;
int getNlines(const char *filename);
void print_disclaimer();
void usage(const int exit_status, const std::string &program_name);

using params::number_of_events;
using params::output_directory;
using params::surface_file;

/**
 * Main program
 * samples hadrons from freezeout hypersurface
 *
 * Does command line parsing
 *
 * \param[in] argc Number of arguments on command-line
 * \param[in] argv List of arguments on command-line
 * \return Either 0 or EXIT_FAILURE.
 */
int main(int argc, char *argv[]) {
  constexpr option long_options[] = {{"config", required_argument, 0, 'c'},
                                     {"help", no_argument, 0, 'h'},
                                     {"num", required_argument, 0, 'n'},
                                     {"output", required_argument, 0, 'o'},
                                     {"surface", required_argument, 0, 's'},
                                     {"quiet", no_argument, 0, 'q'},
                                     {"version", no_argument, 0, 0},
                                     {nullptr, 0, 0, 0}};

  // strip any path to program_name
  const std::string program_name =
      std::filesystem::path(argv[0]).filename().native();

  std::string num{"0"};
  std::string configuration, command_line_output_dir{""},
      command_line_surface_file{""};
  bool suppress_disclaimer_and_parameter_printout = false;

  // parse command-line arguments
  int option;
  while ((option = getopt_long(argc, argv, "c:hn:o:s:q", long_options,
                               nullptr)) != -1) {
    switch (option) {
    case 'c':
      configuration = optarg;
      break;
    case 'h':
      usage(EXIT_SUCCESS, program_name);
      break;
    case 'n':
      num = optarg;
      break;
    case 'o':
      command_line_output_dir = optarg;
      break;
    case 's':
      command_line_surface_file = optarg;
      break;
    case 'q':
      suppress_disclaimer_and_parameter_printout = true;
      break;
    case 0: // --version case
      std::printf("%s\n"
#ifdef GIT_BRANCH
                  "Branch   : %s\n"
#endif
                  "System   : %s\nCompiler : %s %s\n"
                  "Date     : %s\n",
                  SAMPLER_VERSION,
#ifdef GIT_BRANCH
                  GIT_BRANCH,
#endif
                  CMAKE_SYSTEM, CMAKE_CXX_COMPILER_ID,
                  CMAKE_CXX_COMPILER_VERSION, BUILD_DATE);
      std::exit(EXIT_SUCCESS);
    default:
      usage(EXIT_FAILURE, program_name);
    }
  }

  // Abort if there are unhandled arguments left.
  if (optind < argc) {
    std::cout << "\n"
              << argv[0] << ": invalid argument -- '" << argv[optind] << "'\n";
    usage(EXIT_FAILURE, program_name);
  }

  params::read_configuration_file(configuration);
  if (!command_line_surface_file.empty()) {
    params::surface_file = command_line_surface_file;
  }
  if (!command_line_output_dir.empty()) {
    params::output_directory = command_line_output_dir;
  }

  if (!suppress_disclaimer_and_parameter_printout) {
    print_disclaimer();
    params::print_config_parameters();
  }

  ROOT::EnableThreadSafety();
  time_t time0;
  time(&time0);
  int ranseed = time0 + std::stoi(num) * 16;

  TRandom3 *random3 = new TRandom3();
  random3->SetSeed(ranseed);
  std::cout << "Random seed:  " << ranseed << std::endl;
  gen::rnd = random3;

  // ========== generator init
  gen::load(surface_file.c_str(), getNlines(surface_file.c_str()));

  // ========== trees & files
  time_t start_time, end_time;
  time(&start_time);

  //============= main task
  std::string make_output_directory = "mkdir -p " + output_directory;
  system(make_output_directory.c_str());

  gen::generate(); // one call for number_of_events

  // ROOT output disabled by default
  if (params::create_root_output) {

    // Initialize ROOT output
    std::string root_output_file = output_directory + "/" + num + ".root";
    TFile *outputFile = new TFile(root_output_file.c_str(), "RECREATE");
    outputFile->cd();
    MyTree *treeIni = new MyTree(static_cast<const char *>("treeini"));

    // Write ROOT output
    for (int iev = 0; iev < number_of_events; iev++) {
      treeIni->fill(iev);
    } // end events loop
    outputFile->Write();
    outputFile->Close();
  }

  // Write Oscar output
  write_oscar_output();

  time(&end_time);
  float execution_time = difftime(end_time, start_time);
  std::cout << "Event generation done (execution time: " << execution_time
            << " [sec])." << std::endl;

  return 0;
}

/// Function to get the number of lines in the freezeout data file
int getNlines(const char *filename) {
  std::ifstream fin(filename);
  if (!fin) {
    std::cerr << "ERROR: getNlines function cannot open freezeout file "
              << filename << std::endl;
    exit(1);
  }
  std::string line;
  int number_of_lines = 0;
  while (getline(fin, line)) {
    number_of_lines++;
  };
  fin.close();
  return number_of_lines;
}

/// Print the disclaimer.
void print_disclaimer() {
  std::cout
      << "###################################################################"
      << "#############"
      << "\n"
      << "\n"
      << " This is SMASH-hadron-sampler version: " << SAMPLER_VERSION << "\n"
      << "\n"
      << " Distributed under the GNU General Public License 3.0"
      << " (GPLv3 or later)."
      << "\n"
      << " See LICENSE file for details."
      << "\n"
      << "\n"
      << " Please cite the following papers when using the SMASH-hadron-sampler"
      << "\n"
      << "      [1] I. Karpenko et al., Phys. Rev. C 91 (2015) 6, 064901"
      << "\n"
      << "      [2] A. Schäfer et al., arXiv:2112.08724"
      << "\n"
      << "\n"
      << " Report issues via GitHub at"
      << "\n"
      << " https://github.com/smash-transport/smash-hadron-sampler"
      << "\n"
      << "\n"
      << "###################################################################"
      << "#############"
      << "\n"
      << "\n";
}

/**
 * Prints usage information and exits the program
 *
 * \param[out] exit_status Exit status to return
 * \param[in] program_name Name of the program
 *
 * usage() is called when either the `--help` or `-h` command line
 * options are given to the program; in this case, the exit status is
 * EXIT_SUCCESS, or when an unknown option is given; in this case,
 * the exit status is EXIT_FAIL.
 */
void usage(const int exit_status, const std::string &program_name) {
  std::printf("\nUsage: %s [option]\n\n", program_name.c_str());
  std::printf(
      "  -h, --help              print this help message and exit\n"
      "\n"
      "  -c, --config <file>     path to input configuration file\n"
      "  -n, --num <int>         specify integer to create random seed "
      "(default: 0)\n"
      "  -o, --output <dir>      override output directory config value\n"
      "  -s, --surface <dir>     override hypersurface freezeout config value\n"
      "\n"
      "  -q, --quiet             suppress disclaimer print-out\n"
      "  --version               print version of sampler executable\n\n");
  std::exit(exit_status);
}
