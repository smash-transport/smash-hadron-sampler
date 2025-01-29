#include <TFile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <fstream>
#include <limits>
#include <string>

#include "gen.h"
#include "oscaroutput.h"
#include "params.h"
#include "tree.h"

using namespace std;
int get_N_lines_and_set_cell_sizes(const char *filename);
int readCommandLine(int argc, char **argv);

using params::NEVENTS;
using params::output_directory;
using params::surface_file;

// ########## MAIN block ##################

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();
  // command-line parameters
  int prefix = readCommandLine(argc, argv);
  params::printParameters();
  time_t time0;
  time(&time0);
  int ranseed = time0 + prefix * 16;

  TRandom3 *random3 = new TRandom3();
  random3->SetSeed(ranseed);
  cout << "Random seed = " << ranseed << endl;
  gen::rnd = random3;

  // ========== generator init
  gen::load(surface_file.c_str(),
            get_N_lines_and_set_cell_sizes(surface_file.c_str()));

  // ========== trees & files
  time_t start, end;
  time(&start);

  //============= main task
  std::string make_output_directory = "mkdir -p " + output_directory;
  system(make_output_directory.c_str());

  gen::generate(); // one call for NEVENTS

  // ROOT output disabled by default
  if (params::createRootOutput) {

    // Initialize ROOT output
    std::string root_output_file =
        output_directory + "/" + std::to_string(prefix) + ".root";
    TFile *outputFile = new TFile(root_output_file.c_str(), "RECREATE");
    outputFile->cd();
    MyTree *treeIni = new MyTree(static_cast<const char *>("treeini"));

    // Write ROOT output
    for (int iev = 0; iev < NEVENTS; iev++) {
      treeIni->fill(iev);
    } // end events loop
    outputFile->Write();
    outputFile->Close();
  }

  // Write Oscar output
  write_oscar_output();

  cout << "Event generation done\n";
  time(&end);
  float diff2 = difftime(end, start);
  cout << "Execution time = " << diff2 << " [sec]" << endl;
  return 0;
}

int readCommandLine(int argc, char **argv) {
  if (argc == 1) {
    cout << "ERROR: Missing command line parameters!" << endl;
    exit(1);
  }
  bool is_config_given = false;
  int prefix = 0;
  int iarg = 1;
  while (iarg < argc - 1) {
    if (strcmp(argv[iarg], "--config") == 0 || strcmp(argv[iarg], "-c") == 0) {
      params::readParams(argv[iarg + 1]);
      is_config_given = true;
      iarg += 2;
    } else if (strcmp(argv[iarg], "--num") == 0 ||
               strcmp(argv[iarg], "-n") == 0) {
      prefix = atoi(argv[iarg + 1]);
      iarg += 2;
    } else if (strcmp(argv[iarg], "--output") == 0 ||
               strcmp(argv[iarg], "-o") == 0) {
      output_directory = argv[iarg + 1];
      iarg += 2;
    } else if (strcmp(argv[iarg], "--surface") == 0 ||
               strcmp(argv[iarg], "-s") == 0) {
      surface_file = argv[iarg + 1];
      iarg += 2;
    } else {
      cout << "Unknown command line parameter: " << argv[iarg] << endl;
      iarg++;
    }
  }
  if (!is_config_given) {
    cout << "ERROR: No config file provided." << endl;
    exit(126);
  }
  return prefix;
}

// auxiliary function to get the number of lines of the freezeout data file
// and to read out the cell sizes of the hydro evolution
int get_N_lines_and_set_cell_sizes(const char *filename) {
  ifstream fin(filename);
  if (!fin) {
    cout << "get_N_lines_and_set_cell_sizes function error: Cannot open file "
         << filename << endl;
    exit(1);
  }
  std::string line;
  double x_value_prior_line, y_value_prior_line, eta_z_value_prior_line;
  double time_value_current_line, x_value_current_line, y_value_current_line,
      eta_z_value_current_line;
  double machine_epsilon = std::numeric_limits<double>::epsilon();
  int number_of_lines = 0;
  while (getline(fin, line)) {
    std::istringstream sline(line);
    sline >> time_value_current_line >> x_value_current_line >>
        y_value_current_line >> eta_z_value_current_line;
    if (params::transversal_smearing_enabled) {
      if (number_of_lines == 0) {
        x_value_prior_line = x_value_current_line;
        y_value_prior_line = y_value_current_line;
      }
      if (params::dx < machine_epsilon) {
        params::dx = std::abs(x_value_current_line - x_value_prior_line);
        x_value_prior_line = x_value_current_line;
      }
      if (params::dy < machine_epsilon) {
        params::dy = std::abs(y_value_current_line - y_value_prior_line);
        y_value_prior_line = y_value_current_line;
      }
    }
    if (params::deta_dz < machine_epsilon) {
      if (number_of_lines == 0) {
        eta_z_value_prior_line = eta_z_value_current_line;
      }
      params::deta_dz =
          std::abs(eta_z_value_current_line - eta_z_value_prior_line);
      eta_z_value_prior_line = eta_z_value_current_line;
    }
    number_of_lines++;
  };
  fin.close();
  return number_of_lines;
}
