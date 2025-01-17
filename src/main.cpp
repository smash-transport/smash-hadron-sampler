#include <TFile.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <fstream>

#include "gen.h"
#include "oscaroutput.h"
#include "params.h"
#include "tree.h"

using namespace std;
int getNlines(char *filename);
int readCommandLine(int argc, char **argv);

using params::NEVENTS;
using params::output_directory;
using params::surface_file;

int ranseed;

extern "C" {
void getranseedcpp_(int *seed) { *seed = ranseed; }
}

// ########## MAIN block ##################

int main(int argc, char **argv) {
  ROOT::EnableThreadSafety();
  // command-line parameters
  int prefix = readCommandLine(argc, argv);
  params::printParameters();
  time_t time0;
  time(&time0);
  ranseed = time0 + prefix * 16;

  TRandom3 *random3 = new TRandom3();
  random3->SetSeed(ranseed);
  cout << "Random seed = " << ranseed << endl;
  gen::rnd = random3;

  // ========== generator init
  gen::load(surface_file, getNlines(surface_file));

  // ========== trees & files
  time_t start, end;
  time(&start);

  //============= main task
  char sbuffer[255];
  sprintf(sbuffer, "mkdir -p %s", output_directory);
  system(sbuffer);

  gen::generate(); // one call for NEVENTS

  // ROOT output disabled by default
  if (params::createRootOutput) {

    // Initialize ROOT output
    sprintf(sbuffer, "%s/%i.root", output_directory, prefix);
    TFile *outputFile = new TFile(sbuffer, "RECREATE");
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

  cout << "event generation done\n";
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
      strcpy(output_directory, argv[iarg + 1]);
      iarg += 2;
    } else if (strcmp(argv[iarg], "--surface") == 0 ||
               strcmp(argv[iarg], "-s") == 0) {
      strcpy(surface_file, argv[iarg + 1]);
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

// auxiliary function to get the number of lines
int getNlines(char *filename) {
  ifstream fin(filename);
  if (!fin) {
    cout << "getNlines funtion: Cannot open file " << filename << endl;
    exit(1);
  }
  string line;
  int nlines = 0;
  while (fin.good()) {
    getline(fin, line);
    nlines++;
  };
  fin.close();
  return nlines - 1;
}
