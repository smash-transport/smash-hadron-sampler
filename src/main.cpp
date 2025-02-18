#include <TFile.h>
#include <TRandom3.h>

#include <fstream>
#include <chrono> // Include the chrono library

#include "gen.h"
#include "oscaroutput.h"
#include "params.h"
#include "tree.h"
#include "build_metadata.h"

using namespace std;
int getNlines(char *filename);
int readCommandLine(int argc, char** argv);

using params::sSpectraDir;
using params::sSurface;
using params::NEVENTS;

int ranseed;

extern "C" {
    void getranseedcpp_(int *seed) {
        *seed = ranseed;
    }
}

// ########## MAIN block ##################

int main(int argc, char **argv) {
    // command-line parameters
    int prefix = readCommandLine(argc, argv);
    params::printParameters();

    // Use std::chrono to get the current time for the seed
    auto time0 = std::chrono::system_clock::now();
    auto duration_since_epoch = time0.time_since_epoch();
    auto seconds_since_epoch = std::chrono::duration_cast<std::chrono::seconds>(duration_since_epoch);
    ranseed = seconds_since_epoch.count() + prefix * 16;

    TRandom3* random3 = new TRandom3();
    random3->SetSeed(ranseed);
    cout << "Random seed = " << ranseed << endl;
    gen::rnd = random3;

    // ========== generator init
    gen::load(sSurface, getNlines(sSurface));

    // ========== trees & files
    auto start = std::chrono::high_resolution_clock::now(); // Start timing

    // ============= main task
    char sbuffer[255];
    sprintf(sbuffer, "mkdir -p %s", sSpectraDir);
    system(sbuffer);

    // Initialize ROOT output
    sprintf(sbuffer, "%s/%i.root", sSpectraDir, prefix);
    TFile *outputFile = new TFile(sbuffer, "RECREATE");
    outputFile->cd();
    MyTree *treeIni = new MyTree(static_cast<const char*>("treeini"));

    gen::generate(); // one call for NEVENTS

    // Write ROOT output
    for (int iev = 0; iev < NEVENTS; iev++) {
        treeIni->fill(iev);
    } // end events loop
    outputFile->Write();
    outputFile->Close();

    // Write Oscar output
    write_oscar_output();

    // End timing
    auto end = std::chrono::high_resolution_clock::now();
    auto execution_time = std::chrono::duration_cast<std::chrono::duration<float>>(end - start);

    cout << "Event generation done\n";
    cout << "Execution time = " << execution_time.count() << " [sec]" << endl;

    return 0;
}

int readCommandLine(int argc, char** argv) {
    if (argc == 1) { cout << "NO PARAMETERS, exit" << endl; exit(1); }
    int prefix = 0;
    if (strcmp(argv[1], "events") == 0) {
        prefix = atoi(argv[2]);
        cout << "events mode, prefix = " << prefix << endl;
        params::readParams(argv[3]);
    } else if (strcmp(argv[1], "fmax") == 0) {
        if (static_cast<int>(argv[2][0] < 58)) {
            prefix = atoi(argv[2]);
            cout << "fmax mode, prefix = " << prefix << endl;
            params::readParams(argv[3]);
        } else
            params::readParams(argv[2]);
    } else if (strcmp(argv[1], "--version") == 0) {
        std::printf(
            "%s\n"
#ifdef GIT_BRANCH
            "Branch   : %s\n"
#endif
            "System   : %s\nCompiler : %s %s\n"
            "Date     : %s\n",
            SAMPLER_VERSION,
#ifdef GIT_BRANCH
            GIT_BRANCH,
#endif
            CMAKE_SYSTEM, CMAKE_CXX_COMPILER_ID, CMAKE_CXX_COMPILER_VERSION,
            BUILD_DATE);
        std::exit(EXIT_SUCCESS);
    } else { cout << "unknown command-line switch: " << argv[1] << endl; exit(1); }
    return prefix;
}

// auxiliary function to get the number of lines
int getNlines(char *filename) {
    ifstream fin(filename);
    if (!fin) { cout << "getNlines: cannot open file: " << filename << endl; exit(1); }
    string line;
    int nlines = 0;
    while (fin.good()) {
        getline(fin, line); nlines++;
    };
    fin.close();
    return nlines - 1;
}