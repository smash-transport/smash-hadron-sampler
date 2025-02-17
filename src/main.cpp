#include <TFile.h>
#include <TRandom3.h>

#include <fstream>

#include "gen.h"
#include "oscaroutput.h"
#include "params.h"
#include "tree.h"
#include "build_metadata.h"

using namespace std ;
int getNlines(char *filename) ;
int readCommandLine(int argc, char** argv) ;

using params::sSpectraDir ;
using params::sSurface ;
using params::NEVENTS ;

int ranseed ;

extern "C"{
  void getranseedcpp_(int *seed)
  {
    *seed = ranseed ;
  }
}

// ########## MAIN block ##################

int main(int argc, char **argv)
{
  // command-line parameters
  int prefix = readCommandLine(argc, argv) ;
  params::printParameters() ;
  time_t time0 ;
  time(&time0) ;
  ranseed = time0+prefix*16 ;

  TRandom3* random3 = new TRandom3();
	random3->SetSeed(ranseed);
  cout<<"Random seed = "<<ranseed<<endl ;
  gen::rnd = random3 ;

 // ========== generator init
 gen::load(sSurface,getNlines(sSurface)) ;

 // ========== trees & files
 time_t start, end ;
 time(&start);

//============= main task
 char sbuffer [255] ;
 sprintf(sbuffer,"mkdir -p %s",sSpectraDir) ;
 system(sbuffer) ;

 // Initialize ROOT output
 sprintf(sbuffer, "%s/%i.root",sSpectraDir,prefix) ;
 TFile *outputFile = new TFile(sbuffer, "RECREATE");
 outputFile->cd();
 MyTree *treeIni = new MyTree(static_cast<const char*>("treeini")) ;

 gen::generate() ; // one call for NEVENTS

 // Write ROOT output
 for(int iev=0; iev<NEVENTS; iev++){
 treeIni->fill(iev) ;
 } // end events loop
 outputFile->Write() ;
 outputFile->Close() ;

 // Write Oscar output
 write_oscar_output();

 cout << "event generation done\n" ;
 time(&end); float diff2 = difftime(end, start);
 cout<<"Execution time = "<<diff2<< " [sec]" << endl;
 return 0;
}


int readCommandLine(int argc, char** argv)
{
	if(argc==1){cout << "NO PARAMETERS, exit" << endl ; exit(1) ;}
	int prefix = 0 ;
	if(strcmp(argv[1],"events")==0){
	  prefix = atoi(argv[2]) ;
	  cout << "events mode, prefix = " << prefix << endl ;
	  params::readParams(argv[3]) ;
  }else if(strcmp(argv[1],"fmax")==0){
	  if(static_cast<int>(argv[2][0]<58)){
		prefix = atoi(argv[2]) ;
		cout << "fmax mode, prefix = " << prefix << endl ;
		params::readParams(argv[3]) ;
	  }else
	  params::readParams(argv[2]) ;
  } else if (strcmp(argv[1], "--version")==0) {
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
	}else{cout << "unknown command-line switch: " << argv[1] << endl ; exit(1) ;}
	return prefix ;
}


// auxiliary function to get the number of lines
int getNlines(char *filename)
{
  ifstream fin(filename) ;
  if(!fin) {cout<<"getNlines: cannot open file: "<<filename<<endl; exit(1) ; }
  string line ;
  int nlines = 0 ;
  while(fin.good()){
    getline(fin,line) ; nlines++ ;
  } ;
  fin.close() ;
  return nlines-1 ;
}
