#include <TFile.h>
#include <TRandom3.h>

#include <fstream>

#include "gen.h"
#include "oscaroutput.h"
#include "params.h"
#include "tree.h"

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
	int prefix = 0;
	  prefix = atoi(argv[2]) ;
	for(int iarg=1; iarg<argc-1; iarg++){
		/*if(strcmp(argv[iarg],"-Npart")==0) {
			float Nparticipants_float = atof(argv[iarg+1]);
			Nparticipants = (int)round(Nparticipants_float);
		}*/
		if(strcmp(argv[iarg],"--parameters")==0) params::readParams(argv[iarg+1]);
		if(strcmp(argv[iarg],"--surface")==0) strcpy(sSurface, argv[iarg+1]);
		if(strcmp(argv[iarg],"--output")==0) strcpy(sSpectraDir, argv[iarg+1]);
	}
	cout << "hadronSampler: command line parameters are:\n";
	//cout << "Npart  " << Nparticipants << endl;
	cout << "hydro surface:  " << sSurface << endl;
	cout << "output dir:  " << sSpectraDir << endl;
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
