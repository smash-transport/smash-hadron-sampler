#include <omp.h>
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TH1.h>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <TF1.h>
#include <sstream>

#include "DatabasePDG2.h"
#include "gen.h"
#include "cascade.h"
#include "tree.h"
#include "params.h"

using namespace std ;
int getNlines(char *filename) ;
int readCommandLine(int argc, char** argv) ;

using params::sSpectraDir ;
using params::sMultDir ;
using params::sSurface ;
using params::NEVENTS ;
using params::bEventGeneration ;

// ########## MAIN block ##################

int main(int argc, char **argv)
{
  // command-line parameters
  int prefix = readCommandLine(argc, argv) ;
  params::printParameters() ;
  
  TRandom3* random3 = new TRandom3();
	random3->SetSeed();
  gen::rnd = random3 ;
  
//========= particle database init
	DatabasePDG2 *database = new DatabasePDG2("Tb/ptl3.data","Tb/dky3.mar.data");
	database->LoadData();
//	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
//	database->SetWidthRange(0., 10.0);
	database->SortParticlesByMass() ;
	database->CorrectBranching() ;
	database->DumpData() ;
	cout << " pion index = " << database->GetPionIndex() << endl ;
  gen::database = database ;
  
// ========== generator init
 gen::load(sSurface,getNlines(sSurface)) ;
 //cout << "dfMax = " << gen::calcDFMax() << endl ;

 // ========== trees & files
 time_t start, end ;
 time(&start);

//============= main task
if(bEventGeneration){ // ---- generate events
 char sbuffer [255] ;
 sprintf(sbuffer,"mkdir -p %s",sSpectraDir) ;
 system(sbuffer) ;
 sprintf(sbuffer,"%s/all.mult",sMultDir) ;
 gen::loadDFMax(sbuffer,getNlines(sbuffer)) ;
 sprintf(sbuffer, "%s/%i.root",sSpectraDir,prefix) ;
 TFile *outputFile = new TFile(sbuffer, "RECREATE"); 
 outputFile->cd();
 MyTree *treeIni = new MyTree("treeini") ;
 MyTree *treeFin = new MyTree("treefin") ;
 for(int iev=0; iev<NEVENTS; iev++){
 gen::generate() ;
 if(iev<10) cout<<"event #"<<iev<<" done\n" ;
 treeIni->fill() ;
 gen::urqmd(iev) ;
 treeFin->fill() ;
 } // end events loop
 outputFile->Write() ;
 outputFile->Close() ;
}else{ // ---- calculate fMax
char sbuffer [255] ;
 sprintf(sbuffer,"mkdir -p %s",sMultDir) ;
 system(sbuffer) ;
 sprintf(sbuffer, "%s/%i.mult", sMultDir, prefix) ; // here prefix=id of particle
 gen::calcDFMax(prefix, sbuffer) ;
}
 


 cout << "event generation done\n" ;
 time(&end); float diff2 = difftime(end, start);
 cout<<"Execution time = "<<diff2<< " [sec]" << endl;
 return 0;
}


int readCommandLine(int argc, char** argv)
{
	if(argc==1){cout << "NO PARAMETERS, exit" << endl ; exit(1) ;}
	bEventGeneration = 0 ;
	int prefix = 0 ;
	if(strcmp(argv[1],"events")==0){
	  bEventGeneration = 1 ;
	  prefix = atoi(argv[2]) ;
	  cout << "events mode, prefix = " << prefix << endl ;
	  params::readParams(argv[3]) ;
    }else if(strcmp(argv[1],"fmax")==0){
	  if((int)argv[2][0]<58){
		prefix = atoi(argv[2]) ;
		cout << "fmax mode, prefix = " << prefix << endl ;
		params::readParams(argv[3]) ;
	  }else
	  params::readParams(argv[2]) ;
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
