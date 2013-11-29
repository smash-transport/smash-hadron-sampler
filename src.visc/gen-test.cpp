#include <omp.h>
#include <TGraph.h>
#include <TMath.h>
#include <TGraph.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TF1.h>

#include "DatabasePDG2.h"
#include "cascade.h"
#include "params.h"
#include "const.h"

using namespace std ;

#define _USE_MATH_DEFINES
const double C_Feq = (pow(0.5/M_PI/hbarC,3)) ;

// ##########################################################
// #  this version works with arbitrary T/mu distribution   #
// #  on freezeout hypersurface (May'2012)                  #
// #  also, pre-BOOSTED dsigma is used                      #
// ##########################################################


// active Lorentz boost
void fillBoostMatrix(double vx, double vy, double vz, double boostMatrix [4][4])
// here in boostMatrix [0]=t, [1]=x, [2]=y, [3]=z
{
  const double vv [3] = {vx, vy, vz} ;
  const double v2 = vx*vx+vy*vy+vz*vz ;
  const double gamma = 1.0/sqrt(1.0-v2) ;
  if(isinf(gamma)||isnan(gamma)){ cout<<"boost vector invalid; exiting\n" ; exit(1) ; }
  boostMatrix[0][0] = gamma ;
  boostMatrix[0][1] = boostMatrix[1][0] = vx*gamma ;
  boostMatrix[0][2] = boostMatrix[2][0] = vy*gamma ;
  boostMatrix[0][3] = boostMatrix[3][0] = vz*gamma ;
  if(v2>0.0){
  for(int i=1; i<4; i++)
  for(int j=1; j<4; j++)
   boostMatrix[i][j] = (gamma-1.0)*vv[i-1]*vv[j-1]/v2 ;
  }else{
  for(int i=1; i<4; i++)
  for(int j=1; j<4; j++)
   boostMatrix[i][j] = 0.0 ;
  }
  for(int i=1; i<4; i++) boostMatrix[i][i] += 1.0 ;
}


// index44: returns an index of pi^{mu nu} mu,nu component in a plain 1D array
int index44(const int &i, const int &j){
  if(i>3 || j>3 || i<0 || j<0) {std::cout<<"index44: i j " <<i<<" "<<j<<endl ; exit(1) ; }
  if(j<i) return (i*(i+1))/2 + j ;
  else return (j*(j+1))/2 + i ;
}


namespace gen{


int Nelem ;
double *ntherm, dvMax, dsigmaMax ;
TRandom3 *rnd ;
DatabasePDG2 *database ;
int NPART ;
const int npartbuf = 10000 ;

struct element {
 double tau, x, y, eta ;
 double u[4] ;
 double dsigma[4] ;
 double T, mub, muq, mus ;
 double pi[10] ;
 double Pi ;
} ;

element *surf ;
float **Px, **Py, **Pz, **E, **X, **Y, **Z, **T ; // particle arrays
int **Acc, **Id, **MId ;
char **Charge ;
int *npart ;               // number of particles in arrays Px...E, X...T

const double pmax = 10.0 ;
const double rapmax = 6.0 ;

const double c1 = pow(1./2./hbarC/TMath::Pi(),3.0) ;
double *cumulantDensity ; // particle densities (thermal). Seems to be redundant, but needed for fast generation
double totalDensity ; // sum of all thermal densities


// ######## load the elements
void load(char *filename, int N)
{
 double dV, vEff=0.0, vEffOld=0.0, dvEff, dvEffOld ;
 int nfail=0 ;
 TLorentzVector dsigma ;
 Nelem = N ;
 surf = new element [Nelem] ;

 X  = new Float_t* [params::NEVENTS] ;
 Y  = new Float_t* [params::NEVENTS] ;
 Z  = new Float_t* [params::NEVENTS] ;
 T  = new Float_t* [params::NEVENTS] ;
 Px = new Float_t* [params::NEVENTS] ;
 Py = new Float_t* [params::NEVENTS] ;
 Pz = new Float_t* [params::NEVENTS] ;
 E  = new Float_t* [params::NEVENTS] ;
 Acc  = new Int_t* [params::NEVENTS] ;
 Id   = new Int_t* [params::NEVENTS] ;
 MId  = new Int_t* [params::NEVENTS] ;
 Charge = new Char_t* [params::NEVENTS] ;
 for(int i=0; i<params::NEVENTS; i++){
   X[i]  = new Float_t [npartbuf] ;
   Y[i]  = new Float_t [npartbuf] ;
   Z[i]  = new Float_t [npartbuf] ;
   T[i]  = new Float_t [npartbuf] ;
   Px[i] = new Float_t [npartbuf] ;
   Py[i] = new Float_t [npartbuf] ;
   Pz[i] = new Float_t [npartbuf] ;
   E[i]  = new Float_t [npartbuf] ;
   Acc[i]  = new Int_t [npartbuf] ;
   Id[i]   = new Int_t [npartbuf] ;
   MId[i]  = new Int_t [npartbuf] ;
   Charge[i]  = new Char_t [npartbuf] ;
 }
 npart = new Int_t [params::NEVENTS] ;

 cout<<"reading "<<N<<" lines from  "<<filename<<"\n" ;
 ifstream fin(filename) ;
 if(!fin){ cout << "cannot read file " << filename << endl ; exit(1) ; }
 dvMax=0. ;
 dsigmaMax=0. ;
 // ---- reading loop
 string line ;
 istringstream instream ;
 cout<<"1?: failbit="<<instream.fail()<<endl ;
 for(int n=0; n<Nelem; n++){
   getline(fin, line) ;
   instream.str(line) ;
   instream.seekg(0) ;
   if(params::shear){ // read all including viscous terms
    instream>>surf[n].tau>>surf[n].x>>surf[n].y>>surf[n].eta
      >>surf[n].dsigma[0]>>surf[n].dsigma[1]>>surf[n].dsigma[2]>>surf[n].dsigma[3]
      >>surf[n].u[0]>>surf[n].u[1]>>surf[n].u[2]>>surf[n].u[3]
      >>surf[n].T>>surf[n].mub>>surf[n].muq>>surf[n].mus ;
      for(int i=0; i<10; i++) instream>>surf[n].pi[i] ;
      instream>>surf[n].Pi ;
   }else{ // read only ideal part
    instream>>surf[n].tau>>surf[n].x>>surf[n].y>>surf[n].eta
      >>surf[n].dsigma[0]>>surf[n].dsigma[1]>>surf[n].dsigma[2]>>surf[n].dsigma[3]
      >>surf[n].u[0]>>surf[n].u[1]>>surf[n].u[2]>>surf[n].u[3]
      >>surf[n].T>>surf[n].mub>>surf[n].muq>>surf[n].mus>>dV ;
   }
   if(instream.fail()){ cout<<"reading failed at line "<<n<<"; exiting\n" ; exit(1) ; }
   // calculate in the old way
   dvEffOld = surf[n].dsigma[0]*surf[n].u[0]+surf[n].dsigma[1]*surf[n].u[1]+
   surf[n].dsigma[2]*surf[n].u[2]+surf[n].dsigma[3]*surf[n].u[3] ;
   vEffOld += dvEffOld ;
   if(dvEffOld<0.0){
     //cout<<"!!! dvOld!=dV " << dvEffOld <<"  " << dV << "  " << surf[n].tau <<endl ;
     nfail++ ;
   }
   //if(nfail==100) exit(1) ;
   // ---- boost
   dsigma.SetXYZT(-surf[n].dsigma[1],-surf[n].dsigma[2],-surf[n].dsigma[3],surf[n].dsigma[0]) ;
   dsigma.Boost(-surf[n].u[1]/surf[n].u[0],-surf[n].u[2]/surf[n].u[0],-surf[n].u[3]/surf[n].u[0]) ;
   // ######################################################################
   // ###  now and later surf.dsigma are boosted to the fluid rest frame  ##
   // ######################################################################
   surf[n].dsigma[0] = dsigma.T() ;
   surf[n].dsigma[1] = -dsigma.X() ;
   surf[n].dsigma[2] = -dsigma.Y() ;
   surf[n].dsigma[3] = -dsigma.Z() ;
   dvEff = surf[n].dsigma[0] ;
   vEff += dvEff ;
   if(dvMax<dvEff) dvMax = dvEff ;
   // maximal value of the weight max(W) = max(dsigma_0+|\vec dsigma_i|)   for equilibrium DFs
   if(dsigma.T()+dsigma.Rho()>dsigmaMax) dsigmaMax = dsigma.T()+dsigma.Rho() ;
   // ########################
   // pi^{mu nu} boost to fluid rest frame
   // ########################
   if(params::shear){
   double _pi[10], boostMatrix[4][4] ;
   fillBoostMatrix(-surf[n].u[1]/surf[n].u[0],-surf[n].u[2]/surf[n].u[0],-surf[n].u[3]/surf[n].u[0], boostMatrix) ;
   for(int i=0; i<4; i++)
   for(int j=i; j<4; j++){
     _pi[index44(i,j)] = 0.0 ;
     for(int k=0; k<4; k++)
     for(int l=0; l<4; l++)
      _pi[index44(i,j)] += surf[n].pi[index44(k,l)]*boostMatrix[i][k]*boostMatrix[j][l] ;
   }
   for(int i=0; i<10; i++) surf[n].pi[i] = _pi[i] ;
   } // end pi boost
 }
 if(params::shear) dsigmaMax *= 1.3*1.5 ;
 else dsigmaMax *= 1.3 ;

 cout<<"..done.\n" ;
 cout<<"Veff = "<<vEff<<"  dvMax = "<<dvMax<<endl ;
 cout<<"Veff(old) = "<<vEffOld<<endl ;
 cout<<"failed elements = "<<nfail<<endl ;
// ---- prepare some stuff to calculate thermal densities
 NPART=database->GetNParticles() ;
 cout<<"NPART="<<NPART<<endl ;
 cout<<"dsigmaMax="<<dsigmaMax<<endl ;
 cumulantDensity = new double [NPART] ;
}


void loadDFMax(char *filename, int N)
{
  ifstream fin(filename) ;
  int pid ;
  double fmax ;
  if(!fin){ cout<<"readFMax: cannot open: "<<filename<<endl; exit(1); } 
  for(int i=0; i<N; i++){
    fin>>pid>>fmax ;
    ParticlePDG2 *part = database->GetPDGParticle(pid) ;
    if(part!=0) part->SetFMax(fmax) ;
    else{ cout<<"loadFMax: unrecognized particle id = "<<pid<<endl ; exit(1) ; }
  }
}


// ######## Monte-Carlo
double calcDFMax(int pindex, char *fileout)
{
 const double gmumu [4] = {1., -1., -1., -1.} ;
 ParticlePDG2 *part = database->GetPDGParticleByIndex(pindex) ;
 const double mass = part->GetMass() ;
 const double J = part->GetSpin() ;
 const double stat = int(2.*J) & 1 ? -1. : 1. ;
 double dfMax = 0. ;
 for(int iel=0; iel<Nelem; iel++){ // loop over all elements
   for(int imc=0; imc<1000; imc++){ // momentum generation loop
   double p =  10.0*rnd->Rndm() ;
   double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   double sinth = -1.0 + 2.0*rnd->Rndm() ;
   const double px = p*sqrt(1.0-sinth*sinth)*cos(phi) ;
   const double py = p*sqrt(1.0-sinth*sinth)*sin(phi) ;
   const double pz = p*sinth ;
   const double E = sqrt(p*p+mass*mass) ;
   double W = ( surf[iel].dsigma[0]*E + surf[iel].dsigma[1]*px +
        surf[iel].dsigma[2]*py + surf[iel].dsigma[3]*pz ) / E ; 
   if(params::shear && p<2.0){ // cut viscous corrections for p>2 GeV
    const double muf = part->GetBaryonNumber()*surf[iel].mub + part->GetStrangeness()*surf[iel].mus +
               part->GetElectricCharge()*surf[iel].muq ;
    const double feq = C_Feq/( exp((E-muf)/surf[iel].T) - stat ) ;
    double pipp = 0 ;
    double momArray [4] = {E,px,py,pz} ;  // TLorentzVector use different layout of components
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
     pipp += momArray[i]*momArray[j]*gmumu[i]*gmumu[j]*surf[iel].pi[index44(i,j)] ;
    W *= (1.0 + (1.0+stat*feq)*pipp/(2.*surf[iel].T*surf[iel].T*(params::ecrit*1.15))) ;
   }
   if(W>dfMax){
     dfMax = W ;
   }
   } // end momentum generation
 } // loop over all elements
 dfMax *= 1.3 ;
 //part->SetFMax(dfMax) ;
 ofstream fout(fileout) ;
 if(!fout){ cout<<"readFMax: cannot open(for writing): "<<fileout<<endl; exit(1); }
 fout<<setw(14)<<part->GetPDG()<<setw(14)<<dfMax<<endl ;
 fout.close() ;
 return dfMax ;
}


void acceptParticle(int event, int lid, double lx, double ly, double lz, double lt, double lpx, double lpy, double lpz, double lE, int charge) ;


double ffthermal(double *x, double *par)
{
  double &T = par[0] ;
  double &mu = par[1] ;
  double &mass = par[2] ;
  double &stat = par[3] ;
  return x[0]*x[0]/( exp((sqrt(x[0]*x[0]+mass*mass)-mu)/T) - stat ) ;
}

int generate()
{
 const double gmumu [4] = {1., -1., -1., -1.} ;
 TF1 *fthermal = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom ;
 for(int iev=0; iev<params::NEVENTS; iev++) npart[iev] = 0 ;
 int nmaxiter = 0 ;
 double n_pi = 0.0 ;

 for(int iel=0; iel<Nelem; iel++){ // loop over all elements
  // ---> thermal densities, for each surface element
   totalDensity = 0.0 ;
   for(int ip=328; ip<NPART; ip++){
    double density = 0. ;
    ParticlePDG2 *particle = database->GetPDGParticleByIndex(ip) ;
    const double mass = particle->GetMass() ;
    const double J = particle->GetSpin() ;
    const double stat = int(2.*J) & 1 ? -1. : 1. ;
    const double muf = particle->GetBaryonNumber()*surf[iel].mub + particle->GetStrangeness()*surf[iel].mus +
               particle->GetElectricCharge()*surf[iel].muq ;
    for(int i=1; i<11; i++)
    density += (2.*J+1.)*pow(gevtofm,3)/(2.*pow(TMath::Pi(),2))*mass*mass*surf[iel].T*pow(stat,i+1)*TMath::BesselK(2,i*mass/surf[iel].T)*exp(i*muf/surf[iel].T)/i ;
    if(ip>0) cumulantDensity[ip] = cumulantDensity[ip-1] + density ;
        else cumulantDensity[ip] = density ;
    totalDensity += density ;
   }
   //cout<<"thermal densities calculated.\n" ;
   //cout<<cumulantDensity[NPART-1]<<" = "<<totalDensity<<endl ;
 // ---< end thermal densities calc
  double rval, dvEff = 0., W ;
  // dvEff = dsigma_mu * u^mu
  dvEff = surf[iel].dsigma[0] ;
  n_pi += dvEff*totalDensity ;
  for(int ievent=0; ievent<params::NEVENTS; ievent++){
  // ---- number of particles to generate
  int nToGen = 0 ;
  if(dvEff*totalDensity<0.01){
    double x = rnd->Rndm() ; // throw dice
    if(x<dvEff*totalDensity) nToGen = 1 ;
  }else{
    nToGen = rnd->Poisson(dvEff*totalDensity) ;
  }
   // ---- we generate a particle!
  for(int ip=0; ip<nToGen; ip++){
  int isort = 0 ;
  double xsort = rnd->Rndm()*totalDensity ; // throw dice, particle sort
  while(cumulantDensity[isort]<xsort) isort++ ;
   ParticlePDG2 *part = database->GetPDGParticleByIndex(isort) ;
   const double J = part->GetSpin() ;
   const double mass = part->GetMass() ;
   const double stat = int(2.*J) & 1 ? -1. : 1. ;
   const double muf = part->GetBaryonNumber()*surf[iel].mub + part->GetStrangeness()*surf[iel].mus +
               part->GetElectricCharge()*surf[iel].muq ;
   fthermal->SetParameters(surf[iel].T,muf,mass,stat) ;
   //const double dfMax = part->GetFMax() ;
   int niter = 0 ; // number of iterations, for debug purposes
   do{ // fast momentum generation loop
   const double p = fthermal->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   mom.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi), p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   W = ( surf[iel].dsigma[0]*mom.E() + surf[iel].dsigma[1]*mom.Px() +
        surf[iel].dsigma[2]*mom.Py() + surf[iel].dsigma[3]*mom.Pz() ) / mom.E() ;
   double WviscFactor = 1.0 ;
   if(params::shear){ // cut viscous corrections for p>2 GeV
    const double feq = C_Feq/( exp((sqrt(p*p+mass*mass)-muf)/surf[iel].T) - stat ) ;
    double pipp = 0 ;
    double momArray [4] = {mom[3],mom[0],mom[1],mom[2]} ;
    for(int i=0; i<4; i++)
    for(int j=0; j<4; j++)
     pipp += momArray[i]*momArray[j]*gmumu[i]*gmumu[j]*surf[iel].pi[index44(i,j)] ;
    WviscFactor = (1.0 + (1.0+stat*feq)*pipp/(2.*surf[iel].T*surf[iel].T*(params::ecrit*1.15))) ;
    if(WviscFactor<0.5) WviscFactor = 0.5 ;
    if(WviscFactor>1.5) WviscFactor = 1.5 ;
   }
   W *= WviscFactor ;
   rval = rnd->Rndm()*dsigmaMax ;
   niter++ ;
   }while(rval>W) ; // end fast momentum generation
   if(niter>nmaxiter) nmaxiter = niter ;
   // additional random smearing over eta
   const double etaF = 0.5*log((surf[iel].u[0]+surf[iel].u[3])/(surf[iel].u[0]-surf[iel].u[3])) ;
   const double etaShift = params::deta*(-0.5+rnd->Rndm()) ;
   const double vx = surf[iel].u[1]/surf[iel].u[0]*cosh(etaF)/cosh(etaF+etaShift) ;
   const double vy = surf[iel].u[2]/surf[iel].u[0]*cosh(etaF)/cosh(etaF+etaShift) ;
   const double vz = tanh(etaF+etaShift) ;
   mom.Boost(vx,vy,vz) ;
   acceptParticle(ievent,part->GetPDG(), surf[iel].x, surf[iel].y, surf[iel].tau*sinh(surf[iel].eta+etaShift),
     surf[iel].tau*cosh(surf[iel].eta+etaShift), mom.Px(), mom.Py(), mom.Pz(), mom.E(), part->GetElectricCharge()) ;
  } // coordinate accepted
  } // events loop
  if(iel%(Nelem/50)==0) cout<<(iel*100)/Nelem<<" percents done, maxiter= "<<nmaxiter<<endl ;
 } // loop over all elements
 delete fthermal ;
 cout << "n_pi = "<<n_pi<<endl ;
 return npart[0] ;
}



void acceptParticle(int ievent, int lid, double lx, double ly, double lz, double lt, double lpx, double lpy, double lpz, double lE, int charge)
{
  int& npart1 = npart[ievent] ;
 int urqmdid, urqmdiso3 ;
 pdg2id_(&urqmdid, &urqmdiso3, &lid) ;
 if(geteposcode_(&lid)!=0 && abs(urqmdid)<1000){  // particle known to UrQMD
// if(true){ // TEST!! for thermal mult's
    Id[ievent][npart1] = lid ;
    MId[ievent][npart1] = 0 ;
    Charge[ievent][npart1] = charge ;
    X[ievent][npart1] = lx ;
    Y[ievent][npart1] = ly ;
    Z[ievent][npart1] = lz ;
    T[ievent][npart1] = lt ;
    Px[ievent][npart1] = lpx ;
    Py[ievent][npart1] = lpy ;
    Pz[ievent][npart1] = lpz ;
    E[ievent][npart1] = lE ;
#ifdef DEBUG2
    cout << "accepted" << setw(5) << lid << setw(14) << x[jtotpart] << setw(14) << y[jtotpart] << setw(14) << z[jtotpart] << endl ;
#endif
    npart1++ ;
 }else{ // decay particles unknown to UrQMD
  double mom [4] = {lpx,lpy,lpz,lE} ;
#ifdef DEBUG2
  cout << "------ unstable particle decay (Cooper-Frye isotherm) " << lid << endl ;
  cout << setw(14) << "px" << setw(14) << "py" << setw(14) << "pz" << setw(14) << "E" << endl ;
  cout << setw(14) << mom[0] << setw(14) << mom[1] << setw(14) << mom[2] << setw(14) << mom[3] << endl ;
#endif
  int nprod, ppid [4], pchrg [4] ;
  double pmom [3][4] ;
  decay(lid, mom, nprod, ppid, pmom, pchrg) ;
  for(int iprod=0; iprod<nprod; iprod++){ // add decay products to UrQMD input
    pdg2id_(&urqmdid, &urqmdiso3, &ppid[iprod]) ;
    if(geteposcode_(&ppid[iprod])!=0 && abs(urqmdid)<1000){  // particle known to UrQMD
    Id[ievent][npart1] = ppid[iprod] ;
    Charge[ievent][npart1] = pchrg[iprod] ;
    MId[ievent][npart1] = lid ;
    X[ievent][npart1] = lx ;
    Y[ievent][npart1] = ly ;
    Z[ievent][npart1] = lz ;
    T[ievent][npart1] = lt ;
    Px[ievent][npart1] = pmom[iprod][0] ;
    Py[ievent][npart1] = pmom[iprod][1] ;
    Pz[ievent][npart1] = pmom[iprod][2] ;
    E[ievent][npart1] = pmom[iprod][3] ;
    npart1++ ; 
  }}
 }
 if(npart1>npartbuf){ cout<<"Error. Please increase gen::npartbuf\n"; exit(1);}
}

// ################### end #################
} // end namespace gen
