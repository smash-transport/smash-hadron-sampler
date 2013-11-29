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
using params::Temp ;
using params::mu_b ;
using params::mu_q ;
using params::mu_s ;
using params::weakContribution ;
using params::rescatter ;

// #######################
// ###   BOOSTED dsigma  #
// #######################


// active Lorentz boost
void fillBoostMatrix(double vx, double vy, double vz, double boostMatrix [4][4])
// here in boostMatrix [0]=x, [1]=y, [2]=z, [3]=t, which is the same notation as used in TLorentzVector
{
  const double vv [3] = {vx, vy, vz} ;
  const double v2 = vx*vx+vy*vy+vz*vz ;
  const double gamma = 1.0/sqrt(1.0-v2) ;
  if(isinf(gamma)||isnan(gamma)){ cout<<"boost vector invalid; exiting\n" ; exit(1) ; }
  boostMatrix[3][3] = gamma ;
  boostMatrix[3][0] = boostMatrix[0][3] = vx*gamma ;
  boostMatrix[3][1] = boostMatrix[1][3] = vy*gamma ;
  boostMatrix[3][2] = boostMatrix[2][3] = vz*gamma ;
  if(v2>0.0){
  for(int i=0; i<3; i++)
  for(int j=0; j<3; j++)
   boostMatrix[i][j] = (gamma-1.0)*vv[i]*vv[j]/v2 ;
  }else{
  for(int i=0; i<3; i++)
  for(int j=0; j<3; j++)
   boostMatrix[i][j] = 0.0 ;
  }
  for(int i=0; i<3; i++) boostMatrix[i][i] += 1.0 ;
}


// index44: returns an index of pi^{mu nu} mu,nu component in a plain 1D array
int index44(const int &i, const int &j){
  if(i>3 || j>3 || i<0 || j<0) {std::cout<<"index44: i j " <<i<<" "<<j<<endl ; exit(1) ; }
  if(j<i) return (i*(i+1))/2 + j ;
  else return (j*(j+1))/2 + i ;
}


namespace gen{


int Nelem ;
double *ntherm, dvMax ;
TRandom3 *rnd ;
DatabasePDG2 *database ;
int NPART ;

struct element {
 double tau, x, y, eta ;
 double u[4] ;
 double dsigma[4] ;
 double pi[10] ;
} ;

element *surf ;
Float_t *Px, *Py, *Pz, *E, *X, *Y, *Z, *T ; // particle arrays
Int_t *Acc, *Id, *MId ;
Int_t npart ;               // number of particles in arrays Px...E, X...T

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
 X  = new Float_t [10000] ;
 Y  = new Float_t [10000] ;
 Z  = new Float_t [10000] ;
 T  = new Float_t [10000] ;
 Px = new Float_t [10000] ;
 Py = new Float_t [10000] ;
 Pz = new Float_t [10000] ;
 E  = new Float_t [10000] ;
 Acc  = new Int_t [10000] ;
 Id   = new Int_t [10000] ;
 MId  = new Int_t [10000] ;

 cout<<"reading "<<N<<" lines from  "<<filename<<"\n" ;
 ifstream fin(filename) ;
 if(!fin){ cout << "cannot read file " << filename << endl ; exit(1) ; }
 dvMax=0. ;
 // ---- reading loop
 string line ;
 istringstream instream ;
 cout<<"1?: failbit="<<instream.fail()<<endl ;
 for(int n=0; n<Nelem; n++){
   getline(fin, line) ;
   instream.str(line) ;
   instream.seekg(0) ;
    instream>>surf[n].tau>>surf[n].x>>surf[n].y>>surf[n].eta
   >>surf[n].dsigma[0]>>surf[n].dsigma[1]>>surf[n].dsigma[2]>>surf[n].dsigma[3]>>surf[n].u[0]>>surf[n].u[1]>>surf[n].u[2]>>surf[n].u[3]>>dV ;
   if(instream.fail()){ cout<<"reading failed at line "<<n<<"; exiting\n" ; exit(1) ; }
   // calculate in the old way
   dvEffOld = surf[n].dsigma[0]*surf[n].u[0]+surf[n].dsigma[1]*surf[n].u[1]+
   surf[n].dsigma[2]*surf[n].u[2]+surf[n].dsigma[3]*surf[n].u[3] ;
   vEffOld += dvEffOld ;
   if(dV<0.0){
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
 }
 cout<<"..done.\n" ;
 cout<<"Veff = "<<vEff<<"  dvMax = "<<dvMax<<endl ;
 cout<<"Veff(old) = "<<vEffOld<<endl ;
//========= thermal densities
 NPART=database->GetNParticles() ;
 cout<<"NPART="<<NPART<<endl ;
 cumulantDensity = new double [NPART] ;
 totalDensity = 0.0 ;
 for(int ip=0; ip<NPART; ip++){
  double density = 0. ;
  ParticlePDG2 *particle = database->GetPDGParticleByIndex(ip) ;
  const double mass = particle->GetMass() ;
  const double J = particle->GetSpin() ;
  const double stat = int(2.*J) & 1 ? -1. : 1. ;
  const double muf = particle->GetBaryonNumber()*mu_b + particle->GetStrangeness()*mu_s +
               particle->GetElectricCharge()*mu_q ;
  for(int i=1; i<11; i++)
  density += (2.*J+1.)*pow(gevtofm,3)/(2.*pow(TMath::Pi(),2))*mass*mass*Temp*pow(stat,i+1)*TMath::BesselK(2,i*mass/Temp)*exp(i*muf/Temp)/i ;
  if(ip>0) cumulantDensity[ip] = cumulantDensity[ip-1] + density ;
      else cumulantDensity[ip] = density ;
  totalDensity += density ;
  particle->SetDensity(density) ;
 }
 cout<<"thermal densities calculated.\n" ;
 cout<<cumulantDensity[NPART-1]<<" = "<<totalDensity<<endl ;
 cout<<"failed elements = "<<nfail<<endl ;
 exit(1) ;
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
 ParticlePDG2 *part = database->GetPDGParticleByIndex(pindex) ;
 const double mass = part->GetMass() ;
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
   if(W>dfMax) dfMax = W ;
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


void acceptParticle(int lid, double lx, double ly, double lz, double lt, double lpx, double lpy, double lpz, double lE) ;


double ffthermal(double *x, double *par)
{
  double &T = par[0] ;
  double &mu = par[1] ;
  double &mass = par[2] ;
  double &stat = par[3] ;
  return x[0]*x[0]/( exp((sqrt(x[0]*x[0]+mass*mass)-mu)/T) - stat ) ;
}

int generate(void)
{
 TF1 *fthermal = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom ;
 npart = 0 ;
 for(int iel=0; iel<Nelem; iel++){ // loop over all elements
  double rval, dvEff = 0., W ;
  // dvEff = dsigma_mu * u^mu
  dvEff = surf[iel].dsigma[0] ;
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
   const double muf = part->GetBaryonNumber()*mu_b + part->GetStrangeness()*mu_s +
              part->GetElectricCharge()*mu_q ;
   fthermal->SetParameters(Temp,muf,mass,stat) ;
   const double dfMax = part->GetFMax() ;
   do{ // fast momentum generation loop
   const double p = fthermal->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   mom.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi), p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   W = ( surf[iel].dsigma[0]*mom.E() + surf[iel].dsigma[1]*mom.Px() +
        surf[iel].dsigma[2]*mom.Py() + surf[iel].dsigma[3]*mom.Pz() ) / mom.E() ; 
   rval = rnd->Rndm()*dfMax ;
   }while(rval>W) ; // end fast momentum generation
   // additional random smearing over eta
   const double etaF = 0.5*log((surf[iel].u[0]+surf[iel].u[3])/(surf[iel].u[0]-surf[iel].u[3])) ;
   const double etaShift = params::deta*(-0.5+rnd->Rndm()) ;
   const double vx = surf[iel].u[1]/surf[iel].u[0]*cosh(etaF)/cosh(etaF+etaShift) ;
   const double vy = surf[iel].u[2]/surf[iel].u[0]*cosh(etaF)/cosh(etaF+etaShift) ;
   const double vz = tanh(etaF+etaShift) ;
   mom.Boost(vx,vy,vz) ;
   acceptParticle(part->GetPDG(), surf[iel].x, surf[iel].y, surf[iel].tau*sinh(surf[iel].eta+etaShift),
     surf[iel].tau*cosh(surf[iel].eta+etaShift), mom.Px(), mom.Py(), mom.Pz(), mom.E()) ;
  } // coordinate accepted
 } // loop over all elements
 return npart ;
}



void acceptParticle(int lid, double lx, double ly, double lz, double lt, double lpx, double lpy, double lpz, double lE)
{
 int urqmdid, urqmdiso3 ;
 pdg2id_(&urqmdid, &urqmdiso3, &lid) ;
 if(geteposcode_(&lid)!=0 && abs(urqmdid)<1000){  // particle known in UrQMD
// if(true){ // TEST!! for thermal mult's
    Id[npart] = lid ;
    MId[npart] = 0 ;
    X[npart] = lx ;
    Y[npart] = ly ;
    Z[npart] = lz ;
    T[npart] = lt ;
    Px[npart] = lpx ;
    Py[npart] = lpy ;
    Pz[npart] = lpz ;
    E[npart] = lE ;
#ifdef DEBUG2
    cout << "accepted" << setw(5) << lid << setw(14) << x[jtotpart] << setw(14) << y[jtotpart] << setw(14) << z[jtotpart] << endl ;
#endif
    npart++ ;
 }else{ // decay particles unknown to UrQMD
  double mom [4] = {lpx,lpy,lpz,lE} ;
#ifdef DEBUG2
  cout << "------ unstable particle decay (Cooper-Frye isotherm) " << lid << endl ;
  cout << setw(14) << "px" << setw(14) << "py" << setw(14) << "pz" << setw(14) << "E" << endl ;
  cout << setw(14) << mom[0] << setw(14) << mom[1] << setw(14) << mom[2] << setw(14) << mom[3] << endl ;
#endif
  int nprod, ppid [4] ;
  double pmom [3][4] ;
  decay(lid, mom, nprod, ppid, pmom) ;
  for(int iprod=0; iprod<nprod; iprod++){ // add decay products to UrQMD input
    pdg2id_(&urqmdid, &urqmdiso3, &ppid[iprod]) ;
    if(geteposcode_(&ppid[iprod])!=0 && abs(urqmdid)<1000){  // particle known in UrQMD
    Id[npart] = ppid[iprod] ;
    MId[npart] = lid ;
    X[npart] = lx ;
    Y[npart] = ly ;
    Z[npart] = lz ;
    T[npart] = lt ;
    Px[npart] = pmom[iprod][0] ;
    Py[npart] = pmom[iprod][1] ;
    Pz[npart] = pmom[iprod][2] ;
    E[npart] = pmom[iprod][3] ;
    npart++ ; 
  }}
 }
}

// ################### end #################
} // end namespace gen
