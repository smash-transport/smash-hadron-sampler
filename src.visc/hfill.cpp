#include <TString.h>
#include <TLorentzVector.h>
#include "hfill.h"
#include "params.h"
#include "const.h"

//#######################################
//##   this module is NOT used so far  ##
//#######################################
using params::NEVENTS ;
int NPART ;

//----------- constants for correlation analysis
 const int PID = -211 ;  // particle id for correlations (default: pi-, -211)
 const int NPT = 6 ;
 char* binnames [NPT] = {"01-02","02-03","03-04","04-05","055-075","075-095"} ;
 double pt_low [NPT] = {0.1, 0.2, 0.3, 0.4, 0.55, 0.75} ;
 double pt_high [NPT] = {0.2, 0.3, 0.4, 0.5, 0.75, 0.95} ;
 extern double QMAX ;
 extern int NBINS ;
//------- other constants
 int pid_pt [] = {-211,211, -321,321, 2212, -2212, 311, 3122, -3122, /* XiOm */ 3312, -3312, 3322, -3322, 3334, -3334} ;
 const int npid_pt = sizeof(pid_pt)/sizeof(int) ;
 int pid_v2 [] = {211, 321, 2212} ;
 const int npid_v2 = sizeof(pid_v2)/sizeof(int) ;
 char* name_pt [npid_pt] = {"ptpion","ptpion+","ptkaon","ptkaon+","ptproton","ptaproton","ptkaon0","ptLambda","ptaLambda",
	"ptXi-","ptaXi-","ptXi0","ptaXi0","ptOmega","ptaOmega"} ;
 char* name_v2 [npid_v2] = {"pions","kaons","protons"} ;
 double kPtBins [] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.6, 2.0, 2.45, 3.05, 3.8, 4.5, 5.3} ;
//------ variable bin size for (all charged) spectra 
 int nSpectrumBins ;
 double kSpectrumBinsLow[100] ;

 
void makeSpectrumBins()
{
 double pt=0. ;
 nSpectrumBins=0 ;
 while(pt<9.0){
  kSpectrumBinsLow[nSpectrumBins] = pt ;
  if(pt<0.99) pt=pt+0.1 ;
  else if(pt<1.199) pt=pt+0.2 ;
  else if(pt<2.399) pt=pt+0.4 ;
  else pt=pt+0.8 ;
  nSpectrumBins++ ;
  cout << "[" << pt << "]" ;
 }
 cout << endl ;
}


Hfill::Hfill(DatabasePDG2 *d, int *pid, int *mpid, float *px, float *py, float *pz, float *e, float *x, float *y, float *z, float *t, TFile* outputFile, char *suffix)
{
  NPART = database->GetNParticles() ;
  NMIX = NEVENTS ;
  database = d ;
  Px = px; Py = py; Pz = pz; E = e ;
  X = x; Y = y; Z = z; T = t;
  id = pid ; mid = mpid ;
  nevents = 0 ;
  const int bufdim = 1000*100 ;
  Pxbuf = new float [bufdim] ;
  Pybuf = new float [bufdim] ;
  Pzbuf = new float [bufdim] ;
  Ebuf = new float [bufdim] ;
  Xbuf = new float [bufdim] ;
  Ybuf = new float [bufdim] ;
  Zbuf = new float [bufdim] ;
  Tbuf = new float [bufdim] ;
  midbuf = new int [bufdim] ;
  mult_id = new Int_t [NPART] ;
  
  outputFile->cd();
  
   hCFnom = new TH3F* [NPT];
   hCFden = new TH3F* [NPT];

   char hname [NCHARFILENAME] ;
   for(int ipt=0; ipt<NPT; ipt++){
     sprintf(hname,"hCFnum%s%s",suffix,binnames[ipt]) ;
     cout << hname << endl ;
     hCFnom [ipt] = new TH3F(hname,hname,NBINS, -QMAX, QMAX,NBINS, -QMAX, QMAX,NBINS, -QMAX, QMAX);
     sprintf(hname,"hCFden%s%s",suffix,binnames[ipt]) ;
     cout << hname << endl ;
     hCFden [ipt] = new TH3F(hname,hname,NBINS, -QMAX, QMAX,NBINS, -QMAX, QMAX,NBINS, -QMAX, QMAX);
   }
   
   hpt = new TH1F* [npid_pt] ;
   havcos2ident = new TH2F* [npid_v2] ;
   for(int i=0; i<npid_pt; i++){
   sprintf(hname,"%s%s",name_pt[i],suffix) ;
   cout << hname << endl ;
   hpt[i] = new TH1F(hname,hname,100,0.,6.0) ;
   }
   for(int i=0; i<npid_v2; i++){
   sprintf(hname,"havcos2%s%s",name_v2[i],suffix) ;
   havcos2ident[i] = new TH2F(hname,hname,sizeof(kPtBins)/sizeof(double)-1,kPtBins,40, -1.0, 1.0) ;
   }
   sprintf(hname,"tdist%s",suffix) ;
   htimedist = new TH2F(hname,hname,100,0.,20.,100,1.,16.) ;
   sprintf(hname,"hptcharged%s",suffix) ;
   makeSpectrumBins() ; // make variable pt bins for all charged spectrum
   hptch = new TH1F(hname,hname,nSpectrumBins-1,kSpectrumBinsLow) ;
   sprintf(hname,"hptcharged%s_phi",suffix) ;
   hptphich = new TH2F(hname,hname,sizeof(kPtBins)/sizeof(double)-1,kPtBins,40, -TMath::Pi(), TMath::Pi()) ;
   sprintf(hname,"havcos2%s",suffix) ;
   havcos2 = new TH2F(hname,hname,sizeof(kPtBins)/sizeof(double)-1,kPtBins,40, -1.0, 1.0) ;
   
   sprintf(hname,"tree%s",suffix) ;
   tree = new TTree(hname,hname) ;
   tree->Branch("dndeta",&dndeta,"dndeta/I");
   tree->Branch("dndy",&dndy,"dndy/I");
   tree->Branch("NPART",&NPART,"NPART/I") ;
//   tree->Branch("mult_id",&mult_id[0],"mult_id[NPART]/I") ;
}


void Hfill::fillCF3(int pcount, int ievent)
{
 TLorentzVector m1 ;
 TLorentzVector m2 ;
 TVector3 vlong ;
 
 if(ievent%NMIX==0) bufcount = 0 ; // begin filling
 for (Int_t i=0;i<pcount;i++){
 float Pi = sqrt(Px[i]*Px[i]+Py[i]*Py[i]+Pz[i]*Pz[i]) ;
 if(id[i]==PID && (Pi+fabs(Pz[i]))/(Pi-fabs(Pz[i]))<2.7){ // |eta|<0.5
  Pxbuf[bufcount] = Px[i] ;
  Pybuf[bufcount] = Py[i] ;
  Pzbuf[bufcount] = Pz[i] ;
  Ebuf[bufcount] = E[i] ;
  Xbuf[bufcount] = X[i] ;
  Ybuf[bufcount] = Y[i] ;
  Zbuf[bufcount] = Z[i] ;
  Tbuf[bufcount] = T[i] ;
  midbuf[bufcount] = mid[i] ;
  bufcount++ ;
 }
 }

     if(ievent%NMIX==NMIX-1) 
     for (Int_t i=0;i<bufcount;i++)
     for (Int_t j=0; j<bufcount; j++){ // previously: j=i+1
         if(i==j) continue ;
		m1.SetXYZT(Pxbuf[i],Pybuf[i],Pzbuf[i],Ebuf[i]);
		m2.SetXYZT(Pxbuf[j],Pybuf[j],Pzbuf[j],Ebuf[j]);
		vlong.SetXYZ(0., 0., (Pzbuf[i] + Pzbuf[j]) / (Ebuf[i] + Ebuf[j]));
		m1.Boost(-vlong);
		m2.Boost(-vlong);
              Float_t kt =  0.5 * TMath::Hypot(m1.X() + m2.X(), m1.Y() + m2.Y());
               for(int ipt=0; ipt<NPT; ipt++) // pt loop
                 if (kt>pt_low[ipt] && kt<pt_high[ipt]) {

               Float_t pXsum = m1.X() + m2.X();
               Float_t pYsum = m1.Y() + m2.Y();
               Float_t pXdif = m1.X() - m2.X();
               Float_t pYdif = m1.Y() - m2.Y();
               Float_t qOut = 0.5 * (pXsum * pXdif + pYsum * pYdif) / kt;

               Float_t qSide = (m1.X() * m2.Y() - m1.Y() * m2.X()) / kt;

               Float_t qLong = m1.Z() - m2.Z() ;

		Float_t wCos ;
		if((Tbuf[i]*Tbuf[i]-Zbuf[i]*Zbuf[i])<50.*50. && (Tbuf[j]*Tbuf[j]-Zbuf[j]*Zbuf[j])<50.*50.)
		 wCos = 1 + TMath::Cos(5.068423*((m1.X()-m2.X())*(Xbuf[i]-Xbuf[j])
                  + (m1.Y()-m2.Y())*(Ybuf[i]-Ybuf[j]) + (m1.Z()-m2.Z())*(Zbuf[i]-Zbuf[j])
                  - (Tbuf[i]-Tbuf[j])*(m1.T()-m2.T())));
                  else wCos = 1. ;
                       //for pions
                        if (qOut > -QMAX && qOut < QMAX && qSide > -QMAX && qSide <QMAX && qLong > -QMAX && qLong <QMAX) {
                            hCFnom[ipt]->Fill(qOut,qSide,qLong, wCos);
                            hCFden[ipt]->Fill(qOut,qSide,qLong, 1.);
                          }
                    } // pt loop
                } //for part1, part2
}


void Hfill::fillpt(int pcount)
{
//  int n1=0, n2=0 ;
  dndeta = dndy = 0 ;
  for(int i=0; i<NPART; i++) mult_id[i]=0 ;
  int index ;
  double QnSin=0, QnCos=0 ;
  for(int i=0; i<pcount; i++){ // determine reaction plane angle first
   ParticlePDG2 *part = database->GetPDGParticle(id[i]) ;
   double charge = 0 ;
   if(part!=0) charge = part->GetElectricCharge() ;
   if(charge!=0){
      double phi = atan2(Py[i],Px[i]) ;
      double Pt = sqrt(Px[i]*Px[i]+Py[i]*Py[i]) ;
      double w = Pt<2.0 ? Pt : 2.0 ;
      //if(fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))>0.5) w=0.0 ;
      QnSin += w*sin(2.0*phi) ;
      QnCos += w*cos(2.0*phi) ;
   }
  }
  double Psi2 = 0.5*atan(QnSin/QnCos) ;

  for(int i=0; i<pcount; i++){
	  double Pt = sqrt(Px[i]*Px[i]+Py[i]*Py[i]) ;
	  double phi = atan2(Py[i],Px[i]) ;
      for(int ipid=0; ipid<npid_pt; ipid++){
      if(id[i]==pid_pt[ipid] && Pt>0.) hpt[ipid]->Fill(Pt, // exclude feed-down from weak decays
      1./Pt*( fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.5)) ;
      }
      for(int ipid=0; ipid<npid_v2; ipid++){
      havcos2ident[ipid]->Fill(Pt,cos(2.0*(phi)), abs(id[i])==abs(pid_v2[ipid]) && fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.35) ;
      // TEST, before: Fill(Pt,cos(2.0*(phi-Psi2)),...
	  }
   ParticlePDG2 *part = database->GetPDGParticle(id[i],index) ;
   if(part!=0 && index<NPART && fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.5) mult_id[index]=mult_id[index]+1 ;
   double charge = 0 ;
   if(part!=0) charge = part->GetElectricCharge() ;
   double P = sqrt(Px[i]*Px[i]+Py[i]*Py[i]+Pz[i]*Pz[i]) ;
   if(charge!=0 && fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.5) dndy++ ;
   if(charge!=0 && fabs(0.5*log((P+Pz[i])/(P-Pz[i])))<0.5) dndeta++ ;
   if(charge!=0){
	   hptch->Fill(Pt,
      1./Pt*(fabs(0.5*log((P+Pz[i])/(P-Pz[i])))<0.8)) ; // n_ch versus pt, |eta|<0.8
      double phi = atan2(Py[i],Px[i]) ;
      hptphich->Fill(Pt,phi, fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.35) ;
      havcos2->Fill(Pt,cos(2.0*(phi)), fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.35) ;
      //TEST, before: Fill(Pt,cos(2.0*(phi-Psi2)), ...
   }
   htimedist->Fill(sqrt(X[i]*X[i]+Y[i]*Y[i]),sqrt(T[i]*T[i]-Z[i]*Z[i]), fabs(0.5*log((E[i]+Pz[i])/(E[i]-Pz[i])))<0.35) ;
  }
  tree->Fill() ;
//  cout << "debug 21  " << n1 << "  " << n2 << endl ;
//  cout << "2nd-order event plane tg = " << QnSin/QnCos << endl ;
 nevents++ ;
}


void Hfill::write()
{
     for(int ipt=0; ipt<NPT; ipt++)
 {
   hCFnom[ipt]->Sumw2();
   hCFden[ipt]->Sumw2();
 }
 for(int i=0; i<npid_pt; i++){
  if(hpt[i]->GetEntries()>0) hpt[i]->Scale(1./(2.*TMath::Pi()*hpt[i]->GetBinWidth(1)*2.*0.5)) ; // /nevents before
 }
 // old normalization - for constant pt bins
 // hptch->Scale(1./(2.*TMath::Pi()*hptch->GetBinWidth(1)*2.*0.8)) ; // /nevents before
 for(int i=0; i<hptch->GetNbinsX()+2; i++){
  hptch->SetBinContent(i,hptch->GetBinContent(i)*1./(2.*TMath::Pi()*hptch->GetBinWidth(i)*2.*0.8)) ;
 }
}
