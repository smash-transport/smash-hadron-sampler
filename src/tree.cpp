#include <TTree.h>
#include "tree.h"
#include "ParticlePDG2.h"
#include "gen.h"
#include "particle.h"

MyTree::MyTree(char *name)
{
  tree = new TTree(name,name);
  X  = new Float_t [gen::NPartBuf] ;
  Y  = new Float_t [gen::NPartBuf] ;
  Z  = new Float_t [gen::NPartBuf] ;
  T  = new Float_t [gen::NPartBuf] ;
  Px = new Float_t [gen::NPartBuf] ;
  Py = new Float_t [gen::NPartBuf] ;
  Pz = new Float_t [gen::NPartBuf] ;
  E  = new Float_t [gen::NPartBuf] ;
  Id   = new Int_t [gen::NPartBuf] ;
  MId  = new Int_t [gen::NPartBuf] ;
  LastColl  = new Int_t [gen::NPartBuf] ;
  NColl  = new Int_t [gen::NPartBuf] ;
  Origin  = new Int_t [gen::NPartBuf] ;
  Chrg= new Short_t [gen::NPartBuf] ; // particle's electric charge
  Bar = new Short_t [gen::NPartBuf] ; // baryon charge
  Strg= new Short_t [gen::NPartBuf] ; // strangeness

  tree->Branch("npart",&nfill,"npart/I");
//  treefin->Branch("nev",&nev,"nev/I");
  tree->Branch("x",&X[0],"x[npart]/F");
  tree->Branch("y",&Y[0],"y[npart]/F");
  tree->Branch("z",&Z[0],"z[npart]/F");
  tree->Branch("t",&T[0],"t[npart]/F");
  tree->Branch("px",&Px[0],"px[npart]/F");
  tree->Branch("py",&Py[0],"py[npart]/F");
  tree->Branch("pz",&Pz[0],"pz[npart]/F");
  tree->Branch("E",&E[0],"E[npart]/F");
  tree->Branch("id",&Id[0],"id[npart]/I");
  tree->Branch("mid",&MId[0],"mid[npart]/I");
  tree->Branch("lastcoll",&LastColl[0],"lastcoll[npart]/I");
  tree->Branch("ncoll",&NColl[0],"ncoll[npart]/I");
  tree->Branch("origin",&Origin[0],"orig[npart]/I");
  tree->Branch("ele",&Chrg[0],"ele[npart]/S");
  tree->Branch("bar",&Bar[0],"bar[npart]/S");
  tree->Branch("str",&Strg[0],"str[npart]/S");
}


void MyTree::fill(int iev)
{
//   nfill = 0 ;
//  for(int ipart=0; ipart<gen::npart[iev]; ipart++)
//  if(gen::pList[iev][ipart]->def!=0){
//    X[nfill] = gen::pList[iev][ipart]->x ;
//    Y[nfill] = gen::pList[iev][ipart]->y ;
//    Z[nfill] = gen::pList[iev][ipart]->z ;
//    T[nfill] = gen::pList[iev][ipart]->t ;
//   Px[nfill] = gen::pList[iev][ipart]->px ;
//   Py[nfill] = gen::pList[iev][ipart]->py ;
//   Pz[nfill] = gen::pList[iev][ipart]->pz ;
//    E[nfill] = gen::pList[iev][ipart]->e ;
//   Id[nfill] = gen::pList[iev][ipart]->def->GetPDG() ;
//  MId[nfill] = gen::pList[iev][ipart]->mid ;
//  LastColl[nfill] = gen::pList[iev][ipart]->lastcoll ;
//  NColl[nfill] = gen::pList[iev][ipart]->ncoll ;
//  Origin[nfill] = gen::pList[iev][ipart]->origin ;
// Chrg[nfill] = (Char_t)(gen::pList[iev][ipart]->def->GetElectricCharge()) ;
//  Bar[nfill] = (Char_t)(gen::pList[iev][ipart]->def->GetBaryonNumber()) ;
// Strg[nfill] = (Char_t)(gen::pList[iev][ipart]->def->GetStrangeness()) ;
//    nfill++ ;
//  }
//  tree->Fill() ;
}
