#include <TTree.h>
#include "tree.h"
#include "gen.h"

MyTree::MyTree(char *name)
{
  tree = new TTree(name,name);
  tree->Branch("npart",&gen::npart[0],"npart/I");
//  treefin->Branch("nev",&nev,"nev/I");
  tree->Branch("x",&gen::X[0][0],"x[npart]/F");
  tree->Branch("y",&gen::Y[0][0],"y[npart]/F");
  tree->Branch("z",&gen::Z[0][0],"z[npart]/F");
  tree->Branch("t",&gen::T[0][0],"t[npart]/F");
  tree->Branch("px",&gen::Px[0][0],"px[npart]/F");
  tree->Branch("py",&gen::Py[0][0],"py[npart]/F");
  tree->Branch("pz",&gen::Pz[0][0],"pz[npart]/F");
  tree->Branch("E",&gen::E[0][0],"E[npart]/F");
  tree->Branch("acc",&gen::Acc[0][0],"acc[npart]/I");
  tree->Branch("id",&gen::Id[0][0],"id[npart]/I");
  tree->Branch("mid",&gen::MId[0][0],"mid[npart]/I");
}


void MyTree::setEventAddr(int iev)
{
  tree->SetBranchAddress("npart",&gen::npart[iev]);
  tree->SetBranchAddress("x",&gen::X[iev][0]);
  tree->SetBranchAddress("y",&gen::Y[iev][0]);
  tree->SetBranchAddress("z",&gen::Z[iev][0]);
  tree->SetBranchAddress("t",&gen::T[iev][0]);
  tree->SetBranchAddress("px",&gen::Px[iev][0]);
  tree->SetBranchAddress("py",&gen::Py[iev][0]);
  tree->SetBranchAddress("pz",&gen::Pz[iev][0]);
  tree->SetBranchAddress("E",&gen::E[iev][0]);
  tree->SetBranchAddress("acc",&gen::Acc[iev][0]);
  tree->SetBranchAddress("id",&gen::Id[iev][0]);
  tree->SetBranchAddress("mid",&gen::MId[iev][0]);
}


void MyTree::fill()
{
 tree->Fill() ;
}
