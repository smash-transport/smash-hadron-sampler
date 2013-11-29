#include <TTree.h>
#include "tree.h"
#include "gen.h"

MyTree::MyTree(char *name)
{
  tree = new TTree(name,name);
  tree->Branch("npart",&gen::npart,"npart/I");
//  treefin->Branch("nev",&nev,"nev/I");
  tree->Branch("x",&gen::X[0],"x[npart]/F");
  tree->Branch("y",&gen::Y[0],"y[npart]/F");
  tree->Branch("z",&gen::Z[0],"z[npart]/F");
  tree->Branch("t",&gen::T[0],"t[npart]/F");
  tree->Branch("px",&gen::Px[0],"px[npart]/F");
  tree->Branch("py",&gen::Py[0],"py[npart]/F");
  tree->Branch("pz",&gen::Pz[0],"pz[npart]/F");
  tree->Branch("E",&gen::E[0],"E[npart]/F");
  tree->Branch("acc",&gen::Acc[0],"acc[npart]/I");
  tree->Branch("id",&gen::Id[0],"id[npart]/I");
  tree->Branch("mid",&gen::MId[0],"mid[npart]/I");
}


void MyTree::fill()
{
 tree->Fill() ;
}
