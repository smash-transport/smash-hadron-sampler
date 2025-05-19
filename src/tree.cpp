#include "tree.h"

#include <TTree.h>

#include "gen.h"

MyTree::MyTree(const char *name) {
  tree = new TTree(name, name);
  X = new Float_t[gen::NPartBuf];
  Y = new Float_t[gen::NPartBuf];
  Z = new Float_t[gen::NPartBuf];
  T = new Float_t[gen::NPartBuf];
  Px = new Float_t[gen::NPartBuf];
  Py = new Float_t[gen::NPartBuf];
  Pz = new Float_t[gen::NPartBuf];
  E = new Float_t[gen::NPartBuf];
  Id = new Int_t[gen::NPartBuf];
  MId = new Int_t[gen::NPartBuf];
  LastColl = new Int_t[gen::NPartBuf];
  NColl = new Int_t[gen::NPartBuf];
  Origin = new Int_t[gen::NPartBuf];
  Chrg = new Short_t[gen::NPartBuf];  // particle's electric charge
  Bar = new Short_t[gen::NPartBuf];   // baryon charge
  Strg = new Short_t[gen::NPartBuf];  // strangeness

  tree->Branch("npart", &nfill, "npart/I");
  //  treefin->Branch("nev",&nev,"nev/I");
  tree->Branch("x", &X[0], "x[npart]/F");
  tree->Branch("y", &Y[0], "y[npart]/F");
  tree->Branch("z", &Z[0], "z[npart]/F");
  tree->Branch("t", &T[0], "t[npart]/F");
  tree->Branch("px", &Px[0], "px[npart]/F");
  tree->Branch("py", &Py[0], "py[npart]/F");
  tree->Branch("pz", &Pz[0], "pz[npart]/F");
  tree->Branch("E", &E[0], "E[npart]/F");
  tree->Branch("id", &Id[0], "id[npart]/I");
  tree->Branch("mid", &MId[0], "mid[npart]/I");
  tree->Branch("lastcoll", &LastColl[0], "lastcoll[npart]/I");
  tree->Branch("ncoll", &NColl[0], "ncoll[npart]/I");
  tree->Branch("origin", &Origin[0], "orig[npart]/I");
  tree->Branch("ele", &Chrg[0], "ele[npart]/S");
  tree->Branch("bar", &Bar[0], "bar[npart]/S");
  tree->Branch("str", &Strg[0], "str[npart]/S");
}

void MyTree::fill(int iev) {
  nfill = 0;
  for (int ipart = 0; ipart < gen::npart[iev]; ipart++)
    if (gen::pList[iev][ipart]->is_hadron()) {
      X[nfill] = gen::pList[iev][ipart]->position().x1();
      Y[nfill] = gen::pList[iev][ipart]->position().x2();
      Z[nfill] = gen::pList[iev][ipart]->position().x3();
      T[nfill] = gen::pList[iev][ipart]->position().x0();
      Px[nfill] = gen::pList[iev][ipart]->momentum().x1();
      Py[nfill] = gen::pList[iev][ipart]->momentum().x2();
      Pz[nfill] = gen::pList[iev][ipart]->momentum().x3();
      E[nfill] = gen::pList[iev][ipart]->momentum().x0();
      Id[nfill] = gen::pList[iev][ipart]->pdgcode().get_decimal();
      MId[nfill] = 0;       // Used to be set to 0 in the UrQMD version
      LastColl[nfill] = 0;  // Used to be set to 0 in the UrQMD version
      NColl[nfill] = 0;     // Used to be set to 0 in the UrQMD version
      Origin[nfill] = 0;    // Used to be set to 0 in the UrQMD version
      Chrg[nfill] =
          static_cast<Char_t>(gen::pList[iev][ipart]->type().charge());
      Bar[nfill] =
          static_cast<Char_t>(gen::pList[iev][ipart]->type().baryon_number());
      Strg[nfill] =
          static_cast<Char_t>(gen::pList[iev][ipart]->type().strangeness());
      nfill++;
    }
  tree->Fill();
}
