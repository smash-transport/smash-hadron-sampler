class TTree ;

class MyTree{
 TTree *tree ;
 Float_t *X, *Y, *Z, *T, *Px, *Py, *Pz, *E ;
 Int_t *Id, *MId, *LastColl, *NColl, *Origin ;
 Short_t *Chrg, *Bar, *Strg ;
 Int_t nfill ;
public:
 MyTree(const char *name) ;
 void fill(int iev) ;
} ;
