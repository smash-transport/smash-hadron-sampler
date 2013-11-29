#include <TH3F.h>
#include <TH2F.h>
#include <TTree.h>
#include <TFile.h>

#include "DatabasePDG2.h"

class Hfill
{
  float *Px, *Py, *Pz, *E, *X, *Y, *Z, *T ;
  float *Pxbuf, *Pybuf, *Pzbuf, *Ebuf, *Xbuf, *Ybuf, *Zbuf, *Tbuf ; // for event mixing
  int NMIX ;
  int *midbuf ; // for event mixing
  int bufcount ;
  int *id, *mid ;
  int dndeta, dndy ;
  int *mult_id ;
  TH3F **hCFnom, **hCFden ;
  TH2F *htimedist ;
  TH1F **hpt, *hptch ;
  TH2F **havcos2ident, *hptphich, *havcos2 ; // pt,phi distribution, for v_2 etc, _2 for canonical v_2 determination <cos(2phi)>
  TTree *tree ;
  TFile* outputFile ;
  int nevents ;
  
  DatabasePDG2 *database ;
  public:
  Hfill(DatabasePDG2 *d, int *pid, int *mpid, float *px, float *py, float *pz, float *e, float *x, float *y, float *z, float *t, TFile* outputFile, char *suffix) ;
  void fillCF3(int pcount, int ievent) ;
  void fillpt(int pcount) ;
  void write() ;
};
