class TRandom3 ;
class DatabasePDG2;

namespace gen{
// data
extern DatabasePDG2 *database ;
extern TRandom3 *rnd ;
extern Float_t *Px, *Py, *Pz, *E, *X, *Y, *Z, *T ; // particle arrays
extern Int_t *Acc, *Id, *MId ;
extern Int_t npart ;

// functions
void load(char *filename, int N) ;
double calcDFMax(int pindex, char *fileout) ;
void loadDFMax(char *filename, int N) ;
int generate(void) ;
}

