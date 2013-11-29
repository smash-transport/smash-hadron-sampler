class TRandom3 ;
class DatabasePDG2;

namespace gen{
// data
extern DatabasePDG2 *database ;
extern TRandom3 *rnd ;
extern float **Px, **Py, **Pz, **E, **X, **Y, **Z, **T ; // particle arrays
extern int **Acc, **Id, **MId ;
extern char **Charge ;
extern int *npart ;

// functions
void load(char *filename, int N) ;
double calcDFMax(int pindex, char *fileout) ;
void loadDFMax(char *filename, int N) ;
int generate() ;
}

