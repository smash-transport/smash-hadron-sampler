class TTree ;

class MyTree{
 TTree *tree ;
public:
 MyTree(char *name) ;
 void setEventAddr(int iev) ;
 void fill(void) ;
} ;
