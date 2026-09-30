
struct cField
{
  double ro;
  double v;
  double Bt;
  double p;
  double C;
  double Bn;
};

void PrintErg( double *L, double *R, double *pVar );

void Residuum( cField U0, double alpha0, cField U1, double alpha1, double *vt,
               double* res, double* pVar );
int Newton( void fkt(double*,double*,int), double *x, int n );

double fRareMax( double A, double B );
double sShockMax( double A, double B );
cField Rarefaction( cField In, double pVar, char sfFlag, int posneg );
cField sShock( cField In, double pVar, int posneg, double *Ms = NULL );
cField fShock( cField In, double pVar, int posneg );

int getch();
double arctan( double x, double y );
 


