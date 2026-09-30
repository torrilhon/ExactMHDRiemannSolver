//==============================================================================
//  Implementation des klassischen Newton-Verfahrens.
//  Die Jacobi-Matrix wird mit Vorwärts-Differenzen bereitgestellt.
//==============================================================================

//#include <conio.h>
#include <stdio.h>
#include <math.h>

#define EPS_Nwtn 1e-5
#define EPS_Jac 1e-7
#define EPS_Inv 1e-10

void GetJacobi( double *(*jac), double *x, double *f1, double *f2,
                void fkt(double*,double*,int), int n );
void LinSolve( double *(*A), double *x, double *b, int n );
int Inv( double *(*inverse), int n );
double NormInf( double *f, int n );

// Implementation:

int Newton( void fkt(double*,double*,int), double *x, int n )
{
  double *f  = new double[n],
         *f2 = new double[n],
         *dx = new double[n];
  double *(*Df);
  int N=0, Nmax = 100;

  Df = new double*[n];
  for( int i=0; i<n; i++ )
    Df[i] = new double[n];

  do {
    N++;
    GetJacobi(Df,x,f,f2,fkt,n);
    LinSolve(Df,dx,f,n);
    for( int i=0; i<n; i++ ) {
      x[i] -= dx[i];
      printf( "%14.10f", f[i] );
    };
    printf( "\n" ); //getch();
  }
  while( (NormInf(f,n) > EPS_Nwtn) && (N < Nmax) );
  
  return( N );
};

void GetJacobi( double *(*jac), double *x, double *f1, double *f2,
                void fkt(double*,double*,int), int n )
{
  double tmp, h;
  
  fkt( f1, x, n );
  for( int j=0; j<n; j++ ) {
    tmp = x[j];
    h = EPS_Jac*fabs(tmp);
    if( h == 0 ) h = EPS_Jac;
    x[j] = tmp+h;
    fkt( f2, x, n );
    x[j] = tmp;
    for( int i=0; i<n; i++ )
      jac[i][j] = (f2[i]-f1[i])/h;
  };
};

inline double max( double z1, double z2 )
{
  if( z1 > z2 ) return( z1 );
  else return( z2 );
};

double NormInf( double *f, int n )
{
  double erg = 0;
  for( int i=0; i<n; i++ )
    erg = max(fabs(f[i]),erg);
  return( erg );
};

void LinSolve( double *(*A), double *x, double *b, int n )
{
  if( Inv(A,n) == 3 ) {
    printf( " Vorsicht: Singuläre Matrix in LinSolve()!\n" );
    //    getch();
  };
  for( int i=0; i<n; i++ ) {
    x[i] = 0;
    for( int j=0; j<n; j++ )
      x[i] += A[i][j]*b[j];
  };
};

// Alte Matrix-Inversion-Funktion: inverse enthält die Matrix als Input
// und ihre Inverse als Output. Inv liefert 3 bei Singularität.
int Inv( double *(*inverse), int n )
{
  int    *permx,*permy;                   /*  Zeilen/Spaltenpermutationen   */
  int register k, j, i, ix, iy;           /*  Schleifenindizes      */
  int      nx, ny;                        /*  Pivotindizes          */
  char     *malloc();
  double   piv, temp, faktor, *temp1;     /*  Hilfsgroessen         */
  void     free();

  if ( n < 1 ) return(1);          	/* Unzulaessige Eingabeparameter   */
  for (k = 0; k < n; k++)
    if ( inverse[k] == NULL ) return(1);

  if ( n == 1 ) {
    inverse[0][0] = 1/inverse[0][0];
    return( 0 );
  };

  permx = new int [n]; 		 	   /*  Speicher allokieren  */
  if ( permx == NULL ) return(2);
  permy = new int [n];            	   /*  Speicher allokieren  */
  if ( permx == NULL ) return(2);

  for (i = 0; i < n; i++)
    permx[i] = permy[i] = -1;              /*  permx, permy initial.*/

  for (i = 0; i < n; i++)
    { for (piv = 0.0, ix = 0; ix < n; ix++)/*  Suche aktuelles      */
      if ( permx[ix] == -1 )               /*  Pivotelement         */
	{ for (iy = 0; iy < n; iy++)
	  if ( permy[iy] == -1 && fabs(piv) < fabs(inverse[ix][iy]) ) {
	    piv = inverse[ix][iy];  /* merke aktuelle Pivotpos.   */
	    nx = ix; ny = iy;       /* u. deren Indizes           */
	  }
	}
    if ( fabs(piv) < EPS_Inv ) {      /* Wenn piv zu klein, so ist  */
      delete permx;                   /* matrix nahezu singulaer    */
      delete permy;
      return(3);
    }
    
    permx[nx] = ny; permy[ny] = nx;   /* Tausche Pivotpositionen    */
    
    temp = 1.0 / piv;                 /* Pivotschritt ...           */
    for (j = 0; j < n; j++)
      if ( j != nx )
	{ faktor = inverse[j][ny] * temp;
	for (k = 0; k < n; k++)       /* ... ausserhalb von Pivot- */
	                              /*     zeile u. -spalte      */
	  inverse[j][k] -= inverse[nx][k] * faktor;
	
	inverse[j][ny] = faktor;      /* ... in der Pivotspalte    */
	}
    for (k = 0; k < n; k++)
      inverse[nx][k] *= -temp;        /* ... in der Pivotzeile     */
    inverse[nx][ny] = temp;           /* ... fuers Pivotelement    */
    
    }   /*  end i */
		     /* Zeilen- u. Spaltenvertauschungen rueckgaengig machen */
  for (i = 0; i < n; i++) {           /* Bestimme j mit permx[j] = i */
      for (j = i; j < n; j++) if (permx[j] == i) break;
      if ( j != i) {
	  temp1 = inverse[i];            /* Zeilenvertauschung   */
	  inverse[i] = inverse[j];       /* nur Zeilenzeiger     */
	  inverse[j] = temp1;            /* tauschen !           */
	  permx[j] = permx[i]; permx[i] = i;
	}
      /* Bestimme j mit permy[j] = i */
      for (j = i; j < n; j++) if (permy[j] == i) break;
      if ( j != i) {
	  for (k = 0; k < n; k++ ) {
	    temp = inverse[k][i];
	    inverse[k][i] = inverse[k][j];    /* Spaltenvertauschung */
	    inverse[k][j] = temp;
	    }
	  permy[j] = permy[i]; permy[i] = i;
	}
    }  /* end i */
  
  if ( permy != NULL ) delete permy;          /* Speicher freigeben  */
  if ( permx != NULL ) delete permx;          /* Speicher freigeben  */
   
  return(0);
};


