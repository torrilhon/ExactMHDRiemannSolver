//==============================================================================
//  Residuum(L,R,res,pVar) liefert das Residuum res[5] für das Riemannproblem
//  mit Feldwerten L[8] und R[8] links und rechts und den Pfadvariablen pVar[5].
//  pVar[0]-pVar[4] gehören zu den Wellen: Schnell links, Langsam links,
//  Rotation, Langsam rechts, Schnell rechts.
//  Verfahren: Es werden die Feldwerte in der Mitte jeweils von rechts
//  (posneg = -1) und von links (posneg = +1) ausgehend mit Half() berechnet.
//  Das Residuum ist die Differenz der Felder Druck, senkr. Geschwindigkeit,
//  Betrag des tang. B-Feldes und beide Komponenten der tang. Geschwindigkeit
//  der Rechnungen von links und von rechts.
//==============================================================================

//#include <conio.h>
#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include "MHD1.h"

extern double gam;
extern double kap;

static double* Half( cField U, double alpha, double pVarS, double pVarF,
                               double alpha0, int posneg );
static double cf( double A, double B );

// Implementation:

void Residuum( cField U0, double alpha0, cField U1, double alpha1, double *vt,
               double* res, double* pVar )
{
  double *W0, *W1;

  /* ----reduziertes Riemann-Problem---- *

  W0 = Half(U0,M_PI+alpha1,pVar[1],pVar[0],alpha0,-1);
  if( W0[2] < 0.0 ) {
    W0[2]= fabs(W0[2]);
    res[4] = pVar[2]-(M_PI+alpha0);
  };
  W1 = Half(U1,M_PI+alpha0,pVar[3],pVar[4],alpha1,+1);
  if( W1[2] < 0.0 ) {
    W1[2]= fabs(W1[2]);
    res[4] = pVar[2]-(M_PI+alpha1);
  };

  W0[3] += vt[0];

  for( int i=0; i<4; i++ )
    res[i] = W1[i]-W0[i];

  *----------------------------------*/
  /* ----Standard Riemann-Problem---- */

  W0 = Half(U0,pVar[2],pVar[1],pVar[0],alpha0,-1);
  W1 = Half(U1,pVar[2],pVar[3],pVar[4],alpha1,+1);

  W0[3] += vt[0];
  W0[4] += vt[1];

  for( int i=0; i<5; i++ )
    res[i] = W1[i]-W0[i];

  /*----------------------------------*/

  delete(W0);
  delete(W1);
};

static double* Half( cField U, double alpha, double pVarS, double pVarF,
                               double alpha0, int posneg )
{
  cField Fast, Slow;
  double* W = new double[5];
  double ca  = cos(alpha),  sa  = sin(alpha),
         ca0 = cos(alpha0), sa0 = sin(alpha0), srtp = sqrt( U.p );

  if( pVarF > 0 ) {
    Fast = fShock(U,cf(U.Bt/srtp,U.Bn/srtp)+pVarF,posneg);
    //printf( " ro = %.10f, v = %.10f, p = %.10f\n Bt = %.10f, C = %.10f\n",Fast.ro, Fast.v,Fast.p, Fast.Bt, Fast.C ); getch();
  }
  else {
    Fast = Rarefaction(U,fRareMax(U.Bt/srtp,U.Bn/srtp)*tanh(-pVarF),'f',posneg);
    //printf( " ro = %.10f, v = %.10f, p = %.10f\n Bt = %.10f, C = %.10f\n",Fast.ro, Fast.v,Fast.p, Fast.Bt, Fast.C ); getch();
  };

  if( pVarS > 0 ) {
    Slow = sShock(Fast,pVarS,posneg);
    //printf( " ro = %.10f, v = %.10f, p = %.10f\n Bt = %.10f, C = %.10f\n",Slow.ro, Slow.v,Slow.p, Slow.Bt, Slow.C ); getch();
  }
  else {
    Slow = Rarefaction(Fast,-pVarS,'s',posneg );
    //printf( " ro = %.10f, v = %.10f, p = %.10f\n Bt = %.10f, C = %.10f\n",Slow.ro, Slow.v,Slow.p, Slow.Bt, Slow.C ); getch();
  };

  W[0] = Slow.p;
  W[1] = Slow.v;
  W[2] = Slow.Bt; //fabs(Slow.Bt); 
  W[3] = Fast.C*(Fast.Bt-   U.Bt)*ca0 -posneg*Fast.Bt*(ca-ca0)/sqrt(Fast.ro)
           + Slow.C*(Slow.Bt-Fast.Bt)*ca;
  W[4] = Fast.C*(Fast.Bt-   U.Bt)*sa0 -posneg*Fast.Bt*(sa-sa0)/sqrt(Fast.ro)
           + Slow.C*(Slow.Bt-Fast.Bt)*sa;

  if( Slow.Bt < 0 )
    W[3] = Fast.C*(Fast.Bt-   U.Bt)*ca0 + Slow.C*(Slow.Bt-Fast.Bt)*ca0;

  //printf( " v2 = %.10f v3 = %.10f\n", W[3], W[4] ); getch();

  return( W );
};

static double cf( double A, double B )
{
  double tmp = 0.5*(1+(A*A+B*B)/gam);               // lokal normiert
  return( sqrt(tmp+sqrt(tmp*tmp-B*B/gam)) );
};


