//=============================================================================
//  PrintErg(L,R,pVar) schreibt die Felder fuer die Pfadvariablen pVar[] 
//  in die Datei Erg.dat. L und R sind die Felder rechts und links.
//  Verfahren: Die Felder werden zunaechst in die Arrays ErgR[] und ErgL[]
//  zwischengespeichert. Die Array werden beschrieben indem jeweils von rechts 
//  und von links in ErgHalf() durchgerechnet wird. 
//  Faecher werden mit imax Zwischenpunkten ausgegeben.
//=============================================================================

//#include <conio.h>
//#include <dir.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MHD1.h"

extern double gam;

static void ErgHalf( double *F, double pVarF, double pVarS, double alpha, 
                     int posneg );
static void PutErg( double M, double ro, double v1, double v2, double v3,
                    double B1, double B2, double B3, double p, 
                    int n, int posneg, int flag );
static double cf( double A, double B );
static double cs( double A, double B );
static double ca( double B );

static int lmax;
static int rmax;
static double ErgL[21000][9];
static double ErgR[21000][9];

//double Bfac = sqrt(4*M_PI);
double Bfac = 1.0;;

// Implementation:

void PrintErg( double *L, double *R, double *pVar )
{
  lmax = 0;
  rmax = 0;

  ErgHalf(R,pVar[4],pVar[3],pVar[2],+1);
  ErgHalf(L,pVar[0],pVar[1],pVar[2],-1);

  //  chdir( "~manuel/DATEN/Riemann1D/Path" );
  FILE *fp = fopen( "Erg.dat", "w" );

  for( int i=0; i<lmax; i++ ) {
    for( int j=0; j<9; j++ )
      fprintf( fp, "%15.10f ", ErgL[i][j] );
    fprintf( fp, "\n" );
  };
  for( int i=1; i<=rmax; i++ ) {
    for( int j=0; j<9; j++ )
      fprintf( fp, "%15.10f ", ErgR[rmax-i][j] );
    fprintf( fp, "\n" );
  };
  fclose(fp);
};

static void ErgHalf( double *F, double pVarF, double pVarS, double alpha, 
                     int posneg )
{
  int i, imax = 1, count = 0;
  cField U, UTmp;
  double alpha0, pTmp, Btc, eps = 1.0e-6;
  double v2, v3, M, BtTmp, srtp, srta;

  U.ro = F[0];
  U.v  = F[1];
  U.p  = F[7];
  U.Bt = sqrt(F[5]*F[5]+F[6]*F[6]);
  U.Bn = F[4];
  v2 = F[2];
  v3 = F[3];
  alpha0 = arctan( F[5], F[6] );
  
  PutErg(15*posneg,
     U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha0),U.Bt*sin(alpha0),U.p,count++,posneg,0);

  if( pVarF > 0 ) {
    srtp = sqrt( U.p );
    M = U.v+posneg*(cf(U.Bt/srtp,U.Bn/srtp)+pVarF)*srtp/sqrt(U.ro/gam);
    BtTmp = U.Bt;

    PutErg(M+posneg*eps,
       U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha0),U.Bt*sin(alpha0),U.p,count++,posneg,0);

    U = fShock(U,cf(U.Bt/srtp,U.Bn/srtp)+pVarF,posneg);

    v2 = v2 + U.C*(U.Bt-BtTmp)*cos(alpha0);
    v3 = v3 + U.C*(U.Bt-BtTmp)*sin(alpha0);

    PutErg(M-posneg*eps,
       U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha0),U.Bt*sin(alpha0),U.p,count++,posneg,0);
  }
  else {
    srtp = sqrt( U.p );
    M = U.v+posneg*cf(U.Bt/srtp,U.Bn/srtp)*srtp/sqrt(U.ro/gam);
    BtTmp = U.Bt;
    
    PutErg(M+posneg*eps,
       U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha0),U.Bt*sin(alpha0),U.p,count++,posneg,0);

    pTmp = fRareMax(U.Bt/srtp,U.Bn/srtp)*tanh(-pVarF);
    for( i=0; i<imax; i++ ) {
      U = Rarefaction(U,pTmp/imax,'f',posneg );

      v2 = v2 + U.C*(U.Bt-BtTmp)*cos(alpha0);
      v3 = v3 + U.C*(U.Bt-BtTmp)*sin(alpha0);
      srtp = sqrt( U.p );
      M = U.v+posneg*cf(U.Bt/srtp,U.Bn/srtp)*srtp/sqrt(U.ro/gam);
      BtTmp = U.Bt;

      PutErg(M-posneg*eps,
         U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha0),U.Bt*sin(alpha0),U.p,count++,posneg,i-imax+1);
    };
  };

  /* c-wave Kompatibilitaet */
    Btc = sShockMax(fabs(U.Bt)/sqrt(U.p),U.Bn/sqrt(U.p));
    if( pVarS > Btc ) pTmp = Btc;
    else pTmp = pVarS; 
    UTmp = sShock(U,pTmp,posneg);
    if( UTmp.Bt < 0.0 ) alpha = alpha0;  // Falls c-Wave -> keine Rotation
  /*---------------------------*/

  srtp = sqrt( U.p );
  M = U.v+posneg*ca(U.Bn/srtp)*srtp/sqrt(U.ro/gam);

  PutErg(M+posneg*eps,
     U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha0),U.Bt*sin(alpha0),U.p,count++,posneg,0);

  v2 = v2 -posneg*U.Bt*(cos(alpha)-cos(alpha0))/sqrt(U.ro);
  v3 = v3 -posneg*U.Bt*(sin(alpha)-sin(alpha0))/sqrt(U.ro);

  PutErg(M-posneg*eps,
  U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,0);

  if( pVarS > 0 ) {
    srtp = sqrt( U.p );
    srta = srtp/sqrt(U.ro/gam);
    BtTmp = U.Bt;

    Btc = sShockMax(fabs(U.Bt)/srtp,U.Bn/srtp);
    if( pVarS > Btc ) pTmp = Btc;
    else pTmp = pVarS;

    UTmp = sShock(U,pTmp,posneg,&M);
    M = U.v +posneg*M*srta;

    PutErg(M+posneg*eps,
	   U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,0);

    U = UTmp;
    v2 = v2 + U.C*(U.Bt-BtTmp)*cos(alpha);
    v3 = v3 + U.C*(U.Bt-BtTmp)*sin(alpha);

    PutErg(M-posneg*eps,
	   U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,0);

    BtTmp = U.Bt;
    if( pVarS > Btc )
      for( i=0; i<imax; i++ ) {
	U = Rarefaction(U,(pVarS-Btc)/imax,'s',posneg);

	v2 = v2 + U.C*(U.Bt-BtTmp)*cos(alpha);
	v3 = v3 + U.C*(U.Bt-BtTmp)*sin(alpha);
	srtp = sqrt( U.p );
	M = U.v+posneg*cs(U.Bt/srtp,U.Bn/srtp)*srtp/sqrt(U.ro/gam);
	BtTmp = U.Bt;

	PutErg(M-posneg*eps,
	     U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,i-imax+1);
      };
  }
  else {
    srtp = sqrt( U.p );
    M = U.v+posneg*cs(U.Bt/srtp,U.Bn/srtp)*srtp/sqrt(U.ro/gam);
    BtTmp = U.Bt;

    PutErg(M+posneg*eps,
       U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,0);

    pTmp = -pVarS;
    for( i=0; i<imax; i++ ) {
      U = Rarefaction(U,pTmp/imax,'s',posneg);

      v2 = v2 + U.C*(U.Bt-BtTmp)*cos(alpha);
      v3 = v3 + U.C*(U.Bt-BtTmp)*sin(alpha);
      srtp = sqrt( U.p );
      M = U.v+posneg*cs(U.Bt/srtp,U.Bn/srtp)*srtp/sqrt(U.ro/gam);
      BtTmp = U.Bt;

      PutErg(M-posneg*eps,
       U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,i-imax+1);
    };
  };

  M = U.v;

  PutErg(M+posneg*eps,
     U.ro,U.v,v2,v3,U.Bn,U.Bt*cos(alpha),U.Bt*sin(alpha),U.p,count++,posneg,0);
};

static void PutErg( double M, double ro, double v1, double v2, double v3,
                    double B1, double B2, double B3, double p, 
                    int n, int posneg, int flag )
{
  double tend = 0.08;
  if( posneg > 0 ) {
    rmax++;
    ErgR[n][0] = M*tend;
    ErgR[n][1] = ro;
    ErgR[n][2] = v1;
    ErgR[n][3] = v2;
    ErgR[n][4] = v3;
    ErgR[n][5] = Bfac*B1;
    ErgR[n][6] = Bfac*B2;
    ErgR[n][7] = Bfac*B3;
    ErgR[n][8] = p;
  }
  else {
    lmax++; 
    ErgL[n][0] = M*tend;
    ErgL[n][1] = ro;
    ErgL[n][2] = v1;
    ErgL[n][3] = v2;
    ErgL[n][4] = v3;
    ErgL[n][5] = Bfac*B1;
    ErgL[n][6] = Bfac*B2;
    ErgL[n][7] = Bfac*B3;
    ErgL[n][8] = p;
  };

  if( !flag )
    printf(  "%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f%12.8f\n",M,
	     ro,v1,v2,v3,B1,B2,B3,p);
};


static double cf( double A, double B )
{
  double tmp = 0.5*(1+(A*A+B*B)/gam);
  return( sqrt(tmp+sqrt(tmp*tmp-B*B/gam)) );
};

static double ca( double B )
{
  return( fabs(B)/sqrt(gam) );
};

static double cs( double A, double B )
{
  double tmp = 0.5*(1+(A*A+B*B)/gam);
  return( sqrt(tmp-sqrt(tmp*tmp-B*B/gam)) );
};

//static double cf( double ro, double p, double A, double B )
//{
//  double tmp = 0.5*(A*A+B*B+gam*p)/ro;
//  return( sqrt(tmp+sqrt(tmp*tmp-gam*p*B*B/(ro*ro))) );
//};

