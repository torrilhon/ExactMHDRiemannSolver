//==============================================================================
//  sShock(In,pVar,Btc,posneg) liefert die Feldwerte nach einem langsamen
//  Shock mit Stärke pVar in Richtung posneg. Btc ist der größte per Shock
//  erlaubte Sprung in Bt. Ist pVar kleiner als Btc beschreibt es den Sprung
//  in Bt, ist es größer beschreibt pVar-Btc die Stärke des sich an den
//  Shock bis Btc anschließenden Fächers (Compound-Wave).
//  Verfahren: Für den Shock lassen sich die Felder mit Hilfe von vv(Bt,.)
//  direkt ausrechnen. Für den Fächer wird rarefaction() aufgerufen.
//==============================================================================

//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MHD1.h"

extern double gam;
extern double kap;

static double vv( double Bt, double A, double B2 );
static double gM2( double v, double Bt, double A, double B2 );
static double pH( double v, double Bt, double A );
static double sign( double a, double b );

// Implementation:

cField sShock( cField In, double pVar, int posneg, double *Ms )
{
  double srtp = sqrt( In.p ), srta = srtp/sqrt( In.ro ), srtm;
  double A, B, B2, v, M, Bt, Btc;
  cField Out1, Out2;

  B  = In.Bn / srtp;
  B2 = B*B;
  A  = fabs(In.Bt) / srtp;
  Btc = sShockMax(A,B);
  if( pVar > Btc ) Bt = A - Btc;
  else Bt = A - pVar;
  if( A < 0.00001 ) {
    printf( "Vorsicht: Langsamer Shock mit A=0, (%f)\n", A );
    //getch();
  };

  v = vv(Bt,A,B2);
  M = gM2(v,Bt,A,B2);
  srtm = sqrt(M);
  if( Ms ) *Ms = srtm/sqrt(gam);

  Out1.ro = In.ro/v;
  Out1.p  = In.p*pH(v,Bt,A);
  Out1.Bt = Bt*srtp*sign(1,In.Bt);
  Out1.v  = In.v + posneg*srta*srtm*(1-In.ro/Out1.ro);
  Out1.C  = -posneg*In.Bn/(srta*In.ro*srtm); 
  Out1.Bn = In.Bn;

  if( (Out1.Bt < 0) && Ms ) printf( "  Intermediate Wave...\n" );

  if( pVar > Btc ) {
    //printf( "  Compound Wave...\n" );
    Out2 = Rarefaction( Out1, pVar-Btc, 's', posneg );
    Out2.C = (Out2.C*Out2.Bt+(Out1.C-Out2.C)*Out1.Bt-Out1.C*In.Bt)
      /(Out2.Bt-In.Bt);
    return( Out2 );
    };

  return( Out1 );
};

static double vv( double Bt, double A, double B2 )
{
  double aa, bb, cc;
  
  aa = 0.5*Bt*(2*(kap+1)+(Bt-A)*(Bt-A)+kap*(A*A-Bt*Bt))-kap*B2*(Bt-A);
  bb = (1+kap)*(0.5*A*(Bt*Bt-A*A)+B2*(Bt-A)-(Bt+A));
  cc = (1+kap)*A-(A*A+B2)*(Bt-A);
  
  return( -2*cc/(bb-sqrt(bb*bb-4*aa*cc)) );
};

static double gM2( double v, double Bt, double A, double B2 )
{
  return( B2*(Bt-A)/(Bt*v-A) );
};

static double pH( double v, double Bt, double A )
{
  return( (v-kap+0.5*(Bt-A)*(Bt-A)*(v-1))/(1-kap*v) );
};

static double sign( double a, double b )
{
  return (b > 0.0) ? fabs(a) : -fabs(a);
};







