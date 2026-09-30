//==============================================================================
//  fShock(In,pVar,posneg) berechnet Felder nach einem schnellen Shock
//  der Stärke pVar (Machzahl). posneg legt mit +1 oder -1 die Richtung
//  des Shocks (rechts oder links) fest. pVar, d.h. die Machzahl muß
//  größer c_f sein.
//  Verfahren: Fallunterscheidung zwischen Anfangs-Bt ungleich oder gleich
//  null. Schnittpunktsuche zwischen Hugoniot- und Rayleigh-Kurve.
//==============================================================================

//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MHD1.h"

#define TOL 1e-15
#define TOL_Bisek 1e-4
#define EPS 1e-3
#define EPS1 1e-9

extern double gam;
extern double kap;

static double loc_B;
static double loc_B2;
static double loc_A;
static double loc_gM2;   // gam*M^2

static double pH( double v, double Bt, double A );
static double pR( double v, double Bt, double A, double gM2 );
static double Btt( double v, double gM2, double A, double B2 );
static double ff( double v );
static double NST( double f(double), double x1, double x2 );
static double sign( double a, double b );

// Implementation:

cField fShock( cField In, double pVar, int posneg )
{
  double srtp = sqrt( In.p ), srta = srtp/sqrt( In.ro ), srtm;
  double vv, Bt, vvLinks;
  cField Out;

  loc_B  = In.Bn / srtp;
  loc_B2 = loc_B*loc_B;
  loc_A  = fabs(In.Bt) / srtp;
  loc_gM2= gam*pVar*pVar;
  srtm   = sqrt(loc_gM2);

  if( loc_A > EPS ) {
    if( kap > loc_gM2/loc_B2 ) vvLinks = loc_B2/loc_gM2;
    else vvLinks = 1/kap;
    vv = NST( ff, vvLinks+EPS1, 1-EPS1 );
    Bt = Btt(vv,loc_gM2,loc_A,loc_B2);
  }
  else {
    if( loc_gM2 < kap*loc_B2-kap-1 ) {
      vv = loc_B2/loc_gM2;
      Bt = sqrt( 2*(loc_gM2-loc_B2)*((kap*loc_B2-loc_gM2)/(kap-1)-gam)/loc_B2 );
    }
    else {
      vv = (1+(1+kap)/loc_gM2)/kap;
      Bt = 0.0;
    };
  };

  Out.ro = In.ro/vv;
  Out.p  = In.p*pR(vv,Bt,loc_A,loc_gM2);
  Out.Bt = sign(Bt*srtp,In.Bt);
  Out.v  = In.v + posneg*srta*srtm*(1-In.ro/Out.ro);
  Out.C  = -posneg*In.Bn/(srta*In.ro*srtm);     //-posneg*In.Bn/(srta*srtm);    
  Out.Bn = In.Bn;

  return( Out );
};

static double pH( double v, double Bt, double A )
{
  return( (v-kap+0.5*(Bt-A)*(Bt-A)*(v-1))/(1-kap*v) );
};

static double pR( double v, double Bt, double A, double gM2 )
{
  return( 1-gM2*(v-1)+0.5*(A*A-Bt*Bt) );
};

static double Btt( double v, double gM2, double A, double B2 )
{
  return( A*(gM2-B2)/(gM2*v-B2) );
};

static double ff( double v )
{
  return( pH(v,Btt(v,loc_gM2,loc_A,loc_B2),loc_A)
          -pR(v,Btt(v,loc_gM2,loc_A,loc_B2),loc_A,loc_gM2) );
};

static double NST( double f(double), double x1, double x2 )
{
  double xm = 0.5*(x1+x2), dx, f1 = f(x1), f2 = f(x2), fm = f(xm);

  if( f1*f2 > 0 ) printf( "NST Problem" );

  while( fabs( x2-x1 ) > TOL_Bisek*fabs(x2)+TOL_Bisek ) {
    if( f1*fm > 0 ) {             // Bisektion
      x1 = xm;
      f1 = fm;
    }
    else x2 = xm;
    
    xm = 0.5*(x1+x2);
    fm = f(xm);
  };

  while( fabs( x2-x1 ) > TOL*fabs(x2)+TOL ) {
    f2 = f(x2);
    dx = f2*(x2-x1)/(f2-f1);   // Sekanten-Verfahren
    f1 = f2;
    
    x1 = x2;
    x2 = x2 - dx;
  };
  
  return( x2 );
};

static double sign( double a, double b )
{
  return (b > 0.0) ? fabs(a) : -fabs(a);
};




