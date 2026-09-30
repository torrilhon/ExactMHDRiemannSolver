//==============================================================================
//  sShockMax(A,B) liefert den maximalen Sprung in Bt für einen langsamen
//  Shock mit Startwerten A und B. Der Sprung ist positiv, Bt wird allerdings
//  verkleinert. A und B sind lokal normiert.
//  Verfahren: vv(Bt,.) liefert das Volumen für einen Shock nach Bt. Ab
//  einer bestimmten Stärke ist die charakteristische Geschw. nach dem Shock
//  nicht mehr kleiner als die Shockgeschw. Das maximale Bt liefert die
//  Nullstelle in der Differenz zwischen den Geschwindigkeiten (ff(.)).
//==============================================================================

//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MHD1.h"

#define TOL 1e-9
#define TOL_Bisek 1e-4
#define EPS 1e-8

extern double gam;
extern double kap;

static double loc_B2;
static double loc_A;

static double BtLinks( double A, double B2 );
static double vv( double Bt, double A, double B2 );
static double pH( double v, double Bt, double A );
static double gM2( double v, double Bt, double A, double B2 );
static double ff( double Bt );
static double NST( double f(double), double x1, double x2 );

// Implementation:

double sShockMax( double A, double B )
{
  loc_A  = fabs(A);
  loc_B2 = B*B;

  return( loc_A-NST( ff, BtLinks(loc_A,loc_B2)+EPS, 0 ) );
};

static double BtLinks( double A, double B2 )
{
  double tmp1, tmp2, tmp3, bb, cc, tmp4, Bt;

  tmp1 = A*A*(kap-3)*(kap-3)-8*B2*(kap-1);
  tmp2 = 2*A*B2*(1-kap*kap)-A*(2+A*A)*(kap+1)*(kap-3);
  tmp3 = B2*(kap-1)*(2*(1+kap-B2*(kap-1))*(1+kap-B2*(kap-1))
	       +A*A*( 5+kap*(3+kap*(-1+kap))+2*(2*B2+A*A)*(kap-1)*(kap-1) ));
  //  tmp3 = B2*(kap-1)*(4*(1+kap-B2*(kap-1))*(1+kap-B2*(kap-1))
  //	       +A*A*((5+kap*(3+kap*(-1+kap)))+2*(2*B2+A*A)*(kap-1)*(kap-1)));

  Bt = (tmp2+4*sqrt(tmp3))/tmp1;

  bb = 0.5*A*(Bt*Bt-A*A)+B2*(Bt-A)-(Bt+A);
  cc = A-(A*A+B2)*(Bt-A)/(1+kap);
  tmp4 = -2*cc/bb;
  if( (tmp4 < 0) || (tmp4 > 1) ) Bt = -A;

  return( Bt );
};

static double vv( double Bt, double A, double B2 )
{
  double aa, bb, cc;

  aa = 0.5*Bt*(2*(kap+1)+(Bt-A)*(Bt-A)+kap*(A*A-Bt*Bt))-kap*B2*(Bt-A);
  bb = (1+kap)*(0.5*A*(Bt*Bt-A*A)+B2*(Bt-A)-(Bt+A));
  cc = (1+kap)*A-(A*A+B2)*(Bt-A);

  return( -2*cc/(bb-sqrt(bb*bb-4*aa*cc)) );
};

static double pH( double v, double Bt, double A )
{
  return( (v-kap+0.5*(Bt-A)*(Bt-A)*(v-1))/(1-kap*v) );
};

static double gM2( double v, double Bt, double A, double B2 )
{
  return( B2*(Bt-A)/(Bt*v-A) );
};

static double ff( double Bt )
{
  double v, tmp2, tmp3;
  v = vv(Bt,loc_A,loc_B2);
  tmp2 = gam*pH(v,Bt,loc_A);
  tmp3 = 0.5*(loc_B2+Bt*Bt+tmp2);

  return( -v*gM2(v,Bt,loc_A,loc_B2)+tmp3-sqrt(tmp3*tmp3-tmp2*loc_B2) );
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
