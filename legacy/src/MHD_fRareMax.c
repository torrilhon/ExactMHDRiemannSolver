//==============================================================================
//  fRareMax( B, A ) findet den Maximalwert der Pfadvariable für einen
//  schnellen Fächer zum normalen M-feld B und Anfangswert A für Bt.
//  B und A sind lokal normiert.
//  Verfahren: Die DGL für Bt wird mit abgespeckten dopri5 integriert
//  bis die relative Schrittweite unter den Wert hTOL sinkt und damit
//  die Singularität markiert.
//==============================================================================

//#include <conio.h>
#include <math.h>
#include <stdio.h>
#include "MHD1.h"

extern double gam;

static long nstep;
static double loc_B2;

static double rhsBt( double x, double y );
static double sign( double a, double b );
static double min_d( double a, double b );
static double max_d( double a, double b );
static double hinit( double fcn(double,double), double x, double* y, 
                     double posneg, int iord, double hmax, 
                     double atol, double rtol, double* k1);
static double FindMax( double fcn(double,double), double x, double* y,
                       double rtol, double atol );

// Implementation:

double fRareMax( double A, double B )
{
  double t0 = 0.0, tMax, y = A, rtol = 1e-8, atol = 1e-8;

  loc_B2 = B*B;

  if( fabs(A) < 0.0001 ) {
    printf("\n Vorsicht: Probleme in fRareMax(), A = 0 !" );
    //getch();
  };

  tMax = FindMax(rhsBt,t0,&y,rtol,atol);

  if( tMax < 0.0 /*|| y < 0.0*/ ) {
    //printf("  Vorsicht: Probleme in fRareMax(), nstep > nmax !" );
    //printf("  y = %.10f, t = %.10f, n = %ld\n", y, tMax, nstep );
    tMax = fabs(tMax);
  };

  return( tMax );
};

static double rhsBt( double x, double y )
{
  double tmp1, tmp2, tmp3;

  tmp1 = gam*exp(-gam*x)/loc_B2;
  tmp2 = 0.5*(tmp1+y*y/loc_B2+1);
  tmp3 = tmp2+sqrt(tmp2*tmp2-tmp1);

  return( y/(1/tmp3-1) );
};


static double sign (double a, double b)
{
  return (b > 0.0) ? fabs(a) : -fabs(a);
} /* sign */

static double min_d (double a, double b)
{
  return (a < b)?a:b;
} /* min_d */

static double max_d (double a, double b)
{
  return (a > b)?a:b;
} /* max_d */


static double hinit( double fcn(double,double), double x, double* y, 
                     double posneg, int iord, double hmax, 
                     double atol, double rtol, double* k1)
{
  double dnf, dny, sk, h, h1, der2, der12, sqr, yy1, k2;

  dnf = 0.0;
  dny = 0.0;

  sk = atol + rtol * fabs(y[0]);
  sqr = k1[0] / sk;
  dnf += sqr*sqr;
  sqr = y[0] / sk;
  dny += sqr*sqr;

  if ((dnf <= 1.0E-10) || (dny <= 1.0E-10)) h = 1.0E-6;
  else h = sqrt (dny/dnf) * 0.01;

  h = min_d (h, hmax);
  h = sign (h, posneg);

  /* perform an explicit Euler step */
  yy1 = y[0] + h * k1[0];
  k2 = fcn( x+h, yy1 );

  /* estimate the second derivative of the solution */
  der2 = 0.0;
  sqr = (k2 - k1[0]) / sk;
  der2 += sqr*sqr;
  der2 = sqrt (der2) / h;

  /* step size is computed such that h**iord*max_d(norm(f0),norm(der2))=0.01 */
  der12 = max_d (fabs(der2), sqrt(dnf));
  if (der12 <= 1.0E-15) h1 = max_d (1.0E-6, fabs(h)*1.0E-3);
  else h1 = pow (0.01/der12, 1.0/(double)iord);
  h = min_d (100.0 * h, min_d (h1, hmax));

  return sign (h, posneg);
} /* hinit */


/* front-end  und  core integrator von dopri5 */
// modified for present purpose
// return -x if nstep > nmax
static double FindMax( double fcn(double,double), double x, double* y,
                       double rtol, double atol )
{
  double yy1, k1, k2, k3, k4, k5, k6, ysti;

  /* initialisations */
  nstep = 0;

  int nmax = 10000;
  double hTOL = 2.3E-12;
  double safe = 0.9;
  double fac1 = 0.2;
  double fac2 = 10.0;
  double beta = 0.04;
  double hmax = 10.0;

  /* Aufruf dopcor */
  double   facold, expo1, fac, facc1, facc2, fac11, posneg, xph;
  double   h, err, sk, hnew;
  double   sqr;
  int      iord, reject;
  double   c2, c3, c4, c5, e1, e3, e4, e5, e6, e7;
  double   a21, a31, a32, a41, a42, a43, a51, a52, a53, a54;
  double   a61, a62, a63, a64, a65, a71, a73, a74, a75, a76;

  /* initialisations */
  c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;
  a21=0.2, a31=3.0/40.0, a32=9.0/40.0;
  a41=44.0/45.0, a42=-56.0/15.0; a43=32.0/9.0;
  a51=19372.0/6561.0, a52=-25360.0/2187.0;
  a53=64448.0/6561.0, a54=-212.0/729.0;
  a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0;
  a64=49.0/176.0, a65=-5103.0/18656.0;
  a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0;
  a75=-2187.0/6784.0, a76=11.0/84.0;
  e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0;
  e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;

  facold = 1.0E-4;
  expo1 = 0.2 - beta * 0.75;
  facc1 = 1.0 / fac1;
  facc2 = 1.0 / fac2;
  posneg = 1.0;

  /* initial preparations */
  k1 = fcn( x, y[0] );
  hmax = fabs (hmax);
  iord = 5;
  h = hinit (fcn, x, y, posneg, iord, hmax, atol, rtol, &k1);
  reject = 0;

  /* basic integration step */
  while (1)
  {
    if( nstep > nmax ) return( -x );
    if( 0.1*fabs(h) <= fabs(x)*hTOL ) return( x );

    nstep++;

    yy1 = y[0] + h * a21 * k1;                    /* the first 6 stages */
    k2 = fcn( x+c2*h, yy1 );
    yy1 = y[0] + h * (a31*k1 + a32*k2);
    k3 = fcn( x+c3*h, yy1 );
    yy1 = y[0] + h * (a41*k1 + a42*k2 + a43*k3);
    k4 = fcn( x+c4*h, yy1 );
    yy1 = y[0] + h * (a51*k1 + a52*k2 + a53*k3 + a54*k4);
    k5 = fcn( x+c5*h, yy1 );
    ysti = y[0] + h * (a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5);
    xph = x + h;
    k6 = fcn( xph, ysti );
    yy1 = y[0] + h * (a71*k1 + a73*k3 + a74*k4 + a75*k5 + a76*k6);
    k2 = fcn( xph, yy1 );
    k4 = h * (e1*k1 + e3*k3 + e4*k4 + e5*k5 + e6*k6 + e7*k2);

    err = 0.0;                                     /* error estimation */
    sk = atol + rtol * max_d (fabs(y[0]), fabs(yy1));
    sqr = k4 / sk;
    err += sqr*sqr;
    err = sqrt (err);

    fac11 = pow (err, expo1);                       /* computation of hnew */
    fac = fac11 / pow(facold,beta);                 /* Lund-stabilization */
    fac = max_d (facc2, min_d (facc1, fac/safe));   /* fac1 <= hnew/h <= fac2 */
    hnew = h / fac;

    if (err <= 1.0) {                               /* step accepted */
      facold = max_d (err, 1.0E-4);

      k1 = k2;
      y[0] = yy1;
      x = xph;

      if (fabs(hnew) > hmax) hnew = posneg * hmax;
      if (reject) hnew = posneg * min_d (fabs(hnew), fabs(h));

      reject = 0;
    }
    else {                                           /* step rejected */
      hnew = h / min_d (facc1, fac11/safe);
      reject = 1;
    };

    h = hnew;
  };
};




