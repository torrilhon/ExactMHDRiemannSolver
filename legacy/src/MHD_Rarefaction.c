//==============================================================================
//  Rarefaction(In,pVar,sfFlag,posneg) berechnet Feldwerte nach einem
//  Fächer der Stärke pVar. sfFlag legt mit 's' oder 'f' fest, ob es
//  sich um einen schnellen (fast) oder langsamen (slow) Fächer handelt.
//  posneg ist entweder -1 oder +1 falls der Fächer nach links (-1) oder
//  nach rechts (+1) läuft.
//  Die Pfadvariable wird keiner Gültigkeitsprüfung mehr unterzogen !
//  Verfahren: Die Glgn für Bt, v und C (d.h. vt) werden mit dopri5
//  hochintegriert. ro und p werden explizit aus pVar ausgerechnet.
//==============================================================================

//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dopri5.h"
#include "MHD1.h"

extern double gam;

static double loc_B;
static double loc_B2;

static void fRHS( unsigned n, double x, double* y, double* f );
static void sRHS( unsigned n, double x, double* y, double* f );
static void check(long nr, double xold, double x, double* y, 
                  unsigned n, int* irtrn);

// Implementation:

cField Rarefaction( cField In, double pVar, char sfFlag, int posneg )
{
  double srtp = sqrt( In.p ), srta = srtp/sqrt( In.ro );
  double A;
  cField Out;

  loc_B  = In.Bn / srtp;
  loc_B2 = loc_B*loc_B;
  A      = In.Bt / srtp;

  double y[3] = { A, 0.0, 0.0 };
  double te = pVar;

  FcnEqDiff rhs = NULL;
  if( sfFlag == 's' ) rhs = sRHS;
  if( sfFlag == 'f' ) rhs = fRHS;

  unsigned n = 3;
  double t0 = 0.0;
  double rtol = 1e-8;
  double atol = 1e-8;
// Steuerungsparameter dopri5:
  int info;                     int itol = 0;
  int iout = 0;                 int meth = 0;
  double uround = 0.0;          double safe = 0.0;
  double fac1 = 0.0;            double fac2 = 0.0;
  double beta = 0.0;            double hmax = 0.0;
  long nmax = 500000;                long nstiff = 0;
  unsigned nrdens = 0;          unsigned licont = 0;
  double h = 0.0;

  info = dopri5(n,rhs,t0,y,te,&rtol,&atol,itol,check,iout,NULL,uround,safe,
		fac1,fac2,beta,hmax,h,nmax,meth,nstiff,nrdens,NULL,licont);

  if( info < 0 ) {
    printf( "\n Vorsicht: dopri5 liefert %d in Rarefaction()", info );
    getch();
  };

  Out.ro = In.ro*exp(-pVar);
  Out.p  = In.p*exp(-gam*pVar);
  Out.v  = In.v - posneg*srta*y[1];
  Out.Bt = y[0]*srtp;
  Out.C  = -posneg*y[2]*srta/(Out.Bt-In.Bt);
  Out.Bn = In.Bn;

  return( Out );
};

static void fRHS( unsigned n, double x, double* y, double* f )
{
  double tmp1, tmp2, tmp3, tmp4;

  tmp1 = gam*exp(-gam*x)/loc_B2;
  tmp2 = 0.5*(tmp1+y[0]*y[0]/loc_B2+1);
  tmp3 = tmp2+sqrt(tmp2*tmp2-tmp1);      // vf^2/va^2
  tmp4 = sqrt(loc_B2*exp(x)*tmp3);       // vf

  f[0] = y[0]/(1/tmp3-1);                // Glg für Bt
  f[1] = tmp4;                           // Glg für v
  f[2] = tmp4*y[0]/loc_B/(1-tmp3);       // Glg für C
};

static void sRHS( unsigned n, double x, double* y, double* f )
{
  double tmp1, tmp2, tmp3, tmp4;

  tmp1 = gam*exp(-gam*x)/loc_B2;
  tmp2 = 0.5*(tmp1+y[0]*y[0]/loc_B2+1);
  tmp3 = tmp2-sqrt(tmp2*tmp2-tmp1);      // vs^2/va^2
  tmp4 = sqrt(loc_B2*exp(x)*tmp3);       // vs

  f[0] = y[0]/(1/tmp3-1);                // Glg für Bt
  f[1] = tmp4;                           // Glg für v
  f[2] = tmp4*y[0]/loc_B/(1-tmp3);       // Glg für C
};

static void check(long nr, double xold, double x, double* y, 
                  unsigned n, int* irtrn)
{
  printf( "h = %.10f  y = %.10f\n", x-xold, y[0] );
};


