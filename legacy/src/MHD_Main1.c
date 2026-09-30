//=============================================================================
//  Hauptprogramm
//=============================================================================

//#include <vcl\condefs.h>
//#include <conio.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "MHD1.h"

#include "tools.h"

//#pragma hdrstop
//---------------------------------------------------------------------------
//USERES("MHD_Main1.res");
//USEUNIT("MHD_Residuum.cpp");
//USEUNIT("MHD_Rarefaction.cpp");
//USEUNIT("MHD_fShock.cpp");
//USEUNIT("MHD_sShock.cpp");
//USEUNIT("MHD_fRareMax.cpp");
//USEUNIT("MHD_sShockMax.cpp");
//USEUNIT("DiffGlgn.cpp");
//USEUNIT("Newton.cpp");
//---------------------------------------------------------------------------


double gam = 5./3;
//double gam = 2.0;
double kap = (gam+1)/(gam-1);

cField U0, U1;
double alpha0, alpha1;  
double vt[2];

void Continuation( double *L0, double *R0, double *L1, double *R1, double *pVar, int nstep );
void Load_Data( char *fname, double *L, double *R, double *pVar );
void Write_Data( char *fname, double *L, double *R, double *pVar );
void Load_DataBase( char *fname, double *L, double *R, double *L1, double *R1, 
                    double *pVar, int n );
double arctan( double x, double y );
double norm( double *x, int n );

void MHD_fkt( double *f, double *x, int n )
{
  Residuum(U0,alpha0,U1,alpha1,vt,f,x);
};

//  Vorgehen bei der Berechnung ohne Rotationen:
//   In Half() Absolutwert von Bt betrachten:  W[2] = fabs(Slow.Bt);
//   sowie die Rotation in  W[3] auskommentieren: //-posneg*Fast.Bt*...
//   Ausgabe:
//   In ErgHalf die Rotationen auskommentieren und
//   am Anfang alpha = alpha0 setzen. 
//   (Unklar was tangentialen Geschw. und U.C machen...) 
 

int main(int argc, const char **argv)
{
  double fac = sqrt(4*M_PI);
  char StartFile[300] = "DataBase.dat";
  
  //  double alph = 3.0, r = 1.0;
  //double r1 = 0.9, r0  = 0.9, B1 = 1.25, alph = M_PI;
  //double r1 = 1.0, r0  = r1, B1 = 1.0, alph = M_PI;

//    Felder:   {  ro,  v1,  v2,  v3,  B1,  B2,  B3,  p  }
  //double R[8] = { 1.0, 0.0, 0.0, 0.0, 1.5, r1*cos(alph), r1*sin(alph), 1.0 };
  //double L[8] = { 3.0, 0.0, 0.0, 0.0, 1.5, r1*cos(0.0), r1*sin(0.0), 3.0 };

  //double R[8] = { 0.2, 0.0, 0.0, 0.0, B1, r1*cos(alph), r1*sin(alph), 0.2 };
  //double L[8] = { 1.0, 0.0, 0.0, 0.0, B1, r0*cos(0.0), r0*sin(0.0), 1.0 };

  // Coplanares Problem (c-Wave nach links):
  //double r1 = 0.9, r0  = 0.9, B1 = 1.25, alph = M_PI;
  //double R[8] = { 0.5, 0.0, 0.0, 0.0, B1, r1*cos(alph), r1*sin(alph), 0.5 };
  //double L[8] = { 2.0, 0.0, 0.0, 0.0, B1, r0*cos(0.0), r0*sin(0.0), 2.0 };
  //double pVar[5] = { -0.86, 0.64, M_PI, 0.63, -0.28 };//Compound
  //double pVar[5] = { -0.79, 0.08, M_PI, 0.63, -0.28 };//regular

  // Nonplanares Problem (c-Wave nach links):
  //double r1 = 1.0, r0  = r1, B1 = 1.1, alph = 2.4;
  //double R[8] = { 0.2, 0.0, 0.0, 0.0,          B1, r1*cos(alph), r1*sin(alph), 0.2 };
  //double L[8] = { 1.7, 0.0, 0.0, 1.4968909366, B1, r0*cos(0.0), r0*sin(0.0), 1.7 };
  //double pVar[5] = { -1.1, 0.7, M_PI, 0.9, -0.5 };//Compound
  //double pVar[5] = { -1.2, 0.1, M_PI, 0.7, -0.8 };//Regular

  // Nonplanares Problem (c-Wave nach rechts):
  //double r1 = 0.9, r0  = 0.9, B1 = 1.0, alph = 1.75;
  //double R[8] = { 0.75, 0.0, 0.0, 0.0,          B1, r1*cos(alph), r1*sin(alph), 0.75 };
  //double L[8] = {  3.0, 0.0, 0.0, 2.7269731686, B1, r0*cos(0.0), r0*sin(0.0), 3.0 };
  //double pVar[5] = { -0.31, -0.16, 4.89, 1.95, 0.22 };//Compound
  //double pVar[5] = { -0.31, -0.16, 4.89, 0.21, 0.23 };//regular

  // Funktionierend, klassisch mit Rotation:
  //double r1 = 0.9, r0  = 0.9, B1 = 1.25, alph = M_PI;
  //double R[8] = { 0.2, 0.0, 0.0, 0.0, B1, r1*cos(alph), r1*sin(alph), 0.2 };
  //double L[8] = { 1.7, 0.0, 0.0, 0.0, B1, r0*cos(0.0), r0*sin(0.0), 1.7 };
  //double pVar[5] = { -0.79, 0.08, M_PI, 0.63, -0.28 };

  // from R. Keppens, Fusion Sci. Tech. 45 (2004), gamma = 5/3
  //double L[8] = { 0.5, 0.0, 1.0, 0.1, 1.0, 2.5, 0.0, 1.0 };
  //double R[8] = { 0.1, 0.0, 0.0, 0.0, 1.0, 2.0, 0.0, 0.1 };
  //double pVar[5] = { -0.2, -0.2, -0.01, 0.4, 0.3 };

  // from Brio & Wu, JCP 75, (1988) -> with gam = 2.0 (!!!)
  //double L[8] = { 1.0, 0.0, 0.0, 0.0, 0.75, 1.0, 0.0, 1.0 };
  //double R[8] = { 0.125, 0.0, 0.0, 0.0, 0.75, -1.0, 0.0, 0.1 };
  //double pVar[5] = { -0.5, 0.1, M_PI, 1.2, -0.1 }; //R-Solution
  //double pVar[5] = { -0.5, 1.4, M_PI, 1.2, -0.1 }; //C-Solution

  // konstante Dichte, Druck, magnetisches Shocktube:
  //double L[8] = { 1.0, 0.0, 0.0, 0.0, 1.5, 7*cos(0.5), 7*sin(0.5), 1.0 };
  //double R[8] = { 1.0, 0.0, 0.0, 0.0, 1.5, cos(1.5), sin(1.5), 1.0 };
  //double pVar[5] = { -0.25, 0.8, 0.6, -0.4, 1.7 };

  //double L[8] = { 1.7, 0.0, 0.0, 0.0, 1.25, 0.9, 0.0, 1.7 };
  //double R[8] = { 0.2, 0.0, 0.0, 0.0, 1.25, -0.4, 0.8, 0.2 };
  //double pVar[5] = { -0.79, 0.08, 2.0, 0.63, -0.28 };
  
  double L[8] = { 4.5, 0.0, 0.0, 0.0, 2.00, 1.0, 0.0, 2.8 };
  double R[8] = { 1.0, 1.0, 0.0, 0.0, 2.00, -0.4, 0.1, 0.5 };
  double pVar[5] = { -3.4, -0.75, 0.98, 0.0003, -2.1 };  

  //double L[8] = { 0.5, 0.0, 0.0, 1.0, 1.0, 0.5, 0.0, 4.8 };
  //double R[8] = { 2.9, 1.0, 0.0, 0.0, 1.0, -1.0, 0.9, 0.5 };
  //double pVar[5] = { -0.32, 0.06, 1.7, 1.81, 0.30 };
 
  int nstep = 10;
  if (exists_argument(argc, argv, "-n"))
      nstep = (int)get_float_argument(argc, argv, "-n");

  int write = 0;
  if (exists_argument(argc, argv, "-w"))
      write = (int)get_float_argument(argc, argv, "-w");

  int pVarFlag = 0;
  if (exists_argument(argc, argv, "-p"))
      pVarFlag = (int)get_float_argument(argc, argv, "-p");
      
  double L1[8], R1[8], pVar1[5]; 
  
  int DBproblem = 1;
  if (exists_argument(argc, argv, "-db")) {
    DBproblem = (int)get_float_argument(argc, argv, "-db");
    Load_DataBase( StartFile, L, R, L1, R1, pVar, DBproblem );
    Continuation( L, R, L1, R1, pVar, nstep ); 
    for( int i=0; i<8; i++ ) {
      L[i] = L1[i];  R[i] = R1[i];
    };
  };        
      
  if (exists_argument(argc, argv, "-f")) {
    get_string_argument(argc, argv, "-f",StartFile);
    Load_Data( StartFile, L1, R1, pVar1 );
    if( pVarFlag ) 
      for( int i=0; i<5; i++ ) pVar[i] = pVar1[i];
    else
      Continuation( L, R, L1, R1, pVar, nstep ); 
    for( int i=0; i<8; i++ ) {
      L[i] = L1[i];  R[i] = R1[i];
    };
  };

    
  U0.ro = L[0];
  U0.v  = L[1];
  U0.p  = L[7];
  U0.Bt = sqrt(L[5]*L[5]+L[6]*L[6]);
  U0.Bn = L[4];
  alpha0 = arctan( L[5], L[6] );

  U1.ro = R[0];
  U1.v  = R[1];
  U1.p  = R[7];
  U1.Bt = sqrt(R[5]*R[5]+R[6]*R[6]);
  U1.Bn = R[4];
  alpha1 = arctan( R[5], R[6] );

  vt[0] = L[2]-R[2];
  vt[1] = L[3]-R[3];

   //U0 = fShock(U1,cf(U1.Bt/sqrt(U1.p),U1.Bn/sqrt(U1.p))+0.0001,1);
   //printf( " ro = %.10f, v = %.10f, p = %.10f, Bt = %.10f\n",U0.ro,U0.v,U0.p,U0.Bt); getch();

   //Residuum(U0,alpha0,U1,alpha1,f,pVar);

  printf( "\n MHD-Riemann-Loeser:\n\n  Residuum:%f, %f\n", alpha0, alpha1 );

  Newton(MHD_fkt,pVar,5);

  printf( "\n Loesung (Pfadvariablen): " );

  printf( " %12.8f,  %12.8f,  %12.8f,  %12.8f,  %12.8f\n", pVar[0],pVar[1],pVar[2],pVar[3],pVar[4]);



  printf( "\n Ergebnisse: x, rho, v1, v2, v3, B1, B2, B3, p\n\n" );
  PrintErg(L,R,pVar);

  double f[5];
  MHD_fkt( f, pVar, 5 );
  if( (norm(f,5) < 0.0001) && write ) 
      Write_Data( "DataBase.dat", L, R, pVar );
  
  return( 0 );
};

double arctan( double x, double y )
{
  if( fabs(x) < 1e-7 )
    if( y > 0 ) return( M_PI/2 );
    else return( -M_PI/2 );
  if( x > 0 ) return( atan(y/x) );
  else return( M_PI+atan(y/x) );
};

double norm( double *x, int n )
{
  double sum = 0;
  
  for( int i=0; i<n; i++ )
    sum += x[i]*x[i];
    
  return( sqrt(sum) );
};
    

int getch()
{
  return( getchar() );
};


/* Continuation Method */

void Continuation( double *L0, double *R0, double *L1, double *R1, double *pVar, int nstep )
{
  int n;
  cField Data0_L, Data0_R, Data1_L, Data1_R;
  double beta0_L, beta0_R, beta1_L, beta1_R;
  double vt0[2], vt1[2];
  
  Data0_L.ro = L0[0];
  Data0_L.v  = L0[1];
  Data0_L.p  = L0[7];
  Data0_L.Bt = sqrt(L0[5]*L0[5]+L0[6]*L0[6]);
  Data0_L.Bn = L0[4];
  beta0_L = arctan( L0[5], L0[6] );
    
  Data0_R.ro = R0[0];
  Data0_R.v  = R0[1];
  Data0_R.p  = R0[7];
  Data0_R.Bt = sqrt(R0[5]*R0[5]+R0[6]*R0[6]);
  Data0_R.Bn = R0[4];
  beta0_R = arctan( R0[5], R0[6] );
    
  vt0[0] = L0[2]-R0[2];
  vt0[1] = L0[3]-R0[3];
  
  Data1_L.ro = L1[0];
  Data1_L.v  = L1[1];
  Data1_L.p  = L1[7];
  Data1_L.Bt = sqrt(L1[5]*L1[5]+L1[6]*L1[6]);
  Data1_L.Bn = L1[4];
  beta1_L = arctan( L1[5], L1[6] );
    
  Data1_R.ro = R1[0];
  Data1_R.v  = R1[1];
  Data1_R.p  = R1[7];
  Data1_R.Bt = sqrt(R1[5]*R1[5]+R1[6]*R1[6]);
  Data1_R.Bn = R1[4];
  beta1_R = arctan( R1[5], R1[6] );
    
  vt1[0] = L1[2]-R1[2];
  vt1[1] = L1[3]-R1[3];

  for( n=0; n<=nstep; n++ ) {
    
    U0.ro = Data0_L.ro + n*(Data1_L.ro-Data0_L.ro)/nstep;
    U0.v  = Data0_L.v  + n*(Data1_L.v -Data0_L.v)/nstep;
    U0.p  = Data0_L.p  + n*(Data1_L.p -Data0_L.p)/nstep;
    U0.Bt = Data0_L.Bt + n*(Data1_L.Bt-Data0_L.Bt)/nstep;
    U0.Bn = Data0_L.Bn + n*(Data1_L.Bn-Data0_L.Bn)/nstep;
    alpha0 = beta0_L  +  n*(beta1_L  - beta0_L)/nstep;
    
    U1.ro = Data0_R.ro + n*(Data1_R.ro-Data0_R.ro)/nstep;
    U1.v  = Data0_R.v  + n*(Data1_R.v -Data0_R.v)/nstep;
    U1.p  = Data0_R.p  + n*(Data1_R.p -Data0_R.p)/nstep;
    U1.Bt = Data0_R.Bt + n*(Data1_R.Bt-Data0_R.Bt)/nstep;
    U1.Bn = Data0_R.Bn + n*(Data1_R.Bn-Data0_R.Bn)/nstep;
    alpha1 = beta0_R  +  n*(beta1_R  - beta0_R)/nstep;
    
    vt[0] = vt0[0] + n/nstep*(vt1[0]-vt0[0]);
    vt[1] = vt0[1] + n/nstep*(vt1[1]-vt0[1]);

    printf( "  || Continuation Step %d of %d\n  ||\n  || With Data:\n", n, nstep );
    printf( "  ||       rho,    v1,    v2,    v3,    B1,    B2,    B3,    p\n" );
    printf(  "  ||  L: %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
	     U0.ro,U0.v,U0.p,U0.Bt,U0.Bn,alpha0);
    printf(  "  ||  R: %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
	     U1.ro,U1.v,U1.p,U1.Bt,U1.Bn,alpha1);
    
    printf( "  ||\n  || Residuum:\n" );
    
    Newton(MHD_fkt,pVar,5);
    
    printf( "  ||\n  || Loesung (Pfadvariablen):\n  ||" );
    
    printf( "  %12.8f,  %12.8f,  %12.8f,  %12.8f,  %12.8f\n", pVar[0],pVar[1],pVar[2],pVar[3],pVar[4]);
    printf( "  ||TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n  ||\n" );
  };
};

void Write_Data( char *fname, double *L, double *R, double *pVar )
{
  double rho, vx, vy, vz, Bx, By, Bz, p;
  FILE *fp;

  printf( " \n  Write Solution to %s\n", fname );
  
  fp = fopen( fname, "a" );
  
  rho = L[0]; 
  vx = L[1]; vy = L[2]; vz = L[3];
  Bx = L[4]; By = L[5]; Bz = L[6];
  p = L[7];

//    Felder:   {  ro,  v1,  v2,  v3,  B1,  B2,  B3,  p  }
  fprintf( fp, " %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n", 
      rho, vx, vy, vz, Bx, By, Bz, p );
  
  rho = R[0]; 
  vx = R[1]; vy = R[2]; vz = R[3];
  Bx = R[4]; By = R[5]; Bz = R[6];
  p = R[7];
  
  fprintf( fp, " %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n", 
       rho, vx, vy, vz, Bx, By, Bz, p );
   
  fprintf( fp, " %10.5lf %10.5lf %10.5lf %10.5lf %10.5lf\n", 
      pVar[0], pVar[1], pVar[2], pVar[3], pVar[4] );
      
  fprintf( fp, "...\n" );
      
  fclose(fp);
};

void Load_DataBase( char *fname, double *L, double *R, double *L1, double *R1, 
                    double *pVar, int n )
{
  double rho, vx, vy, vz, Bx, By, Bz, p;
  double p1, p2, p3, p4, p5;
  char tmp[300];
  FILE *fp;

  fp = fopen( fname, "r" );

  for( int i=1; i<=n-1; i++ ) {
    fgets(tmp,300,fp);
    fgets(tmp,300,fp);
    fgets(tmp,300,fp);
    fgets(tmp,300,fp);
  }; 
  
  fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      &rho, &vx, &vy, &vz, &Bx, &By, &Bz, &p );
  
  L[0] = rho; 
  L[1] = vx; L[2] = vy; L[3] = vz;
  L[4] = Bx; L[5] = By; L[6] = Bz;
  L[7] = p;
  
  fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      &rho, &vx, &vy, &vz, &Bx, &By, &Bz, &p );
  
  R[0] = rho; 
  R[1] = vx; R[2] = vy; R[3] = vz;
  R[4] = Bx; R[5] = By; R[6] = Bz;
  R[7] = p;

  fscanf( fp, " %lf %lf %lf %lf %lf\n",
      &p1, &p2, &p3, &p4, &p5 );
      
  pVar[0] = p1; pVar[1] = p2; pVar[2] = p3; 
  pVar[3] = p4; pVar[4] = p5;     

  fgets(tmp,300,fp);
  
//    Felder:   {  ro,  v1,  v2,  v3,  B1,  B2,  B3,  p  }
  fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      &rho, &vx, &vy, &vz, &Bx, &By, &Bz, &p );
  
  L1[0] = rho; 
  L1[1] = vx; L1[2] = vy; L1[3] = vz;
  L1[4] = Bx; L1[5] = By; L1[6] = Bz;
  L1[7] = p;
  
  fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      &rho, &vx, &vy, &vz, &Bx, &By, &Bz, &p );
  
  R1[0] = rho; 
  R1[1] = vx; R1[2] = vy; R1[3] = vz;
  R1[4] = Bx; R1[5] = By; R1[6] = Bz;
  R1[7] = p;

  printf( " Data from File:\n" );
  printf(  "  %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
       L1[0],L1[1],L1[2],L1[3],L1[4],L1[5],L1[6],L1[7]);
  printf(  "  %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
       R1[0],R1[1],R1[2],R1[3],R1[4],R1[5],R1[6],R1[7]);
  printf(  "  %8.4f%8.4f%8.4f%8.4f%8.4f\n",
       pVar[0],pVar[1],pVar[2],pVar[3],pVar[4]);
      
  fclose(fp);
};


void Load_Data( char *fname, double *L, double *R, double *pVar )
{
  double rho, vx, vy, vz, Bx, By, Bz, p;
  double p1, p2, p3, p4, p5;
  FILE *fp;

  fp = fopen( fname, "r" );

//    Felder:   {  ro,  v1,  v2,  v3,  B1,  B2,  B3,  p  }
  fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      &rho, &vx, &vy, &vz, &Bx, &By, &Bz, &p );
  
  L[0] = rho; 
  L[1] = vx; L[2] = vy; L[3] = vz;
  L[4] = Bx; L[5] = By; L[6] = Bz;
  L[7] = p;
  
  fscanf( fp, " %lf %lf %lf %lf %lf %lf %lf %lf\n", 
      &rho, &vx, &vy, &vz, &Bx, &By, &Bz, &p );
  
  R[0] = rho; 
  R[1] = vx; R[2] = vy; R[3] = vz;
  R[4] = Bx; R[5] = By; R[6] = Bz;
  R[7] = p;
  
  fscanf( fp, " %lf %lf %lf %lf %lf\n",
      &p1, &p2, &p3, &p4, &p5 );
      
  pVar[0] = p1; pVar[1] = p2; pVar[2] = p3; 
  pVar[3] = p4; pVar[4] = p5;     

  printf( " Data from File:\n" );
  printf(  "  %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
       L[0],L[1],L[2],L[3],L[4],L[5],L[6],L[7]);
  printf(  "  %8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f%8.4f\n",
       R[0],R[1],R[2],R[3],R[4],R[5],R[6],R[7]);
  printf(  "  %8.4f%8.4f%8.4f%8.4f%8.4f\n",
       pVar[0],pVar[1],pVar[2],pVar[3],pVar[4]);
      
  fclose(fp);
};


