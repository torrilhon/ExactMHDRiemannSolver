
// some tools for convenience, command line options and time measurements

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

#include "tools.h"


int exists_argument(int argc,const char **argv,const char *keystr)
{
  int i;

  for (i=1;i<argc;i++) if (strstr(argv[i],keystr)) return 1;
  return 0;
}

float get_float_argument(int argc,const char **argv,const char *keystr)
{
  int i;
  float result=0.0;

  for (i=1;i<argc;i++) {
    if (strstr(argv[i],keystr)) {
      result=atof(argv[i+1]);
      break;
    }
  }
  if (i == argc) {
    printf("Error: get_float_argument(%s)\n",keystr);
    return (float) -999999999;
  }
  else return result;
}

int get_string_argument(int argc,const char **argv,const char *keystr,char *result)
{
  int i;

  for (i=1;i<argc;i++) {
    if (strstr(argv[i],keystr)) {
      sprintf(result,"%s",argv[i+1]);
      break;
    }
  }
  if (i == argc) {
    fprintf(stderr,"Error: get_string_argument(%s)\n",keystr);
    return -1;
  }
  else return 0;
}


char *Convert( double a, char *s )
{
  int i = (int)fabs(a+1e-7);
  int j = (int)fabs(10*(a-i+1e-7));
  int k = (int)fabs(100*(a-i-0.1*j+1e-7));
  int l = (int)fabs(1000*(a-i-0.1*j-0.01*k+1e-7));

  if( l != 0 ) sprintf( s, "%dp%d%d%d", i, j, k, l );
  else if( k != 0 ) sprintf( s, "%dp%d%d", i, j, k );
  else sprintf( s, "%dp%d", i, j );
  return( s );
};


/******** functional implementation of time measurements*********************/
#ifndef WIN32

timeval time0;

unsigned int GetmSec()
{
  struct timeval time;

  gettimeofday(&time, NULL );
  
  return( (unsigned int)(1.0e3*(time.tv_sec-time0.tv_sec)+1.0e-3*(time.tv_usec-time0.tv_usec)) );
}
  
void TickTimer()
{
  struct timeval time;

  gettimeofday(&time, NULL );
  
  if( time.tv_usec > time0.tv_usec ) {
    time0.tv_sec = time.tv_sec-time0.tv_sec;
    time0.tv_usec = time.tv_usec-time0.tv_usec;
  } else {
    time0.tv_sec = time.tv_sec-time0.tv_sec-1;
    time0.tv_usec = time0.tv_usec-time.tv_usec;
  };
};
#endif


//LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
/****** dummy implementation of time measurements for windows systems ******/
/******   reason: microsec-times and gettimeofday() are not available ******/
#ifdef WIN32

unsigned int GetmSec()
{ return( 0 ); };

int gettimeofday(struct timeval *tp, void *tzp)
{ return( 0 ); }; 

#endif
//TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

