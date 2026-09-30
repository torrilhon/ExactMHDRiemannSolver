
// some tools for convenience, command line options and time measurements

#include <unistd.h>
#include <string.h>
#include <sys/time.h>


/*
inline double max( double a, double b ) 
{ 
  return( a > b ? a : b );
}

inline double min( double a, double b ) 
{ 
  return( a < b ? a : b );
}

inline double sgn( double a ) 
{ 
  if( a == 0 ) return( 0 );
  return( a < 0 ? -1 : 1 );
}

inline double abs( double a ) 
{ 
  return( a < 0 ? -a : a );
}
*/

int exists_argument(int argc,const char **argv,const char *keystr);
float get_float_argument(int argc,const char **argv,const char *keystr);
int get_string_argument(int argc,const char **argv,const char *keystr,char *result);
