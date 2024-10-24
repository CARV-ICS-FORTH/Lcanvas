#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>


int main ( int argc, char **argv )
{
  if ( argc == 1 ) {
    printf ( "I expected the names of the tw files to be comparend" );
    exit(1); }

  FILE  *files[2];
  double *P[2];
  int     N[2];
  double  tolerance;

  files[0]  = fopen( *(argv+1), "r" );
  files[1]  = fopen( *(argv+2), "r" );
  tolerance = (argc>3 ? atof(*(argv+3)) : 1e-5 );

  if ( (files[0] == NULL) || (files[1] == NULL) ) {
    for ( int i = 0; i < 2; i++ )
      if ( files[i] == NULL )
	printf ( "something went wrong in opening the file %s\n", *(argv+i+1) );
    exit(2); }
  
  for ( int i = 0; i < 2; i++ ) 
    fread( &N[i], sizeof(int), 1, files[i] );
  
  if ( N[0] != N[1] )
    printf("file %s contains %d particles, whereas file %s contains %d particles, "
	   "which makes them different\n", *(argv+1), N[0], *(argv+2), N[1] );
  
  P[0] = (double*) malloc( N[0]*2*sizeof(double) );
  P[1] = P[0] + N[0];

  fread( P[0], sizeof(double), N[0], files[0] );
  fread( P[1], sizeof(double), N[1], files[1] );
  
  int    faults   = 0;
  int    processed = 0;
  double min_diff = 1e10, max_diff = 0, avg_diff = 0;
  
  for ( int i = 0; i < N[0]; i++ )
    {
      if ( P[0][i] == P[1][i] )
	continue;

      processed++;
      
      double diff = fabs(P[0][i] - P[1][i]);
      if ( P[0][i] != 0 )
	diff /= P[0][i];

      if (isnan(diff ) )
	printf(".");
      faults += ( diff > tolerance );

      min_diff  = ( diff < min_diff ? diff : min_diff );
      max_diff  = ( diff > max_diff ? diff : max_diff );
      avg_diff += diff;
    }

  if ( faults )
    printf ( "   »»» %d particles differ more tha the tolerance (%g)\n",
	     faults, tolerance );
  else
    printf("   :) OK !!\n");
  
  if ( processed > 0 )
    printf ( "min, avg, max difference of %d compared particles is %g, %g, %g\n", processed, 
	     min_diff, avg_diff / processed, max_diff );
  else
    printf ( "no particles pair was any different\n" );

  free( P[0] );
  fclose( files[1] );
  fclose( files[0] );
  
  return 0;
}
