/* ────────────────────────────────────────────────────────────────────────── *
 │                                                                            │
 │ This file is part of the nbody-vectorization-toy-lab in the fram of the    │
 │ SPACE Center of Excellence                                                 │
 │ © 2024                                                                     │
 │ The work has been started by Luca Tornatore @ INAF (Italian National       │
 │ Institute for Astrophysics) and continued in collaboration with            │
 │ Vassilis Flouris @ FORTH (Foundation for Research and Technology - Hellas  │
 │                                                                            │
 │ Consult the accopanying file changes.log for details                       │
 │                                                                            │
 │ contact:                                                                   │
 │  luca.tornatore@inaf.it                                                    │
 │  vflouris@ics.forth.gr                                                     │
 │                                                                            │
 │                                                                            │
 │ ×××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××  │
 │                                                                            │
 │     This is open source software; you can redistribute it and/or modify    │
 │     it under the terms of the "3 clauses BSD license"                      │
 │       https://opensource.org/license/bsd-3-clause                          │
 │                                                                            │
 │     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    │
 │     “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      │
 │     LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS      │
 │     FOR A PARTICULAR PURPOSE ARE DISCLAIMED.                               │
 │                                                                            │
 * ────────────────────────────────────────────────────────────────────────── */


//starting point vect1

//#if defined(__STDC__)
//#  if (__STDC_VERSION__ >= 199901L)
//#     define _XOPEN_SOURCE 700
//#  endif
//#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include "defaults.h"
#include "util_macros.h"
#include "kerntypes.h"

//for now this is the only supported DATASIZE
//we might want to abstract this away later
//but for now it will have to do
#define DRIVER_PDATASIZE 4

#if defined(DRIVER_PDATASIZE) &&  DRIVER_PDATASIZE == 4

#include "dtypes4.h"

#else

#error("unsupported pdata size\n")

#endif



// kernel declarations
#if TYPE == COLDSINGLE

double FUNCTION_NAME(TYPEPREFIX,IMPLEMENTATION)(struct lcanvas_driver_context * restrict lc);

#elif TYPE == WARMSINGLE

double FUNCTION_NAME(TYPEPREFIX,IMPLEMENTATION)(struct lcanvas_driver_context * restrict lc);

#else

#endif


int main( int argc, char **argv )
{

  unsigned int N;       
  long     int seed;
           int Nrepetitions;
  P_t    * restrict P;
  v4df   * restrict vpos;
  double * restrict pos_x;
  double * restrict pos_y;
  double * restrict pos_z;
  double * restrict mass ;               
  double * restrict force; 
  int    * restrict ngb_list;
  int Nactive;
  int AvgNneighbours;
  double average_time0 = 0.0;
  double average_time1 = 0.0;
  double average_time2 = 0.0;
  double average_time3 = 0.0;
  
  FILE   * outfile = NULL;
  FILE   * infile  = NULL;
  lcanvas_driver_context ld_context;

  N           = (argc > 1 ? atoi(*(argv+1)) : NPART );
  seed        = (argc > 2 ? atoi(*(argv+2)) : -1 );     // set the seed for repeatible runs
  Nrepetitions= (argc > 3 ? atoi(*(argv+3)) : NREPETITIONS );
  
  force = (double *)aligned_alloc(getpagesize(), N*sizeof(double));
  memset( (void*)force, 0, N*sizeof(double)); //maybe change this a bit

  printf("> Running case: \"%s\" with %u points and %d repetitions\n", "" TYPENAME "::" STRINGIZE (IMPLEMENTATION), N, Nrepetitions);

  
  // allocate array of big structures
  //
  P = (P_t    *)aligned_alloc(getpagesize(),   N * sizeof(P_t));
  
  // allocate array of small structures (vectors)
  //
  vpos  = (v4df   *)aligned_alloc(getpagesize(),   N * sizeof(v4df));
  // allocate arrays
  //
  pos_x  = (double *)aligned_alloc(getpagesize(), N*sizeof(double));
  pos_y  = (double *)aligned_alloc(getpagesize(), N*sizeof(double));
  pos_z  = (double *)aligned_alloc(getpagesize(), N*sizeof(double));
  mass   = (double *)aligned_alloc(getpagesize(), N*sizeof(double));

  
  if ( seed == -1 )
      srand48(time(NULL));
  else
      srand48( seed );

  // determine how many particles will be active, from
  // a minimum of 1/100-ed to a maximum of 1/100+1/10
  //
  Nactive  = N/100 + (int)lrand48() % (N/10);  // 11% of particle will be processed at max, 1% as minimum
  AvgNneighbours = ( NNEIGHBOURS > N*0.1 ? (int)(N*0.1) : NNEIGHBOURS );
  ngb_list = (int*)calloc( AvgNneighbours * 1.2, sizeof(int) );
  
  for ( int i = 0; i < N; i++ )
  {
	  double data[DATASIZE] __attribute__((aligned(CACHE_LINE_BYTES))); 
	  for ( int j = 0; j < DATASIZE; j++ )
      {
          data[j] = 0.01 + fabs(drand48()); // dividing by zero-ish is bad, mkay?
      }
#if DATASIZE == 4      
      P[i].pos[0] = data[0];
      P[i].pos[1] = data[1];
      P[i].pos[2] = data[2];
      P[i].mass   = data[3];
	  v4df_u * _vpos_ = (v4df_u*)vpos;
	  for ( int j = 0; j < DATASIZE; j++ )
      {
		  _vpos_[i].p[j] = data[j];
      }
      pos_x[i] = data[0];
      pos_y[i] = data[1];
      pos_z[i] = data[2];
      mass [i] = data[3];
  }
#else
#error("unsupported DATASIZE")
#endif


  /* ------------------------------------------------------
   *
   * populate context
   *
   * ------------------------------------------------------ */

  //common fields

  ld_context.P = P;
  ld_context.Vec = (double *)vpos;
  ld_context.x_pos = pos_x;
  ld_context.y_pos = pos_y;
  ld_context.z_pos = pos_z;
  ld_context.mass =  mass;


  //type specific fields

  #if TYPE == COLDSINGLE

  #elif TYPE == WARMSINGLE

  #else

  #endif
  
  /* ------------------------------------------------------
   *
   * simulate interactions
   *
   * ------------------------------------------------------ */


  
  printf ( "> %d particles will be active\n", Nactive );

  // loop over active particles
  //
  for ( int j = 0; j < Nactive; j++ )
  {
      // determine which one will be active
      // of course it is not important which one will be,
      // nor is important that we do not pick the same
      // particle more than once in the loop
      // (shouldn't happen, though, because the
      // random repetition lenght is much larger)
      //
      int target, Nneighbours;
      double myforce;
      double now;
      double time = 0.0;

	  
	  // determine how many will be the neighbours
	  // ( same comment than before apply)
	  Nneighbours = AvgNneighbours + (int)mrand48()%(AvgNneighbours/10) ;  // with how many particles it will interact
	  //   Nneighbours +- 10%
	  
      target      = (int)lrand48() % N; // the index of the to-be-processed particle
	  // decide which will be the neighbours
	  for ( int i = 0; i < Nneighbours; i++ )
	  {
          int draw = (int)lrand48() % N;
          ngb_list[i] = draw;
	  }
      ld_context.target = target;
      ld_context.ngb_list = ngb_list;
      ld_context.Nngb = Nneighbours;

#if TYPE == COLDSINGLE 

      time= 0.0;
      now = CPU_TIME_DSEC; 
	  myforce =  FUNCTION_NAME(TYPEPREFIX,IMPLEMENTATION) (&ld_context);
      time = CPU_TIME_DSEC - now;
      force[target] = myforce / (Nrepetitions+1);
      average_time0 += time; //dividing once (maybe) reduces error so for now this is sum not avg
	}
    average_time0 = average_time0 / Nactive;
    printf("> average time: \"%s\"  is  %lf sec\n", ""TYPENAME"::"STRINGIZE(IMPLEMENTATION), average_time0);
#elif TYPE == WARMSINGLE

      time= 0.0;
      now = CPU_TIME_DSEC;

	  myforce =  FUNCTION_NAME(TYPEPREFIX,IMPLEMENTATION) (&ld_context);
      
      time = CPU_TIME_DSEC - now;
      average_time0 += time; //dividing once (maybe) reduces error so for now this is sum not avg

      time= 0.0;
      now = CPU_TIME_DSEC;

	  myforce +=  FUNCTION_NAME(TYPEPREFIX,IMPLEMENTATION) (&ld_context);
      
      time = CPU_TIME_DSEC - now;
      average_time1 += time; //dividing once (maybe) reduces error so for now this is sum not avg
      force[target] = myforce / 2;
	}
    average_time0 = average_time0 / Nactive;
    average_time1 = average_time1 / Nactive;
    printf("> average time: \"%s\"  is  cold:%lf hot:%lf sec\n", ""TYPENAME"::"STRINGIZE(IMPLEMENTATION), 
            average_time0, average_time1);
#else

#error("type is specified in header file but not implemented in the driver")

#endif //TYPE

 
      // just to trick the compiler, so that it does not optimize out
      // pointless calculations
      //
    FILE *output = fopen( FORCE_OUTPUT_FILE_NAME, "w" );
    if( output != NULL ) {
        fwrite( &N, sizeof(int), 1, output);
        fwrite( force, sizeof(double), N, output );
        fclose(output); 
    }
    else
    {
        printf(">>> wow, I was unable to create stupid file\n");
        exit(-1);
    }
    free ( force    );
    free ( ngb_list );
    free ( P );
    free ( vpos );
    free ( pos_x ), free ( pos_y ), free ( pos_z ), free ( mass );
    return 0;
}

