
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
#include <assert.h>

#include <vector_pragmas.h>
#include <vector_types.h>

#include "mypapi.h"
//#include <immintrin.h>

#define CPU_TIME ({struct timespec ts; clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
					 (double)ts.tv_nsec * 1e-9;})

// with how many neighbours each particle interacts, in avg
//
#define NNEIGHBOURS 200     

#define NREPETITIONS 5

#define DATASIZE 4              // 3 positions + the mass

// how many different implementations
//
#define NSHOTS 5
const char *implementation_labels[NSHOTS+1] = {"generate particles",
  "Using AoS",
  "using AoV (array of vectors)",
  "using AoV with intrinsics",
  "using SoA",
  "using SoA with intrinsics" };

// define a vector type
//

#define DV4SIZE DATASIZE
typedef double v4df __attribute__ ((vector_size (DV4SIZE*sizeof(double))));
typedef union
{
  v4df   P;
  double p[DV4SIZE];
} v4df_u;


// define a structure
typedef struct {
  double pos[3];
  double mass;
 #if defined(FILL)
  double array[FILL];
 #endif
} P_t;




// ···································································
//
/*
 * Processing the particles as structures
 * 
 */

double process_with_structure( const P_t * restrict P, void * restrict mem,
			       const int target, const int *ngb_list, const int Nngb,
			       double *timing)
{

  P_t * restrict memP = __builtin_assume_aligned((P_t*)mem, VALIGN);
  for ( int i = 0; i < Nngb; i++ )
    memP[i] = P[ ngb_list[i] ];
  
  *timing = CPU_TIME;
  
  P_t register Ptarget = P[target];
  double register force = 0;            // this will accumulate the total force on target from all the neigbhours

  IVDEP
  LOOP_VECTORIZE
  VECTOR_ALIGNED
  for ( int i = 0; i < Nngb; i++ )
    {

      // 1D distances
      //
      double register dx = Ptarget.pos[0] - memP[i].pos[0];
      double register dy = Ptarget.pos[1] - memP[i].pos[1];
      double register dz = Ptarget.pos[2] - memP[i].pos[2];

      // square the distance
      //
      double register dist2 = dx*dx + dy*dy + dz*dz;
      
      // accumulate the force
      //
      force += Ptarget.mass * memP[i].mass / dist2;

      /*
      double check = P[target].mass * memP[i].mass / dist2;
      printf("T%d n %d->%d %g %g %g %g\n",
	     target, i, ngb_list[i], dx*dx, dy*dy, dz*dz, check );
      */

    }

  *timing = CPU_TIME - *timing;
  return force;
}


// ···································································
//
/*
 * Here we trat the particles explicitly as vectors
 * The array if structures P is accessed as v4df vectors
 *
 */



double process_with_vectors( const v4df * restrict V,
			           void * restrict mem,
			     const int target,
			     const int *ngb_list,
			     const int Nngb,
			     double *timing )
{
  #define BUNCH  4
  
  alignas(VALIGN) const v4df   Vtarget   = V[target];    // save the target vector in a local vector (much probably a register)
  alignas(VALIGN) const v4df   mask_pos  = {1,1,1,0};    // used to mask out the mass
  alignas(VALIGN) const v4df   mask_mass = {0,0,0,1};    // used to mask out the position

       #if defined(VECTORS_V2)
	double vforce[BUNCH] = {0};
       #endif
	double force     = 0;                // this will accumulate the final total force on target
	                                     // from all the neigbhours

  ASSUME_ALIGNED(mem, VALIGN);
  ASSUME_ALIGNED(V, VALIGN);
  ASSUME_ALIGNED(ngb_list, VALIGN);
  v4df *memV = ASSUME_ALIGNED((v4df*)mem, VALIGN);
  
  IVDEP
  LOOP_UNROLL_N(4)
  VECTOR_ALIGNED
  for ( int i = 0; i < Nngb; i++ )
    memV[i] = V[ ngb_list[i] ];

  *timing = CPU_TIME;
  
  // set-up the new iteration spaces
  //
  int _Nngb_ = Nngb & (-BUNCH);

  IVDEP
  VECTOR_ALIGNED
  LOOP_UNROLL_N(BUNCH)
  for ( int i = 0; i < _Nngb_; i++ )
    {
      // the distance^2 of the i-th neighbour
      //     
      v4df register dist2 = (Vtarget - memV[i  ]) * (Vtarget - memV[i  ]) * mask_pos;
      // the mass product with the i-th neighbour
      //
      v4df register prod  = (Vtarget * memV[i  ]) * mask_mass;
      // get the mass product as a single double
      //      
      double mm = ((v4df_u)prod).p[3];
      // get the force from the i-th neighbour
      //
      force += mm /( ((v4df_u)dist2).p[0] + ((v4df_u)dist2).p[1] + ((v4df_u)dist2).p[2] );
    }
  
  // process the reamining particles
  //
  for ( int i = _Nngb_; i < Nngb; i++ )
    {
      int k = ngb_list[i];
      
      v4df dist2 = (Vtarget - V[k]) * (Vtarget - V[k]) * mask_pos;
      
      // the mass product with the i-th neighbour
      //
      v4df prod  = (Vtarget * V[k]) * mask_mass;
      
      // get the mass product as a single double
      //
      double mm = ((v4df_u)prod).p[3];
      
      // get the force from the i-th neighbour
      //
      force += mm /( ((v4df_u)dist2).p[0] + ((v4df_u)dist2).p[1] + ((v4df_u)dist2).p[2] );
    }

  *timing = CPU_TIME - *timing;
  return force;
 #undef BUNCH
}




// ···································································
//
/*
 *
 */

double process_with_vectors_intrinsics( const double * restrict V,
					      void   * restrict mem,
					const int      target,
					const int    * ngb_list,
					const int      Nngb,
					      double * timing )
  
{

  double * restrict memV = ASSUME_ALIGNED((double*)mem, VALIGN);
  __m256d * restrict memVV = ASSUME_ALIGNED ((__m256d*)mem, VALIGN);

  IVDEP
  LOOP_VECTORIZE
  LOOP_UNROLL
    for ( int i = 0; i < Nngb; i++ ) {
      int off0 = i*4;
      int off1 = ngb_list[i]*4;
      for ( int j = 0; j < 4; j++ )
	memV[off0+j] = V[off1+j]; }

  *timing = CPU_TIME;
  
  const __m256d register vtarget = _mm256_load_pd(&V[target*4]);                       // load the target data    
  long long bb       = -1;                                                             // used for masking entries of a vector
  __m256d register mask_pos   = _mm256_set_pd(0, *(double*)&bb, *(double*)&bb, *(double*)&bb);  // used to mask out the mass
  __m256d register mask_mass  = _mm256_set_pd(*(double*)&bb, 0, 0, 0);                          // used to mask out the position 
  double force = 0;
  
  IVDEP
  LOOP_VECTORIZE
  LOOP_UNROLL
  for ( int i = 0; i < Nngb; i++ )
    {      
      __m256d register vdist;

      // note: the 64bits fields in the following
      //       comments are reported so that at
      //       th beginning there are the first
      //       bits in the vector field, i.e.
      //       those at the least significant bits
      //       >> however, when you consider the
      //       _mm256_set_pd instruction, the arguments
      //       are in the opposite direction: the first
      //       one sets the most significant bits and
      //       so on
      //
      //       least --->  most
      
      // get [ dx, dy, dz dm ]
      vdist = _mm256_sub_pd         ( vtarget, memVV[i] );
      //
      // get [ dx^2, dy^2, dz^2, dm^2 ]
      vdist = _mm256_mul_pd         ( vdist, vdist      );
      
      //
      // get [ dx^2, dy^2, dz^2, 0 ]
      vdist = _mm256_and_pd         ( vdist, mask_pos   );      
      //
      // get [ dx^2+dy^2, dx^2+dy^2, dz^2, dz^2]
      vdist = _mm256_hadd_pd        ( vdist, vdist      );
      //
      // get [ dx^2+dy^2, dz^2, dx^2+dy^2, dz^2 ]
      vdist = _mm256_permute4x64_pd ( vdist, 0b00100111 );
      //
      // get [ dx^2+dy^2+dz^2, .., .., .. ]
      vdist = _mm256_hadd_pd        ( vdist, vdist      );
      
      __m256d mprod;
      //
      // get [ x0*x1, y0*y1, z0*z1, m0*m1 ]
      mprod = _mm256_mul_pd( vtarget, memVV[i] );
      //
      // get [ 0, 0, 0, m0*m1 ]
      mprod = _mm256_and_pd( mprod, mask_mass );      
      //
      // get [ 0, 0, 0, m0*m1 / (dx^2+dy^2+dz^2) ]
      mprod = _mm256_div_pd( mprod, vdist );
      
      // get the uppermost element from the vector
      // --> here we should use the reshuffle+extraction instruction
      force += (((double*)&mprod)[3]);
      
    }

  
  
  *timing = CPU_TIME - *timing;
  return force;
 #undef BUNCH
}



// ···································································
//
/*
 *
 */

double process_with_arrays( const double * restrict pos_x,
			    const double * restrict pos_y,
			    const double * restrict pos_z,
			    const double * restrict mass,
			          void   * restrict mems[],
			    const int               target,
			    const int    * restrict ngb_list,
			    const int               Nngb,
			          double *          timing )
  
{
  #define BUNCH DVSIZE

  ASSUME_ALIGNED(pos_x, VALIGN);
  ASSUME_ALIGNED(pos_y, VALIGN);
  ASSUME_ALIGNED(pos_z, VALIGN);
  ASSUME_ALIGNED( mass, VALIGN );
  ASSUME_ALIGNED( ngb_list, VALIGN );

  double *_memx = ASSUME_ALIGNED((double*)mems[0], VALIGN);
  double *_memy = ASSUME_ALIGNED((double*)mems[1], VALIGN);
  double *_memz = ASSUME_ALIGNED((double*)mems[2], VALIGN);
  double *_memm = ASSUME_ALIGNED((double*)mems[3], VALIGN);

  IVDEP
  LOOP_VECTORIZE
  LOOP_UNROLL_N(DVSIZE)
  for ( int i = 0; i < Nngb; i++ ) {
    int k = ngb_list[i];
    _memx[i] = pos_x[k];
    _memy[i] = pos_y[k];
    _memz[i] = pos_z[k];
    _memm[i] = mass [k]; }
   

  *timing = CPU_TIME;
  
  const double x[BUNCH] = { [0 ... BUNCH-1] = pos_x[target]};
  const double y[BUNCH] = { [0 ... BUNCH-1] = pos_y[target]};
  const double z[BUNCH] = { [0 ... BUNCH-1] = pos_z[target]};
  const double m[BUNCH] = { [0 ... BUNCH-1] = mass [target]};

  double vforce[BUNCH] ATTRIBUTE_ALIGNED((VALIGN))  = { [0 ... BUNCH-1] = 0 };
  double force = 0;                // this will accumulate the total force on target from all the neigbhours

  int _Nngb_ = Nngb & (-(int)BUNCH);
  
  IVDEP
  LOOP_VECTORIZE
  LOOP_UNROLL_N(1)
  for ( int i = 0; i < _Nngb_; i+=BUNCH )
    {

      int register off= i;
      
      // get distances along x, y and z for 4 neighbours
      //
      double dx[BUNCH] ATTRIBUTE_ALIGNED((VALIGN));
      double dy[BUNCH] ATTRIBUTE_ALIGNED((VALIGN));
      double dz[BUNCH] ATTRIBUTE_ALIGNED((VALIGN));
      double mm[BUNCH] ATTRIBUTE_ALIGNED((VALIGN));
      double dist2[BUNCH] ATTRIBUTE_ALIGNED((VALIGN));

      IVDEP
      LOOP_VECTORIZE
      LOOP_UNROLL_N(BUNCH)
      for ( int j = 0; j < BUNCH; j++ ) {
	dx[j] = x[j] - _memx[off+j];
	dy[j] = y[j] - _memy[off+j];
	dz[j] = z[j] - _memz[off+j];
	mm[j] = m[j] * _memm[off+j]; }
      
      // calculate the distance^2 and the forcefor BUNCH neighbours
      //
      IVDEP
      LOOP_UNROLL_N(BUNCH)
      for ( int j = 0; j < BUNCH; j++ )
	{
	  // distance
	  dist2[j] = dx[j]*dx[j]+dy[j]*dy[j]+dz[j]*dz[j];	  
	  // get the force 
	  vforce[j] += mm[j] / dist2[j];
	}
            
    }

  // accumulate on force
  //
  force += (vforce[0] + vforce[1]) + (vforce[2] + vforce[3]);	  
  
  // process the reminder of the iteration space
  //
  for ( int i = _Nngb_; i < Nngb; i++ )
    {
      int k = ngb_list[i];
      double dist2 = (x[0] - pos_x[k])*(x[0] - pos_x[k]) +
	(y[0] - pos_y[k])*(y[0] - pos_y[k]) +
	(z[0] - pos_z[k])*(z[0] - pos_z[k]);
      
      force += mass[target]*mass[ngb_list[i]] / dist2;					      
    }

  *timing = CPU_TIME - *timing;
  return force;
 #undef BUNCH
}



// ···································································
//
/*
 *
 */

double process_with_arrays_intrinsics( const double * restrict pos_x,
				       const double * restrict pos_y,
				       const double * restrict pos_z,
				       const double * restrict mass,
				             void   *          mems[],
				       const int               target,
				       const int    * restrict ngb_list,
				       const int               Nngb,
				             double *          timing)
{

  *timing = CPU_TIME;
  
  ASSUME_ALIGNED(pos_x, VALIGN);
  ASSUME_ALIGNED(pos_y, VALIGN);
  ASSUME_ALIGNED(pos_z, VALIGN);
  ASSUME_ALIGNED( mass, VALIGN );
  ASSUME_ALIGNED( ngb_list, VALIGN );
  
  double *_memx = ASSUME_ALIGNED((double*)mems[0], VALIGN);
  double *_memy = ASSUME_ALIGNED((double*)mems[1], VALIGN);
  double *_memz = ASSUME_ALIGNED((double*)mems[2], VALIGN);
  double *_memm = ASSUME_ALIGNED((double*)mems[3], VALIGN);

  IVDEP
  LOOP_VECTORIZE
  LOOP_UNROLL_N(DVSIZE)
  for ( int i = 0; i < Nngb; i++ ) {
    int k = ngb_list[i];
    _memx[i] = pos_x[k];
    _memy[i] = pos_y[k];
    _memz[i] = pos_z[k];
    _memm[i] = mass [k]; }

  dvector_t *_vmemx = ASSUME_ALIGNED((dvector_t*)mems[0], VALIGN);
  dvector_t *_vmemy = ASSUME_ALIGNED((dvector_t*)mems[1], VALIGN);
  dvector_t *_vmemz = ASSUME_ALIGNED((dvector_t*)mems[2], VALIGN);
  dvector_t *_vmemm = ASSUME_ALIGNED((dvector_t*)mems[3], VALIGN);

 #if defined(__clang__) || defined(__INTEL_LLVM_COMPILER)
  dvector_t register vtargetx ATTRIBUTE_ALIGNED(VALIGN) = (dvector_t)(pos_x[target]);
  dvector_t register vtargety ATTRIBUTE_ALIGNED(VALIGN) = (dvector_t)(pos_y[target]);
  dvector_t register vtargetz ATTRIBUTE_ALIGNED(VALIGN) = (dvector_t)(pos_z[target]);
  dvector_t register vtargetm ATTRIBUTE_ALIGNED(VALIGN) = (dvector_t)(mass[target]);
  dvector_t          one      ATTRIBUTE_ALIGNED(VALIGN) = (dvector_t)(1.0);
 #else
  dvector_t register vtargetx ATTRIBUTE_ALIGNED(VALIGN) = {pos_x[target]};
  dvector_t register vtargety ATTRIBUTE_ALIGNED(VALIGN) = {pos_y[target]};
  dvector_t register vtargetz ATTRIBUTE_ALIGNED(VALIGN) = {pos_z[target]};
  dvector_t register vtargetm ATTRIBUTE_ALIGNED(VALIGN) = {mass[target]};
  dvector_t          one      ATTRIBUTE_ALIGNED(VALIGN) = {1.0};
 #endif
  dvector_t register vforce   ATTRIBUTE_ALIGNED(VALIGN) = {0};
  double             force = 0;           // this will accumulate the total force on target from all the neigbhours

  int _Nngb_ = Nngb & (-DVSIZE);
  _Nngb_ = _Nngb_ / (DVSIZE);
  
  IVDEP
  LOOP_VECTORIZE
  for ( int i = 0; i < _Nngb_; i++ )
    {
      
      dvector_t register dist2x = (vtargetx - _vmemx[i]);      
      dvector_t register dist2y = (vtargety - _vmemy[i]);      
      dvector_t register dist2z = (vtargetz - _vmemz[i]);

      dist2x *= dist2x;
      dist2y *= dist2y;
      dist2z *= dist2z;
      
      dvector_t register inv_dist = one / (dist2x+dist2y+dist2z);
      vforce += (vtargetm * _vmemm[i]) * inv_dist;
      
    }

  // process the reminder of the iteration space
  //
  for ( int i = _Nngb_*(DVSIZE); i < Nngb; i++ )
    {
      double dist2 = (
		      (pos_x[target] - _memx[i])* (pos_x[target] - _memx[i]) +
		      (pos_y[target] - _memy[i])* (pos_y[target] - _memy[i]) +
		      (pos_z[target] - _memz[i])* (pos_z[target] - _memz[i]) );
      
      force += mass[target]*_memm[i] / dist2;		      
    }

  // reduce to a unique double the vector force
  // calculated before
  //
  dvector_u _vforce_;
  _vforce_.V = vforce;
  for ( int i = 0; i < 4; i++ )
    force += _vforce_.v[i];

  *timing = CPU_TIME - *timing;

  return force;
  
 #undef BUNCH
}


int main( int argc, char **argv )
{

           int case_to_run = (argc > 1 ? atoi(*(argv+1)) : 0 );
  unsigned int N           = (argc > 2 ? atoi(*(argv+2)) : 1000000 );
  long     int seed        = (argc > 3 ? atoi(*(argv+3)) : -1 );     // set the seed for repeatible runs
           int Nrepetitions= (argc > 4 ? atoi(*(argv+4)) : NREPETITIONS );
           int dry_run     = (argc > 5 ? atoi(*(argv+5)) : 0 );      // 1 to estimate floats for initialization
           int from_file   = ( case_to_run < 0 );	   
	   
  case_to_run = (case_to_run < 0 ? ~case_to_run + 1 : case_to_run);  // make it positive
	   
  if ( (case_to_run < 0) || (case_to_run > 5) ) {
    printf("unknown case %d\n", case_to_run );
    return 0; }
  

  P_t    * restrict P;
  v4df   * restrict vpos;
  double * restrict pos_x;
  double * restrict pos_y;
  double * restrict pos_z;
  double * restrict mass ;               
  double * restrict force; 
  int    * restrict ngb_list;
  
  FILE   * outfile = NULL;
  FILE   * infile  = NULL;
  

  if ( case_to_run == 0 )
    {
      outfile = fopen( "particle.data", "w" );
      fwrite( &N, sizeof(int), 1, outfile );
      fwrite( &Nrepetitions, sizeof(int), 1, outfile );
    }
  else 
    {
      if ( from_file ) {
	size_t ret;
	int    Nrep;
	infile = fopen( "particle.data", "r" );
	ret = fread( &N, sizeof(int), 1, infile );
	ret = fread( &Nrep, sizeof(int), 1, infile );
	if ( Nrep < Nrepetitions ) {
	  printf("Nrepetitions is %d, larger than that "
		 "used to generate the particle file, which was %d.\n"
		 "This may lead to inf or NaN, better stop.\n",
		 Nrepetitions, Nrep);
	  exit(3); }	  
      }

      force = (double *)aligned_alloc(VALIGN, N*sizeof(double));
    }


  printf(" >  Running case %d -> \"%s\" with %u points and %d repetitions\n"
	 " >  Vectors are %lu bytes long, %lu double long\n",
	 case_to_run, implementation_labels[case_to_run], N, Nrepetitions,
	 VSIZE, DVSIZE);
  if ( dry_run )
    printf("  >>> DRY RUN :: use to estimate ops outside the actual calculations\n" );

  
  if ( case_to_run > 0 )
    {
      switch( case_to_run )
	{
	case 1: {
	  // allocate array of big structures
	  //
	  P = (P_t    *)aligned_alloc(VALIGN,   N * sizeof(P_t)); } break;
	  
	case 2:
	case 3: {      
	  // allocate array of small structures (vectors)
	  //
	  vpos  = (v4df   *)aligned_alloc(VALIGN,   N * sizeof(v4df)); } break;
	  
	case 4:
	case 5: {
	  // allocate arrays
	  //
	  pos_x  = (double *)aligned_alloc(VALIGN, N*sizeof(double));
	  pos_y  = (double *)aligned_alloc(VALIGN, N*sizeof(double));
	  pos_z  = (double *)aligned_alloc(VALIGN, N*sizeof(double));
	  mass   = (double *)aligned_alloc(VALIGN, N*sizeof(double));
	}break;
	}

      memset( (void*)force, 0, N*sizeof(double));
    }
  
  // initialize particle data
  //
    {
      if ( ! from_file ) {
	if ( seed == -1 )
	  srand48(time(NULL));
	else
	  srand48( seed ); }

      for ( int i = 0; i < N; i++ )
	{
	  double data[4] __attribute__((aligned(VALIGN))); 
	  if ( ! from_file )
	    for ( int j = 0; j < 4; j++ )
	      data[j] = drand48();
	  else
	    {
	      size_t ret = fread ( data, 4, sizeof(double), infile );
	    }

	  if ( case_to_run >= 0 )	    
	    switch( case_to_run )
	      {
	      case 0:{
		size_t ret = fwrite( data, sizeof(double), 4, outfile ); } break;
		
	      case 1: {		
		P[i].pos[0] = data[0];
		P[i].pos[1] = data[1];
		P[i].pos[2] = data[2];
		P[i].mass   = data[3]; } break;
		
	      case 2:
	      case 3: {
		v4df_u * _vpos_ = (v4df_u*)vpos;
		for ( int j = 0; j < 4; j++ )
		  _vpos_[i].p[j] = data[j]; } break;
		
	      case 4:
	      case 5:{ 
		pos_x[i] = data[0];
		pos_y[i] = data[1];
		pos_z[i] = data[2];
		mass [i] = data[3]; } break;
	      }
	}
    }

  if ( dry_run ) {
    printf("dry run: going to clean up\n");
    goto clean; }
  
  /* ------------------------------------------------------
   *
   * simulate interactions
   *
   * ------------------------------------------------------ */

  // determine how many particles will be active, from
  // a minimum of 1/100-ed to a maximum of 1/100+1/10
  //

  int Nactive;
  int AvgNneighbours;
  
  if ( case_to_run >= 0 && !from_file )
    {
      Nactive  = N/100 + (int)lrand48() % (N/10);  // 11% of particle will be processed at max, 1% as minimum
      if ( case_to_run == 0 ) {
	size_t ret = fwrite( &Nactive, sizeof(int), 1, outfile );}
    }
  else {
    size_t ret = fread( &Nactive, sizeof(int), 1, infile ); }

  AvgNneighbours = ( NNEIGHBOURS > N*0.1 ? (int)(N*0.1) : NNEIGHBOURS );
  ngb_list = (int*)calloc( AvgNneighbours * 1.2, sizeof(int) );

  PAPI_INIT;
  
  printf ( " *  %d particles will be active\n", Nactive );

  double timing   = 0;
  double intiming = 0;
  
  void *mem = NULL;
  void *mems[4] = {NULL};
  switch ( case_to_run )
    {
    case 1: mem = aligned_alloc(VALIGN, AvgNneighbours*1.2 * sizeof(P_t)); break;
    case 2: 
    case 3: mem = aligned_alloc(VALIGN, AvgNneighbours*1.2 * sizeof(v4df)); break;
    case 4:
    case 5: {
      mems[0] = aligned_alloc(VALIGN, AvgNneighbours*1.2 * sizeof(double));
      mems[1] = aligned_alloc(VALIGN, AvgNneighbours*1.2 * sizeof(double));
      mems[2] = aligned_alloc(VALIGN, AvgNneighbours*1.2 * sizeof(double));
      mems[3] = aligned_alloc(VALIGN, AvgNneighbours*1.2 * sizeof(double)); } break;
    }
      
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

      if ( !from_file )
	{
	  target      = (int)lrand48() % N;                              // the index of the to-be-processed particle
	  
	  // determine how many will be the neighbours
	  // ( same comment than before apply)
	  Nneighbours = AvgNneighbours + (int)mrand48()%(AvgNneighbours/10) ;  // with how many particles it will interact
	  //   Nneighbours +- 10%
	  
	  
	  // decide which will be the neighbours
	  int upperlimit = target+Nrepetitions;
	  for ( int i = 0; i < Nneighbours; i++ )
	    {
	      int draw = (int)lrand48() % N;
	      if ( (draw >= target) && (draw <= upperlimit ) )
		draw = upperlimit+1;
	      ngb_list[i] = draw;
	    }

	  if ( case_to_run == 0 )
	    {
	      printf("iter %d target %d has %u neighbours\n", j,
		     target, (unsigned)Nneighbours );
	      size_t ret;
	      ret = fwrite( &target, sizeof(int), 1, outfile );
	      ret = fwrite( &Nneighbours, sizeof(int), 1, outfile );
	      ret = fwrite( &ngb_list[0], sizeof(int), Nneighbours, outfile );
	    }
	  
	}
      else 
	{
	  size_t ret;
	  ret = fread( &target, sizeof(int), 1, infile );
	  ret = fread( &Nneighbours, sizeof(int), 1, infile );
	  ret = fread( ngb_list, sizeof(int), Nneighbours, infile );
	}

      double myforce;
      double now;
      
      // now, call the force evaluation with the different implementations
      //
      
      if ( case_to_run > 0 ) {

	double this_timing       = 1e10;
	double this_inner_timing = 1e10;
	for ( int shot = 0; shot <= Nrepetitions; shot++ ) {

	  // shot = = iteration is used as a warm-up for the neighbours walk
	  int mytarget = (target+shot)%N;
	 
	  if ( shot >= 1 )
	    now = CPU_TIME;
	  if ( shot == 1 )
	    PAPI_START_CNTR;

	  double inner_timing;
	  switch ( case_to_run )
	    {	  
	      // ------------------------------------------------
	      //
	    case 1: {
	      myforce += process_with_structure( P, mem, mytarget, ngb_list, Nneighbours, &inner_timing );
	    } break;
	    
	      // ------------------------------------------------
	      //     
	    case 2: {
	      myforce += process_with_vectors( vpos, mem, mytarget, ngb_list, Nneighbours, &inner_timing );
	    } break;
	    
	      // ------------------------------------------------
	      //
	    case 3: {
	      myforce += process_with_vectors_intrinsics( (double*)vpos, mem, mytarget, ngb_list, Nneighbours, &inner_timing );
	    } break;
	    
	      // ------------------------------------------------
	      //
	    case 4: {
	      myforce += process_with_arrays( pos_x, pos_y, pos_z, mass, mems,
					      mytarget, ngb_list, Nneighbours, &inner_timing );
	    } break;

	      // ------------------------------------------------
	      //
	    
	    case 5: {
	      myforce +=  process_with_arrays_intrinsics( pos_x, pos_y, pos_z, mass, mems,
							  mytarget, ngb_list, Nneighbours, &inner_timing );
	    } break; 
	    }

	  if (shot >= 1 ) {
	    double chrono = CPU_TIME - now;
	    if ( chrono < this_timing )
	      this_timing = chrono;
	    if ( inner_timing < this_inner_timing )
	      this_inner_timing = inner_timing;
	  }

	  if ( shot == 1 )
	    PAPI_STOP_CNTR;	    

	  /* if ( j == 0 ) */
	  /*   printf("shot %d timing is %g\n", shot, CPU_TIME-now); */
	  
	}

	timing += this_timing;
	intiming += this_inner_timing;
      }

      if ( case_to_run > 0 ) {
       #if !defined(NDEBUG)
       #warning "you're compiling with assert"
       #endif
	assert( myforce != 0 );
	assert( !isnan(myforce) );
	assert( !isinf(myforce) );
	force[target] = myforce / (Nrepetitions+1); }
    }

  if ( infile != NULL )
    fclose(infile);
  if ( outfile != NULL )
    fclose(outfile);

  if ( case_to_run > 0 )
    {
      
      printf(" +  timing and inner_timing for case %d \"%s\" : %g %g\n",
	     case_to_run, implementation_labels[case_to_run],
	     timing, intiming );
      
      // just to trick the compiler, so that it does not optimize out
      // pointless calculations
      //

      char  filename[100];
      sprintf( filename, "force.%d.out", case_to_run );
      FILE *output = fopen( filename, "w" );
      if( output != NULL ) {
	fwrite( &N, sizeof(int), 1, output);
	fwrite( force, sizeof(double), N, output );
	fclose(output); }
      else
	printf(">>> wow, I was unable to create stupid file\n");

      free ( force    );
      if ( !dry_run)
	free ( ngb_list );
  
    }

 clean:
  if ( mem != NULL )
    free ( mem );
  if ( mems[0] != NULL )
    for ( int i = 0; i < 4; i++ )
      free ( mems[i] );
  switch ( case_to_run )
    {
    case 1: {
      free ( P ); } break;
    case 2:
    case 3: {
      free ( vpos ); } break;
    case 4:
    case 5: {
      free ( pos_x ), free ( pos_y ), free ( pos_z ), free ( mass ); } break;
    }

 #if defined(USE_PAPI)
  if ( case_to_run > 0 )
    for ( int i = 0; i < PAPI_EVENTS_NUM; i++ )
      printf("PAPI event %d: %llu\n",
	     i, (unsigned long long)papi_values[i]);
 #endif

  
  
  return 0;
}
