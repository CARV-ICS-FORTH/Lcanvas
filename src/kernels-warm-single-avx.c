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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <immintrin.h>
#include <stdint.h>
#include "kerntypes.h"

//IMPLEMENTATION: structs_baseline
double lcanvas_warm_single_structs_baseline(lcanvas_driver_context * restrict lc)
{

  double force = 0;            // this will accumulate the total force on target from all the neigbhours
  const P_t *P = lc->P;
  const int target = lc->target;
  const int *ngb_list = lc->ngb_list; 
  const int Nngb = lc->Nngb;

  for ( int i = 0; i < Nngb; i++ )
    {
      int k = ngb_list[i];    // the index of the i-th neighbour

      // get the distance of the i-th neighbour
      // along all the coordinates
      double dx = P[target].pos[0] - P[k].pos[0];
      double dy = P[target].pos[1] - P[k].pos[1];
      double dz = P[target].pos[2] - P[k].pos[2];

      // square the distance
      //
      double dist2 = dx*dx + dy*dy + dz*dz;

      // accumulate the force
      //
      force += P[target].mass * P[k].mass / dist2;	  
    }

  return force;
}


double lcanvas_warm_single_vectors_baseline(lcanvas_driver_context * restrict lc)
{
  double force     = 0;            // this will accumulate the total force on target from all the neigbhours
  const v4df * restrict V = (const v4df *)lc->Vec;
  const int target = lc->target;
  const int *ngb_list = lc->ngb_list; 
  const int Nngb = lc->Nngb;
  v4df   mask_pos  = {1,1,1,0};    // used to mask out the mass
  v4df   mask_mass = {0,0,0,1};    // used to mask out the position
  
  for ( int i = 0; i < Nngb; i++ )
    {
      int k = ngb_list[i];         // the index of the i-th neighbour

      // the distance^2f the i-th neighbour
      // 
      v4df dist2 = (V[target] - V[k]) * (V[target] - V[k]) * mask_pos;

      // the mass product with the i-th neighbour
      //
      v4df prod  = (V[target] * V[k]) * mask_mass;

      // get the mass product as a single double
      //
      double mm = ((v4df_u)prod).p[3];

      // get the force from the i-th neighbour
      //
      force += mm /( ((v4df_u)dist2).p[0] + ((v4df_u)dist2).p[1] + ((v4df_u)dist2).p[2] );
    }

  return force;
}

double lcanvas_warm_single_vectors_avx(lcanvas_driver_context * restrict lc)
{
  double    force    = 0;      // this will accumulate the total force on target from all the neigbhours
  const v4df * restrict V = (const v4df *)lc->Vec;
  const int target = lc->target;
  const int *ngb_list = lc->ngb_list; 
  const int Nngb = lc->Nngb;
  long long bb       = -1;     // used for masking entries of a vector
  
  __m256d mask_pos   = _mm256_set_pd(0, *(double*)&bb, *(double*)&bb, *(double*)&bb);  // used to mask out the mass
  __m256d mask_mass  = _mm256_set_pd(*(double*)&bb, 0, 0, 0);                          // used to mask out the position
  
  const __m256d register vtarget = _mm256_load_pd(&V[target*4]);                       // load the target data


  for ( int i = 0; i < Nngb; i++ )
    {      
      int k = ngb_list[i];

      // load the i-th neighbour data
      //
      const __m256d register vobj = _mm256_load_pd(&V[k*4]);       
      
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
      
      // get [ dx, dy, ddz dm ]
      vdist = _mm256_sub_pd         ( vtarget, vobj     );
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
      mprod = _mm256_mul_pd( vtarget, vobj );
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

  return force;
}


double lcanvas_warm_single_arrays_baseline(lcanvas_driver_context * restrict lc)
{

  double force = 0;                // this will accumulate the total force on target from all the neigbhours
  const double * restrict pos_x = lc->x_pos;
  const double * restrict pos_y = lc->x_pos;
  const double * restrict pos_z = lc->x_pos;
  const double * restrict mass = lc->mass;
  const int target = lc->target;
  const int *ngb_list = lc->ngb_list; 
  const int Nngb = lc->Nngb;
  double x     = pos_x[target];    // get the target's position and mass
  double y     = pos_y[target];
  double z     = pos_z[target];
  double m     = mass[target];

  int ngb_4 = (Nngb / 4)*4;       // prepare to unroll by 4
  
  for ( int i = 0; i < ngb_4; i+=4 )
    {

      // get distances along x, y and z for 4 neighbours
      //
      double dx[4] = { x - pos_x[ngb_list[i]], x - pos_x[ngb_list[i+1]], x - pos_x[ngb_list[i+2]], x - pos_x[ngb_list[i+3]] };
      double dy[4] = { y - pos_y[ngb_list[i]], y - pos_y[ngb_list[i+1]], y - pos_y[ngb_list[i+2]], y - pos_y[ngb_list[i+3]] };
      double dz[4] = { z - pos_z[ngb_list[i]], z - pos_z[ngb_list[i+1]], z - pos_z[ngb_list[i+2]], z - pos_z[ngb_list[i+3]] };

      // get mass product for 4 neighbours
      //
      double mm[4] = { m  * mass[ngb_list[i]], m  * mass[ngb_list[i+1]], m  * mass[ngb_list[i+2]], m  * mass[ngb_list[i+3]] };

      // calculate the distance^2 for 4 neighbours
      //
      double dist2[4] = { dx[0]*dx[0]+dy[0]*dy[0]+dz[0]*dz[0],
	dx[1]*dx[1]+dy[1]*dy[1]+dz[1]*dz[1],
	dx[2]*dx[2]+dy[2]*dy[2]+dz[2]*dz[2],
	    dx[3]*dx[3]+dy[3]*dy[3]+dz[3]*dz[3] };

      // get the force from 4 neighbrous
      //
      double _force_[4] = { mm[0] / dist2[0],
	mm[1] / dist2[1],
	mm[2] / dist2[2],
	mm[3] / dist2[3] };

      // accumulate on force
      //
      force += (_force_[0] + _force_[1]) + (_force_[2] + _force_[3]);	  
    }

  // process the reminder of the iteration space
  //
  for ( int i = ngb_4; i < Nngb; i++ )
    {
      double dist2 = (
		      (x - pos_x[ngb_list[i]])*(x - pos_x[ngb_list[i]]) +
		      (y - pos_y[ngb_list[i]])*(y - pos_y[ngb_list[i]]) ) +
	(z - pos_z[ngb_list[i]])*(z - pos_z[ngb_list[i]]);
      
      force += mass[target]*mass[ngb_list[i]] / dist2;					      
    }
  
  return force;
}

double lcanvas_warm_single_arrays_avx(lcanvas_driver_context * restrict lc)
{

  double force = 0;           // this will accumulate the total force on target from all the neigbhours  
  const double * restrict pos_x = lc->x_pos;
  const double * restrict pos_y = lc->x_pos;
  const double * restrict pos_z = lc->x_pos;
  const double * restrict mass = lc->mass;
  const int target = lc->target;
  const int *ngb_list = lc->ngb_list; 
  const int Nngb = lc->Nngb;

                              // make an array of 4 vectors with the x,y,z and mass of the target
  __m256d vtarget[4] = {
    _mm256_set_pd( pos_x[target], pos_x[target], pos_x[target], pos_x[target] ),
    _mm256_set_pd( pos_y[target], pos_y[target], pos_y[target], pos_y[target] ),
    _mm256_set_pd( pos_z[target], pos_z[target], pos_z[target], pos_z[target] ),
    _mm256_set_pd( mass [target], mass [target], mass [target], mass [target] )};
  __m256d vforce = _mm256_set_pd( 0, 0, 0, 0 );
  
  int ngb_4 = (Nngb / 4)*4;   // prepare the unroll by 4
  
  for ( int i = 0; i < ngb_4; i+=4 )
    {

      // get the indexes of the next 4 neighbours
      //
      int k[4] = { ngb_list[i], ngb_list[i+1], ngb_list[i+2], ngb_list[i+3] };

      // make an array of 4 vectors with the x,y,z and mass of the next 4 neighbours
      // each row contains either x, y, z, mass; each column refer to a neighbour
      __m256d v[4] = {
	_mm256_set_pd( pos_x[k[3]], pos_x[k[2]], pos_x[k[1]], pos_x[k[0]]),
	_mm256_set_pd( pos_y[k[3]], pos_y[k[2]], pos_y[k[1]], pos_y[k[0]]),
	_mm256_set_pd( pos_z[k[3]], pos_z[k[2]], pos_z[k[1]], pos_z[k[0]]),
	_mm256_set_pd( mass [k[3]], mass [k[2]], mass [k[1]], mass [k[0]]) };

      // calculate dx^2, dy^2, dz^2 for all the 4 neighbours
      //
      for ( int j = 0; j < 3; j++ ) {
	v[j] = _mm256_sub_pd( vtarget[j], v[j] );
	v[j] = _mm256_mul_pd( v[j], v[j] ); }
      // calculate mass product
      //
      v[3] = _mm256_mul_pd( vtarget[3], v[3] );

      //
      // dx^2 + dy^2 for all the 4 neighbours      
      v[0]   = _mm256_add_pd( v[0], v[1] );
      //
      // add dz^2 for all the 4 neighbours
      v[0]   = _mm256_add_pd( v[0], v[2] );
      //
      // get mm / (dx^2 + dy^2 + dz^2) for all 4 neighbours
      v[0]   = _mm256_div_pd( v[3], v[0] );
      //
      // accumulate the force on 4 lanes
      vforce = _mm256_add_pd( vforce, v[0] );
    }

  // process the reminder of the iteration space
  //
  for ( int i = ngb_4; i < Nngb; i++ )
    {
      double dist2 = (
		      (pos_x[target] - pos_x[ngb_list[i]])*(pos_x[target] - pos_x[ngb_list[i]]) +
		      (pos_y[target] - pos_y[ngb_list[i]])*(pos_y[target] - pos_y[ngb_list[i]]) +
		      (pos_z[target] - pos_z[ngb_list[i]])*(pos_z[target] - pos_z[ngb_list[i]]) );
      
      force += mass[target]*mass[ngb_list[i]] / dist2;					      
    }

  // reduce to a unique double the vector force
  // calculated before
  //
  for ( int i = 0; i < 4; i++ )
    force += ((double*)&vforce)[i];
  
  return force;
}

