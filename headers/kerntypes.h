
/* ────────────────────────────────────────────────────────────────────────── *
 │                                                                            │
 │ This file is part of the nbody-vectorization-toy-lab in the from of the    │
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

/*
 * Defines the types of the kernels and type specific parameters
 * a type may have many implementations but the types interaction with
 * the driver shoould allways be the same (expect same parameters etc)
 */
#ifndef __LCANVAS_KERNTYPES_
#define __LCANVAS_KERNTYPES__

//perhaps we need to make this more general at some point
#include "dtypes4.h"

#define COLDSINGLE 1
#define WARMSINGLE 2


#if TYPE == COLDSINGLE

#define TYPENAME "cold-single"
#define TYPEPREFIX lcanvas_cold_single
typedef struct lcanvas_driver_context
{
    P_t *P;
    double* Vec;
    double* x_pos;
    double* y_pos;
    double* z_pos;
    double* mass;
    int target;
    int *ngb_list;
    int Nngb;
}lcanvas_driver_context;

#elif TYPE == WARMSINGLE

#define TYPENAME "warm-single"
#define TYPEPREFIX lcanvas_warm_single
typedef struct lcanvas_driver_context
{
    P_t *P;
    double* Vec;
    double* x_pos;
    double* y_pos;
    double* z_pos;
    double* mass;
    int target;
    int *ngb_list;
    int Nngb;
}lcanvas_driver_context;

#else

#error("unsupported type")

#endif



#endif /*__LCANVAS_KERNTYPES__ */
