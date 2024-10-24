
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
 * This file defines common marcros 
 */

#ifndef __LCANVAS_UTILMACROS__
#define __LCANVAS_UTILMACROS__
#include <time.h>

#define STRINGIZE(X) INNER_STRINGIZE(X)
#define INNER_STRINGIZE(X) #X

//if not sure, overestimate
#define CACHE_LINE_BYTES 64 


#define CPU_TIME_DSEC ({struct timespec ts; clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ),\
        (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;})


#define CPU_TIME_LLNSEC (long long)({struct timespec ts; clock_gettime( CLOCK_PROCESS_CPUTIME_ID,\
            &ts ), (long long)(ts.tv_sec * 1000000000ULL) + (long long)ts.tv_nsec;})

#endif /* __LCANVAS_UTILMACROS__ */
