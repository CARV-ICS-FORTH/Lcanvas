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

/* @ ······················································· @
   :  Here we collect some macros to express vectorization   :
   :  Pragmas indpendently on the compiler used.             :
   :  This is intended to have only didactial purposes and   :
   :  by no means conveys an exhaustive list and/or          :
   :  combination, nor the maximum tuning for every compiler :
   @ ······················································· @
 */


// the two following macros are useful only to
// produce pragma strings in the source files
//
#define STRINGIFY(X) #X
#define _DO_PRAGMA(x) _Pragma (#x)
//


/* ============================================================================
 * ==                                                                        ==
 * ==                                                                        ==
 * == D I C T I O N A R Y                                                    ==
 * ==                                                                        ==
 * ==                                                                        ==
 * ============================================================================

 * ××××××××××××××××××××
 * × VECTORIZATION
 
 - IVDEP
 - LOOP_VECTORIZE              defined for clang and icx, empty for gcc
 - LOOP_VECTOR_LENGTH(N)       defined for clang, empty for others
 - VECTOR_ALWAYS         |
 - VECTOR_ALIGNED        | ->  defined for intel compiler, empty for others
 - VECTOR_UNALIGNED      |

 * ××××××××××××××××××××
 * × LOOPS

 - LOOP_UNROLL                 generic directive for unrollinf
 - LOOP_UNROLL_N(n)            directive for unrolling n times
 
 * ××××××××××××××××××××
 * × ALIGMENT

 - ASSUME_ALIGNED              assume an alignment for an array
 - ATTRIBUTE_ALIGNED           instruct to aligne a static variable
 
 * ============================================================================
 * ==                                                                        ==
 * == end of Dictionary                                                      ==
 * ============================================================================
 */



/* ············································································
 *
 *  INTEL COMPILER
 */

#if defined(__INTEL_LLVM_COMPILER)
#pragma message "using Intel LLVM Compiler"

#define IVDEP            _Pragma("ivdep")
#define LOOP_VECTORIZE   _Pragma("vector")
#define LOOP_VECTOR_LENGTH(N)
#define VECTOR_ALWAYS    _Pragma("vector always")
#define VECTOR_ALIGNED   _Pragma("vector aligned")
#define VECTOR_UNALIGNED _Pragma("vector unaligned")


#define LOOP_UNROLL      _Pragma("unroll")
#define LOOP_UNROLL_N(N) _DO_PRAGMA(unroll N)

#define ASSUME_ALIGNED(V,A)  __builtin_assume_aligned((V), (A))
#define ATTRIBUTE_ALIGNED(A) __attribute__((aligned((A))))

/* ············································································
 *
 *  CLANG
 */

#elif defined(__clang__)
#pragma message "using clang"

#define IVDEP                 _Pragma("clang ivdep")
#define LOOP_VECTORIZE        _DO_PRAGMA(clang loop vectorize( enable ))
#define LOOP_VECTOR_LENGTH(N) _DO_PRAGMA(clang vectorize_width( N ))

#define LOOP_UNROLL           _DO_PRAGMA(clang loop interleave( enable ))
#define LOOP_UNROLL_N(N)      _DO_PRAGMA(clang loop interleave_count( N ))

#define ASSUME_ALIGNED(V,A) __builtin_assume_aligned((V), (A))
#define ATTRIBUTE_ALIGNED(A) __attribute__((__aligned__((A))))

/* ············································································
 *  
 *  GCC
 */

#elif defined(__GNUC__)
#pragma message "using GCC"

#define IVDEP              _Pragma("GCC ivdep")
#define LOOP_VECTORIZE     _Pragma("GVV ivdep")
#define LOOP_VECTOR_LENGTH(N)

#define LOOP_UNROLL        _Pragma("GCC unroll 4")
#define LOOP_UNROLL_N(N)   _DO_PRAGMA(GCC unroll N)

#define ASSUME_ALIGNED(V,A) __builtin_assume_aligned((V), (A))
#define ATTRIBUTE_ALIGNED(A) __attribute__((aligned((A))))

/* ············································································
 *
 *  ARM COMPILER
 */
#elif defined(__CC_ARM)

#error "ARM compilers are not supported yet")
// TBD, sorry

#else

#error "UKNOWN COMPILER USED"

#endif


#if !defined(__INTEL_LLVM_COMPILER)
#define VECTOR_ALWAYS      //_Pragma("message \"vector always not defined\"")
#define VECTOR_ALIGNED     //_Pragma("message \"vector aligned not defined\"")
#define VECTOR_UNALIGNED   //_Pragma("message \"vector unaligned not defined\"")
#endif