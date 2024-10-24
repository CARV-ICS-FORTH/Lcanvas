
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
   :  We attempt here to uniquely define "our" vector types  :
   :  not using intrinsics but the custom vector types.      :
   :  We exploit the fact that the                           :
   :      __attribute__((vector size (N)))                   :
   :  formalism is widely adopted.                           :
   :                                                         :
   :  This is intended to have only didactical purposes and  :
   :  by no means conveys an exhaustive list and/or          :
   :  combination, nor the maximum tuning for every compiler :
   @ ······················································· @
 */


#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif


// --------------------------------------------------------
//  discover the vector size
// --------------------------------------------------------

#if !defined(VSIZE)

#if !defined(_IMMINTRIN_H_INCLUDED) && !defined(__IMMINTRIN_H)
#include <immintrin.h>
#endif

#define VSIZE_FROM_INTRINSICS
#if defined(__AVX512__)

#warning "found AVX512"
#define VSIZE (sizeof(__m512))
#define vdtype  __m512d
#define vftype  vdtype
#define vitype  __m512
#define vlltype vitype

#elif defined ( __AVX__ ) || defined ( __AVX2__ )

#warning "found AVX/AVX2"

#define VSIZE (sizeof(__m256))
#define vdtype  __m256d
#define vftype  vdtype
#define vitype  __m256
#define vlltype vitype

#elif defined ( __SSE4__ ) || defined ( __SSE3__ )

#warning "found SSE >= 3"
#define VSIZE (sizeof(__m128))
#define vdtype  __m128d
#define vftype  vdtype
#define vitype  __m128
#define vlltype vitype

#else

#define VSIZE (sizeof(double))
#undef VSIZE_FROM_INTRINSICS
#define NO_VECTOR
#endif

#endif

// --------------------------------------------------------
//  VSIZE has been either given or discovered
// --------------------------------------------------------


#if !defined(VSIZE_FROM_INTRINSICS)
#if ( ((VSIZE-1)&VSIZE) > 0 )
#error "the defined vector size is not a power of 2"
#endif

#if (VSIZE <= 8)
#define NO_VECTOR
#endif

#endif

#if defined(NO_VECTOR)
#warning "no vector capability found"
typedef double dvector_t;
typedef float fvector_t;
typedef int ivector_t;

#define DVSIZE 1
#define FVSIZE 1
#define ILLVSIZE 1
#define VALIGN 32

#else


#define DVSIZE  (VSIZE / sizeof(double))
#define FVSIZE  (VSIZE / sizeof(float))
#define IVSIZE  (VSIZE / sizeof(int))
#define LLVSIZE (VSIZE / sizeof(int))
#define VALIGN  (VSIZE)


#if defined(__GNUC__) && !defined(__clang__)

typedef double    dvector_t  __attribute__((vector_size (VSIZE)));
typedef float     fvector_t  __attribute__((vector_size (VSIZE)));
typedef int       ivector_t  __attribute__((vector_size (VSIZE)));
typedef long long llvector_t __attribute__((vector_size (VSIZE)));

#elif defined(__clang__) || defined(__INTEL_LLVM_COMPILER)

typedef double    dvector_t  __attribute__((ext_vector_type (DVSIZE)));
typedef float     fvector_t  __attribute__((ext_vector_type (FVSIZE)));
typedef int       ivector_t  __attribute__((ext_vector_type (IVSIZE)));
typedef long long llvector_t __attribute__((ext_vector_type (LLVSIZE)));
#endif



#if !defined(vdtype)
#define vdtype  dvector_t
#define vftype  fvector_t
#define vitype  ivector_t
#define vlltype llvector_t
#endif

typedef union { dvector_t V;  double    v[DVSIZE]; } dvector_u;
typedef union { fvector_t V;  float     v[FVSIZE]; } fvector_u;
typedef union { ivector_t V;  int       v[IVSIZE]; } ivector_u;
typedef union { llvector_t V; long long v[LLVSIZE];} llvector_u;
#endif
