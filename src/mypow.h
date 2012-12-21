#ifndef _MYPOW_H_
#define _MYPOW_H_

#include "emmintrin.h"

union f4 {
   int i[4];
   unsigned int u[4];
   float f[4];
   __m128 m;
   __m128i mi;
};

union d2 {
   long l[2];
   int i[4];
   double d[2];
   __m128 m;
   __m128i mi;
   __m128d md;
};

__m128 exp2f4(__m128 x);

__m128 log2f4(__m128 x);

static inline __m128 powf4(__m128 x, __m128 y) {
   return exp2f4(_mm_mul_ps(log2f4(x), y));
}

__m128d exp2d2(__m128d x);

__m128d log2d2(__m128d x);

static inline __m128d powd2(__m128d x, __m128d y) {
   return exp2d2(_mm_mul_pd(log2d2(x), y));
}


#endif
