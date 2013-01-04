#ifndef _SSE_H_
#define _SSE_H_

#include "emmintrin.h"

namespace sse {

	const size_t EXP_POLY_DEGREE = 5;
	const size_t LOG_POLY_DEGREE = 5;

	union f4 {
		float f[4];
		unsigned int u[4];
		__m128 m;
		__m128i mi;
	};
	typedef union f4 f4_t;

	union d2 {
		double d[2];
		unsigned long long l[2];
		__m128 m;
		__m128i mi;
		__m128d md;
	};
	typedef union d2 d2_t;

	inline __m128 exp(__m128 x);
	inline __m128 log(__m128 x);
	inline __m128 pow(__m128 x, __m128 y);

	inline __m128d exp(__m128d x);
	inline __m128d log(__m128d x);
	inline __m128d pow(__m128d x, __m128d y);
	inline d2_t pow(d2_t x, d2_t y);

}

#include "sse_impl.hpp"

#endif
