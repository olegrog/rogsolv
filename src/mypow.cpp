#include <math.h>
#include "mypow.h"

#define EXP_POLY_DEGREE 5
#define LOG_POLY_DEGREE 5

#define POLY0(x, c0) _mm_set1_ps(c0)
#define POLY1(x, c0, c1) _mm_add_ps(_mm_mul_ps(POLY0(x, c1), x), _mm_set1_ps(c0))
#define POLY2(x, c0, c1, c2) _mm_add_ps(_mm_mul_ps(POLY1(x, c1, c2), x), _mm_set1_ps(c0))
#define POLY3(x, c0, c1, c2, c3) _mm_add_ps(_mm_mul_ps(POLY2(x, c1, c2, c3), x), _mm_set1_ps(c0))
#define POLY4(x, c0, c1, c2, c3, c4) _mm_add_ps(_mm_mul_ps(POLY3(x, c1, c2, c3, c4), x), _mm_set1_ps(c0))
#define POLY5(x, c0, c1, c2, c3, c4, c5) _mm_add_ps(_mm_mul_ps(POLY4(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0))

typedef union {
	__m128i i;
	__m128  x;
	__m128d d;
} un;

__m128 __mm_castsi128_ps(__m128i n) {
	un u;
	u.i = n;
	return u.x;
}
__m128i __mm_castps_si128(__m128 n) {
	un u;
	u.x = n;
	return u.i;
}
__m128i __mm_castpd_si128(__m128d n) {
	un u;
	u.d = n;
	return u.i;
}
__m128d __mm_castsi128_pd(__m128i n) {
	un u;
	u.i = n;
	return u.d;
}

#define POLY0d(x, c0) _mm_set1_pd(c0)
#define POLY1d(x, c0, c1) _mm_add_pd(_mm_mul_pd(POLY0d(x, c1), x), _mm_set1_pd(c0))
#define POLY2d(x, c0, c1, c2) _mm_add_pd(_mm_mul_pd(POLY1d(x, c1, c2), x), _mm_set1_pd(c0))
#define POLY3d(x, c0, c1, c2, c3) _mm_add_pd(_mm_mul_pd(POLY2d(x, c1, c2, c3), x), _mm_set1_pd(c0))
#define POLY4d(x, c0, c1, c2, c3, c4) _mm_add_pd(_mm_mul_pd(POLY3d(x, c1, c2, c3, c4), x), _mm_set1_pd(c0))
#define POLY5d(x, c0, c1, c2, c3, c4, c5) _mm_add_pd(_mm_mul_pd(POLY4d(x, c1, c2, c3, c4, c5), x), _mm_set1_pd(c0))


__m128d exp2d2(__m128d x)
{
   __m128i ipart;
   __m128d fpart, expipart, expfpart;

   x = _mm_min_pd(x, _mm_set1_pd( 1025.0));
   x = _mm_max_pd(x, _mm_set1_pd(-1022.99999999999999));

   // ipart = int(x - 0.5) 
   ipart = _mm_cvtpd_epi32(_mm_sub_pd(x, _mm_set1_pd(0.5)));

   // fpart = x - ipart 
   fpart = _mm_sub_pd(x, _mm_cvtepi32_pd(ipart));

   // expipart = (float) (1 << ipart) 
   expipart = __mm_castsi128_pd(_mm_slli_epi64(_mm_unpacklo_epi32(_mm_add_epi32(ipart, _mm_set1_epi32(1023)), _mm_set1_epi32(0)), 52));

   // minimax polynomial fit of 2**x, in range [-0.5, 0.5[ 
#if EXP_POLY_DEGREE == 5
   expfpart = POLY5d(fpart, 9.9999994e-1, 6.9315308e-1, 2.4015361e-1, 5.5826318e-2, 8.9893397e-3, 1.8775767e-3);
#elif EXP_POLY_DEGREE == 4
   expfpart = POLY4d(fpart, 1.0000026, 6.9300383e-1, 2.4144275e-1, 5.2011464e-2, 1.3534167e-2);
#elif EXP_POLY_DEGREE == 3
   expfpart = POLY3d(fpart, 9.9992520e-1, 6.9583356e-1, 2.2606716e-1, 7.8024521e-2);
#elif EXP_POLY_DEGREE == 2
   expfpart = POLY2d(fpart, 1.0017247, 6.5763628e-1, 3.3718944e-1);
#else
#error
#endif

   return _mm_mul_pd(expipart, expfpart);
}


__m128d log2d2(__m128d x)
{
   __m128i exp = _mm_set1_epi64(_mm_set_pi32(0x7FF00000, 0x00000000));
   __m128i mant = _mm_set1_epi64(_mm_set_pi32(0x000FFFFF, 0x0FFFFFFFF));

   __m128d one = _mm_set1_pd(1.0);

	__m128i i = __mm_castpd_si128(x);

   	__m128 m1  = __mm_castsi128_ps(_mm_srli_epi64(_mm_and_si128(i, exp), 52));

   __m128d e = _mm_cvtepi32_pd(_mm_sub_epi32(__mm_castps_si128(_mm_shuffle_ps(m1, m1, _MM_SHUFFLE(2, 0, 2, 0))), _mm_set1_epi32(1023)));

   __m128d m = _mm_or_pd(_mm_and_pd(x, __mm_castsi128_pd(mant)), one);
   __m128d p;

   // Minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ 
#if LOG_POLY_DEGREE == 6
   p = POLY5d( m, 3.1157899, -3.3241990, 2.5988452, -1.2315303,  3.1821337e-1, -3.4436006e-2);
#elif LOG_POLY_DEGREE == 5
   p = POLY4d(m, 2.8882704548164776201, -2.52074962577807006663, 1.48116647521213171641, -0.465725644288844778798, 0.0596515482674574969533);
#elif LOG_POLY_DEGREE == 4
   p = POLY3d(m, 2.61761038894603480148, -1.75647175389045657003, 0.688243882994381274313, -0.107254423828329604454);
#elif LOG_POLY_DEGREE == 3
   p = POLY2d(m, 2.28330284476918490682, -1.04913055217340124191, 0.204446009836232697516);
#else
#error
#endif


   // This effectively increases the polynomial degree by one, but ensures that log2(1) == 0
   p = _mm_mul_pd(p, _mm_sub_pd(m, one));

   return _mm_add_pd(p, e);
}


