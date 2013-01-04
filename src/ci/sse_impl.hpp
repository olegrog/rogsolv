#ifndef _SSE_IMPL_H_
#define _SSE_IMPL_H_

#include "sse.hpp"

/* inspired by jrfonseca.blogspot.com/2008/09/fast-sse2-pow-tables-or-polynomials.html */ 

namespace sse {

	/* for SSE instructions */
	inline __m128 poly(__m128 , float c0) {
		return _mm_set1_ps(c0);
	}
	inline __m128 poly(__m128 x, float c0, float c1) {
		return _mm_add_ps(_mm_mul_ps(poly(x, c1), x), _mm_set1_ps(c0));
	}
	inline __m128 poly(__m128 x, float c0, float c1, float c2) {
		return _mm_add_ps(_mm_mul_ps(poly(x, c1, c2), x), _mm_set1_ps(c0));
	}
	inline __m128 poly(__m128 x, float c0, float c1, float c2, float c3) {
		return _mm_add_ps(_mm_mul_ps(poly(x, c1, c2, c3), x), _mm_set1_ps(c0));
	}
	inline __m128 poly(__m128 x, float c0, float c1, float c2, float c3, float c4) {
		return _mm_add_ps(_mm_mul_ps(poly(x, c1, c2, c3, c4), x), _mm_set1_ps(c0));
	}
	inline __m128 poly(__m128 x, float c0, float c1, float c2, float c3, float c4, float c5) {
		return _mm_add_ps(_mm_mul_ps(poly(x, c1, c2, c3, c4, c5), x), _mm_set1_ps(c0));
	}
	/* for SSE2 instructions */
	inline __m128d poly(__m128d , double c0) {
		return _mm_set1_pd(c0);
	}
	inline __m128d poly(__m128d x, double c0, double c1) {
		return _mm_add_pd(_mm_mul_pd(poly(x, c1), x), _mm_set1_pd(c0));
	}
	inline __m128d poly(__m128d x, double c0, double c1, double c2) {
		return _mm_add_pd(_mm_mul_pd(poly(x, c1, c2), x), _mm_set1_pd(c0));
	}
	inline __m128d poly(__m128d x, double c0, double c1, double c2, double c3) {
		return _mm_add_pd(_mm_mul_pd(poly(x, c1, c2, c3), x), _mm_set1_pd(c0));
	}
	inline __m128d poly(__m128d x, double c0, double c1, double c2, double c3, double c4) {
		return _mm_add_pd(_mm_mul_pd(poly(x, c1, c2, c3, c4), x), _mm_set1_pd(c0));
	}
	inline __m128d poly(__m128d x, double c0, double c1, double c2, double c3, double c4,
			double c5) {
		return _mm_add_pd(_mm_mul_pd(poly(x, c1, c2, c3, c4, c5), x), _mm_set1_pd(c0));
	}
	/* minimax polynomial fit of 2**x, in range [-0.5, 0.5[ */
	template <size_t i> inline __m128 polyExp(__m128 x); // SSE
	template <> inline __m128 polyExp<2>(__m128 x) {
		return poly(x, 1.0017247f, 6.5763628e-1f, 3.3718944e-1f);
	}
	template <> inline __m128 polyExp<3>(__m128 x) {
		return poly(x, 9.9992520e-1f, 6.9583356e-1f, 2.2606716e-1f, 7.8024521e-2f);
	}
	template <> inline __m128 polyExp<4>(__m128 x) {
		return poly(x, 1.0000026f, 6.9300383e-1f, 2.4144275e-1f, 5.2011464e-2f,
				1.3534167e-2f);
	}
	template <> inline __m128 polyExp<5>(__m128 x) {
		return poly(x, 9.9999994e-1f, 6.9315308e-1f, 2.4015361e-1f, 5.5826318e-2f,
				8.9893397e-3f, 1.8775767e-3f);
	}
	template <size_t i> inline __m128d polyExp(__m128d x); // SSE2
	template <> inline __m128d polyExp<2>(__m128d x) {
		return poly(x, 1.0017247, 6.5763628e-1, 3.3718944e-1);
	}
	template <> inline __m128d polyExp<3>(__m128d x) {
		return poly(x, 9.9992520e-1, 6.9583356e-1, 2.2606716e-1, 7.8024521e-2);
	}
	template <> inline __m128d polyExp<4>(__m128d x) {
		return poly(x, 1.0000026, 6.9300383e-1, 2.4144275e-1, 5.2011464e-2, 1.3534167e-2);
	}
	template <> inline __m128d polyExp<5>(__m128d x) {
		return poly(x, 9.9999994e-1, 6.9315308e-1, 2.4015361e-1,
				5.5826318e-2, 8.9893397e-3, 1.8775767e-3);
	}
	/* minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ */
	template <size_t i> inline __m128 polyLog(__m128 x); // SSE
	template <> inline __m128 polyLog<2>(__m128 x) {
		return poly(x, 2.28330284476918490682f, -1.04913055217340124191f,
				0.204446009836232697516f);
	}
	template <> inline __m128 polyLog<3>(__m128 x) {
		return poly(x, 2.61761038894603480148f, -1.75647175389045657003f,
				0.688243882994381274313f, -0.107254423828329604454f);
	}
	template <> inline __m128 polyLog<4>(__m128 x) {
		return poly(x, 2.8882704548164776201f, -2.52074962577807006663f,
				1.48116647521213171641f, -0.465725644288844778798f,
				0.0596515482674574969533f);
	}
	template <> inline __m128 polyLog<5>(__m128 x) {
		return poly(x, 3.1157899f, -3.3241990f, 2.5988452f, -1.2315303f,
				3.1821337e-1f, -3.4436006e-2f);
	}
	template <size_t i> inline __m128d polyLog(__m128d x); // SSE2
	template <> inline __m128d polyLog<2>(__m128d x) {
		return poly(x, 2.28330284476918490682, -1.04913055217340124191,
				0.204446009836232697516);
	}
	template <> inline __m128d polyLog<3>(__m128d x) {
		return poly(x, 2.61761038894603480148, -1.75647175389045657003,
				0.688243882994381274313, -0.107254423828329604454);
	}
	template <> inline __m128d polyLog<4>(__m128d x) {
		return poly(x, 2.8882704548164776201, -2.52074962577807006663, 1.48116647521213171641,
				-0.465725644288844778798, 0.0596515482674574969533);
	}
	template <> inline __m128d polyLog<5>(__m128d x) {
		return poly(x, 3.1157899, -3.3241990, 2.5988452, -1.2315303,
				3.1821337e-1, -3.4436006e-2);
	}

	inline __m128 exp(__m128 x)
	{
		__m128i ipart;
		__m128 fpart, expipart, expfpart;

		x = _mm_min_ps(x, _mm_set1_ps( 129.00000f));
		x = _mm_max_ps(x, _mm_set1_ps(-126.99999f));

		/* ipart = int(x - 0.5) */
		ipart = _mm_cvtps_epi32(_mm_sub_ps(x, _mm_set1_ps(0.5f)));

		/* fpart = x - ipart */
		fpart = _mm_sub_ps(x, _mm_cvtepi32_ps(ipart));

		/* expipart = (float) (1 << ipart) */
		expipart = _mm_castsi128_ps(_mm_slli_epi32(
				_mm_add_epi32(ipart, _mm_set1_epi32(127)), 23));

		/* minimax polynomial fit of 2**x, in range [-0.5, 0.5[ */
		expfpart = polyExp<EXP_POLY_DEGREE>(fpart);

		return _mm_mul_ps(expipart, expfpart);
	}

	inline __m128 log(__m128 x)
	{
		__m128i exponent = _mm_set1_epi32(0x7F800000);
		__m128i mant = _mm_set1_epi32(0x007FFFFF);

		__m128 one = _mm_set1_ps( 1.0f);

		__m128i i = _mm_castps_si128(x);

		__m128 e = _mm_cvtepi32_ps(_mm_sub_epi32(
				_mm_srli_epi32(_mm_and_si128(i, exponent), 23), _mm_set1_epi32(127)));

		__m128 m = _mm_or_ps(_mm_castsi128_ps(_mm_and_si128(i, mant)), one);

		/* minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ */
		__m128 p = polyLog<LOG_POLY_DEGREE>(m);

		/* this effectively increases the polynomial degree by one,
		 * but ensures that log2(1) == 0 */
		p = _mm_mul_ps(p, _mm_sub_ps(m, one));

		return _mm_add_ps(p, e);
	}

	inline __m128d exp(__m128d x)
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
		expipart = _mm_castsi128_pd(_mm_slli_epi64(_mm_unpacklo_epi32(
				_mm_add_epi32(ipart, _mm_set1_epi32(1023)), _mm_set1_epi32(0)), 52));

		// minimax polynomial fit of 2**x, in range [-0.5, 0.5[ 
		expfpart = polyExp<EXP_POLY_DEGREE>(fpart);

		return _mm_mul_pd(expipart, expfpart);
	}


	inline __m128d log(__m128d x)
	{
		__m128i exponent = _mm_set1_epi64(_mm_set_pi32(0x7FF00000, 0x00000000));
		__m128i mant = _mm_set1_epi64(_mm_set_pi32(0x000FFFFF, 0x0FFFFFFFF));

		__m128d one = _mm_set1_pd(1.0);

		__m128i i = _mm_castpd_si128(x);

		__m128 m1  = _mm_castsi128_ps(_mm_srli_epi64(_mm_and_si128(i, exponent), 52));

		__m128d e = _mm_cvtepi32_pd(_mm_sub_epi32(_mm_castps_si128(
				_mm_shuffle_ps(m1, m1, _MM_SHUFFLE(2, 0, 2, 0))), _mm_set1_epi32(1023)));

		__m128d m = _mm_or_pd(_mm_and_pd(x, _mm_castsi128_pd(mant)), one);
		__m128d p;

		/* Minimax polynomial fit of log2(x)/(x - 1), for x in range [1, 2[ */
		p = polyLog<LOG_POLY_DEGREE>(m);

		/* this effectively increases the polynomial degree by one,
		 * but ensures that log2(1) == 0 */
		p = _mm_mul_pd(p, _mm_sub_pd(m, one));

		return _mm_add_pd(p, e);
	}

	/* +, -, *, / operations */

	inline __m128 mul(__m128 x, __m128 y) {
		return _mm_mul_ps(x, y);
	}
	inline __m128d mul(__m128d x, __m128d y) {
		return _mm_mul_pd(x, y);
	}
	inline d2_t mul(d2_t x, d2_t y) {
		d2_t z;
		z.md = mul(x.md, y.md);
		return z;
	}
	
	/* pow operaton */

	inline __m128 pow(__m128 x, __m128 y) {
		return exp(mul(log(x), y));
	}

	inline __m128d pow(__m128d x, __m128d y) {
		return exp(mul(log(x), y));
	}
	inline d2_t pow(d2_t x, d2_t y) {
		d2_t z;
		z.md = pow(x.md, y.md);
		return z;
	}

}

#endif
