#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <boost/math/constants/constants.hpp>				// for pi constant

#include "sse.hpp"
#include "v.hpp"

namespace ci {

	using namespace dod_vector;

	struct node_calc {
		int i1,  i2;
		int i1m, i1l, i2m, i2l;
		double r, c;
	};

	extern int N_nu;
	extern std::vector<node_calc> nc;

	extern int ss[9];
	extern const double exact_hit;

	inline int toInt(double x) { return static_cast<int>(x); }

	inline V3d i2xi(V3i i, int xi_rad) { return V3d(i) + 0.5 - xi_rad; }
	inline V3i xi2i(V3d x, int xi_rad) { return V3i(x  - 0.5 + xi_rad); }

	inline int out_of_sphere_i(V3i i, int xi_rad) {
		return (sqr(i2xi(i, xi_rad)) > sqr(xi_rad));
	}

	inline int out_of_sphere_r(V3d xi, int xi_rad) {
		return (sqr(xi) > sqr(xi_rad));
	}

	const V3d scatter(const V3d& x, double theta, double e);
} // namespace

#include "select_xilm.hpp"

namespace ci {

	using namespace dod_vector;
	
	//процедура вычисляет по начальным скоростям, прицельному растоянию и углу
	//конечные параметры, которые нужны для вычисления интеграла столкновений
	template <typename Map>
	void calc_int_node(V3i xi1, V3i xi2, double b, double e, int nk_rad1, int nk_rad2,
			Map& xyz2i1, Map& xyz2i2, double m1, double m2, double,
			const Particle& p1, const Particle& p2)
	{

		V3d rxi1 = i2xi(xi1, nk_rad1);
		V3d rxi2 = i2xi(xi2, nk_rad2);

//		std::cout << "r " << rxi1 << ' ' << rxi2 << std::endl;

		// первая проверка не выходит ли скорость за пределы сферы
		if ( out_of_sphere_r(rxi1, nk_rad1) || out_of_sphere_r(rxi2, nk_rad2) )	{
			ss[0]++; return;
		}

		N_nu++;

		V3d u = (rxi1+rxi2)/(m1+m2);
		V3d g = rxi2 - m2*u;
			
		/** (I) счет разлетных скоростей **/

		double root = m1 * m2 / (m1 + m2) * std::sqrt(sqr(rxi1/m1 - rxi2/m2));
		double teta = potential->theta(p1, p2, b, root); 

		V3d g1 = scatter(g, teta, e);

		V3d wxi1, wxi2;	// абсолютные скорости

		wxi1 = m1*u - g1; 
		wxi2 = m2*u + g1; 

//		std::cout << "w " << rxi1 << ' ' << rxi2 << std::endl;

		// основное выкидывание из-за того, что разлетные скорости больше скорости обрезания
		if ( out_of_sphere_r(wxi1, nk_rad1) || out_of_sphere_r(wxi2, nk_rad2) ) {
			ss[2]++; return;
		}
		/** (II) подгонка разлетных скоростей к узлам сетки **/
				
		// (x, y, z) скорость, от нее мы будем перебирать все скорости в кубе 
		V3i xi = xi2i(wxi2, nk_rad2); 
		
		// стандартный шаблон (8 точек)
		std::vector<V3i> stencil8;
		for (int s1 = 0; s1 < 2; ++s1)
			for (int s2 = 0; s2 < 2; ++s2)
				for (int s3 = 0; s3 < 2; ++s3)
					stencil8.push_back(xi + V3i(s1, s2, s3));

		/** выбор оптимальной пары точек из шаблона, используя соответствующий метод **/
		V3d w = wxi2 - m2*u;
		V3i xi2l, xi2m;
		double r;
		Select_xilm* method = new Min_delta_p (nk_rad2, m2, sqr(g), u, w);
		if (!(*method) (stencil8, r, xi2l, xi2m)) {
			/** если пара xi2l, xi2m не найдена, то используем расширенный шаблон **/
			// расширенный шаблон (32 точки)
			std::vector<V3i> stencil32;
			for (int s1 = -1; s1 < 3; ++s1)
				for (int s2 = -1; s2 < 3; ++s2)
					for (int s3 = -1; s3 < 3; ++s3)
						if (sqr(s1-0.5)+sqr(s2-0.5)+sqr(s3-0.5) < 3)
							stencil32.push_back(xi + V3i(s1, s2, s3));
			if (!(*method) (stencil32, r, xi2l, xi2m)) {
				ss[6]++; N_nu--; return;
			}
		}
		delete method;

		V3i xi1l = xi1 + xi2 - xi2l;
		V3i xi1m = xi1 + xi2 - xi2m;

		if ( ( (xi1 == xi1l) && (xi2 == xi2l) ) || ( (xi1 == xi1m) && (xi2 == xi2m) ) ) {
			ss[7]++; return;
		}

		if (out_of_sphere_i(xi1l, nk_rad1) || out_of_sphere_i(xi1m, nk_rad1)  ||
				out_of_sphere_i(xi2l, nk_rad2) || out_of_sphere_i(xi2m, nk_rad2)) {
			ss[6]++; return;
		}

		if (std::abs(r-1) < exact_hit) 
			ss[8]++;

		node_calc node; 
		
		node.r = r;

#ifdef DODSOLVER
		node.i1  = xyz2i1[ xi1[0]][ xi1[1]][ xi1[2]];
		node.i1l = xyz2i1[xi1l[0]][xi1l[1]][xi1l[2]];
		node.i1m = xyz2i1[xi1m[0]][xi1m[1]][xi1m[2]];

		node.i2  = xyz2i2[ xi2[0]][ xi2[1]][ xi2[2]];
		node.i2l = xyz2i2[xi2l[0]][xi2l[1]][xi2l[2]];
		node.i2m = xyz2i2[xi2m[0]][xi2m[1]][xi2m[2]];
#else
		node.i1  = xyz2i1(xi1[0], xi1[1], xi1[2]);
		node.i1l = xyz2i1(xi1l[0], xi1l[1], xi1l[2]);
		node.i1m = xyz2i1(xi1m[0], xi1m[1], xi1m[2]);

		node.i2  = xyz2i2(xi2[0], xi2[1], xi2[2]);
		node.i2l = xyz2i2(xi2l[0], xi2l[1], xi2l[2]);
		node.i2m = xyz2i2(xi2m[0], xi2m[1], xi2m[2]);
#endif

		node.c = std::sqrt(sqr(rxi1/m1 - rxi2/m2)) * b;

		nc.push_back(node);

	}

	template <class Grid, class Func>
	void integrate(int rad1, int rad2, double bmax, Grid igrid, Func ifunc)
	{
		for (auto point : igrid) {
			V3i xi1(toInt(point[2] * rad1 * 2),
					toInt(point[3] * rad1 * 2),
					toInt(point[4] * rad1 * 2));

			V3i xi2(toInt(point[5] * rad2 * 2),
					toInt(point[6] * rad2 * 2),
					toInt(point[7] * rad2 * 2));

			double b = point[8] * bmax;
			double e = point[9] * 2*boost::math::constants::pi<double>();

			ifunc(xi1, xi2, b, e);
		}
	}

	template <typename Map, typename IGrid>
	int gen(
		double tt,
		int nk_rad1, int nk_rad2,
		const Map& xyz2i1, const Map& xyz2i2,
		double a, double m1, double m2,
		const Particle& p1, const Particle& p2,
		const IGrid& igrid
	)
	{
		size_t nk1 = 0;
		for (int i1 = 0; i1 < 2*nk_rad1; ++i1)
			for (int i2 = 0; i2 < 2*nk_rad1; ++i2)
				for (int i3 = 0; i3 < 2*nk_rad1; ++i3)
#ifdef DODSOLVER
					if (xyz2i1[i1][i2][i3] >= 0)
#else
					if (xyz2i1(i1, i2, i3) >= 0)
#endif
						++nk1;
		size_t nk2 = 0;
		for (int i1 = 0; i1 < 2*nk_rad2; ++i1)
			for (int i2 = 0; i2 < 2*nk_rad2; ++i2)
				for (int i3 = 0; i3 < 2*nk_rad2; ++i3)
#ifdef DODSOLVER
					if (xyz2i2[i1][i2][i3] >= 0)
#else
					if (xyz2i2(i1, i2, i3) >= 0)
#endif
						++nk2;

		N_nu = 0;
		memset(ss, 0, sizeof(int) * 9);

		nc.clear();

		/** integrate at igrid **/
		integrate (nk_rad1, nk_rad2, potential->bMax(p1, p2), igrid,
			[&] (V3i xi1, V3i xi2, double b, double e) {
				calc_int_node(xi1, xi2, b, e, nk_rad1, nk_rad2, xyz2i1, xyz2i2, m1, m2, a, p1, p2);
			});

// 		std::cout << "n_calc = " << nc.size() << " N_nu = " << N_nu << std::endl;
// 		for (int j = 0; j < 9; j++) 
// 			std::cout << ss[j] << ' ';
// 		std::cout << std::endl;

		double B = sqrt(2) * potential->bMax(p1, p2) * static_cast<double>(nk1 * nk2)
			* std::pow(a, 3) * a / static_cast<double>(N_nu) / 4 * tt;
/*
		std::cout << "m = " << m1 << ' ' << m2 << " d = " << p1.d << ' ' << p2.d << 
					" nk = " << nk1 << ' ' << nk2 << 
					" rad = " << nk_rad1 << ' ' << nk_rad2 << 
					" B = " << B << std::endl;
*/
		double r;
		if (symm == YZ_SYMM) r = 4.0;
		else if (symm == Z_SYMM) r = 2.0;
		else r = 1.0;

		for (std::vector<node_calc>::iterator p = nc.begin(); p != nc.end(); ++p) 
			p->c *= B / r;

		std::random_shuffle(nc.begin(), nc.end());

		return 0;
	}

	template <typename F>
	void iter(F& f1, F& f2) {
		int kneg = 0;
		for (std::vector<node_calc>::iterator p = nc.begin(); p != nc.end(); ++p) {
			if (std::abs(p->r-1) > exact_hit) {
				sse::d2_t x, y, z, w, v;

				x.d[0] = f1[p->i1l];
				x.d[1] = f1[p->i1m];

				z.d[0] = f2[p->i2l];
				z.d[1] = f2[p->i2m];

				w = sse::mul(x, z);
			
				y.d[0] = 1. - p->r;
				y.d[1] = p->r;
			
				v = sse::pow(w, y);

				double rr5 = f1[p->i1];
				double rr6 = f2[p->i2];
				double d = ( - v.d[0] * v.d[1] + rr5 * rr6) * p->c; 	

				double dl = (1. - p->r) * d;
				double dm = p->r * d;

				f1[p->i1l] += dl;
				f2[p->i2l] += dl;
				f1[p->i1m] += dm;
				f2[p->i2m] += dm;
				f1[p->i1] -= d;
				f2[p->i2] -= d;

				if ((f1[p->i1l] < 0) ||
					(f1[p->i1m] < 0) ||
					(f2[p->i2l] < 0) ||
					(f2[p->i2m] < 0) ||
					(f1[p->i1] < 0) || 
					(f2[p->i2] < 0))  {

					f1[p->i1l] = x.d[0];
					f1[p->i1m] = x.d[1];
					f2[p->i2l] = z.d[0];
					f2[p->i2m] = z.d[1];
					f1[p->i1 ] = rr5;
					f2[p->i2 ] = rr6;
					kneg++;
				}
			}
			else {
				double g1 = f1[p->i1];
				double g2 = f2[p->i2];
				double g3 = f1[p->i1m];
				double g4 = f2[p->i2m];

				double d = (- g3*g4 + g1*g2) * p->c;

				f1[p->i1]  -= d;
				f2[p->i2]  -= d;
				f1[p->i1m] += d;
				f2[p->i2m] += d;

				if ((f1[p->i1m] < 0) || 
					(f2[p->i2m] < 0) || 
					(f1[p->i1] < 0) || 
					(f2[p->i2] < 0))  {

					f1[p->i1] = g1;
					f2[p->i2] = g2;
					f1[p->i1m] = g3;
					f2[p->i2m] = g4;
					kneg++;
				}
			}
		}
		if (kneg /*> (0.001*N_nu)*/) {
			std::cout.precision (2);
			std::cout << "There was negative f: " << kneg << ", %N_nu = " << 100.*kneg/N_nu << std::endl;
		}
	}
} // namespace



