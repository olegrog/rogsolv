#ifndef _KOROBOV_H_
#define _KOROBOV_H_

#include <cstdlib>
#include <climits>

namespace korobov {

	const int dimension = 10;
	const int coefficients[][1 + dimension] = {
		{    10007,	1, 1, 1,    1905,    6491,    6710,    3611,    4146,    2607,    2863 },
		{    20011,	1, 1, 1,    7191,    2057,    3758,    8928,    5960,   14809,   12988 },
		{    30011,	1, 1, 1,   13167,   26353,    2769,   26069,   14716,   14556,    8606 },
		{    50021,	1, 1, 1,   11281,    7537,   39218,   32534,   11977,    5816,   32765 }, 
		{   100003,	1, 1, 1,   20285,   68883,   49739,   25348,   68757,   93907,   46351 },
		{   200003,	1, 1, 1,   47369,  188507,   54145,  156036,  158419,   37051,   42494 },
		{   300017,	1, 1, 1,   81575,  103565,  136172,  101475,   54078,  262899,  170731 },
		{   500009,	1, 1, 1,   42535,  193663,  307439,  182488,  487373,   37415,  418387 },
		{  1000003,	1, 1, 1,  417564,  171019,  163483,  410620,  615303,  611111,  188073 },
		{  1500007,	1, 1, 1,  413996,  388189,  943278, 1496508,  434758,  733031,  985685 },
		{  2000003,	1, 1, 1,  832685, 1269182, 1228431,  532894,  174792,  458201,  527381 },
		{  3000017,	1, 1, 1,  368334,  166765, 2871452,  407635,  979274, 1865572, 2703215 },
		{  4000037,	1, 1, 1,   72362,  210611,   92212,  583028,  681897, 2974319, 1680656 },
		{  6000011,	1, 1, 1, 1323844, 1723313, 1392620,  251332, 5750225,  908859, 5328166 }, 
		{  8000009,	1, 1, 1,   93973, 6914802, 3957321,  907968, 4380879, 3879127, 4791477 }, 
		{ 10000019,	1, 1, 1, 1833663, 3609180, 7252140, 5522715, 8914182, 4652083, 6262402 },
		{ 1000000000, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
	};

	inline double frac(double x) {
		if (x > 0) return (x - (int)x);
		else if (x == 0) return 0;
		else return (x - (int)x - 1);
	} 

	class Point {
		public:
			Point(const int* line, const double* shift, int s) :
					line(line), shift(shift), s(s) {}
			double operator[](size_t i) const {
				return frac(shift[i] + 
						static_cast<double>(line[i+1]) / 
						line[0] * (s + 1));
			}

		private:
			friend class Iterator;
			const int* line;
			const double* shift;
			int s;
	};

	class Iterator {
		public:
			Iterator(const int* line, const double* shift, int s = 0) : 
					point(line, shift, s) {}

			void operator++() { ++point.s; }
			bool operator!=(const Iterator& other) {
				return point.s != other.point.s;
			}
		
			Point& operator*() { return point; }

		private:
			Point point;
	};

	class Grid {
		public:
			Grid(int size = 0) { resize(size); update(); }

			void resize(int size); 
			int size() const { return sz; }
			typedef Iterator iterator;
			iterator begin() {
				return iterator(coefficients[line], random_shift);
			}
			iterator end() {
				return iterator(coefficients[line], random_shift, coefficients[line][0]);
			}
			void update() {
				for (int i = 0; i < dimension; ++i) 
					random_shift[i] = static_cast<double>(std::rand())/RAND_MAX;
			}
		private:
			int sz, line;
			double random_shift[dimension];
	};

	inline void Grid::resize(int size) {
//		std::cout << "korobov_size = " << size << ' ';
		for (line = 0; coefficients[line][0] < size; ++line) {}
		sz = coefficients[line][0];
//		std::cout << "sz = " << sz << ' ' << std::endl;
		update();
	}

}

#endif
