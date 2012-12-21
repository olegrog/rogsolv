#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cstring>
#include "ci.h"
#include "mypow.h"

#define POW_FAST

// сами коэффициенты в сетке Коробова
#define NETS 16
static int ka[NETS][9] = {
	{ 10007,	1, 1905, 6491, 6710, 3611, 4146, 2607, 2863 },
	{ 20011,	1, 7191, 2057, 3758, 8928, 5960, 14809, 12988 },
	{ 30011,	1, 13167, 26353, 2769, 26069, 14716, 14556, 8606},
	{ 50021,	1, 11281, 7537, 39218, 32534, 11977, 5816, 32765 },
	{ 100003,	1, 20285, 68883, 49739, 25348, 68757, 93907, 46351 },
	{ 200003,	1, 47369, 188507, 54145, 156036, 158419, 37051, 42494 },
	{ 300017,	1, 81575, 103565, 136172, 101475, 54078, 262899, 170731 },
	{ 500009,	1, 42535, 193663, 307439, 182488, 487373, 37415, 418387 },
	{ 1000003,	1, 417564, 171019, 163483, 410620, 615303, 611111, 188073 },
	{ 1500007,	1, 413996, 388189, 943278, 1496508, 434758, 733031, 985685 },
	{ 2000003,	1, 832685, 1269182, 1228431, 532894, 174792, 458201, 527381 },
	{ 3000017,	1, 368334, 166765, 2871452, 407635, 979274, 1865572, 2703215 },
	{ 4000037,	1, 72362, 210611, 92212, 583028, 681897, 2974319, 1680656 },
	{ 6000011,	1, 1323844, 1723313, 1392620, 251332, 5750225, 908859, 5328166 },
	{ 8000009,	1, 93973, 6914802, 3957321, 907968, 4380879, 3879127, 4791477 },
	{ 10000019,	1, 1833663, 3609180, 7252140, 5522715, 8914182, 4652083, 6262402 }
};

static int POWER;

void set_power (int pow)
{
	for (POWER=0; POWER<NETS; POWER++)
		if (ka[POWER][0] > pow) break;
}
int get_power () { return ka[POWER][0]; }


#define PARTGEN 0.25
//	для вычисления угла разлета
#define NG 300
#define NB 100
#define MAXG 35.0
#define MAXB 2.5
#define hhy 0.00001										// шаг при вычислении интеграла
// структура используется при вычислении аппроксимирующих узлов
struct _node_app {
	real x, y, z;
	real d;
};

typedef struct _node_app node_app;

// структура узлов сетки Коробова именно нужные
struct _node_calc {
	int i1, i2;
	int im, il;
	int im_, il_;
	real r, c;
};
typedef struct _node_calc node_calc;

static int radius;
static real cut;
static real BMAX;
static int N_nu;
static int n_calc, n_calc_alloced;
static node_calc *nc;

static int inc_fl, dec_fl;
static int symm, pot;

static int ss[9];
static real **teta_m;
static real Te_d_T0, OMEGA22;

static inline real frac(real x) {
	if (x > 0) return (x - (int)x);

	else return (x - (int)x - 1);
} 
static inline real sqr(real x) { return x*x; }
// процедура которая будет заранее расчитывать углы разлета в зависимости от g и b
// teta_m[i][j] - массив углов разлета, i - индекс по скоростям g, j - индекс по b
// здесь же задается потенциал
//FILE *fdt;
static void teta_init()
{
	int i;
	real hy = hhy;
	real sqrtOm1 = 1/sqrt(OMEGA22);

	teta_m = (real **)malloc(NG * sizeof(real *));
//	fdt = fopen("teta.out", "w+");
	for (i = 0; i < NG; i++) {

		int j;
		real I, y, b, g2, f;

		teta_m[i] = (real *)malloc(NB * sizeof(real));
		for (j = 0; j < NB; j++) {
			b = (real)j/NB*MAXB;
			g2 = sqr((real)i/NG*MAXG);
			I = 0;
			y = 0;
			for(;;) {
				y += hy/2;
				f = sqr(y/b*sqrtOm1);
				f = sqr(f)*f;
				f = 1 - sqr(y) - 4*4*Te_d_T0/g2*(sqr(f)-f);
				y += hy/2;
				if (f > 0) {
					real dI = hy/sqrt(f);
					I += dI;
					if ((dI * hy) > (M_PI * hhy)) {
						hy = hy / 2;
			//			printf("hy = %f\n", hy);
					}
				}
				else break;
			}
			teta_m[i][j] = M_PI - 2*I;
//			fprintf(fdt, "%f ",  (*teta_m)[i][j]);
		}
//		fprintf(fdt, "\n");
	}
//	fclose(fdt);
}

static real get_teta(real g, real b)
{
	real ig = NG*(g/MAXG);
	real ib = NB*(b/MAXB);
	int i = (int)(ig+0.5);
	int j = (int)(ib+0.5);
	ig = ig - i;
	ib = ib - j;

	real f;
	if ((i>=0) && (i<NG) && (j>=0) && (j<NB)) f = teta_m[i][j];				// линейная интерполяция
	else return 0;

	real d1;
	real d2;
	if (ig > 0)
		if ((i+1) < NG) d1 = teta_m[i+1][j]-f;
		else return 0;
	else 
		if ((i-1) >= 0) d1 = f-teta_m[i-1][j];
		else return 0;
	if (ib > 0) 
		if ((j+1) < NB) d2 = teta_m[i][j+1]-f;
		else { return 0; } 
	else 
		if ((j-1) >= 0) d2 = f-teta_m[i][j-1];
		else return 0;

	return (f+d1*ig+d2*ib);
}
// процедура инициализации
// n - число скоростных узлов по одной координате
// v - массив преобразующий индекс скорости в абсолютную скорость
// xcut - радиус обрезания
// s - тип симметрии; NOSYMM нет симметрии, ZSYMM симметрия по оси Z, YZSYMM - симметрия по осям Y и Z
void ci_init (int rad, real cut_, int t, int s) {

	radius = rad;
	cut = cut_;
	
	symm = s;
	
	pot = t;
	BMAX = 1;

	n_calc_alloced = 0;
	dec_fl = 0;
	inc_fl = 0;
}

void ci_finalize() {
	free(nc);
}

void ci_lj_init (real e, real s)
{
	Te_d_T0 = e;
	OMEGA22 = s;
	BMAX = MAXB;
	teta_init();
}

static inline real i2k(int i, int nk_rad) { return 0.5 + i - nk_rad; }
static inline real k2i(real x, int nk_rad) { return (int)(x - 0.5 + nk_rad); }
static inline int out_of_sphere_i(int i, int j, int k, int nk_rad) {
	real r2 = sqr(i2k(i, nk_rad)) + sqr(i2k(j, nk_rad)) + sqr(i2k(k, nk_rad));
	return (r2 > sqr(nk_rad));
}
static inline int out_of_sphere_r(real x, real y, real z, int nk_rad) {
	real r2 = sqr(x) + sqr(y) + sqr(z);
	return (r2 > sqr(nk_rad));
}

void mixer(int n, node_calc* nc) {
	int i, j;
	node_calc temp;
	for (i = n-1; i > 0; i--) {
   		temp = nc[i];
		j = (int)(((real)rand()/RAND_MAX)*(i+1)); // целое в диапазоне 0..i включительно
//		printf("%d %d\n", j, n);
		nc[i] = nc[j];
		nc[j] = temp;
	}
}

//процедура вычисляет по начальным скоростям, прицельному растоянию и углу
//конечные параметры, которые нужны для вычисления интеграла столкновений
static void calc_int_node(int ix1, int iy1, int iz1, int ix2, int iy2, int iz2, real _b, real _e, int nk_rad1, int nk_rad2,
							const Mapper& xyz2i1, const Mapper& xyz2i2, real m1, real m2, real e1, real e2) {
	int ixl, iyl, izl, ixm, iym, izm;
	int ixl_, iyl_, izl_, ixm_, iym_, izm_;
	real _r, _g;

	real rx1 = i2k(ix1, nk_rad1);
	real ry1 = i2k(iy1, nk_rad1);
	real rz1 = i2k(iz1, nk_rad1);
	real rx2 = i2k(ix2, nk_rad2);
	real ry2 = i2k(iy2, nk_rad2);
	real rz2 = i2k(iz2, nk_rad2);

	// первая проверка не выходит ли скорость за пределы сферы
	if ( out_of_sphere_r(rx1, ry1, rz1, nk_rad1) || out_of_sphere_r(rx2, ry2, rz2, nk_rad2) )
		{ ss[0]++; return; }
	N_nu++;

	real ux = (rx1+rx2)/(m1+m2);
	real uy = (ry1+ry2)/(m1+m2);
	real uz = (rz1+rz2)/(m1+m2);

	real g[3]; 		// начальная относительная скорость
	g[0] = rx2 - m2*ux;
	g[1] = ry2 - m2*uy;
	g[2] = rz2 - m2*uz;
		
	// (I) счет разлетных скоростей
	real gxy = sqr(g[0]) + sqr(g[1]);		// модуль проекции отн. скорости на плосскость ОХУ
	real gg = sqrt(gxy + sqr(g[2]));		// модуль относительной скорости
	gxy = sqrt(gxy);
	_g = sqrt(sqr(rx1/m1-rx2/m2)+sqr(ry1/m1-ry2/m2)+sqr(rz1/m1-rz2/m2));
	real teta;
	if (pot == HS_POTENTIAL) teta = 2*acos(_b);			// угол разлета
	else if (pot == LJ_POTENTIAL) teta = get_teta(_g*sqrt(2*m1*m2/(m1+m2))/sqrt(sqrt(e1*e2)), _b); 	
//	printf("%f\n",_g*sqrt(2*m1*m2/(m1+m2))); 	
	
	real se = sin(_e);				// чтобы не вычислять синусы и косинусы много раз
	real ce = cos(_e);		
	real st = sin(teta);
	real ct = cos(teta);		
	
	real g1[3];
	if (fabs(gxy) > 1E-12) {		// условный ноль
		// относительная разлетная скорость
		g1[0] = g[0]*ct - g[0]*g[2]/gxy*ce*st + g[1]/gxy*gg*se*st;
		g1[1] = g[1]*ct - g[1]*g[2]/gxy*ce*st - g[0]/gxy*gg*se*st;
		g1[2] = g[2]*ct + gxy*ce*st;
	}
	else {
		g1[0] = gg*se*st;
		g1[1] = gg*ce*st;
		g1[2] = gg*ct;
	}

	real w1[3], w2[3];	// абсолютные скорости
	w1[0] = m1*ux - g1[0]; 
	w1[1] = m1*uy - g1[1]; 
	w1[2] = m1*uz - g1[2]; 		
	w2[0] = m2*ux + g1[0]; 
	w2[1] = m2*uy + g1[1]; 
	w2[2] = m2*uz + g1[2]; 		
	// проверка не явл. ли отклонение слишком маленьким
	if ((sqr(w1[0]-rx1)+sqr(w1[1]-ry1)+sqr(w1[2]-rz1)) < 0.5) { ss[1]++; return; }

	// основное выкидывание из-за того, что разлетные скорости больше скорости обрезания
	if ( out_of_sphere_r(w1[0], w1[1], w1[2], nk_rad1) || out_of_sphere_r(w2[0], w2[1], w2[2], nk_rad2) ) {
		ss[2]++; return;
	}
	// (II) подгонка разлетных скоростей к узлам сетки					
			
	// (x, y, z) скорость, от нее мы будем перебирать все скорости в кубе 
	int x = k2i(w2[0], nk_rad2); 
	int y = k2i(w2[1], nk_rad2); // округление идет в меньшую сторону
	int z = k2i(w2[2], nk_rad2);
	// аппроксимирующие узлы
	node_app n1; n1.d = -1E10;			// "бесконечность"
	node_app n2; n2.d = 1E10;
//	node_app n2; n2.d = 0;

	int k = 0, j = 0; 					// флажки которые нужны для определения были ли выбраны узлы 
				
	gg = sqr(gg);	// теперь это наше E0 с ним и будем сравнивать
	int i1, j1, k1; 
	real d;
	for (i1 = x; i1 < (x+2); i1++) 
		for (j1 = y; j1 < (y+2); j1++) 	
			for (k1 = z; k1 < (z+2); k1++) {
				if (!out_of_sphere_i(i1, j1, k1, nk_rad2)) {
					d = sqr(i2k(i1, nk_rad2)-m2*ux) + sqr(i2k(j1, nk_rad2)-m2*uy) + sqr(i2k(k1, nk_rad2)-m2*uz);	// растояние до середины
		//			n1 ближайщий узел подходящий снизу, n2 - сверху
					if (((d - gg) < 0) && ((d - gg) > n1.d)) {
						n1.d = d - gg;
						n1.x = i1;
						n1.y = j1;
						n1.z = k1;
						k = 1;
					}
					else if (((d - gg) > 0) && ((d - gg) < n2.d)) {
						n2.d = d - gg;
						n2.x = i1;
						n2.y = j1;
						n2.z = k1;
						j = 1;
					}
				}
				else { ss[3]++; return; }
			}
				
	if ((k == 0) || (j == 0)) { ss[4]++; N_nu--; return; } // это особые выброшенные узлы 				

	ixl = n1.x;
	iyl = n1.y;
	izl = n1.z;
	ixm = ix1 + ix2 - n1.x;		// координаты точек апроксимирующие вторую разлетную скорость
	iym = iy1 + iy2 - n1.y;
	izm = iz1 + iz2 - n1.z;
				
	if ( out_of_sphere_i(ixm, iym, izm, nk_rad1) )	{ ss[5]++; return; }

	ixl_ = n2.x;
	iyl_ = n2.y;
	izl_ = n2.z;	

	ixm_ = ix1 + ix2 - n2.x;
	iym_ = iy1 + iy2 - n2.y;
	izm_ = iz1 + iz2 - n2.z;

	if ( out_of_sphere_i(ixm_, iym_, izm_, nk_rad1) )	{ ss[6]++; return; }

	_r = - n1.d/(n2.d-n1.d);			// коэффициент характеризующий 

	if ((_r > 1) || (_r < 0)) { ss[7]++; return; }	// на всякий случай 

	nc[n_calc].r = _r;
	nc[n_calc].il  = xyz2i2(ixl,  iyl,  izl);
	nc[n_calc].il_ = xyz2i2(ixl_, iyl_, izl_);
	nc[n_calc].im  = xyz2i1(ixm , iym , izm );
	nc[n_calc].im_ = xyz2i1(ixm_, iym_, izm_);
	nc[n_calc].i1  = xyz2i1(ix1,  iy1,  iz1);
	nc[n_calc].i2  = xyz2i2(ix2,  iy2,  iz2);
	nc[n_calc].c = _b * _g;
	n_calc++;
}

void ci_gen (real tt, const Mapper& xyz2i) {
//void ci_gen (real tt, int K, int nk_rad1, int nk_rad2, int ***xyz2i1, int ***xyz2i2, int nk1, int nk2, real k_vol,
//					real m1, real m2, real d1, real d2, real e1, real e2) {

	
	int j = ka[POWER][0] * PARTGEN;
	if (n_calc_alloced < j) {
		if (n_calc_alloced != 0)
			free(nc);
		n_calc_alloced = j;
		nc = (node_calc *)malloc(sizeof(node_calc) * n_calc_alloced);		
	}

	int l;
	real x1, y1, z1, x2, y2, z2, b, e;
	real rr[8];
	n_calc = 0;
	N_nu = 0;
	memset(ss, 0, sizeof(int) * 9);
	for (l = 0; l < ka[POWER][0]; l++) {
		// генерируем сетку Коробова для расчета интеграла
		rr[0] = (real)rand()/RAND_MAX;  // произвольный вектор смещения
		rr[1] = (real)rand()/RAND_MAX;
		rr[2] = (real)rand()/RAND_MAX;
		rr[3] = (real)rand()/RAND_MAX;
		rr[4] = (real)rand()/RAND_MAX;
		rr[5] = (real)rand()/RAND_MAX;
		rr[6] = (real)rand()/RAND_MAX;
		rr[7] = (real)rand()/RAND_MAX;
		// генерируем скорости
		x1 = (int)(2 * radius * frac(rr[0]+(real)ka[POWER][1]*j/ka[POWER][0]) );
		y1 = (int)(2 * radius * frac(rr[1]+(real)ka[POWER][2]*j/ka[POWER][0]) );
		z1 = (int)(2 * radius * frac(rr[2]+(real)ka[POWER][3]*j/ka[POWER][0]) );
		x2 = (int)(2 * radius * frac(rr[3]+(real)ka[POWER][4]*j/ka[POWER][0]) );
		y2 = (int)(2 * radius * frac(rr[4]+(real)ka[POWER][5]*j/ka[POWER][0]) );
		z2 = (int)(2 * radius * frac(rr[5]+(real)ka[POWER][6]*j/ka[POWER][0]) );
		b = BMAX * frac(rr[6]+(real)ka[POWER][7]*j/ka[POWER][0]);		// прицельное растояние
		e = 2 * M_PI * frac(rr[7]+(real)ka[POWER][8]*j/ka[POWER][0]);	// угол
		calc_int_node(x1, y1, z1, x2, y2, z2, b, e, radius, radius, xyz2i, xyz2i, 1., 1., 1., 1.);		// вычиляем остальные параметры 
	}

/*	printf("n_calc = %d, N_nu = %d ", n_calc, N_nu);
	for (j = 0; j < 9; j++) 
		printf("%d ", ss[j]);
	printf("\n");
*/
	real B = (1/sqrt(2)/M_PI) * 2*M_PI * BMAX * std::pow(cut/radius, 3) * sqr(mapper().volume()) / N_nu / 4 * tt;

	real r = 0;
	if (symm == YZ_SYMM) r = 4.0;
	else if (symm == Z_SYMM) r = 2.0;
	else if (symm == NO_SYMM) r = 1.0;
	assert (r);
	
	real q = cut / radius;
	for (j = 0; j < n_calc; j++) 
		nc[j].c = B * r * nc[j].c * q;
	if (n_calc % 2 != 0) {
		nc[n_calc].r = 0;
		nc[n_calc].i1 = 0;
		nc[n_calc].i2 = 0;
		nc[n_calc].il = 0;
		nc[n_calc].il_ = 0;
		nc[n_calc].im = 0;
		nc[n_calc].im_ = 0;
		nc[n_calc].c = 0;
	}
	mixer(n_calc, nc);
}

#ifdef POW_FAST
#ifdef USE_DOUBLE
void ci_iter(Vel_grid& f1, Vel_grid& f2) {
	
	int i, kneg = 0;
	for (i = 0; i < n_calc; i++) {

			union d2 x, y, z, w, v;

			x.d[0] = f2[nc[i].il ];
			x.d[1] = f2[nc[i].il_];

			z.d[0] = f1[nc[i].im ];
			z.d[1] = f1[nc[i].im_];

			w.md = _mm_mul_pd(x.md, z.md);
		
			y.d[0] = 1-nc[i].r;
			y.d[1] = nc[i].r;
		
			v.md = powd2(w.md, y.md);

			real rr5 = f1[nc[i].i1];
			real rr6 = f2[nc[i].i2];
			real d = ( - v.d[0] * v.d[1] + rr5 * rr6) * nc[i].c; 	
/*
				if ((f2[nc[i].il] < 0) ||
					(f1[nc[i].im] < 0) || 
					(f2[nc[i].il_] < 0) || 
					(f1[nc[i].im_] < 0) || 
					(f1[nc[i].i1] < 0) || 
					(f2[nc[i].i2] < 0)) 
 					std::cout << "!!!f<0: " << f2[nc[i].il] << ' ' << f1[nc[i].im] << ' ' << f2[nc[i].il_] << ' ' <<
 						f1[nc[i].im_] << ' ' << f1[nc[i].i1] << ' ' << f2[nc[i].i2] << '\n';*/

			if (d < 1E10) {				// может случится так pow(0, 0) = nan 
			
				real d1 = (1 - nc[i].r) * d;
				real d2 = nc[i].r * d;

				f2[nc[i].il] += d1;
				f1[nc[i].im] += d1;
				f2[nc[i].il_] += d2;
				f1[nc[i].im_] += d2;
				f1[nc[i].i1] -= d;
				f2[nc[i].i2] -= d;

				if ((f2[nc[i].il] < 0) ||
					(f1[nc[i].im] < 0) || 
					(f2[nc[i].il_] < 0) || 
					(f1[nc[i].im_] < 0) || 
					(f1[nc[i].i1] < 0) || 
					(f2[nc[i].i2] < 0))  {

					f2[nc[i].il ] = x.d[0];
					f1[nc[i].im ] = z.d[0];
					f2[nc[i].il_] = x.d[1];
					f1[nc[i].im_] = z.d[1];
					f1[nc[i].i1 ] = rr5;
					f2[nc[i].i2 ] = rr6;
					kneg++;
				}
			}		
	}
	if (kneg > (0.001*ka[POWER][0])) {
		std::cout.precision (2);
		std::cout << "A lot of Korobov errors: " << kneg << ", % = " << static_cast<real>(kneg)/ka[POWER][0]*100 << std::endl;
		// необходимо увеличить число выбрасываемых сеток Коробова
		// inc_fl = 1;
	}
}
#endif
#endif

