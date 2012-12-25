
#ifndef _CI_H_
#define _CI_H_
#include "vel_grid.h"
// типы симметрии
#define NO_SYMM 0	// нет симметрии
#define Z_SYMM  1	// симметрия по оси z
#define YZ_SYMM 2	// симметрия по осям y, z 

// типы потенциалов
#define LJ_POTENTIAL 0	// потенциал Леннарда-Джонса
#define HS_POTENTIAL 1	// потенциал твердых сфер

#define USE_DOUBLE


// устанавливает приблизительный объём сетки Коробова
void set_power (int pow);
// возвращает точный объём сетки Коробова
int get_power ();

// _n_k_rad радицс скоростной сетки
void ci_init (int _n_k_rad, real xcut, int t, int s);

// генерирует расчетную сетку
// xyz2i трехмерный массив который создает соответствие между тремя
// индексами x, y, z и индексом в выпрямленной функции распределения
// т.е. xyz2i[x][y][z] = i
void ci_gen (real tt, const Mapper& xyz2i);

// задает параметры потенциала Леннарда-Джонса
// e = epsilon/(k*T0)
// s = sigma/sigma_eff
void ci_lj_init (real e, real s);

// выполняет шаг протяженностью tt
// f - скоростная функция распределения
void ci_iter (Vel_grid& f1, Vel_grid& f2);

// освобождает всю использованную память
void ci_finalize();

#endif
