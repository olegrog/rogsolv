#ifndef VEL_GRID_H
#define VEL_GRID_H

#include <valarray>

#include "../base/auxiliary.h"
#include "../base/vector3.h"
#include "mapper.h"

class Scalar {															// scalar type for vel_grid arithmetic
	const real value;
public:
	Scalar (real v) : value (v) { }									// use as implicit constructor
	real operator[] (int) const { return value; }
	operator real () { return value; }
};

class Vector {															// vector type for vel_grid arithmetic
	const std::valarray<real>& v;
	const Axis axis;
	const Mapper& m;
public:
	Vector (Axis ax, const std::valarray<real>& ar) : v (ar), axis (ax), m (mapper ()) { }
	real operator[] (unsigned int i) const { assert (m[i][axis]<(int)v.size ()); return v[m[i][axis]]; }
};

class Half_grid;
template<class> class Expr;												// expression for vel_grid arithmetic
template<class> class All_grid;										// container for copy all vel_grid
class Vel_grid {														// velocity grid
	friend class Half_grid;
	static real cut, dV;												// cutting velocity and volume of velocity cell
	static std::valarray<real> vel;									// 1D array of velocity grid
	const Mapper& m;													// reference to mapper
	std::valarray<real>* v;												// data container
public:
	explicit Vel_grid () : 
		m (mapper ()), v (new std::valarray<real> (m.volume ())) { }	// just allocating memory
	Vel_grid (real temp, 
			  real dens, 
			  Real_vect speed = 0,
			  Real_vect qflow = 0,
			  Real_vect shear = 0
			 );															// creating vel_grid with Grad13 approximation 
	template<class T> const Vel_grid& operator= (const Expr<T>&);		// all vel_grid operation do by this call
	template<class T> const Vel_grid& operator= (const All_grid<T>&);
	const Vel_grid& operator= (const Vel_grid&);
	const Vel_grid& operator= (const Half_grid&);
	~Vel_grid ();
	// access to Vel_grid elements
	real operator[] (unsigned int i) const { assert (i<v->size ()); return (*v)[i]; }
	real& operator[] (unsigned int i) { return (*v)[i]; }				// for CI module
	char* raw_data () const { return reinterpret_cast<char*> (&(*v)[0]); }			// for MPI operations, cache_f
	// all macroparameters
	real dens () const;
	Real_vect flow () const;
	Real_vect press (Real_vect speed) const;
	Real_vect qflow (Real_vect speed) const;
	Real_vect shear (Real_vect speed) const;
	// static functions
	static void set_cut (real);
	static Vector velocity (Axis axis) { return Vector (axis, vel); }
	static real cut_vel () { return cut; }
	static const std::valarray<real>& vel_array () { return vel; }
private:
	Vel_grid (const Vel_grid&);											// prohibited copy constructor
	friend void swap (Vel_grid&, Vel_grid&);
};
inline void swap (Vel_grid& gr1, Vel_grid& gr2) { std::swap (gr1.v, gr2.v); }

class Half_grid {														// half of velocity grid with constant axis and direction
	friend class Vel_grid;
	std::valarray<real>* v;
	const Mapper& m;
public:
	explicit Half_grid () : m (mapper ()) { v = new std::valarray<real> (m.volume ()/2); }
	~Half_grid () { delete v; }
	template<class T> const Half_grid&  operator= (const Expr<T>&);		// all half_grid operation do by this call
	const Half_grid& operator= (const Vel_grid&);
	const Half_grid& operator= (const Half_grid&);
	real operator[] (int i) const { assert (m.rindex (i)>=0 && m.rindex (i)<(int)v->size()); return (*v)[m.rindex (i)]; }
	char* raw_data () { return reinterpret_cast<char*> (&(*v)[0]); }	// for MPI operations, cache_f
private:
	Half_grid (const Half_grid&);										// prohibited constructor
};

template<class T> 
const Vel_grid& Vel_grid::operator= (const Expr<T>& expr)
{
	unsigned int size = m.size ();
	assert (size == v->size () || size == v->size ()/2);
	for (unsigned int i=0; i<size; i++) {
		int j = m.index (i);
		(*v)[j] = expr [j];
	}
	return *this;
}

template<class T> 
const Half_grid& Half_grid::operator= (const Expr<T>& expr)
{
	unsigned int size = m.size ();
	assert (size == v->size ());
	for (unsigned int i=0; i<size; i++)
		(*v)[i] = expr [m.index (i)];
	return *this;
}

template<class F>
void for_each (const Vel_grid& grid, F func)
{
	const Mapper& m = mapper();
	unsigned int size = m.size ();
	for (unsigned int i=0; i<size; i++) {
		int j = m.index (i);
		Int_vect p = m[j];
		func (p, grid[j]);
	}
}

// template<class Arg>
// class All_grid {
// 	const Arg& arg;
// public:
// 	All_grid (const Arg& a) : arg (a) { }
// 	real operator[] (int i) const { return arg[i]; }
// };
// template<class T> 
// const Vel_grid& Vel_grid::operator= (const All_grid<T>& expr)
// {
// 	unsigned int volume = m.volume ();
// 	for (unsigned int i=0; i<volume; i++)
// 		(*v)[i] = expr [i];
// 	return *this;
// }
// 
// inline All_grid<Vel_grid> copy_all (const Vel_grid& gr) { return All_grid<Vel_grid> (gr); }
// template<class T> inline All_grid<T> copy_all (const Expr<T>& e) { return All_grid<T> (e ()); }

// ------------ for debug -----------
void print_vel_grid (const Vel_grid & grid);
void print_half_grid (const Half_grid & grid, Side side);

#include "vel_grid_impl.h"
#endif
