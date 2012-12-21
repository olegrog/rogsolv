#ifndef MATRIX_H
#define MATRIX_H

#include <valarray>

#include "auxiliary.h"
#include "vector3.h"


#include "vel_grid.h"

// ----------------------------- Slice ---------------------------------
struct Slice {
	int start;
	Int_vect size, stride;
	Slice () : start (), size (), stride () { }
	Slice (const int s, const Int_vect& sz, const Int_vect& st) : start (s), size (sz), stride (st) { }
};

template<class T> class Slice_iter;

// ----------------------------- Matrix ---------------------------------
template<class T> class Matrix {
	std::valarray<T>* v;
	const Int_vect size;
	Matrix (const Matrix&);
	Matrix& operator= (const Matrix&);
public:
	Matrix (Int_vect s) : size (s) { v = new std::valarray<T> (size.vol ()); }
	~Matrix () { delete v; }
	const T& operator[] (int index) const { return (*v)[index]; }
	operator T* () { return &((*v)[0]); }			
	Slice_iter<T> layers (Side, int, int) const;
	Slice_iter<T> layers_from (Side, int) const;
	Slice_iter<T> layers_from_to (Side, int, int) const;
	Slice_iter<T> layer (Side, int) const;
	Slice_iter<T> all () const;
	Slice_iter<T> bound_rect (Side, Int_vect, Int_vect) const;
	Slice_iter<T> volume (Int_vect, Int_vect) const;
	friend class Slice_iter<T>;
};

template<class T> class Slice_iter {
	const Matrix<T>& m;
	const Slice sl;
	Slice_iter (const Matrix<T>& mm, const Slice s) : m (mm), sl (s) { }
public:
	Slice_iter& operator= (const Slice_iter&);
	Slice_iter& operator= (const T&);
	T& operator[] (unsigned int index) { assert (index < m.v->size ()); return (*m.v)[index]; }
	const T& operator[] (unsigned int index) const { assert (index < m.v->size ()); return (*m.v)[index]; }
	const T& operator() (Int_vect p) { return (*m.v)[sl.start+p.x*sl.stride[0]+p.y*sl.stride[1]+p.z*sl.stride[2]]; } // for bkviewer and write_f
	friend class Matrix<T>;
	template<class C1, class C2, class F>
		friend void for_each (Slice_iter<C1>, Slice_iter<C2>, F);
	template<class C1, class C2, class C3, class F>
		friend void for_each (Slice_iter<C1>, Slice_iter<C2>, Slice_iter<C3>,  F);
	template<class C1, class C2, class C3, class C4, class F>
		friend void for_each (Slice_iter<C1>, Slice_iter<C2>, Slice_iter<C3>, Slice_iter<C4>, F);
	template<class C, class F>
		friend void for_each (Slice_iter<C>, F);
	template<class C, class F>
		friend void for_each_index (Slice_iter<C>, F);
	template<class C, class S>
		friend void update (Slice_iter<C>, Matrix<C>*[], real, S);
};

#define IND(slice) (slice.start+i*slice.stride[0]+j*slice.stride[1]+k*slice.stride[2])
#define FOR_ALL							\
	for (k=0; k<sl.size[2]; k++)				\
		for (j=0; j<sl.size[1]; j++)			\
			for (i=0; i<sl.size[0]; i++)

template<class T> Slice_iter<T>& Slice_iter<T>::operator= (const Slice_iter<T>& sl_it)
{
	int i, j, k;
	const Slice sl2 = sl_it.sl;
	FOR_ALL
		(*m.v)[IND(sl)] = sl_it[IND(sl2)];
	return *this;
}

template<class T> Slice_iter<T>& Slice_iter<T>::operator= (const T& value)
{
	int i, j, k;
	FOR_ALL {
		(*m.v)[IND(sl)] = value;
	}
	return *this;
}

template<class T> Slice_iter<T> Matrix<T>::layers (Side side, int begin, int end) const
{
	const Int_vect stride (1, size.x, size.x*size.y);
	const int start = stride[axis (side)]*(direction (side) == BACKWARD ? begin : size[axis (side)]-1-end);
	Int_vect length (size);
	length[axis (side)] = end-begin+1;
	return Slice_iter<T> (*this, Slice (start, length, stride));
}

template<class T> inline Slice_iter<T> Matrix<T>::layers_from (Side side, int begin) const
{
	return layers (side, begin, size[axis (side)]-1);
}

template<class T> inline Slice_iter<T> Matrix<T>::layers_from_to (Side side, int begin, int end) const
{
	return layers (side, begin, size[axis (side)]-1-end);
}

template<class T> inline Slice_iter<T> Matrix<T>::all () const
{
	return layers_from (Side (0), 0);
}

template<class T> inline Slice_iter<T> Matrix<T>::layer (Side side, int num) const
{
	return layers (side, num, num);
}

template<class T> Slice_iter<T> Matrix<T>::bound_rect (Side side, Int_vect begin, Int_vect end) const
{
	const Int_vect stride (1, size.x, size.x*size.y);
	Int_vect length (end-begin);
	int start = (begin*stride).sum ();
	length[axis (side)] = 1;
	start -= begin[axis (side)]*stride[axis (side)];
	if (direction (side) == FORWARD) start += (size[axis (side)]-1)*stride[axis (side)];
	return Slice_iter<T> (*this, Slice (start, length, stride));
}

template<class T> inline Slice_iter<T> Matrix<T>::volume (Int_vect begin, Int_vect dim) const
{
	const Int_vect stride (1, size.x, size.x*size.y);
	const int start = (begin*stride).sum ();
	const Int_vect length (dim);
	return Slice_iter<T> (*this, Slice (start, length, stride));
}

template<class C1, class C2, class F>
void for_each (Slice_iter<C1> sl_it1, Slice_iter<C2> sl_it2, F func)
{
	int i, j, k;
	const Slice sl = sl_it1.sl;
	const Slice sl2 = sl_it2.sl;
	FOR_ALL
		func (sl_it1[IND(sl)], sl_it2[IND(sl2)]);
}

template<class C1, class C2, class C3, class F>
void for_each (Slice_iter<C1> sl_it1, Slice_iter<C2> sl_it2, Slice_iter<C3> sl_it3, F func)
{
	int i, j, k;
	const Slice sl = sl_it1.sl;
	const Slice sl2 = sl_it2.sl;
	const Slice sl3 = sl_it3.sl;
	FOR_ALL
		func (sl_it1[IND(sl)], sl_it2[IND(sl2)], sl_it3[IND(sl3)]);
}

template<class C1, class C2, class C3, class C4, class F>
void for_each (Slice_iter<C1> sl_it1, Slice_iter<C2> sl_it2, Slice_iter<C3> sl_it3, Slice_iter<C4> sl_it4, F func)
{
	int i, j, k;
	const Slice sl = sl_it1.sl;
	const Slice sl2 = sl_it2.sl;
	const Slice sl3 = sl_it3.sl;
	const Slice sl4 = sl_it4.sl;
	FOR_ALL
		func (sl_it1[IND(sl)], sl_it2[IND(sl2)], sl_it3[IND(sl3)], sl_it4[IND(sl4)]);
}

template<class C, class F>
void for_each (Slice_iter<C> sl_it, F func)
{
	int i, j, k;
	const Slice sl = sl_it.sl;
	FOR_ALL
		func (sl_it[IND(sl)]);
}

template<class C, class F>
void for_each_index (Slice_iter<C> sl_it, F func)
{
	int i, j, k;
	const Slice sl = sl_it.sl;
	FOR_ALL
		func (sl_it[IND(sl)], Int_vect (i, j, k));
}

template<class C, class S>
void update (Slice_iter<C> f, Matrix<C>* F[], real coeff, S space)
{
	int i, j, k;
	Vec3<Slice> p, n;
	for (int ax=0; ax<3; ax++) {
		if (!space.count (Axis(ax))) continue;
		p[ax] = F[ax]->layers_from (side (ax, BACKWARD), 1).sl;
		n[ax] = F[ax]->layers_from (side (ax, FORWARD), 1).sl;
	}
	const Slice sl = f.sl;
	if (space.size () == 1) {
		const Axis ax = *space.begin ();
		FOR_ALL
			f[IND(sl)] = f[IND(sl)] - Scalar (coeff)*((*F[ax])[IND(p[ax])] - (*F[ax])[IND(n[ax])]);
	}
	if (space.size () == 2) {
		const Axis ax1 = *space.begin ();
		const Axis ax2 = *space.rbegin ();
		FOR_ALL
			f[IND(sl)] = f[IND(sl)] - Scalar (coeff)*((*F[ax1])[IND(p[ax1])] - (*F[ax1])[IND(n[ax1])] + (*F[ax2])[IND(p[ax2])] - (*F[ax2])[IND(n[ax2])]);
	}
	if (space.size () == 3) {
		FOR_ALL
			f[IND(sl)] = f[IND(sl)] - Scalar (coeff)*((*F[0])[IND(p[0])] - (*F[0])[IND(n[0])] + (*F[1])[IND(p[1])] - (*F[1])[IND(n[1])] + (*F[2])[IND(p[2])] - (*F[2])[IND(n[2])]);
	}
}

#undef IND
#undef FOR_ALL

#endif
