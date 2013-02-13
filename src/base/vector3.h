#ifndef VECTOR3_H
#define VECTOR3_H

#include <iostream>
#include <iomanip>
#include <boost/math/special_functions/next.hpp>
 
#include "auxiliary.h"

// ----------------------------- Real_vects ---------------------------------
template<class T> struct Vec3 {
	T x, y, z;
	Vec3 (const T& _x, const T& _y, const T& _z) : x (_x), y (_y), z (_z) { }
	Vec3 (const T& t) { x = y = z = t; }
	Vec3 () { x = y = z = T (); }
	Vec3 (const Vec3& v) { __builtin_memcpy (this, &v, 3*sizeof (T)); }
	template<class C> Vec3 (const Vec3<C>& v) : x (static_cast<T> (v.x)), y (static_cast<T> (v.y)), z (static_cast<T> (v.z)) { }
	Vec3& operator= (const Vec3& v) { x = v.x; y = v.y; z = v.z; return *this; }
	Vec3& operator/= (const T& t) { x /= t; y /= t; z /= t; return *this; }
	Vec3& operator*= (const T& t) { x *= t; y *= t; z *= t; return *this; }
	Vec3& operator*= (const Vec3& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
	Vec3& operator+= (const Vec3& v) { x += v.x; y += v.y; z += v.z; return *this; }
	Vec3& operator+= (const T& t) { x += t; y += t; z += t; return *this; }
	Vec3& operator-= (const Vec3& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
	Vec3& operator-= (const T& t) { x -= t; y -= t; z -= t; return *this; }
	T& operator[] (const int i) { return (&x)[i]; }
	const T operator[] (const int i) const { return (&x)[i]; }
	const T sum () const { return x+y+z; }
	const T vol () const { return x*y*z; }
	const T area (int i) const { int a=1; for (int j=0; j<3; j++) if (j!=i) a*=(&x)[j]; return a; }
	Axis max_axis () const
	{
		if (x >= y && x >= z) return XX;
		if (y >= x && y >= z) return YY;
		return ZZ;
	}
};
typedef Vec3<real> Real_vect;
typedef Vec3<int> Int_vect;

template<class T> inline Vec3<T> reflect (const Vec3<T>& vec, Axis ax,  int diameter=0)
{
	Vec3<T> q (vec); q[ax] = diameter-q[ax]; return q;
}

template<class T> inline Vec3<T> abs (const Vec3<T>& vec)
{
	return Vec3<T> (abs (vec.x), abs (vec.y), abs (vec.z));
}

template<class T> inline Vec3<T> rotate (const Vec3<T>& vec, Side s1, Side s2, int diameter=0)
{
	Axis ax1 = axis (s1);
	Axis ax2 = axis (s2);
	if (ax1 == ax2) return vec;
	if ((direction (s1)+direction(s2))%2) std::swap (ax1, ax2);
	Vec3<T> q (vec);
	q[ax1] = diameter - vec[ax2];
	q[ax2] = vec[ax1];
	return q;
}

template<class T> Vec3<T> operator+ (const Vec3<T>& v1, const Vec3<T>& v2)
{
	Vec3<T> v = v1;
	return v += v2;
}

template<class T> Vec3<T> operator+ (const Vec3<T>& v, const T& t)
{
	Vec3<T> vec = v;
	return vec += t;
}

template<class T> Vec3<T> operator+ (const T& t, const Vec3<T>& v)
{
	Vec3<T> vec = v;
	return vec += t;
}

template<class T> Vec3<T> operator- (const Vec3<T>& v1, const Vec3<T>& v2)
{
	Vec3<T> v = v1;
	return v -= v2;
}

template<class T> Vec3<T> operator- (const Vec3<T>& v, const T& t)
{
	Vec3<T> vec = v;
	return vec -= t;
}

template<class T> Vec3<T> operator- (const T& t, const Vec3<T>& v)
{
	Vec3<T> vec = v;
	return vec -= t;
}

template<class T> Vec3<T> operator/ (const Vec3<T>& v, const T& t)
{
	Vec3<T> vec = v;
	return vec /= t;
}

template<class T> Vec3<T> operator* (const Vec3<T>& v, const T& t)
{
	Vec3<T> vec = v;
	return vec *= t;
}

template<class T> Vec3<T> operator* (const T& t, const Vec3<T>& v)
{
	Vec3<T> vec = v;
	return vec *= t;
}

template<class T> Vec3<T> operator* (const Vec3<T>& v1, const Vec3<T>& v2)
{
	Vec3<T> vec = v1;
	return vec *= v2;
}

template<class T> Vec3<T> max (const Vec3<T>& v1, const Vec3<T>& v2)
{
	Vec3<T> v (std::max (v1.x, v2.x), std::max (v1.y, v2.y), std::max (v1.z, v2.z));
	return v;
}

template<class T> Vec3<T> min (const Vec3<T>& v1, const Vec3<T>& v2)
{
	Vec3<T> v (std::min (v1.x, v2.x), std::min (v1.y, v2.y), std::min (v1.z, v2.z));
	return v;
}

template<class T> std::ostream& operator<< (std::ostream& str, const Vec3<T>& v)
{
	return str << '(' << std::setw (3) << v.x << ';' << std::setw (3) << v.y << ';' << std::setw (3) << v.z << ')';
}

template<class T> bool operator== (const Vec3<T>& v1, const Vec3<T>& v2)
{
	for (int i=0; i<3; i++)
		if (v1[i] != v2[i]) return false;
	return true;
}

template<class T> bool operator!= (const Vec3<T>& v1, const Vec3<T>& v2)
{
	for (int i=0; i<3; i++)
		if (v1[i] != v2[i]) return true;
	return false;
}

template<class T> bool operator< (const Vec3<T>& v1, const Vec3<T>& v2)
{
	for (int i=0; i<3; i++)
		if (v1[i] > v2[i]) return false;
	if (std::abs (boost::math::float_distance (v1.z, v2.z)) <= 1) return false; else return true;
}

template<class T> bool operator> (const Vec3<T>& v1, const Vec3<T>& v2)
{
	for (int i=0; i<3; i++)
		if (v1[i] < v2[i]) return false;
	if (std::abs (boost::math::float_distance (v1.z, v2.z)) <= 1) return false; else return true;
}

template<class T> T shear_matrix_product (const Vec3<T>& tau, const Vec3<T>& c)
{
	T result = T ();
	for (int i=0; i<3; i++)
		result += 2 * tau[i] * c[(i+1)%3] * c[(i+2)%3];
	return result;
}

template<class T> T dot (const Vec3<T>& v1, const Vec3<T>& v2)
{
	T result = T ();
	for (int i=0; i<3; i++)
		result += v1[i]*v2[i];
	return result;
}

template<class T> T sqr (const Vec3<T>& v)
{ 
	return dot (v, v);
}

#endif
