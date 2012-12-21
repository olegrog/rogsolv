#ifndef VEL_GRID_AFTER_H
#define VEL_GRID_AFTER_H

// this is fucking unit that nobody must read
// all functions will be inlined when compiled with -O3

#include "auxiliary.h"
#include "vector3.h"
#include "mapper.h"
#include "vel_grid.h"

// the base class Expression, which store used closure
template<class Closure> class Expr {
	Closure clos;
public:
	Expr (const Closure&& cl) : clos (cl) { }
	const Closure& operator () () const { return clos; }
	real operator[] (int i) const { return clos [i]; }
	Expr (Expr&& e) : clos (static_cast<Closure&&>(e.clos)) { }
private:
	Expr (const Closure&); // only move semantic
};

// set of used closures
template<class Oper, class Arg>
class Unary_closure {
	const Arg& arg;
public:
	Unary_closure (const Arg& a) : arg (a) { }
	real operator [] (int i) const { return Oper () (arg[i]); }
};

template<class Oper>
class Unary_closure<Oper, Half_grid> {
	const Half_grid& arg;
public:
	Unary_closure (const Half_grid& a) : arg (a) { }
	real operator [] (int i) const { return Oper () (arg[i]); }
};

template<class Oper>
class Unary_closure<Oper, Vel_grid> {
	const Vel_grid& arg;
public:
	Unary_closure (const Vel_grid& a) : arg (a) { }
	real operator [] (int i) const { return Oper () (arg[i]); }
};

template<class Oper>
class Unary_closure<Oper, Vector> {
	const Vector& arg;
public:
	Unary_closure (const Vector& a) : arg (a) { }
	real operator [] (int i) const { return Oper () (arg[i]); }
};

template<class Oper, class Arg1, class Arg2>
class Binary_closure {
	const Arg1& arg1;
	const Arg2& arg2;
public:
	Binary_closure (const Arg1& a1, const Arg2& a2) : arg1 (a1), arg2 (a2) { }
	real operator [] (int i) const { return Oper () (arg1[i], arg2[i]); }
};

template<class Arg>
class Rotate {
	const Arg& arg;
	const Side side1, side2;
	const Mapper& m;
	const int diameter;
public:
	Rotate (const Arg& a, Side s1, Side s2) : arg (a), side1 (s1), side2 (s2),
		m (mapper ()), diameter (2*mapper ().radius ()-1) { }
	// convert index to vector, rotate, convert back vector to index
	inline real operator[] (int i) const { return arg [m(rotate (m[i], side1, side2, diameter))]; }
};

template<class Arg>
class Reflect {
	const Arg& arg;
	const Axis axis;
	const Mapper& m;
	const int diameter;
public:
	Reflect (const Arg& a, Axis ax) : arg (a), axis (ax), 
		m (mapper ()), diameter (2*mapper ().radius ()-1) { }
	inline real operator[] (int i) const { return arg [m(reflect (m[i], axis, diameter))]; }
};

// unary object-functions
struct Negate {
	real operator () (real a) { return -a; }
};

// binary object-functions
struct Sum {
	real operator () (real a, real b) { return a+b; }
};

struct Difference {
	real operator () (real a, real b) { return a-b; }
};

struct Product {
	real operator () (real a, real b) { return a*b; }
};

struct Max {
	real operator () (real a, real b) { return std::max (a, b); }
};

// object-functions which return scalar
template<class Object>
class Sum_all {
	real sum;
public:
	Sum_all (const Object& object) : sum (0)
	{ 
		const int volume = mapper ().volume ();
		for (int i=0; i<volume; i++)
			sum += object[i];	
	}
	real operator () () { return sum; }
};

template<class Object>
class Sum_half {
	real sum;
public:
	Sum_half (const Object& object, Side side) : sum (0)
	{
		const int volume = mapper ().volume ()/2;
		const std::valarray<int>& index = mapper ().half_ball (side);
		for (int i=0; i<volume; i++)
			sum += object[index[i]];
	}
	real operator () () { return sum; }
};

// define unary operators

#define DEFINE_EXPR_UNARY(Op, Name) 								\
template<class T>													\
inline Expr<Unary_closure<Name, T> > Op (const Expr<T>& e)			\
{																	\
	typedef Unary_closure<Name, T> Clos;							\
	return Expr<Clos> (Clos (e ()));								\
}
#define DEFINE_OBJECT_UNARY(Op, Name, Object) 						\
inline Expr<Unary_closure<Name, Object> > Op (const Object& a)		\
{																	\
	typedef Unary_closure<Name, Object> Clos;						\
	return Expr<Clos> (Clos (a));									\
}
DEFINE_EXPR_UNARY (operator-, Negate)
#undef DEFINE_EXPR_UNARY
DEFINE_OBJECT_UNARY (operator-, Negate, Vel_grid)
DEFINE_OBJECT_UNARY (operator-, Negate, Half_grid)
DEFINE_OBJECT_UNARY (operator-, Negate, Vector)
#undef DEFINE_OBJECT_UNARY

// define binary operators

#define DEFINE_EXPR_EXPR_BINARY(Op, Name) 								\
template<class T1, class T2>													\
inline Expr<Binary_closure<Name, T1, T2> > Op (const Expr<T1>& e1, const Expr<T2>& e2)		\
{																	\
	typedef Binary_closure<Name, T1, T2> Clos;							\
	return Expr<Clos> (Clos (e1 (), e2 ()));										\
}	
#define DEFINE_EXPR_OBJECT_BINARY(Op, Name, Object) 								\
template<class T>													\
inline Expr<Binary_closure<Name, T, Object> > Op (const Expr<T>& e, const Object& a)		\
{																	\
	typedef Binary_closure<Name, T, Object> Clos;							\
	return Expr<Clos> (Clos (e (), a));										\
}																	\
template<class T>													\
inline Expr<Binary_closure<Name, Object, T> > Op (const Object& a, const Expr<T>& e)		\
{																	\
	typedef Binary_closure<Name, Object, T> Clos;							\
	return Expr<Clos> (Clos (a, e ()));										\
}																		\
inline Expr<Binary_closure<Name, Object, Object> > Op (const Object& a1, const Object& a2)		\
{																	\
	typedef Binary_closure<Name, Object, Object> Clos;							\
	return Expr<Clos> (Clos (a1, a2));										\
}
#define DEFINE_OBJECT_OBJECT_BINARY(Op, Name, Object1, Object2) 								\
inline Expr<Binary_closure<Name, Object1, Object2> > Op (const Object1& a1, const Object2& a2)		\
{																	\
	typedef Binary_closure<Name, Object1, Object2> Clos;							\
	return Expr<Clos> (Clos (a1, a2));										\
}																	\
inline Expr<Binary_closure<Name, Object2, Object1> > Op (const Object2& a2, const Object1& a1)		\
{																	\
	typedef Binary_closure<Name, Object2, Object1> Clos;							\
	return Expr<Clos> (Clos (a2, a1));										\
}											

DEFINE_EXPR_EXPR_BINARY (operator+, Sum)
DEFINE_EXPR_EXPR_BINARY (operator-, Difference)
DEFINE_EXPR_EXPR_BINARY (operator*, Product)
#undef DEFINE_EXPR_EXPR_BINARY
DEFINE_EXPR_OBJECT_BINARY (operator+, Sum, Vel_grid)
DEFINE_EXPR_OBJECT_BINARY (operator-, Difference, Vel_grid)
DEFINE_EXPR_OBJECT_BINARY (operator*, Product, Vel_grid)
DEFINE_EXPR_OBJECT_BINARY (operator+, Sum, Half_grid)
DEFINE_EXPR_OBJECT_BINARY (operator-, Difference, Half_grid)
DEFINE_EXPR_OBJECT_BINARY (operator*, Product, Half_grid)
DEFINE_EXPR_OBJECT_BINARY (operator*, Product, Vector)
DEFINE_EXPR_OBJECT_BINARY (operator+, Sum, Vector)
DEFINE_EXPR_OBJECT_BINARY (operator*, Product, Scalar)
DEFINE_EXPR_OBJECT_BINARY (operator-, Difference, Scalar)
DEFINE_EXPR_OBJECT_BINARY (operator-, Difference, Vector)
DEFINE_EXPR_OBJECT_BINARY (operator+, Sum, Scalar)
DEFINE_EXPR_OBJECT_BINARY (max, Max, Scalar)
#undef DEFINE_EXPR_OBJECT_BINARY
DEFINE_OBJECT_OBJECT_BINARY (operator*, Product, Vel_grid, Vector)
DEFINE_OBJECT_OBJECT_BINARY (operator*, Product, Vel_grid, Scalar)
DEFINE_OBJECT_OBJECT_BINARY (operator*, Product, Half_grid, Vector)
DEFINE_OBJECT_OBJECT_BINARY (operator*, Product, Half_grid, Scalar)
DEFINE_OBJECT_OBJECT_BINARY (operator*, Product, Scalar, Vector)
DEFINE_OBJECT_OBJECT_BINARY (max, Max, Scalar, Vel_grid)
DEFINE_OBJECT_OBJECT_BINARY (max, Max, Scalar, Half_grid)
DEFINE_OBJECT_OBJECT_BINARY (max, Max, Half_grid, Vel_grid)
DEFINE_OBJECT_OBJECT_BINARY (operator-, Difference, Half_grid, Vel_grid)
DEFINE_OBJECT_OBJECT_BINARY (operator+, Sum, Scalar, Vel_grid)
DEFINE_OBJECT_OBJECT_BINARY (operator-, Difference, Scalar, Vel_grid)
DEFINE_OBJECT_OBJECT_BINARY (operator-, Difference, Scalar, Half_grid)
DEFINE_OBJECT_OBJECT_BINARY (operator+, Sum, Scalar, Half_grid)

DEFINE_OBJECT_OBJECT_BINARY (operator-, Difference, Scalar, Vector)
DEFINE_OBJECT_OBJECT_BINARY (operator+, Sum, Scalar, Vector)
#undef DEFINE_OBJECT_OBJECT_BINARY
// define functions which return scalar

template<class T>													
inline Scalar sum (const Expr<T>& e)		
{																	
	return Scalar (Sum_all<T> (e ()) ());										
}																	
inline Scalar sum (const Vel_grid& a)		
{																	
	return Scalar (Sum_all<Vel_grid> (a) ());										
}																	

template<class T>													
inline Scalar sum (const Expr<T>& e, Side side)		
{																	
	return Scalar (Sum_half<T> (e (), side) ());										
}																	
inline Scalar sum (const Vel_grid& a, Side side)		
{																	
	return Scalar (Sum_half<Vel_grid> (a, side) ());										
}																	
inline Scalar sum (const Half_grid& a, Side side)		
{																	
	return Scalar (Sum_half<Half_grid> (a, side) ());										
}																	

template<class T>
inline Expr<Rotate<T> > rotate (const Expr<T>& e, Side s1, Side s2) 
{
	typedef Rotate<T> Clos;
	return Expr<Clos> (Clos (e (), s1, s2));
}
inline Expr<Rotate<Vel_grid> > rotate (const Vel_grid& a, Side s1, Side s2) 
{
	typedef Rotate<Vel_grid> Clos;
	return Expr<Clos> (Clos (a, s1, s2));
}

template<class T>
inline Expr<Reflect<T> > reflect (const Expr<T>& e, Axis ax) 
{
	typedef Reflect<T> Clos;
	return Expr<Clos> (Clos (e (), ax));
}
inline Expr<Reflect<Vel_grid> > reflect (const Vel_grid& a, Axis ax) 
{
	typedef Reflect<Vel_grid> Clos;
	return Expr<Clos> (Clos (a, ax));
}


#endif