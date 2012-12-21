#include "vel_grid.h"

class Init_cond {
public:
	virtual void set_grid (Vel_grid&, Int_vect coord) const = 0;
};

class Const_maxwell : public Init_cond {
	real temp, dens;
	Real_vect speed;
public:
	explicit Const_maxwell (real t = 1, real d = 1, Real_vect s = 0) : temp (t), dens (d), speed (s) { }
	void set_grid (Vel_grid& grid, Int_vect) const {
		grid = Vel_grid (temp, dens, speed);
	};
};

template<class Temp, class Dens, class Speed>
class Arbitrary_maxwell : public Init_cond {
	Temp temp;
	Dens dens;
	Speed speed;
public:
	Arbitrary_maxwell (const Temp& t, const Dens& d, const Speed& s)
		: temp (t), dens (d), speed (s) { }
	void set_grid (Vel_grid& grid, Int_vect r) const {
		grid = Vel_grid (temp (r), dens (r), speed (r));
	};
};

template<class Temp, class Dens, class Speed>
Arbitrary_maxwell<Temp, Dens, Speed>* arbitrary_maxwell (const Temp& t, const Dens& d, const Speed& s)
{
	return new Arbitrary_maxwell<Temp, Dens, Speed> (t, d, s);
}

namespace {
	struct Zero_Vect {
		Real_vect operator () (Int_vect) const { return 0; }
	};
}

template<class Temp, class Dens>
Arbitrary_maxwell<Temp, Dens, Zero_Vect>* arbitrary_maxwell (const Temp& t, const Dens& d)
{
	return arbitrary_maxwell (t, d, Zero_Vect ());
}

template<class Temp, class Dens, class Speed, class QFlow, class Shear>
class Arbitrary_grad13 : public Init_cond {
	Temp temp;
	Dens dens;
	Speed speed;
	QFlow qflow;
	Shear shear;
public:
	Arbitrary_grad13 (const Temp& t, const Dens& d, const Speed& sp, const QFlow& q, const Shear& sh)
		: temp (t), dens (d), speed (sp), qflow (q), shear (sh) { }
	void set_grid (Vel_grid& grid, Int_vect r) const {
		grid = Vel_grid (temp (r), dens (r), speed (r), qflow (r), shear (r));
	};
};

template<class Temp, class Dens, class Speed, class QFlow, class Shear>
Arbitrary_grad13<Temp, Dens, Speed, QFlow, Shear>* arbitrary_grad13 (
	const Temp& t, const Dens& d, const Speed& sp, const QFlow& q, const Shear& sh)
{
	return new Arbitrary_grad13<Temp, Dens, Speed, QFlow, Shear> (t, d, sp, q, sh);
}


