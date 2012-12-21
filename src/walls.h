#ifndef WALLS_H
#define WALLS_H

#include <vector>
#include <tuple>

#include "buffer.h"

struct Wall {																// abstact class that is responsible for the box walls
	virtual ~Wall () { } 
};

struct Wall_real : public Wall {											// abstact class that is responsible for the real boundary
	virtual ~Wall_real () { } 
};

class Wall_grid : public Wall_real {										// abstact class that is responsible for the boundary returning the Vel_grid
protected:
	Vel_grid* grid_;
public:
	const Vel_grid& grid () const { return *grid_; }
	virtual void update (const Vel_grid&) = 0;
	virtual ~Wall_grid () { if (grid_) delete grid_; }
};

struct Wall_mirror : public Wall_real {										// for symmetry planes
	~Wall_mirror () { }
};

struct Wall_simple : public Wall_real {										// simple reflect from wall
	const real temp;														// wall temperature
	const Real_vect speed;													// wall speed
	std::vector<Half_grid*> maxwells;										// pointer to maxwell grid with the unit of density 
	Wall_simple (real t, Real_vect s) : temp (t), speed (s) { }
	~Wall_simple () { }
	typedef std::tuple <real, Real_vect> Params;
	Params params () const { return Params (temp, speed); }
};

struct Wall_box : public Wall {												// for linking boxes
	Buffer& buffer;															// usable buffer for read/write
	const int read_index, write_index;										// proper package indexes
	std::vector<Vel_grid*> grids;											// additional cells
	Wall_box (Buffer& b, int r, int w) : buffer (b), read_index (r), write_index (w) { }
	~Wall_box () { } 
};

class Wall_maxwell : public Wall_grid {										// constant maxwellian distribution at wall
	real temp, dens;														// gas parameters behind the wall
	bool is_first;
public:
	Wall_maxwell (real t, real d) : temp (t), dens (d), is_first (true) { }
	void update (const Vel_grid& gr)
	{
		if (is_first) {
			grid_ = new Vel_grid (temp, dens, gr.flow () / gr.dens());
			is_first = false;
			return;
		}
		Vel_grid tmp (temp, dens, gr.flow () / gr.dens());
		swap (*grid_, tmp);
	}
};

class Wall_interface : public Wall_grid {									// constant maxwellian distribution at the interface
	real temp, dens;														// gas parameters behind the interface
	Real_vect speed;
	bool is_first;
public:
	Wall_interface (real t, real d, Real_vect s) : temp (t), dens (d), speed (s), is_first (true) { }
	void update (const Vel_grid&)
	{
		if (is_first) {
			grid_ = new Vel_grid (temp, dens, speed);
			is_first = false;
		}
	}
};

struct Simple_bound {														// abstract class for wall temperature assignment
	virtual real temp (Int_vect vec, Int_vect size) const = 0;
	virtual Real_vect speed (Int_vect vec, Int_vect size) const = 0;
};

class Const_bound : public Simple_bound {									// constant temperature
	real temp_;
	Real_vect speed_;
public:
	Const_bound (real t, Real_vect s = 0) : temp_ (t), speed_ (s) { }
	real temp (Int_vect, Int_vect) const { return temp_; }
	Real_vect speed (Int_vect, Int_vect) const { return speed_; }
};

class Linear_temp_bound : public Simple_bound  {								// linear dependence of temperature
	real t0;
	Real_vect temp_;
public:
	Linear_temp_bound (real tt0, Real_vect t) : t0 (tt0), temp_ (t) { }
	real temp (Int_vect vec, Int_vect size) const
	{
		return t0 + (temp_.x-t0)*(vec.x+.5)/size.x + (temp_.y-t0)*(vec.y+.5)/size.y + (temp_.z-t0)*(vec.z+.5)/size.z;
	}
	Real_vect speed (Int_vect, Int_vect) const { return 0; }
};

#endif