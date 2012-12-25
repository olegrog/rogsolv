#include <mpi.h>
#include <algorithm>
#include <stdexcept>

#include "box.h"
#include "ci.h"

real Box::H, Box::tau;
Init_cond* Box::init_cond = 0;

struct Integral {
	void operator() (Vel_grid& grid) { ci_iter (grid, grid); }
};

void Collision_integral::operator () (Box* box)
{
	for_each (box->f->all (), Integral ());
}

struct Parameters {
	void operator() (const Vel_grid& f, Features& feat)
	{
		feat.dens = f.dens ();
		feat.flow = f.flow ();
		Real_vect speed = feat.flow/feat.dens;
		feat.press = f.press (speed);
		feat.qflow = f.qflow (speed);
		feat.shear = f.shear (speed);
		feat.temp = feat.press.sum ()/3/feat.dens;
	}
};

void Macroparameters::operator () (Box* box)
{
	mapper ().set_side ();
	for_each (box->f->all (), box->features->all (), Parameters ());
}

void Box::init_distribution ()
{
	f = new Matrix<Vel_grid> (size ());
	for (int ax=0; ax<3; ax++) {
		if (!space ().count (Axis (ax))) continue;
		Int_vect tmp = 0; tmp[ax] = 1;
		F[ax] = new Matrix<Vel_grid> (size ()+tmp);
	}
	if (!init_cond) throw std::logic_error ("Initial conditions have not been initialized");
	mapper ().set_side ();
	for_each_index (f->all (), [this] (Vel_grid& gr, Int_vect c) { init_cond->set_grid (gr, coord () + c); });	// initial distribution function
	for (MPI_Buffer::iterator pbuf = MPI_buffer ().begin (); pbuf != MPI_buffer ().end (); ++pbuf)
		(*pbuf)->init_buffer ();
}

Box::Box (Int_vect sz) : size_ (sz), MPI_rank_ (-1)
{
	f = 0;
	for (int ax=0; ax<3; ax++) {
		F[ax] = 0;
		Int_vect tmp = size (); tmp[ax] = 1;
		wall[side (ax, 0)] = new Matrix<Wall*> (tmp);
		wall[side (ax, 1)] = new Matrix<Wall*> (tmp);
	}
	features = new Matrix<Features> (size ());
}

struct Delete_wall {
	void operator() (Wall*& wall) { if (wall) delete wall; }
};

Box::~Box()
{
	if (f != 0) delete f;
	for (int i=0; i<3; i++) if (F[i]) delete F[i];
	if (features) delete features;
	for (int s=0; s<6; s++) {
		for_each (wall[s]->all (), Delete_wall ());
		delete wall[s];
		if (!space ().count (axis (s))) continue;
	}
	for (MPI_Buffer::iterator pbuf = MPI_buffer ().begin (); pbuf != MPI_buffer ().end (); ++pbuf)
		delete *pbuf;
}

void Box::set_mirror (Side side)
{
	for_each (wall[side]->all (), [=] (Wall*& wall) {
		wall = new Wall_mirror;
	});
}

void Box::set_maxwell (Side side, real temp, real dens)
{
	for_each (wall[side]->all (), [=] (Wall*& wall) {
		wall = new Wall_maxwell (temp, dens);
	});
}

void Box::set_interface (Side side, real temp, real dens, Real_vect speed)
{
	for_each (wall[side]->all (), [=] (Wall*& wall) {
		wall = new Wall_interface (temp, dens, speed);
	});
}

void Box::set_simple (Side side, const Simple_bound& bound)
{
	for_each_index (wall[side]->all (), [&] (Wall*& wall, Int_vect coord){
		wall = new Wall_simple (bound.temp (coord, size ()), bound.speed (coord, size ()));
	});
}
