#include <mpi.h>
#include <stdexcept>

#include "first_scheme.h"
#include "printer.h"
#include "box.h"

struct First_scheme::Interior {
	const Axis axis;
	Interior (Axis ax) : axis (ax) { }
	void operator() (Vel_grid& flux, const Vel_grid& grid)
	{
		flux = grid*Vel_grid::velocity (axis);
	}
};

struct First_scheme::Boundary {
	const Side side;																
	Boundary (Side s) : side (s) { }
	void operator() (Vel_grid& flux, Vel_grid& grid, Wall* wall)
	{
		if (Wall_simple* wall_s = dynamic_cast<Wall_simple*> (wall))
			flux = std::abs (sum (flux, side)) * (*wall_s->maxwells[0]);
		if (Wall_grid* wall_gr = dynamic_cast<Wall_grid*> (wall)) {
			wall_gr->update (grid);
			flux = wall_gr->grid () * Vel_grid::velocity (axis (side));
		}
		if (dynamic_cast<Wall_mirror*> (wall))
			flux = -reflect (flux, axis (side));
	}
};

void First_scheme::scheme (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (reflect (side));
		for_each (box->F[axis (s)]->layers_from (side, 1), box->f->all (), Interior (axis (s)));
		mapper ().set_side (side);
		for_each (box->F[axis (s)]->layer (reflect (side), 0), box->f->layer (side, 0), box->wall[reflect (side)]->all (), Boundary (reflect (side)));
	}
}

struct First_scheme::Read_buffer_after {
	void operator() (Vel_grid& flux, Wall* wall)
	{
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall)) {
			First_buffer* buf = dynamic_cast<First_buffer*> (&wall_b->buffer);
			flux = buf->flux (wall_b->read_index);									// read flux from MPI_buffer
		}
	}
};

void First_scheme::read_buffer_after (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (reflect (side));											// read fluxes from wall_box
		for_each (box->F[axis (s)]->layer (side, 0), box->wall[side]->all (), Read_buffer_after ());
	}
}

struct First_scheme::Write_buffer_after {
	void operator() (Vel_grid& flux, Wall* wall)
	{
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall)) {
			First_buffer* buf = dynamic_cast<First_buffer*> (&wall_b->buffer);
			if (axis (buf->this_side) == axis (buf->that_side))
				buf->flux (wall_b->write_index) = flux;							// write flux to MPI_buffer
			else {
				int sign = (direction (buf->this_side)+direction (buf->that_side))%2 ? 1 : -1;
				mapper ().set_side (reflect (buf->that_side));						// use fluxes from box
				buf->flux (wall_b->write_index) = 
					sign*rotate (flux, buf->this_side, buf->that_side);				// rotate and negate before
			}
		}
	}
};

void First_scheme::write_buffer_after (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (side);												// write fluxes to wall_box
		for_each (box->F[axis (s)]->layer (side, 0), box->wall[side]->all (), Write_buffer_after ());
	}
}

void First_buffer::init_buffer ()
{ 
	fluxes = new std::valarray<Half_grid> (capacity);
}

void First_buffer::swap_after (Box* box)
{ 
	First_buffer* buf1 = dynamic_cast<First_buffer*> (twin_buffer (box));
	std::swap (fluxes, buf1->fluxes);
}	

void First_buffer::MPI_swap_after (int dest)
{
	MPI_swap (dest, *fluxes, mapper ().volume ()/2);
}

struct First_scheme::Init_walls {
	const Side side;
	Maxwell_fluxes* maxwell;
	Init_walls (Side s, Maxwell_fluxes* maxw) : side (s), maxwell (maxw) { }
	void operator() (Wall* wall)
	{
		if (Wall_simple* wall_s = dynamic_cast<Wall_simple*> (wall)) {
			Maxwell_fluxes::iterator p = maxwell[side].find (wall_s->params ());
			if (p == maxwell[side].end ()) {
				wall_s->maxwells.push_back (maxwell[side][wall_s->params ()] = new Half_grid);
				Vel_grid tmp_grid (wall_s->temp, 1, wall_s->speed);
				mapper ().set_side ();
				tmp_grid = Vel_grid::velocity (axis (side)) * tmp_grid;
				real norm = std::abs (sum (tmp_grid, side));
				mapper ().set_side (reflect (side));
				*maxwell[side][wall_s->params ()] = 1./norm *tmp_grid;
			}
			else
				wall_s->maxwells.push_back (p->second);
		}
	}
};

void First_scheme::init_boundary (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (reflect (side));
		for_each (box->wall[side]->all (), Init_walls (side, maxwell_fluxes));
	}
}

First_scheme::~First_scheme ()
{
	for (int s=0; s<6; s++)
		for (Maxwell_fluxes::iterator p = maxwell_fluxes[s].begin (); p != maxwell_fluxes[s].end (); ++p)
			delete p->second;
}

void First_scheme::info (Printer* printer)
{
	printer->var ("Scheme type", "first order"); 
}
