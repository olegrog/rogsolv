#include "limiters.h"
#include "tvd_scheme.h"
namespace {
	inline real coeff (Side side) { return (direction (side) == BACKWARD ? 1. : -1.)/sqrt (3)/Vel_grid::cut_vel (); }
}

template<class T>
struct TVD_scheme<T>::Interior {
	const Side side;
	Interior (Side s) : side (s) { }
	void operator() (Vel_grid& flux, const Vel_grid& grid1, const Vel_grid& grid2, const Vel_grid& grid3)
	{
		const Axis axis = ::axis (side);
		const real c = coeff (side);
		flux = (grid2 +	Lmtr::Limiter<T> () (grid1, grid2, grid3, c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
	}
};

template<class T>
struct TVD_scheme<T>::Interior_last {
	const Side side;
	Interior_last (Side s) : side (s) { }
	void operator() (Vel_grid& flux, const Vel_grid& grid1, const Vel_grid& grid2, Wall* wall)
	{
		const Axis axis = ::axis (side);
		const real c = coeff (side);
		// virtual N+1 cell: grid3=2*grid2-grid1
		if (dynamic_cast<Wall_real*> (wall))
			flux = (grid2 + Lmtr::Limiter<T> () (grid1, grid2, max (2*grid2-grid1, 0), c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
		// real N+1 cell: wall->grid
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall))
			flux = (grid2 + Lmtr::Limiter<T> () (grid1, grid2,
				wall_b->buffer.grid (wall_b->read_index), c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
	}
};

template<class T>
struct TVD_scheme<T>::Interior_first {
	const Side side;
	Interior_first (Side s) : side (s) { }
	void operator() (Vel_grid& flux, const Vel_grid& grid2, const Vel_grid& grid3, Wall* wall)
	{
		const Axis axis = ::axis (side);
		const real c = coeff (side);
		if (dynamic_cast<Wall_mirror*> (wall))
			flux = (grid2 + Lmtr::Limiter<T> () (max (0, reflect (max(0, 2*grid2-grid3)+grid2, axis)-grid2), grid2, grid3, c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
		if (Wall_simple* wall_s = dynamic_cast<Wall_simple*> (wall)) {
			real sum = std::abs (::sum (0.5*(max(0, 2*grid2-grid3)+grid2)*Vel_grid::velocity (axis), side));
			flux = (grid2 + Lmtr::Limiter<T> () (max (0, 2*(sum*(*wall_s->maxwells[1]))-grid2), grid2, grid3, c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
		}
		if (Wall_grid* wall_gr = dynamic_cast<Wall_grid*> (wall)) {
			wall_gr->update (grid2);
			flux = (grid2 + Lmtr::Limiter<T> () (wall_gr->grid (), grid2, grid3, c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
		}
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall))
			flux = (grid2 + Lmtr::Limiter<T> () (wall_b->buffer.grid (wall_b->read_index), grid2, grid3, c*Vel_grid::velocity (axis)))
				* Vel_grid::velocity (axis);
	}
};

template<class T>
struct TVD_scheme<T>::Boundary {
	const Side side;
	Boundary (Side s) : side (s) { }
	void operator() (Vel_grid& flux, Wall* wall)
	{
		if (Wall_simple* wall_s = dynamic_cast<Wall_simple*> (wall))
			flux = std::abs (sum (flux, side)) * (*wall_s->maxwells[0]);
		if (Wall_grid* wall_gr = dynamic_cast<Wall_grid*> (wall))
			flux = wall_gr->grid () * Vel_grid::velocity (axis (side));
		if (dynamic_cast<Wall_mirror*> (wall))
			flux = -reflect (flux, axis (side));
	}
};

template<class T>
void TVD_scheme<T>::scheme (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (reflect (side));
		for_each (box->F[axis (s)]->layers_from_to (side, 2, 1), box->f->layers_from (reflect (side), 2), 
			box->f->layers_from_to (side, 1, 1), box->f->layers_from (side, 2), Interior (side));
		for_each (box->F[axis (s)]->layer (reflect (side), 0), box->f->layer (reflect (side), 1), box->f->layer (reflect (side), 0), 
			box->wall[reflect (side)]->all (), Interior_last (side));
		for_each (box->F[axis (s)]->layer (side, 1), box->f->layer (side, 0), box->f->layer (side, 1), 
			box->wall[side]->all (), Interior_first (side));
	}
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (reflect (side));
		for_each (box->F[axis (s)]->layer (side, 0), box->wall[side]->all (), Boundary (side));
	}
}

#include "printer.h"
template<class T>
void TVD_scheme<T>::info (Printer* printer)
{ 
	printer->var ("Scheme type", "TVD second order");
	printer->var ("Limiter", Lmtr::Limiter<T> ().name ());
}
