#include "tvd_scheme.h"
#include "box.h"

struct TVD_scheme_impl::Read_buffer_after {
	void operator() (Vel_grid& flux, Wall* wall)
	{
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall)) {
			TVD_buffer* buf = static_cast<TVD_buffer*> (&wall_b->buffer);
			flux = buf->flux (wall_b->read_index);									// read flux from MPI_buffer
		}
	}
};

void TVD_scheme_impl::read_buffer_after (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (reflect (side));											// read fluxes from wall_box
		for_each (box->F[axis (s)]->layer (side, 0), box->wall[side]->all (), Read_buffer_after ());
	}
}

struct TVD_scheme_impl::Write_buffer_after {
	void operator() (Vel_grid& flux, Wall* wall)
	{
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall)) {
			TVD_buffer* buf = static_cast<TVD_buffer*> (&wall_b->buffer);
			assert (buf);
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

void TVD_scheme_impl::write_buffer_after (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		mapper ().set_side (side);												// write fluxes to wall_box
		for_each (box->F[axis (s)]->layer (side, 0), box->wall[side]->all (), Write_buffer_after ());
	}
}

struct TVD_scheme_impl::Write_buffer_before {
	void operator() (Vel_grid& grid, Wall* wall)
	{
		if (Wall_box* wall_b = dynamic_cast<Wall_box*> (wall)) {
			TVD_buffer* buf = static_cast<TVD_buffer*> (&wall_b->buffer);
			if (axis (buf->this_side) == axis (buf->that_side))
				buf->grid (wall_b->write_index) = grid;							// write grid to MPI_buffer
			else {
				int sign = (direction (buf->this_side)+direction (buf->that_side))%2 ? 1 : -1;
				buf->grid (wall_b->write_index) = 
					sign*rotate (grid, buf->this_side, buf->that_side);				// rotate and negate before
			}
		}
	}
};

void TVD_scheme_impl::write_buffer_before (Box* box)
{
	mapper ().set_side ();
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		for_each (box->f->layer (side, 0), box->wall[side]->all (), Write_buffer_before ());
	}
}

void TVD_buffer::init_buffer ()
{ 
	fluxes = new std::valarray<Half_grid> (capacity);
	grids = new std::valarray<Vel_grid> (capacity);
}

void TVD_buffer::swap_before (Box* box)
{ 
	TVD_buffer* buf2 = dynamic_cast<TVD_buffer*> (twin_buffer (box));
	std::swap (grids, buf2->grids);
}	

void TVD_buffer::swap_after (Box* box)
{ 
	TVD_buffer* buf2 = dynamic_cast<TVD_buffer*> (twin_buffer (box));
	std::swap (fluxes, buf2->fluxes);
}	

void TVD_buffer::MPI_swap_before (int dest)
{
	MPI_swap (dest, *grids, mapper ().volume ());
}

void TVD_buffer::MPI_swap_after (int dest)
{
	MPI_swap (dest, *fluxes, mapper ().volume ()/2);
}

struct TVD_scheme_impl::Init_walls {
	const Side side;
	Maxwell_fluxes* maxwell;
	Maxwell_fluxes* maxwell2;
	Init_walls (Side s, Maxwell_fluxes* maxw, Maxwell_fluxes* maxw2)
		: side (s), maxwell (maxw), maxwell2 (maxw2) { }
	void operator() (Wall* wall)
	{
		if (Wall_simple* wall_s = dynamic_cast<Wall_simple*> (wall)) {
			Maxwell_fluxes::iterator p = maxwell[side].find (wall_s->params ());
			Maxwell_fluxes::iterator p2 = maxwell2[side].find (wall_s->params ());
			// initializate cache of maxwellians
			if (p == maxwell[side].end ()) {
				wall_s->maxwells.push_back (maxwell[side][wall_s->params ()] = new Half_grid);
				wall_s->maxwells.push_back (maxwell2[side][wall_s->params ()] = new Half_grid);
				Vel_grid tmp_grid (wall_s->temp, 1, wall_s->speed);
				mapper ().set_side (reflect (side));
				*maxwell2[side][wall_s->params ()] = tmp_grid;
				mapper ().set_side ();
				tmp_grid = Vel_grid::velocity (axis (side)) * tmp_grid;
				real norm = std::abs (sum (tmp_grid, side));
				mapper ().set_side (reflect (side));
				*maxwell[side][wall_s->params ()] = 1./norm * tmp_grid;
				*maxwell2[side][wall_s->params ()] = 1./norm * (*maxwell2[side][wall_s->params ()]);
			}
			// use this cache
			else {
				wall_s->maxwells.push_back (p->second);
				wall_s->maxwells.push_back (p2->second);
			}
		}
	}
};

void TVD_scheme_impl::init_boundary (Box* box)
{
	for (int s=0; s<6; s++) {
		Side side = static_cast<Side> (s);
		if (!box->space ().count (axis (s))) continue;
		for_each (box->wall[side]->all (), Init_walls (side, cache1, cache2));
	}
}

TVD_scheme_impl::~TVD_scheme_impl ()
{
	for (int s=0; s<6; s++) {
		for (Maxwell_fluxes::iterator p = cache1[s].begin (); p != cache1[s].end (); ++p)
			delete p->second;
		for (Maxwell_fluxes::iterator p = cache2[s].begin (); p != cache2[s].end (); ++p)
			delete p->second;
	}
}

