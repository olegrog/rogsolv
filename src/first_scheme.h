#ifndef FIRST_SCHEME_H
#define FIRST_SCHEME_H

#include "difference_scheme.h"
#include "walls.h"

// buffer for first order scheme
class First_buffer : public Buffer {
	std::valarray<Half_grid>* fluxes;
public:
	First_buffer (int c, Box* b, Side is, Side at) : Buffer (c, b, is, at), fluxes (0) { }
	void init_buffer ();
	~First_buffer () { if (fluxes) delete fluxes; }
	Half_grid& flux (int i) { assert (fluxes); assert (i>=0 && i<capacity); return (*fluxes)[i]; }
	void swap_before (Box*) { }
	void swap_after (Box*);
	void MPI_swap_before (int) { }
	void MPI_swap_after (int);
private:
	Vel_grid& grid (int) { assert (false); return *(new Vel_grid); }	// prohibited in First_buffer
};

// conservative first order flux scheme
class First_scheme : public Difference_scheme {
public:
	void read_buffer_before (Box*) { }
	void write_buffer_before (Box*) { }
	void read_buffer_after (Box*);
	void write_buffer_after (Box*);
	void scheme (Box*);
	Buffer* create_buffer (int c, Box* b, Side is, Side at) { return new First_buffer (c, b, is, at); }
	void init_boundary (Box*);
	void info (Printer*);
	~First_scheme ();
private:
	Maxwell_fluxes maxwell_fluxes [6];				// cache of maxwell fluxes for all six walls
	class Interior;
	class Boundary;
	class Read_buffer_after;
	class Write_buffer_after;
	class Init_walls;
};

#endif // FIRST_SCHEME_H
