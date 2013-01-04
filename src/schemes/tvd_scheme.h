#ifndef TVD_SCHEME_H
#define TVD_SCHEME_H

#include "difference_scheme.h"
#include "limiters.h"

// buffer for TVD scheme
class TVD_buffer : public Buffer {
	std::valarray<Half_grid>* fluxes;
	std::valarray<Vel_grid>* grids;
public:
	TVD_buffer (int c, Box* b, Side is, Side at) : Buffer (c, b, is, at), fluxes (0), grids (0)  { }
	void init_buffer ();
	~TVD_buffer () { if (fluxes) delete fluxes; if (grids) delete grids; }
	Half_grid& flux (int i) { assert (fluxes); assert (i>=0 && i<capacity); return (*fluxes)[i]; }
	Vel_grid& grid (int i) { assert (grids); assert (i>=0 && i<capacity); return (*grids)[i]; }
	void swap_before (Box*);
	void swap_after (Box*);
	void MPI_swap_before (int);
	void MPI_swap_after (int);
};

// conservative second order monotonous TVD scheme
class TVD_scheme_impl {
public:
	void write_buffer_before (Box*);
	void read_buffer_after (Box*);
	void write_buffer_after (Box*);
	void init_boundary (Box*);
	~TVD_scheme_impl ();
private:
	typedef Difference_scheme::Maxwell_fluxes Maxwell_fluxes;
	Maxwell_fluxes cache1[6], cache2[6];				// cache of maxwell fluxes for all six walls
	class Write_buffer_before;
	class Read_buffer_after;
	class Write_buffer_after;
	class Init_walls;
};

template<class Limiter>
class TVD_scheme : public Difference_scheme {
public:
	void read_buffer_before (Box*) { }
	void write_buffer_before (Box* box) { pimpl->write_buffer_before (box); }
	void read_buffer_after (Box* box) { pimpl->read_buffer_after (box); }
	void write_buffer_after (Box* box)  { pimpl->write_buffer_after (box); }
	void scheme (Box*);
	Buffer* create_buffer (int c, Box* b, Side is, Side at) { return new TVD_buffer (c, b, is, at); }
	void init_boundary (Box* box) { pimpl->init_boundary (box); }
	TVD_scheme () : pimpl (new TVD_scheme_impl) { }
	void info (Printer*);
	~TVD_scheme () { delete pimpl; }
private:
	TVD_scheme_impl* pimpl;
	class Interior;
	class Interior_first;
	class Interior_last;
	class Boundary;
};

#include "tvd_scheme_impl.h"
#endif // TVD_SCHEME_H
