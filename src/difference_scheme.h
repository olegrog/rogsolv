#ifndef DIFFERENCE_SCHEME_H
#define DIFFERENCE_SCHEME_H

#include <map>

#include "walls.h"
//#include "printer.h"
#include "box.h"
#include "buffer.h"

class Printer;
// this class is a part of the third abstract level containing difference scheme
// it also provide exchange data between boxes
class Difference_scheme {
public:	
	typedef std::map<Wall_simple::Params, Half_grid*> Maxwell_fluxes;	// cache of maxwell fluxes with the unit of density for all usable temperatures
	void next_layer (Box*);								// update the distribution function with calculated fluxes
	void MPI_exchange_before (Box*);					// exchanging buffers with other boxes
	void MPI_exchange_after (Box*);
	virtual void scheme (Box*) = 0;
	virtual void write_buffer_before (Box*) = 0;
	virtual void write_buffer_after (Box*) = 0;
	virtual void read_buffer_before (Box*) = 0;
	virtual void read_buffer_after (Box*) = 0;
	virtual Buffer* create_buffer (int c, Box* b, Side is, Side at) = 0;
	virtual void init_boundary (Box*) = 0;
	virtual void info (Printer*) = 0;
	virtual ~Difference_scheme () { }
private:
	void MPI_exchange (Box*, void (Buffer::*) (int), void (Buffer::*) (Box*));	// exchanging buffers with other boxes
};

#endif // DIFFERENCE_SCHEME_H
