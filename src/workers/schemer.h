#ifndef SCHEMER_H
#define SCHEMER_H

#include <map>

#include "../containers/walls.h"
#include "../containers/box.h"

// this class is a part of the third abstract level containing difference scheme
// it also provide exchange data between boxes
class Schemer {
public:	
	typedef std::map<Wall_simple::Params, Half_grid*> Maxwell_fluxes;	// cache of maxwell fluxes with the unit of density for all usable temperatures
	void next_layer (Box*);												// update the distribution function with calculated fluxes
	void MPI_exchange_before (Box*);									// exchanging buffers with other boxes
	void MPI_exchange_after (Box*);
	virtual void scheme (Box*) = 0;
	virtual void write_buffer_before (Box*) = 0;
	virtual void write_buffer_after (Box*) = 0;
	virtual void read_buffer_before (Box*) = 0;
	virtual void read_buffer_after (Box*) = 0;
	virtual Buffer* create_buffer (int c, Box* b, Side is, Side at) = 0;
	virtual void init_boundary (Box*) = 0;
	virtual void info () = 0;
	virtual ~Schemer () { }
private:
	void MPI_exchange (Box*, void (Buffer::*) (int), void (Buffer::*) (Box*));	// exchanging buffers with other boxes
};

#endif // SCHEMER_H
