#ifndef BUFFER_H
#define BUFFER_H

#include <mpi.h>
#include <valarray>

#include "vel_grid.h"

class Box;

class Buffer {															// buffer for exchange boxes boundaries
public:
	Box* const that_box;
	const Side this_side, that_side;
	Buffer (int cap, Box* b, Side this_s, Side that_s)
		: that_box (b), this_side (this_s), that_side (that_s), capacity (cap), twin (0) { }
	virtual ~Buffer () { }
	virtual void init_buffer () = 0;										// delayed allocating memory
	virtual void swap_before (Box* dest_box) = 0;							// exchange packets with swaping pointers
	virtual void swap_after (Box* dest_box) = 0;
	virtual void MPI_swap_before (int dest_rank) = 0;						// exchange packets using MPI
	virtual void MPI_swap_after (int dest_rank) = 0;
	virtual Half_grid& flux (int) = 0;										// return flux from wall
	virtual Vel_grid& grid (int) = 0;										// return adjacent grid (for TVD)
	Buffer* twin_buffer (Box*) const;										// return partner for swaping
protected:
	template<class T> void MPI_swap (int dest_rank, T& packages, size_t);	// realization
	const int capacity;														// capacity of buffer = number of common cells
private:
	Buffer* twin;															// pointer to twin buffer
private:																	// prohibited constructors
	Buffer ();
	Buffer (const Buffer&);
	const Buffer& operator= (const Buffer&);
};

template<class T>
void Buffer::MPI_swap (int dest, T& packages, size_t size)
{
	MPI_Datatype datatype;
	MPI_Status status;
	if (sizeof (real) == sizeof (double)) datatype = MPI_DOUBLE;
	if (sizeof (real) == sizeof (float)) datatype = MPI_FLOAT;
	for (int i=0; i<capacity; i++)
		MPI_Sendrecv_replace (packages[i].raw_data (), size, datatype, dest, 0, dest, 0, MPI_COMM_WORLD, &status);
}

#endif
