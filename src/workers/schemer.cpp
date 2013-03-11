#include "schemer.h"

void Schemer::next_layer (Box* box)
{
	mapper ().set_side ();
	update (box->f->all (), box->F, box->tau/box->H, box->space ());
}

void Schemer::MPI_exchange (Box* this_box, void (Buffer::* MPI_swap) (int), void (Buffer::* swap) (Box*))
{
	for (Box::MPI_Buffer::iterator pbuf = this_box->MPI_buffer ().begin (); pbuf != this_box->MPI_buffer ().end (); ++pbuf) {
		Buffer* buf = *pbuf;
		if (buf->that_box->MPI_rank ()  != this_box->MPI_rank ())
			(buf->*MPI_swap) (buf->that_box->MPI_rank ());							// swaping packets using MPI
		else	{
			if (buf->that_box < this_box) continue;									// against double-swap
			if (buf->that_box == this_box && buf->this_side < buf->that_side) continue;	// against double-swap
			(buf->*swap) (this_box);												// swaping pointers to buffers
		}
	}
}

void Schemer::MPI_exchange_before (Box* box)
{
	MPI_exchange (box, &Buffer::MPI_swap_before, &Buffer::swap_before);	
}

void Schemer::MPI_exchange_after (Box* box)
{
	MPI_exchange (box, &Buffer::MPI_swap_after, &Buffer::swap_after);	
}
