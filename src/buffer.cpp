#include <mpi.h>
#include "box.h"

Buffer* Buffer::twin_buffer (Box* this_box) const
{
	if (twin) return twin;
	for (Box::MPI_Buffer::iterator pbuf = that_box->MPI_buffer ().begin (); pbuf != that_box->MPI_buffer ().end (); pbuf++)
		if ((*pbuf)->that_box == this_box && (*pbuf)->that_side == this_side)
			return *pbuf;
	assert (false); return 0;
}
