#ifndef PARAVIEW_H
#define PARAVIEW_H

#include "writer.h"

class Writer_paraview : public Writer {
	int points, cells;
	void prepare_files ();
public:
	Writer_paraview (const Boxes& b, int rank) : Writer (b, rank) { }
	bool write_result (int);
};

namespace Writers {
	class ParaView : public Writer_creator {
	public:
		pWriter create (const Boxes& b, int r) { return pWriter (new Writer_paraview (b,r)); }
	};
}

#endif
