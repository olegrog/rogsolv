#ifndef BKVIEWER_H
#define BKVIEWER_H

#include "writer.h"

class Writer_bkviewer : public Writer {
	Matrix<Map>* map;
	void write_static (std::ofstream& , Int_vect, int);
	void build_map ();
	void prepare_files ();
public:
	Writer_bkviewer (const Boxes& b, int rank) : Writer (b, rank) { }
	bool write_result (int);
	~Writer_bkviewer ();
};

namespace Writers {
	class BKViewer : public Writer_creator {
	public:
		pWriter create (const Boxes& b, int r) { return pWriter (new Writer_bkviewer (b,r)); }
	};
}


#endif
