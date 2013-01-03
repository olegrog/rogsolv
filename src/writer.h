#ifndef WRITER_H
#define WRITER_H

#include <map>
#include <memory>
#include "box.h"

enum Files_format {BIN, TXT, XML};
enum Cell_type {SOLID, GAS};

struct Map {
	const Features* features;
	Cell_type type;
	Map () : type (SOLID) { }
};

typedef std::multimap<Box*, Int_vect> Proj_points;

class Writer {
protected:
	void macroparameters ();
	Int_vect size;
	const Boxes& boxes;
	int MPI_rank;
	Proj_points points;
public:
	bool save_f (int time);
	bool load_f (int& time);
	Writer (const Boxes& b, int rank) : boxes (b), MPI_rank (rank) { }
	void set_size (Int_vect s) { size = s; }
	virtual void prepare_files () = 0;
	virtual bool write_result (int time) = 0;
	void add_point (Box*, const Int_vect&);
	void write_f (int time);
	virtual ~Writer () {}
};

namespace Writers {
	typedef std::unique_ptr<Writer> pWriter;
	struct Writer_creator {
		virtual pWriter create (const Boxes&, int) = 0;
		virtual ~Writer_creator () {}
	};
	void write_param (Files_format, std::ofstream&, real);
	void write_param (Files_format, std::ofstream&, Real_vect);
}

#endif
