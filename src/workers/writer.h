#ifndef WRITER_H
#define WRITER_H

#include <map>
#include <memory>
#include "../containers/box.h"

typedef std::multimap<Box*, Int_vect> Proj_points;

class Writer {
protected:
	const Boxes& boxes;
	const int MPI_rank;
	Int_vect size;
	Proj_points points;
	void macroparameters () const;
public:
	Writer ();
	bool save_f (int time) const;
	bool load_f (int& time) const;
	void set_size (Int_vect s) { size = s; }
	virtual void prepare_files () = 0;
	virtual bool write_result (int time) const = 0;
	void add_point (Box*, const Int_vect&);
	void write_f (int time) const;
	virtual void info () const = 0;
	virtual ~Writer () {}
};

namespace Writers {
	enum Files_format {BIN, TXT, XML};
	enum Cell_type {SOLID, GAS};

	struct Map {
		const Features* features;
		Cell_type type;
		Map () : type (SOLID) {}
	};

	void write_param (Files_format, std::ofstream&, real);
	void write_param (Files_format, std::ofstream&, Real_vect);
}

#endif // WRITER_H
