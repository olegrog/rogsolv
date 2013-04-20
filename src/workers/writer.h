#ifndef WRITER_H
#define WRITER_H

#include <map>
#include <memory>
#include "../containers/box.h"

typedef std::multimap<Box*, Int_vect> Proj_points;

namespace Writers {
	enum Files_format {BIN, TXT, XML};
	enum Cell_type {SOLID, GAS};
	typedef std::set<std::size_t> Vertices;

	struct Map {
		Features* features;
		Cell_type type;
		Vertices vertices;
		Map () : type (SOLID) {}
	};
	
	void write_param (Files_format, std::ofstream&, real);
	void write_param (Files_format, std::ofstream&, const Real_vect&);
	void read_param (Files_format, std::ifstream&, real&);
	void read_param (Files_format, std::ifstream&, Real_vect&);
}

class Writer {
protected:
	const Boxes& boxes;
	const int MPI_rank;
	Int_vect size;
	Proj_points points;
	Matrix<Writers::Map>* map;
	void macroparameters () const;
	void build_map ();
public:
	Writer ();
	bool save_f (int time) const;
	bool load_f (int& time) const;
	bool load_macro () const;
	void set_size (Int_vect s) { size = s; }
	virtual void prepare_files () = 0;
	virtual bool write_result (int time) const = 0;
	void add_point (Box*, const Int_vect&);
	void write_f (int time) const;
	virtual void info () const = 0;
	virtual ~Writer ();
};

#endif // WRITER_H
