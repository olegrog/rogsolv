#ifndef MAPPER_H
#define MAPPER_H

#include <valarray>
#include <cassert>

#include "vector3.h"

// Singleton
class Mapper {
	std::valarray<Int_vect> *ind;											// index->vector
	std::valarray<int> *xyz, *all, *halves[6], *rhalves[6];						// set of indexations
	const std::valarray<int>* index_, *rindex_;								// current indexes of map
	int size_;																// current size of map
	int volume_;															// volume of vel_grid
public:
	int operator() (Int_vect p) const { return (*xyz)[p.x+p.y*(2*R)+p.z*(2*R)*(2*R)]; } // for debug & write_f
	int operator() (int x, int y, int z) const { return (*this) (Int_vect (x, y, z)); }		// for CI module
	const Int_vect& operator[] (int i) const { return (*ind)[i]; }
	const std::valarray<int>& half_ball (Side side) const { return *halves[side]; }	// half grid indexation
	const std::valarray<int>& all_ball () const { return *all; }					// all grid indexation
	int index (int i) const { return (*index_)[i]; }								// indexation
	int rindex (int i) const { return (*rindex_)[i]; }								// anti-indexation
	int size () const { return size_; }
	int volume () const { return volume_; }
	int radius () const { return R; }
	void set_side (Side);											// switch to half grid index
	void set_side ();														// switch to all grid index
	
	static Mapper& instance ();
	static void set_radius (int radius) { assert (!R); R = radius; }
private:
	static int R;															// radius of vel_grid
	Mapper ();
	Mapper (const Mapper&);
	const Mapper& operator= (const Mapper&);
	~Mapper ();
};

inline Mapper& mapper () { return Mapper::instance (); }						// short alias

#endif
