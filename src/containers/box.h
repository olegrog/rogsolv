#ifndef BOX_H
#define BOX_H

#include <map>
#include <vector>
#include <set>
#include <cassert>

#include "vel_grid.h"
#include "matrix.h"
#include "buffer.h"
#include "walls.h"
#include "../base/init_cond.h"

struct Features {															// list of calculating macroparameters
	real dens, temp;														// scalar parameters
	Real_vect flow, qflow, press, shear;									// vector parameters
	Features (real d, real t, Real_vect f, Real_vect qf, Real_vect pr, Real_vect sh) : 
		dens (d), temp (t), flow (f), qflow (qf), press (pr), shear (sh) { }
	Features () { }
};

// Box -- container of the third abstract level
// it works only with matrixes and return macroparameters
class Box {
public:
	struct Box_less {														// for sorting set of Boxes
		bool operator() (const Box* box1 , const Box* box2) const
		{																
			if (box1->size ().vol () != box2->size ().vol ()) 
				return box1->size ().vol () < box2->size ().vol ();		// sorting by box volume 
			for (int i=0; i<3; i++)
				if (box1->coord ()[i] != box2->coord ()[i]) 
					return box1->coord ()[i] < box2->coord ()[i];			// sorting by box coordinate
			assert (false);
			return false;
		}
	};
	struct Buffer_less {													// for sorting MPI_buffer
		bool operator() (const Buffer* buf1, const Buffer* buf2) const
		{
			if (buf1->that_box != buf2->that_box)							// this guarantees lack of MPI blocking
				return Box_less () (buf1->that_box, buf2->that_box);		// always send to boxes sorting by Box_less 
			return buf1 < buf2;												// for complex linking
		}
	};
	typedef std::set<Buffer*, Buffer_less> MPI_Buffer;						// list of buffers for exchanging with other boxes
	typedef std::set<Axis> Space;

	static real H;															// size of spatial grid in terms of mean free molecular path
	static real tau;														// time step in terms of mean free molecular path time
	MPI_Buffer MPI_buffer_;
	Matrix<Vel_grid>* f;													// distribution function
	Matrix<Vel_grid>* F[3];													// fluxes in all directions
	Matrix<Features>* features;												// all macroparameters
	Matrix<Wall*>* wall[6];													// all walls
	static Init_cond* init_cond;
private:
	friend class Manager;
	Space space_;															// space axis set
	const Int_vect size_;													// box size 
	Int_vect coord_;														// box position in the construction 
	real init_temp, init_dens;												// initial temperature and density
	int MPI_rank_;															// MPI_rank of process which calculate this box
public:
	void add_buffer (Buffer* buf) { MPI_buffer_.insert (buf); }				// add buffer
	void remove_buffer (Buffer* buf) { MPI_buffer_.erase (buf); delete buf; }		// remove and delete buffer
	void init_distribution ();												// creating initial distribution
	const Int_vect& size () const { return size_; }
	const Int_vect& coord () const { return coord_; }
	int MPI_rank () const { return MPI_rank_; }
	const Space& space () const { return space_; }
	void set_space (const Space& sp) { space_ = sp; }
	const MPI_Buffer& MPI_buffer () const { return MPI_buffer_; }
	Box (Int_vect box_size);
	~Box ();
	void set_mirror (Side);													// set mirror boundary
	void set_maxwell (Side, real temp, real dens);							// set free Maxwell distribution
	void set_interface (Side, real temp, real dens, Real_vect speed);		// set interface boundary
	void set_simple (Side, const Simple_bound&);							// set simple boundary
};

// calculating all macroparameters
class Macroparameters {
public:
	void operator () (Box*);
};

// defining new type of sort, because sort by address {Box*} depend on MPI_node
typedef std::set<Box*, Box::Box_less> Boxes;
typedef std::set<Box*> Set_of_boxes;
typedef Boxes::const_iterator BI;
typedef Set_of_boxes::const_iterator BI_;

#endif
