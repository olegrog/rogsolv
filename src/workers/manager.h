//#include <memory>

#include "../containers/box.h"
#include "../schemes/difference_scheme.h"
#include "../writers/writer.h"
#include "printer.h"
#include "timer.h"
#include "../base/auxiliary.h"

class Manager {
	Box::Space space;													// space axis set
	Int_vect size;														// size of construction
	Boxes boxes;														// set of all sorting boxes
	Boxes_ boxes_;														// raw set of boxes
	int MPI_rank;														// current MPI_rank
	Printer* printer;
	Writers::pWriter writer;
	Timer* timer;
	Difference_scheme* scheme;
	void MPI_distribute ();												// disrtribute boxes between MPI nodes
	void divide_box (Box*, Axis, int parts);							// divide box along axis
	void join_boundaries (Box*, Side, Box*, Side);						// connect boxes boundaries
	template <class T> void for_all_boxes (T);							// do smth for all boxes
	template <class T> void for_all_boxes (T* obj, void (T::*) (Box*));	// do smth for all boxes
	void transfer ();													// left part of Bolzmann equation
	void collision_integral (real);										// right part of Bolzmann equation
public:
	void set_timer (int end,
					int log, 
					int macroparameters=0, 
					int cache_f=0);
	void set_grid (	real cutting_velocity,
					real knudsen_number, 
					int character_size, 
					Box::Space axis_space = Box::Space ({XX, YY, ZZ}));
	void set_scheme (Difference_scheme*);
	void link_boxes (Box*, Box*, Axis, Int_vect offset);				// linking boxes along axis with offset
	void link_boxes (Box*, Side, Box*, Side);							// complex linking
	void alone_box (Box* box) { boxes_.insert (box); }					// used when no one link_boxes called
	void init_model ();													// preparing construction for calculating
	void iterate ();													// iterations of difference equations
	void macroparameters ();											// calculating all macroparameters
	void write_result (int time) { writer->write_result (time); }		// write them into files
	void add_tracing_f (Box*, const Int_vect&);						// adding point for tracing distribution function
	Manager (Writers::Writer_creator*);
	~Manager ();
};

template <class T>
inline void Manager::for_all_boxes (T fun)
{
	for (BI pbox = boxes.begin (); pbox != boxes.end (); ++pbox) {
		Box* box = *pbox;
		if (MPI_rank != box->MPI_rank ()) continue;
		fun (box);
	}
}

template <class T>
inline void Manager::for_all_boxes (T* obj, void (T::* fun) (Box*))
{
	for (BI pbox = boxes.begin (); pbox != boxes.end (); ++pbox) {
		Box* box = *pbox;
		if (MPI_rank != box->MPI_rank ()) continue;
		(obj->*fun) (box);
	}
}

