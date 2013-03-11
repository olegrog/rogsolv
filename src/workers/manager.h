#ifndef MANAGER_H
#define MANAGER_H

#include "../base/auxiliary.h"
#include "../containers/box.h"
#include "schemer.h"
#include "writer.h"
#include "ci_grider.h"
#include "printer.h"
#include "timer.h"

/** Singleton Manager **/
class Manager {
	Writer* writer;
	CI_grider* ci_grider;
	Timer* timer;
	Schemer* schemer;
	const Printer& printer;
	Box::Space space;														// space axis set
	Int_vect size;															// size of construction
	Boxes boxes;															// set of all sorting boxes
	Set_of_boxes boxes_;													// raw set of boxes
	int rank;																// current MPI_rank
	void MPI_distribute ();													// disrtribute boxes between MPI nodes
	void divide_box (Box*, Axis, int parts);								// divide box along axis
	void join_boundaries (Box*, Side, Box*, Side);							// connect boxes boundaries
	template <class T>
		void for_all_boxes (T) const;										// do smth for all boxes
	template <class T>
		void for_all_boxes (T* obj, void (T::*) (Box*)) const;				// do smth for all boxes
	void transfer () const;													// left part of Boltzmann equation
	void collision_integral (real) const;									// right part of Boltzmann equation
public:
	/** initializing procedures **/
	void set_grids (														// set up velocity&space grids parameters
		real cutting_velocity,
		real knudsen_number, 
		int character_size, 
		Box::Space axis_space = Box::Space ({XX, YY, ZZ})
	);
	void set_workers (														// employ special workers
		Writer*,
		CI_grider*,
		Timer*,
		Schemer*
	);
	void add_tracing_f (Box*, const Int_vect&);							// adding point for tracing distribution function

	/** constructing procedures **/
	void link_boxes (Box*, Box*, Axis, Int_vect offset);					// linking boxes along axis with offset
	void link_boxes (Box*, Side, Box*, Side);								// complex linking
	void alone_box (Box* box) { boxes_.insert (box); }						// used when no one link_boxes called

	/** simulating procedures **/
	void init_model ();														// preparing construction for calculating
	void iterate () const;													// iterations of difference equations

	static Manager& instance () { static Manager m; return m; }
public:
	/** procedures for workers **/
	const Writer& get_writer () const { return *writer; }
	Boxes& get_boxes () { return boxes; }
private:
	Manager ();
	Manager (const Manager&);
	const Manager& operator= (const Manager&);
	~Manager ();
};

inline Manager& manager () { return Manager::instance (); }

template <class T>
inline void Manager::for_all_boxes (T fun) const
{
	for (BI pbox = boxes.begin (); pbox != boxes.end (); ++pbox) {
		Box* box = *pbox;
		if (rank != box->MPI_rank ()) continue;
		fun (box);
	}
}

template <class T>
inline void Manager::for_all_boxes (T* obj, void (T::* fun) (Box*)) const
{
	for (BI pbox = boxes.begin (); pbox != boxes.end (); ++pbox) {
		Box* box = *pbox;
		if (rank != box->MPI_rank ()) continue;
		(obj->*fun) (box);
	}
}

#endif // MANAGER_H