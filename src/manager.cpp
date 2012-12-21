#include <functional>
#include <stack>
#include <mpi.h>
#include <stdexcept>

#include "manager.h"
#include "ci.h"

struct Join_boundaries {
	Buffer& buf;
	const Int_vect size;
	const Int_vect offset;
	Join_boundaries (Buffer& b, Int_vect s, Int_vect o) : buf (b), size (s), offset (o) { }
	int index (Int_vect p, Int_vect s, Axis ax) { return p[(ax+1)%3]+s[(ax+1)%3]*p[(ax+2)%3]; }
	void operator () (Wall*& wall, Int_vect coord) {
		Axis this_axis = axis (buf.this_side);
		Axis that_axis = axis (buf.that_side);
		Int_vect coord2 = rotate (offset+coord, buf.that_side, buf.this_side) - rotate (offset, buf.that_side, buf.this_side);
		Int_vect size2 = abs (rotate (size, buf.that_side, buf.this_side));
		if ((direction (buf.this_side)+direction (buf.that_side))%2)
			coord2[this_axis] += size2[this_axis]-1;
		int read_index = index (coord, size, this_axis);
		int write_index = index (coord2, size2, that_axis);
		if (wall) delete wall;
		wall = new Wall_box (buf, read_index, write_index);
	}
};

void Manager::join_boundaries (Box* box1, Side side1, Box* box2, Side side2)
{
	Int_vect beg1, end1;													// begin and end of common area in box1
	Int_vect beg2, end2;													// begin and end of common area in box2
	Int_vect b1, b2;														// need for Join_boundaries
	if (axis (side1) == axis (side2)) {
		const Axis ax = axis (side1);
		if (box1->coord () [ax] > box2->coord () [ax]) { std::swap (box1, box2); std::swap (side1, side2); }
		Int_vect offset = box2->coord () - box1->coord ();
		beg1 = max (Int_vect (), offset);
		end1 = min (box1->size (), box2->size ()+offset);
		beg2 = beg1 - offset;
		end2 = end1 - offset;
		if (box1->coord() [ax]+box1->size ()[ax] != box2->coord ()[ax])  return;	// assert (distance==0)
	} else {
		const Axis ax1 = ::axis (side1);
		const Axis ax2 = ::axis (side2);
		// coordinates of necessary planes
		b1 = box1->coord (); b1[ax1] += direction (side1)*box1->size ()[ax1];
		Int_vect e1 = box1->coord () + box1->size (); e1[ax1] -= !direction (side1)*box1->size ()[ax1];
		b2 = box2->coord (); b2[ax2] += direction (side2)*box2->size ()[ax2];
		Int_vect e2 = box2->coord () + box2->size (); e2[ax2] -= !direction (side2)*box2->size ()[ax2];
		Int_vect centre; centre[ax1] = b1[ax1]; centre[ax2] = b2[ax2];			// rotate fixed point
		bool swaping = false;
		if ((direction (side1)+direction (side2))%2) {
			std::swap (b2[ax1], e2[ax1]);
			std::swap (b2[ax2], e2[ax2]);
			swaping = true;
		}
		beg1 = max (b1, rotate (b2-centre, side1, side2)+centre);
		end1 = min (e1, rotate (e2-centre, side1, side2)+centre);
		beg2 = rotate (beg1-centre, side2, side1) + centre;
		end2 = rotate (end1-centre, side2, side1) + centre;
		if (swaping) { 													// back swaping
			std::swap (b2[ax1], e2[ax1]);
			std::swap (b2[ax2], e2[ax2]);
			std::swap (beg2[ax1], end2[ax1]);
			std::swap (beg2[ax2], end2[ax2]);
		}
		beg1 -= b1; end1 -= b1; beg2 -= b2; end2 -= b2;						// relavive coordinates
	}
	const Int_vect size = end1 - beg1;										// common boundary size
	const int area = size.area (axis (side1));									// common area between boxes
	if (area <= 0 || size.sum () - size[axis (side1)] <= 0)  return;					// assert (area > 0)
	Buffer* buf1 = scheme->create_buffer (area, box2, side1, side2);				// creating MPI_buffers
	Buffer* buf2 = scheme->create_buffer (area, box1, side2, side1);
	box1->add_buffer (buf1); box2->add_buffer (buf2);							// add them to boxes MPI_Buffer
	for_each_index (box1->wall[side1]->bound_rect (side1, beg1, end1),
		Join_boundaries (*buf1, end1-beg1, b1));								// joinind boxes
	for_each_index (box2->wall[side2]->bound_rect (side2, beg2, end2),
		Join_boundaries (*buf2, end2-beg2, b2));
}

// this function link boxes along {axis} with {offset} and add new boxes to the construction
// assert (box1->coord() [axis]+box1->size[axis] == box2->coord ()[axis])
void Manager::link_boxes (Box* box1, Box* box2, Axis axis, Int_vect offset)
{
	if (box1 == box2) throw std::logic_error ("You tried to link the same box.");
	if (box1->size ().vol () <= 0 || box2->size ().vol () <= 0)
		throw std::logic_error ("The box volume must be positive.");
	if (boxes_.empty ()) { boxes_.insert (box1); box1->coord_ = 0; }				// if it's first call of link_boxes
	Int_vect size = min (box1->size (), box2->size ()+offset)
				- max (Int_vect (), offset);									// common boundary size
	int area = size.area (axis);												// common area between boxes
	bool new1 = boxes_.find (box1) == boxes_.end ();							// looking for new_boxes in set<boxes>
	bool new2 = boxes_.find (box2) == boxes_.end ();
	if (new1 && new2)
		throw std::logic_error ("You tried to link two new boxes. You must join all boxes sequentially.");
	if (new1 || new2) {
		if (area <= 0 || size.sum () - size[axis] <= 0)							// (a*b <=0 || (a+b)<=0) <=> (a<=0 || b<=0)
			throw std::logic_error ("Illegal linking of boxes. The boxes must touch each other.");
		Box* new_box = box1; Box* old_box = box2; int sign = 1;				// if box1 is new
		if (new2) {														// if box2 is new
			std::swap (new_box, old_box);
			sign = -1;
		}
		boxes_.insert (new_box);
		offset[axis] = box1->size ()[axis];
		new_box->coord_ = old_box->coord () - sign*offset;						// setting coordinates to new box
	}
	// now the both boxes include in construction and have coordinates
	join_boundaries (box1, side (axis, FORWARD), box2, side (axis, BACKWARD));
}

void Manager::link_boxes (Box* box1, Side side1, Box* box2, Side side2)
{
	if (boxes_.find (box1) == boxes_.end () || boxes_.find (box2) == boxes_.end ())
		throw std::logic_error ("Illegal complex link. Both boxes must be linked before.");
	join_boundaries (box1, side1, box2, side2);	
}

void Manager::transfer ()
{
	timer->nick ("exchange");
	for_all_boxes (scheme, &Difference_scheme::write_buffer_before);
	for_all_boxes (scheme, &Difference_scheme::MPI_exchange_before);
	for_all_boxes (scheme, &Difference_scheme::read_buffer_before);
	timer->nick ("transfer");
	for_all_boxes (scheme, &Difference_scheme::scheme);
	timer->nick ("exchange");
	for_all_boxes (scheme, &Difference_scheme::write_buffer_after);
	for_all_boxes (scheme, &Difference_scheme::MPI_exchange_after);
	for_all_boxes (scheme, &Difference_scheme::read_buffer_after);
	timer->nick ("transfer");
	for_all_boxes (scheme, &Difference_scheme::next_layer);
}

void Manager::collision_integral (real time)
{
	timer->nick ("integral");
	ci_gen (time, mapper ());
	for_all_boxes (Collision_integral ());
}

void Manager::iterate ()
{
	while (true) {
		if (timer->begin ()) break;											// true is the end
		transfer ();
		collision_integral (2*Box::tau);
		transfer ();
		timer->end ();
	}
}

void Manager::set_timer (int finish, int log, int macro, int cache)
{
	timer->set_parameters (finish, log, macro, cache);
}

void Manager::set_grid (real cut, real knud, int ch_size, Box::Space sp)
{
	assert (scheme);
	space = sp;
	int R = mapper ().radius ();
	Box::H = 1./knud/ch_size;
	Box::tau = 0*Box::H/sqrt (3)/cut;
	Vel_grid::set_cut (cut);
	printer->title ("Grid parameters");
	printer->var ("Knudsen number", knud);
	printer->var ("Coordinate step", Box::H);
	printer->var ("Time step", 2*Box::tau);
	printer->title ("Vel_grid parameters");
	printer->var ("Velocity grid radius", R);
	printer->var ("Cutting velocity", cut);
	printer->title ("Collision_integral");
	printer->var ("Number of korobov points", get_power ());
	int potential = HS_POTENTIAL;
	srand (1000);
	ci_init (R, cut, potential, NO_SYMM);
	switch (potential) {
	case HS_POTENTIAL:
		printer->var ("Molecular potential", "Hard sphere"); break;
	case LJ_POTENTIAL:
		printer->var ("Molecular potential", "Lennard-Jones"); break;
	}
}

void Manager::set_scheme (Difference_scheme* s)
{
	assert (mapper ().radius ());
	printer->title ("Difference scheme");
	scheme=s; scheme->info (printer);
}

Manager::Manager (Writers::Writer_creator* wc)
{
	MPI_Init (0, 0);
	MPI_Comm_rank (MPI_COMM_WORLD, &MPI_rank);
	printer = new Printer (MPI_rank);
	writer = std::move (wc->create(boxes, MPI_rank));
	delete wc;
	timer = new Timer (printer, writer.get(), boxes);
}

Manager::~Manager ()
{
	ci_finalize ();
	delete timer; delete printer; delete scheme;
	for (BI pbox = boxes.begin (); pbox != boxes.end (); ++pbox)
		delete *pbox;
	delete Box::init_cond;
	MPI_Finalize ();
}

struct Transfer_wall {
	void operator () (Wall*& sour, Wall*& dest)
	{
		if (!dynamic_cast<Wall_box*> (sour) && sour !=0)								// except Wall_box
			std::swap (sour, dest);														// transfer wall to new box
	}
};

// check if the boxes touched inner way
bool touch_inside (const Box* box1, const Box* box2, Side side) {
	if (direction (side) == BACKWARD)
		return box1->coord ()[axis (side)] == box2->coord()[axis (side)];
	else
		return (box1->coord ()+box1->size()) [axis (side)] == (box2->coord()+box2->size())[axis (side)];
}

void Manager::divide_box (Box* box, Axis axis, int parts)
{
	boxes_.erase (box);													// remove old box from the construction
	std::vector<Box*> new_box (parts);										// array of new boxes
	std::vector<Int_vect> size (parts);
	assert (box->size ()[axis] / parts);
	for (int i=0; i<parts; i++) {
		size[i] = box->size (); size[i][axis] /= parts;								// size of new boxes
		if (i < box->size ()[axis] % parts) size[i][axis]++;
	}
	Int_vect offset;	 													// begin of old box common boundary 
	for (int i=0; i<parts; i++) {
		new_box[i] = new Box (size[i]);			// create new box
		// transfering boundary conditions
		for (int s=0; s<6; s++) {
			Side side = Side (s);
			if (::axis (side) != axis)
				for_each (box->wall[side]->bound_rect (side, offset, offset+size[i]),
					new_box[i]->wall[side]->all (), Transfer_wall ());
		}
		offset[axis] += size[i][axis];
		if (i == 0) {
			for_each (box->wall[side (axis, BACKWARD)]->all (),
				new_box[i]->wall[side (axis, BACKWARD)]->all (), Transfer_wall ());
			new_box[0]->coord_ = box->coord ();
			boxes_.insert (new_box[0]);										// insert first box into construction
		} else {
			if (i == parts-1)
				for_each (box->wall[side (axis, FORWARD)]->all (), 
					new_box[i]->wall[side (axis, FORWARD)]->all (), Transfer_wall ());
			link_boxes (new_box[i-1], new_box[i], axis, 0);						// linking new boxes between each other
		}
		// connect new boxes with old box neighbours
		for (Box::MPI_Buffer::iterator pbuf = box->MPI_buffer ().begin (); pbuf != box->MPI_buffer ().end (); ++pbuf)
			if (box != (*pbuf)->that_box)
				join_boundaries (new_box[i], (*pbuf)->this_side, (*pbuf)->that_box, (*pbuf)->that_side);
	}
	// connect new boxes between each other in case of complex linking
	for (Box::MPI_Buffer::iterator pbuf = box->MPI_buffer ().begin (); pbuf != box->MPI_buffer ().end (); ++pbuf)
		if (box == (*pbuf)->that_box && *pbuf > (*pbuf)->twin_buffer (box))		// buffer comparison against double connecting
			for (int i=0; i<parts; i++)
				for (int j=0; j<parts; j++)
					if (touch_inside (box, new_box[i], (*pbuf)->this_side) && 
						touch_inside (box, new_box[j], (*pbuf)->that_side))
							join_boundaries (new_box[i], (*pbuf)->this_side, new_box[j], (*pbuf)->that_side);
	// remove old MPI_buffers
	for (Box::MPI_Buffer::iterator pbuf = box->MPI_buffer ().begin (); pbuf != box->MPI_buffer ().end (); ++pbuf)
		(*pbuf)->that_box->remove_buffer ((*pbuf)->twin_buffer (box));
	delete box;
}

void Manager::MPI_distribute ()
{
	int num;																// number of nodes
	MPI_Comm_size (MPI_COMM_WORLD, &num);
	printer->title ("MPI");
	printer->var ("Number of nodes", num);
	int cells = 0;															// total number of cells
	for (BI pbox = boxes_.begin (); pbox != boxes_.end (); ++pbox)
		cells += (*pbox)->size ().vol ();
	printer->var ("Total amount of cells", cells);
	const int average = cells/num;
	assert (average);
	printer->var ("Average to MPI_node", average);
	BI_ pbox = boxes_.begin ();
	while (true) {
		Box* box = *pbox;
		int n = box->size ().vol () / average;
		if (n > 1) {
			Axis axis = Axis (box->size ().max_axis ());
			divide_box (box, axis, n);
			pbox = boxes_.begin ();											// continue dividing till possible
			continue;
		}
		if (++pbox == boxes_.end ()) break;
	}
	for (BI_ pbox = boxes_.begin (); pbox != boxes_.end (); ++pbox)
		boxes.insert (*pbox);
	
	Boxes in_box;															// set<Box*> sorting by Box_less
	std::multimap<int, int> in_MPI;											// multimap of MPI_nodes<amount of own cells, MPI_rank>
	// filling multimaps
	for (BI pbox = boxes.begin (); pbox != boxes.end (); ++pbox) in_box.insert (*pbox);
	for (int i=0; i<num; i++) in_MPI.insert (std::make_pair (0, i));
	// distributing boxes from {in_box} to {in_MPI}
	for (Boxes::reverse_iterator pbox = in_box.rbegin (); pbox != in_box.rend (); ++pbox) {
		Box* box = *pbox;
		box->MPI_rank_ = in_MPI.begin ()->second;							// attaching the biggest box from "in_box" to the most capacious MPI_node
		int node_in_MPI = in_MPI.begin ()->first + box->size ().vol ();		// resulting volume in this MPI_node
		in_MPI.erase (in_MPI.begin ());										// removing this MPI_node from "in_MPI"
		if (node_in_MPI < average)
			in_MPI.insert (std::make_pair (node_in_MPI, box->MPI_rank ()));	// if MPI_node has free volume add it to "in_MPI"
	}
	
	printer->boxes (boxes);
	printer->MPI_ranks (boxes);
}

void Manager::init_model ()
{
	Int_vect dimension[2] = {0, 0};
	for (BI pbox = boxes_.begin (); pbox != boxes_.end (); ++pbox)
		dimension[1] = max (dimension[1], (*pbox)->coord ()+(*pbox)->size ()),
		dimension[0] = min (dimension[0], (*pbox)->coord ());					// overblowing size of construction
	for (BI pbox = boxes_.begin (); pbox != boxes_.end (); ++pbox)
		(*pbox)->coord_ -= dimension[0];									// alignment of coordinates
	size = dimension[1] - dimension[0];										// setting size of construction
	
	printer->title ("Model geometry");
	printer->var ("Size of construction", size);
		
	for (BI_ pbox = boxes_.begin (); pbox != boxes_.end (); ++pbox) {
		boxes.insert (*pbox);
		for (int ax=0; ax<3; ax++)
			if (!space.count (Axis (ax)) && (*pbox)->size () [ax] != 1)
				throw std::logic_error ("Illegal boxes size. Boxes must have unit size along free axis.");
	}
	printer->boxes (boxes);												// print preliminary set of boxes
	boxes.clear ();
	
	MPI_distribute ();

	for_all_boxes (std::bind2nd (std::mem_fun (&Box::set_space), space));
	for_all_boxes (std::mem_fun (&Box::init_distribution));
	for_all_boxes (scheme, &Difference_scheme::init_boundary);
	
	writer->set_size (size);
	writer->prepare_files ();

}

void Manager::add_tracing_f (Box* box, const Int_vect& point)
{
	writer->add_point (box, point);
}
