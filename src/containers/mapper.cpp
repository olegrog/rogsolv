#include "mapper.h"

int Mapper::R = 0;

Mapper& Mapper::instance () 
{
	assert (R);
	static Mapper mapper;
	return mapper;
}

inline real sqr (real x) { return x*x; }

#define FOR_BALL						\
	for (p.z=0; p.z<2*R; p.z++)				\
		for (p.y=0; p.y<2*R; p.y++)			\
			for (p.x=0; p.x<2*R; p.x++)		\
				if (sqr (.5+p.x-R) + sqr (.5+p.y-R) + sqr (.5+p.z-R) <= sqr (R))


Mapper::Mapper ()
{
	volume_ = 0;
	xyz = new std::valarray<int> ((2*R)*(2*R)*(2*R));
	Int_vect p;
	FOR_BALL
		(*xyz)[p.x+p.y*(2*R)+p.z*(2*R)*(2*R)] = volume_++;
	else (*xyz)[p.x+p.y*(2*R)+p.z*(2*R)*(2*R)] = -1;
	ind = new std::valarray<Int_vect> (volume());
	all = new std::valarray<int> (volume ());
	for (int s=0; s<6; s++) {
		halves[s] = new std::valarray<int> (volume()/2);
		rhalves[s] = new std::valarray<int> (volume());
	}
	int index[7];
	for (int i=0; i<7; i++) index[i]=0;
	FOR_BALL {
		(*all)[index[6]] = index[6];
		(*ind)[index[6]++] = p;
		for (int s=0; s<6; s++) {
			Axis axis = Axis (s%3); Direction direction = Direction (s%2);
			(*rhalves[s])[(*xyz)[p.x+p.y*(2*R)+p.z*(2*R)*(2*R)]] = -1;			// for debug
			if (direction == FORWARD && p[axis] < -.5+R) continue;
			if (direction == BACKWARD && p[axis] > -.5+R) continue;
			(*rhalves[s])[(*xyz)[p.x+p.y*(2*R)+p.z*(2*R)*(2*R)]] = index[s];
			(*halves[s])[index[s]++] = (*xyz)[p.x+p.y*(2*R)+p.z*(2*R)*(2*R)];
		}
	}
}
#undef FOR_BALL

Mapper::~Mapper ()
{
	delete ind; delete xyz; delete all;
	for (int s=0; s<6; s++) {
		delete halves[s]; delete rhalves[s];
	}
}

// for_debug
std::ostream& operator<< (std::ostream& str, Side side)
{
	switch (side) {
		case RIGHT: str << " RIGHT "; break;
		case LEFT: str << " LEFT "; break;
		case BOTTOM: str << " BOTTOM "; break;
		case TOP: str << " TOP "; break;
		case FRONT: str << " FRONT "; break;
		case BACK: str << " BACK "; break;
	}
	return str;
}

void Mapper::set_side (Side side)
{
//	std::cout << "set_side" << side << '\n';
	size_ = volume ()/2;
	index_ = halves[side];
	rindex_ = rhalves[side];
}

void Mapper::set_side ()
{
	size_= volume ();
	index_ = all;
}
