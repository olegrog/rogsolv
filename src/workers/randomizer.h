#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <random>
#include <functional>

#include "../base/auxiliary.h"

/** Singleton Randomizer **/
class Randomizer {
	std::default_random_engine generator;
	std::uniform_real_distribution<real> distrib;
public:
	real operator () () const {
		static auto dice (std::bind (distrib, generator));
		return dice ();
	}
	static const Randomizer& instance () { static Randomizer r; return r; }
private:
	Randomizer () : generator (1000), distrib (0, 1) {}
	Randomizer (const Randomizer&);
	const Randomizer& operator= (const Randomizer&);
	~Randomizer () {}
};

inline const Randomizer& randomizer () { return Randomizer::instance (); }

#endif // RANDOMIZER_H