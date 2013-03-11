#ifndef CI_GRIDER_H
#define CI_GRIDER_H

#include "../base/auxiliary.h"
#include "../ci/ci.hpp"

class CI_grider {
protected:
	real mass;
	ci::Particle particle;
public:
	virtual void generate (real time_step) = 0;
	virtual void info () = 0;
	virtual ~CI_grider () {}
public:
	void set_molecule (real m, const ci::Particle& p) {
		mass = m; particle = p;
	}
};

namespace CI_griders {
	constexpr int dimension = 10;									// dimension for ci::gen
	using Point = std::array<real, dimension>;						// integrate point
}

#endif // CI_GRIDER_H
