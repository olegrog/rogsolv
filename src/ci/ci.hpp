// нет поддержки потенциала Ленарда-Джонса

#ifndef _CI_H_
#define _CI_H_

#include <vector>
#include <cstddef>

namespace ci {

	enum Symmetry {
		NO_SYMM	= 0,	// нет симметрии
		Z_SYMM	= 1,	// симметрия по оси z
		YZ_SYMM = 2		// симметрия по осям y, z 
	};

	struct Particle {
		double d;
	};
	struct LJParticle : public Particle {
		double e;
	};

	template <typename Map>
		int gen(double tt, int c_nd, int nk_rad1, int nk_rad2, 
			const Map& xyz2i1, const Map& xyz2i2,
			double a, double m1, double m2, const Particle& p1, const Particle& p2);

	template <typename F> void iter(F& f1, F& f2);

	void finalize();

	class Potential {
		public:
			virtual double theta(const Particle& p1, const Particle& p2, double b, double g) const = 0;
			virtual double bMax(const Particle& p1, const Particle& p2) const = 0;
	};

	class HSPotential : public Potential {
		public:
			double theta(const Particle& p1, const Particle& p2, double b, double g) const;
			double bMax(const Particle& p1, const Particle& p2) const;
	};

	class LJPotential : public Potential {
		public:
			LJPotential(double b_extension = 2.5, size_t b_size = 50,
					double g_step = 0.1, double r_step = 2.5e-4);
			double theta(const Particle& p1, const Particle& p2, double b, double g) const;
			double bMax(const Particle& p1, const Particle& p2) const;
		private:
			double b_extension;
			size_t b_size;
			double g_extension, g_step;
			mutable size_t g_size; 
			double r_step;
			mutable std::vector<double> gb2theta;

			double gbToThetaCalc(double g, double b) const;
			void fillGbToTheta(size_t from, size_t to) const;
			void extentGbToTheta(double g) const;
	};

	void init(const Potential* p, Symmetry s); 

}

#include "ci_impl.hpp"

#endif
