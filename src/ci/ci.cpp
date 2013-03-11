#include <cmath>
#include "ci.hpp"

namespace ci {

	int N_nu;
	std::vector<node_calc> nc;

	int ss[9];
	const double exact_hit = 1e-10;

	void finalize() {
		if (potential) delete potential;
	}

	Symmetry symm;
	const Potential* potential;
	void init (const Potential* p, Symmetry s) {
		potential = p;
		symm = s;
	}

	const V3d scatter(const V3d& x, double theta, double e) {

		double rxy = std::sqrt(sqr(x[0]) + sqr(x[1]));
		double r = std::sqrt(sqr(x));
		
		double se = std::sin(e);
		double ce = std::cos(e);
		double st = std::sin(theta);
		double ct = std::cos(theta);
		
		V3d y;
		if (std::fabs(rxy) > 1E-12) {	
			y[0] = x[0]*ct - x[0]*x[2]/rxy*ce*st + x[1]/rxy*r*se*st;
			y[1] = x[1]*ct - x[1]*x[2]/rxy*ce*st - x[0]/rxy*r*se*st;
			y[2] = x[2]*ct + rxy*ce*st;
		}
		else {
			y[0] = r*se*st;
			y[1] = r*ce*st;
			y[2] = r*ct;
		}

		return y;
	}

	double HSPotential::theta(const Particle& p1, const Particle& p2, double b, double) const {
		return 2. * std::acos( b / bMax(p1, p2) );
	}

	double HSPotential::bMax(const Particle& p1, const Particle& p2) const {
		return (p1.d + p2.d) / 2.;
	}

	LJPotential::LJPotential(double b_extension, size_t b_size,
			double g_step, double r_step) :
			b_extension(b_extension), b_size(b_size),
			g_step(g_step), g_size(0),
			r_step(r_step) {
		std::cout << "b_extension = " << b_extension
				<< ", b_size" << b_size
				<< ", g_step" << g_step
				<< ", g_size" << g_size
				<< ", r_step" << r_step << std::endl;
	}

	void LJPotential::fillGbToTheta(size_t from, size_t to) const {
		std::cout << "fill form " << from << ' ' << to << std::endl;
		for (size_t i_g = from; i_g < to; ++i_g) {
			double g = (static_cast<double>(i_g) + 0.5) * g_step;
			for (size_t i_b = 0; i_b < b_size; ++i_b) {
				double b = b_extension * (static_cast<double>(i_b) + 0.5) / static_cast<double> (b_size);
				gb2theta[i_b + i_g*b_size] = gbToThetaCalc(g, b);
	//			std::cout << "gb2theta" << ' ' << i_b << ' ' << i_g << ' ' << i_b + i_g*b_size << ' ' 
	//					<< gb2theta[i_b + i_g*b_size] << std::endl;
			}
		}
	}

	void LJPotential::extentGbToTheta(double g) const {
		int new_size = static_cast<int>(g / g_step) + 1;
		gb2theta.resize(new_size * b_size);
		fillGbToTheta(g_size, new_size);
		g_size = new_size;
	}

	inline const V2d f(const V2d r) {
		double one_div_r2 = 1. / sqr(r);
		double one_div_r6 = one_div_r2 * one_div_r2 * one_div_r2;
		double tmp = 24*(2 * one_div_r6 * one_div_r6 - one_div_r6)*one_div_r2;
		return tmp * r;
	}
	double LJPotential::gbToThetaCalc(double g, double b) const {
		double dt = r_step / g;

		V2d x(-std::sqrt(sqr(b_extension) - sqr(b)), b);
		V2d u(g, 0.);

		double r2 = 0.;

		while (r2 < sqr(b_extension)) {
			V2d k1x(u);
			V2d k1u = f(x);

			V2d y	= x + 1./3.*k1x*dt;
			V2d k2x	= u + 1./3.*k1u*dt;
			V2d k2u = f(y);
			
			y		= x - 1./3.*k1x*dt + k2x*dt;
			V2d k3x	= u - 1./3.*k1u*dt + k2u*dt;
			V2d k3u = f(y);

			y		= x + k1x*dt - k2x*dt + k3x*dt;
			V2d k4x	= u + k1u*dt - k2u*dt + k3u*dt;
			V2d k4u = f(y);
			
			x += 1./8.*k1x*dt + 3./8.*k2x*dt + 3./8.*k3x*dt + 1./8.*k4x*dt; 
			u += 1./8.*k1u*dt + 3./8.*k2u*dt + 3./8.*k3u*dt + 1./8.*k4u*dt; 

			r2 = sqr(x);
		}
		return arg(u);
	}

	double LJPotential::bMax(const Particle& p1, const Particle& p2) const {
		return b_extension * (p1.d + p2.d) / 2.;
	}

	double LJPotential::theta(const Particle& q1, const Particle& q2, double b, double g) const {

		const LJParticle& p1 = static_cast<const LJParticle&>(q1);
		const LJParticle& p2 = static_cast<const LJParticle&>(q2);

		double e = std::sqrt(p1.e * p2.e);
		g = g / std::sqrt(e);
		b = 2. * b / (p1.d + p2.d);

		double x_g = g / g_step;
		double x_b = static_cast<double>(b_size) * b / b_extension;
		size_t i_g = static_cast<int>(x_g+0.5);
		size_t i_b = static_cast<int>(x_b+0.5);
		x_g -= static_cast<double>(i_g);
		x_b -= static_cast<double>(i_b);

		double d;
		while (i_g + 1 > g_size) 
			extentGbToTheta(1.2 * g_step * static_cast<double>(i_g + 1));
		if (i_b < b_size)
			d = gb2theta[i_b + i_g*b_size];
		else 
			return 0.;

		double d1, d2;
		if (x_g > 0)
			d1 = gb2theta[i_b + (i_g+1)*b_size] - d;
		else 
			if (i_g > 0) d1 = d - gb2theta[i_b + (i_g-1)*b_size];
			else d1 = 0;
		if (x_b > 0) 
			if (i_b + 1 < b_size) d2 = gb2theta[i_b+1 + i_g*b_size] - d;
			else d2 = 0;
		else 
			if (i_b > 0) d2 = d - gb2theta[i_b-1 + i_g*b_size];
			else d2 = 0;

		return d + d1*x_g + d2*x_b;
	}


}
