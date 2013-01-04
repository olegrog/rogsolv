#include <boost/math/constants/constants.hpp>
#include <stdexcept>

#include "vel_grid.h"

real Vel_grid::cut, Vel_grid::dV;
std::valarray<real> Vel_grid::vel;

inline real sqr (real x) { return x*x; }

const Vel_grid& Vel_grid::operator= (const Vel_grid& expr)
{
	int size = m.size ();
	for (int i=0; i<size; i++) {
		int j = m.index (i);
		(*v)[j] = expr [j];
	}
	return *this;
}

const Vel_grid& Vel_grid::operator= (const Half_grid& expr)
{
	int size = m.size ();
	for (int i=0; i<size; i++)
		(*v)[m.index (i)] = (*expr.v)[i];
	return *this;
}

const Half_grid& Half_grid::operator= (const Half_grid& expr)
{
	int size = m.size ();
	for (int i=0; i<size; i++)
		(*v)[i] = (*expr.v)[i];
	return *this;
}

const Half_grid& Half_grid::operator= (const Vel_grid& expr)
{
	int size = m.size ();
	for (int i=0; i<size; i++)
		(*v)[i] = expr [m.index (i)];
	return *this;
}

void Vel_grid::set_cut (real vel_cut)
{
	int R = mapper ().radius ();			// radius of velocity grid
	cut = vel_cut;							// cutting velocity
	vel.resize (2*R);						// velocity vector
	for (int i=0; i<2*R; i++)
		vel[i] = (.5+i-R)*cut/R;	
	dV = std::pow (cut/R, 3);				// velocity cell volume
}

// ----------- <for debug> --------------
void print_vel_grid (const Vel_grid& grid)
{
	const Mapper& m = mapper ();
	Int_vect p;
	int R = m.radius();
	std::cout<<'\n';
	std::cout.precision (0);
	for (p.z=0; p.z<2*R; p.z++)	 {
		for (p.y=0; p.y<2*R; p.y++) {
			for (p.x=0; p.x<2*R; p.x++)
				if (sqr (.5+p.x-R) + sqr (.5+p.y-R) + sqr (.5+p.z-R) <= sqr (R))
					std::cout<< std::scientific << grid[m (p)] <<'\t';
				else  std::cout<<"  *\t";
			std::cout<<'\n';
		}
		std::cout << '\n';
	}
	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}

void print_half_grid (const Half_grid& grid, Side side)
{
	const Mapper& m = mapper ();
	Int_vect p;
	int R = m.radius();
	std::cout<<'\n';
	std::cout.precision (0);
	Axis ax = axis (side);
	Direction dir = direction (side);
	for (p.z=0; p.z<2*R; p.z++)	 {
		for (p.y=0; p.y<2*R; p.y++) {
			for (p.x=0; p.x<2*R; p.x++)
				if (sqr (.5+p.x-R) + sqr (.5+p.y-R) + sqr (.5+p.z-R) <= sqr (R)) {
					if ((dir==BACKWARD && p[ax]<R) || (dir==FORWARD && p[ax]>=R))
						std::cout<< std::scientific << grid[m (p)] <<'\t';
					else
						std::cout << "  ###\t";
				}
				else  std::cout<<"   *\t";
			std::cout<<'\n';
		}
		std::cout << '\n';
	}
	std::cout<<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
}
// ----------- </for debug> --------------

Vel_grid::Vel_grid (real temp, real dens, Real_vect speed, Real_vect qflow, Real_vect shear) : 
	m (mapper ()), v (new std::valarray<real> (m.volume ()))
{
	// create maxwell distribution function with given temperature, density and speed
	const real coeff = dens * std::pow (boost::math::constants::pi<real> ()*temp, -1.5);
	const real k = -1./temp;
	for (int i=0; i<m.volume (); i++) {
		const Int_vect& p = m[i];
		(*v)[i] = coeff * std::exp (k*(sqr(Real_vect (vel[p.x], vel[p.y], vel[p.z]) - speed)));
	}
	if (sqr (qflow) == 0 && sqr (shear) == 0) return;				// just optimization
	// build the Grad 13-moment approximation from maxwell distribution
	const real k1 = 1./(dens*sqr(temp));
	const real k2 = 4.*k1/(5.*temp);
	for (int i=0; i<m.volume (); i++) {
		const Int_vect& p = m[i];
		Real_vect c = Real_vect (vel[p.x], vel[p.y], vel[p.z]) - speed;
		real factor = 1 + k1*shear_matrix_product (shear, c) + k2*dot (qflow, c)*(sqr (c)-2.5*temp);
		if (factor <= 0) throw std::logic_error ("Initial Grad13 distribution has negative values.");
		(*v)[i] *= 1 + k1*shear_matrix_product (shear, c) + k2*dot (qflow, c)*(sqr (c)-2.5*temp);
	}
}

Vel_grid::~Vel_grid () { delete v; }

real Vel_grid::dens () const
{
	return static_cast<real> (sum (*this)) * dV;
}

Real_vect Vel_grid::flow () const
{
	Real_vect result;
	for (int s=0; s<3; s++)
		result[s] = sum (Vector (Axis (s), vel) * (*this));
	return result*dV;
}

Real_vect Vel_grid::qflow (Real_vect speed) const
{
	Real_vect result;
	const std::valarray<real> c[3] = {vel - speed.x, vel - speed.y, vel - speed.z};
	const std::valarray<real> sqr_c[3] = {c[0]*c[0], c[1]*c[1], c[2]*c[2]};
	Vel_grid sqr_c_grid; sqr_c_grid = (*this) *
		(Vector (XX, sqr_c[XX]) + Vector (YY, sqr_c[YY]) + Vector (ZZ, sqr_c[ZZ]));
	for (int s=0; s<3; s++)
		result[s] = sum (Vector (Axis (s), c[s]) * sqr_c_grid);
	return result*dV;
}

Real_vect Vel_grid::press (Real_vect speed) const
{
	Real_vect result;
	const std::valarray<real> sqr_c[3] = {(vel - speed.x)*(vel - speed.x), 
		(vel - speed.y)*(vel - speed.y), (vel - speed.z)*(vel - speed.z)};
	for (int s=0; s<3; s++)
		result[s] = sum (Scalar (2) * Vector (Axis (s), sqr_c[s]) * (*this));
	return result*dV;
}

Real_vect Vel_grid::shear (Real_vect speed) const
{
	Real_vect result;
	const std::valarray<real> c[3] = {vel - speed.x, vel - speed.y, vel - speed.z};
	for (int s=0; s<3; s++)
		result[s] = sum (Scalar (2) * Vector (Axis ((s+1)%3), c[(s+1)%3]) * Vector (Axis ((s+2)%3), c[(s+2)%3]) * (*this));
	return result*dV;
}

