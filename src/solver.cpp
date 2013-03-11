#include "workers/manager.h"
#include "containers/box.h"
#include "ci/ci.hpp"
#include "schemes/tvd_scheme.h"
#include "writers/paraview.h"
#include "ci_griders/korobov.h"

#include <mpi.h>

void couette ()
{
	const real T = 1, n = 1, Ux = 0.01;
	const real mu = 0.562773, k0 = -1.2540, pi = 3.14159;
	const Real_vect U (0, 0, Ux);
	const int Nx = 50;
	const real Kn = 0.03981, corr = 1 - sqrt(pi)*k0*Kn;
	
	manager ().set_grids (4.3*sqrt (T), Kn, Nx, {XX});
	
	Box::init_cond = arbitrary_grad13 (
		[T] (Int_vect) { return T; },
		[n] (Int_vect) { return n; },
		[=] (Int_vect r) { return U * real ((real (r.x)/(Nx-1)-sqrt(pi)/2*k0*Kn)/corr); },
		[] (Int_vect) { return 0; },
		[=] (Int_vect) { return Real_vect (0, -Ux*2*mu*Kn/corr , 0); }	
  	);
	Box* box = new Box (Int_vect (Nx, 1, 1));
	// ------------------------- box --------------------------
	box->set_simple (LEFT, Const_bound (T));
	box->set_simple (RIGHT, Const_bound (T, U));
	manager ().alone_box (box);
}

void heat_transfer ()
{
	const real dT = 0.01, T = 1, Tmax = T+dT/2, Tmin = T-dT/2;
	const real pmax = 1, lambda = 2.129475, d1 = 2.4001, pi = 3.14159;
	const int Nx = 40;
	const real Kn = 0.05, corr = 1 + sqrt(pi)*d1*Kn;
	
	manager ().set_grids (4.3*sqrt (Tmax), Kn, Nx, {XX});
	
 	Box::init_cond = arbitrary_grad13 (
		[=] (Int_vect r) { return (Tmin + dT*r.x/(Nx-1)); },
		[=] (Int_vect r) { return pmax/(Tmin + dT*r.x/(Nx-1)/corr); },
		[] (Int_vect) { return 0; },
		[=] (Int_vect r) { return Real_vect (-lambda*dT*Kn/corr*sqrt (Tmin + dT*r.x/(Nx-1)), 0, 0); }, //-1./sqrt(pi)*dT
		[] (Int_vect) { return 0; }	
  	);
	Box* box = new Box (Int_vect (Nx, 1, 1));
	// ------------------------- box --------------------------
	box->set_simple (LEFT, Const_bound (Tmin));
	box->set_simple (RIGHT, Const_bound (Tmax));
	manager ().alone_box (box);
}

void shock ()
{
	const real M = 3, g = 5./3;
	const real T1 = 1, rho1 = 1, v1 = sqrt(g/2)*M;
	const real T2 = (2*g*sqr(M)-g+1)*((g-1)*sqr(M)+2)/sqr(g+1)/sqr(M), rho2 = (g+1)*sqr(M)/((g-1)*sqr(M)+2), v2 = rho1*v1/rho2;
	const real H = 0.8;
	const int Nx = 8./H;
	
	manager ().set_grids (3.4*sqrt (T2), 1./H, 1, {XX});
	
	Box::init_cond = arbitrary_maxwell (
		[=] (Int_vect r) { return (r.x<Nx) ? T1 : T2; },
		[=] (Int_vect r) { return (r.x<Nx) ? rho1 : rho2; },
		[=] (Int_vect r) { return (r.x<Nx) ? Real_vect (v1,0,0) : Real_vect (v2,0,0); }
  	);
	Box* box = new Box (Int_vect (2*Nx, 1, 1));
	// ------------------------- box --------------------------
	box->set_interface (LEFT, T1, rho1, Real_vect (v1,0,0));
	box->set_interface (RIGHT, T2, rho2, Real_vect (v2,0,0));
	manager ().alone_box (box);
}

void relax ()
{
	const real T = 1, n = 1, u = 1, q = .01, p = 0*.05;
	manager ().set_grids (4.2*sqrt (T), 100, 1, {});
	
	Box::init_cond = arbitrary_grad13 (
		[T] (Int_vect) { return T; },
		[n] (Int_vect) { return n; },
		[u] (Int_vect) { return Real_vect (u, 0, 0); },
		[q] (Int_vect) { return Real_vect (q, 0, 0); },
		[p] (Int_vect) { return Real_vect (p, 0, 0); }	
  	);
	Box* box = new Box (Int_vect (1, 1, 1));
	manager ().alone_box (box);
	manager ().add_tracing_f (box, Int_vect (0));
}

void poiseuille ()
{
	const real T = 1, n1 = 1.2, n2 = 1;
	int Nx, Nz;
	const int multi = 15;
	Nx = 14*multi; Nz = 1*multi;
	real Kn = 0.85;
	
	manager ().set_grids (4.8, Kn, 2*Nz, {XX, ZZ});
	
	Box::init_cond = arbitrary_maxwell ([T] (Int_vect) { return T; }, 
		[Nx, n1, n2] (Int_vect r) { return n1-(n2-n1)*r.x/(Nx-1); });
	Box* box = new Box (Int_vect (Nx, 1, Nz)); // nx,1,nz
	// ------------------------- box --------------------------
	box->set_simple (TOP, Const_bound (T));
	box->set_mirror (BOTTOM);
	box->set_maxwell (LEFT, T, n1);
	box->set_maxwell (RIGHT, T, n2);
	manager ().alone_box (box);
}

void kaskad_knudsen_2d (int num)
{
	const int multi = 2;			// 2
	const int tank = 1;				// 1
	const real Tmin = 1, Tmax = 2;	// 1 2
	int Nx, Nz; 					// 10
	Nx = 32*multi, Nz = 28*multi;	// 32 28
	real Kn = 0.5;					// 0.5
	manager ().set_grids (4.8*sqrt (Tmax), Kn, Nz/2, {XX, ZZ});
	std::vector<Box*> boxes (2*num);
	
	Box::init_cond = new Const_maxwell (Tmin, 1);
	for (int i=0; i<num; i++) {
		boxes[2*i] = new Box (Int_vect (3*Nx/2, 1, Nz/4));		// Nx, 1, Nz/4
		boxes[2*i+1] = new Box (Int_vect (3*Nx/2, 1, Nz/2));	// Nx, 1, Nz/2
		
		boxes[2*i]->set_mirror (BOTTOM);
		boxes[2*i]->set_simple (TOP, Linear_temp_bound (Tmin, Real_vect (Tmax, Tmin, Tmin)));
		boxes[2*i+1]->set_mirror (BOTTOM);
		boxes[2*i+1]->set_simple (TOP, Linear_temp_bound (Tmax, Real_vect (Tmin, Tmax, Tmax)));
		boxes[2*i+1]->set_simple (LEFT, Const_bound (Tmax));
		boxes[2*i+1]->set_simple (RIGHT, Const_bound (Tmin));
	
		if (i>0) manager ().link_boxes (boxes[2*i-1], boxes[2*i], XX, 0);
		manager ().link_boxes (boxes[2*i], boxes[2*i+1], XX, 0);
	}
	Box* box1 = new Box (Int_vect (3*Nx/2*tank, 1, Nz*tank));
	Box* boxN = new Box (Int_vect (3*Nx/2*tank, 1, Nz*tank));
	// ------------------------- box1 --------------------------
	box1->set_mirror (BOTTOM);
	box1->set_simple (TOP, Const_bound (Tmin));
	box1->set_simple (LEFT, Const_bound (Tmin));
	box1->set_simple (RIGHT, Const_bound (Tmin));
	// ------------------------- boxN --------------------------
	boxN->set_mirror (BOTTOM);
	boxN->set_simple (TOP, Const_bound (Tmin));
	boxN->set_simple (LEFT, Const_bound (Tmin));
	boxN->set_simple (RIGHT, Const_bound (Tmin));
	// joining boxes
	// !!! it's obligatory to joining after setting boundary conditions
	manager ().link_boxes (box1, boxes[0], XX, 0);
	manager ().link_boxes (boxes[2*num-1], boxN, XX, 0);
}

void crookes ()
{
	const int multi = 1;
	const int tank = 1;	
	const real Tmin = 1, Tmax = 2;
	const int len = 6*multi;
	real Kn = 1;
	manager ().set_grids (4.8*sqrt (Tmax), Kn, len, {XX, ZZ});
	
	Box::init_cond = new Const_maxwell (Tmin, 1);
	Box* box1 = new Box (Int_vect (len*tank, 1, len*tank));
	Box* box2 = new Box (Int_vect (len*tank, 1, 3*len*tank));
	Box* box3 = new Box (Int_vect (3*len*tank, 1, 3*len*tank));
	Box* box4 = new Box (Int_vect (3*len*tank, 1, len*tank));

	// ------------------------- box1 --------------------------
	box1->set_simple (BOTTOM, Const_bound(Tmin));
	box1->set_simple (LEFT, Const_bound (Tmax));
	// ------------------------- box2 --------------------------
	box2->set_simple (TOP, Const_bound (Tmin));
	// ------------------------- box3 --------------------------
	box3->set_simple (TOP, Const_bound (Tmin));
	box3->set_simple (RIGHT, Const_bound (Tmin));
	// ------------------------- box4 --------------------------
	box4->set_simple (RIGHT, Const_bound (Tmin));
	// joining boxes
	// !!! it's obligatory to joining after setting boundary conditions
	manager ().link_boxes (box1, box2, ZZ, 0);
	manager ().link_boxes (box1, box4, XX, 0);
	manager ().link_boxes (box2, box3, XX, 0);
	manager ().link_boxes (box4, box3, ZZ, 0);
	manager ().link_boxes (box2, LEFT, box4, BOTTOM);
}

int main (int argc, char *argv[])
{
	/** customize collision integral module **/
	ci::Potential* potential = new ci::HSPotential ();		// molecular potential
	ci::Symmetry symmetry = ci::NO_SYMM;					// collision integral symmetry
	ci::init (potential, symmetry);

	/** customize timestamps **/
	int finish = 2, log = 1, macro = 1, cache = 500;
	if (argc>=2) finish = atoi (argv[1]);					// set finish
	if (argc>=3) log = atoi (argv[2]);						// set log
	if (argc>=4) macro = atoi (argv[3]);					// set macro
	if (argc>=5) cache = atoi (argv[4]);					// set cache

	/** customize workers **/
	Mapper::set_radius (16);								// set Vel_grid radius
	manager ().set_workers (
		new Writers::ParaView,								// set program for visualization
		new CI_griders::Korobov (5e5),						// set integrate grid
		new Timer (finish, log, macro, cache),				// set timer parameters
		new TVD_scheme<Lmtr::wide_third>					// set difference scheme
	);
	
	/** choose problem **/
	//kaskad_knudsen_2d (1);
	//crookes ();
	//knudsen_3d ();
	//tube ();
	//cube2 ();
	//poiseuille ();
	heat_transfer ();
	//couette ();
	//shock ();
	//relax ();
	//knudsen_slip ();
	
	/** start simulation **/
	manager ().init_model ();									// prepare calculations
	manager ().iterate ();										// cycle of iterations
	return EXIT_SUCCESS;
}

