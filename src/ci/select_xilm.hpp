#include <boost/math/special_functions/next.hpp>			// for float_distance
	
namespace ci {

	using namespace dod_vector;

	/** Functor for selecting nodes from stencil **/
	class Select_xilm {
	public:
		virtual bool operator () (const std::vector<V3i>&, double&, V3i&, V3i&) const = 0;
		virtual ~Select_xilm () {}
	};
	
	/** minimize functional:
	 *	J = (1-r)*sqr(xi_l-xi) + r*sqr(xi_m-xi) 
	 **/
	class Min_delta_p : public Select_xilm {
		const int nk_rad2;
		const double m2, E0;
		const V3d u, w;
	public:
		Min_delta_p (int r_, double m2_, double E0_, V3d u_, V3d w_) : nk_rad2 (r_), m2 (m2_), E0 (E0_), u (u_), w (w_) { }
		bool operator () (const std::vector<V3i>& stencil, double& r, V3i& xi2l, V3i& xi2m) const {
			double q = std::numeric_limits<double>::max();
			bool is_found = false;
			for (auto pl : stencil) { 
				V3d xil = i2xi(pl, nk_rad2) - m2*u;
				double El = sqr(xil);
				if (boost::math::float_distance(E0, El) > 1) // El > E0 >= Em
					for (auto pm : stencil) {
						V3d xim = i2xi(pm, nk_rad2) - m2*u;
						double Em = sqr(xim);
						if (Em <= E0) {
							double r_ = (E0-El)/(Em-El);
							double ql = sqr(xil - w);
							double qm = sqr(xim - w);
							double q_ = (1-r_)*ql + r_*qm;
							if (q_ < q) {
								r = r_; q = q_;
								xi2l = pl; xi2m = pm;
								is_found = true;
							}
						}
					}
			}
			return is_found;
		}
	};

	/** minimize functionals:
	 *	J_{l,m} = |sqr(xi_{l,m}) - sqr(xi)|
	 *	sqr(xi_m) < sqr(xi) < sqr(xi_l)
	 **/
	class Min_delta_E : public Select_xilm {
		const int nk_rad2;
		const double m2, E0;
		const V3d u, w;
	public:
		Min_delta_E (int r_, double m2_, double E0_, V3d u_, V3d w_) : nk_rad2 (r_), m2 (m2_), E0 (E0_), u (u_), w (w_) { }
		bool operator () (const std::vector<V3i>& stencil, double& r, V3i& xi2l, V3i& xi2m) const {
			double Em = 0;
			double El = std::numeric_limits<double>::max();
			bool found_l = false;
			bool found_m = false;
			for (auto p : stencil) { 
				V3d xi = i2xi(p, nk_rad2) - m2*u;
				double E = sqr(xi);
				if (boost::math::float_distance(E0, E) > 1) { // E > E0
					if (E < El) {
						El = E;
						xi2l = p;
						found_l = true;
					}
				} else {
					if (E > Em) {
						Em = E;
						xi2m = p;
						found_m = true;
					}
				}
			}
			r = (E0-El)/(Em-El);
			return found_l && found_m;
		}
	};
	
} // namespace
