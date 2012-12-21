#ifndef LIMITERS_H
#define LIMITERS_H

#include <cmath> // for abs, min, max
#include <algorithm> // for min_element
#include <string>

#include "auxiliary.h"

namespace Lmtr {

	/**
	 * "Name" is the name of appropriate flux limiter 
	 * class Limiter<Name> is a functor returns closure of limiter calculations for Vel_grid operations (using Expr<> shell)
	 * First, Second, Third are f(i-1), f(i), f(i+1), CFL is a Courant number
	 * 		#  - quoting (stringizing) operator (#Name = "Name")
	 *		## - concatenation operator (Name_##impl = Name_impl)
	 * class Name is empty class-switch without templates (like metafunction class)
	 * class Name_impl<...> is closure of limiter calculations (defined operator[])
	 */
	
	template<class> class Limiter { };
	
	#define DEFINE_LIMITER(Name)															\
	class Name { }; 																		\
	namespace {																				\
		template<class First, class Second, class Third, class CFL>							\
		class Name##_impl {																	\
			const First& grid1;																\
			const Second& grid2;															\
			const Third& grid3;																\
			const CFL& cfl;																	\
		public:																				\
			Name##_impl (const First& g1, const Second& g2, const Third& g3, const CFL& c)	\
				: grid1 (g1), grid2 (g2), grid3 (g3), cfl (c) { }							\
			real operator[] (int) const;													\
		};																					\
	}																						\
	template<> struct Limiter<Name> {														\
		template<class First, class Second, class Third, class CFL>							\
		Expr<Name##_impl<First, Second, Third, CFL> > operator ()							\
			(const First& g1, const Second& g2, const Third& g3, const CFL& cfl)			\
		{																					\
			typedef Name##_impl<First, Second, Third, CFL> Clos;							\
			return Expr<Clos> (Clos (g1, g2, g3, cfl));										\
		}																					\
		std::string name () const { return std::string (#Name); }							\
	};
	DEFINE_LIMITER (none)
	DEFINE_LIMITER (mc)
	DEFINE_LIMITER (Koren)
	DEFINE_LIMITER (minmod)
	DEFINE_LIMITER (superbee)
	DEFINE_LIMITER (van_Leer)
	DEFINE_LIMITER (van_Albada)
	DEFINE_LIMITER (wide_superbee)
	DEFINE_LIMITER (wide_third)
	#undef DEFINE_LIMITER
	
	namespace {
		struct Abs_comp {
			bool operator () (real x, real y)
			{
				return std::abs (x) < std::abs (y);
			}
		};

		template<class First, class Second, class Third, class CFL>
		inline real none_impl<First, Second, Third, CFL>::operator[] (int) const { return 0; }
		
		template<class First, class Second, class Third, class CFL>
		inline real mc_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[3] = {grid2[i] - grid1[i], grid3[i] - grid2[i], (grid3[i] - grid1[i])*.25};
			if (d[0]*d[1] <= 0) return 0;
			return (1.-cfl[i])**std::min_element (d, d+3, Abs_comp ());
		}
		
		template<class First, class Second, class Third, class CFL>
		inline real Koren_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[3] = {grid2[i] - grid1[i], grid3[i] - grid2[i], (2*grid3[i] - grid2[i] - grid1[i])*(1./6)};
			if (d[0]*d[1] <= 0) return 0;
			return (1.-cfl[i])**std::min_element (d, d+3, Abs_comp ());
		}
		
		template<class First, class Second, class Third, class CFL>
		inline real minmod_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[2] = {grid2[i] - grid1[i], grid3[i] - grid2[i]};
			if (d[0]*d[1] <= 0) return 0;
			return 0.5*(1.-cfl[i])**std::min_element (d, d+2, Abs_comp ());
		}
		
		template<class First, class Second, class Third, class CFL>
		inline real superbee_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[2] = {grid2[i] - grid1[i], grid3[i] - grid2[i]};
			if (d[0]*d[1] <= 0) return 0;
			const real sgn = (d[0] < 0) ? -0.5 : 0.5;
			const real d1[2] = {2*std::abs (grid2[i] - grid1[i]), std::abs (grid3[i] - grid2[i])};
			const real d2[2] = {std::abs (grid2[i] - grid1[i]), 2*std::abs (grid3[i] - grid2[i])};
			return sgn*(1.-cfl[i])*std::max (std::min (d1[0], d1[1]), std::min (d2[0], d2[1]));
		}

		template<class First, class Second, class Third, class CFL>
		inline real wide_superbee_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[2] = {grid2[i] - grid1[i], grid3[i] - grid2[i]};
			if (d[0]*d[1] <= 0) return 0;
			const real sgn = (d[0] < 0) ? -0.5 : 0.5;
			const real d1[2] = {2./cfl[i]*std::abs (grid2[i] - grid1[i]), std::abs (grid3[i] - grid2[i])};
			const real d2[2] = {std::abs (grid2[i] - grid1[i]), 2./(1-cfl[i])*std::abs (grid3[i] - grid2[i])};
			return sgn*(1.-cfl[i])*std::max (std::min (d1[0], d1[1]), std::min (d2[0], d2[1]));
		}

		template<class First, class Second, class Third, class CFL>
		inline real van_Leer_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[2] = {grid2[i] - grid1[i], grid3[i] - grid2[i]};
			if (d[0]*d[1] <= 0) return 0;
			return (1.-cfl[i])*d[0]*d[1]/(d[0]+d[1]);
		}
		
		template<class First, class Second, class Third, class CFL>
		inline real van_Albada_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[2] = {grid2[i] - grid1[i], grid3[i] - grid2[i]};
			if (d[0]*d[1] <= 0) return 0;
			return (1.-cfl[i])*d[0]*d[0]*d[1]/(d[0]*d[0]+d[1]*d[1]);
		}
		
		template<class First, class Second, class Third, class CFL>
		inline real wide_third_impl<First, Second, Third, CFL>::operator[] (int i) const
		{
			const real d[3] = {	(1.-cfl[i]) / cfl[i] * (grid2[i] - grid1[i]), (grid3[i] - grid2[i]),
								(1.-cfl[i])/6*(-(1.+cfl[i])*grid1[i] + (2.*cfl[i]-1.)*grid2[i] + (2.-cfl[i])*grid3[i]) };
			if (d[0]*d[1] <= 0) return 0;
			return *std::min_element (d, d+3, Abs_comp ());
		}
		
	} // namespace
} // namespace Lmtr


#endif // LIMITERS_H