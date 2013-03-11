#ifndef KOROBOV_H
#define KOROBOV_H

#include <array>

#include "../workers/ci_grider.h"
#include "../base/auxiliary.h"

namespace CI_griders {

	/** STL-like iterator returns Korobov points **/
	class Korobov_iterator : public std::iterator<std::forward_iterator_tag, Point> {
	public: /* FIXME */
		const std::size_t line;
		const Point& shift;
		std::size_t iter;
	public:
		Korobov_iterator (
			std::size_t l,												// line in Korobov table
			const Point& sh,											// random shift
			std::size_t it												// Korobov grid counter
		) : line (l), shift (sh), iter (it) {}
		void operator++ () { ++iter; }
		bool operator!= (const Korobov_iterator& other) {
			return iter != other.iter;
		}
		Point&& operator* () const;
	};
	
	/** STL-like container of Korobov points **/
	class Korobov_grid {
		std::size_t size_, line;
		Point shift;
	public:
		Korobov_grid (std::size_t size);
		std::size_t size () const { return size_; }
		typedef Korobov_iterator const_iterator;
		const_iterator begin () const;
		const_iterator end () const;
	};
	
	/** CI_grider using Korobov grid **/
	class Korobov : public CI_grider {
	public: /* FIXME */
		const Korobov_grid grid;
	public:
		Korobov (std::size_t size) : grid (size) {}
		void generate (real time_step);
		void info ();
	};

} // namespace CI_grids

#endif // KOROBOV_H