#ifndef KOROBOV_H
#define KOROBOV_H

#include "../workers/ci_grider.h"
#include "../base/auxiliary.h"

namespace CI_griders {

	/** STL-like iterator returns integrate points **/
	class Korobov_iterator : public std::iterator<std::forward_iterator_tag, Point> {
		const std::size_t line;
		const Point& shift;
		std::size_t iter;
	public:
		Korobov_iterator (
			std::size_t l,												// line in Korobov table
			const Point& sh,											// random shift
			std::size_t it												// grid counter
		) : line (l), shift (sh), iter (it) {}
		void operator++ () { ++iter; }
		bool operator!= (const Korobov_iterator& other) {
			return iter != other.iter;
		}
		Point&& operator* () const;
	};
	
	/** STL-like container of integrate points **/
	class Korobov_grid {
		std::size_t size_, line;
		Point shift;
	public:
		Korobov_grid (std::size_t size);
		std::size_t size () const { return size_; }
		typedef Korobov_iterator const_iterator;
		const_iterator begin () const {
			return const_iterator (line, shift, 0);
		}
		const_iterator end () const {
			return const_iterator (line, shift, size_);
		}
	};
	
	/** CI_grider using Korobov grid **/
	class Korobov : public CI_grider {
		const Korobov_grid grid;
	public:
		Korobov (std::size_t size) : grid (size) {}
		void generate (real time_step);
		void info ();
	};

} // namespace

#endif // KOROBOV_H