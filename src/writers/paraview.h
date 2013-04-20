#ifndef PARAVIEW_H
#define PARAVIEW_H

#include <array>

#include "../workers/writer.h"

namespace Writers {
	class ParaView : public Writer {
		typedef std::array<int,8> Cell;
		std::vector<Real_vect> points;
		std::vector<Cell> cells;
		void prepare_files ();
	public:
		bool write_result (int) const;
		void info () const;
	};
}

#endif
