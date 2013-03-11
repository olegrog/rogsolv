#ifndef PARAVIEW_H
#define PARAVIEW_H

#include "../workers/writer.h"

namespace Writers {
	class ParaView : public Writer {
		int points, cells;
		void prepare_files ();
	public:
		bool write_result (int) const;
		void info () const;
	};
}

#endif
