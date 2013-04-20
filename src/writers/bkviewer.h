#ifndef BKVIEWER_H
#define BKVIEWER_H

#include "../workers/writer.h"

namespace Writers {
	class BKViewer : public Writer {
		void write_static (std::ofstream&, Int_vect, int) const;
		void prepare_files ();
	public:
		bool write_result (int) const;
		void info () const;
		~BKViewer ();
	};
}


#endif
