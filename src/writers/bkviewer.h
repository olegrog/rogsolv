#ifndef BKVIEWER_H
#define BKVIEWER_H

#include "../workers/writer.h"

namespace Writers {
	class BKViewer : public Writer {
		Matrix<Map>* map;
		void write_static (std::ofstream&, Int_vect, int) const;
		void build_map ();
		void prepare_files ();
	public:
		bool write_result (int) const;
		void info () const;
		~BKViewer ();
	};
}


#endif
