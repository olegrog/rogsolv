#ifndef PRINTER_H
#define PRINTER_H

#include <set>

#include "../containers/box.h"

class Printer {
	int MPI_rank;
public:
	Printer ();
	template<class T>
		void var (const std::string&, T) const;							// print special variable
	void title (const std::string&) const;									// chapter name
	void log (const std::string&) const;									// simple log printing
	void task (const std::string&) const;									// begin of task
	void result (bool) const;												// end of task
	void boxes (const Boxes&) const;										// print table of boxes
	void MPI_ranks (const Boxes&) const;									// print map of boxes & MPI_ranks
};

template<class T>
void Printer::var (const std::string& str, T v) const
{
	if (!MPI_rank) std::cout << str << " = " << v << std::endl;
}

#endif // PRINTER_H
