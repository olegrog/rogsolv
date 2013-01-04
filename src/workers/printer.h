#ifndef PRINTER_H
#define PRINTER_H

#include <set>

#include "../containers/box.h"

class Printer {
	int MPI_rank;
public:
	Printer (int rank) : MPI_rank (rank) { }
	template<class T>
		void var (const std::string&, T);
	void title (const std::string&);
	void log (const std::string&);
	void task (const std::string&);
	void result (bool);
	void boxes (const Boxes&);
	void MPI_ranks (const Boxes&);
};

template<class T>
void Printer::var (const std::string& str, T v)
{
	if (!MPI_rank) std::cout << str << " = " << v << std::endl;
}

#endif
