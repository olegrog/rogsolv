#ifndef TIMER_H
#define TIMER_H

#include <map>
#include <string>

#include "box.h"
#include "printer.h"
#include "writer.h"

class Timer {
	int finish, log, macro, cache;
	const Boxes& boxes;
	int counter;
	std::map<std::string, double> nick_cpu, nick_real;
	timespec cpu, real;
	Printer* printer;
public:
	Writer* writer;
	Timer (Printer*, Writer*, const Boxes&);
	bool begin ();
	void end ();
	void nick (const std::string);
	void set_parameters (int finish, int log, int macroparameters, int cache_f);
};

#endif

