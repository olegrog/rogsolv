#ifndef TIMER_H
#define TIMER_H

#include <map>
#include <string>
#include <chrono>
#include <boost/chrono.hpp>

#include "../containers/box.h"
#include "printer.h"
#include "../writers/writer.h"

namespace Real_chrono = std::chrono;
namespace Cpu_chrono = boost::chrono;

class Timer {
	int finish, log, macro, cache;
	const Boxes& boxes;
	int counter;
	std::map<std::string, double> nick_cpu, nick_real;

	typedef Real_chrono::steady_clock Real_clock;
	typedef Cpu_chrono::process_real_cpu_clock Cpu_clock;
	Real_chrono::time_point<Real_clock> curr_real;
	Cpu_chrono::time_point<Cpu_clock> curr_cpu;
	
	Printer* printer;
	void macroparameters ();
public:
	Writer* writer;
	Timer (Printer*, Writer*, const Boxes&);
	bool begin ();
	void end ();
	void nick (const std::string);
	void set_parameters (int finish, int log, int macroparameters, int cache_f);
};

#endif

