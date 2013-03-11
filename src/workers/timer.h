#ifndef TIMER_H
#define TIMER_H

#include <map>
#include <string>
#include <chrono>
#include <boost/chrono.hpp>

#include "../containers/box.h"
#include "printer.h"
#include "writer.h"

namespace Real_chrono = std::chrono;
namespace Cpu_chrono = boost::chrono;

class Timer {
	int finish, log, macro, cache;
	int counter;
	const Printer& printer;
	std::map<std::string, double> nick_cpu, nick_real;

	typedef Real_chrono::system_clock Real_clock;
	typedef Cpu_chrono::process_real_cpu_clock Cpu_clock;
	Real_chrono::time_point<Real_clock> curr_real;
	Cpu_chrono::time_point<Cpu_clock> curr_cpu;
	
	void macroparameters ();
public:
	Timer (
		int finish,												// stop time
		int log,												// interval of logging
		int macro,												// interval of saving macroparameters
		int cache												// interval of caching time sample
	);
	bool begin ();												// timestamp for starting iteration
	void end ();												// timestamp for finishing iteration
	void nick (const std::string);								// special timestamps
	void info ();												// print special information to log file
};

#endif

