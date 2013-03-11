#include <ctime>
#include <algorithm>

#include "timer.h"
#include "manager.h"

#include <mpi.h>

Timer::Timer (int f, int l, int m, int c)
	: finish (f), log (l), macro (m), cache (c), counter (0), 
		printer (manager ().get_printer ())
{
	nick ("init");
}

void Timer::info ()
{
	printer.var ("Total amount of iterations", finish);
	printer.var ("Step of logging", log);
	if (macro > 0) printer.var ("Step of recording macroparameters", macro);
	if (cache > 0) printer.var ("Step of caching distribution function", cache);
}

void Timer::macroparameters ()
{
	nick ("macro");
	printer.result (manager ().get_writer ().write_result (counter));
	manager ().get_writer ().write_f (counter);
}

bool Timer::begin ()
{
	if (counter == 0) {
		nick ("load_f");
		printer.title ("Cached distribution function");
		printer.task ("Loading data");
		printer.result (manager ().get_writer ().load_f (counter));
		printer.var ("Initial time", counter);

		printer.title ("Begin of iterations");
		printer.task ("Saving initial macroparameters");
		macroparameters();
	}
	nick_cpu.clear ();
	nick_real.clear ();
	return (counter++ >= finish);
}

void Timer::end ()
{
	MPI_Barrier (MPI_COMM_WORLD);
	if (counter % macro == 0) {
		printer.task ("Saving macroparameters");
		macroparameters ();
	}
	if (counter % cache == 0) {
		printer.task ("Saving distribution function");
		nick ("save_f");
		printer.result (manager ().get_writer ().save_f (counter));
	}
	nick ("end");
	if (counter % log == 0) {
		printer.title ("LOG");
		printer.var ("Iteration number", counter);
		printer.var ("Transfer spend time", nick_cpu["transfer"]);
		printer.var ("MPI_exchange spend time", nick_real["exchange"]);
		printer.var ("Integral spend time", nick_cpu["integral"]);
		printer.var ("Writing files spend time", nick_real["macro"] + nick_real["save_f"]);
		printer.var ("Total spend cpu time", nick_cpu["all"]);
		printer.var ("Total spend real time", nick_real["all"]);
		printer.title ("/LOG");
	}
}

void Timer::nick (const std::string str)
{
	static std::string phase;
	
	auto prev_real = curr_real;
	auto prev_cpu = curr_cpu;
	
	curr_real = Real_clock::now ();
	curr_cpu = Cpu_clock::now ();
	if (phase == "") { phase = str; return; }
	
	double time;
	const double nano = 1e9;

	time = double (Cpu_chrono::duration_cast<boost::chrono::nanoseconds> (curr_cpu - prev_cpu).count ())/nano;
	nick_cpu[phase] += time;
	nick_cpu["all"] += time;
	
	time = double (Real_chrono::duration_cast<std::chrono::nanoseconds> (curr_real - prev_real).count ())/nano;
	nick_real[phase] += time;
	nick_real["all"] += time;
	
	phase = str;
}
