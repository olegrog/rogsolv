#include <sstream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include "boost/filesystem.hpp"
#include <sys/stat.h>

#include "writer.h"
#include "manager.h"
#include "../containers/vel_grid.h"

#include "mpi.h"

using namespace Writers;

const int precision = 6;	// number of figures in text output (float has 6 precise figures, double - 15)
void Writers::write_param (Files_format fmt, std::ofstream& file, real value)
{
	if (fmt == BIN) {
		float v = static_cast<float> (value);
		file.write (reinterpret_cast<const char*> (&v), sizeof (float));
	} 
	if (fmt == TXT)
		file << std::scientific << std::setprecision (precision) << value << std::endl;
}

void Writers::write_param (Files_format fmt, std::ofstream& file, const Real_vect& value)
{
	if (fmt == BIN) {
		Vec3<float> v = value;
		file.write (reinterpret_cast<const char*> (&v), sizeof (Vec3<float>));
	}
	if (fmt == TXT) 
		file << std::scientific << std::setprecision (precision)
			<< value.x << ' ' << value.y << ' ' << value.z << std::endl;
}

void Writers::read_param (Files_format fmt, std::ifstream& file, real& value)
{
	if (fmt == BIN) {
		float v = static_cast<float> (value);
		file.read (reinterpret_cast<char*> (&v), sizeof (float));
	} 
	if (fmt == TXT)
		file >> value;
}

void Writers::read_param (Files_format fmt, std::ifstream& file, Real_vect& value)
{
	if (fmt == BIN) {
		Vec3<float> v = value;
		file.read (reinterpret_cast<char*> (&v), sizeof (Vec3<float>));
	}
	if (fmt == TXT) 
		file >> value.x >> value.y >> value.z;
}

void Writer::build_map ()
{
	static Features solid (-1, -1, 0, 0, 0, 0);
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Box* box = *pbox;
		for_each (map->volume (box->coord (), box->size ()), box->features->all (), [] (Map& map, Features& feat) {
			map.type = GAS;
			map.features = &feat;										// fill gas cells
		});
	}
	for_each (map->all (), [&] (Map& map) {
		if (map.type == SOLID)
			map.features = &solid;										// fill solid cells
	});
}

// 1) calculate macroparameters from distribution function
// 2) send data to process with rank 0
void Writer::macroparameters () const
{
	MPI_Datatype datatype;
	MPI_Status status;
	if (sizeof (real) == sizeof (double)) datatype = MPI_DOUBLE;
	if (sizeof (real) == sizeof (float)) datatype = MPI_FLOAT;
	
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Box* box = *pbox;
		if (MPI_rank == box->MPI_rank ()) Macroparameters () (box);
		
		int package_size = box->size ().vol () * static_cast<int> (sizeof (Features) / sizeof (real));
		if (MPI_rank && MPI_rank == box->MPI_rank ())					// if (I'm not writer and it's my box)
			MPI_Send (*box->features, package_size, datatype, 0, 0, MPI_COMM_WORLD);
		if (!MPI_rank && MPI_rank != box->MPI_rank ()) 				// if (I'm writer and it's not my box)
			MPI_Recv (*box->features, package_size, datatype, box->MPI_rank (), 0, MPI_COMM_WORLD, &status);
	}
}

void Writer::add_point (Box* box, const Int_vect& pp)
{
	points.insert(std::make_pair (box, pp));
}

void Writer::write_f (int time) const
{
	/** write down macroparameters */
	if (!MPI_rank) {
		const Files_format format = TXT;
		std::string dir = "macro";
		boost::filesystem::create_directory (dir);
		boost::filesystem::path path (dir);
		std::ofstream f_temp ((path/"temp").string ());
		std::ofstream f_dens ((path/"dens").string ());
		std::ofstream f_flow ((path/"flow").string ());
		std::ofstream f_qflow ((path/"qflow").string ());
		std::ofstream f_shear ((path/"shear").string ());
		
		for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, f_temp, m.features->temp); });
		for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, f_dens, m.features->dens); });
		for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, f_flow, m.features->flow); });
		for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, f_qflow, m.features->qflow); });
		for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, f_shear, m.features->shear); });
	}
	
	/** write down distribution function */
	if (points.empty ()) return;
	std::string dir = "dist_fun";
	boost::filesystem::create_directory (dir);
	mapper ().set_side ();
	int num = 0;
	for (Proj_points::const_iterator it = points.begin (); it != points.end (); ++it) {
		Box* box = it->first;
		if (MPI_rank != box->MPI_rank ()) continue;

		std::ostringstream filename;
		filename << time << "_" << num << "_f.txt";
		boost::filesystem::path path (dir); path /= filename.str ();
		std::ofstream file (path.string ());
		if (file.fail ()) throw std::runtime_error ("Cannot write files with distribution function");

		char delim = ' ';
		file << mapper ().radius () << delim <<  Vel_grid::cut_vel () << std::endl;
		const Int_vect& pp = it->second;
		auto& node = box->features->all () (pp);
		file << node.dens << delim << node.temp << std::endl;
		for (int ax=0; ax<3; ax++) file << node.flow[ax] << delim;
		file << std::endl << std::endl;
		for_each (box->f->all () (pp), [&file, delim] (Int_vect r, real dist_fun) {
			for (int ax=0; ax<3; ax++) file << Vel_grid::vel_array () [r[ax]] << delim;
			file << dist_fun << std::endl;
		});
		num++;
	}
}

bool Writer::save_f (int time) const
{
	std::ostringstream filename;
	filename << "f.cache_" << std::setw (2) << std::setfill('0') << MPI_rank;
	std::ofstream file (filename.str ());
	if (file.fail ()) return false;
	int num;
	MPI_Comm_size (MPI_COMM_WORLD, &num);
	file.write (reinterpret_cast<char*> (&time), sizeof (int));
	file.write (reinterpret_cast<char*> (&num), sizeof (int));
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Box* box = *pbox;
		if (MPI_rank != box->MPI_rank ()) continue;
		for_each (box->f->all (), [&file] (const Vel_grid& dist) {
			file.write (dist.raw_data (), mapper ().volume () * sizeof (real));
		});
	}
	return true;
}

bool Writer::load_f (int& time) const
{
	std::ostringstream filename;
	filename << "f.cache_" << std::setw (2) << std::setfill('0') << MPI_rank;
	std::ifstream file (filename.str ());
	if (file.fail ()) return false;
	int num, real_num;
	MPI_Comm_size (MPI_COMM_WORLD, &num);
	file.read (reinterpret_cast<char*> (&time), sizeof (int));
	file.read (reinterpret_cast<char*> (&real_num), sizeof (int));
	if (num != real_num) throw std::runtime_error ("Another MPI_Comm_size");
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Box* box = *pbox;
		if (MPI_rank != box->MPI_rank ()) continue;
		for_each (box->f->all (), [&file] (Vel_grid& dist) {
			file.read (dist.raw_data (), mapper ().volume () * sizeof (real));
		});
	}
	if (file.eof ()) throw std::runtime_error ("Small cache files");
	file.read (reinterpret_cast<char*> (&num), 1);				// just read 1 byte
	if (!file.eof ()) throw std::runtime_error ("Big cache files");
	return true;
}

bool Writer::load_macro () const
{
	const Files_format format = TXT;
	std::string dir = "macro";
	boost::filesystem::path path (dir);
	if (!boost::filesystem::exists(path)) return false;

	MPI_Datatype datatype;
	MPI_Status status;
	if (sizeof (real) == sizeof (double)) datatype = MPI_DOUBLE;
	if (sizeof (real) == sizeof (float)) datatype = MPI_FLOAT;
	
	// load macroparameters from files
	if (!MPI_rank) {
		std::ifstream f_temp ((path/"temp").string ());
		std::ifstream f_dens ((path/"dens").string ());
		std::ifstream f_flow ((path/"flow").string ());
		std::ifstream f_qflow ((path/"qflow").string ());
		std::ifstream f_shear ((path/"shear").string ());
		
		for_each (map->all (), [&] (Map& m) { if (m.type == GAS) read_param (format, f_temp, m.features->temp); });
		for_each (map->all (), [&] (Map& m) { if (m.type == GAS) read_param (format, f_dens, m.features->dens); });
		for_each (map->all (), [&] (Map& m) { if (m.type == GAS) read_param (format, f_flow, m.features->flow); });
		for_each (map->all (), [&] (Map& m) { if (m.type == GAS) read_param (format, f_qflow, m.features->qflow); });
		for_each (map->all (), [&] (Map& m) { if (m.type == GAS) read_param (format, f_shear, m.features->shear); });
		
		if (f_temp.eof ()) throw std::runtime_error ("Small macro files");
		int dummy;
		f_temp.read (reinterpret_cast<char*> (&dummy), 1);				// just read 1 byte
		if (!f_temp.eof ()) throw std::runtime_error ("Big macro files");

	}
	
	// send them to the appropriate process
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Box* box = *pbox;
		int package_size = box->size ().vol () * static_cast<int> (sizeof (Features) / sizeof (real));
		if (MPI_rank && MPI_rank == box->MPI_rank ())					// if (I'm not writer and it's my box)
			MPI_Recv (*box->features, package_size, datatype, 0, 0, MPI_COMM_WORLD, &status);
		if (!MPI_rank && MPI_rank != box->MPI_rank ()) 				// if (I'm writer and it's not my box)
			MPI_Send (*box->features, package_size, datatype, box->MPI_rank (), 0, MPI_COMM_WORLD);
	}
	
	// update Vel_grid
	mapper().set_side ();
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Box* box = *pbox;
		if (MPI_rank != box->MPI_rank ()) continue;
		for_each (box->f->all (), box->features->all (), [&] (Vel_grid& grid, const Features& feat) {
			grid = std::move (Vel_grid (feat.temp, feat.dens, feat.flow/feat.dens));
		});
	}
	
	return true;
}

Writer::Writer () :
	boxes (manager ().get_boxes ()), MPI_rank (printer ().MPI_rank ()), map (nullptr) {}

Writer::~Writer()
{
	if (map) delete map;
}
