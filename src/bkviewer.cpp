#include <sstream>
#include <fstream>
#include <stdexcept>
#include "boost/filesystem.hpp"

#include "bkviewer.h"

using namespace Writers;

void Writer_bkviewer::write_static (std::ofstream& file, Int_vect coord, int type)
{
	Vec3<float> v = Vec3<float> (coord) + 0.5f; v *= Box::H;
	file.write (reinterpret_cast<const char*> (&v), sizeof (Vec3<float>));
	v = Vec3<float> (1);
	file.write (reinterpret_cast<const char*> (&v), sizeof (Vec3<float>));
	file.write (reinterpret_cast<const char*> (&type), sizeof (int));
}

void Writer_bkviewer::build_map ()
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

void Writer_bkviewer::prepare_files ()
{
	if (MPI_rank) return;
	assert (size.vol () > 0);
	map = new Matrix<Map> (size);
	build_map ();
	boost::filesystem::create_directory ("static"); 
	boost::filesystem::create_directory ("dynamic");
	boost::filesystem::path filename ("static/net.cells");
	std::ofstream file (filename.string ());
	if (file.fail ()) throw std::runtime_error ("Cannot write files");
	file.write (reinterpret_cast<const char*> (&size), sizeof (Int_vect));
	
	for (int l=0; l<size.x; l++)
		for (int m=0; m<size.y; m++)
			for (int n=0; n<size.z; n++)
				write_static (file, Int_vect (l,m,n), map->all () (Int_vect (l,m,n)).type);
}

Writer_bkviewer::~Writer_bkviewer ()
{
	if (MPI_rank) return;
	if (map) delete map;
}

bool Writer_bkviewer::write_result (int time)
{
	macroparameters ();
	if (MPI_rank) return true;
	std::ostringstream filename;
	filename << std::setw (5) << std::setfill('0') << time;
	boost::filesystem::path path ("dynamic"); path /= filename.str ();
	const int number_of_files = 5;
	std::ofstream files [number_of_files];
	files[0].open (path.string () + "_0.temp");
	files[1].open (path.string () + "_0.dens");
	files[2].open (path.string () + "_0.flow");
	files[3].open (path.string () + "_0.qflow");
	files[4].open (path.string () + "_0.press");
	//for_each_index (map->all (), Write_scalar (files, std::mem_fun (&Features::temp)));
	for (int l=0; l<size.x; l++) 
		for (int m=0; m<size.y; m++)
			for (int n=0; n<size.z; n++) {
				Int_vect p (l,m,n);
				write_param (BIN, files[0], map->all () (p).features->temp);
				write_param (BIN, files[1], map->all () (p).features->dens);
				write_param (BIN, files[2], map->all () (p).features->flow);
				write_param (BIN, files[3], map->all () (p).features->qflow);
				write_param (BIN, files[4], map->all () (p).features->press);
			}
	return true;
}
