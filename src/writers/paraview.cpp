#include <sstream>
#include <fstream>
#include <stdexcept>
#include "boost/filesystem.hpp"

#include "paraview.h"
#include "../workers/printer.h"

using namespace Writers;

void ParaView::prepare_files ()
{
	if (printer ().MPI_rank ()) return;
	assert (size.vol () > 0);
	map = new Matrix<Map> (size);
	build_map ();

	// fill vectot points
	std::vector<Int_vect> stencil;
	for (int s1 = 0; s1 < 2; ++s1)
		for (int s2 = 0; s2 < 2; ++s2)
			for (int s3 = 0; s3 < 2; ++s3)
				stencil.push_back (Int_vect (s1,s2,s3));
	for (int n=0; n<=map->size ().z; n++)
		for (int m=0; m<=map->size ().y; m++)
			for (int l=0; l<=map->size ().x; l++) {
				bool found = false;
				Int_vect c (l,m,n);
				std::for_each (stencil.begin (), stencil.end (), [&] (const Int_vect& s) {
					Int_vect p = c - s;
					if (p < 0 || (!(p - map->size ()).vol())) return;
					if (map->point(p).type == GAS) {
						map->point(p).vertices.insert (points.size ());
						found = true;
					}
				});
				if (found)
					points.push_back (c);
			}

	// fill vector cells
	for_each (map->all (), [&] (const Map& map) {
		if (map.type == GAS) {
			Cell cell;
			std::copy (map.vertices.begin (), map.vertices.end (), cell.begin ());
			cells.push_back (cell);
		}
	});

	boost::filesystem::create_directory ("VTK");
}

bool ParaView::write_result (int time) const
{
	const Files_format format = TXT;
	macroparameters ();
	if (printer ().MPI_rank ()) return true;
	std::ostringstream filename;
	filename << time << "data.vtk";
	boost::filesystem::path path ("VTK"); path /= filename.str ();
	std::ofstream file (path.string ());
	if (file.fail ()) throw std::runtime_error ("Cannot write VTK files");
	std::string fmt;
	if (format == BIN) fmt = "BINARY";
	if (format == TXT) fmt = "ASCII";
	file << "# vtk DataFile Version 3.0\nrogsolv output for iteration number " << time 
		<< "\n" << fmt << "\nDATASET UNSTRUCTURED_GRID\n";
	
	file << "POINTS " << points.size () << " double\n";
	std::for_each (points.begin (), points.end (), [&] (const Int_vect& p) {
		write_param (format, file, Real_vect (p) * Box::H);
	});

	file << "CELLS " << cells.size () << ' ' << 9*cells.size () << '\n';
	std::for_each (cells.begin (), cells.end (), [&] (const Cell& c) {
		file << 8;
		std::for_each (c.begin (), c.end (), [&] (int p) {
			file << ' ' << p;
		});
		file << std::endl;
	});
	
	file << "CELL_TYPES " << cells.size () << '\n';
	for (std::size_t i=0; i<cells.size (); i++) {
		const int cell_type = 11;
		if (format == BIN) file.write (reinterpret_cast<const char*> (&cell_type), sizeof (int));
		if (format == TXT) file << cell_type << '\n';
	}
	
	file << "CELL_DATA " << cells.size () << '\n';
	
	file << "SCALARS temperature double\n";
	file << "LOOKUP_TABLE default\n";
	for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, file, m.features->temp); });

	file << "SCALARS density double\n";
	file << "LOOKUP_TABLE default\n";
	for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, file, m.features->dens); });
		
	file << "VECTORS mass_flow double\n";
	for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, file, m.features->flow); });
		
	file << "VECTORS heat_flow double\n";
	for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, file, m.features->qflow); });
		
	file << "VECTORS pressure double\n";
	for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, file, m.features->press); });

	file << "VECTORS shear_stress double\n";
	for_each (map->all (), [&] (const Map& m) { if (m.type == GAS) write_param (format, file, m.features->shear); });

	return true;
}

void ParaView::info () const
{
	printer ().var ("Visualization program", "ParaView"); 
	printer ().var ("Format files", "ASCII"); 
}
