#include <sstream>
#include <fstream>
#include <stdexcept>
#include "boost/filesystem.hpp"

#include "paraview.h"

using namespace Writers;

void Writer_paraview::prepare_files ()
{
	if (MPI_rank) return;
	assert (size.vol () > 0);
	cells = 0; points = 0;
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		cells += (*pbox)-> size ().vol ();
		points += ((*pbox)->size ()+1).vol ();
	}
	boost::filesystem::create_directory ("VTK");
}

struct Write_cells {
	const Files_format fmt;
	std::ofstream& file;
	const Int_vect s;
	const int curr;
	Write_cells (Files_format t, std::ofstream& f, Int_vect size, int c) : fmt (t), file (f), s (size), curr (c) { }
	void operator () (const Features&, Int_vect c)
	{
		const int nodes = 8;
		const int p0 = c.z*s.y*s.x + c.y*s.x+c.x + curr;
		const int p4 = p0 + s.x*s.y;
		if (fmt == BIN) {
			const int data[nodes+1] ={nodes, p0, p0+1, p0+s.x, p0+s.x+1, p4, p4+1, p4+s.x, p4+s.x+1};
			file.write (reinterpret_cast<const char*> (&data), (nodes+1)*sizeof(int));
		}
		if (fmt == TXT) file << nodes << ' ' << p0 << ' ' << p0+1 << ' ' << p0+s.x << ' ' << p0+s.x+1 << ' ' 
			<< p4 << ' ' << p4+1 << ' ' << p4+s.x << ' ' << p4+s.x+1 << std::endl;
	}
};

bool Writer_paraview::write_result (int time)
{
	const Files_format format = TXT;
	macroparameters ();
	if (MPI_rank) return true;
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
	
	file << "POINTS " << points << " double\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for (int n=0; n<=(*pbox)->size ().z; n++)
			for (int m=0; m<=(*pbox)->size ().y; m++)
				for (int l=0; l<=(*pbox)->size ().x; l++) 
					write_param (format, file, Real_vect ((*pbox)->coord () + Int_vect (l,m,n)) * Box::H);
					 
	file << "CELLS " << cells << ' ' << 9*cells << '\n';
	int current = 0;
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++) {
		Int_vect full_size = (*pbox)->size () + 1;
		for_each_index ((*pbox)->features->all (), Write_cells (format, file, full_size, current));
		current += full_size.vol ();
	}
				
	file << "CELL_TYPES " << cells << '\n';
	for (int i=0; i<cells; i++) {
		const int cell_type = 11;
		if (format == BIN) file.write (reinterpret_cast<const char*> (&cell_type), sizeof (int));
		if (format == TXT) file << cell_type << '\n';
	}
	
	file << "CELL_DATA " << cells << '\n';
	
	file << "SCALARS temperature double\n";
	file << "LOOKUP_TABLE default\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for_each ((*pbox)->features->all (), [format, &file] (const Features& feat) { write_param (format, file, feat.temp); });	

	file << "SCALARS density double\n";
	file << "LOOKUP_TABLE default\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for_each ((*pbox)->features->all (), [format, &file] (const Features& feat) { write_param (format, file, feat.dens); });	
		
	file << "VECTORS mass_flow double\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for_each ((*pbox)->features->all (), [format, &file] (const Features& feat) { write_param (format, file, feat.flow); });	
		
	file << "VECTORS heat_flow double\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for_each ((*pbox)->features->all (), [format, &file] (const Features& feat) { write_param (format, file, feat.qflow); });	
		
	file << "VECTORS pressure double\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for_each ((*pbox)->features->all (), [format, &file] (const Features& feat) { write_param (format, file, feat.press); });	

	file << "VECTORS shear_stress double\n";
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		for_each ((*pbox)->features->all (), [format, &file] (const Features& feat) { write_param (format, file, feat.shear); });	

	return true;
}
