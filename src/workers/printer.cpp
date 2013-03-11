#include <iomanip>

#include "printer.h"
#include "manager.h"

using namespace std;

Printer::Printer ()
{
	MPI_Comm_rank (MPI_COMM_WORLD, &MPI_rank);
}

const int string_length = 60;

void Printer::title (const string& str) const
{
	if (MPI_rank) return;
	int len = (string_length+static_cast<int> (str.size ()))/2;
	cout << setw (len) << setfill ('-') << str << setw (string_length-len) << "" << setfill (' ') << endl;
}

void Printer::log (const string& str) const
{
	if (MPI_rank) return;
	cout << str << flush;
}

void Printer::task (const string& str) const
{
	if (MPI_rank) return;
	cout << left << setw (string_length-6) << setfill ('.') << str << setfill (' ') << right << flush;
}

void Printer::result (bool res) const
{
	if (MPI_rank) return;
	if (res) cout << "..[OK]"; else cout << "[FAIL]";
	cout << endl << flush;
}

void Printer::boxes (const Boxes& boxes) const
{
	if (MPI_rank) return;
	Boxes sort_boxes;
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		sort_boxes.insert (*pbox);
	const int len = 51;
	const int tab = (string_length - len)/2;
	cout << setw (tab) << "" << setw ((len+14)/2) << "Table of boxes" << std::endl;
	cout << setw (tab) << "" << setw (len) << setfill ('~') << "" << endl << setfill (' ');
	cout << setw (tab) << "" << '|' << setw (3) << "#|" << setw (16) << "size     |" << setw (16) << "coord     |" << setw (8) << "cells |" << setw (7) << "rank |" << endl;
	cout << setw (tab) << "" << setw (len) << setfill ('~') << "" << endl << setfill (' ');
	int c = 1;
	for (Boxes::const_iterator pbox = sort_boxes.begin (); pbox != sort_boxes.end (); pbox++) {
		Box* box = *pbox;
		cout << setw (tab) << "";
		cout << '|' << setw (2) << c++;
		cout << "| " << box->size ();
		cout << " | " << box->coord () << " |";
		cout << setw (7) << box->size ().vol () << '|';
		cout << setw (6);
		if (box->MPI_rank () == -1) cout << ""; else cout << box->MPI_rank ();
		cout << '|' << endl;
	}
	cout << setw (tab) << "" << setw (len) << setfill ('~') << "" << endl << setfill (' ');
}

void Printer::MPI_ranks (const Boxes& boxes) const
{
	if (MPI_rank) return;
	int num;
	MPI_Comm_size (MPI_COMM_WORLD, &num);
	std::vector<int> ranks (num, 0);
	for (BI pbox = boxes.begin (); pbox != boxes.end (); pbox++)
		ranks[(*pbox)->MPI_rank ()] += (*pbox)->size ().vol ();
	const int len = 14;
	const int tab = (string_length - len)/2;
	cout << setw (tab-2) << "" << setw ((len-2)/2) << "Table of MPI_ranks" << std::endl;
	cout << setw (tab) << "" << setw (len) << setfill ('~') << "" << endl << setfill (' ');
	cout << setw (tab) << "" << '|' << setw (3) << "rank|" << setw (8) << "cells |" << endl;
	cout << setw (tab) << "" << setw (len) << setfill ('~') << "" << endl << setfill (' ');
	for (int i=0; i<num; i++) {
		cout << setw (tab) << "";
		cout << '|' << setw (4) << i << '|';
		cout << setw (7) << ranks[i] << '|' << endl;
	}
	cout << setw (tab) << "" << setw (len) << setfill ('~') << "" << endl << setfill (' ');
}

