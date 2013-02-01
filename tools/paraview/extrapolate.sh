#!/bin/bash

usage() {
	echo "This script extrapolates investigating value along time using gnuplot"
	echo "Usage: `basename $0` <script> <path>"
}

[[ $# -ne 2 ]] && { echo "Illegal amount of arguments."; usage; exit 1; }
script=$1
path=$2
[[ -f $script ]] || { echo "Bad script filename."; usage; exit 1; }
[[ -d $path ]] || { echo "Bad dirname."; usage; exit 1; }

tmp=$(mktemp /tmp/$(basename $0).XXXXXX)
for dir in $(find $path -name VTK); do
	echo $dir
	: > $tmp
	for f in $(ls $dir | grep vtk | sort -n); do
		$(dirname $script)/$(basename $script) $dir/$f | tail -1 | awk '{print $3}' >> $tmp		# write to tmp file investigated values
	done
	sed -i '1,10d' $tmp							# remove first 10 values
	num=$(wc -l $tmp | awk '{print $1}')
	echo "$num files were used"
	[[ $num -eq 0 ]] && continue
	# evaluate mean value
	cat $tmp | python -c 'from sys import stdin; import numpy as np; \
		a = [float(i) for i in stdin.read().split()]; print "mean = %.4f" % np.mean(a); print "sigma = %.2e" % (np.std(a)/np.sqrt('$num'))'
done

rm $tmp
