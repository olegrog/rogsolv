#!/bin/bash

usage() {
	echo -e "This script extrapolates investigating value along time"
	echo -e "Usage:\n\t`basename $0` [<options>] <script> <path>"
	echo -e "Options:"
	echo -e "\t--verbose, -v"
	echo -e "\t\tPrint results for each file"
}

while [[ $# -gt 2 ]]; do
	case $1 in
		--verbose|-v)
			verbose=1
			shift
			;;
	esac
done

[[ $# -ne 2 ]] && { echo "Illegal amount of arguments."; usage; exit 1; }
script=$1
path=$2
tools=$(dirname $0)
[[ -f "$tools/$script" ]] || { echo "Bad script filename."; usage; exit 1; }
[[ -d "$path" ]] || { echo "Bad dirname."; usage; exit 1; }

tmp=$(mktemp /tmp/$(basename $0).XXXXXX)		# temporary file with results
for dir in $(find $path -name VTK); do
	echo $dir
	: > $tmp
	i=0
	i0=10						# skip first ((i0)) values
	for f in $(ls $dir | grep vtk | sort -n); do
		i=$((i+1))
		[[ i -le i0 ]] && continue
		value=$("$tools/$script" $dir/$f | tail -1 | awk '{print $3}')
		[[ -n $verbose ]] && { echo "$f $value"; }
		echo $value >> $tmp			# write to tmp file investigated values
	done
	num=$(wc -l $tmp | awk '{print $1}')
	echo "$num files were used"
	[[ $num -eq 0 ]] && continue
	cat $tmp | python -c '
from sys import stdin;
import numpy as np;
a = [float(i) for i in stdin.read().split()];
print "mean = %.4f" % np.mean(a);
print "sigma = %.2e" % (np.std(a)/np.sqrt('$num'))
	'						# evaluate mean value
done

rm $tmp
