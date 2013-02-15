#!/bin/bash

[[ -n $1 ]] && path=$1 || path=`pwd`

for dir in $(ls $path); do
	[[ ! -d $path/$dir ]] && continue
	(	
		cd $path/$dir
		[[ -f f.cache_00 ]] && continue
		pwd
		mpirun -np 4 ./rogsolv 500 10 10 > log 2> err &
	)
done
