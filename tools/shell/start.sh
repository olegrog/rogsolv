#!/bin/bash

path=`pwd`

for dir in $(ls $path); do
	[[ ! -d $path/$dir ]] && continue
	(	
		cd $path/$dir
		[[ -f f.cache_00 ]] && continue
		pwd
		mpirun -np 2 ./rogsolv 500 20 20 > log 2> err &
	)
done
