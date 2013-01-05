#!/bin/bash
# script for replacing xyz2i[][][] to xyz2i(,,)

for f in $(ls | grep pp); do
	sed -i 's/xyz2i\(.\)\[\(.*\)\]\[\(.*\)\]\[\(.*\)\]/xyz2i\1(\2, \3, \4)/g' $f
done
