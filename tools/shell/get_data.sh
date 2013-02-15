#!/bin/bash

[[ -n $1 ]] && path=$1 || path=`pwd`

[[ -f $path/result ]] || extrapolate.sh kappa.py $path > $path/result

dirs="$(ls -p $path | grep /)
$(grep -oE '([A-Z,a-z]+[0-9,\.]+/|)' $path/result)"

for d in $(echo "$dirs" | sort | uniq -u); do
	printf "Adding $d\t..."
	extrapolate.sh kappa.py $path/$d >> result
	echo "done."
done

(
	cd $path
	grep -oE '([A-Z,a-z]+[0-9,\.]+|[0-9]+\..*[0-9]+)' result | xargs -L3 | sed 's/^[A-Z,a-z]*//' | sort -n
)
