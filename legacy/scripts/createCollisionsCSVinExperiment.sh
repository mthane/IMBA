#!/bin/bash
DIR="$1"
for i in $(find "$DIR" -iname stdout.log) ; do 
	echo "IN $(dirname $i)" 
	bash columnizeColisionGraphWModelBorders.sh "$i" > "$(dirname "$i")"/collisions.csv ; 
	if [ -n "$(find "$(dirname $i)" -type f -size 1c -iname "collisions.csv")" ] ; then
		# collisions is empty
		echo "EMPTY COLLISIONS FILE!!"
		rm "$(dirname $i)"/collisions.csv
	else
		echo "Running R analysis"
		arg="$(dirname $(dirname $i))"
		Rscript test.r $arg
	fi
done
