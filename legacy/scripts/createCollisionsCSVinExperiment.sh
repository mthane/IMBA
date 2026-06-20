#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$1"
for i in $(find "$DIR" -iname stdout.log) ; do
	echo "IN $(dirname $i)"
	LOGDIR="$(dirname "$i")"
	LOGFILE="$(basename "$i")"
	( cd "$LOGDIR" && bash "$SCRIPT_DIR/columnizeColisionGraphWModelBorders.sh" "$LOGFILE" ) > "$LOGDIR"/collisions.csv
	if [ -n "$(find "$LOGDIR" -type f -size 1c -iname "collisions.csv")" ] ; then
		echo "EMPTY COLLISIONS FILE!!"
		rm "$LOGDIR"/collisions.csv
	else
		echo "Running R analysis"
		arg="$(dirname "$(dirname "$i")")"
		Rscript "$SCRIPT_DIR/test.r" "$arg"
	fi
done
