#!/bin/bash

DEBUG=""
TMPFILE="/tmp/TMPMP"
TMPFILE2="/tmp/TMPMP2"

debugEcho()
{
	if [ -n "$DEBUG" ];  then
		echo "$*"
	fi
}

if  grep -q 'MODEL: ID:' $1 ; then
	MODEL="1"
	MODELDATA=$(cat "$1" | grep "MODEL: ID:" | sed 's/MODEL: ID: //' | sort -n )
else
	MODEL=""
fi

function process_NEW
{
	debugEcho "-----NEW----"
	debugEcho "$NEW"
	for LRVLINE in $NEW ; do
		debugEcho "Processing $LRVLINE"
		LRV=$(echo $LRVLINE | sed 's/.*NEW.*->.*(\(.*\))$/\1/')
		TO=$(echo $LRVLINE | sed 's/.*NEW,\(.*\)->.*$/\1/')
		debugEcho "LRV $LRV -> $TO"
		if [ ! -s "../${LRV}.csv" ]; then
			#cat $TMPFILE | grep -v "LOST,$LRV" > $TMPFILE2
			#mv $TMPFILE2 $TMPFILE
			debugEcho "Continuing"
			continue
		fi
		if [ "$TO" = "0" ]; then
			ODORALOCATIONX=0
			ODORALOCATIONY=0
			INEQ=">"
			CMPVAL=4450
		else
			ODORALOCATIONX=$(cat metadata.txt | grep OdorALocation | sed 's/.*=\(.*\),.*$/\1/')
			ODORALOCATIONY=$(cat metadata.txt | grep OdorALocation | sed 's/.*=.*, \(.*\)$/\1/')
			INEQ="<"
			CMPVAL=40
		fi
		LRVLOCATIONX=$(head -n 1 ../${LRV}.csv | cut -d',' -f 70)
		LRVLOCATIONY=$(head -n 1 ../${LRV}.csv | cut -d',' -f 71)
		DIFFXSQ="(($LRVLOCATIONX)-($ODORALOCATIONX))*(($LRVLOCATIONX)-($ODORALOCATIONX))"
		DIFFYSQ="(($LRVLOCATIONY)-($ODORALOCATIONY))*(($LRVLOCATIONY)-($ODORALOCATIONY))"
		debugEcho "$DIFFXSQ $DIFFYSQ" 
		DSTSQ=$(echo "$DIFFXSQ + $DIFFYSQ" | bc -l)
		CMP=$(echo "$DSTSQ $INEQ $CMPVAL" | bc -l )
		if [ "$CMP" = "1" ] ; then
			#cat $TMPFILE | grep -v "LOST,$LRV" > $TMPFILE2
			#mv $TMPFILE2 $TMPFILE
			#echo "$LRV->$TO: $DSTSQ"
			DATA_UP=$(echo "$DATA_NEW" | sed "s/\(^$LRV.*\),NA,\(.*\)$/\1,$TO,\2/")
			DATA_NEW="$DATA_UP"
		fi
	done
	debugEcho "-----NEW END-----"
}

function process_LOST
{
	for LRVLINE in $LOST ; do
		LRV=$(echo $LRVLINE | sed 's/.*LOST,\(.*\)->.*$/\1/')
		TO=$(echo $LRVLINE | sed 's/.*LOST,.*->\(.*\)$/\1/')
		if [ ! -s "../${LRV}.csv" ]; then
			continue
		fi
		if [ "$TO" = "0" ]; then
			ODORALOCATIONX=0
			ODORALOCATIONY=0
			INEQ=">"
			CMPVAL=4450
		else
			ODORALOCATIONX=$(cat metadata.txt | grep OdorALocation | sed 's/.*=\(.*\),.*$/\1/')
			ODORALOCATIONY=$(cat metadata.txt | grep OdorALocation | sed 's/.*=.*, \(.*\)$/\1/')
			INEQ="<"
			CMPVAL=40
		fi
		LRVLOCATIONX=$(tail -n 1 ../${LRV}.csv | cut -d',' -f 70)
		LRVLOCATIONY=$(tail -n 1 ../${LRV}.csv | cut -d',' -f 71)
		DIFFXSQ="(($LRVLOCATIONX)-($ODORALOCATIONX))*(($LRVLOCATIONX)-($ODORALOCATIONX))"
		DIFFYSQ="(($LRVLOCATIONY)-($ODORALOCATIONY))*(($LRVLOCATIONY)-($ODORALOCATIONY))"
		DSTSQ=$(echo "$DIFFXSQ + $DIFFYSQ" | bc -l)
		CMP=$(echo "$DSTSQ $INEQ $CMPVAL" | bc -l )
		#echo "$LRV->$TO: $DSTSQ"
		if [ "$CMP" = "1" ] ; then
			DATA_UP=$(echo "$DATA_NEW" | sed "s/^$LRV,,/$LRV,$TO,/")
			DATA_NEW="$DATA_UP"
		fi
	done
}

function process_data
{
	ARRAY=( 0 $(cat "$1"  | grep 'RES: OK' | tail -n 1 | sed 's/.*RES: OK Sizes: //' | sed 's/, / /g'))

	DATA=$(cat "$1"  | grep CNTTRACK | grep -v 'Check' | sed 's/.*CNTTRACK: [^,]*,//')
	DATANEW=$(echo "$DATA" |  sed "s/1->N,\([0-9]*\)->\[.*\](\([^ ]*\) \([^ ]*\) )/\1,\2,\3,,NA/" )
	DATANEW2=$(echo "$DATANEW" |  sed "s/1->N,\([0-9]*\)->\[.*\](\([^ ]*\) \([^ ]*\) \([^ ]*\) )/\1,\2,\3,\4,NA/" )
	#DATANEW3=$(echo "$DATANEW2" | sed "s/LOST,\([0-9]*\)->\(0\|-1\|-2\)$/\1,\2,,,0/")
	#DATANEW3=$(echo "$DATANEW2" | sed "s/NEW,\([0-9]*\)->\(0\|-1\|-2\)$/\1,,,,\2/")
	DATA=$(echo "$DATANEW2" | sed "s/N->1,\[\([^,]*\),\([^,]*\)\]->[0-9]*(\([0-9]*\))/\1,\3,,,NA\n\2,\3,,,NA/")

	for (( i=1; i<${#ARRAY[@]} ; i++ )); do 
		DATALINE=$(echo "$DATA" | grep "^$i," | tail -n 1) 
		if [ -z "$DATALINE" ] ; then
			echo "$i,,,,NA,${ARRAY[$i]}"
			#echo -n ""
		else
			echo "$DATALINE,${ARRAY[$i]}"
			#echo -n ""
		fi
	done
}

function findOther
{

ILINE="$1"
NOT="$2"

F2="$(echo "$ILINE" | cut -d',' -f 2)"
F3="$(echo "$ILINE" | cut -d',' -f 3)"

if [ "$F2" = "$NOT" ] ; then
	echo "$F3"
else
	echo "$F2"
fi

}

LOST=$(cat "$1" | grep CNTTRACK | grep LOST | sed 's/.*CNTTRACK: [^,]*,//')
NEW=$(cat "$1" | grep CNTTRACK | grep NEW | sed 's/.*CNTTRACK: [^,]*,//' | sed 's/NEW ,/NEW,/')
#echo "==================="
#echo "$LOST"
#echo "==================="

DATA=$(process_data $1)
debugEcho ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
debugEcho "$DATA"
debugEcho ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"

if [ -n "$MODEL" ] ; then
	for MODELLINE in $MODELDATA; do
		debugEcho "<<<<"$MODELLINE">>>>"
		#Find out reassigned number
		NEWNUM=$(echo "$MODELLINE" | cut -d">" -f1)
		RENUM=$(echo "$MODELLINE" | cut -d">" -f2)
		debugEcho "NEWNUM: $NEWNUM"
		debugEcho "RENUM: $RENUM"
		#Find line of RENUM
		RENUMLINE=$(echo "$DATA" | grep "^${RENUM}," | tail -n 1) 
		debugEcho "RENUMLINE: $RENUMLINE"
		#Find the cluster that RENUM went in
		CLUSTERNUM=$(echo "$RENUMLINE" | cut -d',' -f2)
		debugEcho "CLUSTERNUM: $CLUSTERNUM"
		#Find line of the cluster that RENUM went in
		CLUSTERLINE=$(echo "$DATA" | grep "^${CLUSTERNUM}," | grep ",${NEWNUM}," ) 
		debugEcho "CLUSTERLINE: $CLUSTERLINE"
		#Find the other member (larva or cluster) that the model collided to
		CLUSTERMATELINE=$(echo "$DATA" | grep "^[^,]*,$CLUSTERNUM," | grep -v "^${RENUM},")
		CLUSTERMATE=$(echo "$CLUSTERMATELINE" | cut -d',' -f1)
		CLUSTERMATESIZE=$(echo "$CLUSTERMATELINE" | cut -d',' -f6)
		debugEcho "CLUSTERMATE: $CLUSTERMATE"
		debugEcho "CLUSTERMATESIZE: $CLUSTERMATESIZE"
		#Find the clustermate number after collision
		OCLUSTERMATE=$(findOther "$CLUSTERLINE" "$NEWNUM")
		debugEcho "OCLUSTERMATE: $OCLUSTERMATE"
		# If cluster line has two larvae in it (for 2-2 collisions)
		if echo $CLUSTERLINE | grep -q ',2$'; then
			#REMOVE RENUMLINE and CLUSTERLINE from DATA
			debugEcho "CLUSTERLINE THERE"
			#REMOVE RENUMLINE from DATA
			debugEcho "========================"
			debugEcho "$DATA"
			debugEcho "========================"
			DATA=$(echo "$DATA" | grep -v "^$RENUMLINE")
			debugEcho "+++++++++++++++++++++++++++++++"
			debugEcho "$DATA"
			debugEcho "+++++++++++++++++++++++++++++++"
			#REMOVE CLUSTERLINE from DATA
			debugEcho "========================"
			debugEcho "$DATA"
			debugEcho "========================"
			DATA=$(echo "$DATA" | grep -v "^$CLUSTERLINE")
			debugEcho "+++++++++++++++++++++++++++++++"
			debugEcho "$DATA"
			debugEcho "+++++++++++++++++++++++++++++++"
			#REPLACE ALL INSTANCES OF NEWNUM with RENUM
			debugEcho "REPLACE ALL INSTANCES OF : $NEWNUM with $RENUM"
			DATA=$(echo "$DATA" | sed "s/\(^\|,\)$NEWNUM,/\1$RENUM,/")
			#UPDATE LOST LIST
			LOST=$(echo "$LOST" | sed "s/LOST,${NEWNUM}->/LOST,${RENUM}->/")
		elif [ -n "$CLUSTERLINE" ] ; then # The clusterline with cluster size>2
			#REMOVE RENUMLINE and CLUSTERLINE from DATA
			debugEcho "CLUSTERLINE THERE > 2"
			#REMOVE RENUMLINE from DATA
			debugEcho "========================"
			debugEcho "$DATA"
			debugEcho "========================"
			DATA=$(echo "$DATA" | grep -v "^$RENUMLINE")
			debugEcho "+++++++++++++++++++++++++++++++"
			debugEcho "$DATA"
			debugEcho "+++++++++++++++++++++++++++++++"
			#REMOVE CLUSTERLINE from DATA
			debugEcho "========================"
			debugEcho "$DATA"
			debugEcho "========================"
			DATA=$(echo "$DATA" | grep -v "^$CLUSTERLINE")
			debugEcho "+++++++++++++++++++++++++++++++"
			debugEcho "$DATA"
			debugEcho "+++++++++++++++++++++++++++++++"
			#REPLACE ALL INSTANCES OF NEWNUM with RENUM
			debugEcho "REPLACE ALL INSTANCES OF : $NEWNUM with $RENUM"
			DATA=$(echo "$DATA" | sed "s/\(^\|,\)$NEWNUM,/\1$RENUM,/")
			LOST=$(echo "$LOST" | sed "s/LOST,${NEWNUM}->/LOST,${RENUM}->/")
			#REPLACE ALL INSTANCES OF OCLUSTERMATE with CLUSTERMATE
			DATA=$(echo "$DATA" | sed "s/\(^\|,\)$OCLUSTERMATE,/\1$CLUSTERMATE,/")
			LOST=$(echo "$LOST" | sed "s/LOST,${OCLUSTERMATE}->/LOST,${CLUSTERMATE}->/")
		elif [ -z "$CLUSTERLINE" ] ; then # The clusterline was removed by the other part of the collision
			debugEcho "CLUSTERLINE NOT THERE"
			debugEcho "NEWNUM: $NEWNUM"
			debugEcho "RENUM: $RENUM"
			debugEcho "RENUMLINE: <$RENUMLINE>"
			#REMOVE RENUMLINE from DATA
			debugEcho "========================"
			debugEcho "$DATA"
			debugEcho "========================"
			DATA=$(echo "$DATA" | grep -v "^$RENUMLINE")
			debugEcho "+++++++++++++++++++++++++++++++"
			debugEcho "$DATA"
			debugEcho "+++++++++++++++++++++++++++++++"
			#REPLACE ALL INSTANCES OF NEWNUM with RENUM
			debugEcho "REPLACE ALL INSTANCES OF : $NEWNUM with $RENUM"
			DATA=$(echo "$DATA" | sed "s/\(^\|,\)$NEWNUM,/\1${RENUM},/")
			LOST=$(echo "$LOST" | sed "s/LOST,${NEWNUM}->/LOST,${RENUM}->/")
		fi
		debugEcho ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
		debugEcho "$DATA"
		debugEcho ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
	done
fi
export DATA_NEW=$(echo "$DATA" | sort -n)
process_NEW
process_LOST
echo "$DATA_NEW"
