#!/bin/bash
#
# save current field output in directory 'fields#'
#
# usage: saveField.sh fieldNum

echo ""
echo "saving field $1 ..."
echo ""

saveDir="field$1"

if [ -d $saveDir ]; then
	#
	# ERROR: don't overwrite
	echo "ERROR: $saveDir directory already present"
	echo "exiting"
	echo ""
	exit
	#
else

	mkdir $saveDir
	
	mv couette0.f*    $saveDir/
	mv divcouette0.f* $saveDir/
	mv stfcouette0.f* $saveDir/
	mv hpts.out       $saveDir/
	mv logfile*       $saveDir/
	mv errfile*       $saveDir/
	mv SESSION.NAME   $saveDir/
	mv couette.sch    $saveDir/

fi
