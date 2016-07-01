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
	
	mv pipe0.f*     $saveDir/
	mv ???pipe0*    $saveDir/
	mv hpts.out     $saveDir/
	mv logfile*     $saveDir/
	mv errfile*     $saveDir/
	mv SESSION.NAME $saveDir/
	mv pipe.sch     $saveDir/

fi
