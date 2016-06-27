#!/bin/bash
#
function createVisitFile() {

   echo "NEK5000"                          > vis.nek3d
   echo "version: 1.0"                    >> vis.nek3d
   echo "filetemplate: couette%01d.f%05d" >> vis.nek3d
   echo "firsttimestep: 1"                >> vis.nek3d
   echo "numtimesteps: $1"                >> vis.nek3d

}
#
#
sN="1"
max="1000000"
one="1"
nine="9"
nnine="99"
nnnine="999"
nnnnine="9999"
nnnnnine="99999"
while [ $sN -le $max ]; do

   if [ $sN -le $nine ]; then

      if [ -f couette0.f0000$sN ]; then
         sN=`expr $sN + $one`
      else
         createVisitFile `expr $sN - $one`
         exit 0
      fi

   elif [ $sN -le $nnine ]; then

      if [ -f couette0.f000$sN ]; then
         sN=`expr $sN + $one`
      else
         createVisitFile `expr $sN - $one`
         exit 0
      fi

   elif [ $sN -le $nnnine ]; then

      if [ -f couette0.f00$sN ]; then
         sN=`expr $sN + $one`
      else
         createVisitFile `expr $sN - $one`
         exit 0
      fi

   elif [ $sN -le $nnnnine ]; then

      if [ -f couette0.f0$sN ]; then
         sN=`expr $sN + $one`
      else
         createVisitFile `expr $sN - $one`
         exit 0
      fi

   elif [ $sN -le $nnnnnine ]; then

      if [ -f couette0.f$sN ]; then
         sN=`expr $sN + $one`
      else
         createVisitFile `expr $sN - $one`
         exit 0
      fi

   else

		echo "**********************"
		echo "* Error in visnek.sh *"
		echo "* add nines.         *"
		echo "**********************"

   fi

done
