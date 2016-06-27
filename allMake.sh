#
genbox << EOF
couette.box
EOF
#
genmap << EOF
box
0.05
EOF
#
mv box.map couette.map
#
reatore2 << EOF
box
couette
EOF
#
rm box.rea box.tmp
#
./makenek couette
