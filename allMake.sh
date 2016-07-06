../pipeMeshNek/pipeMeshNek

../nek5/bin/reatore2 << EOF
base
pipe
EOF

../nek5/bin/genmap << EOF
pipe
0.05
EOF

rm base2d.rea base.rea

./makenek pipe
