../pipeMeshNek/pipeMeshNek

reatore2 << EOF
base
pipe
EOF

genmap << EOF
pipe
0.05
EOF

rm base2d.rea base.rea

./makenek pipe
