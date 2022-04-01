#!/bin/bash

export GBD3=/home/ccr/Source/github/geombd/

$GBD3/bin/Parameterize.py -d $GBD3/Parameters/AA.gdbp -i TMB.pdb -o TMB.pqr
$GBD3/bin/Parameterize.py -d $GBD3/Parameters/AA.gdbp -i HRP.pdb -o HRP.pqr
#$GBD3/bin/Gridder-ES -d $GBD3/Parameters/AA.gdbp -r HRP.pqr -l TMB.pqr -o HRP_TMB- -n 15 -s 0.5 -p 40 -q 0,0.2,0.1,0
#$GBD3/bin/Gridder-LJ -d $GBD3/Parameters/AA.gdbp -r HRP.pqr -l TMB.pqr -o HRP_TMB- -n 1 -p 12 -s 0.25
#$GBD3/bin/Gridder-EX -r HRP.pqr -o HRP_TMB- -n 20 -p 5 -s 0.5
