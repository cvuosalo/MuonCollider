#!/bin/bash

mkdir jobHome$$
cp ilcsoft.sh jobHome$$

singularity exec -B jobHome$$:$HOME  mucoll_1.6_v02-07MC.sif /bin/bash $HOME/ilcsoft.sh $*
