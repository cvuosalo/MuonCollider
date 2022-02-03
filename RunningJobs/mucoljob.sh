#!/bin/bash

mkdir jobHome
mv ilcsoft.sh jobHome
WORKDIR=`pwd`/jobHome

singularity exec -B /cvmfs --contain --home=$WORKDIR --workdir=$WORKDIR /cvmfs/cms.hep.wisc.edu/mucol/reference/mucoll_1.6_v02-07MC.sif /bin/bash ilcsoft.sh $*

mv jobHome/*.out .
