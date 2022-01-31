#!/bin/bash

mkdir jobHome
mv ilcsoft.sh jobHome

singularity exec -B jobHome:/tmp /cvmfs/cms.hep.wisc.edu/mucol/reference/mucoll_1.6_v02-07MC.sif /bin/bash /tmp/ilcsoft.sh $*

mv jobHome/*.out .

### clean dir

rm -rf jobHome

