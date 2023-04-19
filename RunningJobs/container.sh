#!/bin/bash

WORKDIR=`pwd`

apptainer exec -B /cvmfs --contain --home=$WORKDIR --workdir=$WORKDIR /cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel7 /bin/bash $*
