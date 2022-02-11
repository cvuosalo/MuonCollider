#!/bin/bash

source /opt/ilcsoft/muonc/init_ilcsoft.sh

git clone https://github.com/MuonColliderSoft/MuonCutil.git

WORKDIR=`pwd`

cd MuonCutil/SoftCheck

GEO="/opt/ilcsoft/muonc/detector-simulation/geometries/MuColl_v1/MuColl_v1.xml"

ddsim --compactFile ${GEO} --steeringFile sim_steer.py &> $WORKDIR/sim$$.out

Marlin --InitDD4hep_mod4.DD4hepXMLFile=${GEO} reco_steer.xml &> $WORKDIR/reco$$.out

