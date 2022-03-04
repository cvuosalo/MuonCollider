#!/bin/bash

mkdir jobHome
mv madgraph.sh jobHome
WORKDIR='pwd'/jobHome

#Download package and unzip
tar xvf transferInputFiles.tar.gz
cp -r ../Delphes ./
cp -r ../MG5_aMC_v3_3_1 ./
cp ../vvqqz.txt ./
rm -r runMG_job
export mg5dir=$PWD/MG5_aMC_v3_3_1
export Delphesdir=$PWD/Delphes
export workdir=$PWD
export datadir=$workdir/data;
echo "Using ROOTSYS=$ROOTSYS"
echo "Using mg5dir=$mg5dir"
echo "Using Delphesdir=$Delphesdir"
echo "Using workdir=$workdir"
echo "Using datadir=$datadir"
echo "Current directory set to $datadir"
mkdir -p $datadir

#Edit script for each MadGraph job
cd $mg5dir 
cp ../vvqqz.txt $mg5dir/vvqqz_${jobHome}.txt
export textdir=$mg5dir/vvqqz_${jobHome}.txt
export random='date + "%S%M%H%S"'
echo "random" > "$textdir"

#Run MadGraph job
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh 
export PYTHIA8DATA=$PWD/HEPTools/pythia8/share/Pythia8/xmldoc/
./bin/mg5_aMC ./vvqqz_${jobHome}.xt
cd ../bkgvvqqz_10TeV/Events/run_01
gunzip tag_1_pythia8_events.hepmc.gz

#Move the result file under jobHome and clean the MadGraph directory
cp tag_1_pythia8_events.hepmc ../../../../tag_1_pythia8_events.hepmc
cd $mg5dir
cd ..
rm -r MG5_aMC_v3_3_1

#Run Delphes simulation for each job result from MadGraph
cd $Delphesdir
pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_1_6; cmsenv;popd
mkdir ./bkgvvqqz_10TeV
./DelphesHepMC cards/delphes_card_MuonColliderDet_MuonInJetShort.tcl ./bkgvvqqz_10TeV_${jobHome}.root $mg5dir/bkgvvqqz_10TeV/Events/run_01/tag_1_pythia8_events.hepmc
cp ./bkgvvqqz_10TeV_${jobHome}.root ../
rm -r Delphes

mv jobHome/*.out .
