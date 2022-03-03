#!/bin/bash

mkdir jobHome
mv madgraph.sh jobHome
WORKDIR='pwd'/jobHome

#Download MadGraph 5 Package and unzip
wget https://cms-project-generators.web.cern.ch/cms-project-generators/MG5_aMC_v3.3.1.tar.gz
tar xvf MG5_aMC_v3.3.1.tar.gz

#Edit script for each MadGraph job
cd MG5_aMC_v3_3_1
cp ../vvqqz.txt ./vvqqz_${jobHome}.txt
destdir=./vvqqz_${jobHome}.txt
if [ -f "$destdir" ]
then
	echo "18232$jobHome\n0" > "$destdir"
fi

#Install Pythia8 in MadGraph 5
./bin/mg5_ac Pythia8install.txt

#Run MadGraph job
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh 
export PYTHIA8DATA=$PWD/HEPTools/pythia8/share/Pythia8/xmldoc/
./bin/mg5_aMC ./vvqqz_${jobHome}.xt
cd ../bkgvvqqz_10TeV/Events/run_01
gunzip tag_1_pythia8_events.hepmc.gz

#Move the result file under jobHome and clean the MadGraph directory
cp tag_1_pythia8_events.hepmc ../../../../tag_1_pythia8_events.hepmc
cd ../../../../
rm -r MG5_aMC_v3_3_1

#Run Delphes simulation for each job result from MadGraph
mkdir $WORKDIR/jobHome/bkgvvqqz_10TeV
~/Delphes/DelphesHepMC cards/delphes_card_MuonColliderDet_MuonInJetShort.tcl bkgvvqqz_10TeV_${jobHome}.root /nfs_scratch/bkgvvqqz_10TeV/bkgvvqqz_10TeV_${jobHome}/Events/run_01/tag_1_pythia8_events.hepmc

mv jobHome/*.out .
