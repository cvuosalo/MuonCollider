# MadGraph5_aMC@NLO and DELPHES
## Monte Carlo simulation with MadGraph5
### MG5 installation
```bash
wget https://cms-project-generators.web.cern.ch/cms-project-generators/MG5_aMC_v2.7.2.tar.gz
tar xvf MG5_aMC_v2.7.2.tar.gz
```
### Install pythia8 on MG5
```bash
cd MG5_aMC_v2_7_2
./bin/mg5_aMC
```
Do not recommend to upgrade here when it ask you to do so, as old version is less buggy.
```MG5
install pythia8
exit
```
### Simulation
```bash
pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_1_6; cmsenv;popd
cd MG5_aMC_v2_7_2
./bin/mg5_aMC
```
Here is a example for generating the signal for diHiggs channel, similar for background or other channel 
```MG5
import model sm
generate mu+ mu- > vm vm~ h h
output mucol-hpair_sig
launch
```
Then you will see a menu for MG5, type 1 to open up the pythia8 for format transformation to .hepmc, then type enter to next step--setting parameters:
```MG5
set mh 125.0
set wh 0.004
set ebeam1 1500.
set ebeam2 1500.
set nevents 1000
set iseed 1823211
```
The process is expected to take about 5-60 min. 
After generate the evetns, don't forget to gunzip the .hepmc.tar.gz for preparation for delphes simulation:
```bash 
cd <your output directory>/Events/run_01/
gunzip tag_1_pythia8_events.hepmc.gz
```
## Detector simulation on Delphes
### Installation of Delphes
```bash
git clone git://github.com/delphes/delphes.git Delphes
cd Delphes
```
If you are on UW machine, you need to change the 23-26 lines of Makefile from:
```
CXXFLAGS += -std=c++0x -I$(subst :, -I,$(CMSSW_FWLITE_INCLUDE_PATH))
OPT_LIBS += -L$(subst include,lib,$(subst :, -L,$(CMSSW_FWLITE_INCLUDE_PATH)))
ifneq ($(CMSSW_RELEASE_BASE),)
CXXFLAGS += -I$(CMSSW_RELEASE_BASE)/src
```
to:
```
CXXFLAGS += -std=c++17 -I$(subst :, -I,$(CMSSW_FWLITE_INCLUDE_PATH))
OPT_LIBS += -L$(subst include,lib,$(subst :, -L,$(CMSSW_FWLITE_INCLUDE_PATH)))
ifneq ($(CMSSW_BASE),)
CXXFLAGS += -I$(CMSSW_BASE)/src
```
then you could compile
```bash
pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_1_6; cmsenv;popd
make -j 8
```
However if you are on lxplus.cern.ch, you could simply run:
```bash
source  /afs/cern.ch/sw/lcg/external/gcc/4.9.3/x86_64-slc6/setup.sh
source /afs/cern.ch/sw/lcg/app/releases/ROOT/6.06.00/x86_64-slc6-gcc49-opt/root/bin/thisroot.sh
make
```
There will be some warning or errors when compile, to test if it works fine:
```bash
./DelphesHepMC
#(it just prints usage command as we did not give input - next get input and use)
```
If it shows:
```bash
$ ./DelphesHepMC
Usage: DelphesHepMC config_file output_file [input_file(s)]
 config_file - configuration file in Tcl format,
 output_file - output file in ROOT format,
 input_file(s) - input file(s) in HepMC format,
 with no input_file, or when input_file is -, read standard input.
```
Then you are ready to run simulation for your events!
### Simulation
```bash
./DelphesHepMC cards/delphes_card_MuonColliderDet.tcl <name for output root file>.root ~/MG5_aMC_v2_7_2/<your output directory>/Events/run_01/tag_1_pythia8_events.hepmc
```
Here you could also send your simulation job with nohup command to avoid failing by losing connection
```bash
    nohup ./DelphesHepMC cards/delphes_card_MuonColliderDet.tcl nfs_scratch/<your user name>/<name for output root file>.root ~/MG5_aMC_v2_7_2/<your output directory>/Events/run_01/tag_1_pythia8_events.hepmc >& job.log &
```
To check:
```bash
    ps -fu your_username
```
## Event Display
Firstly compile the display library
```bash
cd Delphes/
make display
```
Then we could call the event display macro:
```bash
pushd /cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/cmssw/CMSSW_11_1_6; cmsenv;popd
cd Delphes
root -l examples/EventDisplay.C'("cards/delphes_card_MuonColliderDet.tcl","<your_output_file>.root","ParticlePropagator","ChargedHadronTrackingEfficiency","MuonTrackingEfficiency","Ecal,Hcal")''
