# Running MadGraph-Pythia8-Delphes on cluster
Since UW machine does not allows the built-in running-on-cluster mode, this instruction contains a specific approach to
1) Generate hard events by MadGraph,
2) Shower by Pythia 8,
3) Run fast simulation for ATLAS/CMS/C3/MuonCollider detector's response using Delphes
4) Submit jobs to Condor in UW-Madison HEP group computing environment to produce large number of events.
## Installation
If you install MadGraph from the official site, this instruction will not work (especially for Muon Collider). There are modifications hard-code in this package.

```
git clone https://github.com/cvuosalo/MuonCollider/tree/main/runMGjobs runMGjobs
```
## Configuration
Configuration files includes a card for MadGraph command and a detector card for Delphes. In our example, we have modified Muon Collider card without VLC jets and adding anti-kt jets with flavor association, B-Tagging, and Tau-Tagging modules. You could change with your preference and change the detector card name inside script "runMG\_Del.sh" and also when preparing the package.

MadGraph configuration examples are in directory "runMGjobs/MadGraph\_configure".

Delphes detector card exmaple is "runMGjobs/runMG\_job/delphes\_card\_MuonColliderDet\_HHstudy.tcl".

Prepare package of pre-installed Madgraph with Pythia8 and Delphes:


```
cd runMGjobs/runMG_job
tar -czvf transferInputFiles.tar.gz ./delphes_card_MuonColliderDet_HHstudy.tcl MG5_aMC_v3_3_1 MuonCollider
```
## Submit jobs to Condor

If you are on machines with /cvmfs and CentOS7 (login.hep.wisc.edu), you may use ROOT from there:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh 
```

If you are on machines with /cvmfs and CentOS8 (mucol01.hep.wisc.edu), you may use ROOT from there:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos8-gcc11-opt/setup.sh
```
Create proxy for GRID access:

```
voms-proxy-init -rfc -valid 144:00 -voms cms
```

Example for generating 100k events with a single MadGraph configuration:

```
runWiscJobs.py \
  --WorkFlow MG_sig_3TeV \
  --Executable=runMG_Del.sh \
  --Arguments=mc-sig_3TeV.txt \
  --nJobs=1000 \
  --TransferInputFile=transferInputFiles.tar.gz,mc-sig_3TeV.txt\
  --OutputDir=/nfs_scratch/$USER \
  --HDFSProdDir None  \
  --Experiment mucol \
  --DiskRequirements 50000
```
Examples for submitting several different jobs are "submitMCJobs\*".

To check jobs progress:

```
condor_q -nobatch
```

or accessing "https://www.hep.wisc.edu/cms/comp/jobs/".

Results are a large number of ".root", ".log", and ".html" files, feel free to modify "runMG_Del.sh" to change what will be send back. Use "hadd" command to combine root files (syntax are in the Detailed Description of "https://root.cern/doc/master/hadd_8cxx.html").

## How to Build Your Own
This instruction uses MadGraph version 3.3.1 rather than the latest version. If you would like to change to the newest version. You could build from scratch by installing MadGraph from the official site and install LHAPDF6, Pythia8, Delphes inside MadGraph. 

For non-SM model, since the models are develop in python2 environment, it need to by convert to python3 version inside MadGraph before packing the zip file. Take HEFT as an example:

```
> convert model $PWD/models/heft
```
Then we need to change the 33rd line of the executable calling Delphes inside MadGraph "$mg5dir/Template/LO/bin/internal/run_delphes3" from:

```
gunzip --stdout $file | $delphesdir/DelphesHepMC2 ../Cards/delphes_card.dat ${run}/${tag}_delphes_events.root
```

to:

```
gunzip --stdout $file | $delphesdir/DelphesHepMC ../Cards/delphes_card_default.dat ${run}/${tag}_delphes_events.root
```

also for the 40th line from:

```
$delphesdir/DelphesHepMC2 ../Cards/delphes_card.dat  ${run}/${tag}_delphes_events.root $file
```

to:

```
$delphesdir/DelphesHepMC ../Cards/delphes_card_default.dat  ${run}/${tag}_delphes_events.root $file
```
Since the installation address are hard-coded in MadGraph, hence when sending to the Condor it need to be rewrite. Therefore, when making your own MadGraph package you would need to change the 8-9th lines of runMG_Del.sh from:

```
for file in `grep '/nfs_scratch/hjia38/runMG_job/' $mg5dir -Ilr`; do
    cat $file | sed 's|/nfs_scratch/hjia38/runMG_job|'$PWD'|g' > $file
```

to:

```
for file in `grep '<your installation address>' $mg5dir -Ilr`; do
    cat $file | sed 's|<your installation address>|'$PWD'|g' > $file
```
## BIB Simulation in Delphes
For BIB Simulation in Delphes see "runMG_job/BIBsimulation".

## Debugging and other Issues
The scripts are written specifically for Muon Collider simulation, other detectors (ATLAS/CMS/C3) compatibilities are added later without any test. Please feel free to contact Kenny Jia through email: hjia38@wisc.edu, hjia625@stanford.edu, or haoyi.jia@cern.ch. I would be happy to help with the debugging process or questions!


