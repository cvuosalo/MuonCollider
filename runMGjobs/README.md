# Running MadGraph-Pythia8-Delphes on cluster
Since UW machine does not allows the built-in running-on-cluster mode, this instruction contains a specific approach to
1) Generate hard event by MadGraph,
2) Shower by Pythia 8,
3) Run fast simulation for ATLAS/CMS/C^3/MuonCollider detector's response using Delphes
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
cd runMGjobs/runMG\_job
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

Results are a large number of ".root", ".log", and ".html" files, feel free to modify "runMG\_Del.sh" to change what will be send back. Use "hadd" command to combine root files (syntax are inthe Detailed Description of "https://root.cern/doc/master/hadd\_8cxx.html").

## Debugging and other Issues
Notice that many things are hard-coding in the scripts. Please feel free to contact Kenny through email: hjia38@wisc.edu or haoyi.jia@cern.ch.
