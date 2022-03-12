## Command for running MadGraph and Delphes on the cluster
transferInputFile.tar.gz contains pre-installed package of MadGraph with Pythia8 and Delphes, and also the modified mucol collider detector card for Delphes
Command:
```
/usr/local/bin/runWiscJobs.py \
  --WorkFlow run_MG \
  --Executable=/path/to/runMG_Del.sh \
  --Arguments=10 \
  --nJobs=1000 \
  --TransferInputFile=/path/to/MG5zipDir/transferInputFile.tar.gz \
  --OutputDir=/nfs_scratch/$USER \
  --HDFSProdDir None \
  --Experiment mucol \
  --MemoryRequirements 2048
  --ExtraAttributes +BIG_MEMORY_JOB=True
```
