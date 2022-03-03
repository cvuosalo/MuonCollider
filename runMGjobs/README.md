## Command for running MadGraph and Delphes on the cluster
Command:
```
python /usr/local/bin/runWiscJobs.py \
  --WorkFlow run_MG \
  --Executable=/nfs_scratch/hjia38/MG_job.sh \
  --Arguments=10 \
  --nJobs=1000 \
  --TransferInputFile=/nfs_scratch/hjia38/inputFile_tar.gz \
  --OutputDir=/nfs_scratch/hjia38 \
  --HDFSProdDir None \
  --Experiment mucol \
  --MemoryRequirements 2048
  --ExtraAttributes +BIG_MEMORY_JOB=True
```
