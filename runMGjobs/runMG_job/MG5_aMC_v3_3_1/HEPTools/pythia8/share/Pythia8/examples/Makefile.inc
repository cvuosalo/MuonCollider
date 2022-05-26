# PYTHIA configuration file.
# Generated on Mon Mar 7 16:58:30 CST 2022 with the user supplied options:
# --prefix=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/pythia8/
# --with-hepmc2=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/hepmc
# --with-hepmc2-include=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/hepmc/include
# --with-gzip=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/zlib
# --with-lhapdf6=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/lhapdf6_py3
# --cxx-common='-ldl -fPIC -lstdc++ -std=c++11 -DHEPMC2HACK'

# Install directory prefixes.
PREFIX_BIN=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/pythia8//bin
PREFIX_INCLUDE=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/pythia8//include
PREFIX_LIB=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/pythia8//lib
PREFIX_SHARE=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/pythia8//share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-ldl -fPIC -lstdc++ -std=c++11 -DHEPMC2HACK -DGZIP -I/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/zlib/include
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=true
HEPMC2_CONFIG=
HEPMC2_BIN=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/hepmc/
HEPMC2_INCLUDE=-I/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/hepmc/include
HEPMC2_LIB=-L/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/hepmc/lib -Wl,-rpath,/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/hepmc/lib -lHepMC

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/lhapdf6_py3/bin/
LHAPDF6_INCLUDE=-I/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/lhapdf6_py3/include
LHAPDF6_LIB=-L/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/lhapdf6_py3/lib -Wl,-rpath,/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/lhapdf6_py3/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=false
ROOT_CONFIG=
ROOT_BIN=
ROOT_INCLUDE=
ROOT_LIB=

GZIP_USE=true
GZIP_CONFIG=
GZIP_BIN=/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/zlib/
GZIP_INCLUDE=-I/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/zlib/include
GZIP_LIB=-L/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/zlib/lib -Wl,-rpath,/nfs_scratch/hjia38/runMG_job/MG5_aMC_v3_3_1/HEPTools/zlib/lib -lz

PYTHON_USE=false
PYTHON_CONFIG=
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=

MPICH_USE=false
MPICH_CONFIG=
MPICH_BIN=
MPICH_INCLUDE=
MPICH_LIB=

HDF5_USE=false
HDF5_CONFIG=
HDF5_BIN=
HDF5_INCLUDE=
HDF5_LIB=

HIGHFIVE_USE=false
HIGHFIVE_CONFIG=
HIGHFIVE_BIN=
HIGHFIVE_INCLUDE=
HIGHFIVE_LIB=
