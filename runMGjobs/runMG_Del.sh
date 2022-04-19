#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_101/x86_64-centos7-gcc11-opt/setup.sh

#Download package and unzip
tar xvf transferInputFiles.tar.gz
export mg5dir=$PWD/MG5_aMC_v3_3_1
jobfilename="$(echo $1)"
for file in `grep '/nfs_scratch/hjia38/runMG_job/' $mg5dir -Ilr`; do
    cat $file | sed 's|/nfs_scratch/hjia38/runMG_job|'$PWD'|g' > $file
done
if [[ $1 == cms-*  ]]; then
    cp $mg5dir/Delphes/cards/delphes_card_CMS.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
    jobfilename2="$(echo ${jobfilename#c*-})"
    jobname="$(echo ${jobfilename2%.*t})"
elif [[ $1 == mc-*  ]]; then
    cp $PWD/delphes_card_MuonColliderDet_MuonInJetShort.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
    cp -r $PWD/MuonCollider $mg5dir/Template/Common/Cards/
    jobfilename2="$(echo ${jobfilename#m*-})"
    jobname="$(echo ${jobfilename2%.*t})"
elif [[ $1 == c3-*  ]]; then
    cp $mg5dir/Delphes/cards/delphes_card_ILD.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
    jobfilename2="$(echo ${jobfilename#c*-})"
    jobname="$(echo ${jobfilename2%.*t})"
elif [[ $1 == atlas-*  ]]; then
    cp $mg5dir/Delphes/cards/delphes_card_ATLAS.tcl $mg5dir/Template/Common/Cards/delphes_card_default.dat
    jobfilename2="$(echo ${jobfilename#a*-})"
    jobname="$(echo ${jobfilename2%.*t})"
fi
export workdir=$PWD
echo "Using ROOTSYS=$ROOTSYS"
echo "Using mg5dir=$mg5dir"
echo "Using Delphesdir=$delphesdir"
echo "Using workdir=$workdir"
echo "Using datadir=$datadir"
echo "Current directory set to $datadir"
export PATH=$mg5dir/bin:$mg5dir/HEPTools/bin:$PATH
export LD_LIBRARY_PATH=$mg5dir/HEPTools/lib/:$mg5dir/HEPTools/lhapdf6_py3/lib/:$mg5dir/HEPTools/lhapdf6_py3/lib/python3.9/site-packages/:$mg5dir/HEPTools/hepmc/lib/:$mg5dir/HEPTools/pythia8//lib:$mg5dir/HEPTools/zlib/lib/:$LD_LIBRARY_PATH
export PYTHIA8DATA=$mg5dir/HEPTools/pythia8/share/Pythia8/xmldoc
export LHAPDF_DATA_PATH=$mg5dir/HEPTools/lhapdf6_py3/share/LHAPDF
source $mg5dir/Delphes/DelphesEnv.sh
RANDOM=$$

#Edit script for each MadGraph job
cd $mg5dir 
cp $workdir/$1 $mg5dir/mg5-configure.txt
cat $workdir/$1 | sed 's/set iseed 0/set iseed '"$RANDOM"'/g' > $mg5dir/mg5-configure.txt

#Run MadGraph job
python $mg5dir/bin/mg5_aMC mg5-configure.txt
ls -alh
cd $mg5dir/${jobname}/Events/run_01
ls -alh
echo 'using random seed '"$RANDOM"''

#Move the result file under jobHome and clean the MadGraph directory
cp $mg5dir/${jobname}/index.html $workdir/index_${jobname}_${RANDOM}.html
cp $mg5dir/${jobname}/crossx.html $workdir/crossx_${jobname}_${RANDOM}.html
cp $mg5dir/${jobname}/Events/run_01/*.root $workdir/${jobname}_${RANDOM}.root
cd $workdir
rm -r MG5_aMC_v3_3_1

