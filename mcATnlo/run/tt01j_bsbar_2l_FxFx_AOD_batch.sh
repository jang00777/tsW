#!/bin/bash

date
hostname

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/scratch/iwatson/tsW/mcATnlo/CMSSW_8_0_30/src
eval `scramv1 runtime -sh`
cd -

echo "-- Setup CMS environment"

TDIR=`mktemp -d`

mkdir -p $TDIR
cp /cms/scratch/iwatson/tsW/mcATnlo/run/*py $TDIR
cd $TDIR

echo `pwd` @ `hostname` [ `date` ] > test

xrdcp test root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/LOG/$1.start

echo "-- Finished setup"
pwd
echo "-- Running Step 1"
cmsRun ./tt01j_bsbar_2l_FxFx_LHEGS.py > step1.log 2>&1
echo "-- Running Step 2"
# cmsRun ./step2_DIGI_L1_DIGI2RAW_HLT_2016.py > step2.log 2>&1
cmsRun ./step2_DIGIPREMIX_S2_DATAMIX_L1_DIGI2RAW_HLT.py > step2.log 2>&1
echo "-- Running Step 3"
cmsRun ./step3AOD_RAW2DIGI_L1Reco_RECO_EI_PAT_2016.py > step3.log 2>&1

xrdcp step1.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/LOG/$1_step1.log
xrdcp step2.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/LOG/$1_step2.log
xrdcp step3.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/LOG/$1_step3.log
xrdcp step3.root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/GEN/$1.root
xrdcp step3_inAODSIM.root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/AODSIM/$1.root

rm step1.root
rm step2.root

cd /cms/scratch/iwatson/tsW/cms/CMSSW_9_4_6_patch1/src
eval `scramv1 runtime -sh` 2>&1 > /dev/null # This will fail
eval `scramv1 runtime -sh`
cd -

echo "-- Running NANOAOD"

cp /cms/scratch/iwatson/tsW/cms/CMSSW_9_4_6_patch1/src/nano/nanoAOD/prod/run2_2016MC_NANO.py .
sed -i.bak 's/fileNames/fileNames = cms.untracked.vstring("file:step3_inAODSIM.root"),   # /' run2_2016MC_NANO.py
cmsRun ./run2_2016MC_NANO.py > step4.log 2>&1

rm step3_inAODSIM.root

xrdcp step4.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/LOG/$1_step4.log
xrdcp run2_2016MC_NANO.root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/NANOAOD/$1.root

echo "-- Running HADAOD"

cp /cms/scratch/iwatson/tsW/cms/CMSSW_9_4_6_patch1/src/nano/nanoAOD/prod/hadAOD.py .
sed -i.bak 's/fileNames/fileNames = cms.untracked.vstring("file:step3.root"),   # /' hadAOD.py
cmsRun ./hadAOD.py > step5.log 2>&1

xrdcp step5.log root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/LOG/$1_step5.log
xrdcp nanoAOD.root root://cms-xrdr.sdfarm.kr:1094///xrd/store/user/iawatson/tt01j_bsbar_2l_FxFx_AOD/HADAOD/$1.root

rm step3.root

rm nanoAOD.root
rm run2_2016MC_NANO.root

echo "-- Done."
date