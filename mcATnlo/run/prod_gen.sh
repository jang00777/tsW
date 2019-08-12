#!/bin/bash

date
hostname

BASE=$1
PROCESS=$2
SUFFIX=$3
XRDCP=$4

CP=cp
if [ $XRDCP -eq 1 ]; then
    CP=xrdcp
fi

echo ${CP} ${BASE} ${PROCESS} ${SUFFIX}

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cms/scratch/wjjang/Vts/tsW/mcATnlo/CMSSW_8_0_30/src #/cms/scratch/iwatson/tsW/mcATnlo/CMSSW_8_0_30/src
eval `scramv1 runtime -sh`
cd -

echo "CMSSW BASE: ${CMSSW_BASE}"
echo "SCRAM ARCH: ${SCRAM_ARCH}"

echo "-- Setup CMS environment"

TDIR=`mktemp -d`

mkdir -p $TDIR

cp /cms/scratch/wjjang/Vts/tsW/mcATnlo/run/*py $TDIR
#cp /cms/scratch/iwatson/tsW/mcATnlo/run/*py $TDIR
cd $TDIR

echo `pwd` @ `hostname` [ `date` ] > test

${CP} test ${BASE}/${PROCESS}/LOG/${SUFFIX}.start

echo "-- Finished setup"
pwd
echo "-- Running Step 1"
cmsRun ./LHEGS.py ${PROCESS} ${SUFFIX} > step1.log 2>&1
echo "-- Running Step 2"
# cmsRun ./step2_DIGI_L1_DIGI2RAW_HLT_2016.py > step2.log 2>&1
cmsRun ./step2_DIGIPREMIX_S2_DATAMIX_L1_DIGI2RAW_HLT.py > step2.log 2>&1
echo "-- Running Step 3"
cmsRun ./step3AOD_RAW2DIGI_L1Reco_RECO_EI_PAT_2016.py > step3.log 2>&1

rm step1.root
rm step2.root

${CP} step1.log ${BASE}/${PROCESS}/LOG/${SUFFIX}_step1.log
${CP} step2.log ${BASE}/${PROCESS}/LOG/${SUFFIX}_step2.log
${CP} step3.log ${BASE}/${PROCESS}/LOG/${SUFFIX}_step3.log
${CP} step3.root ${BASE}/${PROCESS}/GEN/${SUFFIX}.root
${CP} step3_inMINIAODSIM.root ${BASE}/${PROCESS}/MINIAODSIM/${SUFFIX}.root
${CP} step3_inAODSIM.root ${BASE}/${PROCESS}/AODSIM/${SUFFIX}.root

rm step3.root
rm step3_inMINIAODSIM.root
rm step3_inAODSIM.root

echo "-- Done."

rm -rf $TDIR
ls

date
