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

echo "CMSSW BASE: ${CMSSW_BASE}"
echo "SCRAM ARCH: ${SCRAM_ARCH}"

echo "-- Setup CMS environment"

TDIR=`mktemp -d`

mkdir -p $TDIR

cp /cms/scratch/wjjang/Vts/tsW/mcATnlo/run/*py $TDIR
cd $TDIR

cd /cms/ldap_home/wjjang/production_nanoAOD_CMSSW_9_4_9/src
eval `scramv1 runtime -sh`
eval `scramv1 runtime -sh`
cd -

echo `pwd`
echo "-- Run NANOAOD production"
cmsRun ./run2_2016MC_NANO.py ${PROCESS}/MINIAODSIM/${SUFFIX}.root > ${SUFFIX}.nano.log #2>&1
${CP} nano*root ${BASE}/${PROCESS}/NANOAOD/${SUFFIX}.root
${CP} ${SUFFIX}.nano.log ${BASE}/${PROCESS}/LOG/${SUFFIX}.nano.log
rm nano*root

#cmsRun /cms/scratch/iwatson/tsW/mcATnlo/run/run2_2016MC_HADAOD.py ${PROCESS}/AODSIM/${SUFFIX}.root > ${SUFFIX}.had.log 2>&1
#cmsRun /cms/ldap_home/wjjang/Vts/tsW/mcATnlo/run/run2_2016MC_HADAOD.py ${PROCESS}/AODSIM/${SUFFIX}.root > ${SUFFIX}.had.log 2>&1
#cmsRun /cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/nanoAOD/prod/run2_2016MC_HADAOD.py ${PROCESS}/AODSIM/${SUFFIX}.root > ${SUFFIX}.had.log 2>&1
#${CP} nanoAOD_AOD.root ${BASE}/${PROCESS}/HADAOD/${SUFFIX}.root
#${CP} ${SUFFIX}.had.log ${BASE}/${PROCESS}/LOG/${SUFFIX}.had.log
#rm nanoAOD_AOD.root

#cmsRun /cms/scratch/iwatson/tsW/mcATnlo/run/hadAOD.py ${PROCESS}/GEN/${SUFFIX}.root > ${SUFFIX}.hadt.log 2>&1
#cmsRun /cms/ldap_home/wjjang/Vts/tsW/mcATnlo/run/hadAOD.py ${PROCESS}/GEN/${SUFFIX}.root > ${SUFFIX}.hadt.log 2>&1
#cmsRun /cms/ldap_home/wjjang/wj_nanoAOD_CMSSW_9_4_4/src/nano/nanoAOD/prod/hadAOD.py ${PROCESS}/GEN/${SUFFIX}.root > ${SUFFIX}.hadt.log 2>&1
#${CP} hadAOD.root ${BASE}/${PROCESS}/HADTRUTHAOD/${SUFFIX}.root
#${CP} ${SUFFIX}.hadt.log ${BASE}/${PROCESS}/LOG/${SUFFIX}.hadt.log
#rm hadAOD.root

echo "-- Run genHadronAOD production"
cmsRun ./runGenHadronProducer.py ${PROCESS}/GEN/${SUFFIX}.root > ${SUFFIX}.ghad.log #2>&1
${CP} genHadronAOD.root ${BASE}/${PROCESS}/HADTRUTHAOD/${SUFFIX}.root
${CP} ${SUFFIX}.ghad.log ${BASE}/${PROCESS}/LOG/${SUFFIX}.ghad.log
rm genHadronAOD.root

echo "-- Done."

rm -rf $TDIR
ls

date
