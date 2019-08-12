#!/bin/sh

BASEP="/xrootd/store/user/wjjang/ttbar"
BASEX="/xrd/store/user/wjjang/ttbar"

XRDFS="root://cms-xrdr.private.lo:2094/"
BASE=$XRDFS"/"$BASEX
XRDCP=1
#XRDCP=0

PROCESS="tt012j_bsbar_1lm_FxFx"

if [ "$#" -eq  "0" ]
  then
    echo "______________________________________________________________"
    echo "No arguments supplied, if you want to put them, don't forget an order "
    echo "====> PROCESS, IsXRDCP(0 or 1), User, SubDir in XRD"
    echo "______________________________________________________________"
    echo "BASE PATH : ${BASE} | XRDCP : ${XRDCP} | PROCESS : ${PROCESS}"
fi
if [ "$1" != "" ]
  then
    PROCESS=$1
fi
if [ "$2" != "" ]
  then
    XRDCP=$2
fi
if [ "$3" != "" ] 
  then
    BASEP="/xrootd/store/user/${3}/ttbar"
    BASEX="/xrd/store/user/${3}/ttbar"
    if [ "$4" != "" ]
      then
        BASEP="/xrootd/store/user/${3}/${4}"
        BASEX="/xrd/store/user/${3}/${4}"
    fi
    BASE=$XRDFS"/"$BASEX
fi
echo "BASE PATH : ${BASE} | XRDCP : ${XRDCP} | PROCESS : ${PROCESS}"

xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/LOG
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/GEN
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/MINIAODSIM
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/HADAOD
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/HADTRUTHAOD
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/NANOAOD
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/AODSIM
xrdfs $XRDFS mkdir -p $BASEX/$PROCESS/GENHADAOD

if [ "$HOSTNAME" = gate.sscc.uos.ac.kr ]; then
    chmod o+rwx $BASEP/$PROCESS/
    chmod o+rwx $BASEP/$PROCESS/*
fi

#for i in {51..1050}; do
for i in {51..250}; do
    condor_submit prod_gen.jds -batch-name "GEN_${PROCESS}_${i}" -append "output = log/runCMS-out_${PROCESS}_${i}.txt" -append "error = log/runCMS-error_${PROCESS}_${i}.txt" -append "log = log/runCMS-log_${PROCESS}_${i}.txt" -append "arguments = $BASE $PROCESS $i $XRDCP"
done
