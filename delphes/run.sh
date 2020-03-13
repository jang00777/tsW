#!/bin/sh

# Run using:

decay_ch=$1

DELPHES_HOME=/home/wjang/CMSSW_9_3_9_patch1/src/Delphes

export PYTHONPATH=${DELPHES_HOME}/python:$PYTHONPATH
export LD_LIBRARY_PATH=${DELPHES_HOME}:$LD_LIBRARY_PATH
export ROOT_INCLUDE_PATH=${DELPHES_HOME}/external:${ROOT_INCLUDE_PATH}


samples=("tt012j_bbars_2l_FxFx" "tt012j_bsbar_2l_FxFx" "tt012j_bsbar_1lm_FxFx" "tt012j_bsbar_1lp_FxFx" "tt012j_sbbar_1lm_FxFx" "tt012j_sbbar_1lp_FxFx") 
#samples=$2
#("tt012j_bsbar_2l_FxFx")

#for f in /home/scratch/tsW/*t*; do
#for f in /scratch/tsW/*k.root*; do
for s in ${samples[@]}; do
    if [[ ${s} == *"_2l_"* ]]; then
        decay_ch="di" 
    else
        decay_ch="semi"
    fi
    sample_full=${s}_"elIso03_muIso04_newResForm"
    for f in /home/wjang/CMSSW_9_3_9_patch1/src/DelphesRun/delphes_ttbar/${sample_full}/*.root; do
        OUT="Out_`basename $f`"
        if [ ! -d delphes_result/${sample_full} ]; then
            mkdir delphes_result/${sample_full}
        else 
            echo "Folder already exists"
        fi
        if [[ $f -nt delphes_result/${sample_full}/`basename $f` ]]; then
            echo "Processing " $f
            ./wj_analysis $f delphes_result/${sample_full}/$OUT ${decay_ch} 
        else
            echo "Not processing" $f
        fi
    done
done

#./wj_analysis /home/wjang/CMSSW_9_3_9_patch1/src/DelphesRun/delphes_ttbar/tt012j_bsbar_2l_FxFx_elIso03_muIso04/tt012j_bsbar_2l_FxFx_elIso03_muIso04_0.root delphes_result/out_tt012j_bsbar_2l_FxFx_elIso03_muIso04_0.root ${decay_ch}
#./wj_analysis /home/wjang/CMSSW_9_3_9_patch1/src/DelphesRun/delphes_ttbar/tt012j_bsbar_2l_FxFx_elIso03_muIso04/tt012j_bsbar_2l_FxFx_elIso03_muIso04_0.root delphes_result/out_test.root ${decay_ch}

#./wj_analysis /home/wjang/CMSSW_9_3_9_patch1/src/DelphesRun/test/decay_filter.root delphes_result/out_test.root ${decay_ch}
