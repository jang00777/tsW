In submit_prod_gen_pythia_limitCylinder_noVtxSmearing.sh
you can change a output path from L3 to L8 
The numbers in L56 are seeds used for generation (default : 51 to 1050)

The script is assumed to run on KISTI (so I'm not sure the script works on gate also)

For running the script, just do :
```
./run_submit_prod_gen_pythia_limitCylinder_noVtxSmearing.sh
```

If there's no tarball in default tarball path (/cms/scratch/wjjang/Vts/tsW/mcATnlo/run/tarball/), change a line below  
```
In step1_pythia_limitCylinder_noVtxSmearing.py
    ====>
    args = cms.vstring('/cms/scratch/wjjang/Vts/tsW/mcATnlo/run/tarball/%s_slc6_amd64_gcc481_CMSSW_7_1_30_tarball.tar.xz' % sname), # L163
```



