# Use singularity image to run make, i.e.
# singularity ~/Images/Madgraph.img make run_tsW_dilep

wj_jetAnalyser_cms : wj_jetAnalyser_cms.cxx
	g++ -g -std=c++1y -o wj_jetAnalyser_cms wj_jetAnalyser_cms.cxx -L/code/MG5_aMC_v2_6_0/Delphes/ -I../../install/include `root-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/Delphes/ -I/code/MG5_aMC_v2_6_0/Delphes/external -lEG -lDelphes

wj_jetAnalyser : wj_jetAnalyser.cxx
	g++ -g -std=c++1y -o wj_jetAnalyser wj_jetAnalyser.cxx -L/code/MG5_aMC_v2_6_0/Delphes/ -I../../install/include `root-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/Delphes/ -I/code/MG5_aMC_v2_6_0/Delphes/external -lEG -lDelphes

wj_analysis : wj_analysis.cxx
	g++ -g -std=c++1y -o wj_analysis wj_analysis.cxx -L/code/MG5_aMC_v2_6_0/Delphes/ -I../../install/include `root-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/Delphes/ -I/code/MG5_aMC_v2_6_0/Delphes/external -lEG -lDelphes

wj_analysis_track : wj_analysis_track.cxx
	g++ -g -std=c++1y -o wj_analysis_track wj_analysis_track.cxx -L/code/MG5_aMC_v2_6_0/Delphes/ -I../../install/include `root-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/Delphes/ -I/code/MG5_aMC_v2_6_0/Delphes/external -lEG -lDelphes


pythia_tqW : pythia_tqW.cc
#	g++ pythia_tqW.cc -o pythia_tqW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl
	g++ pythia_tqW.cc -std=c++11 -o pythia_tqW  -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include  -I/code/MG5_aMC_v2_6_0/HEPTools/pythia8//include -L/code/MG5_aMC_v2_6_0/HEPTools/pythia8//lib -lpythia8 -L/usr/lib -lboost_iostreams -L/code/MG5_aMC_v2_6_0/HEPTools/zlib/lib -lz -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run_s : pythia_tqW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tqW -s -n 100000 s.hepmc
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl s.root s.hepmc
run_b : pythia_tqW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tqW -b -n 100000 b.hepmc
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl b.root b.hepmc
run_h : pythia_tqW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tqW -h -n 20000 h.hepmc
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl h.root h.hepmc

run_l : pythia_tqW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tqW -l -n 10000 l.hepmc
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl l.root l.hepmc

run_ls : pythia_tqW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tqW -l -n 10 ls.hepmc
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl ls.root ls.hepmc

pythia_tsW : pythia_tsW.cc
#	g++ pythia_tsW.cc -o pythia_tsW `/code/MG5_aMC_v2_6_0/HEPTools/pythia8/bin/pythia8-config --cflags --libs` -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl
	g++ pythia_tsW.cc -o pythia_tsW  -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include  -I/code/MG5_aMC_v2_6_0/HEPTools/pythia8//include -L/code/MG5_aMC_v2_6_0/HEPTools/pythia8//lib -lpythia8 -L/usr/lib -lboost_iostreams -L/code/MG5_aMC_v2_6_0/HEPTools/zlib/lib -lz -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

pythia_tbW : pythia_tbW.cc
	g++ pythia_tbW.cc -o pythia_tbW  -I/code/MG5_aMC_v2_6_0/HEPTools/pythia8//include -L/code/MG5_aMC_v2_6_0/HEPTools/pythia8//lib -lpythia8 -L/usr/lib -lboost_iostreams -L/code/MG5_aMC_v2_6_0/HEPTools/zlib/lib -lz -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run_tsW : pythia_tsW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tsW
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl tsW.root tsW.hepmc

run_tbW : pythia_tbW
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tbW
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl tbW.root tbW.hepmc

pythia_tsW_dilep : pythia_tsW_dilep.cc
	g++ pythia_tsW_dilep.cc -o pythia_tsW_dilep -I/code/MG5_aMC_v2_6_0/HEPTools/pythia8//include -L/code/MG5_aMC_v2_6_0/HEPTools/pythia8//lib -lpythia8 -L/usr/lib -lboost_iostreams -L/code/MG5_aMC_v2_6_0/HEPTools/zlib/lib -lz -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

pythia_tbW_dilep : pythia_tbW_dilep.cc
	g++ pythia_tbW_dilep.cc -o pythia_tbW_dilep -I/code/MG5_aMC_v2_6_0/HEPTools/pythia8//include -L/code/MG5_aMC_v2_6_0/HEPTools/pythia8//lib -lpythia8 -L/usr/lib -lboost_iostreams -L/code/MG5_aMC_v2_6_0/HEPTools/zlib/lib -lz -I/code/MG5_aMC_v2_6_0/HEPTools/hepmc/include -L/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib -lHepMC -ldl

run_tsW_dilep : pythia_tsW_dilep
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tsW_dilep
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl tsW_dilep.root tsW_dilep.hepmc

run_tbW_dilep : pythia_tbW_dilep
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:$LD_LIBRARY_PATH ./pythia_tbW_dilep
	/code/MG5_aMC_v2_6_0/Delphes/DelphesHepMC delphes_card_CMS.tcl tbW_dilep.root tbW_dilep.hepmc

analyse :
	singularity exec ~/Images/Madgraph.img python analyse.py

WcsAnalysis : WcsAnalysis.cc
	g++ -g -std=c++1y -o WcsAnalysis WcsAnalysis.cc `root-config --cflags --libs` -I/usr/include/root -L/code/MG5_aMC_v2_6_0/Delphes/ -I../../install/include  -I/code/MG5_aMC_v2_6_0/Delphes/ -I/code/MG5_aMC_v2_6_0/Delphes/external -lEG -lDelphes

run_WcsAnalysis : WcsAnalysis
	LD_LIBRARY_PATH=/code/MG5_aMC_v2_6_0/HEPTools/hepmc/lib:/code/MG5_aMC_v2_6_0/Delphes:$LD_LIBRARY_PATH ./WcsAnalysis
