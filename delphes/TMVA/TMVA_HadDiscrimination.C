#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void addJetVariable(TMVA::DataLoader *dataloader) {
  dataloader->AddVariable("pt",          'F');
  dataloader->AddVariable("eta",         'F');
  dataloader->AddVariable("phi",         'F');
  dataloader->AddVariable("mass",        'F');
  dataloader->AddVariable("cMult",       'I');
  dataloader->AddVariable("nMult",       'I');
  dataloader->AddVariable("radiusInEta", 'F');
  dataloader->AddVariable("radiusInPhi", 'F');
  dataloader->AddVariable("area := 3.14159*radiusInEta*radiusInPhi", 'F');
  //dataloader->AddVariable("ptD",         'F');
  dataloader->AddVariable("cpt1",        'F');
  dataloader->AddVariable("cpt2",        'F');
  dataloader->AddVariable("cpt3",        'F');
  dataloader->AddVariable("cptA",        'F');
  dataloader->AddVariable("npt1",        'F');
  dataloader->AddVariable("npt2",        'F');
  dataloader->AddVariable("npt3",        'F');
  dataloader->AddVariable("nptA",        'F');

  dataloader->AddVariable("c_x1 := cpt1/pt",        'F');
  dataloader->AddVariable("c_x2 := cpt2/pt",        'F');
  dataloader->AddVariable("c_x3 := cpt3/pt",        'F');
  dataloader->AddVariable("c_xA := cptA/pt",        'F');
  dataloader->AddVariable("n_x1 := npt1/pt",        'F');
  dataloader->AddVariable("n_x2 := npt2/pt",        'F');
  dataloader->AddVariable("n_x3 := npt3/pt",        'F');
  dataloader->AddVariable("n_xA := nptA/pt",        'F');


}

void addHadVariablePP(TMVA::DataLoader *dataloader) {
  //dataloader->AddVariable("had_x",            'F');
  //dataloader->AddVariable("had_dr",           'F');
  dataloader->AddVariable("had_pt",           'F');
  dataloader->AddVariable("had_eta",          'F');
  dataloader->AddVariable("had_phi",          'F');
  //dataloader->AddVariable("had_mass",         'F');
  dataloader->AddVariable("dau1_pt",          'F');
  dataloader->AddVariable("dau1_eta",         'F');
  dataloader->AddVariable("dau1_phi",         'F');
  dataloader->AddVariable("dau1_ctgTheta",    'F');
  //dataloader->AddVariable("fabs(dau1_D0)",    'F');
  //dataloader->AddVariable("fabs(dau1_DZ)",    'F');
  //dataloader->AddVariable("fabs(dau1_D0Sig)", 'F');
  //dataloader->AddVariable("fabs(dau1_DZSig)", 'F');
  dataloader->AddVariable("dau1_D0",    'F');
  dataloader->AddVariable("dau1_DZ",    'F');
  dataloader->AddVariable("dau1_D0Sig", 'F');
  dataloader->AddVariable("dau1_DZSig", 'F');
  dataloader->AddVariable("dau2_pt",          'F');
  dataloader->AddVariable("dau2_eta",         'F');
  dataloader->AddVariable("dau2_phi",         'F');
  dataloader->AddVariable("dau2_ctgTheta",    'F');
  //dataloader->AddVariable("fabs(dau2_D0)",    'F');
  //dataloader->AddVariable("fabs(dau2_DZ)",    'F');
  //dataloader->AddVariable("fabs(dau2_D0Sig)", 'F');
  //dataloader->AddVariable("fabs(dau2_DZSig)", 'F');
  dataloader->AddVariable("dau2_D0",    'F');
  dataloader->AddVariable("dau2_DZ",    'F');
  dataloader->AddVariable("dau2_D0Sig", 'F');
  dataloader->AddVariable("dau2_DZSig", 'F');

}

void addTree(TMVA::DataLoader *dataloader, TTree* sig, TTree* bkg, Double_t sigW = 1.0, Double_t bkgW = 1.0) {
  dataloader->AddSignalTree    ( sig, sigW );
  dataloader->AddBackgroundTree( bkg, bkgW );
}
void addTree(TMVA::DataLoader *dataloader, TTree* sigTrain, TTree* sigTest, TTree* bkgTrain, TTree* bkgTest, Double_t sigWTrain = 1.0, Double_t sigWTest = 1.0, Double_t bkgWTrain = 1.0, Double_t bkgWTest = 1.0) {
  dataloader->AddSignalTree    ( sigTrain, sigWTrain, "Training");
  dataloader->AddSignalTree    ( sigTest,  sigWTest,  "Test");
  dataloader->AddBackgroundTree( bkgTrain, bkgWTrain, "Training");
  dataloader->AddBackgroundTree( bkgTest,  bkgWTest,  "Test");
}

void addMethod(std::map<std::string,int> Use, TMVA::Factory *factory, TMVA::DataLoader *dataloader) {
  if (Use["MLP"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", 
      "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" ); // tanh and back-propagation
  if (Use["MLPSigmoid"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPSigmoid",
      "H:!V:NeuronType=sigmoid:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" ); // sigmoid and back-propagation
  if (Use["MLPBNN"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", 
      "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators
  if (Use["MLPBFGS"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", 
      "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" ); // BFGS without bayesian regulators
  if (Use["TMlpANN"])
    factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", 
      "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // sigmoid and BFGS

  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG", 
      "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT", 
      "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  if (Use["KNN"]) // K-Nearest Neighbour classifier (KNN)
    factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN", 
      "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

}

void setModel(std::map<std::string,int> Use,
            TMVA::Factory *factory,
            TMVA::DataLoader *dataloader,
            TTree* sigTrain,
            TTree* bkgTrain,
            float signalWeight = 1,
            float backgroundWeight = 1,
            TCut cut_sig = "abs(isFrom) == 3",
            TCut cut_bkg = "abs(isFrom) != 3",
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7,
            bool flag_jet=true,
            bool flag_had=true) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  if (flag_jet) addJetVariable(dataloader);
  if (flag_had) addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, bkgTrain, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}
void setModel(std::map<std::string,int> Use,
            TMVA::Factory *factory,
            TMVA::DataLoader *dataloader,
            TTree* sigTrain,
            TTree* sigTest,
            TTree* bkgTrain,
            TTree* bkgTest,
            float signalWeight = 1,
            float backgroundWeight = 1,
            TCut cut_sig = "abs(isFrom) == 3",
            TCut cut_bkg = "abs(isFrom) != 3",
            TString opt = "SplitMode=Random:NormMode=NumEvents:!V",
            float sig_ratio=0.7,
            float bkg_ratio=0.7,
            bool flag_jet=true,
            bool flag_had=true) {
  TString opt_t = "nTrain_Signal="+std::to_string(int(sigTrain->Draw("",cut_sig,"")*sig_ratio))+":nTrain_Background="+std::to_string(int(bkgTrain->Draw("",cut_bkg,"")*bkg_ratio))+":"+opt;
  if (flag_jet) addJetVariable(dataloader);
  if (flag_had) addHadVariablePP(dataloader);
  addTree(dataloader, sigTrain, sigTest, bkgTrain, bkgTest, signalWeight, backgroundWeight);
  dataloader->PrepareTrainingAndTestTree(cut_sig, cut_bkg, opt_t );
  addMethod(Use, factory, dataloader);
}


int TMVA_HadDiscrimination( TString myMethodList = "" )
{
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPSigmoid"]      = 0;
  Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator <-- very slow...
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["TMlpANN"]         = 0; // ROOT's own ANN

  // Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost

  Use["KNN"]             = 0;
  std::map<std::string, int> Opt;

  //
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVA_HadDiscrimination" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);
      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
        std::cout << std::endl;
        return 1;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  // Here the preparation phase begins

  // Read training and test data
  // (it is also possible to use ASCII format as s_vs_b_highest -> see TMVA Users Guide)
  TFile *dilep_bbars(0);       
  TFile *dilep_bsbar(0);       
  TFile *dilep_bbars_bsbar(0); 

  //TString base_path = "/hdfs/store/user/wjang/delphes_result/";
  TString base_path = "/home/wjang/CMSSW_9_3_9_patch1/src/tsW/delphes/delphes_result/";

  TString sample_dilep_bbars       = base_path+"sum_tt012j_bbars_2l_FxFx_elIso03_muIso04_newResForm_trackSmearing.root";
  TString sample_dilep_bsbar       = base_path+"sum_tt012j_bsbar_2l_FxFx_elIso03_muIso04_newResForm_trackSmearing.root";
  TString sample_dilep_bbars_bsbar = base_path+"sum_tt012j_bbars_bsbar_2l_FxFx_elIso03_muIso04_newResForm_trackSmearing.root";

  dilep_bbars       = TFile::Open( sample_dilep_bbars ); 
  dilep_bsbar       = TFile::Open( sample_dilep_bsbar ); 
  dilep_bbars_bsbar = TFile::Open( sample_dilep_bbars_bsbar );

  std::cout << "--- Delphes Jet Discrimination : Using input file 1 : " << dilep_bbars->GetName() << std::endl;
  std::cout << "--- Delphes Jet Discrimination : Using input file 2 : " << dilep_bsbar->GetName() << std::endl;
  std::cout << "--- Delphes Jet Discrimination : Using input file 3 : " << dilep_bbars_bsbar->GetName() << std::endl;

  // Register the training and test trees
  TTree *dilep_bbars_tree       = (TTree*)dilep_bbars->Get("MVA_had");
  TTree *dilep_bsbar_tree       = (TTree*)dilep_bsbar->Get("MVA_had");
  TTree *dilep_bbars_bsbar_tree = (TTree*)dilep_bbars_bsbar->Get("MVA_had");

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString outfileName( "./output/TMVA_HadDiscrimination.root" );
  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  TMVA::Factory *factory = new TMVA::Factory( "TMVA_HadDiscrimination", outputFile,
//                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );
  /* Hadron Events */
  TMVA::DataLoader *DL_KS_real_vs_fake     = new TMVA::DataLoader("DL_KS_real_vs_fake");
  TMVA::DataLoader *DL_Lambda_real_vs_fake = new TMVA::DataLoader("DL_Lambda_real_vs_fake");


  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 99999;
  (TMVA::gConfig().GetVariablePlotting()).fNbinsXOfROCCurve = 1000;

  // global event weights per tree (see below for setting event-wise weights) 
  Double_t signalWeight     = 1.0; Double_t backgroundWeight = 1.0;

  /* Cuts */
  TCut cut_recoKs      = "abs(had_pdgId) == 310";
  TCut cut_trueKs      = "abs(genHad_pdgId) == 310";
  TCut cut_recoLambda  = "abs(had_pdgId) == 3122";
  TCut cut_trueLambda  = "abs(genHad_pdgId) == 3122";

  TCut cut_noHad       = "had_pdgId == -99";

  TCut avoid_overflow  = "fabs(dau1_D0Sig) < 1000000 && fabs(dau2_D0Sig) < 1000000 && fabs(dau1_DZSig) < 1000000 && fabs(dau1_DZSig) < 1000000";

  /* Hadron */
  auto cut_had_real_KS = cut_recoKs + cut_trueKs + avoid_overflow;
  auto cut_had_fake_KS = cut_recoKs + !cut_trueKs + avoid_overflow;

  auto cut_had_real_Lambda = cut_recoLambda + cut_trueLambda + avoid_overflow;
  auto cut_had_fake_Lambda = cut_recoLambda + !cut_trueLambda + avoid_overflow;


  /* Args for Jet w/o hadron */
  TString opt = "SplitMode=Random:NormMode=NumEvents:!V";
  float sig_ratio=0.7;
  float bkg_ratio=0.7;
  bool flag_jet=false;
  bool flag_had=true;

  /* Hadron discrimination */
  setModel(Use, factory, DL_KS_real_vs_fake,     dilep_bbars_bsbar_tree, dilep_bbars_bsbar_tree, signalWeight, backgroundWeight, cut_had_real_KS,     cut_had_fake_KS,     opt, sig_ratio, bkg_ratio, flag_jet, flag_had);
  setModel(Use, factory, DL_Lambda_real_vs_fake, dilep_bbars_bsbar_tree, dilep_bbars_bsbar_tree, signalWeight, backgroundWeight, cut_had_real_Lambda, cut_had_fake_Lambda, opt, sig_ratio, bkg_ratio, flag_jet, flag_had);

  // For an example of the category classifier usage, see: TMVA_HadDiscriminationCategory
  //
  // --------------------------------------------------------------------------------------------------
  //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
  // STILL EXPERIMENTAL and only implemented for BDT's !
  //
  //     factory->OptimizeAllMethods("SigEffAt001","Scan");
  //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
  //
  // --------------------------------------------------------------------------------------------------

  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVA_HadDiscrimination is done!" << std::endl;

  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

  return 0;
}

int main( int argc, char** argv )
{
  // Select methods (don't look at this code - not of interest)
  TString methodList;
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(",");
    methodList += regMethod;
  }
  return TMVA_HadDiscrimination(methodList);
}










