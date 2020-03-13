//define
#define Branch_(type, name, suffix) outtr->Branch(#name, &(b_##name), #name "/" #suffix);
#define BranchI(name) Branch_(Int_t, name, I)
#define BranchF(name) Branch_(Float_t, name, F)
#define BranchO(name) Branch_(Bool_t, name, O)
#define BranchA_(type, name, size, suffix) outtr->Branch(#name, &(b_##name), #name"["#size"]/"#suffix);
#define BranchAI(name, size) BranchA_(Int_t, name, size, I);
#define BranchAF(name, size) BranchA_(Float_t, name, size, F);
#define BranchAO(name, size) BranchA_(Bool_t, name, size, O);
#define BranchVF(name) outtr->Branch(#name, "vector<float>", &(b_##name));
#define BranchVI(name) outtr->Branch(#name, "vector<int>", &(b_##name));
#define BranchVO(name) outtr->Branch(#name, "vector<bool>", &(b_##name));
#define BranchTLV(name) outtr->Branch(#name, "TLorentzVector", &(b_##name));

#define BranchP_(type, br, name, suffix) TBranch *br =  outtr->Branch(#name, &name, #name "/" #suffix);
#define BranchPI(br,name) BranchP_(Int_t, br,name, I);
#define BranchPF(br,name) BranchP_(Float_t,br, name, F);
#define BranchPO(br,name) BranchP_(Bool_t,br, name, O);

//struct for selected leptons
struct Lepton {
  TLorentzVector tlv;
  int charge;
  int pdgid;
};

struct RecoJet {
  int idx;
  Jet*  j;
  int   pdgid;
  int   gen_pdgid;
  float dr;
  float drl1;
  float drl2;
  float x;
  bool  isOverlap;
  int   isFrom;
  bool  hasHighestPt;
  bool  hasClosestLep;
};

struct RecoHad {
  TLorentzVector tlv;
  TLorentzVector dau1_tlv;
  TLorentzVector dau2_tlv;
  int  idx;
  int  dau1_idx;
  int  dau2_idx;
  int  pdgid;
  int  dau1_pdgid;
  int  dau2_pdgid;
  int  isFrom;
  bool isSelected;
  bool isJetMatched;
};

struct TruthHad {
  TLorentzVector tlv;
  TLorentzVector dau1_tlv;
  TLorentzVector dau2_tlv;
  float x;
  int   idx;
  int   pdgid;
  int   dau1_idx;
  int   dau1_pdgid;
  int   dau2_idx;
  int   dau2_pdgid;
  int   isFrom;
  bool  isFromTop;
  bool  isFromW;
};

//make out file
//auto out   = TFile::Open(outf.c_str(), "RECREATE");
//auto outtr = new TTree("event", "event");
//auto jettr = new TTree("MVA_jet", "MVA_jet");
//auto hadtr = new TTree("MVA_had", "MVA_had");

// const values

float ksMass_  = 0.49761; float lambda0Mass_  = 1.11568; float pionMass_  = 0.13967; float protonMass_  = 0.938272;
int   ksPdgId_ = 310;     int   lambda0PdgId_ = 3122;    int   pionPdgId_ = 211;     int   protonPdgId_ = 2212;


//read data
TClonesArray *gen_jets = 0, *particles = 0;
TClonesArray *electrons = 0, *muons = 0, *jets = 0, *missingET = 0;
TClonesArray *tracks = 0;


//declare global variables
std::vector<const GenParticle*> m_genQuark;
std::map<int, struct TruthHad> m_genHadron;
std::vector<struct Lepton> m_genLepton;
std::map<int, Jet*> m_selectedJet;
std::vector<RecoJet> m_matchedJet;
std::vector<RecoHad> m_recoHad;
bool m_hadPreCut = false;
std::string m_decayChannel;

//declare variable for branch
//  EventSelection()
TLorentzVector b_recoLep1, b_recoLep2;
int b_recoLep1_pdgId, b_recoLep2_pdgId;

//  MatchingGenJet()
float b_dilepton_mass;
int   b_dilepton_ch, b_channel, b_step;

float b_Ks_x,   b_Ks_x_S,   b_Ks_x_B,   b_Ks_rho,   b_Ks_rho_S,   b_Ks_rho_B,   b_Ks_d,   b_Ks_d_S,   b_Ks_d_B;
float b_lamb_x, b_lamb_x_S, b_lamb_x_B, b_lamb_rho, b_lamb_rho_S, b_lamb_rho_B, b_lamb_d, b_lamb_d_S, b_lamb_d_B;
float b_significance_S, b_significance_B;

bool b_gen_step0, b_gen_step1;
std::vector<float> b_jet_pt, b_jet_eta, b_jet_phi, b_jet_energy;
std::vector<int>   b_jet_pdgId;

std::vector<float> b_KsInJet_pt,   b_KsInJet_eta,   b_KsInJet_phi,   b_KsInJet_energy,   b_KsInJet_R,   b_KsInJet_outR,   b_KsInJet_rho, b_KsInJet_d;
std::vector<float> b_lambInJet_pt, b_lambInJet_eta, b_lambInJet_phi, b_lambInJet_energy, b_lambInJet_R, b_lambInJet_outR, b_lambInJet_rho, b_lambInJet_d;
std::vector<float> b_lepInJet_pt,  b_lepInJet_eta,  b_lepInJet_phi,  b_lepInJet_energy,  b_lepInJet_R,  b_lepInJet_outR;

std::vector<int> b_nKsInJet, b_nLambInJet, b_nLepInJet;

std::vector<float> b_jet1_diHadron_mass, b_jet2_diHadron_mass;

// FillJetTree() and FillHadTree()
bool  b_isSelectedJet,    b_hasHighestPt,     b_hasClosestLep;
float b_pt,               b_eta,              b_phi,              b_mass,         b_ptD,              b_radiusInEta, b_radiusInPhi;
int   b_pdgId,            b_cMult,            b_nMult,            b_flavorAlgo,   b_flavorPhys,       b_bTag,        b_bTagAlgo,    b_bTagPhys, b_tauTag, b_isFrom;
float b_tauWeight;

float b_had_pt,           b_had_eta,          b_had_phi,          b_had_mass,     b_had_x,            b_had_dr;
int   b_had_pdgId,        b_had_isFrom,       b_had_nMatched;
bool  b_had_isSelected,   b_had_isJetMatched;
int   b_dau1_pdgId;
float b_dau1_pt,          b_dau1_eta,         b_dau1_phi,         b_dau1_mass;
float b_dau1_D0,          b_dau1_DZ,          b_dau1_ctgTheta;    
float b_dau1_ptErr,       b_dau1_phiErr,      b_dau1_D0Err,       b_dau1_DZErr,   b_dau1_ctgThetaErr;
float b_dau1_D0Sig,       b_dau1_DZSig;
int   b_dau2_pdgId;
float b_dau2_pt,          b_dau2_eta,         b_dau2_phi,         b_dau2_mass;
float b_dau2_D0,          b_dau2_DZ,          b_dau2_ctgTheta;
float b_dau2_ptErr,       b_dau2_phiErr,      b_dau2_D0Err,       b_dau2_DZErr,   b_dau2_ctgThetaErr;
float b_dau2_D0Sig,       b_dau2_DZSig;

float b_genHad_pt,        b_genHad_eta,       b_genHad_phi,       b_genHad_mass,  b_genHad_x;
int   b_genHad_pdgId,     b_genHad_isFrom;    
bool  b_genHad_isFromTop, b_genHad_isFromW,   b_genHad_isPU;
float b_genHad_D0,        b_genHad_DZ,        b_genHad_ctgTheta;   
int   b_genDau1_pdgId;   
float b_genDau1_pt,       b_genDau1_eta,      b_genDau1_phi,      b_genDau1_mass;
float b_genDau1_D0,       b_genDau1_DZ,       b_genDau1_ctgTheta;
int   b_genDau2_pdgId;   
float b_genDau2_pt,       b_genDau2_eta,      b_genDau2_phi,      b_genDau2_mass;
float b_genDau2_D0,       b_genDau2_DZ,       b_genDau2_ctgTheta;


void ResetHadValues() {
  b_had_pt           = -99;   b_had_eta          = -99;   b_had_phi          = -99;   b_had_mass     = -99; b_had_x            = -99; b_had_dr     = -99;
  b_had_pdgId        = -99;   b_had_isFrom       = -99;   b_had_nMatched     = -99;
  b_had_isSelected   = false; b_had_isJetMatched = false;
  b_dau1_pdgId       = -99;
  b_dau1_pt          = -99;   b_dau1_eta         = -99;   b_dau1_phi         = -99;   b_dau1_mass    = -99;
  b_dau1_D0          = -99;   b_dau1_DZ          = -99;   b_dau1_ctgTheta    = -99;   
  b_dau1_ptErr       = -99;   b_dau1_phiErr      = -99;   b_dau1_D0Err       = -99;   b_dau1_DZErr   = -99; b_dau1_ctgThetaErr = -99;
  b_dau1_D0Sig       = -99;   b_dau1_DZSig       = -99;
  b_dau2_pdgId       = -99;
  b_dau2_pt          = -99;   b_dau2_eta         = -99;   b_dau2_phi         = -99;   b_dau2_mass    = -99;
  b_dau2_D0          = -99;   b_dau2_DZ          = -99;   b_dau2_ctgTheta    = -99;   
  b_dau2_ptErr       = -99;   b_dau2_phiErr      = -99;   b_dau2_D0Err       = -99;   b_dau2_DZErr   = -99; b_dau2_ctgThetaErr = -99;
  b_dau2_D0Sig       = -99;   b_dau2_DZSig       = -99;

  b_genHad_pt        = -99;   b_genHad_eta       = -99;   b_genHad_phi       = -99;   b_genHad_mass  = -99; b_genHad_x         = -99;
  b_genHad_pdgId     = -99;   b_genHad_isFrom    = -99;
  b_genHad_isFromTop = false; b_genHad_isFromW   = false; b_genHad_isPU      = false;
  b_genHad_D0        = -99;   b_genHad_DZ        = -99;   b_genHad_ctgTheta  = -99;
  b_genDau1_pdgId    = -99;
  b_genDau1_pt       = -99;   b_genDau1_eta      = -99;   b_genDau1_phi      = -99;   b_genDau1_mass = -99;
  b_genDau1_D0       = -99;   b_genDau1_DZ       = -99;   b_genDau1_ctgTheta = -99;
  b_genDau2_pdgId    = -99;
  b_genDau2_pt       = -99;   b_genDau2_eta      = -99;   b_genDau2_phi      = -99;   b_genDau2_mass = -99;
  b_genDau2_D0       = -99;   b_genDau2_DZ       = -99;   b_genDau2_ctgTheta = -99;

}

void ResetJetValues() {
  b_isSelectedJet  = false; b_hasHighestPt     = false;  b_hasClosestLep = false;
  b_pt             = -99;   b_eta              = -99;    b_phi           = -99;   b_mass        = -99; b_ptD         = -99; b_radiusInEta = -99; b_radiusInPhi = -99;
  b_pdgId          = -99;   b_cMult            = -99;    b_nMult         = -99;   b_flavorAlgo  = -99; b_flavorPhys  = -99; b_bTag        = -99; b_bTagAlgo    = -99; b_bTagPhys = -99; b_tauTag = -99; b_isFrom = -99;
  b_tauWeight      = -99;
  ResetHadValues();
}

// object cut variables
float cut_ElectronPt, cut_ElectronEta, cut_ElectronRelIso03All;
float cut_MuonPt,     cut_MuonEta,     cut_MuonRelIso04All;

float cut_VetoElectronPt, cut_VetoElectronEta, cut_VetoElectronRelIso03All;
float cut_VetoMuonPt,     cut_VetoMuonEta,     cut_VetoMuonRelIso04All;
  
float cut_GenJetPt, cut_GenJetEta, cut_GenJetConeSizeOverlap;
float cut_JetPt,    cut_JetEta,    cut_JetConeSizeOverlap;   
float cut_BJetPt,   cut_BJetEta,   cut_BJetConeSizeOverlap;  

float cut_Dilep_M, cut_MET_Pt;
int   cut_nJet, cut_nBJet;

// hadron reconstruction cut
// Most of them aren't used here (just copy & paste from nano_cff)
  // Track normalized Chi2 <
float tkChi2Cut_      = 100.;
  // Number of valid hits on track >=
float tkNHitsCut_     = 0.;
  // Pt of track >
float tkPtCut_        = 0.;
  // Track impact parameter significance >
float tkIPSigXYCut_   = 2.;
float tkIPSigZCut_    = 100000;
  // -- miscellaneous cuts --
  // POCA distance between tracks <
float tkDCACut_       = 2.;
  // cos(angleXY) between x and p of V0 candidate >
float cosThetaXYCut_  = -2.;
  // cos(angleXYZ) between x and p of V0 candidate >
float cosThetaXYZCut_ = -2.;

//declare functions
struct Lepton toLepton(Muon* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=11*(p->Charge); return l;}
struct Lepton toLepton(Electron* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=13*(p->Charge); return l;}
struct Lepton toLepton(const GenParticle* p){ struct Lepton l; l.tlv = p->P4(); l.charge = p->Charge; l.pdgid=p->PID; return l;}

void DefBranch(TTree* outtr);
void SetJetBranch(TTree* tr);
void SetHadBranch(TTree* tr);

const GenParticle* getLast(TClonesArray * particles, const GenParticle* p);
std::vector<const GenParticle*> getMlist(TClonesArray * particles, const GenParticle* p);
std::vector<float> collectHadron(std::vector<GenParticle> hadInJet, bool motherCheck);
std::vector<Double_t> cross3D(std::vector<Double_t> & a, std::vector<Double_t> & b);
Double_t DeltaPhi(Double_t phi1, Double_t phi2);
Double_t DeltaR(Double_t deta, Double_t dphi);

void ResetBranch();
void EventSelection(TH1F*, UInt_t);
void MatchingGenJet();
void HadronReconstruction();
void HadronPreselection(std::vector<RecoHad>);
void FindTruthMatchedHadron();
int  FindMatchedHadron(TLorentzVector, TString);
void FillJetTree(TClonesArray* jets, TTree*, TString matchingMethod = "pT");
void FillHadTree(TTree*);
void SetJetValues(Jet*);
void SetHadValues(TTree*, int,TLorentzVector jet_tlv = TLorentzVector(0, 0, 0, 0));
//std::vector<Jet*>
std::map<int, Jet*> JetSelection(TClonesArray* jets, std::vector<struct Lepton> recoLep);

//define functions
void DefBranch(TTree* outtr){

  BranchTLV(recoLep1);     BranchTLV(recoLep2);
  BranchI(recoLep1_pdgId); BranchI(recoLep2_pdgId);

  BranchI(step);
  BranchF(dilepton_mass);
  BranchI(dilepton_ch);
  BranchI(channel);

  BranchF(Ks_x);   BranchF(Ks_x_S);   BranchF(Ks_x_B);   BranchF(Ks_rho);   BranchF(Ks_rho_S);   BranchF(Ks_rho_B);   BranchF(Ks_d);   BranchF(Ks_d_S);   BranchF(Ks_d_B);
  BranchF(lamb_x); BranchF(lamb_x_S); BranchF(lamb_x_B); BranchF(lamb_rho); BranchF(lamb_rho_S); BranchF(lamb_rho_B); BranchF(lamb_d); BranchF(lamb_d_S); BranchF(lamb_d_B);

  BranchF(significance_S); BranchF(significance_B);

  BranchO(gen_step0);
  BranchO(gen_step1);

  BranchVF(jet_pt);       BranchVF(jet_eta);       BranchVF(jet_phi);       BranchVF(jet_energy);       BranchVI(jet_pdgId);
  BranchVF(KsInJet_pt);   BranchVF(KsInJet_eta);   BranchVF(KsInJet_phi);   BranchVF(KsInJet_energy);   BranchVF(KsInJet_R);   BranchVF(KsInJet_outR);   BranchVF(KsInJet_rho);   BranchVF(KsInJet_d); 
  BranchVF(lambInJet_pt); BranchVF(lambInJet_eta); BranchVF(lambInJet_phi); BranchVF(lambInJet_energy); BranchVF(lambInJet_R); BranchVF(lambInJet_outR); BranchVF(lambInJet_rho); BranchVF(lambInJet_d);
  BranchVF(lepInJet_pt);  BranchVF(lepInJet_eta);  BranchVF(lepInJet_phi);  BranchVF(lepInJet_energy);  BranchVF(lepInJet_R);  BranchVF(lepInJet_outR);
  BranchVI(nKsInJet);     BranchVI(nLambInJet);    BranchVI(nLepInJet);

  BranchVF(jet1_diHadron_mass);
  BranchVF(jet2_diHadron_mass);
}

void SetJetBranch(TTree* tr) {
  tr->Branch("isSelectedJet",    &b_isSelectedJet,    "isSelectedJet/O");
  tr->Branch("hasHighestPt",     &b_hasHighestPt,     "hasHighestPt/O");
  tr->Branch("hasClosestLep",    &b_hasClosestLep,    "hasClosestLep/O");
  tr->Branch("isFrom",           &b_isFrom,           "isFrom/I");
  tr->Branch("pt",               &b_pt,               "pt/F");
  tr->Branch("eta",              &b_eta,              "eta/F");
  tr->Branch("phi",              &b_phi,              "phi/F");
  tr->Branch("mass",             &b_mass,             "mass/F");
  tr->Branch("ptD",              &b_ptD ,             "ptD/F");
  tr->Branch("radiusInEta",      &b_radiusInEta,      "radiusInEta/F");
  tr->Branch("radiusInPhi",      &b_radiusInPhi,      "radiusInPhi/F");
  tr->Branch("pdgId",            &b_pdgId,            "pdgId/I");
  tr->Branch("cMult",            &b_cMult,            "cMult/I");
  tr->Branch("nMult",            &b_nMult,            "nMult/I");
  tr->Branch("flavorAlgo",       &b_flavorAlgo,       "flavorAlgo/I");
  tr->Branch("flavorPhys",       &b_flavorPhys,       "flavorPhys/I");
  tr->Branch("bTag",             &b_bTag,             "bTag/I");
  tr->Branch("bTagAlgo",         &b_bTagAlgo,         "bTagAlgo/I");
  tr->Branch("bTagPhys",         &b_bTagPhys,         "bTagPhys/I");
  tr->Branch("tauTag",           &b_tauTag,           "tauTag/I");
  tr->Branch("tauWeight",        &b_tauWeight,        "tauWeight/F");
  SetHadBranch(tr);
}

void SetHadBranch(TTree* tr){
  if (std::string(tr->GetName()).find("jet") != std::string::npos) {
    tr->Branch("had_x",            &b_had_x,            "had_x/F");
    tr->Branch("had_dr",           &b_had_dr,           "had_dr/F");
  }
  tr->Branch("had_pt",           &b_had_pt,           "had_pt/F");
  tr->Branch("had_eta",          &b_had_eta,          "had_eta/F");
  tr->Branch("had_phi",          &b_had_phi,          "had_phi/F");
  tr->Branch("had_mass",         &b_had_mass,         "had_mass/F");
  tr->Branch("had_pdgId",        &b_had_pdgId,        "had_pdgId/I");
  tr->Branch("had_isFrom",       &b_had_isFrom,       "had_isFrom/I");   
  tr->Branch("had_nMatched",     &b_had_nMatched,     "had_nMatched/I");
  tr->Branch("had_isSelected",   &b_had_isSelected,   "had_isSelected/O");
  tr->Branch("had_isJetMatched", &b_had_isJetMatched, "had_isJetMatched/O");
  tr->Branch("dau1_pdgId",       &b_dau1_pdgId,       "dau1_pdgId/I");
  tr->Branch("dau1_pt",          &b_dau1_pt,          "dau1_pt/F");
  tr->Branch("dau1_eta",         &b_dau1_eta,         "dau1_eta/F");
  tr->Branch("dau1_phi",         &b_dau1_phi,         "dau1_phi/F");
  tr->Branch("dau1_mass",        &b_dau1_mass,        "dau1_mass/F");
  tr->Branch("dau1_D0",          &b_dau1_D0,          "dau1_D0/F");   
  tr->Branch("dau1_DZ",          &b_dau1_DZ,          "dau1_DZ/F"); 
  tr->Branch("dau1_ctgTheta",    &b_dau1_ctgTheta,    "dau1_ctgTheta/F");
  tr->Branch("dau1_ptErr",       &b_dau1_ptErr,       "dau1_ptErr/F");
  tr->Branch("dau1_phiErr",      &b_dau1_phiErr,      "dau1_phiErr/F");
  tr->Branch("dau1_D0Err",       &b_dau1_D0Err,       "dau1_D0Err/F");
  tr->Branch("dau1_DZErr",       &b_dau1_DZErr,       "dau1_DZErr/F"); 
  tr->Branch("dau1_ctgThetaErr", &b_dau1_ctgThetaErr, "dau1_ctgThetaErr/F");
  tr->Branch("dau1_D0Sig",       &b_dau1_D0Sig,       "dau1_D0Sig/F");
  tr->Branch("dau1_DZSig",       &b_dau1_DZSig,       "dau1_DZSig/F");
  tr->Branch("dau2_pdgId",       &b_dau2_pdgId,       "dau2_pdgId/I");
  tr->Branch("dau2_pt",          &b_dau2_pt,          "dau2_pt/F");
  tr->Branch("dau2_eta",         &b_dau2_eta,         "dau2_eta/F");
  tr->Branch("dau2_phi",         &b_dau2_phi,         "dau2_phi/F");
  tr->Branch("dau2_mass",        &b_dau2_mass,        "dau2_mass/F");
  tr->Branch("dau2_D0",          &b_dau2_D0,          "dau2_D0/F"); 
  tr->Branch("dau2_DZ",          &b_dau2_DZ,          "dau2_DZ/F"); 
  tr->Branch("dau2_ctgTheta",    &b_dau2_ctgTheta,    "dau2_ctgTheta/F"); 
  tr->Branch("dau2_ptErr",       &b_dau2_ptErr,       "dau2_ptErr/F");
  tr->Branch("dau2_phiErr",      &b_dau2_phiErr,      "dau2_phiErr/F");
  tr->Branch("dau2_D0Err",       &b_dau2_D0Err,       "dau2_D0Err/F");
  tr->Branch("dau2_DZErr",       &b_dau2_DZErr,       "dau2_DZErr/F"); 
  tr->Branch("dau2_ctgThetaErr", &b_dau2_ctgThetaErr, "dau2_ctgThetaErr/F"); 
  tr->Branch("dau2_D0Sig",       &b_dau2_D0Sig,       "dau2_D0Sig/F");
  tr->Branch("dau2_DZSig",       &b_dau2_DZSig,       "dau2_DZSig/F");

  tr->Branch("genHad_x",         &b_genHad_x,         "genHad_x/F");
  tr->Branch("genHad_pt",        &b_genHad_pt,        "genHad_pt/F");
  tr->Branch("genHad_eta",       &b_genHad_eta,       "genHad_eta/F");
  tr->Branch("genHad_phi",       &b_genHad_phi,       "genHad_phi/F");
  tr->Branch("genHad_mass",      &b_genHad_mass,      "genHad_mass/F");
  tr->Branch("genHad_pdgId",     &b_genHad_pdgId,     "genHad_pdgId/I");
  tr->Branch("genHad_isFrom",    &b_genHad_isFrom,    "genHad_isFrom/I");
  tr->Branch("genHad_isFromTop", &b_genHad_isFromTop, "genHad_isFromTop/O");
  tr->Branch("genHad_isFromW",   &b_genHad_isFromW,   "genHad_isFromW/O");
  tr->Branch("genHad_isPU",      &b_genHad_isPU,      "genHad_isPU/O");
  tr->Branch("genHad_D0",        &b_genHad_D0,        "genHad_D0/F");
  tr->Branch("genHad_DZ",        &b_genHad_DZ,        "genHad_DZ/F"); 
  tr->Branch("genHad_ctgTheta",  &b_genHad_ctgTheta,  "genHad_ctgTheta/F");
  tr->Branch("genDau1_pdgId",    &b_genDau1_pdgId,    "genDau1_pdgId/I");
  tr->Branch("genDau1_pt",       &b_genDau1_pt,       "genDau1_pt/F");
  tr->Branch("genDau1_eta",      &b_genDau1_eta,      "genDau1_eta/F");
  tr->Branch("genDau1_phi",      &b_genDau1_phi,      "genDau1_phi/F");
  tr->Branch("genDau1_mass",     &b_genDau1_mass,     "genDau1_mass/F");
  tr->Branch("genDau1_D0",       &b_genDau1_D0,       "genDau1_D0/F");
  tr->Branch("genDau1_DZ",       &b_genDau1_DZ,       "genDau1_DZ/F"); 
  tr->Branch("genDau1_ctgTheta", &b_genDau1_ctgTheta, "genDau1_ctgTheta/F");
  tr->Branch("genDau2_pdgId",    &b_genDau2_pdgId,    "genDau2_pdgId/I");
  tr->Branch("genDau2_pt",       &b_genDau2_pt,       "genDau2_pt/F");
  tr->Branch("genDau2_eta",      &b_genDau2_eta,      "genDau2_eta/F");
  tr->Branch("genDau2_phi",      &b_genDau2_phi,      "genDau2_phi/F");
  tr->Branch("genDau2_mass",     &b_genDau2_mass,     "genDau2_mass/F");
  tr->Branch("genDau2_D0",       &b_genDau2_D0,       "genDau2_D0/F"); 
  tr->Branch("genDau2_DZ",       &b_genDau2_DZ,       "genDau2_DZ/F"); 
  tr->Branch("genDau2_ctgTheta", &b_genDau2_ctgTheta, "genDau2_ctgTheta/F");

}

void SetCutValues(UInt_t channel){
  if (channel == 1) {
    cut_ElectronPt     = 35; cut_ElectronEta     = 2.1; cut_ElectronRelIso03All     = 0.1;
    cut_MuonPt         = 26; cut_MuonEta         = 2.4; cut_MuonRelIso04All         = 0.15;
    
    cut_VetoElectronPt = 15; cut_VetoElectronEta = 2.4; cut_VetoElectronRelIso03All = 0.25;
    cut_VetoMuonPt     = 15; cut_VetoMuonEta     = 2.4; cut_VetoMuonRelIso04All     = 0.25;
    
    cut_GenJetPt       = 30; cut_GenJetEta       = 2.4; cut_GenJetConeSizeOverlap   = 0.4;
    cut_JetPt          = 30; cut_JetEta          = 2.4; cut_JetConeSizeOverlap      = 0.4;
    cut_BJetPt         = 30; cut_BJetEta         = 2.4; cut_BJetConeSizeOverlap     = 0.4;
    
    cut_Dilep_M        = 20.; cut_MET_Pt         = 40.; cut_nJet                    = 4;           cut_nBJet = 1;
  } else if (channel == 2) {
    cut_ElectronPt     = 20; cut_ElectronEta     = 2.4; cut_ElectronRelIso03All     = 0.12;
    cut_MuonPt         = 20; cut_MuonEta         = 2.4; cut_MuonRelIso04All         = 0.15; 
    
    cut_VetoElectronPt = 20; cut_VetoElectronEta = 2.4; cut_VetoElectronRelIso03All = 10000000000; 
    cut_VetoMuonPt     = 10; cut_VetoMuonEta     = 2.4; cut_VetoMuonRelIso04All     = 0.25; 
    
    cut_GenJetPt       = 30; cut_GenJetEta       = 2.4; cut_GenJetConeSizeOverlap   = 0.4;
    cut_JetPt          = 30; cut_JetEta          = 2.4; cut_JetConeSizeOverlap      = 0.4;
    cut_BJetPt         = 30; cut_BJetEta         = 2.4; cut_BJetConeSizeOverlap     = 0.4;
    
    cut_Dilep_M        = 20.; cut_MET_Pt         = 40.; cut_nJet                    = 2;           cut_nBJet = 1;
  } else {
    std::cout << " Set semiletpon (1) or dilepton (2) !!!! " << std::endl;
    return;
  }

}

const GenParticle* getLast(TClonesArray * particles, const GenParticle* p){
  auto mom = p;
  while(true){
    auto dau = (const GenParticle*)particles->At(mom->D1);
    if( abs(p->PID) != abs(dau->PID) ) break;
    mom = dau;
  }
  return mom;
}

std::vector<const GenParticle*> getMlist(TClonesArray * particles, const GenParticle* p){
  std::vector<const GenParticle*> mlst;
  auto idx = p->M1;
  if (idx == -1) return mlst;
  while(true){
    auto m = (const GenParticle*)particles->At(idx);
    mlst.push_back(m);
    idx = m->M1;
    if ( idx == -1) break;
  }
  return mlst;
}


std::vector<Double_t> cross3D(std::vector<Double_t> & a, std::vector<Double_t> & b){
  std::vector<Double_t> c = { a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] };
  return c;
}

Double_t DeltaPhi(Double_t phi1, Double_t phi2) {
  static const Double_t kPI = TMath::Pi();
  static const Double_t kTWOPI = 2*TMath::Pi();
  Double_t x = phi1 - phi2;
  if(TMath::IsNaN(x)){
    std::cerr << "DeltaPhi function called with NaN" << std::endl;
    return x;
  }
  while (x >= kPI) x -= kTWOPI;
  while (x < -kPI) x += kTWOPI;
  return x;
}

Double_t DeltaR(Double_t deta, Double_t dphi) {
  return TMath::Sqrt(deta*deta + dphi*dphi);
}

