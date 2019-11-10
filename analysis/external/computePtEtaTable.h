#ifndef computePtEtaTable_H
#define computePtEtaTable_H


#include <TFile.h>
#include <TH2F.h>
#include <TParticle.h>

class computePtEtaTable {
public:
  bool m_bValid;
  
  std::vector<Float_t> m_listVal;
  std::vector<Float_t> m_listErr;
  
  std::vector<std::vector<Float_t>> m_listBinMulDim;
  
  bool m_bUseAbsEta;
  
public:
  computePtEtaTable(std::string strPath) { LoadDataFromCSV(strPath); };
  ~computePtEtaTable() {};
  
  bool isValid() {return m_bValid;};
  
  int LoadDataFromCSV(std::string strPath);
  
  double getScaleFactor(std::vector<Double_t> listInput, int direction=0);
  double getScaleFactor(const TParticle& p, const int pid, const double shift = 0);
  double getScaleFactor(const TParticle& lep1, const TParticle& lep2, int direction=0);
};


#endif // computePtEtaTable_H


