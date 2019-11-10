#include "computePtEtaTable.h"

#include <fstream>
#include <sstream>
#include<vector>
#include<string>

int computePtEtaTable::LoadDataFromCSV(std::string strPath) {
  Int_t i;
  
  Int_t nDim;
  
  std::ifstream ifsCSV(strPath.c_str());
  std::string strLine, strItem;
  std::vector<std::string> listLines;
  
  m_bValid = false;
  
  if ( !ifsCSV.is_open() ) return 1; // Oops!
  
  // Because we need the number of lines at the first, this keeping is needed
  while ( std::getline(ifsCSV, strLine) ) {
    listLines.push_back(strLine);
  }
  
  ifsCSV.close(); // have finished to read
  
  nDim = (Int_t)listLines.size() - 2; // The last two lines are for the SF values and SF errors, respectively
  m_listBinMulDim.clear();
  
  // Each of lines, except the list two lines, contains bins, separated by ','
  for ( i = 0 ; i < nDim ; i++ ) {
    std::istringstream issItems(listLines[ i ]);
    
    m_listBinMulDim.emplace_back();
    while ( std::getline(issItems, strItem, ',') ) m_listBinMulDim[ i ].push_back(std::stof(strItem));
  }
  
  // The second of the last lines contains SF values, separated by ','
  std::istringstream issItemsVal(listLines[ nDim ]);
  m_listVal.clear();
  while ( std::getline(issItemsVal, strItem, ',') ) m_listVal.push_back(std::stof(strItem));
  
  // The second of the last lines contains SF errors, separated by ','
  std::istringstream issItemsErr(listLines[ nDim + 1 ]);
  m_listErr.clear();
  while ( std::getline(issItemsErr, strItem, ',') ) m_listErr.push_back(std::stof(strItem));
  
  m_bValid = true;
  
  return 0;
}

double computePtEtaTable::getScaleFactor(std::vector<Double_t> listInput, int direction) {
  if (!m_bValid) return -1;
  
  // Seeking the bin in which the value is contained
  auto GetIdxRange = [](Float_t fX, std::vector<Float_t> listEdges) {
    if (listEdges[0] >= 0) fX = std::abs(fX); // E.g., pT bins, |eta| bins
    
    if (fX < listEdges[0]) return 0;
    else if (fX >= listEdges[listEdges.size()-1]) return ((Int_t)listEdges.size())-2;
    else return (int)(std::lower_bound(listEdges.begin(), listEdges.end(), fX)-listEdges.begin())-1;
  };
  
  Int_t i;
  Int_t nIdx;
  
  //if (m_bUseAbsEta) fEta = std::abs(fEta);
  
  nIdx = 0;
  
  // Finding the index in multi-dimensional binning
  for (i = 0; i < (Int_t)m_listBinMulDim.size(); i++) {
    nIdx *= (Int_t)m_listBinMulDim[i].size()-1;
    nIdx += GetIdxRange(listInput[i], m_listBinMulDim[i]);
  }
  //nIdxPt  = GetIdxRange(fPt,  m_listPtBin);
  //nIdxEta = GetIdxRange(fEta, m_listEtaBin);
  
  return m_listVal[nIdx]+direction*m_listErr[nIdx];
}

double computePtEtaTable::getScaleFactor(const TParticle& p, const int pid, const double shift) {
  int32_t nIdLep = abs(p.GetPdgCode());
  if ( nIdLep != pid ) return 1.0; // For dilepton channel team
  
  // Note that the eta_SC is saved in the 'weight'
  // see topObjectSelection::elecSelection()
  Double_t dEta = ( nIdLep != 11 ? p.Eta() : p.GetWeight() );
  
  return getScaleFactor({dEta, p.Pt()}, shift);
}

double computePtEtaTable::getScaleFactor(const TParticle& lep1, const TParticle& lep2, int direction) {
  // Note that the eta_SC is saved in the 'weight'
  // see topObjectSelection::elecSelection()
  Double_t dEta1 = ( abs(lep1.GetPdgCode()) != 11 ? lep1.Eta() : lep1.GetWeight() );
  Double_t dEta2 = ( abs(lep2.GetPdgCode()) != 11 ? lep2.Eta() : lep2.GetWeight() );
  return getScaleFactor({dEta1, lep1.Pt(), dEta2, lep2.Pt()}, direction);
}
