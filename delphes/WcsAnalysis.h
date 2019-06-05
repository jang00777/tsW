#include "TClonesArray.h"
#include "TRefArray.h"
#include "TRef.h"
#include "TTree.h"

const int ks_pdg = 310;
const int kl_pdg = 130;
const int k0_pdg = 311;

const int pi_pdg = 211;
const int pi0_pdg = 111;

template <class T>
class TClonesIter {
public:
  TClonesArray *a;
  int i;

  TClonesIter(TClonesArray *_a, int _i) : a(_a), i(_i) {}
  
  T& operator*() { return *((T*) a->At(i)); }
  TClonesIter& operator++() { ++i; return *this; }
  bool operator!=(TClonesIter oth) { return oth.i != i; }
};

template <class T>
class TClones {
public:
  TClonesArray *a;
  TClonesIter<T> begin() { return TClonesIter<T>(a, 0); }
  TClonesIter<T> end() { return TClonesIter<T>(a, a->GetEntries()); }
  T* operator[](size_t i) { return (T*) a->At(i); }
  TClones(TClonesArray *_a=0) : a(_a) {}
  void FromTree(TTree *t, std::string s) { t->SetBranchAddress(s.c_str(), &a); }
};


template <class T>
class TRefIter {
public:
  TRefArray *a;
  int i;

  TRefIter(TRefArray *_a, int _i) : a(_a), i(_i) {}

  T* operator*() {
    auto obj = dynamic_cast<T*>(a->At(i));
    return obj;
  }
  TRefIter& operator++() { ++i; return *this; }
  bool operator!=(TRefIter oth) { return oth.i != i; }
};

template <class T>
class TRefs {
public:
  TRefArray *a;
  TRefIter<T> begin() { return TRefIter<T>(a, 0); }
  TRefIter<T> end() { return TRefIter<T>(a, a->GetEntries()); }
  T* operator[](size_t i) { return (T*) a->At(i); }

  template <class U>
  U* as(size_t i) { return (U*) a->At(i); }
  
  TRefs(TRefArray *_a=0) : a(_a) {}
  size_t size() { return a->GetEntries(); }
  //  void FromTree(TTree *t, std::string s) { t->SetBranchAddress(s.c_str(), &a); }
};
