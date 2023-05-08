#ifndef BIBNeutralHadrons_h
#define BIBNeutralHadrons_h

/** \class BIBModule
 * 
 *  Psudo-simulation of BIB effect with given distribution of energy and position.
 * 
 *  \author H. Jia - UW-Madison
 *
 */

#include "classes/DelphesModule.h"
#include "classes/DelphesClasses.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TString.h"
#include "TFile.h"
#include "TIterator.h"
#include "TObjArray.h"

class TIterator;
class TObjArray;
class TLorentzVector;

class BIBNeutralHadrons: public DelphesModule
{
public:
  BIBNeutralHadrons();
  ~BIBNeutralHadrons();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!
  TIterator *fItOutputArray; //!

  const TObjArray *fInputArray; //!

  TObjArray *fOutputArray; //!

  TString fxHistName; //!
  TString fPositionHistName; //!
  TString fPhiHistName; //!
  TString fThetaHistName; //!
  TString fMomentumHistName; //!
  TString fPdgIDHistName; //!
  //TString fPdgEnergyHistName; //!
/*
  TH1D* fxHist; //!
  TH2D* fPositionHist; //!
  //TH1D* fPhiHist; //!
  //TH1D* fThetaHist; //!
  TH3D* fMomentumHist; //!
  TH2D* fPdgEnergyHist; //!
*/
  TString fileName; //!
  //TFile* file; //!
  //Candidate* bib; //!
   
  Int_t fNumParticles; //! 

  Double_t fDeltaR;

  //DelphesFactory *factory; //!
  //TLorentzVector bibPosition, bibMomentum; //!
  //Double_t px, py, pz, energy, mass, charge, pdgid, x, y, z, r;

  ClassDef(BIBNeutralHadrons, 1)
						  
};

#endif
