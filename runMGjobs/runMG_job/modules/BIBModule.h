/** \class BIBModule
 * 
 *  Pseudo-simulation of Beam-Induced-Backgorund effect in multi-TeV Muon Collider with given distribution of momentum, position and particle type.
 *  THIS VERSION IS NOT RECOMMENDED!!! Due to Delphes's modular structure, memory use of Delphes involving BIB simulation for multi-TeV Muon Collider is considerable. DO NOT RUN WITH MORE THAN 10^5 PARICLES WITH LESS THAN 10GB RAM!!! See BIBECal and BIBNeutralHadrons for new version of BIB pseudo-simulatioon.
 * 
 *  \author H. Jia - UW-Madison
 *
 */

#ifndef BIBModule_h
#define BIBModule_h

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

class BIBModule: public DelphesModule
{
public:
  BIBModule();
  ~BIBModule();

  void Init();
  void Process();
  void Finish();

private:

  TIterator *fItInputArray; //!
  TIterator *fItStableInputArray; //!

  const TObjArray *fInputArray; //!
  const TObjArray *fStableInputArray; //!

  TObjArray *fOutputArray; //!
  TObjArray *fStableOutputArray; //!

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

  //DelphesFactory *factory; //!
  //TLorentzVector bibPosition, bibMomentum; //!
  //Double_t px, py, pz, energy, mass, charge, pdgid, x, y, z, r;

  ClassDef(BIBModule, 1)
						  
};

#endif
