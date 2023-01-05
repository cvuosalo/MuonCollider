/** \class BIBModule
 
 * 
 
 *  Psudo-simulation of BIB effect with given distribution of energy and position.
 
 * 
 
 *  \author H. Jia - UW-Madison
 
 *
 
 */

#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesFormula.h"

#include "ExRootAnalysis/ExRootClassifier.h"
#include "ExRootAnalysis/ExRootFilter.h"
#include "ExRootAnalysis/ExRootResult.h"

#include "TRandom.h"
#include "TRandom2.h"
#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TObjArray.h"
#include "TString.h"
#include "TClonesArray.h"
#include "modules/Delphes.h"
#include "TDatabasePDG.h"

#include "TLorentzVector.h"
#include "TMath.h"

#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>

using namespace std;

//------------------------------------------------------------------------------

BIBModule::BIBModule() :
  fItInputArray(0), fItStableInputArray(0)
{
}

//------------------------------------------------------------------------------

BIBModule::~BIBModule() 
{
}

//------------------------------------------------------------------------------

void BIBModule::Init()
{
  fNumParticles = GetInt("NumParticles", 100);

  fPositionHistName = GetString("PositionHistName", "position_histogram");
  fEnergyHistName = GetString("EnergyHistName", "energy_histogram");
  
  fileName = GetString("FileName", "file");

  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fStableInputArray = ImportArray(GetString("StableInputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();
  fItStableInputArray = fStableInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "allParticles"));
  fStableOutputArray = ExportArray(GetString("StableOutputArray", "stableParticles"));

  file = TFile::Open(fileName);
  
  fEnergyHist = (TH1D*) file->Get(fEnergyHistName);

  fPosistionHist = (TH2D*) file->Get(fPositionHistName);
}

void BIBModule::Finish()
{
  if(fItInputArray) delete fItInputArray; 
  if(fItStableInputArray) delete fItStableInputArray; 
}

void Process()
{
  Candidate *original, *bib;
  Double_t pt, energy, eta, phi, m;

  if (!fEnergyHist) return;
  if (!fPositionHist) return;

  fItInputArray->Reset();
  fItStableInputArray->Reset();

  while((original = static_cast<Candidate *>(fItInputArray->Next())))
  {
    original = static_cast<Candidate *>(original->Clone());
    fOutputArray->Add(original);
  }


  while((original = static_cast<Candidate *>(fItStableInputArray->Next())))
  {
    original = static_cast<Candidate *>(original->Clone());
    fStableOutputArray->Add(original);
  }

  for (UInt i = 0; i < fNumParticles; i++)
  {
    fEnergyHist->GetRandom(energy);
    fPositionHist->GetRandom2(eta,phi);
    
    if(energy <= 0.0) continue;
    pt = TMath::Sqrt(energy*energy - m*m)/TMath::CosH(eta);
    bib->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    bib->Charge = ?;
    bib->PdgCode = ?;

    fOutputArray->Add(bib);
    fStableOutputArray->Add(bib);
  }
}
