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
  fItInputArray(0)
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

  fHistogramName = GetString("HistogramName", "energy_histogram");

  fInputArray = ImportArray(GetString("InputArray", "ParticlePropagator/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "stableParticles"));

  file = TFile::Open("histograms.root");
  
  fEHistogram = (TH1D*) file->Get("EHistogram");

  fPositionHistogram = (TH2D*) file->Get("PositionHistogram");
}

void BIBModule::Finish()
{
  if(fItInputArray) delete fItInputArray; 
}

void Process()
{
  Candidate *original, *bib;
  Double_t pt, energy, eta, phi, m;

  if (!fEHistogram) return;
  if (!fPositionHistogram) return;

  fItInputArray->Reset();

  while((original = static_cast<Candidate *>(fItInputArray->Next())))
  {
    original = static_cast<Candidate *>(original->Clone());
    fOutputArray->Add(original);
  }

  for (UInt i = 0; i < fNumParticles; i++)
  {
    fEHistogram->GetRandom(energy);
    fPositionHistogram->GetRandom2(eta,phi);
    
    if(energy <= 0.0) continue;
    pt = TMath::Sqrt(energy*energy - m*m)/TMath::CosH(eta);
    bib->Momentum.SetPtEtaPhiE(pt, eta, phi, energy);
    bib->Charge = ?;
    bib->PdgCode = ?;

    fOutputArray->Add(bib);
  }
}
