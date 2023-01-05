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
  fNumParticles = GetInt("NumParticles", 30000000);

  fPositionHistName = GetString("PositionHistName", "xyz");
  fAngPosHistName = GetString("AndPosHistName", "thetaphi");
  fEnergyHistName = GetString("EnergyHistName", "energy");
  fMomentumHistName = GetString("EnergyHistName", "pxpypz");
  fPdgIDHistName = GetString("PdgIDHistName", "pdgid");
  
  fileName = GetString("FileName", "file");

  fInputArray = ImportArray(GetString("InputArray", "Delphes/allParticles"));
  fStableInputArray = ImportArray(GetString("StableInputArray", "Delphes/stableParticles"));
  fItInputArray = fInputArray->MakeIterator();
  fItStableInputArray = fStableInputArray->MakeIterator();

  fOutputArray = ExportArray(GetString("OutputArray", "allParticles"));
  fStableOutputArray = ExportArray(GetString("StableOutputArray", "stableParticles"));

  file = TFile::Open(fileName);
  
  fPosistionHist = (TH3F*) file->Get(fPositionHistName);
  fAngPosHist = (TH2F*) file->Get(fAngPosHistName);
  fEnergyHist = (TH1F*) file->Get(fEnergyHistName);
  fMomentumHist = (TH3F*) file->Get(fPositionHistName);
  fPdgIDHist = (TH1F*) file->Get(fPdgIDHistName);
}

void BIBModule::Finish()
{
  if(fItInputArray) delete fItInputArray; 
  if(fItStableInputArray) delete fItStableInputArray; 
}

void Process()
{
  GenParticle *original, *bib;
  Float_t px, py, pz, energy, eta, phi, m, charge, pdgid, x, y, z;

  if (!fPositionHist) return;
  if (!fAngPosHist) return;
  if (!fEnergyHist) return;
  if (!fMomentumHist) return;
  if (!fPdgIDHist) return;

  fItInputArray->Reset();
  fItStableInputArray->Reset();

  while((original = static_cast<GenParticle *>(fItInputArray->Next())))
  {
    original = static_cast<GenParticle *>(original->Clone());
    fOutputArray->Add(original);
  }


  while((original = static_cast<GenParticle *>(fItStableInputArray->Next())))
  {
    original = static_cast<GenParticle *>(original->Clone());
    fStableOutputArray->Add(original);
  }

  for (UInt i = 0; i < fNumParticles; i++)
  {
    fEnergyHist -> GetRandom(energy);
    fAngPosHist -> GetRandom2(eta, phi);
    fPositionHist -> GetRandom3(x, y, z);
    fMomentumHist -> GetRandom3(px, py, pz);
    fPdgIDHist -> GetRandom(pdgid);
    
    if(energy <= 0.0) continue;
    x = bib.X;
    y = bib.Y;
    z = bib.Z;
    eta = bib.Eta;
    phi = bib.Phi;
    energy = bib.E;
    px = bib.Px;
    py = bib.Py;
    pz = bib.Pz;

    pdgid = bib.PID; 

    //assign mass and charge by pdgID here

    fOutputArray->Add(bib);
    fStableOutputArray->Add(bib);
  }
}
